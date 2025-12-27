## =======================================================
##  PACOTES
## =======================================================
suppressPackageStartupMessages({
  library(AER)
  library(systemfit)
  library(sandwich)
  library(lmtest)
  library(car)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(MASS)   # ginv
  library(openxlsx)
})

## =======================================================
##  HELPERS BÁSICOS
## =======================================================

# escolhe variável de cluster automaticamente
pick_cluster <- function(data, candidates = c("psu","uf"), fallback = NULL) {
  nm <- names(data); hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else fallback
}

# alinha vetor de cluster às observações do modelo
.align_cluster <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)
  mf  <- stats::model.frame(model); rid <- rownames(mf)
  if (!is.null(names(cluster))) return(cluster[rid])
  idx <- suppressWarnings(as.integer(rid))
  if (all(is.finite(idx))) return(cluster[idx])
  cluster
}

# model matrix segura (remove constantes/singulares)
.mm <- function(formula_or_terms, data, add_intercept = TRUE) {
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm)
  if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}

# extrai X, Z e resíduos de um ivreg
.get_XZ_u <- function(ivm) {
  X <- try(model.matrix(ivm, "regressors"),  silent = TRUE)
  Z <- try(model.matrix(ivm, "instruments"), silent = TRUE)
  if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NULL)
  list(X = X, Z = Z, u = residuals(ivm))
}

# Hansen-J robusto/cluster (Arellano) eq-a-eq
.j_overid_manual <- function(ivm, cluster = NULL, ridge = 1e-8) {
  mats <- .get_XZ_u(ivm); if (is.null(mats)) return(c(J = NA_real_, df = NA_real_, p = NA_real_))
  X <- mats$X; Z <- mats$Z; u <- as.numeric(mats$u)
  n <- NROW(Z); kz <- ncol(Z); kx <- ncol(X); df <- kz - kx
  if (!is.finite(df) || df <= 0) return(c(J = NA_real_, df = df, p = NA_real_))
  cl <- .align_cluster(ivm, cluster)
  Zu <- Z * u
  if (!is.null(cl)) {
    G  <- rowsum(Zu, group = as.factor(cl))
    S  <- crossprod(as.matrix(G)) / n
  } else {
    S  <- crossprod(Zu) / n
  }
  S <- (S + t(S)) / 2
  d <- mean(diag(S)); if (!is.finite(d) || d <= 0) d <- 1
  S_r <- S + diag(ridge * d, ncol(S))
  gbar <- colMeans(Zu)
  J    <- as.numeric(n * crossprod(gbar, solve(S_r, gbar)))
  p    <- stats::pchisq(J, df = df, lower.tail = FALSE)
  c(J = J, df = df, p = p)
}

# constrói instrumentos = exógenos do RHS + IVs excluídos (keep_excl)
.build_inst_from_fit <- function(fit_obj, keep_excl = NULL) {
  rhs_all  <- unique(fit_obj$rhs_terms)
  base_exog <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  inst <- unique(c(base_exog, keep_excl))
  inst <- inst[!is.na(inst) & nzchar(inst)]
  setdiff(inst, c("1","-1","0","+","~","|"))
}

# pool de IVs excluídos possíveis (inst - exógenos do RHS)
.get_excluded_pool <- function(fit_obj) {
  rhs_all  <- unique(fit_obj$rhs_terms)
  base_exog <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  inst_all <- unique(fit_obj$inst_terms)
  out <- setdiff(inst_all, base_exog)
  out <- out[!grepl("^(1|-1|0|\\+|~|\\|)$", out)]
  unique(out)
}

# df_J mínimo dado um conjunto de instrumentos (por ranks)
.min_dfJ_given_Z <- function(fit_obj, data, inst_terms){
  Xr <- vapply(fit_obj$eq, function(m) qr(.mm(delete.response(terms(formula(m))), data, TRUE))$rank, 1L)
  Zr <- vapply(fit_obj$eq, function(m) qr(.mm(reformulate(inst_terms, intercept = TRUE), data, TRUE))$rank, 1L)
  min(Zr - Xr, na.rm = TRUE)
}

# J por equação para um conjunto arbitrário de IVs excluídos
.j_by_eq_with_excl <- function(fit_obj, data, inst_excl, cluster_var = NULL, verbose = FALSE) {
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_terms <- .build_inst_from_fit(fit_obj, keep_excl = inst_excl)
  rows <- lapply(seq_along(fit_obj$eq), function(i){
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm <- try(AER::ivreg(f_yx, instruments = reformulate(inst_terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) {
      if (isTRUE(verbose)) message("ivreg falhou em ", eqn, ": ", attr(ivm, "condition")$message)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats)) {
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    s_h <- try(AER::sargan(ivm), silent = TRUE)
    J_s <- if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$statistic)) else NA_real_
    p_s <- if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$p.value))   else NA_real_
    Jr  <- .j_overid_manual(ivm, cluster = cl)
    data.frame(eq = eqn, rank_X = rX, rank_Z = rZ, df_J = dfJ,
               J_sargan = J_s, p_sargan = p_s,
               J_hansen = unname(Jr["J"]), p_hansen = unname(Jr["p"]),
               stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(rows)
}

## =======================================================
##  DIAGNÓSTICOS PADRÃO
## =======================================================

jtests_by_eq_manual <- function(fit_obj, data, cluster_var = NULL, verbose = TRUE) {
  inst_excl <- .get_excluded_pool(fit_obj)
  .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var = cluster_var, verbose = verbose)
}

# F em bloco (1º estágio) robusto/cluster e tolerante a singularidade
fs_block_F <- function(fit_obj, data, cluster_var = NULL) {
  data <- as.data.frame(data)
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
  safe_F_one <- function(y) {
    f <- reformulate(c(controls, iv_excl), response = y)
    m <- stats::lm(f, data = data)
    X <- stats::model.matrix(m); u <- stats::residuals(m); n <- NROW(X)
    XtX <- crossprod(X); mu <- mean(diag(XtX)); if (!is.finite(mu) || mu <= 0) mu <- 1
    XtX <- XtX + diag(1e-10 * mu, ncol(XtX)); XtXi <- MASS::ginv(XtX)
    Xu <- X * as.numeric(u)
    meat <- if (!is.null(cluster_var) && cluster_var %in% names(data)) {
      g <- rowsum(Xu, group = as.factor(data[[cluster_var]])); crossprod(as.matrix(g))
    } else crossprod(Xu)
    V <- XtXi %*% meat %*% XtXi
    b  <- coef(m); cn <- names(b); iv_in <- intersect(iv_excl, cn)
    if (!length(iv_in)) return(data.frame(endog=y, F=NA_real_, df1=NA_integer_, df2=NA_integer_, p=NA_real_))
    R <- matrix(0, nrow=length(iv_in), ncol=length(cn), dimnames=list(iv_in, cn)); R[cbind(iv_in,iv_in)] <- 1
    rb <- as.matrix(R %*% b); RVRT <- R %*% V %*% t(R); RVRT <- (RVRT + t(RVRT))/2
    mu2 <- mean(diag(RVRT)); if (!is.finite(mu2) || mu2 <= 0) mu2 <- 1
    RVRT <- RVRT + diag(1e-10 * mu2, nrow(RVRT)); RVRTi <- MASS::ginv(RVRT)
    Wald <- drop(t(rb) %*% RVRTi %*% rb); r <- length(iv_in); Fst <- Wald / r
    rX <- qr(X)$rank; df2 <- n - rX
    if (!is.finite(df2) || df2 <= 0) {
      p <- stats::pchisq(Wald, df = r, lower.tail = FALSE)
      data.frame(endog=y, F=Fst, df1=r, df2=NA_integer_, p=p)
    } else {
      p <- stats::pf(Fst, r, df2, lower.tail = FALSE)
      data.frame(endog=y, F=Fst, df1=r, df2=df2, p=p)
    }
  }
  do.call(rbind, lapply(endo_ln, safe_F_one))
}

# partial R^2 dos excluídos
partial_R2_exclIV <- function(fit_obj, data){
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
  res <- lapply(endo_ln, function(y){
    mfY <- lm(reformulate(controls, response=y), data = data)
    y_t <- residuals(mfY)
    Zt  <- model.matrix(reformulate(iv_excl, intercept = FALSE), data)
    if (length(controls)) {
      Ct  <- model.matrix(reformulate(controls, intercept = TRUE), data)
      Zt  <- resid(lm(Zt ~ Ct))
    }
    fit <- lm(y_t ~ Zt)
    data.frame(endog=y, partial_R2=as.numeric(summary(fit)$r.squared))
  })
  do.call(rbind, res)
}

# Wu–Hausman (control function) robusto/cluster
wu_hausman_cf <- function(fit_obj, data, cluster_var = NULL){
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (!length(endo_ln)) stop("Nenhuma variável endógena detectada (ln_...).")
  Vhat <- lapply(endo_ln, function(y){
    f1 <- reformulate(c(controls, iv_excl), response = y)
    resid(lm(f1, data = data))
  })
  names(Vhat) <- paste0(endo_ln, "_res")
  df2 <- cbind(data, as.data.frame(Vhat))
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  out <- lapply(seq_along(eqs), function(i){
    f <- eqs[[i]]; y <- as.character(f[[2]]); xrhs <- attr(terms(f), "term.labels")
    f2 <- reformulate(c(xrhs, names(Vhat)), response = y)
    ols <- lm(f2, data = df2)
    cl  <- if (!is.null(cluster_var) && cluster_var %in% names(df2)) df2[[cluster_var]] else NULL
    V   <- if (!is.null(cl)) sandwich::vcovCL(ols, cluster = cl, type = "HC3") else sandwich::vcovHC(ols, type = "HC3")
    coef_nms <- names(coef(ols)); test_nms <- intersect(names(Vhat), coef_nms)
    if (!length(test_nms)) return(data.frame(eq=y, k_endog=length(endo_ln), stat=NA_real_, df=NA_integer_, p_wu=NA_real_))
    K <- paste0(test_nms, " = 0")
    lh <- try(car::linearHypothesis(ols, hypothesis.matrix = K, vcov = V, test = "Chisq"), silent = TRUE)
    if (inherits(lh, "try-error")) return(data.frame(eq=y, k_endog=length(endo_ln), stat=NA_real_, df=NA_integer_, p_wu=NA_real_))
    data.frame(eq=y, k_endog=length(endo_ln), stat=as.numeric(lh$Chisq[2]), df=as.integer(lh$Df[2]), p_wu=as.numeric(lh$`Pr(>Chisq)`[2]))
  })
  dplyr::bind_rows(out)
}

## =======================================================
##  C-TESTS por BLOCOS
## =======================================================

make_blocks_explicit <- function(fit_obj) {
  pool <- .get_excluded_pool(fit_obj)
  out <- list()
  out$IV_d_uf_X <- grep("^IV_d_uf_X", pool, value = TRUE)
  out$iv_op     <- grep("^iv_op\\d+",  pool, value = TRUE)
  out <- out[vapply(out, length, integer(1)) > 0]
  out
}

c_tests_blocks_adaptive <- function(fit_obj, data, blocks, cluster_var = NULL, max_k_drop = 2, verbose = FALSE){
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_full <- .get_excluded_pool(fit_obj)
  out <- list()
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]]); eqn <- as.character(f_yx[[2]])
    iv_full <- try(ivreg(f_yx, instruments = reformulate(.build_inst_from_fit(fit_obj, inst_full)), data = data), silent = TRUE)
    if (inherits(iv_full, "try-error")) next
    Jf <- .j_overid_manual(iv_full, cluster = cl)
    if (!is.finite(Jf["df"]) || Jf["df"] <= 0) { if (verbose) message("Eq ", eqn, ": just-ID"); next }
    for (bk in names(blocks)) {
      B <- intersect(blocks[[bk]], inst_full); if (!length(B)) next
      kmax <- min(max_k_drop, length(B))
      for (k in 1:kmax) {
        combs <- utils::combn(B, k, simplify = FALSE)
        for (drop_set in combs) {
          keep <- setdiff(inst_full, drop_set)
          iv_r <- try(ivreg(f_yx, instruments = reformulate(.build_inst_from_fit(fit_obj, keep)), data = data), silent = TRUE)
          if (inherits(iv_r, "try-error")) next
          Jr <- .j_overid_manual(iv_r, cluster = cl)
          C  <- as.numeric(Jr["J"] - Jf["J"])
          pC <- stats::pchisq(C, df = k, lower.tail = FALSE)
          out[[length(out)+1L]] <- data.frame(
            eq = eqn, block = bk, drop_set = paste(drop_set, collapse = " + "),
            k_removed = k, J_full = unname(Jf["J"]), dfJ_full = unname(Jf["df"]),
            J_restr = unname(Jr["J"]), dfJ_restr = unname(Jr["df"]),
            C_stat = C, df_C = k, p_C = pC, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(out)) return(tibble::tibble(eq=character(), block=character(), drop_set=character(), k_removed=integer(),
                                          J_full=double(), dfJ_full=integer(), J_restr=double(), dfJ_restr=integer(),
                                          C_stat=double(), df_C=integer(), p_C=double()))
  dplyr::bind_rows(out)
}

suggest_prune <- function(ctab, alpha = 0.05) {
  if (!nrow(ctab)) return(character(0))
  drops <- unique(unlist(strsplit(ctab$drop_set[which(ctab$p_C < alpha)], "\\s*\\+\\s*")))
  drops[ nzchar(drops) ]
}

## =======================================================
##  REFIT & RELATÓRIO
## =======================================================

refit_quaids_try <- function(fit_obj, data, inst_final_excl) {
  data <- as.data.frame(data)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  inst_terms_new <- .build_inst_from_fit(fit_obj, keep_excl = inst_final_excl)
  fit_new <- try(systemfit(eqs, data = data, method = "3SLS",
                           inst = stats::reformulate(inst_terms_new, intercept = FALSE)),
                 silent = TRUE)
  if (inherits(fit_new, "try-error")) return(NULL)
  shareNames <- fit_obj$shareNames
  keep_sh <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  yhat_mat <- sapply(fit_new$eq, stats::fitted); if (is.list(yhat_mat)) yhat_mat <- do.call(cbind, yhat_mat)
  fitted_included_df <- as.data.frame(yhat_mat, optional = TRUE); colnames(fitted_included_df) <- keep_sh
  fitted_full <- fitted_included_df
  omit_share <- setdiff(shareNames, keep_sh)
  if (length(omit_share)==1) fitted_full[[omit_share]] <- 1 - rowSums(fitted_included_df)
  fitted_full <- fitted_full[, shareNames, drop = FALSE]
  observed_full <- as.data.frame(data[, shareNames, drop = FALSE])
  res <- list(
    method      = "3SLS",
    priceIndex  = fit_obj$priceIndex,
    omit_share  = fit_obj$omit_share,
    drop_price  = fit_obj$drop_price,
    priceNames  = fit_obj$priceNames,
    shareNames  = fit_obj$shareNames,
    rhs_terms   = fit_obj$rhs_terms,
    inst_terms  = inst_terms_new,
    fit         = fit_new,
    eq          = fit_new$eq,
    fitted_shares   = fitted_full,
    observed_shares = observed_full
  )
  class(res) <- "quaids_km1_fit"
  res
}

post_refit_report <- function(fit_new, data, cluster_var = NULL, R = 2000, level = 0.95, seed = 123){
  cat("\n== Over-ID (Hansen J) por equação ==\n")
  jt <- jtests_by_eq_manual(fit_new, data, cluster_var = cluster_var, verbose = FALSE)
  print(jt)
  if (any(jt$df_J <= 0, na.rm = TRUE)) {
    message("Aviso: há equações just/under-ID (df_J <= 0). Hansen J não é aplicável nessas eqs.")
  }
  cat("\n== Força dos IVs (F-HC3/CL) ==\n")
  print(fs_block_F(fit_new, data, cluster_var = cluster_var))
  cat("\n== partial R^2 dos IVs excluídos ==\n")
  print(partial_R2_exclIV(fit_new, data))
  cat("\n== Wu-Hausman (control function) por equação ==\n")
  print(wu_hausman_cf(fit_new, data, cluster_var = cluster_var))
  if (exists("quaids_elasticities_ci", mode = "function")) {
    cat("\n== Elasticidades com IC (ponto: median/observed) ==\n")
    res_ci <- try(quaids_elasticities_ci(fit_new, at="medians", w_source="observed", R=R, level=level, seed=seed), silent = TRUE)
    if (!inherits(res_ci, "try-error") && is.list(res_ci) && !is.null(res_ci$expenditure)) {
      print(res_ci$expenditure); invisible(res_ci)
    } else {
      message("(Falha ao computar IC; reportando ponto).")
      if (exists("quaids_elasticities", mode = "function")) {
        res_pt <- try(quaids_elasticities(fit_new, at="medians", w_source="observed"), silent = TRUE)
        if (!inherits(res_pt, "try-error") && is.list(res_pt) && !is.null(res_pt$expenditure)) print(res_pt$expenditure)
      }
      invisible(NULL)
    }
  } else {
    message("\n(Elasticidades/IC puladas: 'quaids_elasticities_ci' não encontrada.)")
    invisible(NULL)
  }
}

# Moment diagnostics (com opção de ver apenas excluídos, sem under-ID)
iv_moment_diag <- function(fit_obj_or_fitnew, data, cluster_var = NULL, only_excluded = FALSE){
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_obj_or_fitnew$eq, function(m) formula(m))
  inst_terms <- fit_obj_or_fitnew$inst_terms
  rhs_all    <- fit_obj_or_fitnew$rhs_terms
  base_exog  <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  excl       <- setdiff(inst_terms, base_exog)
  out <- lapply(eqs, function(f){
    eqn <- as.character(f[[2]])
    ivm <- AER::ivreg(f, instruments = reformulate(inst_terms), data = data)
    mats <- .get_XZ_u(ivm); if (is.null(mats)) return(NULL)
    Z <- mats$Z; u <- as.numeric(mats$u); n <- NROW(Z)
    Zu <- Z * u
    S  <- if (!is.null(cl)) { G <- rowsum(Zu, group = as.factor(.align_cluster(ivm, cl))); crossprod(as.matrix(G))/n
    } else crossprod(Zu)/n
    gbar <- colMeans(Zu); se <- sqrt(diag(S)/n); t <- gbar / se
    tab  <- data.frame(eq = eqn, inst = colnames(Z), gbar = gbar, se = se, t = t, abs_t = abs(t))
    if (isTRUE(only_excluded)) tab <- dplyr::filter(tab, inst %in% excl)
    tab
  })
  dplyr::bind_rows(out) |> dplyr::filter(inst != "(Intercept)")
}

## =======================================================
##  NOVO: PRUNE FOCADO NA PIOR EQUAÇÃO DO J (com guarda de df_J)
## =======================================================

prune_focus_on_eq <- function(fit_base, data, fit_current, focus_eq,
                              cluster_var = NULL, alpha = 0.05,
                              max_k_drop = 2, target_min_dfJ = 1,
                              max_steps = 3, verbose = TRUE){
  current_excl <- .get_excluded_pool(fit_current)
  for (step in seq_len(max_steps)) {
    if (isTRUE(verbose)) message("\n[focus] passo ", step, " | eq alvo: ", focus_eq)
    # C-tests por blocos, focando a eq alvo
    blocks <- make_blocks_explicit(fit_current)
    ctab <- c_tests_blocks_adaptive(fit_current, data, blocks,
                                    cluster_var = cluster_var, max_k_drop = max_k_drop, verbose = FALSE)
    ctab_eq <- ctab |> dplyr::filter(eq == !!focus_eq) |> dplyr::arrange(p_C)
    if (!nrow(ctab_eq) || all(!is.finite(ctab_eq$p_C))) { if (verbose) message("  (sem C-tests viáveis)"); break }
    to_drop <- suggest_prune(ctab_eq, alpha = alpha)
    if (!length(to_drop)) { if (verbose) message("  (nenhum drop com p_C<", alpha, ")"); break }
    if (isTRUE(verbose)) message("  candidatos: {", paste(to_drop, collapse=", "), "}")
    # tenta remover em lote
    try_batch <- setdiff(current_excl, to_drop)
    inst_terms_try <- .build_inst_from_fit(fit_base, keep_excl = try_batch)
    mindf <- .min_dfJ_given_Z(fit_base, data, inst_terms_try)
    if (mindf >= target_min_dfJ) {
      if (isTRUE(verbose)) message("  aceita lote (min df_J=", mindf, ")")
      current_excl <- try_batch
    } else {
      if (isTRUE(verbose)) message("  lote derrubaria df_J < ", target_min_dfJ, " — testando remoções unitárias")
      # ordena pelos menores p_C (mais "suspeitos")
      ord <- unique(unlist(strsplit(ctab_eq$drop_set[ctab_eq$p_C < alpha], "\\s*\\+\\s*")))
      ord <- ord[nzchar(ord)]
      accepted_any <- FALSE
      for (z in ord) {
        try_one <- setdiff(current_excl, z)
        inst_terms_try2 <- .build_inst_from_fit(fit_base, keep_excl = try_one)
        mindf2 <- .min_dfJ_given_Z(fit_base, data, inst_terms_try2)
        if (mindf2 >= target_min_dfJ) {
          if (isTRUE(verbose)) message("   - removendo {", z, "} (min df_J=", mindf2, ")")
          current_excl <- try_one; accepted_any <- TRUE; break
        } else if (isTRUE(verbose)) {
          message("   - NÃO removi {", z, "} (ficaria min df_J=", mindf2, ")")
        }
      }
      if (!accepted_any) { if (verbose) message("  (nenhuma remoção unitária aceita)"); break }
    }
    # refit e atualiza fit_current
    fit_current <- refit_quaids_try(fit_base, data, inst_final_excl = current_excl)
    if (is.null(fit_current)) { warning("Refit falhou após prune focado."); break }
    # se a eq alvo já ficou ok (p>alpha), pode parar cedo
    jt_now <- jtests_by_eq_manual(fit_current, data, cluster_var = cluster_var, verbose = FALSE)
    p_focus <- jt_now$p_hansen[jt_now$eq == focus_eq]
    if (isTRUE(verbose)) message("  p_hansen(", focus_eq, ") = ", signif(p_focus, 3))
    if (is.finite(p_focus) && p_focus > alpha) break
  }
  list(new_fit = fit_current, inst_excl_final = current_excl)
}

## =======================================================
##  DRIVER — PRÓXIMOS PASSOS “TOP TIER”
##  Pré-requisito: objetos `fit_q` (quaids_km1_fit) e `df_iv` no ambiente
## =======================================================

stopifnot(exists("fit_q"), exists("df_iv"))
rhs_all <- fit_q$rhs_terms
cat("\nEndógenos no RHS (esperado fora de Z): ",
    paste(intersect(rhs_all, grep('^ln_', rhs_all, value=TRUE)), collapse=', '), "\n", sep="")

# escolhe fit de partida (priorize o mais recente se existir)
fit_start <- if (exists("fit_q_dfJ2")) fit_q_dfJ2 else if (exists("fit_q_new")) fit_q_new else fit_q
cl_var <- pick_cluster(df_iv)

cat("\n== Diagnóstico base (fit de partida) ==\n")
jt0 <- jtests_by_eq_manual(fit_start, df_iv, cluster_var = cl_var, verbose = FALSE)
print(jt0)
eq_worst <- jt0$eq[ which.min(jt0$p_hansen) ]
cat("\nEquação com menor p(Hansen J): ", eq_worst, "\n", sep="")

# prune focado na pior eq (guarda df_J >= 1), no máx. 3 passos
focus_out <- prune_focus_on_eq(
  fit_base      = fit_q,
  data          = df_iv,
  fit_current   = fit_start,
  focus_eq      = eq_worst,
  cluster_var   = cl_var,
  alpha         = 0.05,
  max_k_drop    = 2,
  target_min_dfJ= 1,
  max_steps     = 3,
  verbose       = TRUE
)

fit_top <- focus_out$new_fit

cat("\n== Relatório pós-prune focado (TOP TIER snapshot) ==\n")
invisible(post_refit_report(fit_top, df_iv, cluster_var = cl_var, R = 2000, level = 0.95, seed = 123))

cat("\n== Moment diagnostics (apenas IVs excluídos) no set final ==\n")
diag_tab <- iv_moment_diag(fit_top, df_iv, cluster_var = cl_var, only_excluded = TRUE)
print(dplyr::arrange(diag_tab, dplyr::desc(abs_t)) %>%
        dplyr::group_by(eq) %>% dplyr::slice_head(n = 8), n = 36)

cat("\n[INFO] Conjunto final de IVs EXCLUÍDOS:\n  ",
    paste(.get_excluded_pool(fit_top), collapse = " + "), "\n", sep="")
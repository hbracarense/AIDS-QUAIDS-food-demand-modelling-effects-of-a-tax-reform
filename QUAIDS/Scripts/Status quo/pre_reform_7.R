## =======================================================
## TOP-TIER (plug-and-play): poda adaptativa + refit QUAIDS + checklist
## =======================================================

suppressPackageStartupMessages({
  library(AER)
  library(systemfit)
  library(sandwich)
  library(dplyr)
  library(tidyr)
  library(car)
  library(openxlsx)
})

## =======================================================
## SUA FUNÇÃO DE ESTIMAÇÃO (fornecida)  -------------------
## =======================================================
fit_quaids_manual_km1 <- function(
    prices,                      # data.frame/matrix com preços (colunas = priceNames)
    shares,                      # data.frame/matrix com shares (colunas = shareNames)
    x,                           # gasto total (mantido por compatibilidade)
    priceIndex = "Ls",
    estMethod  = "3SLS",
    omit_share = 1,              # índice do share omitido (w1 por padrão)
    drop_price = 1,              # índice do preço/base omitido no RHS (p1 por padrão)
    instNames,                   # vetor com IVs "core"
    instData,                    # data.frame com todas as variáveis (df_iv)
    use_z2     = FALSE,
    maxiter    = 500,
    # opções:
    include_cap     = TRUE,
    demo_pref       = c("p_n_0esc","p_n_8esc"),
    center_sex      = TRUE,
    add_interaction = TRUE,
    compute_diag    = TRUE       # apenas kappa / R2_corr / adding-up
) {
  
  if (!requireNamespace("systemfit", quietly = TRUE)) {
    stop("Pacote 'systemfit' é necessário.")
  }
  
  # ---- Nomes e checagens básicas ----
  priceNames <- colnames(as.data.frame(prices))
  shareNames <- colnames(as.data.frame(shares))
  if (is.null(priceNames) || is.null(shareNames)) {
    stop("As matrizes 'prices' e 'shares' precisam ter nomes de colunas.")
  }
  if (nrow(prices) != nrow(shares) || nrow(prices) != nrow(instData)) {
    stop("As quantidades de linhas de 'prices', 'shares' e 'instData' devem coincidir.")
  }
  if (!all(instNames %in% names(instData))) {
    faltam <- setdiff(instNames, names(instData))
    stop("IVs ausentes em instData: ", paste(faltam, collapse = ", "))
  }
  
  df <- instData
  
  # ---- (A) f_cap (se faltar) ----
  if (include_cap && !"f_cap" %in% names(df)) {
    if (!"capital_interior" %in% names(df)) {
      stop("Nem 'f_cap' nem 'capital_interior' encontrados em instData.")
    }
    df$f_cap <- factor(df$capital_interior)
  }
  
  # Garante fatores existentes
  if (!"f_reg"  %in% names(df)) stop("'f_reg' ausente em instData.")
  if (!"f_area" %in% names(df)) stop("'f_area' ausente em instData.")
  if (!is.factor(df$f_reg))  df$f_reg  <- factor(df$f_reg)
  if (!is.factor(df$f_area)) df$f_area <- factor(df$f_area)
  if (include_cap && !is.factor(df$f_cap)) df$f_cap <- factor(df$f_cap)
  
  # ---- (B) demográfico extra ----
  demo_pref  <- match.arg(demo_pref)
  demo_extra <- if (demo_pref %in% names(df)) {
    demo_pref
  } else if (identical(demo_pref, "p_n_0esc") && "p_n_8esc" %in% names(df)) {
    "p_n_8esc"
  } else if (identical(demo_pref, "p_n_8esc") && "p_n_0esc" %in% names(df)) {
    "p_n_0esc"
  } else {
    stop("Nem 'p_n_0esc' nem 'p_n_8esc' encontrados em instData.")
  }
  
  # ---- (C) centraliza p_sexofem e cria interação única com f_reg ----
  if (!"p_sexofem" %in% names(df)) stop("'p_sexofem' ausente em instData.")
  sex_var <- "p_sexofem"
  if (isTRUE(center_sex)) {
    df$p_sexofem_c <- as.numeric(scale(df$p_sexofem, center = TRUE, scale = FALSE))
    sex_var <- "p_sexofem_c"
  }
  inter_term <- if (isTRUE(add_interaction)) paste0(sex_var, ":f_reg") else NULL
  
  # ---- ln-preços no RHS (K-1) ----
  ln_keep <- setdiff(paste0("ln_", priceNames), paste0("ln_", priceNames[drop_price]))
  if (!all(ln_keep %in% names(df))) {
    faltam <- setdiff(ln_keep, names(df))
    stop("Variáveis de ln-preço ausentes em instData: ", paste(faltam, collapse = ", "))
  }
  
  # ---- RHS e instrumentos (limpando NAs/vazios) ----
  base_exog <- c("z", "f_reg", "f_area", if (include_cap) "f_cap",
                 "p_n_mais65anos", sex_var, demo_extra)
  rhs_terms  <- c(ln_keep, base_exog, inter_term)
  inst_terms <- c("z", "f_reg", "f_area", if (include_cap) "f_cap",
                  "p_n_mais65anos", sex_var, demo_extra, inter_term, instNames)
  if (isTRUE(use_z2) && "z2" %in% names(df)) {
    rhs_terms  <- c(rhs_terms,  "z2")
    inst_terms <- c(inst_terms, "z2")
  }
  rhs_terms  <- unique(rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)])
  inst_terms <- unique(inst_terms[!is.na(inst_terms) & nzchar(inst_terms)])
  
  # ---- Equações K-1 (omite share omitido) ----
  keep_sh <- setdiff(shareNames, shareNames[omit_share])
  if (!all(keep_sh %in% names(df))) {
    faltam <- setdiff(keep_sh, names(df))
    stop("Shares ausentes em instData: ", paste(faltam, collapse = ", "))
  }
  eqs <- lapply(keep_sh, function(y) stats::reformulate(rhs_terms, response = y))
  names(eqs) <- gsub("[^[:alnum:]]", "", keep_sh)  # systemfit não aceita '_' nos rótulos
  
  # ---- Estimação ----
  if (!identical(estMethod, "3SLS")) {
    warning("Estimation method alterado para '3SLS' (único implementado aqui).")
  }
  fit <- systemfit::systemfit(
    eqs, data = df,
    method = "3SLS",
    inst   = stats::reformulate(inst_terms, intercept = FALSE),
    maxit  = maxiter
  )
  
  # ---- Recria share omitido (observado e previsto) ----
  shares_df <- as.data.frame(shares)
  # Observado via adding-up das K-1 observadas
  obs_included <- as.matrix(shares_df[, keep_sh, drop = FALSE])
  w_omit_obs_from_add <- 1 - rowSums(obs_included)
  
  # Previsto via soma dos fitted das K-1
  yhat_mat <- sapply(fit$eq, stats::fitted)
  if (is.list(yhat_mat)) yhat_mat <- do.call(cbind, yhat_mat)
  w_omit_hat <- 1 - rowSums(yhat_mat)
  
  # ---- Diagnósticos leves (opcional) ----
  diag <- NULL
  if (isTRUE(compute_diag)) {
    kappa_eq <- sapply(fit$eq, function(m){
      X <- stats::model.matrix(m)
      X <- scale(X[, colnames(X) != "(Intercept)"], TRUE, TRUE)
      kappa(X)
    })
    r2_corr <- sapply(fit$eq, function(m){
      y <- stats::model.response(stats::model.frame(m))
      stats::cor(y, stats::fitted(m))^2
    })
    add_up_max_abs <- max(abs(w_omit_hat + rowSums(yhat_mat) - 1))
    names(kappa_eq) <- names(fit$eq)
    names(r2_corr)  <- names(fit$eq)
    diag <- list(
      kappa_eq       = kappa_eq,
      r2_corr        = r2_corr,
      add_up_max_abs = add_up_max_abs
    )
  }
  
  # ---- Monta fitted/observed completos (incluindo a omitida) ----
  fitted_included_df <- as.data.frame(yhat_mat, optional = TRUE)
  colnames(fitted_included_df) <- keep_sh
  fitted_full <- fitted_included_df
  fitted_full[[shareNames[omit_share]]] <- w_omit_hat
  fitted_full <- fitted_full[, shareNames, drop = FALSE]
  
  observed_included_df <- as.data.frame(obs_included, optional = TRUE)
  colnames(observed_included_df) <- keep_sh
  observed_full <- observed_included_df
  observed_full[[shareNames[omit_share]]] <- w_omit_obs_from_add
  observed_full <- observed_full[, shareNames, drop = FALSE]
  
  # ---- Retorno ----
  res <- list(
    call        = match.call(),
    method      = "3SLS",
    priceIndex  = priceIndex,
    omit_share  = shareNames[omit_share],
    drop_price  = priceNames[drop_price],
    priceNames  = priceNames,
    shareNames  = shareNames,
    rhs_terms   = rhs_terms,
    inst_terms  = inst_terms,
    fit         = fit,
    eq          = fit$eq,
    fitted_shares   = fitted_full,
    observed_shares = observed_full,
    omitted_recreated = list(
      name     = shareNames[omit_share],
      observed = w_omit_obs_from_add,
      fitted   = w_omit_hat
    ),
    diagnostics = diag
  )
  class(res) <- "quaids_km1_fit"
  return(res)
}

summary.quaids_km1_fit <- function(object, ...) {
  cat("Demand analysis with the QUADRATIC Almost Ideal Demand System (QUAIDS)\n")
  cat("Estimation Method:", object$method, "\n")
  cat("Price Index:", object$priceIndex, "\n")
  cat("Omitted share equation:", object$omit_share, "\n")
  cat("Dropped price from RHS:", object$drop_price, "\n\n")
  print(summary(object$fit))
  if (!is.null(object$diagnostics)) {
    cat("\nQuick diagnostics (no tests)\n")
    with(object$diagnostics, {
      cat("  kappa (std. X) by equation:\n"); print(kappa_eq)
      cat("  R2 (corr^2) by equation:\n");   print(r2_corr)
      cat("  Adding-up |max abs deviation| (fitted):",
          format(add_up_max_abs, digits = 6), "\n")
    })
  }
  invisible(object)
}

## =======================================================
## HELPERs + DIAGNÓSTICOS (J/C-tests, força IV, poda) ----
## =======================================================

.align_cluster <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)
  mf  <- stats::model.frame(model)
  rid <- rownames(mf)
  if (!is.null(names(cluster))) return(cluster[rid])
  idx <- suppressWarnings(as.integer(rid))
  if (all(is.finite(idx))) return(cluster[idx])
  cluster
}

.get_XZ_u <- function(ivm) {
  X <- try(model.matrix(ivm, "regressors"),  silent = TRUE)
  Z <- try(model.matrix(ivm, "instruments"), silent = TRUE)
  if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NULL)
  list(X = X, Z = Z, u = residuals(ivm))
}

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
  S <- (S + t(S))/2
  d <- mean(diag(S)); if (!is.finite(d) || d <= 0) d <- 1
  S_r <- S + diag(ridge * d, ncol(S))
  gbar <- colMeans(Zu)
  J    <- as.numeric(n * crossprod(gbar, solve(S_r, gbar)))
  p    <- stats::pchisq(J, df = df, lower.tail = FALSE)
  c(J = J, df = df, p = p)
}

.clean_terms <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(x)]
  setdiff(x, c("1","-1","0","+","~","|"))
}

.build_inst_terms <- function(fit_obj) {
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  unique(c(exogs, inst))
}

.get_excluded_IVs <- function(fit_obj) {
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  setdiff(inst, exogs)
}

.fit_iv_safe <- function(f_yx, inst_terms, data) {
  try(AER::ivreg(f_yx, instruments = reformulate(inst_terms), data = data), silent = TRUE)
}

jtests_by_eq_manual <- function(fit_obj, data, cluster_var = NULL, verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_terms <- .build_inst_terms(fit_obj)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  rows <- lapply(seq_along(fit_obj$eq), function(i){
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm  <- .fit_iv_safe(f_yx, inst_terms, data)
    if (inherits(ivm, "try-error")) {
      if (verbose) message("ivreg falhou na eq ", eqn, ": ", attr(ivm, "condition")$message)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats)) {
      if (verbose) message("Matrizes X/Z indisponíveis para ", eqn)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    J_h <- try(AER::sargan(ivm), silent = TRUE)
    J_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$statistic)) else NA_real_
    p_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$p.value))   else NA_real_
    Jrob <- .j_overid_manual(ivm, cluster = cl)
    data.frame(
      eq = eqn, rank_X = rX, rank_Z = rZ, df_J = dfJ,
      J_sargan = J_s, p_sargan = p_s,
      J_hansen = unname(Jrob["J"]), p_hansen = unname(Jrob["p"]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

make_blocks_explicit <- function(fit_obj) {
  excl <- .get_excluded_IVs(fit_obj)
  if (!length(excl)) return(list())
  base <- gsub("\\d+$", "", excl)
  base <- gsub("\\s+", "", base)
  sp   <- split(excl, base)
  names(sp) <- make.unique(names(sp))
  sp
}

make_blocks_from_vector <- function(excl_vec) {
  if (!length(excl_vec)) return(list())
  base <- gsub("\\d+$", "", excl_vec)
  base <- gsub("\\s+", "", base)
  sp   <- split(excl_vec, base)
  names(sp) <- make.unique(names(sp))
  sp
}

c_tests_blocks_adaptive <- function(fit_obj, data, blocks,
                                    cluster_var = NULL,
                                    max_k_drop = 2,
                                    verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_base <- .build_inst_terms(fit_obj)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  out <- list()
  
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    
    iv_full <- .fit_iv_safe(f_yx, inst_base, data)
    if (inherits(iv_full, "try-error")) { if (verbose) message("ivreg (full) falhou em ", eqn); next }
    Jf <- .j_overid_manual(iv_full, cluster = cl)
    fm <- .get_XZ_u(iv_full); if (is.null(fm)) next
    df_full <- ncol(fm$Z) - ncol(fm$X)
    if (!is.finite(df_full) || df_full <= 0) { if (verbose) message("Eq ", eqn, ": full não sobre-ID."); next }
    k_cap <- max(1, min(max_k_drop, df_full - 1))
    
    for (bk in names(blocks)) {
      in_block <- intersect(inst_base, blocks[[bk]])
      if (!length(in_block)) next
      for (k in 1:k_cap) {
        if (length(in_block) < k) next
        cmb <- utils::combn(in_block, k, simplify = FALSE)
        for (drop_set in cmb) {
          inst_restr <- setdiff(inst_base, drop_set)
          iv_r <- .fit_iv_safe(f_yx, inst_restr, data)
          if (inherits(iv_r, "try-error")) next
          rm <- .get_XZ_u(iv_r); if (is.null(rm)) next
          df_restr <- ncol(rm$Z) - ncol(rm$X)
          if (!is.finite(df_restr) || df_restr <= 0) next
          Jr <- .j_overid_manual(iv_r, cluster = cl)
          C  <- as.numeric(Jr["J"] - Jf["J"])
          dfC <- length(drop_set)
          pC <- stats::pchisq(C, df = dfC, lower.tail = FALSE)
          out[[length(out) + 1L]] <- data.frame(
            eq = eqn, block = bk, drop_set = paste(drop_set, collapse=" + "),
            k_removed = dfC,
            J_full = unname(Jf["J"]), dfJ_full = unname(Jf["df"]),
            J_restr = unname(Jr["J"]), dfJ_restr = unname(Jr["df"]),
            C_stat = C, df_C = dfC, p_C = pC, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(out)) return(data.frame(
    eq=character(), block=character(), drop_set=character(), k_removed=integer(),
    J_full=double(), dfJ_full=integer(), J_restr=double(), dfJ_restr=integer(),
    C_stat=double(), df_C=integer(), p_C=double()
  ))
  dplyr::bind_rows(out)
}

suggest_prune <- function(ctab, alpha = 0.05) {
  if (!NROW(ctab)) return(character(0))
  bad <- subset(ctab, is.finite(p_C) & p_C < alpha)
  unique(bad$drop_set)
}

new_inst_after_prune <- function(fit_obj, to_drop_sets, blocks = NULL) {
  inst_base <- .build_inst_terms(fit_obj)
  rhs_all   <- .clean_terms(fit_obj$rhs_terms)
  endo_ln   <- grep("^ln_", rhs_all, value = TRUE)
  exogs     <- setdiff(rhs_all, endo_ln)
  drop_vec  <- character(0)
  if (length(to_drop_sets)) drop_vec <- unique(unlist(strsplit(to_drop_sets, "\\s\\+\\s")))
  excl_ivs  <- setdiff(inst_base, exogs)
  excl_new  <- setdiff(excl_ivs, drop_vec)
  unique(c(exogs, excl_new))
}

auto_prune_until_ok <- function(fit_obj, data, cluster_var = NULL,
                                alpha = 0.05, max_k_drop = 2, max_rounds = 5,
                                verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_excl <- .get_excluded_IVs(fit_obj)
  exogs     <- setdiff(.build_inst_terms(fit_obj), inst_excl)
  
  history <- list()
  for (r in seq_len(max_rounds)) {
    fit_tmp <- fit_obj
    fit_tmp$inst_terms <- inst_excl   # apenas excluídos; exógenos serão recolocados por .build_inst_terms()
    blocks  <- make_blocks_from_vector(inst_excl)
    if (!length(blocks)) {
      if (verbose) message("Sem IVs excluídos; nada para podar.")
      break
    }
    ctab <- c_tests_blocks_adaptive(
      fit_tmp, data, blocks,
      cluster_var = cluster_var, max_k_drop = max_k_drop, verbose = verbose
    )
    if (!nrow(ctab)) { if (verbose) message("Nenhum C-test viável nesta rodada."); break }
    bad_sets <- suggest_prune(ctab, alpha = alpha)
    history[[r]] <- list(ctab = ctab, bad_sets = bad_sets, inst_excl = inst_excl)
    if (!length(bad_sets)) { if (verbose) message("Nenhum subconjunto com p_C < ", alpha, "."); break }
    worst <- ctab %>% filter(is.finite(p_C)) %>% arrange(p_C) %>% slice(1) %>% pull(drop_set)
    drop_vec <- unique(unlist(strsplit(worst, "\\s\\+\\s")))
    inst_excl <- setdiff(inst_excl, drop_vec)
    if (verbose) message("Rodada ", r, ": removendo {", paste(drop_vec, collapse=", "), "}.")
  }
  
  inst_final <- unique(c(exogs, inst_excl))
  fit_final  <- fit_obj; fit_final$inst_terms <- inst_excl  # para diagnósticos (apenas excluídos)
  list(inst_final = inst_final, inst_excl_final = inst_excl,
       history = history, fit_for_diagnostics = fit_final)
}

fs_block_F2 <- function(fit_obj, data, robust = TRUE){
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  inst_all <- unique(c(.get_excluded_IVs(fit_obj), setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))))
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
  form_fs <- function(y) reformulate(c(controls, iv_excl), response = y)
  get_one <- function(y){
    m  <- stats::lm(form_fs(y), data = data)
    K  <- paste(iv_excl, "= 0")
    lh <- try(
      car::linearHypothesis(m, K,
                            vcov = if (robust) sandwich::vcovHC(m, type="HC3") else NULL,
                            test = "F"), silent = TRUE)
    if (inherits(lh, "try-error")) {
      return(data.frame(endog=y, F=NA_real_, df1=NA_real_, df2=NA_real_, p=NA_real_))
    }
    data.frame(endog=y, F=unname(lh$F[2]), df1=unname(lh$Df[2]),
               df2=unname(lh$Res.Df[2]), p=unname(lh$`Pr(>F)`[2]))
  }
  do.call(rbind, lapply(endo_ln, get_one))
}

partial_R2_exclIV2 <- function(fit_obj, data){
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  inst_all <- unique(c(.get_excluded_IVs(fit_obj), setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))))
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
  res <- lapply(endo_ln, function(y){
    mfY <- lm(reformulate(controls, response=y), data = data)
    y_t <- residuals(mfY)
    Zt  <- model.matrix(reformulate(iv_excl, intercept = FALSE), data = data)
    if (length(controls)) {
      Ct  <- model.matrix(reformulate(controls, intercept = TRUE), data = data)
      Zt  <- resid(lm(Zt ~ Ct))
    }
    fit <- lm(y_t ~ Zt)
    R2  <- summary(fit)$r.squared
    data.frame(endog=y, partial_R2=as.numeric(R2))
  })
  do.call(rbind, res)
}

## =======================================================
## REFIT AUTOMÁTICO USANDO SUA FUNÇÃO  -------------------
## =======================================================

refit_quaids_try <- function(fit_obj, data, inst_final_excl) {
  # prepara entradas
  priceNames <- fit_obj$priceNames
  shareNames <- fit_obj$shareNames
  prices <- as.data.frame(data[, priceNames, drop = FALSE])
  shares <- as.data.frame(data[, shareNames, drop = FALSE])
  
  # x: tenta achar gasto total usado no seu fluxo
  x_vec <- NULL
  for (cand in c("gasto_total_atualhat","gasto_total","x","despesa_total")) {
    if (cand %in% names(data)) { x_vec <- data[[cand]]; break }
  }
  if (is.null(x_vec)) x_vec <- rowSums(prices)  # fallback inócuo
  
  # índices omitidos (se no objeto vierem como nomes)
  omit_share_idx <- match(fit_obj$omit_share, shareNames)
  if (is.na(omit_share_idx)) omit_share_idx <- 1L
  drop_price_idx <- match(fit_obj$drop_price, priceNames)
  if (is.na(drop_price_idx)) drop_price_idx <- 1L
  
  # demográfico preferido
  demo_pref <- if ("p_n_0esc" %in% names(data)) "p_n_0esc" else "p_n_8esc"
  
  # usar z2 se existir
  use_z2 <- "z2" %in% names(data)
  
  fit_new <- fit_quaids_manual_km1(
    prices      = prices,
    shares      = shares,
    x           = x_vec,
    priceIndex  = "Ls",
    estMethod   = "3SLS",
    omit_share  = omit_share_idx,
    drop_price  = drop_price_idx,
    instNames   = inst_final_excl,     # <<< só os IVs EXCLUÍDOS
    instData    = as.data.frame(data),
    use_z2      = use_z2,
    maxiter     = 600,
    include_cap     = "f_cap" %in% names(data) || "capital_interior" %in% names(data),
    demo_pref       = demo_pref,
    center_sex      = TRUE,
    add_interaction = TRUE,
    compute_diag    = TRUE
  )
  fit_new
}

post_refit_report <- function(fit_new, data, cluster_var = NULL, R = 2000, level = 0.95, seed = 123) {
  message("\n== Over-ID (Hansen J) por equação ==")
  print(jtests_by_eq_manual(fit_new, data, cluster_var = cluster_var, verbose = FALSE))
  message("\n== Força dos IVs (F-HC3) ==")
  print(fs_block_F2(fit_new, data, robust = TRUE))
  message("\n== partial R^2 dos IVs excluídos ==")
  print(partial_R2_exclIV2(fit_new, data))
  if (exists("wu_por_eq_manual2")) {
    message("\n== Wu-Hausman (control function) por equação ==")
    print(wu_por_eq_manual2(fit_new, data, robust = TRUE, verbose = FALSE))
  }
  if (exists("quaids_elasticities_ci")) {
    message("\n== Elasticidades com IC (ponto: median/observed) ==")
    res_ci <- quaids_elasticities_ci(fit_new, at="medians", w_source="observed",
                                     R=R, level=level, seed=seed)
    print(res_ci$expenditure)
    return(invisible(res_ci))
  }
  invisible(NULL)
}

## =======================================================
## COMO USAR (roda direto)  ------------------------------
## =======================================================

# Supõe que você já tem 'fit_q' (objeto quaids_km1_fit) e 'df_iv' no ambiente.

# 0) Sanity: endógenos no RHS
rhs_all <- fit_q$rhs_terms
cat("\nEndógenos no RHS (esperado fora de Z): ",
    paste(intersect(rhs_all, grep('^ln_', rhs_all, value=TRUE)), collapse=', '), "\n", sep="")

# 1) Poda adaptativa até estabilizar (ou max_rounds)
prune_out <- auto_prune_until_ok(
  fit_obj     = fit_q,
  data        = df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  alpha       = 0.05,
  max_k_drop  = 2,     # df_J=3 => testa 1 e 2 por vez
  max_rounds  = 5,
  verbose     = TRUE
)

cat("\nInst terms (FINAIS p/ diagnóstico/refit):\n",
    paste(prune_out$inst_final, collapse = " + "), "\n")

# 2) Diagnóstico com inst_finais (ainda no fit antigo)
message("\n== Diagnóstico com inst_finais (pré-refit) ==")
print(jtests_by_eq_manual(prune_out$fit_for_diagnostics, df_iv,
                          cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"))

# 3) Refit do QUAIDS com IVs podados (usa SUA função acima)
inst_excl_final <- prune_out$inst_excl_final
fit_q_new <- refit_quaids_try(fit_q, df_iv, inst_final_excl = inst_excl_final)

# 4) Checklist pós-refit (+ elasticidades/IC se sua função existir)
if (!is.null(fit_q_new)) {
  res_ci_new <- post_refit_report(
    fit_new    = fit_q_new,
    data       = df_iv,
    cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
    R          = 2000, level = 0.95, seed = 123
  )
  # (opcional) daqui você chama seus heatmaps/export se já tiver tidy_elas_safe/plot_heatmap/export_elas
}

## FIM

## =======================================================
## EXTRA TOP-TIER: momentos por IV + poda gulosa (greedy)
## =======================================================

# --- usa helpers já definidos no seu script:
# .build_inst_terms(), .get_excluded_IVs(), .get_XZ_u(), .align_cluster(),
# .fit_iv_safe(), .j_overid_manual(), jtests_by_eq_manual()

# 1) Tabela de momentos por instrumento (t robusto/cluster)
iv_moment_diag <- function(fit_obj, data, cluster_var = NULL) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  
  inst_base <- .build_inst_terms(fit_obj)
  excl_only <- .get_excluded_IVs(fit_obj)         # só IVs excluídos para ranquear
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  rows <- lapply(seq_along(fit_obj$eq), function(i){
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm  <- .fit_iv_safe(f_yx, inst_base, data)
    mm   <- .get_XZ_u(ivm); if (is.null(mm)) return(NULL)
    
    # momentos e cov robusta (Arellano/cluster) como no seu J manual
    u  <- as.numeric(mm$u); Z <- mm$Z
    Zu <- Z * u
    if (!is.null(cl)) {
      G  <- rowsum(Zu, group = as.factor(.align_cluster(ivm, cl)))
      S  <- crossprod(as.matrix(G)) / nrow(Z)      # var de sqrt(n)*gbar
    } else {
      S  <- crossprod(Zu) / nrow(Z)
    }
    gbar <- colMeans(Zu)                           # (1/n) Z'u
    # t para gbar: gbar_j / sqrt(Var(gbar_j)) = gbar_j / sqrt(S_jj / n)
    n <- nrow(Z)
    v_gbar <- diag(S) / n
    v_gbar[!is.finite(v_gbar) | v_gbar <= 0] <- NA_real_
    t_j <- gbar / sqrt(v_gbar)
    
    tibble::tibble(
      eq = eqn,
      inst = colnames(Z),
      t_moment = as.numeric(t_j),
      abs_t = abs(t_moment),
      is_excluded = inst %in% excl_only
    )
  })
  out <- dplyr::bind_rows(rows)
  # mantém só IVs excluídos (os que podemos podar)
  dplyr::filter(out, is_excluded) %>%
    dplyr::arrange(eq, dplyr::desc(abs_t))
}

# 2) Pequeno utilitário: J por equação com um conjunto específico de IVs excluídos
.j_by_eq_with_excl <- function(fit_obj, data, inst_excl, cluster_var = NULL) {
  fit_tmp <- fit_obj
  fit_tmp$inst_terms <- inst_excl        # só excluídos; exógenos entram por .build_inst_terms()
  jtests_by_eq_manual(fit_tmp, data, cluster_var = cluster_var, verbose = FALSE)
}

# 3) Poda gulosa orientada por momentos (até p_Hansen >= alpha ou df_J == 0)
# === Greedy SAFE: preserva sobre-ID (df_J >= dfJ_target) e exige melhora no J ===
auto_prune_greedy <- function(fit_obj, data, start_excl = NULL,
                              cluster_var = NULL, alpha = 0.05,
                              max_drops = 5, dfJ_target = 1,
                              improve_tol = 1e-6, verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  
  # ponto de partida = IVs excluídos de referência
  inst_excl <- if (is.null(start_excl)) .get_excluded_IVs(fit_obj) else unique(start_excl)
  history <- list()
  
  for (r in seq_len(max_drops)) {
    jtab <- .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var)
    min_p_old <- suppressWarnings(min(jtab$p_hansen, na.rm = TRUE))
    
    # parada 1: já OK
    if (is.finite(min_p_old) && min_p_old >= alpha) {
      if (verbose) message("Greedy: min p_Hansen = ", signif(min_p_old,3),
                           " ≥ ", alpha, ". Parando.")
      break
    }
    # parada 2: já não há sobre-ID suficiente
    if (any(!is.finite(jtab$df_J) | jtab$df_J < dfJ_target)) {
      if (verbose) message("Greedy: alguma eq tem df_J < ", dfJ_target,
                           " (sobre-ID insuficiente). Parando.")
      break
    }
    
    # equação mais crítica
    eq_bad <- jtab$eq[which.min(jtab$p_hansen)]
    
    # momentos na eq crítica com o conjunto atual
    fit_tmp <- fit_obj; fit_tmp$inst_terms <- inst_excl
    mt <- iv_moment_diag(fit_tmp, data, cluster_var)
    if (is.null(mt) || !nrow(mt)) {
      if (verbose) message("Greedy: não consegui calcular momentos. Parando.")
      break
    }
    cand_list <- mt |>
      dplyr::filter(eq == eq_bad) |>
      dplyr::arrange(dplyr::desc(abs_t)) |>
      dplyr::pull(inst)
    
    accepted <- FALSE
    for (cand in cand_list) {
      new_excl <- setdiff(inst_excl, cand)
      jnew <- .j_by_eq_with_excl(fit_obj, data, new_excl, cluster_var)
      
      # regra 1: manter df_J alvo em todas
      if (any(!is.finite(jnew$df_J) | jnew$df_J < dfJ_target)) next
      
      # regra 2: só aceita se melhora o pior p_Hansen ou já atinge alpha
      min_p_new <- suppressWarnings(min(jnew$p_hansen, na.rm = TRUE))
      if (!is.finite(min_p_new)) next
      if (min_p_new >= alpha || min_p_new > (min_p_old + improve_tol)) {
        if (verbose) message("Greedy rodada ", r, ": removendo {", cand,
                             "} (eq crítica = ", eq_bad, ").  min p: ",
                             signif(min_p_old,3), " -> ", signif(min_p_new,3))
        history[[length(history) + 1L]] <- list(drop = cand,
                                                eq = eq_bad,
                                                jtab_before = jtab,
                                                jtab_after  = jnew)
        inst_excl <- new_excl
        accepted <- TRUE
        break
      }
    }
    
    if (!accepted) {
      if (verbose) message("Greedy: nenhum candidato melhora min p mantendo df_J ≥ ",
                           dfJ_target, ". Parando.")
      break
    }
  }
  
  list(inst_excl_final = inst_excl,
       history   = history,
       jtab_final= .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var))
}


## =======================================================
## COMO USAR (continuação do seu fluxo) ------------------
## =======================================================

# Você já rodou a poda por C-tests (auto_prune_until_ok) e obteve 'prune_out'.
# Se ainda ficou J significativo (como no seu output), rode o greedy:

# 1) Checa situação atual (pós-C-tests)
j_pre <- jtests_by_eq_manual(
  prune_out$fit_for_diagnostics, df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  verbose = FALSE
)
cat("\nHansen-J (pré-greedy) -- min p:",
    signif(suppressWarnings(min(j_pre$p_hansen, na.rm=TRUE)),3), "\n")

# 2) Se precisar, aplica a poda gulosa iniciando dos excluídos finais do C-test
if (suppressWarnings(min(j_pre$p_hansen, na.rm = TRUE)) < 0.05) {
  greedy_out <- auto_prune_greedy(
    fit_obj     = fit_q,
    data        = df_iv,
    start_excl  = prune_out$inst_excl_final,
    cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
    alpha       = 0.05,
    max_drops   = 5,
    verbose     = TRUE
  )
  
  cat("\nIVs EXCLUÍDOS após greedy:\n  ",
      paste(greedy_out$inst_excl_final, collapse = " + "), "\n", sep = "")
  
  cat("\nJ (pós-greedy):\n"); print(greedy_out$jtab_final)
  
  # 3) Refit QUAIDS com o conjunto resultante da poda gulosa
  fit_q_new2 <- refit_quaids_try(fit_q, df_iv, inst_final_excl = greedy_out$inst_excl_final)
  
  if (!is.null(fit_q_new2)) {
    # Checklist pós-refit (o mesmo do seu post_refit_report)
    res_ci_new2 <- post_refit_report(
      fit_new    = fit_q_new2,
      data       = df_iv,
      cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
      R          = 2000, level = 0.95, seed = 123
    )
    # (opcional) heatmaps/export usando suas funções já criadas
  }
} else {
  message("Todos os Hansen-J já estão OK após C-tests; greedy não necessário.")
}

## =======================================================
## BUSCA COM TROCAS (swap) MANTENDO SOBRE-ID (df_J >= alvo)
## =======================================================

# --- util: pega IVs excluídos que existem no fit (pool total) ---
.get_excluded_pool <- function(fit_obj){
  rhs_all <- fit_obj$rhs_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst_all <- unique(fit_obj$inst_terms)
  inst_all <- inst_all[!is.na(inst_all) & nzchar(inst_all)]
  inst_all <- setdiff(inst_all, c("1","-1","0","+","~","|"))
  sort(unique(setdiff(inst_all, exogs)))
}

# --- util: constrói a lista de instrumentos = exógenos + (inst_excl escolhidos) ---
.compose_inst_from_excl <- function(fit_obj, inst_excl){
  rhs_all <- fit_obj$rhs_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  unique(c(exogs, inst_excl))
}

# --- J por equação dado um conjunto de IVs excluídos (usa seu .j_overid_manual) ---
.j_by_eq_with_excl <- function(fit_obj, data, inst_excl, cluster_var = NULL){
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  inst_terms_use <- .compose_inst_from_excl(fit_obj, inst_excl)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  rows <- lapply(seq_along(fit_obj$eq), function(i){
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm  <- try(AER::ivreg(f_yx, instruments = reformulate(inst_terms_use), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) {
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats)) {
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    
    # Sargan homoc (opcional) + Hansen robusto/cluster
    J_h <- try(AER::sargan(ivm), silent = TRUE)
    J_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$statistic)) else NA_real_
    p_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$p.value))   else NA_real_
    Jrob <- .j_overid_manual(ivm, cluster = cl)
    
    data.frame(eq = eqn, rank_X = rX, rank_Z = rZ, df_J = dfJ,
               J_sargan = J_s, p_sargan = p_s,
               J_hansen = unname(Jrob["J"]), p_hansen = unname(Jrob["p"]))
  })
  dplyr::bind_rows(rows)
}

# --- diagnóstico de momentos por instrumento (t approx. robusto/cluster) ---
iv_moment_diag <- function(fit_obj, data, cluster_var = NULL){
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_terms <- unique(fit_obj$inst_terms)
  inst_terms <- inst_terms[!is.na(inst_terms) & nzchar(inst_terms)]
  inst_terms <- setdiff(inst_terms, c("1","-1","0","+","~","|"))
  
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  out <- list()
  
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm  <- try(AER::ivreg(f_yx, instruments = reformulate(inst_terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) next
    XZu <- .get_XZ_u(ivm); if (is.null(XZu)) next
    Z <- XZu$Z; u <- as.numeric(XZu$u)
    n <- NROW(Z)
    
    # nomes das colunas de Z
    zn <- colnames(Z); if (is.null(zn)) zn <- paste0("z", seq_len(ncol(Z)))
    
    # g_i (n x k): cada coluna j é z_ij * u_i
    G_i <- Z * u
    # cluster robust S ≈ (1/n) Σ_c (Σ_i∈c g_i)(Σ_i∈c g_i)'
    if (!is.null(cl)) {
      Gc <- rowsum(G_i, group = as.factor(.align_cluster(ivm, cl)))
      S  <- crossprod(as.matrix(Gc)) / n
    } else {
      S  <- crossprod(G_i) / n
    }
    S <- (S + t(S))/2
    
    gbar <- colMeans(G_i)                  # k x 1
    se   <- sqrt(pmax(diag(S), 0) / n)     # se(gbar)
    tstat<- gbar / se
    
    out[[length(out)+1L]] <- data.frame(
      eq = eqn, inst = zn, gbar = as.numeric(gbar),
      se = as.numeric(se), t = as.numeric(tstat),
      abs_t = abs(as.numeric(tstat)),
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(data.frame(eq=character(), inst=character(),
                                      gbar=double(), se=double(), t=double(), abs_t=double()))
  dplyr::bind_rows(out)
}

# --- Greedy com trocas 1x1 (mais "drop puro") preservando df_J >= dfJ_target ---
auto_prune_swap_search <- function(fit_obj, data,
                                   start_excl,
                                   cluster_var = NULL,
                                   alpha = 0.05,
                                   dfJ_target = 1,
                                   max_rounds = 10,
                                   try_pure_drop = TRUE,
                                   verbose = TRUE){
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  pool_all <- .get_excluded_pool(fit_obj)
  
  # estado inicial
  inst_excl <- sort(unique(start_excl))
  best_j    <- .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var)
  best_minp <- suppressWarnings(min(best_j$p_hansen, na.rm = TRUE))
  
  if (verbose) message("Swap-search: min p inicial = ", signif(best_minp,3))
  
  for (r in seq_len(max_rounds)) {
    improved <- FALSE
    
    # equação crítica
    eq_bad <- best_j$eq[which.min(best_j$p_hansen)]
    
    # candidatos a sair/entrar
    drop_pool <- inst_excl
    add_pool  <- setdiff(pool_all, inst_excl)
    
    # tenta "drop puro" (se permitido)
    if (try_pure_drop) {
      for (cand_drop in drop_pool) {
        new_set <- setdiff(inst_excl, cand_drop)
        jnew <- .j_by_eq_with_excl(fit_obj, data, new_set, cluster_var)
        if (any(!is.finite(jnew$df_J) | jnew$df_J < dfJ_target)) next
        minp_new <- suppressWarnings(min(jnew$p_hansen, na.rm = TRUE))
        if (is.finite(minp_new) && (minp_new > best_minp || minp_new >= alpha)) {
          if (verbose) message("Round ", r, ": DROP {", cand_drop, "}  ",
                               signif(best_minp,3), " -> ", signif(minp_new,3))
          inst_excl <- sort(new_set); best_j <- jnew; best_minp <- minp_new
          improved <- TRUE; break
        }
      }
      if (improved) next
    }
    
    # tenta trocas 1x1, priorizando inst com maior |t| na eq crítica
    fit_tmp <- fit_obj; fit_tmp$inst_terms <- .compose_inst_from_excl(fit_obj, inst_excl)
    mt <- iv_moment_diag(fit_tmp, data, cluster_var)
    if (!nrow(mt)) break
    out_eq <- mt |> dplyr::filter(eq == eq_bad, inst %in% inst_excl) |>
      dplyr::arrange(dplyr::desc(abs_t)) |>
      dplyr::pull(inst)
    
    tried <- 0L
    for (cand_drop in out_eq) {
      for (cand_add in add_pool) {
        new_set <- sort(unique(c(setdiff(inst_excl, cand_drop), cand_add)))
        jnew <- .j_by_eq_with_excl(fit_obj, data, new_set, cluster_var)
        if (any(!is.finite(jnew$df_J) | jnew$df_J < dfJ_target)) next
        minp_new <- suppressWarnings(min(jnew$p_hansen, na.rm = TRUE))
        tried <- tried + 1L
        if (is.finite(minp_new) && (minp_new > best_minp || minp_new >= alpha)) {
          if (verbose) message("Round ", r, ": SWAP {", cand_drop, " -> ", cand_add, "}  ",
                               signif(best_minp,3), " -> ", signif(minp_new,3))
          inst_excl <- new_set; best_j <- jnew; best_minp <- minp_new
          improved <- TRUE; break
        }
      }
      if (improved) break
    }
    
    if (!improved) {
      if (verbose) message("Swap-search: nenhuma melhoria (tentativas=", tried, "). Parando.")
      break
    }
    if (best_minp >= alpha) {
      if (verbose) message("Swap-search: min p ≥ ", alpha, ". Parando.")
      break
    }
  }
  
  list(inst_excl_final = inst_excl,
       jtab_final      = best_j,
       min_p_final     = best_minp)
}

# ponto de partida = IVs EXCLUÍDOS finais pós C-tests
start_set <- prune_out$inst_excl_final

swap_out <- auto_prune_swap_search(
  fit_obj     = fit_q,
  data        = df_iv,
  start_excl  = start_set,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  alpha       = 0.05,
  dfJ_target  = 1,      # pode subir p/ 2 se quiser mais folga de sobre-ID
  max_rounds  = 10,
  try_pure_drop = TRUE, # tenta dropar 1 IV antes de procurar trocas
  verbose     = TRUE
)

cat("\nIVs EXCLUÍDOS (swap-search):\n  ",
    paste(swap_out$inst_excl_final, collapse = " + "), "\n", sep = "")
cat("\nJ (pós-swap-search):\n"); print(swap_out$jtab_final)

# Refit (se quiser seguir o seu fluxo)
fit_q_swap <- refit_quaids_try(fit_q, df_iv, inst_final_excl = swap_out$inst_excl_final)
if (!is.null(fit_q_swap)) {
  res_ci_swap <- post_refit_report(
    fit_new    = fit_q_swap,
    data       = df_iv,
    cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
    R          = 2000, level = 0.95, seed = 123
  )
}

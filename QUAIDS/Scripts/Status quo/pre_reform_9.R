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
  library(openxlsx)
})

## =======================================================
##  HELPERS BÁSICOS (robustos a cluster)
## =======================================================

# Alinha vetor de cluster às observações usadas em um modelo
.align_cluster <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)
  mf  <- stats::model.frame(model)
  rid <- rownames(mf)
  if (!is.null(names(cluster))) {
    return(cluster[rid])
  } else {
    idx <- suppressWarnings(as.integer(rid))
    if (all(is.finite(idx))) return(cluster[idx])
  }
  cluster
}

# Model matrix seguro e com limpeza de colunas singulares/constantes
.mm <- function(formula_or_terms, data, add_intercept = TRUE) {
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm  <- mm[, ok, drop = FALSE]
  q   <- qr(mm)
  if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}

# Extrai X, Z e resíduos de um ivreg
.get_XZ_u <- function(ivm) {
  X <- try(model.matrix(ivm, "regressors"),  silent = TRUE)
  Z <- try(model.matrix(ivm, "instruments"), silent = TRUE)
  if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NULL)
  list(X = X, Z = Z, u = residuals(ivm))
}

# J robusto/cluster (Arellano). Se cluster=NULL => HC0.
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

# Constrói lista de instrumentos: EXÓGENOS do RHS + IVs a manter (excluídos)
# - base_exog = RHS sem os ln_ (endógenos)
.build_inst_from_fit <- function(fit_obj, keep_excl = NULL) {
  rhs_all <- unique(fit_obj$rhs_terms)
  base_exog <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  inst <- unique(c(base_exog, keep_excl))
  inst <- inst[!is.na(inst) & nzchar(inst)]
  setdiff(inst, c("1","-1","0","+","~","|"))
}

# Pool de IVs "excluídos" possíveis = instrumentos do fit menos os exógenos do RHS
.get_excluded_pool <- function(fit_obj) {
  rhs_all <- unique(fit_obj$rhs_terms)
  base_exog <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  inst_all <- unique(fit_obj$inst_terms)
  out <- setdiff(inst_all, base_exog)
  out <- out[!grepl("^(1|-1|0|\\+|~|\\|)$", out)]
  unique(out)
}

# df_J mínimo para um conjunto de instrumentos (ranks de X e Z)
.min_dfJ_given_Z <- function(fit_obj, data, inst_terms){
  Xr <- vapply(fit_obj$eq, function(m) qr(.mm(delete.response(terms(formula(m))), data, TRUE))$rank, 1L)
  Zr <- vapply(fit_obj$eq, function(m) qr(.mm(reformulate(inst_terms, intercept = TRUE), data, TRUE))$rank, 1L)
  min(Zr - Xr, na.rm = TRUE)
}

# Calcula J eq-a-eq para um conjunto de IVs excluídos especificado
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
##  DIAGNÓSTICO GERAL
## =======================================================

# J por equação com o Z atual do fit (exógenos + IVs excluídos do fit)
jtests_by_eq_manual <- function(fit_obj, data, cluster_var = NULL, verbose = TRUE) {
  inst_excl <- .get_excluded_pool(fit_obj)
  .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var = cluster_var, verbose = verbose)
}

# Força dos IVs (teste em bloco dos IVs excluídos)
fs_block_F <- function(fit_obj, data, robust = TRUE){
  data <- as.data.frame(data)
  rhs_all <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)   # exógenos no Z
  iv_excl  <- setdiff(inst_all, controls)    # IVs excluídos
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

# partial R^2 dos IVs excluídos
partial_R2_exclIV <- function(fit_obj, data){
  rhs_all <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
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

# Wu-Hausman por equação (via control function) robusto/cluster
wu_hausman_cf <- function(fit_obj, data, cluster_var = NULL){
  rhs_all <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (!length(endo_ln)) stop("Nenhuma variável endógena detectada (ln_...).")
  
  ## 1ª etapa: resíduos das regressões de cada ln_preco em Z (controles + IVs excl.)
  Vhat <- lapply(endo_ln, function(y){
    f1 <- reformulate(c(controls, iv_excl), response = y)
    resid(lm(f1, data = data))
  })
  names(Vhat) <- paste0(endo_ln, "_res")
  df2 <- cbind(data, as.data.frame(Vhat))
  
  ## 2ª etapa: para cada equação K-1, inclui os resíduos e testa H0: coef(residuos)=0
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  out <- lapply(seq_along(eqs), function(i){
    f <- eqs[[i]]
    y <- as.character(f[[2]])
    xrhs <- attr(terms(f), "term.labels")
    f2 <- reformulate(c(xrhs, names(Vhat)), response = y)
    ols <- lm(f2, data = df2)
    
    ## VCOV robusto/cluster
    cl  <- if (!is.null(cluster_var) && cluster_var %in% names(df2)) df2[[cluster_var]] else NULL
    V   <- if (!is.null(cl)) sandwich::vcovCL(ols, cluster = cl, type = "HC3") else sandwich::vcovHC(ols, type = "HC3")
    
    ## Restrições como strings (evita mismatches de dimensão)
    coef_nms <- names(coef(ols))
    test_nms <- intersect(names(Vhat), coef_nms)  # só os resíduos que realmente entraram
    if (!length(test_nms)) {
      return(data.frame(eq = y, k_endog = length(endo_ln),
                        stat = NA_real_, df = NA_integer_, p_wu = NA_real_))
    }
    K <- paste0(test_nms, " = 0")
    lh <- try(car::linearHypothesis(ols, hypothesis.matrix = K, vcov = V, test = "Chisq"),
              silent = TRUE)
    if (inherits(lh, "try-error")) {
      return(data.frame(eq = y, k_endog = length(endo_ln),
                        stat = NA_real_, df = NA_integer_, p_wu = NA_real_))
    }
    data.frame(eq = y, k_endog = length(endo_ln),
               stat = as.numeric(lh$Chisq[2]),
               df   = as.integer(lh$Df[2]),
               p_wu = as.numeric(lh$`Pr(>Chisq)`[2]))
  })
  dplyr::bind_rows(out)
}

## =======================================================
##  C-TESTS ADAPTATIVOS + BLOCOS
## =======================================================

# Constroi blocos a partir de prefixos comuns
make_blocks_explicit <- function(fit_obj) {
  pool <- .get_excluded_pool(fit_obj)
  out <- list()
  # blocos padrão do seu caso:
  out$IV_d_uf_X <- grep("^IV_d_uf_X", pool, value = TRUE)
  out$iv_op     <- grep("^iv_op\\d+",  pool, value = TRUE)
  # remove blocos vazios
  out <- out[vapply(out, length, integer(1)) > 0]
  out
}

# C-test (diferença de Hansen) permitindo k removidos por bloco (até max_k_drop)
c_tests_blocks_adaptive <- function(fit_obj, data, blocks, cluster_var = NULL, max_k_drop = 2, verbose = FALSE){
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_full <- .get_excluded_pool(fit_obj)
  out <- list()
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    # modelo FULL (com todos os excluídos atuais do fit)
    iv_full <- try(ivreg(f_yx, instruments = reformulate(.build_inst_from_fit(fit_obj, inst_full)),
                         data = data), silent = TRUE)
    if (inherits(iv_full, "try-error")) next
    Jf <- .j_overid_manual(iv_full, cluster = cl)
    if (!is.finite(Jf["df"]) || Jf["df"] <= 0) { if (verbose) message("Eq ", eqn, ": just-ID"); next }
    # por bloco
    for (bk in names(blocks)) {
      B <- intersect(blocks[[bk]], inst_full)
      if (!length(B)) next
      kmax <- min(max_k_drop, length(B))
      for (k in 1:kmax) {
        combs <- utils::combn(B, k, simplify = FALSE)
        for (drop_set in combs) {
          keep <- setdiff(inst_full, drop_set)
          iv_r <- try(ivreg(f_yx, instruments = reformulate(.build_inst_from_fit(fit_obj, keep)),
                            data = data), silent = TRUE)
          if (inherits(iv_r, "try-error")) next
          Jr <- .j_overid_manual(iv_r, cluster = cl)
          C  <- as.numeric(Jr["J"] - Jf["J"])
          pC <- stats::pchisq(C, df = k, lower.tail = FALSE)
          out[[length(out)+1L]] <- data.frame(
            eq = eqn,
            block = bk,
            drop_set = paste(drop_set, collapse = " + "),
            k_removed = k,
            J_full = unname(Jf["J"]), dfJ_full = unname(Jf["df"]),
            J_restr = unname(Jr["J"]), dfJ_restr = unname(Jr["df"]),
            C_stat = C, df_C = k, p_C = pC, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(out)) return(tibble::tibble(
    eq=character(), block=character(), drop_set=character(), k_removed=integer(),
    J_full=double(), dfJ_full=integer(), J_restr=double(), dfJ_restr=integer(),
    C_stat=double(), df_C=integer(), p_C=double()
  ))
  dplyr::bind_rows(out)
}

# Sugerir poda a partir do ctab (p_C < alpha)
suggest_prune <- function(ctab, alpha = 0.05) {
  if (!nrow(ctab)) return(character(0))
  drops <- unique(unlist(strsplit(ctab$drop_set[which(ctab$p_C < alpha)], "\\s*\\+\\s*")))
  drops[ nzchar(drops) ]
}

# Após decidir o que dropar, devolve novos instrumentos (texto)
new_inst_after_prune <- function(fit_obj, to_drop, blocks) {
  inst_full <- .get_excluded_pool(fit_obj)
  keep <- setdiff(inst_full, to_drop)
  .build_inst_from_fit(fit_obj, keep_excl = keep)
}

## =======================================================
##  ROTINAS DE BUSCA: poda, greedy, swap, elevar df_J
## =======================================================

# Poda automática via C-tests até estabilizar (com guarda de df_J)
auto_prune_until_ok <- function(fit_obj, data, cluster_var = NULL,
                                alpha = 0.05, max_k_drop = 2, max_rounds = 5,
                                target_min_dfJ = 1, verbose = TRUE){
  current_excl <- .get_excluded_pool(fit_obj)
  fit_diag_obj <- fit_obj
  for (r in 1:max_rounds) {
    blocks <- make_blocks_explicit(fit_diag_obj)
    ctab <- c_tests_blocks_adaptive(fit_diag_obj, data, blocks, cluster_var, max_k_drop, verbose = FALSE)
    ctab <- dplyr::arrange(ctab, eq, p_C)
    to_drop_sets <- suggest_prune(ctab, alpha = alpha)
    if (!length(to_drop_sets)) { if (verbose) message("Nenhum C-test viável nesta rodada."); break }
    # tenta remover tudo; se violar df_J alvo, recua removendo 1 a 1
    drop_batch <- unique(to_drop_sets)
    inst_try <- setdiff(current_excl, drop_batch)
    inst_terms_try <- .build_inst_from_fit(fit_obj, keep_excl = inst_try)
    mindf <- .min_dfJ_given_Z(fit_obj, data, inst_terms_try)
    if (mindf >= target_min_dfJ) {
      if (verbose) message("Rodada ", r, ": removendo {", paste(drop_batch, collapse = ", "), "} (min df_J=", mindf, ").")
      current_excl <- inst_try
    } else {
      if (verbose) message("Rodada ", r, ": remoção em lote reduziria df_J abaixo de ", target_min_dfJ, ". Testando remoções unitárias.")
      for (z in drop_batch) {
        inst_try2 <- setdiff(current_excl, z)
        inst_terms_try2 <- .build_inst_from_fit(fit_obj, keep_excl = inst_try2)
        mindf2 <- .min_dfJ_given_Z(fit_obj, data, inst_terms_try2)
        if (mindf2 >= target_min_dfJ) {
          if (verbose) message("Rodada ", r, ": removendo {", z, "} (min df_J=", mindf2, ").")
          current_excl <- inst_try2
        } else if (verbose) {
          message("Rodada ", r, ": NÃO removi {", z, "} para não cair abaixo de df_J=", target_min_dfJ, ".")
        }
      }
    }
    # atualiza fit para diagnóstico (mesmo RHS; muda Z)
    inst_terms_new <- .build_inst_from_fit(fit_obj, keep_excl = current_excl)
    eqs <- lapply(fit_obj$eq, function(m) formula(m))
    fit_new <- systemfit(eqs, data = as.data.frame(data),
                         method = "3SLS",
                         inst = stats::reformulate(inst_terms_new, intercept = FALSE))
    fit_diag_obj$inst_terms <- inst_terms_new
    fit_diag_obj$fit <- fit_new
    fit_diag_obj$eq  <- fit_new$eq
  }
  list(inst_final = .build_inst_from_fit(fit_obj, keep_excl = current_excl),
       inst_excl_final = current_excl,
       fit_for_diagnostics = fit_diag_obj)
}

# Greedy: tenta remover 1 IV por vez (a partir de um set), max_drops, mantendo df_J>=1
auto_prune_greedy <- function(fit_obj, data, start_excl, cluster_var = NULL,
                              alpha = 0.05, max_drops = 5, verbose = TRUE){
  cur <- start_excl
  best <- .j_by_eq_with_excl(fit_obj, data, cur, cluster_var)
  best_minp <- suppressWarnings(min(best$p_hansen, na.rm = TRUE))
  for (t in 1:max_drops) {
    pool <- cur
    improved <- FALSE
    for (z in pool) {
      trial <- setdiff(cur, z)
      Jt <- .j_by_eq_with_excl(fit_obj, data, trial, cluster_var)
      if (any(Jt$df_J <= 0, na.rm = TRUE)) next
      minp <- suppressWarnings(min(Jt$p_hansen, na.rm = TRUE))
      if (isTRUE(verbose)) message("Greedy rodada ", t, ": removendo {", z, "} => ", signif(best_minp,3), " -> ", signif(minp,3))
      if (is.finite(minp) && minp > best_minp) {
        cur <- trial; best_minp <- minp; best <- Jt; improved <- TRUE
        break
      }
    }
    if (!improved) { if (verbose) message("Greedy: nenhum candidato melhora min p mantendo df_J ≥ 1. Parando."); break }
  }
  list(inst_excl_final = cur, jtab_final = best)
}

# Swap-search: tenta trocas 1-a-1 mantendo df_J alvo
auto_prune_swap_search <- function(fit_obj, data, start_excl, cluster_var = NULL,
                                   alpha = 0.05, dfJ_target = 1, max_rounds = 12,
                                   try_pure_drop = TRUE, verbose = TRUE){
  cur <- start_excl
  J0 <- .j_by_eq_with_excl(fit_obj, data, cur, cluster_var)
  best_minp <- suppressWarnings(min(J0$p_hansen, na.rm = TRUE))
  if (isTRUE(verbose)) message("Swap-search: min p inicial = ", signif(best_minp,3))
  pool_all <- .get_excluded_pool(fit_obj)
  for (r in 1:max_rounds) {
    improved <- FALSE
    # opcional: tentar um drop puro sem violar dfJ_target
    if (try_pure_drop) {
      for (z in cur) {
        trial <- setdiff(cur, z)
        Jt <- .j_by_eq_with_excl(fit_obj, data, trial, cluster_var)
        if (any(Jt$df_J < dfJ_target, na.rm = TRUE)) next
        minp <- suppressWarnings(min(Jt$p_hansen, na.rm = TRUE))
        if (is.finite(minp) && minp > best_minp) {
          if (isTRUE(verbose)) message("Round ", r, ": DROP {", z, "}  ", signif(best_minp,3), " -> ", signif(minp,3))
          cur <- trial; best_minp <- minp; improved <- TRUE; break
        }
      }
      if (improved) next
    }
    # tenta swaps 1-por-1
    add_cand <- setdiff(pool_all, cur)
    for (outz in cur) for (inz in add_cand) {
      trial <- sort(unique(c(setdiff(cur, outz), inz)))
      Jt <- .j_by_eq_with_excl(fit_obj, data, trial, cluster_var)
      if (any(Jt$df_J < dfJ_target, na.rm = TRUE)) next
      minp <- suppressWarnings(min(Jt$p_hansen, na.rm = TRUE))
      if (is.finite(minp) && minp > best_minp) {
        if (isTRUE(verbose)) message("Round ", r, ": SWAP {", outz, " -> ", inz, "}  ", signif(best_minp,3), " -> ", signif(minp,3))
        cur <- trial; best_minp <- minp; improved <- TRUE; break
      }
    }
    if (!improved) { if (isTRUE(verbose)) message("Swap-search: nenhuma melhoria (tentativas=", r, "). Parando."); break }
  }
  list(inst_excl_final = cur, jtab_final = .j_by_eq_with_excl(fit_obj, data, cur, cluster_var))
}

# Elevar df_J mínimo para um alvo adicionando poucos IVs e monitorando p_min
raise_dfJ_search <- function(fit_obj, data, start_excl, cluster_var = NULL,
                             target_dfJ = 2, max_add = 3, alpha = 0.05, verbose = TRUE){
  cur <- start_excl
  J0  <- .j_by_eq_with_excl(fit_obj, data, cur, cluster_var)
  min_dfJ <- suppressWarnings(min(J0$df_J, na.rm = TRUE))
  min_p   <- suppressWarnings(min(J0$p_hansen, na.rm = TRUE))
  if (isTRUE(verbose)) message("df_J inicial(min): ", min_dfJ, " | min p: ", signif(min_p,3))
  adds <- character(0)
  pool <- setdiff(.get_excluded_pool(fit_obj), cur)
  while (min_dfJ < target_dfJ && length(adds) < max_add && length(pool)) {
    # escolhe o add que maximiza o min p, respeitando dfJ_target
    cand_tab <- lapply(pool, function(z){
      trial <- sort(unique(c(cur, z)))
      Jt <- .j_by_eq_with_excl(fit_obj, data, trial, cluster_var)
      data.frame(add = z,
                 dfJ_min = suppressWarnings(min(Jt$df_J, na.rm = TRUE)),
                 p_min   = suppressWarnings(min(Jt$p_hansen, na.rm = TRUE)))
    }) %>% bind_rows() %>% arrange(desc(dfJ_min), desc(p_min))
    if (!nrow(cand_tab) || cand_tab$dfJ_min[1] < min(target_dfJ, 1)) break
    zstar <- cand_tab$add[1]
    cur <- sort(unique(c(cur, zstar))); adds <- c(adds, zstar)
    pool <- setdiff(pool, zstar)
    J0  <- .j_by_eq_with_excl(fit_obj, data, cur, cluster_var)
    min_dfJ <- suppressWarnings(min(J0$df_J, na.rm = TRUE))
    min_p   <- suppressWarnings(min(J0$p_hansen, na.rm = TRUE))
    if (isTRUE(verbose)) message("Add ", zstar, " => df_J(min)=", min_dfJ, " | min p=", signif(min_p,3))
  }
  list(inst_excl_final = cur, jtab_final = J0, adds = adds,
       min_dfJ_final = min_dfJ, min_p_final = min_p)
}

# Sensibilidade a especificação de cluster
compare_cluster_specs <- function(fit_obj, data, inst_excl, cluster_vars = c("psu","uf",NA)){
  tabs <- lapply(cluster_vars, function(cv){
    tab <- .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var = cv)
    tab$cluster <- if (is.na(cv)) "none" else as.character(cv)
    tab
  })
  bind_rows(tabs)
}

## =======================================================
##  MOMENT DIAGNOSTICS (quais IVs puxam o J)
## =======================================================
iv_moment_diag <- function(fit_obj_or_fitnew, data, cluster_var = NULL){
  # aceita objeto quaids_km1_fit (com eqs) retornado por refit_quaids_try
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_obj_or_fitnew$eq, function(m) formula(m))
  inst_terms <- fit_obj_or_fitnew$inst_terms
  out <- list()
  for (f in eqs) {
    eqn <- as.character(f[[2]])
    ivm <- AER::ivreg(f, instruments = reformulate(inst_terms), data = data)
    mats <- .get_XZ_u(ivm); if (is.null(mats)) next
    Z <- mats$Z; u <- as.numeric(mats$u); n <- NROW(Z)
    Zu <- Z * u
    if (!is.null(cl)) {
      G  <- rowsum(Zu, group = as.factor(.align_cluster(ivm, cl)))
      S  <- crossprod(as.matrix(G)) / n
    } else {
      S  <- crossprod(Zu) / n
    }
    gbar <- colMeans(Zu)
    se   <- sqrt(diag(S) / n)
    t    <- gbar / se
    tab  <- data.frame(eq = eqn, inst = colnames(Z), gbar = gbar, se = se, t = t, abs_t = abs(t))
    out[[length(out)+1L]] <- tab
  }
  bind_rows(out)
}

## =======================================================
##  REFIT & RELATÓRIO PÓS-REFIT
## =======================================================

# refit mantendo as equações (RHS) e trocando apenas instrumentos
refit_quaids_try <- function(fit_obj, data, inst_final_excl) {
  data <- as.data.frame(data)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  inst_terms_new <- .build_inst_from_fit(fit_obj, keep_excl = inst_final_excl)
  fit_new <- try(systemfit(eqs, data = data, method = "3SLS",
                           inst = stats::reformulate(inst_terms_new, intercept = FALSE)),
                 silent = TRUE)
  if (inherits(fit_new, "try-error")) return(NULL)
  # Recria fitted/observed completos (inclui omitida por adding-up dos K-1)
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

# Relatório pós-refit (J, F, pR2, Wu-Hausman, elasticidades/IC com fallback)
post_refit_report <- function(fit_new, data, cluster_var = NULL, R = 2000, level = 0.95, seed = 123){
  cat("\n== Over-ID (Hansen J) por equação ==\n")
  jt <- jtests_by_eq_manual(fit_new, data, cluster_var = cluster_var, verbose = FALSE)
  print(jt)
  if (any(jt$df_J <= 0, na.rm = TRUE)) {
    message("Aviso: há equações just/under-ID (df_J <= 0). Hansen J não é aplicável nessas eqs.")
  }
  cat("\n== Força dos IVs (F-HC3) ==\n")
  print(fs_block_F(fit_new, data, robust = TRUE))
  cat("\n== partial R^2 dos IVs excluídos ==\n")
  print(partial_R2_exclIV(fit_new, data))
  cat("\n== Wu-Hausman (control function) por equação ==\n")
  print(wu_hausman_cf(fit_new, data, cluster_var = cluster_var))
  
  # Elasticidades com IC: fallback elegante se a vcov não for P.D.
  if (exists("quaids_elasticities_ci", mode = "function")) {
    cat("\n== Elasticidades com IC (ponto: median/observed) ==\n")
    res_ci <- try(quaids_elasticities_ci(fit_new, at="medians", w_source="observed",
                                         R=R, level=level, seed=seed), silent = TRUE)
    if (!inherits(res_ci, "try-error") && is.list(res_ci) && !is.null(res_ci$expenditure)) {
      print(res_ci$expenditure)
      invisible(res_ci)
    } else {
      message("(Falha ao computar IC: ", if (!inherits(res_ci, "try-error")) "desconhecida" else attr(res_ci, "condition")$message,
              "). Reportando apenas ponto.")
      if (exists("quaids_elasticities", mode = "function")) {
        res_pt <- try(quaids_elasticities(fit_new, at="medians", w_source="observed"), silent = TRUE)
        if (!inherits(res_pt, "try-error") && is.list(res_pt) && !is.null(res_pt$expenditure)) {
          print(res_pt$expenditure)
          invisible(res_pt)
        } else invisible(NULL)
      } else invisible(NULL)
    }
  } else {
    message("\n(Elasticidades/IC puladas: função 'quaids_elasticities_ci' não encontrada.)")
    invisible(NULL)
  }
}

## =======================================================
##  HEATMAPS (opcional)
## =======================================================
tidy_elas_safe <- function(elas, share_names = NULL, price_names = NULL,
                           stat_set = c("est","se","z","p","lwr","upr")) {
  stopifnot(is.list(elas))
  first_mat <- elas[[ stat_set[ stat_set %in% names(elas) ][1] ]]
  if (is.null(first_mat)) stop("Nenhuma matriz encontrada entre: ", paste(stat_set, collapse=", "))
  rn <- rownames(first_mat); cn <- colnames(first_mat)
  if (is.null(rn) || is.null(cn)) stop("As matrizes precisam ter rownames e colnames.")
  if (is.null(share_names)) share_names <- rn
  if (is.null(price_names)) price_names <- cn
  rn <- trimws(rn); cn <- trimws(cn)
  share_names <- trimws(share_names); price_names <- trimws(price_names)
  use_rows <- intersect(share_names, rn)
  use_cols <- intersect(price_names,  cn)
  if (!length(use_cols)) stop("Nenhum dos 'price_names' aparece nas matrizes. Colunas: ", paste(cn, collapse=", "))
  if (!length(use_rows)) use_rows <- rn
  make_long <- function(stat){
    if (!stat %in% names(elas)) return(NULL)
    M <- as.matrix(elas[[stat]])[use_rows, use_cols, drop = FALSE]
    df <- as.data.frame(M, stringsAsFactors = FALSE)
    df$eq <- rownames(df)
    tidyr::pivot_longer(df, cols = dplyr::all_of(use_cols),
                        names_to = "price", values_to = stat)
  }
  pieces <- lapply(stat_set, make_long)
  pieces <- Filter(Negate(is.null), pieces)
  out <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("eq","price")), pieces)
  out$eq    <- factor(out$eq,    levels = share_names)
  out$price <- factor(out$price, levels = price_names)
  out[order(out$eq, out$price), ]
}

plot_heatmap <- function(df, title = "", limits = NULL, label_digits = 2,
                         use_ci_when_p_missing = TRUE) {
  stopifnot(all(c("eq","price","est") %in% names(df)))
  df2 <- df
  if (!"sig" %in% names(df2)) {
    if ("p" %in% names(df2) && any(!is.na(df2$p))) {
      df2$sig <- dplyr::case_when(
        is.na(df2$p)  ~ "",
        df2$p < 0.01  ~ "***",
        df2$p < 0.05  ~ "**",
        df2$p < 0.10  ~ "*",
        TRUE          ~ ""
      )
    } else if (use_ci_when_p_missing && all(c("lwr","upr") %in% names(df2))) {
      df2$sig <- dplyr::case_when(
        is.na(df2$lwr) | is.na(df2$upr) ~ "",
        df2$lwr > 0 | df2$upr < 0       ~ "*",
        TRUE                            ~ ""
      )
    } else df2$sig <- ""
  }
  df2$label <- sprintf(paste0("%.", label_digits, "f%s"), df2$est, df2$sig)
  ggplot(df2, aes(x = price, y = eq, fill = est)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(limits = limits, midpoint = 0, na.value = "grey85") +
    geom_text(aes(label = label), size = 3) +
    labs(title = title, x = "Price", y = "Good", fill = "Elasticity") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## =======================================================
##  COMO USAR — FASE FINAL “TOP TIER”
## =======================================================
## Pré-requisito: você já tem `fit_q` (quaids_km1_fit) e `df_iv`.

## (0) Sanity: quais ln_ (endógenos) estão no RHS
rhs_all <- fit_q$rhs_terms
cat("\nEndógenos no RHS (esperado fora de Z): ",
    paste(intersect(rhs_all, grep('^ln_', rhs_all, value=TRUE)), collapse=', '), "\n", sep="")

## (1) Poda adaptativa por C-tests até estabilizar (com guarda de df_J)
prune_out <- auto_prune_until_ok(
  fit_obj     = fit_q,
  data        = df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  alpha       = 0.05,
  max_k_drop  = 2,     # df_J=3 => testa 1 e 2
  max_rounds  = 5,
  target_min_dfJ = 1,
  verbose     = TRUE
)
cat("\nInst terms (FINAIS p/ diagnóstico/refit):\n",
    paste(prune_out$inst_final, collapse = " + "), "\n")

## (2) Diagnóstico com inst_finais (pré-refit)
cat("\n== Diagnóstico com inst_finais (pré-refit) ==\n")
print(jtests_by_eq_manual(prune_out$fit_for_diagnostics, df_iv,
                          cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"))

## (3) Refit com IVs podados
fit_q_new <- refit_quaids_try(fit_q, df_iv, inst_final_excl = prune_out$inst_excl_final)

## (4) Checklist pós-refit (tolerante a vcov não P.D.)
if (!is.null(fit_q_new)) {
  res_ci_new <- post_refit_report(
    fit_new    = fit_q_new,
    data       = df_iv,
    cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
    R          = 2000, level = 0.95, seed = 123
  )
}

## (5) Se ainda houver J significativo, rode greedy (mantendo df_J>=1)
j_pre <- jtests_by_eq_manual(
  prune_out$fit_for_diagnostics, df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  verbose = FALSE
)
if (suppressWarnings(min(j_pre$p_hansen, na.rm=TRUE)) < 0.05) {
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
}

## (6) Swap-search mantendo df_J alvo (use o set do seu último swap_out se já existir)
start_excl <- if (exists("greedy_out")) greedy_out$inst_excl_final else prune_out$inst_excl_final
swap_out <- auto_prune_swap_search(
  fit_obj      = fit_q,
  data         = df_iv,
  start_excl   = start_excl,
  cluster_var  = if ("psu" %in% names(df_iv)) "psu" else "uf",
  alpha        = 0.05,
  dfJ_target   = 1,
  max_rounds   = 10,
  try_pure_drop= TRUE,
  verbose      = TRUE
)
cat("\nIVs EXCLUÍDOS (swap-search):\n  ",
    paste(swap_out$inst_excl_final, collapse = " + "), "\n", sep = "")
cat("\nJ (pós-swap-search):\n"); print(swap_out$jtab_final)

## (7) (Opcional) Elevar df_J para 2 e tentar recuperar p com novo swap
lift <- raise_dfJ_search(
  fit_obj     = fit_q,
  data        = df_iv,
  start_excl  = swap_out$inst_excl_final,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  target_dfJ  = 2,
  max_add     = 3,
  alpha       = 0.05,
  verbose     = TRUE
)
cat("\n[raise_dfJ] adds: ", paste(lift$adds, collapse=", "), "\n", sep="")
cat("[raise_dfJ] min df_J=", lift$min_dfJ_final, " | min p=", signif(lift$min_p_final,3), "\n", sep="")
print(lift$jtab_final)

# ==== Helper: escolhe automaticamente a variável de cluster ====
pick_cluster <- function(data, candidates = c("psu","uf"), fallback = NULL) {
  nm <- names(data)
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else fallback   # retorna "psu", "uf" ou NULL
}

# Força dos IVs (teste em bloco dos IVs excluídos) com opção de cluster
# ========= PATCH 1: fs_block_F backward-compatible (aceita ... e usa cluster_var) =========
# ===== DROP-IN PATCH: fs_block_F robusto, tolerante a singularidade =====
# ===== DROP-IN PATCH: fs_block_F robusto, tolerante a singularidade =====
fs_block_F <- function(fit_obj, data, cluster_var = NULL) {
  data <- as.data.frame(data)
  
  # identifica endógenas (ln_), controles (exógenos no RHS) e IVs excluídos
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)       # exógenos em Z (também no RHS)
  iv_excl  <- setdiff(inst_all, controls)        # IVs excluídos = alvo do teste
  
  if (length(iv_excl) == 0L)
    stop("Não há IVs excluídos detectados.")
  
  safe_F_one <- function(y) {
    # 1) Regressão do 1º estágio "ampliada": y ~ controles + IVs excluídos
    f <- reformulate(c(controls, iv_excl), response = y)
    m <- stats::lm(f, data = data)
    
    # Extrai X, resíduos, dims
    X <- stats::model.matrix(m)
    u <- stats::residuals(m)
    n <- NROW(X)
    
    # Bread ~ (X'X)^-1 com pseudoinversa + ridge leve
    XtX <- crossprod(X)
    mu  <- mean(diag(XtX)); if (!is.finite(mu) || mu <= 0) mu <- 1
    XtX <- XtX + diag(1e-10 * mu, ncol(XtX))
    XtX_inv <- MASS::ginv(XtX)
    
    # Meat robusto:
    #   sem cluster -> X'u^2 X  == crossprod(X * u)
    #   com cluster -> sum_g (X_g' u_g)(X_g' u_g)'  == crossprod(rowsum(X * u, g))
    Xu <- X * as.numeric(u)
    if (!is.null(cluster_var) && cluster_var %in% names(data)) {
      g     <- rowsum(Xu, group = as.factor(data[[cluster_var]]))
      meat  <- crossprod(as.matrix(g))
    } else {
      meat  <- crossprod(Xu)
    }
    
    # VCOV robusto
    V <- XtX_inv %*% meat %*% XtX_inv
    
    # Monta H0: betas dos IVs excluídos == 0 (apenas os que ficaram no modelo)
    b     <- stats::coef(m)
    cn    <- names(b)
    iv_in <- intersect(iv_excl, cn)
    if (!length(iv_in))
      return(data.frame(endog = y, F = NA_real_, df1 = NA_integer_, df2 = NA_integer_, p = NA_real_))
    
    R <- matrix(0, nrow = length(iv_in), ncol = length(cn),
                dimnames = list(iv_in, cn))
    R[cbind(iv_in, iv_in)] <- 1
    
    rb    <- as.matrix(R %*% b)
    RVRT  <- R %*% V %*% t(R)
    RVRT  <- (RVRT + t(RVRT)) / 2
    mu2   <- mean(diag(RVRT)); if (!is.finite(mu2) || mu2 <= 0) mu2 <- 1
    RVRT  <- RVRT + diag(1e-10 * mu2, nrow(RVRT))
    RVRTi <- MASS::ginv(RVRT)
    
    Wald <- drop(t(rb) %*% RVRTi %*% rb)
    r    <- length(iv_in)
    Fst  <- Wald / r
    
    # df2 ~ n - rank(X); se não for válido, cai pra qui-quadrado assintótico
    rX   <- qr(X)$rank
    df2  <- n - rX
    if (!is.finite(df2) || df2 <= 0) {
      p <- stats::pchisq(Wald, df = r, lower.tail = FALSE)
      data.frame(endog = y, F = Fst, df1 = r, df2 = NA_integer_, p = p)
    } else {
      p <- stats::pf(Fst, r, df2, lower.tail = FALSE)
      data.frame(endog = y, F = Fst, df1 = r, df2 = df2, p = p)
    }
  }
  
  do.call(rbind, lapply(endo_ln, safe_F_one))
}

# ========= PATCH 2: post_refit_report chamando fs_block_F com cluster_var =========
post_refit_report <- function(fit_new, data, cluster_var = NULL, R = 2000, level = 0.95, seed = 123){
  cat("\n== Over-ID (Hansen J) por equação ==\n")
  jt <- jtests_by_eq_manual(fit_new, data, cluster_var = cluster_var, verbose = FALSE)
  print(jt)
  if (any(jt$df_J <= 0, na.rm = TRUE)) {
    message("Aviso: há equações just/under-ID (df_J <= 0). Hansen J não é aplicável nessas eqs.")
  }
  
  cat("\n== Força dos IVs (F-HC3/CL) ==\n")
  print(fs_block_F(fit_new, data, cluster_var = cluster_var))  # << sem robust
  
  cat("\n== partial R^2 dos IVs excluídos ==\n")
  print(partial_R2_exclIV(fit_new, data))
  
  cat("\n== Wu-Hausman (control function) por equação ==\n")
  print(wu_hausman_cf(fit_new, data, cluster_var = cluster_var))
  
  # Elasticidades com IC (fallback elegante)
  if (exists("quaids_elasticities_ci", mode = "function")) {
    cat("\n== Elasticidades com IC (ponto: median/observed) ==\n")
    res_ci <- try(quaids_elasticities_ci(fit_new, at="medians", w_source="observed",
                                         R=R, level=level, seed=seed), silent = TRUE)
    if (!inherits(res_ci, "try-error") && is.list(res_ci) && !is.null(res_ci$expenditure)) {
      print(res_ci$expenditure)
      invisible(res_ci)
    } else {
      message("(Falha ao computar IC: ",
              if (!inherits(res_ci, "try-error")) "desconhecida" else attr(res_ci, "condition")$message,
              "). Reportando apenas ponto.")
      if (exists("quaids_elasticities", mode = "function")) {
        res_pt <- try(quaids_elasticities(fit_new, at="medians", w_source="observed"), silent = TRUE)
        if (!inherits(res_pt, "try-error") && is.list(res_pt) && !is.null(res_pt$expenditure)) {
          print(res_pt$expenditure)
          invisible(res_pt)
        } else invisible(NULL)
      } else invisible(NULL)
    }
  } else {
    message("\n(Elasticidades/IC puladas: função 'quaids_elasticities_ci' não encontrada.)")
    invisible(NULL)
  }
}

cat("\n== Força dos IVs (F-HC3/CL) ==\n")
# F-HC3 para o set final (df_J ≥ 2)
swap_dfJ2 <- auto_prune_swap_search(
  fit_obj      = fit_q,
  data         = df_iv,
  start_excl   = lift$inst_excl_final,
  cluster_var  = if ("psu" %in% names(df_iv)) "psu" else "uf",
  alpha        = 0.05,
  dfJ_target   = 2,
  max_rounds   = 12,
  try_pure_drop= FALSE,
  verbose      = TRUE
)
cat("\nIVs EXCLUÍDOS (df_J≥2 pós-swap):\n  ",
    paste(swap_dfJ2$inst_excl_final, collapse=" + "), "\n", sep="")
cat("\nJ (df_J≥2 pós-swap):\n"); print(swap_dfJ2$jtab_final)

fit_q_dfJ2 <- refit_quaids_try(fit_q, df_iv, inst_final_excl = swap_dfJ2$inst_excl_final)
print(fs_block_F(fit_q_dfJ2, df_iv, cluster_var = pick_cluster(df_iv)))

print(fs_block_F(fit_q_new, df_iv, cluster_var = pick_cluster(df_iv)))
print(fs_block_F(prune_out$fit_for_diagnostics, df_iv, cluster_var = pick_cluster(df_iv)))
summarize_overid <- function(jtab){
  data.frame(
    min_dfJ = suppressWarnings(min(jtab$df_J, na.rm = TRUE)),
    min_p   = suppressWarnings(min(jtab$p_hansen, na.rm = TRUE))
  )
}
jt <- jtests_by_eq_manual(fit_q_dfJ2, df_iv, cluster_var = pick_cluster(df_iv), verbose = FALSE)
print(jt); print(summarize_overid(jt))
iv_moment_diag <- function(fit_obj_or_fitnew, data, cluster_var = NULL, only_excluded = FALSE){
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_obj_or_fitnew$eq, function(m) formula(m))
  inst_terms <- fit_obj_or_fitnew$inst_terms
  if (isTRUE(only_excluded)) {
    rhs_all <- fit_obj_or_fitnew$rhs_terms
    base_exog <- setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
    inst_terms <- setdiff(inst_terms, base_exog)
  }
  out <- list()
  for (f in eqs) {
    eqn <- as.character(f[[2]])
    ivm <- AER::ivreg(f, instruments = reformulate(inst_terms), data = data)
    mats <- .get_XZ_u(ivm); if (is.null(mats)) next
    Z <- mats$Z; u <- as.numeric(mats$u); n <- NROW(Z)
    Zu <- Z * u
    S  <- if (!is.null(cl)) {
      G <- rowsum(Zu, group = as.factor(.align_cluster(ivm, cl))); crossprod(as.matrix(G))/n
    } else crossprod(Zu)/n
    gbar <- colMeans(Zu); se <- sqrt(diag(S)/n); t <- gbar / se
    tab  <- data.frame(eq = eqn, inst = colnames(Z), gbar = gbar, se = se, t = t, abs_t = abs(t))
    out[[length(out)+1L]] <- tab
  }
  res <- dplyr::bind_rows(out)
  res <- dplyr::filter(res, inst != "(Intercept)")
  res
}
diag_tab <- iv_moment_diag(fit_q_dfJ2, df_iv, cluster_var = pick_cluster(df_iv), only_excluded = TRUE)


## (8) Moment diagnostics (quais IVs mais “tensionam” os momentos)

cat("\n== Moment diagnostics (abs t) no set final ==\n")
diag_tab <- iv_moment_diag(fit_q_dfJ2, df_iv, cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf")
print(dplyr::arrange(diag_tab, dplyr::desc(abs_t)) %>% dplyr::group_by(eq) %>% dplyr::slice_head(n = 8))

## (9) Relatório final (J/F/pR2/Wu + elasticidades/IC)
if (!is.null(fit_q_dfJ2)) {
  res_ci_final <- post_refit_report(
    fit_new    = fit_q_dfJ2,
    data       = df_iv,
    cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
    R=2000, level=0.95, seed=123
  )
}
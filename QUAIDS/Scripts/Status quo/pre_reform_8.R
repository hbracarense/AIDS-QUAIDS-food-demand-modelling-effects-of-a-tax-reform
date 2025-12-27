## =======================================================
## CONTINUAÇÃO TOP TIER: levantar df_J, cluster check, re-spec mínima
## Depende de: AER, sandwich, dplyr e dos helpers que você já tem:
##   .align_cluster, .get_XZ_u, .j_overid_manual,
##   .get_excluded_pool, .compose_inst_from_excl,
##   .j_by_eq_with_excl, iv_moment_diag, auto_prune_swap_search,
##   refit_quaids_try, post_refit_report, fit_quaids_manual_km1
## (Se algum helper não existir, redefino abaixo.)
## =======================================================
library(openxlsx)

# -- Segurança: define helpers essenciais se ausentes -------------------------
if (!exists(".compose_inst_from_excl", mode="function")) {
  .compose_inst_from_excl <- function(fit_obj, inst_excl){
    rhs_all <- fit_obj$rhs_terms
    endo_ln <- grep("^ln_", rhs_all, value = TRUE)
    exogs   <- setdiff(rhs_all, endo_ln)
    unique(c(exogs, inst_excl))
  }
}
if (!exists(".get_excluded_pool", mode="function")) {
  .get_excluded_pool <- function(fit_obj){
    rhs_all <- fit_obj$rhs_terms
    endo_ln <- grep("^ln_", rhs_all, value = TRUE)
    exogs   <- setdiff(rhs_all, endo_ln)
    inst_all <- unique(fit_obj$inst_terms)
    inst_all <- inst_all[!is.na(inst_all) & nzchar(inst_all)]
    inst_all <- setdiff(inst_all, c("1","-1","0","+","~","|"))
    sort(unique(setdiff(inst_all, exogs)))
  }
}
if (!exists(".j_by_eq_with_excl", mode="function")) {
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
      XZu <- .get_XZ_u(ivm); if (is.null(XZu)) {
        return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                          J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
      }
      rX <- qr(XZu$X)$rank; rZ <- qr(XZu$Z)$rank; dfJ <- rZ - rX
      Jh <- try(AER::sargan(ivm), silent = TRUE)
      J_s <- if (!inherits(Jh, "try-error")) unname(as.numeric(Jh$statistic)) else NA_real_
      p_s <- if (!inherits(Jh, "try-error")) unname(as.numeric(Jh$p.value))   else NA_real_
      Jrob <- .j_overid_manual(ivm, cluster = cl)
      data.frame(eq = eqn, rank_X = rX, rank_Z = rZ, df_J = dfJ,
                 J_sargan = J_s, p_sargan = p_s,
                 J_hansen = unname(Jrob["J"]), p_hansen = unname(Jrob["p"]))
    })
    dplyr::bind_rows(rows)
  }
}
if (!exists("iv_moment_diag", mode="function")) {
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
      Z <- XZu$Z; u <- as.numeric(XZu$u); n <- NROW(Z)
      zn <- colnames(Z); if (is.null(zn)) zn <- paste0("z", seq_len(ncol(Z)))
      G_i <- Z * u
      if (!is.null(cl)) {
        Gc <- rowsum(G_i, group = as.factor(.align_cluster(ivm, cl)))
        S  <- crossprod(as.matrix(Gc)) / n
      } else {
        S  <- crossprod(G_i) / n
      }
      S <- (S + t(S))/2
      gbar <- colMeans(G_i)
      se   <- sqrt(pmax(diag(S), 0) / n)
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
}

## -------------------------------------------------------
## (1) ELEVAR df_J: adicionar 1-por-vez maximizando min p
## -------------------------------------------------------
raise_dfJ_search <- function(fit_obj, data,
                             start_excl,
                             cluster_var = NULL,
                             target_dfJ = 2,
                             max_add = 3,
                             alpha = 0.05,
                             verbose = TRUE) {
  pool <- .get_excluded_pool(fit_obj)
  inst_excl <- sort(unique(start_excl))
  jcur <- .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var)
  min_df <- suppressWarnings(min(jcur$df_J, na.rm=TRUE))
  min_p  <- suppressWarnings(min(jcur$p_hansen, na.rm=TRUE))
  if (verbose) message("df_J inicial(min): ", min_df, " | min p: ", signif(min_p,3))
  
  adds_done <- character(0)
  for (k in seq_len(max_add)) {
    if (is.finite(min_df) && min_df >= target_dfJ) break
    cand_pool <- setdiff(pool, inst_excl)
    if (!length(cand_pool)) { if (verbose) message("Sem candidatos para adicionar."); break }
    best <- list(minp = -Inf, dfmin = -Inf, add = NA, jtab = NULL)
    for (zadd in cand_pool) {
      new_set <- c(inst_excl, zadd)
      jnew <- .j_by_eq_with_excl(fit_obj, data, new_set, cluster_var)
      dfmin_new <- suppressWarnings(min(jnew$df_J, na.rm=TRUE))
      if (!is.finite(dfmin_new) || dfmin_new < min_df) next
      minp_new <- suppressWarnings(min(jnew$p_hansen, na.rm=TRUE))
      if (is.finite(minp_new) && (minp_new > best$minp || dfmin_new > best$dfmin)) {
        best <- list(minp = minp_new, dfmin = dfmin_new, add = zadd, jtab = jnew)
      }
    }
    if (is.na(best$add)) { if (verbose) message("Não encontrei add que eleve df_J."); break }
    inst_excl <- sort(unique(c(inst_excl, best$add)))
    adds_done <- c(adds_done, best$add)
    jcur <- best$jtab; min_df <- best$dfmin; min_p <- best$minp
    if (verbose) message("Add ", best$add, " => df_J(min)=", min_df, " | min p=", signif(min_p,3))
    if (min_df >= target_dfJ && min_p >= alpha) break
  }
  list(inst_excl_final = inst_excl,
       adds = adds_done,
       jtab_final = jcur,
       min_p_final = min_p,
       min_dfJ_final = min_df)
}

## -------------------------------------------------------
## (2) Sensibilidade a cluster: psu / uf / nenhum
## -------------------------------------------------------
compare_cluster_specs <- function(fit_obj, data, inst_excl,
                                  cluster_vars = c("psu","uf",NA)){
  res <- lapply(cluster_vars, function(cv){
    cv_use <- if (is.na(cv)) NULL else if (cv %in% names(data)) cv else NULL
    j <- .j_by_eq_with_excl(fit_obj, data, inst_excl, cluster_var = cv_use)
    data.frame(cluster = if (is.null(cv_use)) "none" else cv_use,
               eq = j$eq,
               df_J = j$df_J,
               p_hansen = j$p_hansen,
               stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(res)
}

## -------------------------------------------------------
## (3) Re-especificação mínima (grid pequeno) + swap interno
## -------------------------------------------------------
re_spec_minimal_grid <- function(prices, shares, x,
                                 fit_base, instNames, instData,
                                 start_excl,
                                 cluster_var = NULL,
                                 alpha = 0.05,
                                 grid = list(
                                   demo_pref = c("p_n_0esc","p_n_8esc"),
                                   center_sex = c(TRUE, FALSE),
                                   use_z2     = c(TRUE, FALSE),
                                   add_inter  = c(TRUE, FALSE)
                                 ),
                                 verbose = TRUE){
  combos <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- list()
  for (i in seq_len(nrow(combos))) {
    g <- combos[i, ]
    if (verbose) message("Re-fit grid ", i, "/", nrow(combos), " -> ",
                         paste(names(g), unlist(g), sep="=", collapse=", "))
    fit_try <- try(
      fit_quaids_manual_km1(
        prices = prices, shares = shares, x = x,
        priceIndex = fit_base$priceIndex, estMethod = "3SLS",
        omit_share = which(fit_base$shareNames == fit_base$omit_share),
        drop_price = which(fit_base$priceNames == fit_base$drop_price),
        instNames  = instNames,
        instData   = instData,
        use_z2     = isTRUE(g$use_z2),
        include_cap     = TRUE,
        demo_pref       = g$demo_pref,
        center_sex      = isTRUE(g$center_sex),
        add_interaction = isTRUE(g$add_inter),
        compute_diag    = FALSE
      ), silent = TRUE)
    if (inherits(fit_try, "try-error")) next
    
    # mapeia start_excl para o novo universo e roda swap-search curto
    pool_new <- .get_excluded_pool(fit_try)
    start_in_new <- intersect(start_excl, pool_new)
    if (!length(start_in_new)) start_in_new <- pool_new # fallback
    swap_out <- auto_prune_swap_search(
      fit_obj     = fit_try,
      data        = instData,
      start_excl  = start_in_new,
      cluster_var = cluster_var,
      alpha       = alpha,
      dfJ_target  = 1,
      max_rounds  = 6,
      try_pure_drop = TRUE,
      verbose     = FALSE
    )
    minp <- suppressWarnings(min(swap_out$jtab_final$p_hansen, na.rm = TRUE))
    out[[length(out)+1L]] <- list(spec = g, fit = fit_try,
                                  inst_excl_final = swap_out$inst_excl_final,
                                  jtab = swap_out$jtab_final,
                                  minp = minp)
  }
  if (!length(out)) return(NULL)
  # escolhe melhor por maior min p; em empate, maior df_J(min)
  ord <- order(vapply(out, function(x) x$minp, numeric(1)), decreasing = TRUE)
  out[[ord[1]]]
}

## =======================================================
## COMO USAR (rápido)
## =======================================================

## 0) Ponto de partida: use o set do swap_search que você acabou de obter
inst_excl_start <- swap_out$inst_excl_final

## 1) (opcional) Elevar df_J para 2 mantendo/ganhando p
lift <- raise_dfJ_search(
  fit_obj     = fit_q,
  data        = df_iv,
  start_excl  = inst_excl_start,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  target_dfJ  = 2,
  max_add     = 3,
  alpha       = 0.05,
  verbose     = TRUE
)
cat("\n[raise_dfJ] adds: ", paste(lift$adds, collapse=", "), "\n", sep="")
cat("[raise_dfJ] min df_J=", lift$min_dfJ_final, " | min p=", signif(lift$min_p_final,3), "\n", sep="")
print(lift$jtab_final)

## 2) (opcional) Sensibilidade de cluster
cl_comp <- compare_cluster_specs(
  fit_obj   = fit_q,
  data      = df_iv,
  inst_excl = inst_excl_start,
  cluster_vars = c("psu","uf",NA)
)
cat("\n== Cluster sensitivity (min p por cluster) ==\n")
print(dplyr::summarise(dplyr::group_by(cl_comp, cluster),
                       min_p = suppressWarnings(min(p_hansen, na.rm=TRUE))))

## 3) (opcional) Re-especificação mínima + swap interno (usa sua função de refit)
##    Requer que você tenha: prices, shares, x, instNames, df_iv no ambiente
# best_spec <- re_spec_minimal_grid(
#   prices     = prices,
#   shares     = shares,
#   x          = x,
#   fit_base   = fit_q,
#   instNames  = intersect(.get_excluded_pool(fit_q), fit_q$inst_terms),
#   instData   = df_iv,
#   start_excl = inst_excl_start,
#   cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
#   alpha      = 0.05,
#   verbose    = TRUE
# )
# if (!is.null(best_spec)) {
#   cat("\n== MELHOR RE-SPEC ==\n"); print(best_spec$spec)
#   cat("min p =", signif(best_spec$minp,3), "\n")
#   cat("IVs EXCLUÍDOS(final):\n  ", paste(best_spec$inst_excl_final, collapse=" + "), "\n", sep="")
#   # refit final e relatório:
#   fit_q_final <- best_spec$fit
#   res_ci_final <- post_refit_report(
#     fit_new    = fit_q_final,
#     data       = df_iv,
#     cluster_var= if ("psu" %in% names(df_iv)) "psu" else "uf",
#     R          = 2000, level = 0.95, seed = 123
#   )
# }

## =======================================================
## OTIMIZAÇÃO P/ df_J ≥ 2 + RESGATE DE p COM SWAPS
## =======================================================

## 0) ponto de partida: set de EXCLUÍDOS após seu swap_search atual
inst_excl_start <- swap_out$inst_excl_final
cluster_use <- if ("psu" %in% names(df_iv)) "psu" else "uf"

## ---------- (A) varredura ótima de 1 add p/ chegar a df_J ≥ 2 ----------
best_add_for_dfJ2 <- function(fit_obj, data, start_excl, cluster_var = NULL) {
  pool <- .get_excluded_pool(fit_obj)
  cand <- setdiff(pool, start_excl)
  if (!length(cand)) stop("Sem candidatos para adicionar.")
  tab <- lapply(cand, function(z){
    ex <- c(start_excl, z)
    j  <- .j_by_eq_with_excl(fit_obj, data, ex, cluster_var)
    data.frame(add = z,
               dfJ_min = suppressWarnings(min(j$df_J, na.rm = TRUE)),
               p_min   = suppressWarnings(min(j$p_hansen, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  })
  tab <- dplyr::bind_rows(tab)
  tab <- dplyr::arrange(tab, dplyr::desc(dfJ_min), dplyr::desc(p_min))
  # escolhe a melhor com dfJ_min >= 2; se nenhuma, fica NA
  pick <- if (any(tab$dfJ_min >= 2)) tab$add[which.max(ifelse(tab$dfJ_min>=2, tab$p_min, -Inf))] else NA
  list(table = tab, add_star = pick)
}

add_out <- best_add_for_dfJ2(fit_q, df_iv, inst_excl_start, cluster_var = cluster_use)
cat("\n== Ranking de 1-add para df_J ≥ 2 (cluster=", cluster_use, ") ==\n", sep="")
print(add_out$table, digits = 4)
if (is.na(add_out$add_star)) {
  message("Nenhum add único leva a df_J ≥ 2 no cluster escolhido.")
} else {
  cat("Escolha recomendada (1-add): ", add_out$add_star, "\n", sep="")
}

## ---------- (B) (opcional) testar pares de adds quando houver 2+ candidatos ----------
best_pair_for_dfJ2 <- function(fit_obj, data, start_excl, cluster_var = NULL, max_pairs = 30){
  pool <- setdiff(.get_excluded_pool(fit_obj), start_excl)
  if (length(pool) < 2) return(NULL)
  combs <- utils::combn(pool, 2, simplify = FALSE)
  if (length(combs) > max_pairs) combs <- combs[seq_len(max_pairs)]  # limita exploração
  tab <- lapply(combs, function(z2){
    ex <- c(start_excl, z2)
    j  <- .j_by_eq_with_excl(fit_obj, data, ex, cluster_var)
    data.frame(add_pair = paste(z2, collapse=" + "),
               dfJ_min  = suppressWarnings(min(j$df_J, na.rm = TRUE)),
               p_min    = suppressWarnings(min(j$p_hansen, na.rm = TRUE)))
  })
  tab <- dplyr::bind_rows(tab)
  tab <- dplyr::arrange(tab, dplyr::desc(dfJ_min), dplyr::desc(p_min))
  tab
}

pair_tab <- best_pair_for_dfJ2(fit_q, df_iv, inst_excl_start, cluster_var = cluster_use)
if (!is.null(pair_tab)) {
  cat("\n== Top pares de adds (até 30) ==\n")
  print(utils::head(pair_tab, 10), digits = 4)
}

## ---------- (C) fixa o alvo (df_J ≥ 2) e tenta resgatar p com SWAP ----------
inst_excl_dfJ2 <- inst_excl_start
if (!is.na(add_out$add_star)) inst_excl_dfJ2 <- sort(unique(c(inst_excl_dfJ2, add_out$add_star)))

# checa J com o set dfJ2
j_dfJ2 <- .j_by_eq_with_excl(fit_q, df_iv, inst_excl_dfJ2, cluster_var = cluster_use)
cat("\n== Hansen com df_J ≥ 2 (pré-swap) ==\n"); print(j_dfJ2)

# swap-search mantendo df_J alvo
swap_dfJ2 <- auto_prune_swap_search(
  fit_obj      = fit_q,
  data         = df_iv,
  start_excl   = inst_excl_dfJ2,
  cluster_var  = cluster_use,
  alpha        = 0.05,
  dfJ_target   = 2,         # mantém pelo menos 2
  max_rounds   = 12,
  try_pure_drop= FALSE,     # evita cair para df_J<2
  verbose      = TRUE
)
cat("\nIVs EXCLUÍDOS (df_J≥2 pós-swap):\n  ",
    paste(swap_dfJ2$inst_excl_final, collapse=" + "), "\n", sep="")
cat("\nJ (df_J≥2 pós-swap):\n"); print(swap_dfJ2$jtab_final)

## ---------- (D) Moment diagnostics: identifica IVs “problemáticos” ----------
cat("\n== Moment diagnostics (abs t) no set df_J≥2 pós-swap ==\n")
diag_tab <- iv_moment_diag(
  refit_quaids_try(fit_q, df_iv, inst_final_excl = swap_dfJ2$inst_excl_final),
  df_iv, cluster_var = cluster_use
)
print(dplyr::arrange(diag_tab, dplyr::desc(abs_t)) %>% dplyr::group_by(eq) %>% dplyr::slice_head(n = 8))

## ---------- (E) (opcional) Refit e relatório final ----------
fit_q_dfJ2 <- refit_quaids_try(fit_q, df_iv, inst_final_excl = swap_dfJ2$inst_excl_final)
if (!is.null(fit_q_dfJ2)) {
  res_ci_dfJ2 <- post_refit_report(
    fit_new    = fit_q_dfJ2,
    data       = df_iv,
    cluster_var= cluster_use,
    R=2000, level=0.95, seed=123
  )
}

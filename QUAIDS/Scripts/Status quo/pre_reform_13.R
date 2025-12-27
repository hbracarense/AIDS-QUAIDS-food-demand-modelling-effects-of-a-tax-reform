## ===============================================================
## TOP TIER – ROUND 2/3: busca ampliada + sensibilidade + C-tests
## ===============================================================
library(openxlsx)
## ---------- (0) Setup ----------
stopifnot(exists("fit_q"), exists("df_iv"))
cl_var <- pick_cluster(df_iv)
eq_names <- vapply(fit_q$eq, function(m) as.character(formula(m)[[2]]), character(1))

## seu conjunto base (ponto de partida)
excl_now <- c("IV_d_uf_X11","IV_d_uf_X12","IV_d_uf_X13","IV_d_uf_X15","iv_op01","iv_op02","iv_op03")

## (PATCH) Moment diagnostics: evitar warnings mesmo com only_excluded=TRUE
iv_moment_diag_auto <- (function(orig_fun){
  function(fit_obj, data, cluster_var=NULL, only_excluded=FALSE){
    data <- as.data.frame(data)
    cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
    eqs <- lapply(fit_obj$eq, function(m) formula(m))
    inst_excl_by_eq <- attr(fit_obj, "inst_excl_by_eq")
    base_exog <- get_base_exog(fit_obj)
    use_eqspec <- is.list(inst_excl_by_eq)
    out <- list()
    for (f in eqs) {
      eqn <- as.character(f[[2]])
      terms_full <- if (use_eqspec) unique(c(base_exog, inst_excl_by_eq[[eqn]])) else fit_obj$inst_terms
      ivm <- AER::ivreg(f, instruments = reformulate(terms_full), data = data)
      mats <- .get_XZ_u(ivm); if (is.null(mats)) next
      Z <- mats$Z; u <- as.numeric(mats$u); n <- NROW(Z)
      Zu <- Z * u
      S  <- if (!is.null(cl)) { G <- rowsum(Zu, group = as.factor(.align_cluster(ivm, cl))); crossprod(as.matrix(G))/n } else crossprod(Zu)/n
      gbar <- colMeans(Zu); se <- sqrt(diag(S)/n); t <- gbar / se
      tab  <- data.frame(eq = eqn, inst = colnames(Z), gbar = gbar, se = se, t = t, abs_t = abs(t))
      if (isTRUE(only_excluded)) tab <- dplyr::filter(tab, !(inst %in% base_exog))
      out[[length(out)+1L]] <- tab
    }
    res <- dplyr::bind_rows(out)
    dplyr::filter(res, inst != "(Intercept)")
  }
})(iv_moment_diag_auto)

## ---------- (1) ROUND 2: relaxar “must_keep” e permitir df_J <= 2 ----------
sweep2 <- eqspec_sweep_all(
  fit_base         = fit_q,
  data             = df_iv,
  inst_excl_start  = excl_now,
  must_keep_global = character(0),   # << agora pode dropar iv_op* se ajudar
  target_dfJ_max   = 2,              # << permite mais folga por eq
  alpha_flag       = 0.05,
  cluster_var      = cl_var,
  max_add_size     = 5,
  max_combos       = 1000,
  verbose          = TRUE
)

cat("\n== J (ROUND 2, eq-specific) ==\n")
print(jtests_overid_auto(sweep2$fit, df_iv, cluster_var = cl_var))
cat("\nResumo J (ROUND 2):\n"); print(summarize_overid(jtests_overid_auto(sweep2$fit, df_iv, cluster_var = cl_var)))

## ---------- (2) Sensibilidade a cluster para o mesmo fit eq-specific ----------
compare_cluster_specs_auto <- function(fit_obj, data, cluster_vars = c("psu","uf",NA)){
  tabs <- lapply(cluster_vars, function(cv){
    jt <- jtests_overid_auto(fit_obj, data, cluster_var = if (is.na(cv)) NULL else cv)
    jt$cluster <- if (is.na(cv)) "none" else as.character(cv)
    jt
  })
  dplyr::bind_rows(tabs)
}
cat("\n== Sensibilidade de cluster (ROUND 2 fit) ==\n")
print(compare_cluster_specs_auto(sweep2$fit, df_iv, cluster_vars = c("psu","uf",NA)))

## ---------- (3) C-tests eq-specific por bloco (diferença de Hansen) ----------
make_blocks_explicit <- function(fit_obj){
  # usa todos os candidatos que apareceram no pool global
  pool <- setdiff(unique(fit_obj$inst_terms), get_base_exog(fit_obj))
  list(
    IV_d_uf_X = grep("^IV_d_uf_X", pool, value = TRUE),
    iv_op     = grep("^iv_op\\d+",  pool, value = TRUE)
  ) |> (\(x) x[vapply(x, length, 1L)>0])()
}

c_tests_eqspec_blocks <- function(fit_base, data, fit_eqspec, cluster_var = NULL, max_k_drop = 2){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_base$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  inst_excl_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  blocks <- make_blocks_explicit(fit_base)
  
  out <- list()
  for (eqname in names(eqs)) {
    f <- eqs[[eqname]]
    full_terms <- unique(c(base_exog, inst_excl_by_eq[[eqname]]))
    iv_full <- AER::ivreg(f, instruments = reformulate(full_terms), data = data)
    Jf <- .j_overid_manual(iv_full, cluster = cl); if (!is.finite(Jf["df"]) || Jf["df"]<=0) next
    
    for (bk in names(blocks)) {
      B <- intersect(blocks[[bk]], inst_excl_by_eq[[eqname]])
      if (!length(B)) next
      kmax <- min(max_k_drop, length(B))
      for (k in 1:kmax) {
        for (drop_set in utils::combn(B, k, simplify = FALSE)) {
          keep_excl <- setdiff(inst_excl_by_eq[[eqname]], drop_set)
          terms_r   <- unique(c(base_exog, keep_excl))
          iv_r <- AER::ivreg(f, instruments = reformulate(terms_r), data = data)
          Jr <- .j_overid_manual(iv_r, cluster = cl)
          C  <- as.numeric(Jr["J"] - Jf["J"]); pC <- stats::pchisq(C, df = k, lower.tail = FALSE)
          out[[length(out)+1L]] <- data.frame(
            eq = eqname, block = bk, drop_set = paste(drop_set, collapse = " + "),
            k_removed = k, J_full = unname(Jf["J"]), dfJ_full=unname(Jf["df"]),
            J_restr = unname(Jr["J"]), dfJ_restr=unname(Jr["df"]),
            C_stat = C, df_C = k, p_C = pC, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  dplyr::bind_rows(out)
}

apply_eqspec_prune <- function(fit_base, data, fit_eqspec, drops_by_eq){
  inst_excl_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  for (eqn in names(drops_by_eq)) {
    inst_excl_by_eq[[eqn]] <- setdiff(inst_excl_by_eq[[eqn]], drops_by_eq[[eqn]])
  }
  refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)
}

## roda C-tests (ROUND 2 fit) e sugere drops por eq onde p_C < 0.05
ctab <- c_tests_eqspec_blocks(fit_q, df_iv, sweep2$fit, cluster_var = cl_var, max_k_drop = 2)
cat("\n== C-tests eq-specific (p_C) ==\n");
ctab_ord <- dplyr::arrange(ctab, eq, p_C)
print(utils::head(ctab_ord, 30))


drops_auto <- by(ctab, ctab$eq, \(df) unique(unlist(strsplit(df$drop_set[df$p_C < 0.05], "\\s*\\+\\s*"))))
drops_auto <- lapply(drops_auto, \(v) v[nzchar(v)])
drops_auto <- drops_auto[ vapply(drops_auto, length, 1L) > 0 ]
if (length(drops_auto)) {
  cat("\nDrops sugeridos por eq (p_C<0.05):\n"); print(drops_auto)
  fit_round3 <- apply_eqspec_prune(fit_q, df_iv, sweep2$fit, drops_auto)
} else {
  message("\n[C-tests] Nenhum drop sugerido. Mantendo ROUND 2.")
  fit_round3 <- sweep2$fit
}

## ---------- (4) JUST-ID opcional nas eqs ainda ruins ----------
force_justID_eqs <- function(fit_base, data, fit_eqspec, target_eqs, must_keep = character(0), cluster_var = NULL){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  
  for (eqname in target_eqs) {
    cand <- setdiff(inst_by_eq[[eqname]], must_keep)
    best <- NULL; best_p <- -Inf
    # testa subsets s.t. df_J == 0 (just-ID)
    f <- formula(fit_base$eq[[ match(eqname, vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))) ]])
    for (k in 0:length(cand)) {
      for (S in utils::combn(cand, k, simplify = FALSE)) {
        excl <- unique(c(must_keep, S))
        terms <- unique(c(base_exog, excl))
        ivm <- try(AER::ivreg(f, instruments = reformulate(terms), data = data), silent = TRUE)
        if (inherits(ivm, "try-error")) next
        mats <- .get_XZ_u(ivm); if (is.null(mats)) next
        rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank
        if ((rZ - rX) != 0) next
        Jr <- .j_overid_manual(ivm, cluster = cl)
        p  <- unname(Jr["p"])  # no just-ID, J é inaplicável; guardamos p para consistência (NA)
        if (is.na(p)) p <- 1  # prioriza realização de just-ID
        if (p > best_p) { best_p <- p; best <- excl }
      }
    }
    if (!is.null(best)) inst_by_eq[[eqname]] <- best
  }
  refit_quaids_try_eqlist(fit_base, data, inst_by_eq)
}

## identifica eqs ainda com p(J)<0.05
jt_r3 <- jtests_overid_auto(fit_round3, df_iv, cluster_var = cl_var)
eq_bad_r3 <- jt_r3$eq[ jt_r3$p_hansen < 0.05 ]

## (opcional) força just-ID apenas nas piores
if (length(eq_bad_r3)) {
  fit_final <- force_justID_eqs(
    fit_base  = fit_q,
    data      = df_iv,
    fit_eqspec= fit_round3,
    target_eqs= eq_bad_r3,
    must_keep = character(0),
    cluster_var = cl_var
  )
} else {
  fit_final <- fit_round3
}

## ---------- (5) Relatório final ----------
cat("\n== RELATÓRIO FINAL (eq-specific, pós-ROUND 3/just-ID se aplicável) ==\n")
invisible(post_refit_report_auto(fit_final, df_iv, cluster_var = cl_var, R = 2000, level = 0.95))

cat("\n== Moment diagnostics (only excluded, eq-specific, FINAL) ==\n")
diag_eq_final <- iv_moment_diag_auto(fit_final, df_iv, cluster_var = cl_var, only_excluded = TRUE)
print(dplyr::arrange(diag_eq_final, dplyr::desc(abs_t)) %>% dplyr::group_by(eq) %>% dplyr::slice_head(n=8), n=36)
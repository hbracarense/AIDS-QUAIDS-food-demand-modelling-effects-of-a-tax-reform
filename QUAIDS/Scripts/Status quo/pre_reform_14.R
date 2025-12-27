## =========================
## TOP TIER – FUNÇÕES CORE
## =========================
library(openxlsx)
get_base_exog <- function(fit_obj){
  rhs_all <- unique(fit_obj$rhs_terms)
  setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
}

pick_cluster <- function(data, candidates = c("psu","uf"), fallback = NULL) {
  nm <- names(data); hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else fallback
}

## F-HC3/cluster robusto (singularidade tolerante)
fs_block_F <- function(fit_obj, data, cluster_var = NULL) {
  data <- as.data.frame(data)
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  stopifnot(length(iv_excl) > 0)
  
  safe_F_one <- function(y) {
    f <- reformulate(c(controls, iv_excl), response = y)
    m <- stats::lm(f, data = data)
    X <- stats::model.matrix(m); u <- stats::residuals(m); n <- NROW(X)
    XtX <- crossprod(X); mu <- mean(diag(XtX)); if (!is.finite(mu) || mu <= 0) mu <- 1
    XtX <- XtX + diag(1e-10 * mu, ncol(XtX))
    XtXi <- MASS::ginv(XtX)
    Xu <- X * as.numeric(u)
    meat <- if (!is.null(cluster_var) && cluster_var %in% names(data)) {
      g <- rowsum(Xu, group = as.factor(data[[cluster_var]])); crossprod(as.matrix(g))
    } else crossprod(Xu)
    V <- XtXi %*% meat %*% XtXi
    b <- stats::coef(m); cn <- names(b)
    iv_in <- intersect(iv_excl, cn)
    if (!length(iv_in)) return(data.frame(endog=y,F=NA,df1=NA,df2=NA,p=NA))
    R <- matrix(0, nrow=length(iv_in), ncol=length(cn), dimnames=list(iv_in,cn))
    R[cbind(iv_in, iv_in)] <- 1
    rb <- as.matrix(R %*% b)
    RVRT <- R %*% V %*% t(R); RVRT <- (RVRT + t(RVRT))/2
    mu2 <- mean(diag(RVRT)); if (!is.finite(mu2) || mu2 <= 0) mu2 <- 1
    RVRT <- RVRT + diag(1e-10 * mu2, nrow(RVRT))
    RVRTi <- MASS::ginv(RVRT)
    Wald <- drop(t(rb) %*% RVRTi %*% rb); r <- length(iv_in); Fst <- Wald/r
    rX <- qr(X)$rank; df2 <- n - rX
    if (!is.finite(df2) || df2 <= 0) {
      p <- stats::pchisq(Wald, df=r, lower.tail=FALSE)
      data.frame(endog=y, F=Fst, df1=r, df2=NA, p=p)
    } else {
      p <- stats::pf(Fst, r, df2, lower.tail=FALSE)
      data.frame(endog=y, F=Fst, df1=r, df2=df2, p=p)
    }
  }
  do.call(rbind, lapply(endo_ln, safe_F_one))
}

## ID por-equação (garante dfJ >= 0)
ensure_identification_eqlist <- function(fit_obj, data, inst_excl_by_eq, extra_pool = NULL){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  if (is.null(extra_pool)) extra_pool <- setdiff(unique(unlist(inst_excl_by_eq)), base_exog)
  added <- list()
  for (eqname in names(eqs)) {
    X <- .mm(delete.response(terms(eqs[[eqname]])), data, add_intercept = TRUE)
    Z <- .mm(reformulate(unique(c(base_exog, inst_excl_by_eq[[eqname]])), intercept = TRUE), data, add_intercept = TRUE)
    dfJ <- ncol(Z) - qr(X)$rank
    if (!is.finite(dfJ) || dfJ < 0) {
      cand <- setdiff(extra_pool, inst_excl_by_eq[[eqname]])
      if (!length(cand)) stop("Sem IVs para consertar identificação em ", eqname)
      k_add <- 0L
      for (z in cand) {
        inst_excl_by_eq[[eqname]] <- unique(c(inst_excl_by_eq[[eqname]], z))
        Z <- .mm(reformulate(unique(c(base_exog, inst_excl_by_eq[[eqname]])), intercept = TRUE), data, add_intercept = TRUE)
        dfJ <- ncol(Z) - qr(X)$rank; k_add <- k_add + 1L
        if (dfJ >= 0) break
      }
      added[[eqname]] <- tail(inst_excl_by_eq[[eqname]], k_add)
      if (dfJ < 0) stop("Ainda under-ID em ", eqname)
    }
  }
  attr(inst_excl_by_eq, "added_to_fix_id") <- added
  inst_excl_by_eq
}

## Refit 3SLS eq-specific (labels seguros)
refit_quaids_try_eqlist <- function(fit_obj, data, inst_excl_by_eq) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  lhs <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  inst_excl_by_eq <- ensure_identification_eqlist(
    fit_obj, data, inst_excl_by_eq,
    extra_pool = setdiff(unique(unlist(inst_excl_by_eq)), base_exog)
  )
  inst_list <- lapply(seq_along(eqs), function(i){
    terms <- unique(c(base_exog, inst_excl_by_eq[[lhs[i]]]))
    stats::reformulate(terms, intercept = FALSE)
  })
  eq_labels <- paste0("eq", seq_along(eqs))
  names(eqs) <- eq_labels; names(inst_list) <- eq_labels
  fit_new <- try(systemfit(eqs, data=data, method="3SLS", inst=inst_list), silent=TRUE)
  if (inherits(fit_new, "try-error")) stop("systemfit(3SLS) falhou: ", attr(fit_new,"condition")$message)
  
  shareNames <- fit_obj$shareNames; yhat_mat <- sapply(fit_new$eq, stats::fitted)
  if (is.list(yhat_mat)) yhat_mat <- do.call(cbind, yhat_mat)
  fitted_included_df <- as.data.frame(yhat_mat, optional = TRUE)
  colnames(fitted_included_df) <- lhs
  fitted_full <- fitted_included_df
  omit_share  <- setdiff(shareNames, lhs)
  if (length(omit_share) == 1L) fitted_full[[omit_share]] <- 1 - rowSums(fitted_included_df)
  fitted_full <- fitted_full[, shareNames, drop = FALSE]
  observed_full <- as.data.frame(data[, shareNames, drop = FALSE])
  
  res <- list(
    method="3SLS", priceIndex=fit_obj$priceIndex, omit_share=fit_obj$omit_share,
    drop_price=fit_obj$drop_price, priceNames=fit_obj$priceNames,
    shareNames=fit_obj$shareNames, rhs_terms=fit_obj$rhs_terms,
    inst_terms=unique(c(base_exog, unlist(inst_excl_by_eq))),
    fit=fit_new, eq=fit_new$eq,
    fitted_shares=fitted_full, observed_shares=observed_full
  )
  class(res) <- "quaids_km1_fit"
  attr(res, "inst_excl_by_eq") <- inst_excl_by_eq
  res
}

## Hansen J autodetectando eq-specific
jtests_overid_auto <- function(fit_obj, data, cluster_var=NULL, verbose=FALSE){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_excl_by_eq <- attr(fit_obj, "inst_excl_by_eq")
  use_eqspec <- is.list(inst_excl_by_eq)
  rows <- lapply(names(eqs), function(eqname){
    f <- eqs[[eqname]]
    inst_terms <- if (use_eqspec) unique(c(base_exog, inst_excl_by_eq[[eqname]])) else fit_obj$inst_terms
    ivm <- try(AER::ivreg(f, instruments=reformulate(inst_terms), data=data), silent=TRUE)
    if (inherits(ivm, "try-error")) {
      if (verbose) message("ivreg falhou em ", eqname, ": ", attr(ivm,"condition")$message)
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA, J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats))
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA, J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    s_h <- try(AER::sargan(ivm), silent=TRUE)
    J_s <- if (!inherits(s_h,"try-error")) unname(as.numeric(s_h$statistic)) else NA_real_
    p_s <- if (!inherits(s_h,"try-error")) unname(as.numeric(s_h$p.value))   else NA_real_
    Jr  <- .j_overid_manual(ivm, cluster=cl)
    data.frame(eq=eqname, rank_X=rX, rank_Z=rZ, df_J=dfJ,
               J_sargan=J_s, p_sargan=p_s, J_hansen=unname(Jr["J"]), p_hansen=unname(Jr["p"]))
  })
  dplyr::bind_rows(rows)
}

summarize_overid <- function(jtab){
  data.frame(min_dfJ = suppressWarnings(min(jtab$df_J, na.rm=TRUE)),
             min_p   = suppressWarnings(min(jtab$p_hansen, na.rm=TRUE)))
}

## Moment diagnostics (evita warnings com only_excluded=TRUE)
iv_moment_diag_auto <- function(fit_obj, data, cluster_var=NULL, only_excluded=FALSE){
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
    ivm <- AER::ivreg(f, instruments=reformulate(terms_full), data=data)
    mats <- .get_XZ_u(ivm); if (is.null(mats)) next
    Z <- mats$Z; u <- as.numeric(mats$u); n <- NROW(Z); Zu <- Z*u
    S <- if (!is.null(cl)) { G <- rowsum(Zu, group=as.factor(.align_cluster(ivm, cl))); crossprod(as.matrix(G))/n } else crossprod(Zu)/n
    gbar <- colMeans(Zu); se <- sqrt(diag(S)/n); t <- gbar/se
    tab <- data.frame(eq=eqn, inst=colnames(Z), gbar=gbar, se=se, t=t, abs_t=abs(t))
    if (isTRUE(only_excluded)) tab <- dplyr::filter(tab, !(inst %in% base_exog))
    out[[length(out)+1L]] <- tab
  }
  res <- dplyr::bind_rows(out)
  dplyr::filter(res, inst != "(Intercept)")
}

## Relatório auto (respeita eq-specific)
post_refit_report_auto <- function(fit_new, data, cluster_var=NULL, R=2000, level=0.95, seed=123){
  cat("\n== Over-ID (Hansen J) por equação ==\n")
  jt <- jtests_overid_auto(fit_new, data, cluster_var=cluster_var, verbose=FALSE); print(jt)
  if (any(jt$df_J <= 0, na.rm=TRUE)) message("Aviso: há equações just/under-ID (df_J <= 0). Hansen J não é aplicável nessas eqs.")
  cat("\n== Força dos IVs (F-HC3/CL) ==\n"); print(fs_block_F(fit_new, data, cluster_var=cluster_var))
  cat("\n== partial R^2 dos IVs excluídos ==\n"); print(partial_R2_exclIV(fit_new, data))
  cat("\n== Wu-Hausman (control function) por equação ==\n"); print(wu_hausman_cf(fit_new, data, cluster_var=cluster_var))
  if (exists("quaids_elasticities_ci", mode="function")) {
    cat("\n== Elasticidades com IC (ponto: median/observed) ==\n")
    res_ci <- try(quaids_elasticities_ci(fit_new, at="medians", w_source="observed", R=R, level=level, seed=seed), silent=TRUE)
    if (!inherits(res_ci,"try-error") && is.list(res_ci) && !is.null(res_ci$expenditure)) print(res_ci$expenditure)
  } else message("\n(Elasticidades/IC puladas: 'quaids_elasticities_ci' não encontrada.)")
  invisible(jt)
}

## Busca local por eq foco (df_J <= alvo)
tune_focus_eq_subset <- function(fit_obj, data, current_excl, focus_eq,
                                 must_keep = NULL, choose_from = NULL,
                                 target_dfJ_max = 1, cluster_var = NULL,
                                 max_add_size = 5, max_combos = 600L) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  if (is.null(choose_from)) choose_from <- setdiff(current_excl, must_keep)
  if (is.null(must_keep))    must_keep  <- character(0)
  best <- list(p=-Inf, set=NULL, dfJ=NA); sizes <- 0:min(length(choose_from), max_add_size); tested <- 0L
  fmla_eq <- formula(fit_obj$eq[[focus_eq]])
  for (s in sizes) {
    if (tested > max_combos) break
    cmbs <- utils::combn(choose_from, s, simplify=FALSE)
    for (S in cmbs) {
      tested <- tested + 1L
      excl_focus <- unique(c(must_keep, S))
      inst_terms_focus <- unique(c(base_exog, excl_focus))
      ivm <- try(AER::ivreg(fmla_eq, instruments=reformulate(inst_terms_focus), data=data), silent=TRUE)
      if (inherits(ivm, "try-error")) next
      mats <- .get_XZ_u(ivm); if (is.null(mats)) next
      rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
      if (!is.finite(dfJ) || dfJ < 0 || dfJ > target_dfJ_max) next
      Jr <- .j_overid_manual(ivm, cluster = if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL)
      p  <- unname(Jr["p"])
      if (is.finite(p) && p > best$p) best <- list(p=p, set=excl_focus, dfJ=dfJ)
    }
  }
  if (is.null(best$set)) return(NULL); best
}

## Varredura eq-specific para eqs com p(J) baixo
eqspec_sweep_all <- function(fit_base, data, inst_excl_start,
                             must_keep_global=character(0),
                             target_dfJ_max=1, alpha_flag=0.05,
                             cluster_var=NULL, max_add_size=5, max_combos=600L, verbose=TRUE){
  data <- as.data.frame(data)
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) cluster_var else NULL
  inst_excl_by_eq <- setNames(replicate(length(eq_names), inst_excl_start, simplify=FALSE), eq_names)
  tmp_fit <- refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)
  jt <- jtests_overid_auto(tmp_fit, data, cluster_var=cl, verbose=FALSE)
  bad_eqs <- jt$eq[ jt$p_hansen < alpha_flag ]
  if (verbose) message("[eqspec_sweep_all] Equações com p(J) < ", alpha_flag, ": ",
                       if (length(bad_eqs)) paste(bad_eqs, collapse=", ") else "(nenhuma)")
  if (!length(bad_eqs)) return(list(inst_excl_by_eq=inst_excl_by_eq, fit=tmp_fit, jt=jt, improved=FALSE))
  for (eqname in bad_eqs) {
    choose_from <- setdiff(inst_excl_by_eq[[eqname]], must_keep_global)
    best <- tune_focus_eq_subset(fit_obj=fit_base, data=data,
                                 current_excl=inst_excl_by_eq[[eqname]],
                                 focus_eq=match(eqname, eq_names),
                                 must_keep=must_keep_global,
                                 choose_from=choose_from,
                                 target_dfJ_max=target_dfJ_max,
                                 cluster_var=cl, max_add_size=max_add_size, max_combos=max_combos)
    if (!is.null(best)) {
      inst_excl_by_eq[[eqname]] <- best$set
      if (verbose) message("  [",eqname,"] novo subset (df_J=", best$dfJ, ") com p(J)=", signif(best$p,3))
    } else if (verbose) message("  [",eqname,"] sem subset viável com df_J <=", target_dfJ_max, ". Mantendo.")
  }
  fit_eqspec <- refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)
  jt2 <- jtests_overid_auto(fit_eqspec, data, cluster_var=cl)
  list(inst_excl_by_eq=inst_excl_by_eq, fit=fit_eqspec, jt=jt2, improved=TRUE)
}

## C-tests por bloco (diferença de Hansen)
make_blocks_explicit <- function(fit_obj){
  pool <- setdiff(unique(fit_obj$inst_terms), get_base_exog(fit_obj))
  lst <- list(IV_d_uf_X = grep("^IV_d_uf_X", pool, value=TRUE),
              iv_op     = grep("^iv_op\\d+",  pool, value=TRUE))
  lst[vapply(lst, length, 1L) > 0]
}

c_tests_eqspec_blocks <- function(fit_base, data, fit_eqspec, cluster_var=NULL, max_k_drop=2){
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
    iv_full <- AER::ivreg(f, instruments=reformulate(full_terms), data=data)
    Jf <- .j_overid_manual(iv_full, cluster=cl)
    if (!is.finite(Jf["df"]) || Jf["df"] <= 0) next
    for (bk in names(blocks)) {
      B <- intersect(blocks[[bk]], inst_excl_by_eq[[eqname]])
      if (!length(B)) next
      kmax <- min(max_k_drop, length(B))
      for (k in 1:kmax) for (drop_set in utils::combn(B, k, simplify=FALSE)) {
        keep_excl <- setdiff(inst_excl_by_eq[[eqname]], drop_set)
        terms_r   <- unique(c(base_exog, keep_excl))
        iv_r <- AER::ivreg(f, instruments=reformulate(terms_r), data=data)
        Jr <- .j_overid_manual(iv_r, cluster=cl)
        if (!is.finite(Jr["df"]) || Jr["df"] <= 0) next
        C  <- max(0, as.numeric(Jr["J"] - Jf["J"]))
        pC <- stats::pchisq(C, df=k, lower.tail=FALSE)
        out[[length(out)+1L]] <- data.frame(
          eq=eqname, block=bk, drop_set=paste(drop_set, collapse=" + "),
          k_removed=k, J_full=unname(Jf["J"]), dfJ_full=unname(Jf["df"]),
          J_restr=unname(Jr["J"]), dfJ_restr=unname(Jr["df"]),
          C_stat=C, df_C=k, p_C=pC, stringsAsFactors=FALSE
        )
      }
    }
  }
  dplyr::bind_rows(out)
}

apply_eqspec_prune <- function(fit_base, data, fit_eqspec, drops_by_eq){
  inst_excl_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  for (eqn in names(drops_by_eq)) inst_excl_by_eq[[eqn]] <- setdiff(inst_excl_by_eq[[eqn]], drops_by_eq[[eqn]])
  refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)
}

force_justID_eqs <- function(fit_base, data, fit_eqspec, target_eqs, must_keep=character(0), cluster_var=NULL){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  for (eqname in target_eqs) {
    cand <- setdiff(inst_by_eq[[eqname]], must_keep)
    best <- NULL; best_p <- -Inf
    f <- formula(fit_base$eq[[ match(eqname, eq_names) ]])
    for (k in 0:length(cand)) for (S in utils::combn(cand, k, simplify=FALSE)) {
      excl <- unique(c(must_keep, S)); terms <- unique(c(base_exog, excl))
      ivm <- try(AER::ivreg(f, instruments=reformulate(terms), data=data), silent=TRUE)
      if (inherits(ivm, "try-error")) next
      mats <- .get_XZ_u(ivm); if (is.null(mats)) next
      rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank
      if ((rZ - rX) != 0) next
      Jr <- .j_overid_manual(ivm, cluster=cl); p <- unname(Jr["p"]); if (is.na(p)) p <- 1
      if (p > best_p) { best_p <- p; best <- excl }
    }
    if (!is.null(best)) inst_by_eq[[eqname]] <- best
  }
  refit_quaids_try_eqlist(fit_base, data, inst_by_eq)
}

## =========================
## DRIVER ÚNICO (plug & play)
## =========================
run_top_tier <- function(fit_q, df_iv,
                         excl_now,
                         must_keep_global = c("iv_op01","iv_op02","iv_op03"),
                         cluster_preference = c("psu","uf",NA),
                         alpha_flag = 0.05,
                         sweep_dfJ_max = 1,
                         do_round2_relax = TRUE,
                         round2_dfJ_max = 2,
                         do_ctests = TRUE,
                         enforce_just_id = c("never","if_reject","always")) {
  
  enforce_just_id <- match.arg(enforce_just_id)
  cl_var <- {prefs <- cluster_preference; prefs <- prefs[!is.na(prefs)]
  pick_cluster(df_iv, candidates = c(prefs, "psu","uf"), fallback = NULL)}
  
  cat("\n[RUN] Cluster escolhido para diagnósticos: ", if (is.null(cl_var)) "none" else cl_var, "\n", sep="")
  
  ## ROUND 1: varredura eq-specific com dfJ <= sweep_dfJ_max
  sweep1 <- eqspec_sweep_all(fit_q, df_iv, inst_excl_start=excl_now,
                             must_keep_global=must_keep_global,
                             target_dfJ_max=sweep_dfJ_max, alpha_flag=alpha_flag,
                             cluster_var=cl_var, max_add_size=5, max_combos=600, verbose=TRUE)
  
  fit_curr <- sweep1$fit
  
  ## ROUND 2 (opcional): relaxa must_keep e dfJ max
  if (isTRUE(do_round2_relax)) {
    sweep2 <- eqspec_sweep_all(fit_q, df_iv, inst_excl_start=excl_now,
                               must_keep_global=character(0),
                               target_dfJ_max=round2_dfJ_max, alpha_flag=alpha_flag,
                               cluster_var=cl_var, max_add_size=5, max_combos=1000, verbose=TRUE)
    fit_curr <- sweep2$fit
  }
  
  ## C-tests por bloco (sugestões de drop)
  drops_auto <- NULL; ctab <- NULL
  if (isTRUE(do_ctests)) {
    ctab <- c_tests_eqspec_blocks(fit_q, df_iv, fit_curr, cluster_var=cl_var, max_k_drop=2)
    if (!is.null(ctab) && nrow(ctab)) {
      ctab_ord <- dplyr::arrange(ctab, eq, p_C)
      cat("\n== C-tests eq-specific (top 30 por p_C) ==\n"); print(utils::head(ctab_ord, 30))
      tmp <- by(ctab_ord, ctab_ord$eq, \(df) unique(unlist(strsplit(df$drop_set[df$p_C < alpha_flag], "\\s*\\+\\s*"))))
      tmp <- lapply(tmp, \(v) v[nzchar(v)]); tmp <- tmp[vapply(tmp, length, 1L) > 0]
      if (length(tmp)) {
        drops_auto <- tmp
        cat("\n[Drops sugeridos por eq - p_C <", alpha_flag, "]\n", sep="")
        print(drops_auto)
        fit_curr <- apply_eqspec_prune(fit_q, df_iv, fit_curr, drops_auto)
      }
    }
  }
  
  ## Just-ID seletivo/sempre/never
  if (enforce_just_id != "never") {
    jt_now <- jtests_overid_auto(fit_curr, df_iv, cluster_var=cl_var)
    eq_bad <- if (enforce_just_id == "if_reject") jt_now$eq[jt_now$p_hansen < alpha_flag] else
      vapply(fit_q$eq, function(m) as.character(formula(m)[[2]]), character(1))
    if (length(eq_bad)) {
      fit_curr <- force_justID_eqs(fit_q, df_iv, fit_curr, target_eqs=eq_bad, must_keep=character(0), cluster_var=cl_var)
    }
  }
  
  ## Relatório final
  cat("\n== RELATÓRIO FINAL (TOP TIER) ==\n")
  jt_final <- post_refit_report_auto(fit_curr, df_iv, cluster_var=cl_var, R=2000, level=0.95)
  
  invisible(list(
    fit_final = fit_curr,
    tables = list(
      jt_final = jt_final,
      ctests   = ctab,
      F_block  = fs_block_F(fit_curr, df_iv, cluster_var=cl_var),
      pR2      = partial_R2_exclIV(fit_curr, df_iv),
      wu       = wu_hausman_cf(fit_curr, df_iv, cluster_var=cl_var),
      moments_excl = iv_moment_diag_auto(fit_curr, df_iv, cluster_var=cl_var, only_excluded=TRUE)
    ),
    cluster_used = cl_var
  ))
}

excl_now <- c("IV_d_uf_X11","IV_d_uf_X12","IV_d_uf_X13","IV_d_uf_X15","iv_op01","iv_op02","iv_op03")

out <- run_top_tier(
  fit_q, df_iv, excl_now,
  must_keep_global = c("iv_op01","iv_op02","iv_op03"),
  cluster_preference = c("psu","uf", NA),   # >> use "psu" por padrão
  alpha_flag = 0.05,
  sweep_dfJ_max = 1,
  do_round2_relax = TRUE,
  round2_dfJ_max = 2,
  do_ctests = TRUE,
  enforce_just_id = "if_reject"            # "never" | "if_reject" | "always"
)

out2 <- run_top_tier(
  fit_q, df_iv, excl_now,
  must_keep_global = c("iv_op01","iv_op02","iv_op03"),
  cluster_preference = "uf",
  alpha_flag = 0.05,
  sweep_dfJ_max = 1,
  do_round2_relax = TRUE,
  round2_dfJ_max = 2,
  do_ctests = TRUE,
  enforce_just_id = "if_reject"   # << não "always"
)

## 1) Reconstroi o mapeamento eq->IVs do ROUND 1 (df_J=1) e aplica só os drops com p_C<0.05
s1 <- eqspec_sweep_all(
  fit_base         = fit_q,
  data             = df_iv,
  inst_excl_start  = excl_now,
  must_keep_global = c("iv_op01","iv_op02","iv_op03"),
  target_dfJ_max   = 1,
  alpha_flag       = 0.05,
  cluster_var      = "uf",
  max_add_size     = 5,
  max_combos       = 600,
  verbose          = FALSE
)

inst_map <- s1$inst_excl_by_eq
inst_map$w_despesahat2 <- setdiff(inst_map$w_despesahat2, "IV_d_uf_X15")
inst_map$w_despesahat6 <- setdiff(inst_map$w_despesahat6, "IV_d_uf_X11")

fit_rerun <- refit_quaids_try_eqlist(fit_q, df_iv, inst_map)

## 2) Só força just-ID onde AINDA rejeitar sob cluster=uf
jt2 <- jtests_overid_auto(fit_rerun, df_iv, cluster_var = "uf")
eq_bad <- jt2$eq[ is.finite(jt2$p_hansen) & jt2$p_hansen < 0.05 ]

fit_fix <- if (length(eq_bad)) {
  force_justID_eqs(
    fit_base   = fit_q,
    data       = df_iv,
    fit_eqspec = fit_rerun,
    target_eqs = eq_bad,
    must_keep  = character(0),
    cluster_var = "uf"
  )
} else fit_rerun

## 3) Relatório final — agora mantendo over-ID onde possível
post_refit_report_auto(fit_fix, df_iv, cluster_var = "uf", R = 2000, level = 0.95)

## === Promote over-ID para eqs que ficaram just-ID (cluster = uf) ===

bump_dfJ_to_one <- function(fit_base, data, fit_eqspec, eqn, cluster_var = "uf",
                            avoid = character(0)) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  stopifnot(is.list(inst_by_eq), eqn %in% names(inst_by_eq))
  current <- unique(inst_by_eq[[eqn]])
  
  # pool global de IVs excluídos possíveis
  all_pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  cand <- setdiff(all_pool, union(current, avoid))
  
  # fórmula da equação alvo
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  f <- formula(fit_base$eq[[ match(eqn, eq_names) ]])
  
  best <- NULL; bestp <- -Inf
  for (z in cand) {
    terms <- unique(c(base_exog, current, z))
    ivm <- try(AER::ivreg(f, instruments = reformulate(terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) next
    mats <- .get_XZ_u(ivm); if (is.null(mats)) next
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank
    if ((rZ - rX) != 1) next  # queremos df_J exatamente 1
    
    Jr <- .j_overid_manual(ivm, cluster = if (cluster_var %in% names(data)) data[[cluster_var]] else NULL)
    p  <- unname(Jr["p"])
    if (is.finite(p) && p > bestp) { bestp <- p; best <- z }
  }
  
  if (!is.null(best)) {
    inst_by_eq[[eqn]] <- unique(c(current, best))
    message("[", eqn, "] promovida a df_J=1 com IV adicionado: ", best, " | p(J)≈", signif(bestp,3))
  } else {
    message("[", eqn, "] não encontrei candidato que dê df_J=1 sem quebrar o modelo.")
  }
  refit_quaids_try_eqlist(fit_base, data, inst_by_eq)
}

## 1) Parta do fit atual (quase todo just-ID)
fit_now <- fit_fix  # ex.: o 'fit_eqspec' ou 'fit_final' da sua última corrida
  
  ## 2) Promova over-ID nas eqs onde IV é claramente necessário (pela Wu): 2 e 6
  fit_bumped <- bump_dfJ_to_one(
    fit_base  = fit_q,
    data      = df_iv,
    fit_eqspec= fit_now,
    eqn       = "w_despesahat2",
    cluster_var = "uf",
    avoid     = c("IV_d_uf_X15") # não recolocar o que o C-test condenou
  )
fit_bumped <- bump_dfJ_to_one(
  fit_base  = fit_q,
  data      = df_iv,
  fit_eqspec= fit_bumped,
  eqn       = "w_despesahat6",
  cluster_var = "uf",
  avoid     = c("IV_d_uf_X11")
)

## 3) Relatório
post_refit_report_auto(fit_bumped, df_iv, cluster_var = "uf", R = 2000, level = 0.95)

c_tests_eqspec_blocks_fix <- function(fit_base, data, fit_eqspec, cluster_var=NULL, max_k_drop=2){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_base$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  
  # blocos explícitos (ajuste aos seus nomes):
  pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  blocks <- list(
    IV_d_uf_X = grep("^IV_d_uf_X", pool, value=TRUE),
    iv_op     = grep("^iv_op\\d+",  pool, value=TRUE)
  )
  blocks <- blocks[vapply(blocks, length, 1L)>0]
  
  out <- list()
  for (eqn in names(eqs)) {
    f <- eqs[[eqn]]
    full_terms <- unique(c(base_exog, inst_by_eq[[eqn]]))
    iv_full <- AER::ivreg(f, instruments = reformulate(full_terms), data = data)
    Jf <- .j_overid_manual(iv_full, cluster = cl); if (!is.finite(Jf["df"]) || Jf["df"]<=0) next
    
    for (bk in names(blocks)) {
      B <- intersect(blocks[[bk]], inst_by_eq[[eqn]])
      if (!length(B)) next
      kmax <- min(max_k_drop, length(B))
      for (k in 1:kmax) for (drop_set in utils::combn(B, k, simplify=FALSE)) {
        keep_excl <- setdiff(inst_by_eq[[eqn]], drop_set)
        iv_r <- AER::ivreg(f, instruments = reformulate(unique(c(base_exog, keep_excl))), data = data)
        Jr <- .j_overid_manual(iv_r, cluster = cl)
        C  <- as.numeric(Jf["J"] - Jr["J"])              # <<< diferença corrigida
        pC <- stats::pchisq(C, df = k, lower.tail = FALSE)
        out[[length(out)+1L]] <- data.frame(
          eq=eqn, block=bk, drop_set=paste(drop_set, collapse=" + "),
          k_removed=k, J_full=unname(Jf["J"]), dfJ_full=unname(Jf["df"]),
          J_restr=unname(Jr["J"]), dfJ_restr=unname(Jr["df"]),
          C_stat=C, df_C=k, p_C=pC, stringsAsFactors=FALSE
        )
      }
    }
  }
  dplyr::bind_rows(out)
}

swap_to_dfJ_one <- function(fit_base, data, fit_eqspec, eqn, cluster_var="uf",
                            avoid = character(0), max_drop_each = 2){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq"); stopifnot(is.list(inst_by_eq))
  cur <- unique(inst_by_eq[[eqn]])
  
  all_pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  add_cand <- setdiff(all_pool, union(cur, avoid))
  
  f <- formula(fit_base$eq[[ match(eqn, vapply(fit_base$eq, \(m) as.character(formula(m)[[2]]), character(1))) ]])
  
  best <- list(p=-Inf, add=NA, drop=character(0))
  # tenta: (i) só add (se df_J==0), (ii) drop 1 + add 1 mantendo df_J==1
  comb_drop <- c(list(character(0)), lapply(cur, \(z) z))  # drop 0 ou 1
  for (D in comb_drop) for (A in add_cand) {
    terms <- unique(c(base_exog, setdiff(cur, D), A))
    ivm <- try(AER::ivreg(f, instruments=reformulate(terms), data=data), silent=TRUE)
    if (inherits(ivm,"try-error")) next
    mats <- .get_XZ_u(ivm); if (is.null(mats)) next
    if ((qr(mats$Z)$rank - qr(mats$X)$rank) != 1) next
    Jr <- .j_overid_manual(ivm, cluster = if (cluster_var %in% names(data)) data[[cluster_var]] else NULL)
    p  <- unname(Jr["p"])
    if (is.finite(p) && p > best$p) best <- list(p=p, add=A, drop=D)
  }
  
  if (is.finite(best$p) && best$p > -Inf) {
    inst_by_eq[[eqn]] <- unique(c(setdiff(cur, best$drop), best$add))
    message("[swap_to_dfJ_one:", eqn, "] +", best$add,
            if (length(best$drop)) paste0(" −", paste(best$drop, collapse=",")) else "",
            " | p(J)≈", signif(best$p,3))
    return(refit_quaids_try_eqlist(fit_base, data, inst_by_eq))
  } else {
    message("[swap_to_dfJ_one:", eqn, "] não encontrou swap viável.")
    return(fit_eqspec)
  }
}

fit_sw2 <- swap_to_dfJ_one(fit_q, df_iv, fit_bumped, "w_despesahat2", cluster_var="uf",
                           avoid=c("IV_d_uf_X15"))
fit_sw2 <- swap_to_dfJ_one(fit_q, df_iv, fit_sw2,      "w_despesahat6", cluster_var="uf",
                           avoid=c("IV_d_uf_X11"))
post_refit_report_auto(fit_sw2, df_iv, cluster_var="uf", R=2000, level=0.95)
print(compare_cluster_specs_auto(fit_sw2, df_iv, cluster_vars = c("psu","uf",NA)))

c_tests_eqspec_blocks_fix <- function(fit_base, data, fit_eqspec, cluster_var="uf", max_k_drop=2){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_base$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  
  pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  blocks <- list(IV_d_uf_X = grep("^IV_d_uf_X", pool, value=TRUE),
                 iv_op     = grep("^iv_op\\d+",  pool, value=TRUE))
  blocks <- blocks[vapply(blocks, length, 1L)>0]
  
  out <- list(); eqn <- "w_despesahat6"   # <<< foco
  f <- eqs[[eqn]]
  full_terms <- unique(c(base_exog, inst_by_eq[[eqn]]))
  iv_full <- AER::ivreg(f, instruments=reformulate(full_terms), data=data)
  Jf <- .j_overid_manual(iv_full, cluster=cl); if (!is.finite(Jf["df"]) || Jf["df"]<=0) stop("eq6 não over-ID")
  
  for (bk in names(blocks)) {
    B <- intersect(blocks[[bk]], inst_by_eq[[eqn]]); if (!length(B)) next
    kmax <- min(max_k_drop, length(B))
    for (k in 1:kmax) for (drop_set in utils::combn(B, k, simplify=FALSE)) {
      keep <- setdiff(inst_by_eq[[eqn]], drop_set)
      iv_r <- AER::ivreg(f, instruments=reformulate(unique(c(base_exog, keep))), data=data)
      Jr <- .j_overid_manual(iv_r, cluster=cl)
      C  <- as.numeric(Jf["J"] - Jr["J"])             # <<< diferença correta
      pC <- stats::pchisq(C, df=k, lower.tail=FALSE)
      out[[length(out)+1L]] <- data.frame(eq=eqn, block=bk,
                                          drop_set=paste(drop_set, collapse=" + "), k_removed=k,
                                          J_full=unname(Jf["J"]), J_restr=unname(Jr["J"]),
                                          C_stat=C, df_C=k, p_C=pC)
    }
  }
  dplyr::arrange(dplyr::bind_rows(out), p_C)
}

ct6 <- c_tests_eqspec_blocks_fix(fit_q, df_iv, fit_sw2, cluster_var="uf", max_k_drop=2)
print(head(ct6, 20))

fit_try <- swap_to_dfJ_one(
  fit_q, df_iv, fit_sw2, "w_despesahat6", cluster_var="uf",
  avoid=c("IV_d_uf_X11"),      # mantemos fora (reprovado antes)
  max_drop_each = 1            # drop 0 ou 1 (o helper já faz isso)
)
post_refit_report_auto(fit_try, df_iv, cluster_var="uf", R=2000, level=0.95)
fit_try2 <- swap_to_dfJ_one(
  fit_q, df_iv, fit_sw2, "w_despesahat6", cluster_var="uf",
  avoid=c("IV_d_uf_X11","IV_d_uf_X14")
)
post_refit_report_auto(fit_try2, df_iv, cluster_var="uf", R=2000, level=0.95)

c_tests_eq6_fix <- function(fit_base, data, fit_eqspec, cluster_var="uf", max_k_drop=2){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  cl <- if (cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  # fórmulas e IVs da eq6
  eqs <- lapply(fit_base$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  eqn <- "w_despesahat6"
  f <- eqs[[eqn]]
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  full_terms <- unique(c(base_exog, inst_by_eq[[eqn]]))
  
  iv_full <- AER::ivreg(f, instruments=reformulate(full_terms), data=data)
  Jf <- .j_overid_manual(iv_full, cluster=cl); if (!is.finite(Jf["df"]) || Jf["df"]<=0) stop("eq6 não over-ID")
  
  # blocos
  pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  blocks <- list(
    IV_d_uf_X = grep("^IV_d_uf_X", pool, value=TRUE),
    iv_op     = grep("^iv_op\\d+",  pool, value=TRUE)
  )
  blocks <- blocks[vapply(blocks, length, 1L)>0]
  
  out <- list()
  for (bk in names(blocks)) {
    B <- intersect(blocks[[bk]], inst_by_eq[[eqn]]); if (!length(B)) next
    kmax <- min(max_k_drop, length(B))
    for (k in 1:kmax) for (drop_set in utils::combn(B, k, simplify=FALSE)) {
      keep <- setdiff(inst_by_eq[[eqn]], drop_set)
      iv_r <- AER::ivreg(f, instruments=reformulate(unique(c(base_exog, keep))), data=data)
      Jr <- .j_overid_manual(iv_r, cluster=cl)
      
      dfR <- unname(Jr["df"])
      Jr_val <- if (!is.finite(dfR) || dfR <= 0) 0 else unname(Jr["J"])   # <<< AQUI: just-ID -> J_restr = 0
      C  <- unname(Jf["J"]) - Jr_val
      pC <- stats::pchisq(C, df=length(drop_set), lower.tail=FALSE)
      
      out[[length(out)+1L]] <- data.frame(
        eq=eqn, block=bk, drop_set=paste(drop_set, collapse=" + "),
        k_removed=length(drop_set),
        J_full=unname(Jf["J"]), J_restr=Jr_val,
        C_stat=C, df_C=length(drop_set), p_C=pC
      )
    }
  }
  dplyr::arrange(dplyr::bind_rows(out), dplyr::desc(C_stat))
}

swap_to_dfJ_one_keep <- function(fit_base, data, fit_eqspec, eqn,
                                 cluster_var="uf", must_keep=character(0),
                                 avoid=character(0)) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  cur <- unique(inst_by_eq[[eqn]])
  
  # candidatos globais possíveis (exclui base, atuais e avoids)
  all_pool <- setdiff(unique(fit_base$inst_terms), base_exog)
  add_cand <- setdiff(all_pool, union(cur, avoid))
  
  f <- formula(fit_base$eq[[ match(eqn, vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))) ]])
  
  best <- NULL; bestp <- -Inf; best_move <- NULL
  for (z in add_cand) {
    # tenta: (i) add z sem dropar; (ii) add z e dropar 1 para manter df_J=1
    try_sets <- list(cur_plus=unique(c(cur, z)))
    for (d in cur) try_sets[[length(try_sets)+1L]] <- setdiff(unique(c(cur, z)), d)
    
    for (S in try_sets) {
      # respeita must_keep: não pode ter dropado ninguém do must_keep
      if (any(!(must_keep %in% S))) next
      terms <- unique(c(base_exog, S))
      ivm <- try(AER::ivreg(f, instruments=reformulate(terms), data=data), silent=TRUE)
      if (inherits(ivm, "try-error")) next
      mats <- .get_XZ_u(ivm); if (is.null(mats)) next
      if ((qr(mats$Z)$rank - qr(mats$X)$rank) != 1) next  # precisa df_J == 1
      
      Jr <- .j_overid_manual(ivm, cluster = if (cluster_var %in% names(data)) data[[cluster_var]] else NULL)
      p  <- unname(Jr["p"])
      if (is.finite(p) && p > bestp) { bestp <- p; best <- S; best_move <- list(add=z, new_set=S) }
    }
  }
  if (!is.null(best)) {
    inst_by_eq[[eqn]] <- best
    message("[", eqn, "] df_J=1 com swap keeped: +", best_move$add, " | p(J)≈", signif(bestp,3))
    return(refit_quaids_try_eqlist(fit_base, data, inst_by_eq))
  } else {
    message("[", eqn, "] não achei swap com df_J=1 respeitando must_keep.")
    return(fit_eqspec)
  }
}
fit_tryK <- swap_to_dfJ_one_keep(
  fit_q, df_iv, fit_sw2, "w_despesahat6", cluster_var="uf",
  must_keep = c("IV_d_uf_X12"),
  avoid = c("IV_d_uf_X11","IV_d_uf_X14")
)
post_refit_report_auto(fit_tryK, df_iv, cluster_var="uf", R=2000, level=0.95)
fit_tryK2 <- swap_to_dfJ_one_keep(
  fit_q, df_iv, fit_tryK, "w_despesahat6", cluster_var="uf",
  must_keep = c("IV_d_uf_X12"),
  avoid = c("IV_d_uf_X11")  # agora liberamos X14, caso precise sair
)
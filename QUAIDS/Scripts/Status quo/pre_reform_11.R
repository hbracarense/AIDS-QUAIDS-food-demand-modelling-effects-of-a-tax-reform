## =========================
## PATCHES SEGUROS (DROP-IN)
## =========================

# Usa o mesmo helper .mm do seu script. Se não existir, descomente esta versão mínima:
# .mm <- function(formula_or_terms, data, add_intercept = TRUE) {
#   if (inherits(formula_or_terms, "formula")) {
#     mm <- model.matrix(formula_or_terms, data = data)
#   } else {
#     mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
#   }
#   ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
#   mm  <- mm[, ok, drop = FALSE]
#   q   <- qr(mm)
#   if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
#   mm
# }

library(openxlsx)

get_base_exog <- function(fit_obj){
  rhs_all  <- unique(fit_obj$rhs_terms)
  setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
}

# ---------- Busca subset p/ eq foco (agora exige df_J >= 0 e <= alvo) ----------
tune_focus_eq_subset <- function(fit_obj, data, current_excl,
                                 focus_eq, must_keep = NULL,
                                 choose_from = NULL,
                                 target_dfJ_max = 1,
                                 cluster_var = NULL,
                                 max_add_size = 5,
                                 max_combos   = 600L) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  if (is.null(choose_from)) choose_from <- setdiff(current_excl, must_keep)
  if (is.null(must_keep))    must_keep  <- character(0)
  
  best <- list(p = -Inf, set = NULL, dfJ = NA)
  sizes <- 0:min(length(choose_from), max_add_size)
  tested <- 0L
  
  fmla_eq <- formula(fit_obj$eq[[focus_eq]])
  
  for (s in sizes) {
    if (tested > max_combos) break
    cmbs <- utils::combn(choose_from, s, simplify = FALSE)
    for (S in cmbs) {
      tested <- tested + 1L
      excl_focus <- unique(c(must_keep, S))
      inst_terms_focus <- unique(c(base_exog, excl_focus))
      ivm <- try(AER::ivreg(fmla_eq, instruments = reformulate(inst_terms_focus), data = data),
                 silent = TRUE)
      if (inherits(ivm, "try-error")) next
      mats <- .get_XZ_u(ivm); if (is.null(mats)) next
      rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
      # Agora: exige identificação (dfJ >= 0) e limita dfJ <= alvo
      if (!is.finite(dfJ) || dfJ < 0 || dfJ > target_dfJ_max) next
      Jr <- .j_overid_manual(ivm, cluster = if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL)
      p  <- unname(Jr["p"])
      if (is.finite(p) && p > best$p) best <- list(p = p, set = excl_focus, dfJ = dfJ)
    }
  }
  if (is.null(best$set)) return(NULL)
  best
}

# ---------- Checagem/ajuste de identificação por equação ----------
ensure_identification_eqlist <- function(fit_obj, data, inst_excl_by_eq, extra_pool = NULL){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  
  if (is.null(extra_pool)) {
    # usa a união de todos excluídos informados como pool adicional
    extra_pool <- setdiff(unique(unlist(inst_excl_by_eq)), base_exog)
  }
  
  added <- list()
  for (eqname in names(eqs)) {
    # monta X e Z e checa df_J
    X <- .mm(delete.response(terms(eqs[[eqname]])), data, add_intercept = TRUE)
    Z <- .mm(reformulate(unique(c(base_exog, inst_excl_by_eq[[eqname]])), intercept = TRUE), data, add_intercept = TRUE)
    dfJ <- ncol(Z) - qr(X)$rank
    if (!is.finite(dfJ) || dfJ < 0) {
      # precisa adicionar IVs do pool até dfJ >= 0 (just-ID no limite)
      need <- 0 - dfJ
      cand <- setdiff(extra_pool, inst_excl_by_eq[[eqname]])
      if (length(cand) == 0L) stop("Sem IVs no pool para consertar identificação em ", eqname)
      k_add <- 0L
      for (z in cand) {
        inst_excl_by_eq[[eqname]] <- unique(c(inst_excl_by_eq[[eqname]], z))
        Z <- .mm(reformulate(unique(c(base_exog, inst_excl_by_eq[[eqname]])), intercept = TRUE), data, add_intercept = TRUE)
        dfJ <- ncol(Z) - qr(X)$rank
        k_add <- k_add + 1L
        if (dfJ >= 0) break
      }
      added[[eqname]] <- setdiff(inst_excl_by_eq[[eqname]], setdiff(inst_excl_by_eq[[eqname]], cand))[seq_len(k_add)]
      if (dfJ < 0) stop("Ainda under-ID em ", eqname, " mesmo após tentar adicionar do pool.")
    }
  }
  attr(inst_excl_by_eq, "added_to_fix_id") <- added
  inst_excl_by_eq
}

# ---------- Refit 3SLS eq-specific (com validação & mensagens de erro úteis) ----------
# ============ PATCH: refit 3SLS com instrumentos por equação (labels safe) ============
# requer: get_base_exog(), ensure_identification_eqlist() já definidos no seu script

refit_quaids_try_eqlist <- function(fit_obj, data, inst_excl_by_eq) {
  data <- as.data.frame(data)
  
  # Exógenos de base = RHS sem ln_*
  get_base_exog <- function(fit_obj){
    rhs_all  <- unique(fit_obj$rhs_terms)
    setdiff(rhs_all, grep("^ln_", rhs_all, value = TRUE))
  }
  base_exog <- get_base_exog(fit_obj)
  
  # Fórmulas originais K-1 (mantém LHS com underscores)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  lhs <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  
  # Garante identificação por equação (df_J >= 0), se precisar adiciona IV do pool
  if (!exists("ensure_identification_eqlist", mode = "function")) {
    stop("Defina 'ensure_identification_eqlist' antes de chamar refit_quaids_try_eqlist().")
  }
  inst_excl_by_eq <- ensure_identification_eqlist(
    fit_obj, data, inst_excl_by_eq,
    extra_pool = setdiff(unique(unlist(inst_excl_by_eq)), base_exog)
  )
  
  # Monta lista de instrumentos por equação seguindo a ORDEM das fórmulas
  inst_list <- lapply(seq_along(eqs), function(i){
    eq_lhs <- lhs[i]
    keep_excl <- unique(inst_excl_by_eq[[eq_lhs]])
    terms <- unique(c(base_exog, keep_excl))
    stats::reformulate(terms, intercept = FALSE)
  })
  
  # >>> Rótulos “seguros” para o systemfit (sem espaços/underscores)
  eq_labels <- paste0("eq", seq_along(eqs))
  names(eqs)      <- eq_labels
  names(inst_list) <- eq_labels
  
  # Chama o 3SLS com labels seguros
  fit_new <- try(systemfit(eqs, data = data, method = "3SLS", inst = inst_list),
                 silent = TRUE)
  if (inherits(fit_new, "try-error")) {
    msg <- attr(fit_new, "condition")$message
    stop("systemfit(3SLS) falhou no eq-specific. Motivo: ", msg)
  }
  
  # Reconstrói objeto tipo quaids_km1_fit (colunas = nomes originais do LHS)
  shareNames <- fit_obj$shareNames
  keep_sh    <- lhs
  yhat_mat <- sapply(fit_new$eq, stats::fitted)
  if (is.list(yhat_mat)) yhat_mat <- do.call(cbind, yhat_mat)
  fitted_included_df <- as.data.frame(yhat_mat, optional = TRUE)
  colnames(fitted_included_df) <- keep_sh
  
  fitted_full <- fitted_included_df
  omit_share  <- setdiff(shareNames, keep_sh)
  if (length(omit_share) == 1L) fitted_full[[omit_share]] <- 1 - rowSums(fitted_included_df)
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
    inst_terms  = unique(c(base_exog, unlist(inst_excl_by_eq))), # união para compat.
    fit         = fit_new,
    eq          = fit_new$eq,
    fitted_shares   = fitted_full,
    observed_shares = observed_full
  )
  class(res) <- "quaids_km1_fit"
  
  # Guarda o mapeamento eq->IVs excluídos com CHAVES = LHS originais (com "_")
  attr(res, "inst_excl_by_eq") <- inst_excl_by_eq
  res
}
# ========================= FIM DO PATCH =========================

# ---------- Hansen J por equação (eq-specific instruments) ----------
jtests_by_eq_eqlist <- function(fit_obj, data, inst_excl_by_eq, cluster_var = NULL){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  rows <- lapply(names(eqs), function(eqname){
    f <- eqs[[eqname]]
    inst_terms <- unique(c(base_exog, inst_excl_by_eq[[eqname]]))
    ivm <- try(AER::ivreg(f, instruments = reformulate(inst_terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) {
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA, J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats))
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA, J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    s_h <- try(AER::sargan(ivm), silent = TRUE)
    J_s <- if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$statistic)) else NA_real_
    p_s <- if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$p.value))   else NA_real_
    Jr  <- .j_overid_manual(ivm, cluster = cl)
    data.frame(eq=eqname, rank_X=rX, rank_Z=rZ, df_J=dfJ,
               J_sargan=J_s, p_sargan=p_s, J_hansen=unname(Jr["J"]), p_hansen=unname(Jr["p"]))
  })
  dplyr::bind_rows(rows)
}

# =========================
# DRIVER (TOP TIER, EQ-SPEC)
# =========================

stopifnot(exists("fit_q"), exists("df_iv"))
eq_names <- vapply(fit_q$eq, function(m) as.character(formula(m)[[2]]), character(1))
cl_var <- if (exists("pick_cluster")) pick_cluster(df_iv) else { nm <- names(df_iv); if ("psu" %in% nm) "psu" else if ("uf" %in% nm) "uf" else NULL }

# seu set atual de IVs excluídos (ajuste se precisar)
excl_now <- c("IV_d_uf_X11","IV_d_uf_X12","IV_d_uf_X13","IV_d_uf_X15","iv_op01","iv_op02","iv_op03")

# foco na pior eq do Hansen atual (pior = w_despesahat3)
jt_base <- jtests_by_eq_manual(fit_q, df_iv, cluster_var = cl_var, verbose = FALSE)
focus_eqname <- jt_base$eq[ which.min(jt_base$p_hansen) ]
focus_idx <- match(focus_eqname, eq_names)

# força manter (pode trocar/zerar)
must_keep <- c("iv_op01","iv_op02","iv_op03")
choose_from <- setdiff(excl_now, must_keep)

best_focus <- tune_focus_eq_subset(
  fit_obj = fit_q, data = df_iv,
  current_excl = excl_now,
  focus_eq = focus_idx,
  must_keep = must_keep,
  choose_from = choose_from,
  target_dfJ_max = 1,
  cluster_var = cl_var
)

if (is.null(best_focus)) {
  message("[EQ-SPECIFIC] Nenhum subset viável com df_J <= 1 para ", focus_eqname, ". Mantendo conjunto atual nessa equação.")
}

# lista eq->excluídos: resto fica com excl_now; foco usa subset se existir
inst_excl_by_eq <- setNames(replicate(length(eq_names), excl_now, simplify = FALSE), eq_names)
if (!is.null(best_focus)) inst_excl_by_eq[[focus_eqname]] <- best_focus$set

# Refit 3SLS eq-specific (com validação de identificação)
fit_eqspec <- refit_quaids_try_eqlist(fit_q, df_iv, inst_excl_by_eq)


# Hansen J por equação (eq-specific)
cat("\n== Hansen J (eq-specific instruments) ==\n")
print(jtests_by_eq_eqlist(fit_q, df_iv, inst_excl_by_eq, cluster_var = cl_var))

# Relatório completo
post_refit_report(fit_eqspec, df_iv, cluster_var = pick_cluster(df_iv), R = 2000, level = 0.95)

# Moment diagnostics (apenas IVs excluídos), já no eq-specific
cat("\n== Moment diagnostics (only excluded, eq-specific) ==\n")
diag_tab_eqs <- iv_moment_diag(fit_eqspec, df_iv, cluster_var = cl_var, only_excluded = TRUE)
print(dplyr::arrange(diag_tab_eqs, dplyr::desc(abs_t)) %>% dplyr::group_by(eq) %>% dplyr::slice_head(n = 8), n = 36)
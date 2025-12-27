## =========================
## TOP TIER EQ-SPEC: PLUG & PLAY
## Requer: AER, systemfit, sandwich, lmtest, car, dplyr, tidyr, MASS
## Requer helpers do seu script: .get_XZ_u(), .mm(), .align_cluster(), .j_overid_manual(),
##                               pick_cluster(), fs_block_F() (sua versão robusta), 
##                               refit_quaids_try_eqlist(), ensure_identification_eqlist(),
##                               get_base_exog()
## =========================
library(openxlsx)
# ---------- 1) J por equação: autodetecta eq-specific ----------
jtests_overid_auto <- function(fit_obj, data, cluster_var = NULL, verbose = FALSE){
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_obj)
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  names(eqs) <- vapply(eqs, function(f) as.character(f[[2]]), character(1))
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  inst_excl_by_eq <- attr(fit_obj, "inst_excl_by_eq")
  use_eqspec <- is.list(inst_excl_by_eq)
  
  rows <- lapply(names(eqs), function(eqname){
    f <- eqs[[eqname]]
    inst_terms <- if (use_eqspec) {
      unique(c(base_exog, inst_excl_by_eq[[eqname]]))
    } else {
      # fallback: união (comportamento antigo)
      fit_obj$inst_terms
    }
    ivm <- try(AER::ivreg(f, instruments = reformulate(inst_terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) {
      if (verbose) message("ivreg falhou em ", eqname, ": ", attr(ivm, "condition")$message)
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA,
                        J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
    }
    mats <- .get_XZ_u(ivm)
    if (is.null(mats))
      return(data.frame(eq=eqname, rank_X=NA, rank_Z=NA, df_J=NA,
                        J_sargan=NA, p_sargan=NA, J_hansen=NA, p_hansen=NA))
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

summarize_overid <- function(jtab){
  data.frame(
    min_dfJ = suppressWarnings(min(jtab$df_J, na.rm = TRUE)),
    min_p   = suppressWarnings(min(jtab$p_hansen, na.rm = TRUE))
  )
}

# ---------- 2) Moment diagnostics: autodetecta eq-specific ----------
iv_moment_diag_auto <- function(fit_obj, data, cluster_var = NULL, only_excluded = FALSE){
  data <- as.data.frame(data)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  eqs <- lapply(fit_obj$eq, function(m) formula(m))
  inst_excl_by_eq <- attr(fit_obj, "inst_excl_by_eq")
  base_exog <- get_base_exog(fit_obj)
  use_eqspec <- is.list(inst_excl_by_eq)
  
  out <- list()
  for (f in eqs) {
    eqn <- as.character(f[[2]])
    inst_terms <- if (use_eqspec) {
      terms0 <- unique(c(base_exog, inst_excl_by_eq[[eqn]]))
      if (isTRUE(only_excluded)) setdiff(terms0, base_exog) else terms0
    } else {
      terms0 <- fit_obj$inst_terms
      if (isTRUE(only_excluded)) setdiff(terms0, base_exog) else terms0
    }
    
    # Nota: se only_excluded=TRUE, pode haver warning "more regressors than instruments".
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

# ---------- 3) Relatório que respeita eq-specific ----------
post_refit_report_auto <- function(fit_new, data, cluster_var = NULL, R = 2000, level = 0.95, seed = 123){
  cat("\n== Over-ID (Hansen J) por equação ==\n")
  jt <- jtests_overid_auto(fit_new, data, cluster_var = cluster_var, verbose = FALSE)
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
    res_ci <- try(quaids_elasticities_ci(fit_new, at="medians", w_source="observed",
                                         R=R, level=level, seed=seed), silent = TRUE)
    if (!inherits(res_ci, "try-error") && is.list(res_ci) && !is.null(res_ci$expenditure)) {
      print(res_ci$expenditure); invisible(res_ci)
    } else {
      message("(Falha ao computar IC. Reportando apenas ponto.)")
      if (exists("quaids_elasticities", mode = "function")) {
        res_pt <- try(quaids_elasticities(fit_new, at="medians", w_source="observed"), silent = TRUE)
        if (!inherits(res_pt, "try-error") && is.list(res_pt) && !is.null(res_pt$expenditure)) {
          print(res_pt$expenditure); invisible(res_pt)
        }
      }
    }
  } else {
    message("\n(Elasticidades/IC puladas: 'quaids_elasticities_ci' não encontrada.)")
  }
  invisible(jt)
}

# ---------- 4) Varredor eq-specific para equações com J ruim ----------
# Usa sua tune_focus_eq_subset() internamente para cada eq ruim.
eqspec_sweep_all <- function(fit_base, data,
                             inst_excl_start,                  # vetor para eqs "boas"
                             must_keep_global = character(0),  # IVs que não podem sair (global)
                             target_dfJ_max = 1,
                             alpha_flag = 0.05,
                             cluster_var = NULL,
                             max_add_size = 5,
                             max_combos   = 600L,
                             verbose = TRUE){
  data <- as.data.frame(data)
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) cluster_var else NULL
  
  # J atual com conjunto comum
  inst_excl_by_eq <- setNames(replicate(length(eq_names), inst_excl_start, simplify = FALSE), eq_names)
  tmp_fit <- refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)  # checa consistência
  jt <- jtests_overid_auto(tmp_fit, data, cluster_var = cl, verbose = FALSE)
  
  # Eqs problemáticas
  bad_eqs <- jt$eq[ which(jt$p_hansen < alpha_flag) ]
  if (verbose) {
    message("[eqspec_sweep_all] Equações com p(J) < ", alpha_flag, ": ",
            if (length(bad_eqs)) paste(bad_eqs, collapse=", ") else "(nenhuma)")
  }
  if (!length(bad_eqs)) return(list(inst_excl_by_eq = inst_excl_by_eq,
                                    fit = tmp_fit, jt = jt, improved = FALSE))
  
  # Tenta melhorar cada eq ruim individualmente
  for (eqname in bad_eqs) {
    choose_from <- setdiff(inst_excl_by_eq[[eqname]], must_keep_global)
    best <- tune_focus_eq_subset(
      fit_obj = fit_base, data = data,
      current_excl = inst_excl_by_eq[[eqname]],
      focus_eq = match(eqname, eq_names),
      must_keep = must_keep_global,
      choose_from = choose_from,
      target_dfJ_max = target_dfJ_max,
      cluster_var = cl,
      max_add_size = max_add_size,
      max_combos = max_combos
    )
    if (!is.null(best)) {
      inst_excl_by_eq[[eqname]] <- best$set
      if (verbose) message("  [", eqname, "] novo subset (df_J=", best$dfJ, ") com p(J)=", signif(best$p, 3))
    } else if (verbose) {
      message("  [", eqname, "] sem subset viável com df_J <=", target_dfJ_max, ". Mantendo.")
    }
  }
  
  # Refit final eq-specific
  fit_eqspec <- refit_quaids_try_eqlist(fit_base, data, inst_excl_by_eq)
  jt2 <- jtests_overid_auto(fit_eqspec, data, cluster_var = cl)
  
  list(inst_excl_by_eq = inst_excl_by_eq, fit = fit_eqspec, jt = jt2, improved = TRUE)
}

## =========================
## COMO USAR (3 passos)
## =========================

## 0) Pré-requisito:
##    - objetos fit_q (quaids_km1_fit) e df_iv no ambiente
##    - suas funções auxiliares já carregadas (as do seu script + patches prévios)

## 1) ponto de partida: conjunto comum de IVs excluídos (o seu atual)
 excl_now <- c("IV_d_uf_X11","IV_d_uf_X12","IV_d_uf_X13","IV_d_uf_X15","iv_op01","iv_op02","iv_op03")
 cl_var   <- pick_cluster(df_iv)

## 2) varrer eqs com p(J) baixo e buscar subset eq-specific com df_J <= 1
 sweep <- eqspec_sweep_all(
   fit_base         = fit_q,
   data             = df_iv,
   inst_excl_start  = excl_now,
   must_keep_global = c("iv_op01","iv_op02","iv_op03"),   # relaxe isto se necessário!
   target_dfJ_max   = 1,
   alpha_flag       = 0.05,
   cluster_var      = cl_var,
   max_add_size     = 5,
   max_combos       = 600
 )

## 3) relatório final respeitando eq-specific
 post_refit_report_auto(sweep$fit, df_iv, cluster_var = cl_var, R = 2000, level = 0.95)

## (Opcional) Moment diagnostics só com IVs excluídos — já eq-specific
 diag_eq <- iv_moment_diag_auto(sweep$fit, df_iv, cluster_var = cl_var, only_excluded = TRUE)
 print(dplyr::arrange(diag_eq, dplyr::desc(abs_t)) %>% dplyr::group_by(eq) %>% dplyr::slice_head(n=8), n=36)
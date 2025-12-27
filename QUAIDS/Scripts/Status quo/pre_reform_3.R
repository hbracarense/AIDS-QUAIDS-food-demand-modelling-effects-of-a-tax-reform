## =========================================================
##  ADD-ONS "TOP TIER" — J-tests + Robustez + Export
##  Supõe que você já tem: fit_q (quaids_km1_fit) e df_iv
## =========================================================

suppressPackageStartupMessages({
  library(AER); library(sandwich); library(lmtest); library(car)
  library(ggplot2); library(dplyr); library(tidyr); library(openxlsx)
})

## ---------- Helpers numéricos/robustos ----------
.build_mm <- function(formula_or_terms, data, add_intercept = TRUE) {
  mm <- if (inherits(formula_or_terms, "formula")) {
    model.matrix(formula_or_terms, data = data)
  } else {
    model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}
.solve_safe <- function(M, b) {
  out <- try(solve(M, b), silent = TRUE)
  if (inherits(out, "try-error") || any(!is.finite(out))) MASS::ginv(M) %*% b else out
}
.J_stat <- function(gbar, S, n, ridge_scale = 1e-8, eig_tol = 1e-10) {
  S <- (S + t(S))/2
  dmean <- mean(diag(S)); if (!is.finite(dmean) || dmean <= 0) dmean <- 1
  S_r <- S + diag(ridge_scale * dmean, ncol(S))
  val <- try(as.numeric(n * t(gbar) %*% solve(S_r, gbar)), silent = TRUE)
  if (!inherits(val, "try-error") && is.finite(val)) return(val)
  ee <- eigen(S, symmetric = TRUE)
  pos <- ee$values > (max(ee$values, na.rm = TRUE) * eig_tol)
  if (any(pos)) {
    S_psd <- ee$vectors[, pos, drop = FALSE] %*% diag(ee$values[pos], sum(pos)) %*%
      t(ee$vectors[, pos, drop = FALSE])
    S_psd <- S_psd + diag(ridge_scale * dmean, ncol(S_psd))
    val2 <- try(as.numeric(n * t(gbar) %*% solve(S_psd, gbar)), silent = TRUE)
    if (!inherits(val2, "try-error") && is.finite(val2)) return(val2)
  }
  as.numeric(n * t(gbar) %*% MASS::ginv(S) %*% gbar)
}

## ---------- J por equação: Sargan (homoc) e Hansen-J (HC3) ----------
jtests_by_eq <- function(fit_obj, data) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  inst_terms <- unique(fit_obj$inst_terms)
  eqs <- fit_obj$eq
  out <- lapply(seq_along(eqs), function(i){
    f_yx <- formula(eqs[[i]])
    y  <- model.response(model.frame(f_yx, data = data))
    X  <- .build_mm(delete.response(terms(f_yx)), data, add_intercept = TRUE)
    Z  <- .build_mm(inst_terms, data, add_intercept = TRUE)
    n  <- nrow(X); kx <- qr(X)$rank; kz <- qr(Z)$rank; dfJ <- kz - kx
    if (dfJ <= 0) {
      return(data.frame(eq = as.character(f_yx[[2]]),
                        rank_X = kx, rank_Z = kz, df_J = dfJ,
                        J_sargan = NA_real_, p_sargan = NA_real_,
                        J_hansen = NA_real_, p_hansen = NA_real_))
    }
    # 2SLS rápido
    Pz   <- Z %*% .solve_safe(t(Z) %*% Z, t(Z))
    Xhat <- Pz %*% X
    b2s  <- .solve_safe(t(X) %*% Xhat, t(Xhat) %*% y)
    uhat <- as.numeric(y - X %*% b2s)
    
    # Sargan clássico: e'Z(Z'Z)^{-1}Z'e / s2
    s2   <- sum(uhat^2) / (n - kx)
    J_s  <- as.numeric(t(uhat) %*% Z %*% .solve_safe(t(Z) %*% Z, t(Z) %*% uhat)) / s2
    p_s  <- 1 - pchisq(J_s, df = dfJ)
    
    # Hansen-J (HC3): momentos e ẑ: g_t = z_t * u_t; S = (1/n)Σg_t g_t'
    # HC3 ≈ reponderar por (1 - h_t)^{-2}; usa h dos 2SLS
    h    <- diag(X %*% .solve_safe(t(X) %*% Xhat, t(Xhat) %*% X))
    w3   <- 1 / (1 - pmin(h, 0.999))^2
    # em jtests_by_eq(), troque:
    # G <- (Z * uhat) * w3
    # por:
    G <- sweep(sweep(Z, 1, uhat, `*`), 1, w3, `*`)
    
    gbar <- colMeans(G)
    S    <- crossprod(scale(G, scale = FALSE, center = TRUE)) / n
    J_h  <- .J_stat(gbar, S, n)
    p_h  <- 1 - pchisq(J_h, df = dfJ)
    
    data.frame(eq = as.character(f_yx[[2]]),
               rank_X = kx, rank_Z = kz, df_J = dfJ,
               J_sargan = J_s, p_sargan = p_s,
               J_hansen = J_h, p_hansen = p_h)
  })
  do.call(rbind, out)
}

## ---------- Força dos IVs (bloco) + partial R² (já usados) ----------
fs_block_F <- function(fit_obj, data, robust = TRUE){
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  stopifnot(length(iv_excl) > 0)
  form_fs <- function(y) reformulate(c(controls, iv_excl), response = y)
  get_one <- function(y){
    m  <- lm(form_fs(y), data = data)
    K  <- paste(iv_excl, "= 0")
    lh <- car::linearHypothesis(m, K,
                                vcov = if (robust) sandwich::vcovHC(m, type="HC3") else NULL,
                                test = "F")
    data.frame(endog=y, F=unname(lh$F[2]), df1=unname(lh$Df[2]),
               df2=unname(lh$Res.Df[2]), p=unname(lh$`Pr(>F)`[2]))
  }
  do.call(rbind, lapply(endo_ln, get_one))
}
partial_R2_exclIV <- function(fit_obj, data){
  rhs_all  <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  stopifnot(length(iv_excl) > 0)
  res <- lapply(endo_ln, function(y){
    mfY <- lm(reformulate(controls, response=y), data = data)
    y_t <- residuals(mfY)
    Zt  <- .build_mm(iv_excl, data, add_intercept = FALSE)
    if (length(controls)) {
      Ct  <- .build_mm(controls, data, add_intercept = TRUE)
      Zt  <- resid(lm(Zt ~ Ct))
    }
    R2  <- summary(lm(y_t ~ Zt))$r.squared
    data.frame(endog=y, partial_R2=as.numeric(R2))
  })
  do.call(rbind, res)
}

## ---------- OLS w3: cluster/ponderado (desenho) ----------
ols_w3_design_check <- function(fit_obj, df) {
  eq_lhs <- fit_obj$shareNames[3]
  f_w3 <- try({
    lhs_vec <- vapply(fit_obj$eq, function(m) as.character(formula(m)[[2]]), character(1))
    idx <- which(lhs_vec == eq_lhs)[1]; formula(fit_obj$eq[[idx]])
  }, silent = TRUE)
  if (inherits(f_w3, "try-error")) {
    rhs <- fit_obj$rhs_terms[fit_obj$rhs_terms %in% names(df)]
    f_w3 <- reformulate(rhs, response = eq_lhs)
  }
  m_ols <- lm(f_w3, data = df)
  
  cluster_id <- if ("psu" %in% names(df)) df[["psu"]] else
    if ("uf"  %in% names(df)) df[["uf"]]  else
      if ("f_reg" %in% names(df)) df[["f_reg"]] else NULL
  
  vc_cl <- if (!is.null(cluster_id)) sandwich::vcovCL(m_ols, cluster = cluster_id, type = "HC3")
  else sandwich::vcovHC(m_ols, type = "HC3")
  
  out <- list(ols = m_ols,
              coeftest_cluster = lmtest::coeftest(m_ols, vcov = vc_cl))
  
  if ("peso" %in% names(df)) {
    m_w <- lm(f_w3, data = df, weights = df[["peso"]])
    vcw <- if (!is.null(cluster_id)) sandwich::vcovCL(m_w, cluster = cluster_id, type = "HC3")
    else sandwich::vcovHC(m_w, type = "HC3")
    out$ols_weighted <- m_w
    out$coeftest_weighted_cluster <- lmtest::coeftest(m_w, vcov = vcw)
  }
  out
}

## ---------- Tidy + heatmap (para figuras) ----------
tidy_elas_safe <- function(elas, share_names, price_names,
                           stat_set = c("est","se","z","p","lwr","upr")) {
  first <- elas[[ stat_set[ stat_set %in% names(elas) ][1] ]]
  rn <- rownames(first); cn <- colnames(first)
  use_rows <- intersect(share_names, rn); if (!length(use_rows)) use_rows <- rn
  use_cols <- intersect(price_names,  cn); stopifnot(length(use_cols) > 0)
  make_long <- function(stat){
    if (!stat %in% names(elas)) return(NULL)
    M <- as.matrix(elas[[stat]])[use_rows, use_cols, drop = FALSE]
    df <- as.data.frame(M, stringsAsFactors = FALSE); df$eq <- rownames(df)
    tidyr::pivot_longer(df, cols = dplyr::all_of(use_cols),
                        names_to = "price", values_to = stat)
  }
  pieces <- Filter(Negate(is.null), lapply(stat_set, make_long))
  out <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("eq","price")), pieces)
  out$eq    <- factor(out$eq,    levels = share_names)
  out$price <- factor(out$price, levels = price_names)
  out[order(out$eq, out$price), ]
}
plot_heatmap <- function(df, title = "", limits = NULL, label_digits = 2) {
  df2 <- df
  if (!"sig" %in% names(df2)) {
    if ("p" %in% names(df2) && any(!is.na(df2$p))) {
      df2$sig <- dplyr::case_when(
        is.na(df2$p) ~ "", df2$p < .01 ~ "***", df2$p < .05 ~ "**",
        df2$p < .10 ~ "*", TRUE ~ ""
      )
    } else if (all(c("lwr","upr") %in% names(df2))) {
      df2$sig <- dplyr::case_when(
        is.na(df2$lwr)|is.na(df2$upr) ~ "", df2$lwr > 0 | df2$upr < 0 ~ "*", TRUE ~ ""
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

## ---------- Sensibilidades (pontos de avaliação) ----------
run_sensitivities <- function(fit_obj, data, R = 2000, level = 0.95, seed = 123){
  px_med <- vapply(data[, fit_obj$priceNames, drop=FALSE], median, numeric(1), na.rm=TRUE)
  x_med  <- median(data$gasto_total_atualhat, na.rm=TRUE)
  list(
    means_obs = quaids_elasticities_ci(fit_obj, at="means",   w_source="observed", R=R, level=level, seed=seed),
    med_obs   = quaids_elasticities_ci(fit_obj, at="medians", w_source="observed", R=R, level=level, seed=seed),
    fixed_fit = quaids_elasticities_ci(fit_obj, at="medians", w_source="fitted",
                                       point=list(prices=px_med, x=x_med), R=R, level=level, seed=seed)
  )
}

## ---------- EXPORTA tudo em um xlsx ----------
export_top_tier <- function(path_xlsx, fit_obj, df, res_ci,
                            sens = NULL, jtab = NULL, fs = NULL, pR2 = NULL) {
  wb <- createWorkbook()
  addWorksheet(wb, "Expenditure"); writeData(wb, "Expenditure", res_ci$expenditure)
  
  M_tidy <- tidy_elas_safe(res_ci$marshallian, fit_obj$shareNames, fit_obj$priceNames)
  H_tidy <- tidy_elas_safe(res_ci$hicksian,    fit_obj$shareNames, fit_obj$priceNames)
  addWorksheet(wb, "Marshall_tidy"); writeData(wb, "Marshall_tidy", M_tidy)
  addWorksheet(wb, "Hicks_tidy");    writeData(wb, "Hicks_tidy",    H_tidy)
  
  if (!is.null(jtab)) { addWorksheet(wb, "J_tests"); writeData(wb, "J_tests", jtab) }
  if (!is.null(fs))   { addWorksheet(wb, "FS_F_HC3"); writeData(wb, "FS_F_HC3", fs) }
  if (!is.null(pR2))  { addWorksheet(wb, "FS_partial_R2"); writeData(wb, "FS_partial_R2", pR2) }
  
  if (!is.null(sens)) {
    addWorksheet(wb, "Sens_means_M");  writeData(wb, "Sens_means_M",
                                                 tidy_elas_safe(sens$means_obs$marshallian, fit_obj$shareNames, fit_obj$priceNames))
    addWorksheet(wb, "Sens_means_H");  writeData(wb, "Sens_means_H",
                                                 tidy_elas_safe(sens$means_obs$hicksian,    fit_obj$shareNames, fit_obj$priceNames))
    addWorksheet(wb, "Sens_fixed_M");  writeData(wb, "Sens_fixed_M",
                                                 tidy_elas_safe(sens$fixed_fit$marshallian, fit_obj$shareNames, fit_obj$priceNames))
    addWorksheet(wb, "Sens_fixed_H");  writeData(wb, "Sens_fixed_H",
                                                 tidy_elas_safe(sens$fixed_fit$hicksian,    fit_obj$shareNames, fit_obj$priceNames))
  }
  saveWorkbook(wb, path_xlsx, overwrite = TRUE)
  invisible(list(M_tidy=M_tidy, H_tidy=H_tidy, path=normalizePath(path_xlsx)))
}

## ===================== EXECUTAR =====================

## 1) J por equação
jtab <- jtests_by_eq(fit_q, df_iv)
print(jtab)

## 2) Força dos IVs (para a seção de diagnóstico)
fs   <- fs_block_F(fit_q, df_iv, robust = TRUE); print(fs)
pR2  <- partial_R2_exclIV(fit_q, df_iv);         print(pR2)

## 3) Desenho amostral (OLS w3)
design_out <- ols_w3_design_check(fit_q, df_iv)
print(design_out$coeftest_cluster)
if (!is.null(design_out$coeftest_weighted_cluster)) {
  cat("\n== OLS w3 ponderado + cluster ==\n"); print(design_out$coeftest_weighted_cluster)
}

## 4) Sensibilidades dos pontos de avaliação (opcional, mas ótimo p/ robustez)
sens <- run_sensitivities(fit_q, df_iv, R = 2000, level = 0.95, seed = 123)
print(sens)

## 5) Heatmaps (para Figuras)
M_tidy <- tidy_elas_safe(res_ci$marshallian, fit_q$shareNames, fit_q$priceNames)
H_tidy <- tidy_elas_safe(res_ci$hicksian,    fit_q$shareNames, fit_q$priceNames)
M_tidy_trans <- M_tidy
M_tidy_trans <-  M_tidy_trans %>% mutate(eq = case_when(
  eq == 'w_despesahat1' ~ 'w_expenditure_hat1',
  eq == 'w_despesahat2' ~ 'w_expenditure_hat2',
  eq == 'w_despesahat3' ~ 'w_expenditure_hat3',
  eq == 'w_despesahat4' ~ 'w_expenditure_hat4',
  eq == 'w_despesahat5' ~ 'w_expenditure_hat5',
  eq == 'w_despesahat6' ~ 'w_expenditure_hat6'
),
price = case_when(
  price == 'preco_por_kg1' ~ 'price_per_kg_1',
  price == 'preco_por_kg2' ~ 'price_per_kg_2',
  price == 'preco_por_kg3' ~ 'price_per_kg_3',
  price == 'preco_por_kg4' ~ 'price_per_kg_4',
  price == 'preco_por_kg5' ~ 'price_per_kg_5',
  price == 'preco_por_kg6' ~ 'price_per_kg_6'
))

H_tidy_trans <- H_tidy
H_tidy_trans <-  H_tidy_trans %>% mutate(eq = case_when(
  eq == 'w_despesahat1' ~ 'w_expenditure_hat1',
  eq == 'w_despesahat2' ~ 'w_expenditure_hat2',
  eq == 'w_despesahat3' ~ 'w_expenditure_hat3',
  eq == 'w_despesahat4' ~ 'w_expenditure_hat4',
  eq == 'w_despesahat5' ~ 'w_expenditure_hat5',
  eq == 'w_despesahat6' ~ 'w_expenditure_hat6'
),
price = case_when(
  price == 'preco_por_kg1' ~ 'price_per_kg_1',
  price == 'preco_por_kg2' ~ 'price_per_kg_2',
  price == 'preco_por_kg3' ~ 'price_per_kg_3',
  price == 'preco_por_kg4' ~ 'price_per_kg_4',
  price == 'preco_por_kg5' ~ 'price_per_kg_5',
  price == 'preco_por_kg6' ~ 'price_per_kg_6'
))



print( plot_heatmap(M_tidy_trans, "Elasticities Marshall") )
print( plot_heatmap(H_tidy_trans, "Elasticities Hicks") )

## 6) Exporta tudo em um Excel (anexos do paper)
#ex_out <- export_top_tier("elasticidades_top_tier.xlsx",
#                          fit_obj = fit_q, df = df_iv, res_ci = res_ci,
#                          sens = sens, jtab = jtab, fs = fs, pR2 = pR2)
#cat("\nArquivo salvo em:\n", ex_out$path, "\n")

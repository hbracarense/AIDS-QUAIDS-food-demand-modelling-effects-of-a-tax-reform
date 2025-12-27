## ============================================================
##  PACOTES
## ============================================================
suppressPackageStartupMessages({
  library(AER); library(sandwich); library(lmtest); library(car)
  library(ggplot2); library(tidyr); library(dplyr); library(openxlsx)
})

## ============================================================
## 1) STATUS DE IDENTIFICAÇÃO & SOBRE-ID + J (HC3)
##    - rank(Z) e rank(X) por equação; J de Sargan/Hansen se df_J>0
## ============================================================
.build_mm <- function(formula_or_terms, data, add_intercept = TRUE) {
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  # remove colunas com var ~ 0 / NA; resolve colinearidade por QR
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}

overid_status_by_eq <- function(fit_obj, data) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  inst_terms <- unique(fit_obj$inst_terms)
  eqs <- fit_obj$eq
  out <- lapply(seq_along(eqs), function(i){
    f_yx <- formula(eqs[[i]])
    # X = regressors da equação i
    X <- .build_mm(delete.response(terms(f_yx)), data, add_intercept = TRUE)
    # Z = instrumentos (inclui exógenos + IVs excluídos)
    Z <- .build_mm(inst_terms, data, add_intercept = TRUE)
    rX <- qr(X)$rank; rZ <- qr(Z)$rank
    dfJ <- rZ - rX
    # J (Sargan e Hansen-HC3) via ivreg se sobre-ID
    J_sarg <- J_hans <- p_sarg <- p_hans <- NA_real_
    if (dfJ > 0) {
      ivm <- AER::ivreg(f_yx, instruments = reformulate(inst_terms), data = data)
      s_h <- try(AER::sargan(ivm), silent = TRUE)
      s_r <- try(AER::sargan(ivm, vcov = sandwich::vcovHC(ivm, type = "HC3")), silent = TRUE)
      if (!inherits(s_h, "try-error")) {
        J_sarg <- unname(as.numeric(s_h$statistic)); p_sarg <- unname(as.numeric(s_h$p.value))
      }
      if (!inherits(s_r, "try-error")) {
        J_hans <- unname(as.numeric(s_r$statistic)); p_hans <- unname(as.numeric(s_r$p.value))
      }
    }
    data.frame(
      eq      = as.character(f_yx[[2]]),
      rank_X  = rX, rank_Z = rZ, df_J = dfJ,
      J_sargan = J_sarg, p_sargan = p_sarg,
      J_hansen = J_hans, p_hansen = p_hans,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ============================================================
## 2) FORÇA DOS IVs: F (HC3) já em bloco + PARTIAL R^2 dos IVs excluídos
## ============================================================
fs_block_F <- function(fit_obj, data, robust = TRUE){
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  rhs_all <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)   # instrumentos incluídos (exógenos)
  iv_excl  <- setdiff(inst_all, controls)    # IVs excluídos (relevância)
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

partial_R2_exclIV <- function(fit_obj, data){
  # partial R^2 dos IVs excluídos (residualizando nos exógenos incluídos)
  rhs_all <- fit_obj$rhs_terms
  inst_all <- fit_obj$inst_terms
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  controls <- intersect(inst_all, rhs_all)
  iv_excl  <- setdiff(inst_all, controls)
  if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
  res <- lapply(endo_ln, function(y){
    # residualiza Y e Z_excl nos controles
    mfY <- lm(reformulate(controls, response=y), data = data)
    y_t <- residuals(mfY)
    Zt  <- .build_mm(iv_excl, data, add_intercept = FALSE)
    if (length(controls)) {
      Ct  <- .build_mm(controls, data, add_intercept = TRUE)
      Zt  <- resid(lm(Zt ~ Ct))
    }
    fit <- lm(y_t ~ Zt)
    R2  <- summary(fit)$r.squared
    data.frame(endog=y, partial_R2=as.numeric(R2))
  })
  do.call(rbind, res)
}

## ============================================================
## 3) NEGATIVIDADE/ CURVATURA (Slutsky) a partir das Hicksianas
##     S_ij = e^H_ij * (q_i / p_j), onde q_i = w_i * x / p_i
##     Usa preços ponto/medianas; w_i por default = médias das shares observadas.
## ============================================================
slutsky_from_hicks <- function(H_mat, prices, x, shares) {
  K <- length(shares)
  stopifnot(nrow(H_mat)==K, ncol(H_mat)==K, length(prices)==K)
  q <- (shares * x) / prices
  S <- matrix(NA_real_, K, K, dimnames = dimnames(H_mat))
  for (i in seq_len(K)) for (j in seq_len(K)) S[i,j] <- H_mat[i,j] * (q[i] / prices[j])
  S_sym <- (S + t(S))/2
  eig <- eigen(S_sym, symmetric = TRUE, only.values = TRUE)$values
  list(S = S, S_sym = S_sym, eig = eig,
       neg_semi_def = all(eig <= 1e-8),
       max_eig = max(eig), min_eig = min(eig),
       asymmetry = mean(abs(S - t(S))))
}

run_negativity_test <- function(res_ci, df_iv, fit_obj,
                                x_point = c("median","mean")[1],
                                price_point = c("median","mean")[1],
                                share_source = c("observed","median")[1]) {
  H <- res_ci$hicksian$est
  priceNames <- fit_obj$priceNames
  shareNames <- fit_obj$shareNames
  p <- vapply(df_iv[, priceNames, drop=FALSE],
              if (price_point=="median") median else mean, numeric(1), na.rm=TRUE)
  x <- if (x_point=="median") median(df_iv$gasto_total_atualhat, na.rm=TRUE)
  else mean(df_iv$gasto_total_atualhat, na.rm=TRUE)
  w <- vapply(df_iv[, shareNames, drop=FALSE],
              if (share_source=="median") median else mean, numeric(1), na.rm=TRUE)
  w <- pmax(pmin(w, 0.999), 1e-6)  # sanity
  slutsky_from_hicks(H, prices = p, x = x, shares = w)
}

## ============================================================
## 4) PESOS / DESENHO AMOSTRAL (checagem rápida em OLS da w3)
##    - Se existir 'peso' e 'psu', cluster por psu; senão cluster por 'uf' ou 'f_reg'
## ============================================================
ols_w3_design_check <- function(fit_obj, df) {
  # fórmula segura para w3
  eq_lhs <- fit_obj$shareNames[3]
  f_w3 <- try({
    lhs_vec <- vapply(fit_obj$eq, function(m) as.character(formula(m)[[2]]), character(1))
    idx <- which(lhs_vec == eq_lhs)[1]
    formula(fit_obj$eq[[idx]])
  }, silent = TRUE)
  if (inherits(f_w3, "try-error")) {
    rhs <- fit_obj$rhs_terms[fit_obj$rhs_terms %in% names(df)]
    f_w3 <- reformulate(rhs, response = eq_lhs)
  }
  
  # OLS simples
  m_ols <- lm(f_w3, data = df)
  
  # cluster id (ordem: psu, uf, f_reg) com indexação segura
  cluster_id <- if ("psu" %in% names(df)) df[["psu"]] else
    if ("uf"  %in% names(df)) df[["uf"]]  else
      if ("f_reg" %in% names(df)) df[["f_reg"]] else
        NULL
  
  vc_cl <- if (!is.null(cluster_id)) sandwich::vcovCL(m_ols, cluster = cluster_id, type = "HC3") else
    sandwich::vcovHC(m_ols, type = "HC3")
  
  out <- list(
    ols = m_ols,
    coeftest_cluster = lmtest::coeftest(m_ols, vcov = vc_cl)
  )
  
  # Versão ponderada (se existir coluna 'peso')
  # dentro de ols_w3_design_check(...)
  if ("peso" %in% names(df)) {
    m_w <- lm(f_w3, data = df, weights = peso)  # <<-- aqui
    vcw <- if (!is.null(cluster_id)) sandwich::vcovCL(m_w, cluster = cluster_id, type = "HC3")
    else sandwich::vcovHC(m_w, type = "HC3")
    out$ols_weighted <- m_w
    out$coeftest_weighted_cluster <- lmtest::coeftest(m_w, vcov = vcw)
  }
  
  
  out
}


## ============================================================
## 5) SENSIBILIDADES DE ELASTICIDADES (means/medians/fixed point)
## ============================================================
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

## ============================================================
## 6) TIDY + HEATMAPS (Marshall/Hicks) e EXPORTAÇÃO
## ============================================================
tidy_elas <- function(elas_list, good_names, price_names, stat_set = c("est","p","lwr","upr")) {
  mats <- lapply(stat_set, function(s) as.data.frame(elas_list[[s]], check.names = FALSE))
  names(mats) <- stat_set
  for (nm in names(mats)) { mats[[nm]]$good <- rownames(mats[[nm]]) }
  td <- Reduce(function(a,b) full_join(a, b, by="good"), mats)
  td <- pivot_longer(td, cols = all_of(price_names), names_to = "price", values_to = "est") %>%
    rename(est = est)
  # re-anexar as demais métricas
  for (s in setdiff(stat_set, "est")) {
    td_s <- as.data.frame(elas_list[[s]], check.names = FALSE); td_s$good <- rownames(td_s)
    td <- left_join(td, pivot_longer(td_s, cols = all_of(price_names),
                                     names_to = "price", values_to = s), by = c("good","price"))
  }
  td <- td %>% mutate(sig = ifelse(!is.na(p) & p < 0.05, TRUE, FALSE))
  td
}

plot_heatmap <- function(tidy_df, title = "Elasticity (Marshall)") {
  ggplot(tidy_df, aes(price, good, fill = est)) +
    geom_tile(color = "white") +
    geom_point(data = subset(tidy_df, sig), aes(price, good), shape = 8, size = 1.8, color = "black") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    labs(x = "Price", y = "Good", fill = "Elasticity", title = title,
         subtitle = "* indicates p < 0.05") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

export_elas <- function(res_ci, fit_obj, path_xlsx = "elasticidades_top_tier.xlsx") {
  wb <- createWorkbook()
  # despesa
  addWorksheet(wb, "Expenditure")
  writeData(wb, "Expenditure", res_ci$expenditure)
  # marshall
  M_tidy <- tidy_elas(res_ci$marshallian, fit_obj$shareNames, fit_obj$priceNames,
                      stat_set = c("est","se","z","p","lwr","upr"))
  addWorksheet(wb, "Marshall_tidy")
  writeData(wb, "Marshall_tidy", M_tidy)
  # hicks
  H_tidy <- tidy_elas(res_ci$hicksian, fit_obj$shareNames, fit_obj$priceNames,
                      stat_set = c("est","se","z","p","lwr","upr"))
  addWorksheet(wb, "Hicks_tidy")
  writeData(wb, "Hicks_tidy", H_tidy)
  saveWorkbook(wb, path_xlsx, overwrite = TRUE)
  invisible(list(M_tidy=M_tidy, H_tidy=H_tidy, path=normalizePath(path_xlsx)))
}

## ============================================================
## 7) RODAR TUDO (exemplo com res_ci já estimado)
##    - Ajuste seeds/R conforme desejar
## ============================================================
# 7.1) Elasticidades + IC (se ainda não tiver)
# res_ci <- quaids_elasticities_ci(fit_q, at="medians", w_source="observed", R=2000, level=0.95, seed=123)

# 7.2) Identificação & J por equação
ov_stat <- overid_status_by_eq(fit_q, df_iv)
print(ov_stat)

# 7.3) Força dos IVs (F-HC3) e partial R^2
fs_hc3 <- fs_block_F(fit_q, df_iv, robust = TRUE);  print(fs_hc3)
pR2    <- partial_R2_exclIV(fit_q, df_iv);          print(pR2)

# 7.4) Negatividade/curvatura (no mesmo ponto do res_ci)
neg_chk <- run_negativity_test(res_ci, df_iv, fit_q,
                               x_point = "median", price_point = "median", share_source = "observed")
cat("\nNegatividade Slutsky (compensado): ",
    if (neg_chk$neg_semi_def) "OK (semi-definida negativa)" else "FALHA", "\n",
    "max eigen = ", round(neg_chk$max_eig,6), " | min eigen = ", round(neg_chk$min_eig,6),
    " | assimetria média = ", round(neg_chk$asymmetry,6), "\n", sep="")

# 7.5) Checagem com pesos/desenho (OLS w3)
design_out <- ols_w3_design_check(fit_q, df_iv)

print(design_out$coeftest_cluster)
if (!is.null(design_out$coeftest_weighted_cluster)) {
  cat("\n== OLS w3 ponderado + cluster ==\n"); print(design_out$coeftest_weighted_cluster)
}

# 7.6) Sensibilidades (opcional; muitas já rodadas)
sens <- run_sensitivities(fit_q, df_iv, R = 2000, level = 0.95, seed = 123)

# 7.7) Heatmaps e exportação
tidy_elas_safe <- function(elas, share_names = NULL, price_names = NULL,
                           stat_set = c("est","se","z","p","lwr","upr")) {
  stopifnot(is.list(elas))
  
  # pega dimnames da primeira matriz disponível
  first_mat <- elas[[ stat_set[ stat_set %in% names(elas) ][1] ]]
  if (is.null(first_mat)) stop("Nenhuma matriz encontrada entre: ", paste(stat_set, collapse=", "))
  
  rn <- rownames(first_mat); cn <- colnames(first_mat)
  if (is.null(rn) || is.null(cn)) stop("As matrizes precisam ter rownames e colnames.")
  
  # nomes de referência (se não passados)
  if (is.null(share_names)) share_names <- rn
  if (is.null(price_names)) price_names <- cn
  
  # saneia espaços/acentos acidentais
  rn <- trimws(rn); cn <- trimws(cn)
  share_names <- trimws(share_names); price_names <- trimws(price_names)
  
  use_rows <- intersect(share_names, rn)
  use_cols <- intersect(price_names,  cn)
  if (!length(use_cols)) {
    stop("Nenhum dos 'price_names' aparece nas matrizes. Colunas presentes: ",
         paste(cn, collapse=", "))
  }
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
  # ordena opcionalmente
  out$eq    <- factor(out$eq,    levels = share_names)
  out$price <- factor(out$price, levels = price_names)
  out[order(out$eq, out$price), ]
}

M_tidy <- tidy_elas_safe(
  res_ci$marshallian,
  share_names = fit_q$shareNames,
  price_names = fit_q$priceNames,
  stat_set = c("est","p","lwr","upr")
)

H_tidy <- tidy_elas_safe(
  res_ci$hicksian,
  share_names = fit_q$shareNames,
  price_names = fit_q$priceNames,
  stat_set = c("est","p","lwr","upr")
)


plot_heatmap <- function(df, title = "", limits = NULL, label_digits = 2,
                         use_ci_when_p_missing = TRUE) {
  stopifnot(all(c("eq","price","est") %in% names(df)))
  df2 <- df
  
  # cria 'sig' se não existir
  if (!"sig" %in% names(df2)) {
    if ("p" %in% names(df2) && any(!is.na(df2$p))) {
      df2$sig <- dplyr::case_when(
        is.na(df2$p)        ~ "",
        df2$p < 0.01        ~ "***",
        df2$p < 0.05        ~ "**",
        df2$p < 0.10        ~ "*",
        TRUE                ~ ""
      )
    } else if (use_ci_when_p_missing && all(c("lwr","upr") %in% names(df2))) {
      # significância via IC que não cruza zero
      df2$sig <- dplyr::case_when(
        is.na(df2$lwr) | is.na(df2$upr) ~ "",
        df2$lwr > 0 | df2$upr < 0       ~ "*",  # marque como quiser (poderia mapear níveis)
        TRUE                            ~ ""
      )
    } else {
      df2$sig <- ""
    }
  }
  
  # rótulo com estrelas
  df2$label <- sprintf(paste0("%.", label_digits, "f%s"), df2$est, df2$sig)
  
  # plot
  library(ggplot2)
  p <- ggplot(df2, aes(x = price, y = eq, fill = est)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(limits = limits, midpoint = 0, na.value = "grey85") +
    geom_text(aes(label = label), size = 3) +
    labs(title = title, x = "Price", y = "Good", fill = "Elasticity") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


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


print(plot_heatmap(M_tidy_trans, "Elasticity - Marshall"))
print(plot_heatmap(H_tidy_trans, "Elasticity - Hicks"))

#exp_out <- export_elas(res_ci, fit_q, path_xlsx = "elasticidades_top_tier.xlsx")
#cat("\nArquivo Excel salvo em:\n", exp_out$path, "\n")

## ============================================================
##  FIM
## ============================================================

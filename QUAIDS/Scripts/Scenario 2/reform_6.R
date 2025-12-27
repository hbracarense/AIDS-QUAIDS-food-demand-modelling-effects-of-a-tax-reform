# ========================
# PACOTES
# ========================
suppressPackageStartupMessages({
  library(AER)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(car)
})

# ========================
# (0*) GARANTIA: monta os 2SLS base se não existir 'fits_iv'
# ========================
if (!exists("fits_iv")) {
  fits_iv <- map(share_vars, ~ fit_iv_weighted(
    dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
    iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
    z_var = z_var, z2_var = z2_var,
    seasonal_controls = seasonal_controls_full,
    use_spline = FALSE, spline_df = 4
  ))
  names(fits_iv) <- share_vars
}

# ========================
# (1) DIAGNÓSTICOS DE IV (Weak / Wu-Hausman / Sargan) a partir do 2SLS (AER::ivreg)
# ========================
extract_iv_diags <- function(fits_iv, dat = NULL) {
  # tibble vazio com as colunas certas (evita erro no relocate/select etc.)
  empty <- tibble::tibble(
    share     = character(),
    teste     = character(),
    statistic = double(),
    df1       = double(),
    df2       = double(),
    `p-value` = double()
  )
  
  purrr::map_dfr(names(fits_iv), function(v) {
    dm <- fits_iv[[v]]$summary$diagnostics
    
    # se não veio, tenta reestimar SEM pesos (algumas versões só calculam assim)
    if ((is.null(dm) || inherits(dm, "try-error")) && !is.null(dat)) {
      f  <- formula(fits_iv[[v]]$fit)
      su <- try(summary(AER::ivreg(f, data = dat), diagnostics = TRUE), silent = TRUE)
      dm <- if (inherits(su, "try-error")) NULL else su$diagnostics
    }
    
    if (is.null(dm)) {
      empty[0, ]  # devolve 0 linhas, mas com as colunas corretas
    } else {
      tibble::as_tibble(dm, rownames = "teste") |>
        dplyr::mutate(share = v, .before = 1)
    }
  })
}

iv_diags <- extract_iv_diags(fits_iv)
print(iv_diags)

# ========================
# (2) & (3) 3SLS MÍNIMO — VARIANTE A (z, sem z2) E VARIANTE B (spline(z))
# ========================
fit_3sls_minA <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = 4,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE
)

fit_3sls_minB <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = TRUE, spline_df = 4,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE
)

# Comparação in-sample de RMSE ponderado: 2SLS vs 3SLS_minA vs 3SLS_minB
obs_mat <- as.data.frame(dat[, share_vars])
w_vec   <- dat[[weights_col]]

rmse_w_2sls <- sapply(share_vars, function(v) {
  res <- resid(fits_iv[[v]]$fit)
  sqrt(weighted.mean(res^2, w_vec))
})

pred_minA <- predict_3sls(fit_3sls_minA, dat)
rmse_w_3sls_minA <- sapply(share_vars, function(v) {
  sqrt(weighted.mean((obs_mat[[v]] - pred_minA[[v]])^2, w_vec))
})

pred_minB <- predict_3sls(fit_3sls_minB, dat)
rmse_w_3sls_minB <- sapply(share_vars, function(v) {
  sqrt(weighted.mean((obs_mat[[v]] - pred_minB[[v]])^2, w_vec))
})

comp_in_tbl <- tibble(
  share = share_vars,
  rmse_w_2sls      = as.numeric(rmse_w_2sls),
  rmse_w_3sls_minA = as.numeric(rmse_w_3sls_minA),
  rmse_w_3sls_minB = as.numeric(rmse_w_3sls_minB)
)
print(comp_in_tbl)

# ========================
# (4) TESTES CONJUNTOS PARA PREÇOS (2SLS com VCOV robusto HC1)
# ========================
joint_price_tests <- bind_rows(lapply(share_vars, function(v) {
  fit <- fits_iv[[v]]$fit
  cn  <- names(coef(fit))
  price_cn <- cn[grepl("^ln_preco_com_reforma2", cn)]
  if (length(price_cn) == 0) return(tibble(share = v, F = NA, df = NA, p = NA))
  hyp <- paste0(price_cn, " = 0")
  lh  <- car::linearHypothesis(fit, hyp, white.adjust = "hc1")
  tibble(share = v,
         F  = as.numeric(lh[2, "F"]),
         df = paste0(lh[2, "Df"], collapse = ","),
         p  = as.numeric(lh[2, "Pr(>F)"]))
}))
print(joint_price_tests)

# ========================
# (5) VALIDAÇÃO OUT-OF-SAMPLE (hold-out simples 80/20)
# ========================
set.seed(123)
n   <- nrow(dat)
ix  <- sample.int(n, size = floor(0.8 * n))
trn <- dat[ix, ]
tst <- dat[-ix, ]

# Refit 2SLS no treino
fits_iv_trn <- map(share_vars, ~ fit_iv_weighted(
  dat = trn, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv_trn) <- share_vars

# Refit 3SLS mínimo (Variante A) no treino
fit_3sls_trn <- fit_3sls_system(
  dat = trn, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = 4,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE
)

# Previsões no teste
pred_iv_tst <- sapply(share_vars, function(v) {
  as.numeric(predict(fits_iv_trn[[v]]$fit, newdata = tst))
})
colnames(pred_iv_tst) <- share_vars
pred_3sls_tst <- predict_3sls(fit_3sls_trn, tst)

w_tst   <- tst[[weights_col]]
obs_tst <- as.data.frame(tst[, share_vars])

rmse_w <- function(obs, pred, w) sqrt(weighted.mean((obs - pred)^2, w))

rmse_w_iv_tst <- sapply(share_vars, function(v)
  rmse_w(obs_tst[[v]], pred_iv_tst[, v], w_tst))

rmse_w_3sls_tst <- sapply(share_vars, function(v)
  rmse_w(obs_tst[[v]], pred_3sls_tst[[v]], w_tst))

oos_tbl <- tibble(
  share = share_vars,
  rmse_w_2sls_oos = as.numeric(rmse_w_iv_tst),
  rmse_w_3sls_oos = as.numeric(rmse_w_3sls_tst)
)
print(oos_tbl)

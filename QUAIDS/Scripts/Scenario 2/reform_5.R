# ========================
# PRÉ-REQUISITOS
# ========================
library(dplyr)
library(purrr)
library(tibble)
library(broom)
library(car)

# Assume que os objetos e funções abaixo já existem no ambiente:
# dat, share_vars, price_vars, weights_col
# iv_include_patterns_full, iv_include_patterns_min, iv_exclude_patterns
# z_var, z2_var
# seasonal_controls_full, seasonal_controls_min
# fit_iv_weighted(), diagnostics_iv(), fit_3sls_system(), predict_3sls()

# ========================
# (0) 2SLS POR EQUAÇÃO (BASE PARA OS PASSOS 1, 4 e 5)
# ========================
fits_iv <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv) <- share_vars

# Métricas de acurácia in-sample do 2SLS (inclusive ponderadas)

# --- Versões robustas (sem augment) ---

w_rmse <- function(res, w) sqrt(weighted.mean(res^2, w))
w_mae  <- function(res, w) weighted.mean(abs(res), w)

diagnostics_iv <- function(fit_obj, dep, w = NULL) {
  r <- stats::residuals(fit_obj)
  f <- stats::fitted(fit_obj)
  
  ok <- is.finite(r) & is.finite(f)
  if (!is.null(w)) {
    # alinha pesos aos casos válidos
    w_use <- w[ok & is.finite(w)]
    r_ok  <- r[ok & is.finite(w)]
  } else {
    w_use <- NULL
    r_ok  <- r[ok]
  }
  
  tibble::tibble(
    share  = dep,
    rmse   = sqrt(mean(r_ok^2)),
    mae    = mean(abs(r_ok)),
    rmse_w = if (is.null(w_use)) NA_real_ else w_rmse(r_ok, w_use),
    mae_w  = if (is.null(w_use)) NA_real_ else w_mae(r_ok, w_use)
  )
}


diag_tbl <- map2_dfr(fits_iv, share_vars, ~ diagnostics_iv(.x$fit, .y, .x$weights)) |>
  select(share, rmse, mae, rmse_w, mae_w)
print(diag_tbl)

# ========================
# (1) DIAGNÓSTICOS DE IV (weak / over-id) A PARTIR DO 2SLS
# ========================
extract_iv_diags <- function(fits_iv, dat = NULL) {
  # tibble "vazio, mas tipado"
  empty <- tibble::tibble(
    share     = character(),
    teste     = character(),
    statistic = double(),
    df1       = double(),
    df2       = double(),
    `p-value` = double()
  )
  
  purrr::map_dfr(names(fits_iv), function(v) {
    # tenta pegar do objeto salvo...
    dm <- fits_iv[[v]]$summary$diagnostics
    
    # ...ou refaz sem pesos (às vezes só assim o AER calcula os diagnósticos)
    if ((is.null(dm) || inherits(dm, "try-error")) && !is.null(dat)) {
      f  <- formula(fits_iv[[v]]$fit)
      su <- try(summary(AER::ivreg(f, data = dat), diagnostics = TRUE), silent = TRUE)
      dm <- if (inherits(su, "try-error")) NULL else su$diagnostics
    }
    
    if (is.null(dm)) {
      empty[0, ]  # devolve tibble vazio com as colunas certas
    } else {
      tibble::as_tibble(dm, rownames = "teste") |>
        dplyr::mutate(share = v, .before = 1)
    }
  })
}

iv_diags <- extract_iv_diags(fits_iv)
print(iv_diags)
# # Sugestões de filtro:
# dplyr::filter(iv_diags, grepl("Weak|first-stage", teste, ignore.case = TRUE))
# dplyr::filter(iv_diags, grepl("Sargan|Overid",  teste, ignore.case = TRUE))

# ========================
# (2) & (3) 3SLS MÍNIMO — VARIANTES A (z sem z2) E B (spline(z))
# ========================

# Variante A: z (sem z2), sazonais reduzidos e IVs mínimos
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
print(summary(fit_3sls_minA$fit))

# Variante B: spline(z), sazonais reduzidos e IVs mínimos
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
print(summary(fit_3sls_minB$fit))

# Comparação in-sample (RMSE ponderado): 2SLS vs 3SLS_minA vs 3SLS_minB
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
  rmse_w_2sls     = as.numeric(rmse_w_2sls),
  rmse_w_3sls_minA = as.numeric(rmse_w_3sls_minA),
  rmse_w_3sls_minB = as.numeric(rmse_w_3sls_minB)
)
print(comp_in_tbl)

# ========================
# (4) TESTES CONJUNTOS PARA PREÇOS (2SLS c/ VCOV ROBUSTO HC1)
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
# (5) VALIDAÇÃO OUT-OF-SAMPLE (HOLD-OUT 80/20)
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

w_tst  <- tst[[weights_col]]
obs_tst <- as.data.frame(tst[, share_vars])

rmse_w_iv_tst   <- sapply(share_vars, function(v) sqrt(weighted.mean((obs_tst[[v]] - pred_iv_tst[, v])^2, w_tst)))
rmse_w_3sls_tst <- sapply(share_vars, function(v) sqrt(weighted.mean((obs_tst[[v]] - pred_3sls_tst[[v]])^2, w_tst)))

oos_tbl <- tibble(
  share = share_vars,
  rmse_w_2sls_oos = as.numeric(rmse_w_iv_tst),
  rmse_w_3sls_oos = as.numeric(rmse_w_3sls_tst)
)
print(oos_tbl)

# FIM (1→5)
# ===== (6.1) In-sample: soma das previsões 2SLS =====
pred_iv_in <- sapply(share_vars, function(v) as.numeric(predict(fits_iv[[v]]$fit)))
colnames(pred_iv_in) <- share_vars
summary(rowSums(pred_iv_in))  # deve estar ~1

# ===== (6.2) Out-of-sample: soma e normalização =====
# Usa objetos já criados no (5): pred_iv_tst, pred_3sls_tst, obs_tst, w_tst

# helper de normalização linha-a-linha (com clip opcional para [0,1])
renorm_simplex <- function(M, clip01 = TRUE) {
  if (clip01) M <- pmin(pmax(M, 0), 1)
  s <- rowSums(M)
  s[s == 0] <- 1
  M / s
}

# matrizes (3SLS vem como lista -> colar em matriz)
pred_3sls_tst_mat <- do.call(cbind, pred_3sls_tst[share_vars])
colnames(pred_3sls_tst_mat) <- share_vars

# checar soma antes da normalização
summary(rowSums(pred_iv_tst))
summary(rowSums(pred_3sls_tst_mat))

# normalizar (adding-up garantido)
pred_iv_tst_norm    <- renorm_simplex(pred_iv_tst, clip01 = TRUE)
pred_3sls_tst_norm  <- renorm_simplex(pred_3sls_tst_mat, clip01 = TRUE)

# RMSE ponderado pós-normalização
rmse_w <- function(obs, pred, w) sqrt(weighted.mean((obs - pred)^2, w))
rmse_tab_norm <- tibble::tibble(
  share = share_vars,
  rmse_w_2sls_oos_norm = sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_iv_tst_norm[, v],   w_tst)),
  rmse_w_3sls_oos_norm = sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_3sls_tst_norm[, v], w_tst))
)
print(rmse_tab_norm)

# ===== (6.3) Híbrido: escolhe o melhor (2SLS vs 3SLS) por share e reenforça soma=1 =====
winners <- ifelse(rmse_tab_norm$rmse_w_3sls_oos_norm < rmse_tab_norm$rmse_w_2sls_oos_norm,
                  "3SLS_minA", "2SLS")
pred_hybrid <- pred_iv_tst_norm
for (j in seq_along(share_vars)) {
  v <- share_vars[j]
  if (winners[j] == "3SLS_minA") pred_hybrid[, v] <- pred_3sls_tst_norm[, v]
}
pred_hybrid <- renorm_simplex(pred_hybrid, clip01 = TRUE)

# RMSE do híbrido
rmse_hybrid <- sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_hybrid[, v], w_tst))
print(tibble::tibble(share = share_vars, rmse_w_hybrid_oos = as.numeric(rmse_hybrid)))

# ===== Hansen-J robusto (sem ivpack/gmm), com QR/PCA que não dependem de fórmula =====
library(AER)      # precisa para model.matrix(..., component="regressors"/"instruments")
library(dplyr)
library(purrr)

# --- util: pega y,X,Z,w do ivreg ---
.get_mats <- function(fit){
  mf <- model.frame(fit)
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Zf <- model.matrix(fit, component = "instruments")   # exógenas + IVs
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, NROW(X))
  storage.mode(X)  <- "double"
  storage.mode(Zf) <- "double"
  list(y=y, X=X, Zfull=Zf, w=as.numeric(w))
}

# --- util: remove colunas de variância ~0 ---
.drop_nzv <- function(A, tol = 1e-12){
  if (NCOL(A) == 0) return(A)
  keep <- apply(A, 2, function(x) sd(x) > tol)
  A[, keep, drop = FALSE]
}

# --- 2SLS ponderado em matriz (pré-multiplicação por sqrt(w)) ---
.w2sls <- function(y, X, Z, w, ridge = 1e-8){
  s  <- sqrt(w / mean(w))
  yS <- y * s; XS <- X * s; ZS <- Z * s
  
  ZtZ <- crossprod(ZS)
  if (!is.finite(rcond(ZtZ)) || rcond(ZtZ) < 1e-12)
    diag(ZtZ) <- diag(ZtZ) + ridge * mean(diag(ZtZ))
  
  PZ_X <- ZS %*% qr.solve(ZtZ, crossprod(ZS, XS))
  PZ_y <- ZS %*% qr.solve(ZtZ, crossprod(ZS, yS))
  
  XtPZ_X <- crossprod(XS, PZ_X)
  XtPZ_y <- crossprod(XS, PZ_y)
  if (!is.finite(rcond(XtPZ_X)) || rcond(XtPZ_X) < 1e-12)
    diag(XtPZ_X) <- diag(XtPZ_X) + ridge * mean(diag(XtPZ_X))
  
  beta <- solve(XtPZ_X, XtPZ_y)
  e    <- as.numeric(yS - XS %*% beta)
  list(beta = beta, e = e, X = XS, Z = ZS, n = NROW(XS), p = ncol(XS), k = ncol(ZS))
}

# --- Hansen-J robusto (HC0) a partir de e e ZS ---
.hansenJ_from_res <- function(e, ZS, p, ridge_seq = c(0, 1e-8, 1e-6, 1e-4)){
  n <- NROW(ZS)
  m <- crossprod(ZS, e) / n
  S <- crossprod(ZS * e) / n; S <- (S + t(S))/2
  note <- NA_character_; sol <- NULL; ok <- FALSE
  for (rg in ridge_seq){
    Sreg <- S
    if (rg > 0) diag(Sreg) <- diag(Sreg) + rg * mean(diag(Sreg))
    sol <- try(qr.solve(Sreg, m), silent = TRUE)
    if (!inherits(sol, "try-error") && all(is.finite(sol))) { ok <- TRUE; note <- if (rg>0) paste0("ridge=", rg) else NA_character_; break }
  }
  if (!ok) return(tibble(J = NA_real_, df = NA_real_, p = NA_real_, note = "S_singular"))
  J  <- drop(n * crossprod(m, sol))
  df <- ncol(ZS) - p
  tibble(J = if (df>0) J else NA_real_,
         df = if (df>0) df else 0L,
         p  = if (df>0) pchisq(J, df, lower.tail = FALSE) else NA_real_,
         note = note)
}

# --- QR: escolhe IVs independentes (mantém todas exógenas de X em Z) ---
J_qr_for_fit <- function(fit, tol = 1e-8){
  mats <- .get_mats(fit)
  X <- mats$X; Zf <- mats$Zfull
  Xn <- colnames(X)
  Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE]
  Ziv <- .drop_nzv(Ziv)
  if (NCOL(Ziv) > 0){
    Ziv_s <- scale(Ziv, center = TRUE, scale = TRUE)
    qrz   <- qr(Ziv_s)
    r     <- qrz$rank
    cols  <- if (r > 0) sort(qrz$pivot[seq_len(r)]) else integer(0)
    Zsel  <- if (length(cols)) Ziv[, cols, drop = FALSE] else Ziv[, 0, drop = FALSE]
  } else {
    Zsel <- Ziv
  }
  Z <- cbind(X, Zsel)
  est <- .w2sls(mats$y, X, Z, mats$w)
  out <- .hansenJ_from_res(est$e, est$Z, est$p)
  mutate(out, k_after = ncol(Z))
}

# --- PCA: condensa IVs em q PCs e testa ---
J_pca_for_fit <- function(fit, var_explained = 0.99, max_comp = 15){
  mats <- .get_mats(fit)
  X <- mats$X; Zf <- mats$Zfull; Xn <- colnames(X)
  Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE]
  Ziv <- .drop_nzv(Ziv)
  if (NCOL(Ziv) == 0) return(tibble(J = NA_real_, df = 0L, p = NA_real_, note = "no_iv", k_after = ncol(X)))
  pc <- prcomp(Ziv, center = TRUE, scale. = TRUE)
  varex <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  q <- min(which(varex >= var_explained)); q <- min(q, max_comp)
  Z <- cbind(X, pc$x[, seq_len(q), drop = FALSE])
  colnames(Z) <- c(colnames(X), paste0("PCiv", seq_len(q)))
  est <- .w2sls(mats$y, X, Z, mats$w)
  out <- .hansenJ_from_res(est$e, est$Z, est$p)
  mutate(out, k_after = ncol(Z), note = paste0("PCA(q=", q, ")"))
}

# ---- roda para todas as shares do seu objeto fits_iv ----
J_qr_tbl  <- map_dfr(share_vars, \(dep){
  fit <- fits_iv[[dep]]$fit
  res <- try(J_qr_for_fit(fit), silent = TRUE)
  if (inherits(res, "try-error")) tibble(share = dep, J = NA_real_, df = NA_real_, p = NA_real_, note = "erro_qr")
  else mutate(res, share = dep, .before = 1)
})

J_pca_tbl <- map_dfr(share_vars, \(dep){
  fit <- fits_iv[[dep]]$fit
  res <- try(J_pca_for_fit(fit, var_explained = 0.99, max_comp = 12), silent = TRUE)
  if (inherits(res, "try-error")) tibble(share = dep, J = NA_real_, df = NA_real_, p = NA_real_, note = "erro_pca")
  else mutate(res, share = dep, .before = 1)
})

list(J_qr = J_qr_tbl, J_pca = J_pca_tbl)

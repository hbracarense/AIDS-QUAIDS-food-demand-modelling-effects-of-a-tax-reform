# =========================================================
# INFERÊNCIA com IV: 2SLS (mínimo) + robustez
# Passos:
# (1) 2SLS com IVs mínimos (erros robustos) + diagnósticos (weak / Wu-Hausman / Sargan)
# (2) Hansen-J robusto (HC0) via matriz, com seleções QR e PCA (fallbacks estáveis)
# (3) Teste conjunto do bloco de preços (HC1)
# (4) Sensibilidade: 2SLS "full" vs "mínimo" (coeficientes dos preços)
# (5) Resumo executivo para reporte
# Requer objetos já existentes no ambiente: dat, share_vars, price_vars, weights_col,
#   iv_include_patterns_min, iv_include_patterns_full, iv_exclude_patterns,
#   z_var, z2_var, seasonal_controls_min, seasonal_controls_full,
#   e a função fit_iv_weighted(...) que você já vem usando.
# =========================================================

suppressPackageStartupMessages({
  library(AER)
  library(sandwich)
  library(car)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(broom)
})

# -------------------------
# Helpers genéricos
# -------------------------
# -------- 1) Diagnósticos "padrão" do summary(ivreg), sem quebrar quando vier vazio --------
extract_iv_diags <- function(fits_iv){
  out <- dplyr::bind_rows(lapply(names(fits_iv), function(v){
    dmat <- try(fits_iv[[v]]$summary$diagnostics, silent = TRUE)
    if (inherits(dmat, "try-error") || is.null(dmat)) return(tibble::tibble())
    dd <- as.data.frame(dmat)
    dd$teste <- rownames(dmat)
    dd$share <- v
    rownames(dd) <- NULL
    tibble::as_tibble(dd)
  }))
  if (nrow(out) == 0) {
    # retorna tibble vazio com o "esqueleto" certo
    tibble::tibble(
      share = character(),
      teste = character(),
      statistic = double(),
      df1 = double(),
      df2 = double(),
      `p-value` = double()
    )
  } else {
    dplyr::relocate(out, share, teste)
  }
}

iv_diags_min <- extract_iv_diags(fits_iv_min)
iv_diags_min

# --- util para Hansen-J robusto (HC0) sem ivpack/gmm
.get_mats <- function(fit){
  mf <- model.frame(fit)
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Zf <- model.matrix(fit, component = "instruments") # exógenas + IVs
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, NROW(X))
  storage.mode(X)  <- "double"
  storage.mode(Zf) <- "double"
  list(y=y, X=X, Zfull=Zf, w=as.numeric(w))
}
.drop_nzv <- function(A, tol = 1e-12){
  if (is.null(dim(A)) || NCOL(A) == 0) return(A)
  keep <- apply(A, 2, function(x) sd(x) > tol)
  A[, keep, drop = FALSE]
}
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
J_qr_for_fit <- function(fit, tol = 1e-8){
  mats <- .get_mats(fit)
  X <- mats$X; Zf <- mats$Zfull
  Xn <- colnames(X)
  Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE]
  Ziv <- .drop_nzv(Ziv)
  if (NCOL(Ziv) > 0){
    Ziv_s <- scale(Ziv, center = TRUE, scale = TRUE)
    qrz   <- qr(Ziv_s); r <- qrz$rank
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
J_pca_for_fit <- function(fit, var_explained = 0.99, max_comp = 12){
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

# =========================================================
# (1) 2SLS com IVs MÍNIMOS (erros robustos) + diagnósticos padrão
# =========================================================
# (1) 2SLS com IVs MÍNIMOS (erros robustos) + diagnósticos padrão
fits_iv_min <- purrr::map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x,
  price_vars = price_vars,
  weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_min,
  iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var,           # pode manter z2_var mesmo que não use
  z2_var = z2_var,         # a função simplesmente ignorará se não estiver implementado
  seasonal_controls = seasonal_controls_min,
  use_spline = FALSE,
  spline_df = 4
))
names(fits_iv_min) <- share_vars

iv_diags_min <- extract_iv_diags(fits_iv_min)
cat("\n==== (1) DIAGNÓSTICOS 2SLS (IVs mínimos) ====\n")
print(iv_diags_min)

# Tabela de coeficientes dos preços com SE robusto (HC1)
coef_iv_min <- map_dfr(share_vars, function(v){
  fit <- fits_iv_min[[v]]$fit
  tt  <- tidy(fit, conf.int = TRUE, vcov. = sandwich::vcovHC(fit, type = "HC1")) |>
    filter(grepl("^ln_preco_com_reforma2", term)) |>
    mutate(share = v, .before = 1)
  tt
})
cat("\n==== (1b) COEFICIENTES (preços) – 2SLS mínimo (HC1) ====\n")
print(coef_iv_min)

# =========================================================
# (2) Hansen-J ROBUSTO (HC0) via QR e PCA (sem ivpack/gmm)
#     -> usa a mesma especificação mínima do passo (1)
# =========================================================
J_qr_tbl  <- map_dfr(share_vars, \(dep){
  fit <- fits_iv_min[[dep]]$fit
  res <- try(J_qr_for_fit(fit), silent = TRUE)
  if (inherits(res, "try-error")) tibble(share = dep, J = NA_real_, df = NA_real_, p = NA_real_, note = "erro_qr", k_after = NA_integer_)
  else mutate(res, share = dep, .before = 1)
})
J_pca_tbl <- map_dfr(share_vars, \(dep){
  fit <- fits_iv_min[[dep]]$fit
  res <- try(J_pca_for_fit(fit, var_explained = 0.99, max_comp = 12), silent = TRUE)
  if (inherits(res, "try-error")) tibble(share = dep, J = NA_real_, df = NA_real_, p = NA_real_, note = "erro_pca", k_after = NA_integer_)
  else mutate(res, share = dep, .before = 1)
})
cat("\n==== (2) HANSEN-J robusto (QR) – 2SLS mínimo ====\n")
print(J_qr_tbl)
cat("\n==== (2) HANSEN-J robusto (PCA) – 2SLS mínimo ====\n")
print(J_pca_tbl)

# =========================================================
# (3) TESTE CONJUNTO: bloco de preços = 0 (2SLS mínimo, HC1)
# =========================================================
joint_price_tests_min <- bind_rows(lapply(share_vars, function(v) {
  fit <- fits_iv_min[[v]]$fit
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
cat("\n==== (3) TESTE CONJUNTO (preços) – 2SLS mínimo (HC1) ====\n")
print(joint_price_tests_min)

# =========================================================
# (4) SENSIBILIDADE: 2SLS "FULL" vs "MÍNIMO" (coeficientes de preços)
# =========================================================
fits_iv_full <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv_full) <- share_vars

coef_iv_full <- map_dfr(share_vars, function(v){
  fit <- fits_iv_full[[v]]$fit
  tidy(fit, conf.int = TRUE, vcov. = sandwich::vcovHC(fit, type = "HC1")) |>
    filter(grepl("^ln_preco_com_reforma2", term)) |>
    transmute(share = v, term, estimate_full = estimate, conf.low_full = conf.low, conf.high_full = conf.high)
})

coef_compare <- coef_iv_min |>
  select(share, term, estimate_min = estimate, conf.low_min = conf.low, conf.high_min = conf.high) |>
  left_join(coef_iv_full, by = c("share","term")) |>
  mutate(
    diff = estimate_min - estimate_full,
    same_sign = sign(estimate_min) == sign(estimate_full)
  )

cat("\n==== (4) SENSIBILIDADE: 2SLS FULL vs MÍNIMO (preços, HC1) ====\n")
print(coef_compare)

stab_by_share <- coef_compare |>
  group_by(share) |>
  summarise(
    n_terms = n(),
    n_same_sign = sum(same_sign, na.rm = TRUE),
    frac_same_sign = n_same_sign / n_terms,
    max_abs_diff = max(abs(diff), na.rm = TRUE),
    .groups = "drop"
  )
cat("\n---- Estabilidade por share (sinais e dif. máximas) ----\n")
print(stab_by_share)

# =========================================================
# (5) RESUMO EXECUTIVO para INFERÊNCIA
#    - Use 2SLS com IVs mínimos para inferir elasticidades
#    - Reporte HC1, teste conjunto do bloco de preços,
#      e os Hansen-J robustos (QR/PCA) como evidência de sobre-ID
# =========================================================
cat("\n==== (5) RESUMO ====\n")
cat("* Especificação escolhida para inferência: 2SLS com IVs mínimos + HC1.\n")
cat("* Reporte: (i) coeficientes e ICs (tabela 'coef_iv_min');\n")
cat("           (ii) teste conjunto do bloco de preços (tabela 'joint_price_tests_min');\n")
cat("           (iii) Hansen-J robusto QR/PCA (tabelas 'J_qr_tbl' e 'J_pca_tbl');\n")
cat("           (iv) sensibilidade vs IVs FULL (tabelas 'coef_compare' e 'stab_by_share').\n")

## ============================
## Setup
## ============================
suppressPackageStartupMessages({
  library(AER)
  library(dplyr); library(purrr); library(tibble); library(broom)
  library(sandwich); library(car)
})

# >>>> PARAMS <<<<
cluster_var <- NULL   # ex: "id_regiao" ou c("id_regiao","mes")
alpha <- 0.05

## ============================
## Helpers já compatíveis com seus objetos
## ============================
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


.get_mats <- function(fit){
  mf <- model.frame(fit)
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Zf <- model.matrix(fit, component = "instruments")
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, NROW(X))
  storage.mode(X)  <- "double"; storage.mode(Zf) <- "double"
  list(y=y, X=X, Zfull=Zf, w=as.numeric(w))
}
.drop_nzv <- function(A, tol = 1e-12){
  if (is.null(dim(A)) || NCOL(A) == 0) return(A)
  keep <- apply(A, 2, function(x) sd(x) > tol)
  A[, keep, drop = FALSE]
}
.w2sls <- function(y, X, Z, w, ridge = 1e-8){
  s  <- sqrt(w / mean(w)); yS <- y * s; XS <- X * s; ZS <- Z * s
  ZtZ <- crossprod(ZS); if (!is.finite(rcond(ZtZ)) || rcond(ZtZ) < 1e-12)
    diag(ZtZ) <- diag(ZtZ) + ridge * mean(diag(ZtZ))
  PZ_X <- ZS %*% qr.solve(ZtZ, crossprod(ZS, XS))
  PZ_y <- ZS %*% qr.solve(ZtZ, crossprod(ZS, yS))
  XtPZ_X <- crossprod(XS, PZ_X); XtPZ_y <- crossprod(XS, PZ_y)
  if (!is.finite(rcond(XtPZ_X)) || rcond(XtPZ_X) < 1e-12)
    diag(XtPZ_X) <- diag(XtPZ_X) + ridge * mean(diag(XtPZ_X))
  beta <- solve(XtPZ_X, XtPZ_y); e <- as.numeric(yS - XS %*% beta)
  list(beta = beta, e = e, X = XS, Z = ZS, n = NROW(XS), p = ncol(XS), k = ncol(ZS))
}
.hansenJ_from_res <- function(e, ZS, p, ridge_seq = c(0, 1e-8, 1e-6, 1e-4)){
  n <- NROW(ZS); m <- crossprod(ZS, e) / n
  S <- crossprod(ZS * e) / n; S <- (S + t(S))/2
  note <- NA_character_; sol <- NULL; ok <- FALSE
  for (rg in ridge_seq){
    Sreg <- S; if (rg > 0) diag(Sreg) <- diag(Sreg) + rg * mean(diag(Sreg))
    sol <- try(qr.solve(Sreg, m), silent = TRUE)
    if (!inherits(sol, "try-error") && all(is.finite(sol))) { ok <- TRUE; note <- if (rg>0) paste0("ridge=", rg) else NA_character_; break }
  }
  if (!ok) return(tibble(J = NA_real_, df = NA_real_, p = NA_real_, note = "S_singular"))
  J  <- drop(n * crossprod(m, sol)); df <- ncol(ZS) - p
  tibble(J = if (df>0) J else NA_real_, df = if (df>0) df else 0L,
         p  = if (df>0) pchisq(J, df, lower.tail = FALSE) else NA_real_, note = note)
}
J_qr_for_fit <- function(fit, tol = 1e-8){
  mats <- .get_mats(fit); X <- mats$X; Zf <- mats$Zfull; Xn <- colnames(X)
  Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE] |> .drop_nzv()
  if (NCOL(Ziv) > 0){
    Ziv_s <- scale(Ziv, center = TRUE, scale = TRUE)
    qrz <- qr(Ziv_s); r <- qrz$rank; cols <- if (r > 0) sort(qrz$pivot[seq_len(r)]) else integer(0)
    Zsel <- if (length(cols)) Ziv[, cols, drop = FALSE] else Ziv[, 0, drop = FALSE]
  } else Zsel <- Ziv
  Z <- cbind(X, Zsel)
  est <- .w2sls(mats$y, X, Z, mats$w)
  out <- .hansenJ_from_res(est$e, est$Z, est$p)
  mutate(out, k_after = ncol(Z))
}
J_pca_for_fit <- function(fit, var_explained = 0.99, max_comp = 12){
  mats <- .get_mats(fit); X <- mats$X; Zf <- mats$Zfull; Xn <- colnames(X)
  Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE] |> .drop_nzv()
  if (NCOL(Ziv) == 0) return(tibble(J = NA_real_, df = 0L, p = NA_real_, note = "no_iv", k_after = ncol(X)))
  pc <- prcomp(Ziv, center = TRUE, scale. = TRUE)
  varex <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  q <- min(which(varex >= var_explained)); q <- min(q, max_comp)
  Z <- cbind(X, pc$x[, seq_len(q), drop = FALSE]); colnames(Z) <- c(colnames(X), paste0("PCiv", seq_len(q)))
  est <- .w2sls(mats$y, X, Z, mats$w)
  out <- .hansenJ_from_res(est$e, est$Z, est$p)
  mutate(out, k_after = ncol(Z), note = paste0("PCA(q=", q, ")"))
}

## vcov robusto (HC1 ou cluster, se informado)
.vcov_robust <- function(fit, data, cluster_var = NULL){
  if (is.null(cluster_var)) return(vcovHC(fit, type = "HC1"))
  cl <- data[, cluster_var, drop = FALSE]
  if (ncol(cl) == 1L) return(vcovCL(fit, cluster = cl[[1]]))
  # multiway cluster (adição “à la Cameron-Gelbach-Miller” simples)
  V <- Reduce(`+`, lapply(cl, function(ci) vcovCL(fit, cluster = ci))) -
    (length(cl)-1) * vcovHC(fit, type = "HC1")
  V
}

## ============================
## (1) 2SLS MÍNIMO + HC1/cluster + diagnósticos
## ============================
fits_iv_min <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var, seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv_min) <- share_vars

iv_diags_min <- extract_iv_diags(fits_iv_min); print(iv_diags_min)

coef_iv_min <- map_dfr(share_vars, function(v){
  fit <- fits_iv_min[[v]]$fit
  V   <- .vcov_robust(fit, dat, cluster_var)
  tidy(fit, conf.int = TRUE, vcov. = V) |>
    filter(grepl("^ln_preco_com_reforma2", term)) |>
    mutate(share = v, .before = 1)
})
print(coef_iv_min)

## ============================
## (2) Hansen-J ROBUSTO (QR e PCA)
## ============================
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
print(J_qr_tbl); print(J_pca_tbl)

## ============================
## (3) TESTE CONJUNTO bloco de preços (HC1/cluster)
## ============================
joint_price_tests_min <- bind_rows(lapply(share_vars, function(v) {
  fit <- fits_iv_min[[v]]$fit; V <- .vcov_robust(fit, dat, cluster_var)
  cn  <- names(coef(fit))
  price_cn <- cn[grepl("^ln_preco_com_reforma2", cn)]
  if (length(price_cn) == 0) return(tibble(share = v, F = NA, df = NA, p = NA))
  hyp <- paste0(price_cn, " = 0")
  lh  <- car::linearHypothesis(fit, hyp, vcov. = V)
  tibble(share = v,
         F  = as.numeric(lh[2, "F"]),
         df = paste0(lh[2, "Df"], collapse = ","),
         p  = as.numeric(lh[2, "Pr(>F)"]))
}))
print(joint_price_tests_min)

## ============================
## (4) SENSIBILIDADE: FULL vs MÍNIMO (coef. de preços)
## ============================
fits_iv_full <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var, seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv_full) <- share_vars

coef_iv_full <- map_dfr(share_vars, function(v){
  fit <- fits_iv_full[[v]]$fit; V <- .vcov_robust(fit, dat, cluster_var)
  tidy(fit, conf.int = TRUE, vcov. = V) |>
    filter(grepl("^ln_preco_com_reforma2", term)) |>
    transmute(share = v, term, estimate_full = estimate, conf.low_full = conf.low, conf.high_full = conf.high)
})

coef_compare <- coef_iv_min |>
  select(share, term, estimate_min = estimate, conf.low_min = conf.low, conf.high_min = conf.high) |>
  left_join(coef_iv_full, by = c("share","term")) |>
  mutate(diff = estimate_min - estimate_full,
         same_sign = sign(estimate_min) == sign(estimate_full))
print(coef_compare)

stab_by_share <- coef_compare |>
  group_by(share) |>
  summarise(n_terms = n(),
            n_same_sign = sum(same_sign, na.rm = TRUE),
            frac_same_sign = n_same_sign / n_terms,
            max_abs_diff = max(abs(diff), na.rm = TRUE), .groups = "drop")
print(stab_by_share)

## ============================
## (5) OPCIONAL – Robustez a IV fraco (LIML/Fuller, AR/CLR)
## ============================
# LIML/Fuller (se pacote 'ivreg' novo estiver instalado)
if (requireNamespace("ivreg", quietly = TRUE)) {
  fits_liml_min <- map(share_vars, function(v){
    f <- fits_iv_min[[v]]$fit
    # Reusa a mesma fórmula de f
    fm <- formula(f)
    ivreg::ivreg(fm, data = dat, weights = model.weights(model.frame(f)), method = "LIML")
    # Para Fuller, use: method = "LIML", kappa = 1 (ou 4)
  })
  names(fits_liml_min) <- share_vars
  
  coef_liml_min <- map_dfr(share_vars, function(v){
    fit <- fits_liml_min[[v]]
    V   <- .vcov_robust(fit, dat, cluster_var)
    broom::tidy(fit, conf.int = TRUE, vcov. = V) |>
      filter(grepl("^ln_preco_com_reforma2", term)) |>
      mutate(share = v, .before = 1)
  })
  print(coef_liml_min)
} else {
  message("Opcional LIML/Fuller: instale o pacote 'ivreg' (CRAN) para rodar este bloco.")
}

# Anderson–Rubin/CLR (opcional, se 'ivmodel' disponível)
if (requireNamespace("ivmodel", quietly = TRUE)) {
  # Exemplo raso (o pacote exige especificação de Y, D, Z, X)
  message("Use ivmodel::ivmodel() para AR/CLR no seu share-chave; montar Y,D,Z,X conforme seu caso.")
} else {
  message("Opcional AR/CLR: instale 'ivmodel' para testes robustos a IV fraco.")
}

## ============================
## Saída executiva (para paper)
## ============================
cat("\n==== RESUMO PARA INFERÊNCIA ====\n")
cat("* Especificação base: 2SLS com IVs mínimos + erros HC1", ifelse(is.null(cluster_var),"", paste0(" (cluster: ", paste(cluster_var, collapse = ", "), ")")), ".\n", sep = "")
cat("* Reporte: coeficientes e ICs (coef_iv_min), teste conjunto de preços (joint_price_tests_min),\n")
cat("          Hansen-J robusto QR/PCA (J_qr_tbl, J_pca_tbl),\n")
cat("          e sensibilidade vs FULL (coef_compare, stab_by_share).\n")
cat("          (Opcional) LIML/Fuller e AR/CLR se pacotes disponíveis.\n")

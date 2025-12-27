## =========================================================
## INFERÊNCIA com IV – Roteiro 1–5 (IC via Bootstrap Poisson+1 cluster UF)
## =========================================================
suppressPackageStartupMessages({
  library(AER)
  library(dplyr); library(purrr); library(tibble); library(broom)
  library(sandwich); library(car)
})

## -----------------------------
## PARAMS (ajuste se quiser)
## -----------------------------
cluster_var <- NULL     # ex: "id_regiao" ou c("id_regiao","mes"); NULL = HC1 (fora do bootstrap)
alpha <- 0.05
zcrit <- qnorm(1 - alpha/2)

## ---- PARÂMETROS DO BOOTSTRAP (Poisson+1) ----
cluster_var_boot <- "UF"  # variável de cluster para o bootstrap (troque se necessário)
boot_B    <- 2000          # número de réplicas
winsor_p  <- 0.02         # winsorização bilateral (1% e 99%) sobre os draws
abs_cap   <- 10            # cap absoluto: |e| <= abs_cap
boot_seed <- 123

## =========================================================
## Helpers genéricos (compatíveis com seus objetos)
## =========================================================

## (A) Diagnósticos padrão do AER::ivreg
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

## (B) VCOV robusto: HC1 ou (multiway) cluster
.vcov_robust <- function(fit, data, cluster_var = NULL){
  if (is.null(cluster_var)) return(vcovHC(fit, type = "HC1"))
  cl <- data[, cluster_var, drop = FALSE]
  if (ncol(cl) == 1L) return(vcovCL(fit, cluster = cl[[1]]))
  V <- Reduce(`+`, lapply(cl, function(ci) vcovCL(fit, cluster = ci))) -
    (length(cl)-1) * vcovHC(fit, type = "HC1")
  V
}

## (C) Hansen-J robusto (HC0) em matriz – estável (sem ivpack/gmm)
.get_mats <- function(fit){
  mf <- model.frame(fit)
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Zf <- model.matrix(fit, component = "instruments")  # exógenas + IVs
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
    qrz <- qr(Ziv_s); r <- qrz$rank
    cols <- if (r > 0) sort(qrz$pivot[seq_len(r)]) else integer(0)
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

## (D) Pesos e shares médios
w_bar <- setNames(
  sapply(share_vars, function(v) weighted.mean(dat[[v]], dat[[weights_col]])),
  share_vars
)

## =========================================================
## (1) 2SLS MÍNIMO + diagnósticos + tabela de coef. (HC1/cluster)
## =========================================================
fits_iv_min <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var, seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = 4
))
names(fits_iv_min) <- share_vars

iv_diags_min <- extract_iv_diags(fits_iv_min)

coef_iv_min <- map_dfr(share_vars, function(v){
  fit <- fits_iv_min[[v]]$fit
  V   <- .vcov_robust(fit, dat, cluster_var)
  tidy(fit, conf.int = TRUE, vcov. = V) |>
    filter(grepl("^ln_preco_com_reforma2", term)) |>
    mutate(share = v, .before = 1)
})

cat("\n==== (1) DIAGNÓSTICOS 2SLS (mínimo) ====\n"); print(iv_diags_min)
cat("\n==== (1b) COEFICIENTES DE PREÇO – 2SLS (mínimo, robusto) ====\n"); print(coef_iv_min, n = 36)

## =========================================================
## (2) HANSEN-J ROBUSTO (QR & PCA) – 2SLS (mínimo)
## =========================================================
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

cat("\n==== (2) Hansen-J (QR) – 2SLS mínimo ====\n");  print(J_qr_tbl)
cat("\n==== (2) Hansen-J (PCA) – 2SLS mínimo ====\n"); print(J_pca_tbl)

## =========================================================
## (3) TESTE CONJUNTO – bloco de preços = 0 (HC1/cluster)
## =========================================================
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
cat("\n==== (3) TESTE CONJUNTO (preços) – 2SLS mínimo ====\n"); print(joint_price_tests_min)

## =========================================================
## (4) SENSIBILIDADE: 2SLS FULL vs MÍNIMO (coef. de preços)
## =========================================================
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
    transmute(share = v, term, estimate_full = estimate,
              conf.low_full = conf.low, conf.high_full = conf.high)
})

coef_compare <- coef_iv_min |>
  select(share, term, estimate_min = estimate, conf.low_min = conf.low, conf.high_min = conf.high) |>
  left_join(coef_iv_full, by = c("share","term")) |>
  mutate(diff = estimate_min - estimate_full,
         same_sign = sign(estimate_min) == sign(estimate_full))

stab_by_share <- coef_compare |>
  group_by(share) |>
  summarise(n_terms = n(),
            n_same_sign = sum(same_sign, na.rm = TRUE),
            frac_same_sign = n_same_sign / n_terms,
            max_abs_diff = max(abs(diff), na.rm = TRUE), .groups = "drop")

cat("\n==== (4) SENSIBILIDADE – FULL vs MÍNIMO ====\n"); print(coef_compare, n = 36)
cat("\n---- Estabilidade por share ----\n"); print(stab_by_share)

## =========================================================
## (5) ELASTICIDADES (AIDS) – ponto + IC via Bootstrap Poisson+1 (cluster UF)
##     Marshall:  e_ij = -1{i=j} + (γ_ij - β_i * w̄_j) / w̄_i
##     Hicks:     e^c_ij = -1{i=j} + γ_ij/w̄_i + w̄_j
##     Renda:     η_i  = 1 + β_i / w̄_i
##  >>> Substitui IC delta por IC percentil (winsor + cap) <<<
## =========================================================

## ---- helpers para elasticidades (ponto) ----
.elasticities_point_one <- function(fit, i_idx, share_name, price_vars, z_var, w_bar_vec){
  cn <- names(coef(fit))
  w_i <- w_bar_vec[share_name]
  stopifnot(is.finite(w_i), w_i > 0)
  
  beta_nm <- z_var
  beta_ok <- beta_nm %in% cn
  beta_i  <- if (beta_ok) coef(fit)[beta_nm] else NA_real_
  
  PEs <- map_dfr(seq_along(price_vars), function(j){
    pj <- price_vars[j]
    w_j <- w_bar_vec[j]
    
    gamma_ok <- pj %in% cn
    gamma_ij <- if (gamma_ok) coef(fit)[pj] else NA_real_
    
    e_m <- if (gamma_ok && beta_ok) (-as.numeric(i_idx==j) + (gamma_ij - beta_i*w_j)/w_i) else NA_real_
    e_h <- if (gamma_ok) (-as.numeric(i_idx==j) + gamma_ij/w_i + w_j) else NA_real_
    
    tibble(
      share = share_name, price = pj, own = (i_idx==j),
      e_mar = e_m, e_hix = e_h
    )
  })
  
  eta  <- if (beta_ok) 1 + beta_i / w_i else NA_real_
  list(PEs = PEs, YEl = tibble(share = share_name, eta = eta))
}

## ---- winsor + cap para um vetor de draws ----
.apply_winsor_cap <- function(x, p = winsor_p, cap = abs_cap){
  x <- pmax(pmin(x, cap), -cap)  # cap absoluto
  if (is.finite(p) && p > 0 && sum(is.finite(x)) >= 5){
    qs <- quantile(x, c(p, 1 - p), na.rm = TRUE, type = 7)
    x <- pmin(pmax(x, qs[1], na.rm = TRUE), qs[2], na.rm = TRUE)
  }
  x
}

## ---- ponto (amostra original) ----
elas_lists_point <- map2(seq_along(share_vars), share_vars, function(i, si){
  fit <- fits_iv_min[[si]]$fit
  .elasticities_point_one(fit, i, si, price_vars, z_var, w_bar)
})
elas_point_marhix <- bind_rows(lapply(elas_lists_point, `[[`, "PEs"))
elas_point_income <- bind_rows(lapply(elas_lists_point, `[[`, "YEl"))

## ---- Bootstrap Poisson+1 clusterizado por UF ----
set.seed(boot_seed)

## gera pesos Poisson(1)+1 por cluster e retorna coluna 'w_eff'
.gen_boot_weights <- function(data, cl_var, base_w_col){
  if (!is.null(cl_var) && cl_var %in% names(data)){
    g <- droplevels(as.factor(data[[cl_var]]))
    lam <- setNames(rpois(nlevels(g), 1) + 1, levels(g))  # Poisson+1
    w_boot <- lam[as.character(g)]
  } else {
    w_boot <- rpois(nrow(data), 1) + 1                     # Poisson+1 por linha
  }
  w_eff <- w_boot * data[[base_w_col]]
  w_eff
}

## calcula elasticidades (ponto) numa réplica bootstrap com pesos Poisson+1
.elasticities_point_weighted <- function(data0, base_w_col, cl_var){
  dataB <- data0
  dataB$w_eff <- .gen_boot_weights(dataB, cl_var, base_w_col)
  
  ## w_bar (amostra bootstrap) com peso efetivo
  w_bar_B <- setNames(
    sapply(share_vars, function(v) weighted.mean(dataB[[v]], dataB$w_eff)),
    share_vars
  )
  
  ## reestima fits mínimos usando w_eff como coluna de pesos
  fitsB <- map(share_vars, ~ fit_iv_weighted(
    dat = dataB, dep = .x, price_vars = price_vars, weights_col = "w_eff",
    iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
    z_var = z_var, z2_var = z2_var, seasonal_controls = seasonal_controls_min,
    use_spline = FALSE, spline_df = 4
  ))
  names(fitsB) <- share_vars
  
  ## elasticidades ponto
  out_list <- map2(seq_along(share_vars), share_vars, function(i, si){
    fit <- fitsB[[si]]$fit
    .elasticities_point_one(fit, i, si, price_vars, z_var, w_bar_B)
  })
  mh <- bind_rows(lapply(out_list, `[[`, "PEs"))
  yl <- bind_rows(lapply(out_list, `[[`, "YEl"))
  list(mh = mh, yl = yl)
}

## roda o bootstrap e empilha resultados
boot_draws <- replicate(boot_B, {
  res <- try(.elasticities_point_weighted(dat, weights_col, cluster_var_boot), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res)) {
    list(mh = tibble(), yl = tibble())
  } else res
}, simplify = FALSE)

boot_mh_all <- bind_rows(lapply(seq_along(boot_draws), function(b){
  if (nrow(boot_draws[[b]]$mh) == 0) return(tibble())
  mutate(boot_draws[[b]]$mh, .b = b)
}))
boot_yl_all <- bind_rows(lapply(seq_along(boot_draws), function(b){
  if (nrow(boot_draws[[b]]$yl) == 0) return(tibble())
  mutate(boot_draws[[b]]$yl, .b = b)
}))

## sumariza por (share, price) para Marshall/Hicks
summ_mh <- boot_mh_all |>
  group_by(share, price, own) |>
  summarise(
    ## Marshall
    lo_mar = {v <- .apply_winsor_cap(e_mar); quantile(v, alpha/2, na.rm = TRUE, type = 7)},
    hi_mar = {v <- .apply_winsor_cap(e_mar); quantile(v, 1 - alpha/2, na.rm = TRUE, type = 7)},
    se_mar_boot = {v <- .apply_winsor_cap(e_mar); sd(v, na.rm = TRUE)},
    ## Hicks
    lo_hix = {v <- .apply_winsor_cap(e_hix); quantile(v, alpha/2, na.rm = TRUE, type = 7)},
    hi_hix = {v <- .apply_winsor_cap(e_hix); quantile(v, 1 - alpha/2, na.rm = TRUE, type = 7)},
    se_hix_boot = {v <- .apply_winsor_cap(e_hix); sd(v, na.rm = TRUE)},
    .groups = "drop"
  )

## sumariza por (share) para renda
summ_yl <- boot_yl_all |>
  group_by(share) |>
  summarise(
    lo_eta = {v <- .apply_winsor_cap(eta); quantile(v, alpha/2, na.rm = TRUE, type = 7)},
    hi_eta = {v <- .apply_winsor_cap(eta); quantile(v, 1 - alpha/2, na.rm = TRUE, type = 7)},
    se_eta_boot = {v <- .apply_winsor_cap(eta); sd(v, na.rm = TRUE)},
    .groups = "drop"
  )

## junta ponto (amostra original) + IC bootstrap
elas_marshall <- elas_point_marhix |>
  left_join(summ_mh, by = c("share","price","own")) |>
  relocate(share, price, own, e_mar, lo_mar, hi_mar, se_mar_boot)

elas_hicks <- elas_point_marhix |>
  left_join(summ_mh, by = c("share","price","own")) |>
  transmute(share, price, own,
            e_hix = e_hix,
            lo_hix = lo_hix, hi_hix = hi_hix, se_hix_boot = se_hix_boot)

elas_income <- elas_point_income |>
  left_join(summ_yl, by = "share") |>
  relocate(share, eta, lo_eta, hi_eta, se_eta_boot)

## prints finais
cat("\n==== (5) ELASTICIDADES – Marshall (ponto + IC bootstrap Poisson+1 UF) ====\n")
print(elas_marshall, n = nrow(elas_marshall))

cat("\n==== (5) ELASTICIDADES – Hicks (ponto + IC bootstrap Poisson+1 UF) ====\n")
print(elas_hicks, n = nrow(elas_hicks))

cat("\n==== (5) ELASTICIDADE-RENDA (ponto + IC bootstrap Poisson+1 UF) ====\n")
print(elas_income, n = nrow(elas_income))

## Objetos principais criados:
## - iv_diags_min, coef_iv_min
## - J_qr_tbl, J_pca_tbl
## - joint_price_tests_min
## - coef_compare, stab_by_share
## - elas_marshall, elas_hicks, elas_income   (com ICs por bootstrap Poisson+1 + winsor + cap)

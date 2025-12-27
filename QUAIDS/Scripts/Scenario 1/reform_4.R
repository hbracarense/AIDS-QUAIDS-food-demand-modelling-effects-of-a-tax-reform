## =========================================================
##  BLOCO ADICIONAL — IVs por bem, F_cond e AR/CLR por preço
## =========================================================

## 0) Helpers (não mexa)
`%||%` <- function(a,b) if (is.null(a)) b else a

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Faltam pacotes: ", paste(miss, collapse=", "))
}
need_pkg(c("AER","dplyr","purrr","tibble"))

has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)

first_stage_F <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if (length(X)) X else "1"
  rhs1 <- c(X, Z); if (!length(rhs1)) rhs1 <- "1"
  f0 <- lm(stats::reformulate(rhs0, response = D), data = data)
  f1 <- lm(stats::reformulate(rhs1, response = D), data = data)
  as.numeric(anova(f0, f1)$F[2])
}
partial_R2 <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if (length(X)) X else "1"
  rhs1 <- c(X, Z); if (!length(rhs1)) rhs1 <- "1"
  f0 <- lm(stats::reformulate(rhs0, response = D), data = data)
  f1 <- lm(stats::reformulate(rhs1, response = D), data = data)
  1 - sum(residuals(f1)^2)/sum(residuals(f0)^2)
}

## 1) Gera IVs aumentados (interações com shifters binários)
##    Ex.: with = c("SH_SH_area_1","SH_SH_capital_1"). Se não houver, passa reto.
augment_ivs <- function(df, base_iv, with = c("SH_SH_area_1","SH_SH_capital_1")) {
  with <- intersect(with, names(df))
  out_names <- base_iv
  if (!length(with)) return(list(data=df, iv_names=unique(out_names)))
  for (b in base_iv) if (b %in% names(df)) {
    for (w in with) {
      nm <- paste0(b,"_x_", w)
      if (!nm %in% names(df)) df[[nm]] <- as.numeric(df[[b]]) * as.numeric(df[[w]])
      out_names <- c(out_names, nm)
    }
  }
  list(data=df, iv_names=unique(out_names))
}

## 2) Prepara (X,Z) condicionais corretos para um preço de interesse
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  stopifnot(poi %in% names(df), eq_y %in% names(df))
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  xnames <- c(intersect(c("z","z2"), names(df)), intersect(ln_others, names(df)))
  znames <- intersect(iv_names, names(df))
  use <- unique(c(eq_y, poi, xnames, znames))
  dat <- df[stats::complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  list(dat=dat, 
       Y=as.numeric(dat[[eq_y]]), 
       D=as.numeric(dat[[poi]]),
       X=if(length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL,
       Z=if(length(znames)) as.matrix(dat[, znames, drop = FALSE]) else NULL,
       Xnames=xnames, Znames=znames)
}

## 3) AR/CLR via grid adaptativo (usa ivmodel se disponível; senão retorna NA)
ar_clr_ci_grid_ivmodel <- function(Y, D, Z, X = NULL, alpha = 0.05, grid_center = NULL, halfwidth = 50) {
  if (!has_ivmodel) return(list(ci_ar=c(NA_real_,NA_real_), ci_clr=c(NA_real_,NA_real_),
                                center=NA_real_, halfwidth=NA_real_))
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  if (is.null(grid_center)) {
    b2 <- tryCatch({ unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]) }, error=function(e) NA_real_)
    if (!is.finite(b2)) b2 <- 0
  } else b2 <- grid_center
  grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out = 4001)
  ivm  <- ivmodel::ivmodel(Y=Y, D=D, Z=Z, X=X)
  p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0=b0)$p.value, numeric(1))
  p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0=b0)$p.value, numeric(1))
  keep_ar  <- which(p_ar  >= alpha)
  keep_clr <- which(p_clr >= alpha)
  ci_ar  <- if(length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
  ci_clr <- if(length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
  list(ci_ar=ci_ar, ci_clr=ci_clr, center=b2, halfwidth=halfwidth)
}

## 4) F_cond + AR/CLR para TODOS os preços endógenos (j*)
run_iv_suite_all_prices <- function(df, priceNames, shareNames, omit_share, drop_price, base_iv_set,
                                    shifters_for_iv = c("SH_SH_area_1","SH_SH_capital_1"),
                                    exogs = c("z","z2"),
                                    alpha = 0.05, halfwidth = 50) {
  # IVs aumentados
  aug <- augment_ivs(df, base_iv_set, with = shifters_for_iv)
  dfA <- aug$data
  ivA <- aug$iv_names
  
  # listas
  endogs <- paste0("ln_", setdiff(priceNames, drop_price))
  eqs    <- setdiff(shareNames, omit_share)
  exogs  <- intersect(exogs, names(dfA))
  
  # 4.1 — Tabela de 1º estágio (F_cond, R² parcial) por preço
  sw_tbl <- purrr::map_dfr(endogs, function(Dv){
    use <- unique(c(Dv, exogs, ivA))
    cc  <- stats::complete.cases(dfA[, use, drop = FALSE])
    dat <- dfA[cc, , drop = FALSE]
    tibble::tibble(
      price = Dv,
      n_cc  = nrow(dat),
      F_SW_cond  = first_stage_F(Dv, exogs, ivA, dat),
      R2_partial = partial_R2(Dv, exogs, ivA, dat),
      df1 = length(ivA),
      df2 = nrow(dat) - length(c(exogs, ivA)) - 1
    )
  })
  
  # 4.2 — AR/CLR por (preço j*, equação i)
  arclr <- purrr::map_dfr(endogs, function(poi){
    purrr::map_dfr(eqs, function(eq_y){
      pr <- prep_XZ_cond(eq_y, poi, dfA, ivA, drop_price, priceNames)
      F_cond <- first_stage_F(poi, pr$Xnames, pr$Znames, pr$dat)
      ci     <- ar_clr_ci_grid_ivmodel(pr$Y, pr$D, pr$Z, pr$X, alpha=alpha, halfwidth=halfwidth)
      tibble::tibble(price = poi, eq = eq_y, n = nrow(pr$dat),
                     rankX = qr(cbind(1, pr$X))$rank,
                     rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
                     k_excl = if(is.null(pr$Z)) 0L else ncol(pr$Z),
                     F_cond = F_cond,
                     AR_low = ci$ci_ar[1],  AR_high = ci$ci_ar[2],
                     CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2],
                     grid_center_2SLS = ci$center, grid_halfwidth = ci$halfwidth)
    })
  })
  
  # 4.3 — Resumo por preço (largura média das bandas)
  summary_by_price <- arclr |>
    group_by(price) |>
    summarise(
      F_cond_med = mean(F_cond, na.rm=TRUE),
      CLR_mean_width = mean(CLR_high - CLR_low, na.rm=TRUE),
      AR_mean_width  = mean(AR_high  - AR_low,  na.rm=TRUE),
      .groups="drop"
    )
  
  list(
    data_augmented = dfA,
    iv_names_aug   = ivA,
    first_stage    = sw_tbl,
    arclr_byprice_eq = arclr,
    summary_by_price = summary_by_price
  )
}

## =========== USO ===========
## Ajuste os nomes de shifters binários abaixo se quiser mais interações:
result_all <- run_iv_suite_all_prices(
  df            = df_iv_ok,
  priceNames    = priceNames,
  shareNames    = shareNames,
  omit_share    = best_fit$omit_share,
  drop_price    = best_fit$drop_price,
  base_iv_set   = best_iv_set,
  shifters_for_iv = c("SH_SH_area_1","SH_SH_capital_1"),  # edite/expanda aqui
  exogs         = c("z","z2"),
  alpha         = 0.05,
  halfwidth     = 50
)

cat("\n[ALL] 1º estágio condicional por preço (IVs aumentados):\n")
print(result_all$first_stage, digits=3)

cat("\n[ALL] AR/CLR por (preço, equação):\n")
print(head(result_all$arclr_byprice_eq, 15), digits=3)

cat("\n[ALL] Resumo por preço (largura média das bandas):\n")
print(result_all$summary_by_price, digits=3)

first_stage_conditional_table <- purrr::map_dfr(
  paste0("ln_", setdiff(priceNames, best_fit$drop_price)),
  function(poi){
    eqs <- setdiff(shareNames, best_fit$omit_share)
    vals <- purrr::map_dfr(eqs, function(eq_y){
      pr <- prep_XZ_cond(eq_y, poi, result_all$data_augmented,
                         result_all$iv_names_aug, best_fit$drop_price, priceNames)
      tibble::tibble(
        eq = eq_y,
        F_cond = first_stage_F(poi, pr$Xnames, pr$Znames, pr$dat),
        R2_partial = partial_R2(poi, pr$Xnames, pr$Znames, pr$dat)
      )
    })
    vals |>
      dplyr::summarise(
        F_cond_med = median(F_cond, na.rm=TRUE),
        F_cond_min = min(F_cond, na.rm=TRUE),
        R2p_med    = median(R2_partial, na.rm=TRUE)
      ) |>
      dplyr::mutate(price = poi, .before=1)
  }
)
print(first_stage_conditional_table, digits=3)

result_all_wide <- run_iv_suite_all_prices(
  df = df_iv_ok, priceNames = priceNames, shareNames = shareNames,
  omit_share = best_fit$omit_share, drop_price = best_fit$drop_price,
  base_iv_set = best_iv_set,
  shifters_for_iv = c("SH_SH_area_1","SH_SH_capital_1"),
  exogs = c("z","z2"),
  alpha = 0.05,
  halfwidth = 200   # <= aumente para 200–500
)
print(result_all_wide$summary_by_price, digits=3)

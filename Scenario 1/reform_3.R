## =========================================================
##  TOP-TIER SCRIPT: P1..P6 (Identificação → Robustez total)
##  Autor: você ; Assist: GPT
##  Data: Sys.Date()
## =========================================================
best_iv_set <- iv_set_full
## ---------- 0) Setup ----------
need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(sprintf("Pacotes ausentes: %s.\nInstale-os, ex.: install.packages(c(%s))",
                 paste(miss, collapse=", "),
                 paste(sprintf('"%s"', miss), collapse=", ")), call. = FALSE)
  }
}
need_pkg(c("AER","sandwich","clubSandwich","dplyr","purrr","tibble","MASS"))
## opcionais para IV fraco: ivmodel, ivpack
has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)
has_ivpack  <- requireNamespace("ivpack",  quietly = TRUE)

set.seed(12345)

## Sanidade do ambiente principal
stopifnot(exists("df_iv_ok"), is.data.frame(df_iv_ok))
stopifnot(exists("priceNames"), exists("shareNames"))
stopifnot(length(priceNames) >= 2, length(shareNames) >= 2)
stopifnot(exists("best_iv_set"), is.character(best_iv_set))
stopifnot(exists("best_fit"))

## Convenções
drop_price <- if (!is.null(best_fit$drop_price)) best_fit$drop_price else colnames(df_iv_ok)[1]
omit_share <- if (!is.null(best_fit$omit_share)) best_fit$omit_share else shareNames[1]

## Cluster (ajuste se tiver outro nome)
cluster_var <- c("cluster_id","id_domicilio","id_municipio","id_setor_censitario")
cluster_var <- cluster_var[cluster_var %in% names(df_iv_ok)]
if (!length(cluster_var)) {
  message("Aviso: nenhuma coluna de cluster encontrada; usando linhas como clusters (conservador).")
  df_iv_ok$.__rowid__ <- seq_len(nrow(df_iv_ok))
  cluster_var <- ".__rowid__"
}
cluster_var <- cluster_var[1]

## Exógenos do índice de Stone/QUAIDS na RHS
exogs <- intersect(c("z","z2"), names(df_iv_ok))

## Endógenos (preços log) exceto o dropado
endogs <- paste0("ln_", setdiff(priceNames, drop_price))
stopifnot(all(endogs %in% names(df_iv_ok)))

## ---------------------------------------------------------
##  Utils básicos (primeiro e segundo estágio, F SW condicional)
## ---------------------------------------------------------
## --------- versões à prova de X/Z vazios ---------
first_stage_F <- function(D, X, Z, data) {
  # sem IVs excluídos => F não definido
  if (!length(Z)) return(NA_real_)
  # fórmulas com fallback para intercepto (~ 1) quando X estiver vazio
  rhs0 <- if (length(X)) X else "1"
  rhs1 <- c(X, Z)
  if (!length(rhs1)) rhs1 <- "1"
  
  f0 <- lm(stats::reformulate(rhs0, response = D), data = data)
  f1 <- lm(stats::reformulate(rhs1, response = D), data = data)
  
  a <- anova(f0, f1)
  as.numeric(a$F[2])
}

partial_R2 <- function(D, X, Z, data) {
  if (!length(Z)) return(NA_real_)
  rhs0 <- if (length(X)) X else "1"
  rhs1 <- c(X, Z)
  if (!length(rhs1)) rhs1 <- "1"
  
  f0 <- lm(stats::reformulate(rhs0, response = D), data = data)
  f1 <- lm(stats::reformulate(rhs1, response = D), data = data)
  
  1 - sum(residuals(f1)^2) / sum(residuals(f0)^2)
}

## ---------------------------------------------------------
##  P1 — Reforço de identificação / Weak IV
##   - IVs aumentados (interações com dummies plausíveis)
##   - SW condicional e diagnósticos
##   - AR/CLR (ivmodel/ivpack) + mapeamento para elasticidades (conjuntos)
## ---------------------------------------------------------

## 1.1 Criar IVs aumentados (interações com shifters binários)
augment_ivs <- function(df, base_iv, with = c("SH_SH_area_1","SH_SH_capital_1")) {
  with <- intersect(with, names(df))
  out_names <- base_iv
  if (!length(with)) return(list(data = df, iv_names = unique(out_names)))
  for (b in base_iv) if (b %in% names(df)) {
    for (w in with) {
      nm <- paste0(b,"_x_", w)
      if (!nm %in% names(df)) df[[nm]] <- as.numeric(df[[b]]) * as.numeric(df[[w]])
      out_names <- c(out_names, nm)
    }
  }
  list(data = df, iv_names = unique(out_names))
}

aug <- augment_ivs(df_iv_ok, best_iv_set)
df_aug        <- aug$data
best_iv_aug   <- aug$iv_names

## 1.2 SW condicional e R2 parcial por preço endógeno
exogs  <- intersect(c("z","z2"), names(df_aug))   # pode ser character(0) sem problemas
endogs <- paste0("ln_", setdiff(priceNames, drop_price))

sw_tbl <- purrr::map_dfr(endogs, function(Dv){
  use <- unique(c(Dv, exogs, best_iv_aug))
  cc  <- stats::complete.cases(df_aug[, use, drop = FALSE])
  dat <- df_aug[cc, , drop = FALSE]
  tibble::tibble(
    var        = Dv,
    n_cc       = nrow(dat),
    F_SW       = first_stage_F(Dv, exogs, best_iv_aug, dat),
    R2_partial = partial_R2(Dv, exogs, best_iv_aug, dat),
    df1        = length(best_iv_aug),
    df2        = nrow(dat) - length(c(exogs, best_iv_aug)) - 1
  )
})
print(sw_tbl, digits = 3)

cat("\n[P1] SW condicional (IVs aumentados)\n"); print(sw_tbl, digits = 3)

## 1.3 AR/CLR (intervalos) — para cada equação de share vs. preço de interesse
##     Via 'ivmodel' (robusto) — se disponível — e fallback
ar_clr_ci_grid_ivmodel <- function(Y, D, Z, X = NULL, alpha = 0.05, grid = NULL) {
  if (!has_ivmodel) return(list(ci_ar = c(NA,NA), ci_clr = c(NA,NA)))
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  
  if (is.null(grid)) {
    b2 <- tryCatch({ unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]) }, error = function(e) NA_real_)
    if (!is.finite(b2)) b2 <- 0
    grid <- seq(b2 - 50, b2 + 50, length.out = 4001)
  }
  ivm <- ivmodel::ivmodel(Y = Y, D = D, Z = Z, X = X)
  p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0 = b0)$p.value, numeric(1))
  p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0 = b0)$p.value, numeric(1))
  keep_ar  <- which(p_ar  >= alpha)
  keep_clr <- which(p_clr >= alpha)
  ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
  ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
  list(ci_ar = ci_ar, ci_clr = ci_clr)
}

## monta X/Z condicionais corretos: D ~ exogs + (outros ln_p)  | exogs + IVs
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  stopifnot(poi %in% names(df), eq_y %in% names(df))
  Y <- as.numeric(df[[eq_y]]); D <- as.numeric(df[[poi]])
  ln_all <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  xnames <- c(intersect(c("z","z2"), names(df)), intersect(ln_others, names(df)))
  znames <- intersect(iv_names, names(df))
  use <- unique(c(eq_y, poi, xnames, znames))
  dat <- df[stats::complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  list(dat=dat, Y=as.numeric(dat[[eq_y]]), D=as.numeric(dat[[poi]]),
       X=if(length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL,
       Z=if(length(znames)) as.matrix(dat[, znames, drop = FALSE]) else NULL,
       Xnames=xnames, Znames=znames)
}

## roda AR/CLR para um preço de interesse (por eq de share)
price_of_interest <- "ln_preco_com_reforma13"  ## <-- ajuste se quiser outro
eqs <- setdiff(shareNames, omit_share)

p1_arclr <- purrr::map_dfr(eqs, function(eq_y){
  pr <- prep_XZ_cond(eq_y, price_of_interest, df_aug, best_iv_aug, drop_price, priceNames)
  F_cond <- first_stage_F(price_of_interest, pr$Xnames, pr$Znames, pr$dat)
  ci <- ar_clr_ci_grid_ivmodel(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05)
  tibble::tibble(eq = eq_y, n = nrow(pr$dat),
                 rankX = qr(cbind(1, pr$X))$rank,
                 rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
                 k_excl = if(is.null(pr$Z)) 0L else ncol(pr$Z),
                 F_cond = F_cond,
                 AR_low = ci$ci_ar[1], AR_high = ci$ci_ar[2],
                 CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2])
})
cat("\n[P1] AR/CLR (ivmodel) p/ ", price_of_interest, "\n", sep=""); print(p1_arclr, digits=3)

## 1.4 “Sets” de elasticidades via mapeamento do coeficiente (γ_ij) do preço de interesse
##      — substitui γ_{i, j*} por um valor dentro do CI e recomputa elasticidades
elas_with_gamma_override <- function(fit, j_name, gamma_vec_by_i,
                                     x_star, p_star,
                                     normalize_eval="softmax", normalize_deriv="softmax") {
  f2 <- fit
  G  <- f2$coef$gamma
  j  <- match(sub("^ln_", "", j_name), colnames(G))
  stopifnot(is.finite(j))
  for (i in seq_along(gamma_vec_by_i)) {
    ii <- match(names(gamma_vec_by_i)[i], rownames(G))
    if (is.finite(ii)) G[ii, j] <- gamma_vec_by_i[i]
  }
  f2$coef$gamma <- G
  elas_quaids_manual(f2, x=x_star, p=p_star,
                     normalize_eval=normalize_eval, normalize_deriv=normalize_deriv,
                     dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE)
}

## exemplo: usa extremos CLR por eq (quando existir) para produzir bandas diagonais de Marshall
px_emp <- list(
  x_star = median(df_iv_ok$gasto_total_atualhat, na.rm=TRUE),
  p_star = exp(colMeans(best_fit$data[paste0("ln_", priceNames)], na.rm=TRUE))
)
## monta vetor named para cada extremidade (por eq de share)
clr_low <- setNames(p1_arclr$CLR_low,  p1_arclr$eq)
clr_hi  <- setNames(p1_arclr$CLR_high, p1_arclr$eq)

## baixa
if (all(is.finite(clr_low))) {
  E_low <- elas_with_gamma_override(best_fit, price_of_interest, clr_low, px_emp$x_star, px_emp$p_star)
  cat("\n[P1] diag(M) com γ_{·,j} = CLR_low:\n"); print(round(diag(E_low$marshall), 3))
}
## alta
if (all(is.finite(clr_hi))) {
  E_hi <- elas_with_gamma_override(best_fit, price_of_interest, clr_hi, px_emp$x_star, px_emp$p_star)
  cat("\n[P1] diag(M) com γ_{·,j} = CLR_high:\n"); print(round(diag(E_hi$marshall), 3))
}

## ---------------------------------------------------------
##  P2 — Inferência sobre elasticidades: cluster bootstrap
## ---------------------------------------------------------
boot_elasticities <- function(fit, data, priceNames, B = 200, cluster = cluster_var,
                              x_star = NULL, p_star = NULL,
                              normalize_eval="softmax", normalize_deriv="softmax") {
  if (is.null(x_star)) x_star <- median(data$gasto_total_atualhat, na.rm=TRUE)
  if (is.null(p_star)) p_star <- exp(colMeans(fit$data[paste0("ln_", priceNames)], na.rm=TRUE))
  cl <- data[[cluster]]
  ids <- unique(cl); nC <- length(ids)
  diagM <- matrix(NA_real_, nrow=B, ncol=length(priceNames))
  eta   <- matrix(NA_real_, nrow=B, ncol=length(priceNames))
  
  for (b in seq_len(B)) {
    pick <- sample(ids, nC, replace = TRUE)
    idx  <- unlist(lapply(pick, function(k) which(cl == k)), use.names = FALSE)
    df_b <- data[idx, , drop = TRUE]
    ## re-estima por 3SLS (ou 2SLS eq-a-eq como fallback)
    fit_b <- try(
      fit_quaids_manual_km1(
        prices     = df_b[, priceNames, drop = FALSE],
        shares     = df_b[, fit$shareNames, drop = FALSE],
        x          = df_b[["gasto_total_atualhat"]],
        priceIndex = fit$priceIndex %||% "Ls",
        estMethod  = "3SLS",
        omit_share = match(fit$omit_share, fit$shareNames),
        drop_price = match(fit$drop_price, priceNames),
        instNames  = intersect(best_iv_aug, names(df_b)),
        instData   = df_b, use_z2 = "z2" %in% names(df_b),
        maxiter    = fit$maxiter %||% 500
      ), silent = TRUE)
    if (inherits(fit_b, "try-error")) next
    fit_b <- try(reparse_and_fix_quaids(fit_b), silent=TRUE)
    if (inherits(fit_b, "try-error")) next
    
    Eb <- try(elas_quaids_manual(fit_b, x=x_star, p=p_star,
                                 normalize_eval=normalize_eval, normalize_deriv=normalize_deriv,
                                 dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE), silent=TRUE)
    if (inherits(Eb, "try-error")) next
    diagM[b, ] <- suppressWarnings(diag(Eb$marshall))
    eta[b, ]   <- suppressWarnings(as.numeric(Eb$expenditure))
  }
  
  colnames(diagM) <- priceNames; colnames(eta) <- priceNames
  list(
    diagM_quant = apply(diagM, 2, quantile, na.rm=TRUE, probs=c(0.05,0.5,0.95)),
    eta_quant   = apply(eta,   2, quantile, na.rm=TRUE, probs=c(0.05,0.5,0.95))
  )
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

cat("\n[P2] Bootstrap (cluster) das elasticidades — rodando 200 reps...\n")
boot_out <- boot_elasticities(best_fit, df_aug, priceNames, B=200, cluster = cluster_var,
                              x_star = px_emp$x_star, p_star = px_emp$p_star)
print(boot_out$diagM_quant, digits=3); print(boot_out$eta_quant, digits=3)

## ---------------------------------------------------------
##  P3 — Sistema robusto: IV-GMM eq-a-eq + VCOV cluster
## ---------------------------------------------------------
## ---------------------------------------------------------
##  P3 — Sistema robusto: IV-GMM eq-a-eq + VCOV cluster
## ---------------------------------------------------------

## helpers
.f_rhs <- function(vars) if (length(vars)) vars else "1"

.fit_ivreg_one <- function(dat, y, endogs, exogs, iv_names) {
  rhs_y    <- c(endogs, exogs)                     # regressors na estrutural
  inst_rhs <- unique(c(exogs, iv_names))           # instrumentos = exógenos + IVs
  
  if (!length(inst_rhs))
    stop(sprintf("Eq %s: sem instrumentos disponíveis.", y))
  
  f_y <- stats::reformulate(.f_rhs(rhs_y), response = y)
  f_z <- stats::as.formula(paste("~", paste(.f_rhs(inst_rhs), collapse = " + ")))
  
  # checa IVs excluídos (instrumentos que não estão no RHS)
  excl <- setdiff(inst_rhs, rhs_y)
  
  fit <- AER::ivreg(formula = f_y, instruments = f_z, data = dat)
  list(fit = fit, k_excl = length(excl), excl = excl,
       f_y = f_y, f_z = f_z, rhs_y = rhs_y, inst_rhs = inst_rhs)
}

fit_ivreg_system <- function(df, shareNames, omit_share, endogs, exogs, iv_names, cluster) {
  stopifnot(cluster %in% names(df))
  eqs <- setdiff(shareNames, omit_share)
  
  out <- vector("list", length(eqs))
  names(out) <- eqs
  
  for (i in seq_along(eqs)) {
    y <- eqs[i]
    rhs_y    <- c(endogs, exogs)
    inst_rhs <- unique(c(exogs, iv_names))
    
    use <- unique(c(y, rhs_y, inst_rhs, cluster))
    cc  <- stats::complete.cases(df[, use, drop = FALSE])
    dat <- df[cc, , drop = FALSE]
    if (nrow(dat) == 0L) stop(sprintf("Eq %s: sem casos completos.", y))
    
    fit_info <- .fit_ivreg_one(dat, y, endogs, exogs, iv_names)
    
    ## VCOV cluster-robusto; se falhar, cai para HC1
    vc <- tryCatch(
      sandwich::vcovCL(fit_info$fit, cluster = dat[[cluster]]),
      error = function(e) sandwich::vcovHC(fit_info$fit, type = "HC1")
    )
    
    out[[i]] <- list(
      eq     = y,
      n      = nrow(dat),
      fit    = fit_info$fit,
      vcov   = vc,
      k_excl = fit_info$k_excl,
      excl   = fit_info$excl,
      f_y    = fit_info$f_y,
      f_z    = fit_info$f_z
    )
    
    if (fit_info$k_excl == 0)
      warning(sprintf("Eq %s: nenhum IV excluído (instrumentos = regressors).", y))
  }
  out
}

## Tidy: extrai tabela com erros-padrão cluster-robustos
ivsys_tidy <- function(ivsys_list) {
  do.call(
    rbind,
    lapply(ivsys_list, function(x) {
      ct <- lmtest::coeftest(x$fit, vcov. = x$vcov)
      data.frame(
        eq       = x$eq,
        term     = rownames(ct),
        estimate = ct[, 1],
        std.error= ct[, 2],
        z        = ct[, 3],
        p.value  = ct[, 4],
        n        = x$n,
        k_excl   = x$k_excl,
        row.names = NULL
      )
    })
  )
}

## ===== Uso =====
## endogs: e.g. paste0("ln_", setdiff(priceNames, drop_price))
## exogs : e.g. intersect(c("z","z2"), names(df_aug))  (pode ser character(0))
## ivs   : best_iv_aug  (precisa ter ao menos 1)
## cluster_var: nome da coluna de cluster em df_aug

ivsys <- fit_ivreg_system(
  df       = df_aug,
  shareNames = shareNames,
  omit_share = omit_share,
  endogs     = endogs,
  exogs      = exogs,
  iv_names   = best_iv_aug,
  cluster    = cluster_var
)

tab_ivsys <- ivsys_tidy(ivsys)
print(head(tab_ivsys), digits = 3)

## (opcional) checagens rápidas por equação
for (nm in names(ivsys)) {
  cat("\n--", nm, "--\n")
  m <- ivsys[[nm]]
  cat("n =", m$n, " | k_excl =", m$k_excl, " | excl =", paste(m$excl, collapse=", "), "\n")
  cat("formula: "); print(m$f_y)
  cat("instruments: "); print(m$f_z)
}


## ---------------------------------------------------------
##  P4 — Validade extra: leave-one-IV-out, split-sample, (opcional) plausible exog.
## ---------------------------------------------------------
loo_iv <- function(df, shareNames, omit_share, endogs, exogs, iv_names, cluster, target = price_of_interest) {
  base <- fit_ivreg_system(df, shareNames, omit_share, endogs, exogs, iv_names, cluster)
  b0 <- purrr::map_dbl(base, ~coef(.x$fit)[target])
  res <- list(b_base = b0)
  for (iv in iv_names) {
    iv_s <- setdiff(iv_names, iv)
    fit  <- fit_ivreg_system(df, shareNames, omit_share, endogs, exogs, iv_s, cluster)
    res[[paste0("drop_", iv)]] <- purrr::map_dbl(fit, ~coef(.x$fit)[target])
  }
  as_tibble <- function(x) tibble::tibble(eq = names(x), beta = as.numeric(x))
  out <- dplyr::bind_rows(lapply(names(res), function(k){
    a <- as_tibble(res[[k]]); a$spec <- k; a
  }))
  out
}
cat("\n[P4] Leave-one-IV-out no coeficiente de ", price_of_interest, ":\n", sep="")
print(loo_iv(df_aug, shareNames, omit_share, endogs, exogs, best_iv_aug, cluster_var) |>
        tidyr::pivot_wider(names_from = spec, values_from = beta), digits=3)

split_sample_compare <- function(df, shareNames, omit_share, endogs, exogs, iv_names, cluster) {
  set.seed(123)
  idx <- sample.int(nrow(df))
  S   <- rep(1:2, length.out = nrow(df)); S <- S[order(idx)]
  df$.__split__ <- S
  out <- vector("list", 2)
  for (s in 1:2) {
    fit <- fit_ivreg_system(df[df$.__split__==s, , drop=FALSE],
                            shareNames, omit_share, endogs, exogs, iv_names, cluster)
    out[[s]] <- purrr::map_dbl(fit, ~coef(.x$fit)[price_of_interest])
  }
  tibble::tibble(eq = names(out[[1]]),
                 beta_split1 = as.numeric(out[[1]]),
                 beta_split2 = as.numeric(out[[2]]),
                 diff = beta_split2 - beta_split1)
}
cat("\n[P4] Split-sample stability (metade/metade) — coef de ", price_of_interest, ":\n", sep="")
print(split_sample_compare(df_aug, shareNames, omit_share, endogs, exogs, best_iv_aug, cluster_var), digits=3)

## (Opcional) Plausible exogeneity (Conley) — se pacote disponível
if (requireNamespace("plausexog", quietly = TRUE)) {
  cat("\n[P4] Plausible exogeneity (Conley) — exemplo para uma eq:\n")
  y <- setdiff(shareNames, omit_share)[1]
  f_y <- stats::as.formula(paste(y, "~", paste(c(endogs, exogs), collapse = " + ")))
  f_z <- stats::as.formula(paste(y, "~", paste(c(exogs, best_iv_aug), collapse = " + ")))
  fit <- AER::ivreg(f_y | f_z, data = df_aug)
  pe  <- plausexog::plausexog(fit, param = price_of_interest, delta = c(-0.1, 0.1)) ## ajuste banda
  print(pe)
} else {
  message("[P4] 'plausexog' não disponível — pulando Conley bands (opcional).")
}

## ---------------------------------------------------------
##  P5 — Zeros e seleção: Shonkwiler–Yen (se função estiver disponível)
## ---------------------------------------------------------
sh_cols <- grep("^SH_", names(df_iv_ok), value = TRUE)
if (length(sh_cols) && exists("shonkwiler_yen_fit")) {
  cat("\n[P5] Shonkwiler–Yen — rodando (ponto empírico):\n")
  shy_fit <- try(shonkwiler_yen_fit(df_iv_ok, priceNames, shareNames, "gasto_total_atualhat",
                                    best_iv_aug, shifter_names = sh_cols), silent=TRUE)
  if (!inherits(shy_fit, "try-error")) {
    if (exists("elas_at")) {
      E_shy <- elas_at(shy_fit, x=px_emp$x_star, p=px_emp$p_star)
      cat("diag(M) [SHY]:\n"); print(round(diag(E_shy$marshall), 3))
      cat("eta [SHY]:\n"); print(round(E_shy$expenditure, 3))
    }
  } else {
    message("[P5] SHY falhou — verifique a função 'shonkwiler_yen_fit'.")
  }
} else {
  message("[P5] Sem shifters 'SH_*' ou 'shonkwiler_yen_fit' indisponível — pulando SHY.")
}

## ---------------------------------------------------------
##  P6 — Checagens adicionais (normalizações, subamostras, outliers, 1º estágio)
## ---------------------------------------------------------

## 6.1 Normalizações alternativas (muda omit_share/drop_price) no sistema IV eq-a-eq
alt_norm_results <- function(drop_price_alt = NULL, omit_share_alt = NULL) {
  dp <- drop_price_alt %||% drop_price
  os <- omit_share_alt %||% omit_share
  end <- paste0("ln_", setdiff(priceNames, dp))
  eqs <- setdiff(shareNames, os)
  fit <- fit_ivreg_system(df_aug, shareNames, os, end, exogs, best_iv_aug, cluster_var)
  tibble::tibble(
    eq = names(fit),
    coef_poi = purrr::map_dbl(fit, ~coef(.x$fit)[price_of_interest])
  )
}
cat("\n[P6] Normalização alternativa (troca drop/omit) — exemplo:\n")
dp_alt <- setdiff(priceNames, drop_price)[1]
os_alt <- setdiff(shareNames, omit_share)[1]
print(alt_norm_results(dp_alt, os_alt), digits=3)

## 6.2 Subamostras (ex.: capital vs interior, ou mediana de renda)
subsample_compare <- function(var, cutfun) {
  stopifnot(var %in% names(df_aug))
  g <- cutfun(df_aug[[var]])
  out <- list()
  for (lv in unique(g)) {
    dat <- df_aug[g == lv, , drop=FALSE]
    if (nrow(dat) < 50) next
    fit <- fit_ivreg_system(dat, shareNames, omit_share, endogs, exogs, best_iv_aug, cluster_var)
    out[[as.character(lv)]] <- purrr::map_dbl(fit, ~coef(.x$fit)[price_of_interest])
  }
  tibble::tibble(eq = names(out[[1]]),
                 across = names(out),
                 !!!lapply(out, as.numeric))
}
if ("SH_SH_area_1" %in% names(df_aug)) {
  cat("\n[P6] Subamostra por área (0=interior, 1=capital/urbano):\n")
  print(subsample_compare("SH_SH_area_1", function(x) x), digits=3)
}

## 6.3 Influência/outliers — DFBETAs do primeiro estágio (para o preço de interesse)
first_stage_influence <- function(D, X, Z, data) {
  f1 <- lm(stats::reformulate(c(X, Z), response = D), data = data)
  inf <- as.data.frame(dfbetas(f1))
  big <- which(apply(abs(inf), 1, max, na.rm=TRUE) > 0.2)  ## limiar exemplificativo
  list(n = nrow(data), n_big = length(big), idx = head(big, 10))
}
cat("\n[P6] Influência no 1º estágio de ", price_of_interest, ":\n", sep="")
print(first_stage_influence(price_of_interest, exogs, best_iv_aug, df_aug))

## 6.4 Tabela padrão de 1º estágio (R2 parcial, F condicional)
first_stage_table <- purrr::map_dfr(endogs, function(Dv){
  tibble::tibble(
    var = Dv,
    R2_partial = partial_R2(Dv, exogs, best_iv_aug, df_aug),
    F_SW       = first_stage_F(Dv, exogs, best_iv_aug, df_aug),
    n          = nrow(df_aug), k_IV = length(best_iv_aug)
  )
})
cat("\n[P6] Tabela de 1º estágio (padrão):\n"); print(first_stage_table, digits=3)

## 6.5 Reprodutibilidade
cat("\n[P6] Reprodutibilidade — sessionInfo():\n")
print(utils::sessionInfo())

cat("\n=== FIM: P1..P6 executados. Veja os blocos acima para diagnósticos e tabelas. ===\n")

## === KP rk F e Hansen–J (robusto) por equação no sistema P3 ===
kp_hj_table <- do.call(rbind, lapply(ivsys, function(m) {
  ss <- try(summary(m$fit, diagnostics = TRUE), silent = TRUE)
  kpF <- NA; kpP <- NA
  if (!inherits(ss, "try-error") && !is.null(ss$diagnostics)) {
    rn  <- rownames(ss$diagnostics)
    iKP <- grep("Kleibergen|Weak instruments", rn, ignore.case = TRUE)
    if (length(iKP)) {
      kpF <- unname(ss$diagnostics[iKP[1], "statistic"])
      kpP <- unname(ss$diagnostics[iKP[1], "p-value"])
    }
  }
  ## Hansen J robusto via momentos (g' W g) usando vcovCL residual-based:
  u   <- resid(m$fit)
  Zm  <- model.matrix(m$f_z, data = model.frame(m$fit))
  g   <- colMeans(Zm * u)
  ## HC sandwich para momentos (usar meat(B) = Z'u u'Z / n):
  meat <- crossprod(Zm * u) / nrow(Zm)
  Winv <- try(solve(meat), silent = TRUE)
  J    <- if (inherits(Winv,"try-error")) NA_real_ else nrow(Zm) * drop(t(g) %*% Winv %*% g)
  dfJ  <- max(0, ncol(Zm) - length(coef(m$fit)))  # (IVs - params)
  pJ   <- if (is.finite(J) && dfJ > 0) pchisq(J, df = dfJ, lower.tail = FALSE) else NA
  
  data.frame(eq = m$eq, KP_rkF = kpF, KP_p = kpP, HansenJ = J, dfJ = dfJ, HansenJ_p = pJ)
}))
print(kp_hj_table, digits = 3)

## === AR/CLR com grid adaptativo (mais informativo) ===
ar_clr_ci_adaptive <- function(Y,D,Z,X=NULL, alpha=0.05, start_w=5, max_w=200){
  if (!has_ivmodel) return(list(ci_ar=c(NA,NA), ci_clr=c(NA,NA)))
  b2 <- tryCatch(unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]), error=function(e) 0)
  ivm <- ivmodel::ivmodel(Y=Y, D=D, Z=Z, X=X)
  w <- start_w; ci_ar <- c(NA,NA); ci_clr <- c(NA,NA)
  while (w <= max_w) {
    grid <- seq(b2 - w, b2 + w, length.out = 4001)
    p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0=b0)$p.value,  numeric(1))
    p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0=b0)$p.value,  numeric(1))
    keep_ar  <- which(p_ar  >= alpha); keep_clr <- which(p_clr >= alpha)
    if (length(keep_ar))  ci_ar  <- c(min(grid[keep_ar]),  max(grid[keep_ar]))
    if (length(keep_clr)) ci_clr <- c(min(grid[keep_clr]), max(grid[keep_clr]))
    if (all(is.finite(ci_ar)) && all(is.finite(ci_clr))) break
    w <- w * 2
  }
  list(ci_ar=ci_ar, ci_clr=ci_clr, center=b2, width=w)
}

## ===== AR/CLR ADAPTATIVO (ivmodel) =====
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("AER","ivmodel","tibble","purrr","dplyr"))

# se df_aug/best_iv_aug não existirem, caia pro original
if(!exists("df_aug"))      df_aug    <- df_iv_ok
if(!exists("best_iv_aug")) best_iv_aug <- best_iv_set

# helper: X/Z "condicionais" corretos (D ~ exogs + outros ln_p | exogs + IVs)
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  stopifnot(poi %in% names(df), eq_y %in% names(df))
  Y <- as.numeric(df[[eq_y]]); D <- as.numeric(df[[poi]])
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  xnames <- c(intersect(c("z","z2"), names(df)), intersect(ln_others, names(df)))
  znames <- intersect(iv_names, names(df))
  use <- unique(c(eq_y, poi, xnames, znames))
  dat <- df[stats::complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  list(dat=dat, Y=as.numeric(dat[[eq_y]]), D=as.numeric(dat[[poi]]),
       X=if(length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL,
       Z=if(length(znames)) as.matrix(dat[, znames, drop = FALSE]) else NULL,
       Xnames=xnames, Znames=znames)
}

# AR/CLR com grade adaptativa (expande até fechar os intervalos ou atingir max)
ar_clr_ci_adaptive <- function(Y, D, Z, X=NULL, alpha=0.05,
                               start=50, step=50, max_width=2000) {
  if(!is.null(X)) X <- as.matrix(X)
  if(!is.null(Z)) Z <- as.matrix(Z)
  # centro da grade = 2SLS (fallback 0)
  b2 <- tryCatch({ unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]) }, error=function(e) NA_real_)
  if(!is.finite(b2)) b2 <- 0
  width <- start
  ivm <- ivmodel::ivmodel(Y=Y, D=D, Z=Z, X=X)
  
  one_pass <- function(width){
    grid <- seq(b2 - width, b2 + width, length.out = 4001)
    p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0=b0)$p.value,  numeric(1))
    p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0=b0)$p.value,  numeric(1))
    keep_ar  <- which(p_ar  >= alpha)
    keep_clr <- which(p_clr >= alpha)
    ci_ar  <- if(length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA,NA)
    ci_clr <- if(length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA,NA)
    list(ci_ar=ci_ar, ci_clr=ci_clr)
  }
  
  out <- one_pass(width)
  while( (any(!is.finite(out$ci_ar)) || any(!is.finite(out$ci_clr))) && width < max_width ){
    width <- width + step
    out <- one_pass(width)
  }
  c(out, list(center=b2, width=width))
}

# roda para TODAS as equações de share, para 1 preço de interesse
run_arclr_all <- function(price_of_interest, df, shareNames, omit_share,
                          drop_price, priceNames, iv_names){
  eqs <- setdiff(shareNames, omit_share)
  purrr::map_dfr(eqs, function(eq_y){
    pr <- prep_XZ_cond(eq_y, price_of_interest, df, iv_names, drop_price, priceNames)
    F_cond <- {
      f0 <- lm(pr$D ~ pr$X, data = pr$dat)
      f1 <- lm(pr$D ~ pr$X + pr$Z, data = pr$dat)
      as.numeric(anova(f0, f1)$F[2])
    }
    ci <- ar_clr_ci_adaptive(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05,
                             start = 50, step = 50, max_width = 2000)
    tibble::tibble(
      eq      = eq_y, n = nrow(pr$dat),
      rankX   = qr(cbind(1, pr$X))$rank,
      rankZ   = qr(cbind(1, pr$X, pr$Z))$rank,
      k_excl  = if(is.null(pr$Z)) 0L else ncol(pr$Z),
      F_cond  = F_cond,
      AR_low  = ci$ci_ar[1],  AR_high  = ci$ci_ar[2],
      CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2],
      grid_center_2SLS = ci$center, grid_halfwidth = ci$width
    )
  })
}

## === USO ===
price_of_interest <- "ln_preco_com_reforma13"  # ajuste se quiser
p1_arclr_adapt <- run_arclr_all(price_of_interest, df_aug, shareNames, best_fit$omit_share,
                                best_fit$drop_price, priceNames, best_iv_aug)
cat("\n[AR/CLR adaptativo] Coef. do preço de interesse:\n")
print(p1_arclr_adapt, digits = 3)

## ===== LIML / FULLER(1) vs 2SLS =====
need_pkg(c("AER","sandwich","tibble","purrr","dplyr"))

# kappa LIML estável + k-class interno (mesmo desenho que você já usou)
.kappa_liml_stable <- function(y, X, Z,
                               tau_grid = c(0,1e-12,1e-10,1e-8,1e-6)){
  n <- NROW(X); p <- ncol(X) + 1L
  Y <- cbind(y, X)
  ZZ <- crossprod(Z); ZZinv <- tryCatch(solve(ZZ), error=function(e) MASS::ginv(ZZ))
  PZy <- Z %*% (ZZinv %*% crossprod(Z, y))
  PZX <- Z %*% (ZZinv %*% crossprod(Z, X))
  MZy <- y - PZy; MZX <- X - PZX
  A0 <- crossprod(cbind(MZy, MZX))
  PX <- X %*% solve(crossprod(X), t(X))
  MXy <- y - PX %*% y; MXX <- X - PX %*% X
  B0 <- crossprod(cbind(MXy, MXX))
  A0 <- (A0 + t(A0))/2; B0 <- (B0 + t(B0))/2
  for(tau in tau_grid){
    A <- A0 + diag(tau, p)
    ch <- try(chol(A), silent=TRUE); if(inherits(ch,"try-error")) next
    RiT <- forwardsolve(t(ch), diag(p))
    S   <- RiT %*% B0 %*% t(RiT); S <- (S + t(S))/2
    ev  <- eigen(S, symmetric=TRUE, only.values=TRUE)$values
    k_raw <- min(Re(ev))
    return(min(max(k_raw, 1e-10), 1 - 1e-10))
  }
  stop("não estabilizou κ")
}
.kclass_coef <- function(y, X, Z, k){
  ZZ <- crossprod(Z); ZZinv <- tryCatch(solve(ZZ), error=function(e) MASS::ginv(ZZ))
  ZtX <- crossprod(Z, X); Zty <- crossprod(Z, y)
  XPZX <- t(ZtX) %*% ZZinv %*% ZtX
  XPZy <- t(ZtX) %*% ZZinv %*% Zty
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  Ak <- (1-k) * XtX + k * XPZX
  bk <- (1-k) * Xty + k * XPZy
  drop(solve(Ak, bk))
}

make_table_liml_fuller_2sls <- function(df, shareNames, omit_share,
                                        endogs, exogs, iv_names,
                                        price_of_interest, cluster){
  eqs <- setdiff(shareNames, omit_share)
  purrr::map_dfr(eqs, function(y){
    Xn <- c(endogs, exogs); Zn <- unique(c(exogs, iv_names))
    use <- unique(c(y, Xn, Zn, cluster))
    dat <- df[complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
    # 2SLS com SE cluster
    f_y <- stats::reformulate(Xn, response=y)
    f_z <- stats::as.formula(paste("~", paste(Zn, collapse=" + ")))
    fit2 <- AER::ivreg(formula=f_y, instruments=f_z, data=dat)
    vc2  <- tryCatch(sandwich::vcovCL(fit2, cluster=dat[[cluster]]),
                     error = function(e) sandwich::vcovHC(fit2, type="HC1"))
    b2   <- unname(coef(fit2)[price_of_interest])
    se2  <- unname(sqrt(diag(vc2))[price_of_interest])
    
    # LIML / Fuller(1) coef (sem SE)
    yv <- dat[[y]]
    Xm <- model.matrix(f_y, dat)
    Zm <- model.matrix(f_z, dat)
    k  <- .kappa_liml_stable(yv, Xm, Zm)
    bL <- .kclass_coef(yv, Xm, Zm, k = k)
    bF <- .kclass_coef(yv, Xm, Zm, k = 1 - 1/(nrow(Xm) - qr(Zm)$rank))  # Fuller(1)
    
    tibble::tibble(
      eq = y,
      b_2SLS = b2, se_2SLS = se2,
      b_LIML = unname(bL[colnames(Xm)==price_of_interest]),
      b_Fuller1 = unname(bF[colnames(Xm)==price_of_interest])
    )
  })
}

## === USO ===
endogs <- paste0("ln_", setdiff(priceNames, best_fit$drop_price))
exogs  <- intersect(c("z","z2"), names(df_aug))
cluster_var <- if(exists("cluster_var")) cluster_var else {
  df_aug$.__rowid__ <- seq_len(nrow(df_aug)); ".__rowid__"
}
tab_LF2 <- make_table_liml_fuller_2sls(df_aug, shareNames, best_fit$omit_share,
                                       endogs, exogs, best_iv_aug,
                                       price_of_interest, cluster_var)
cat("\n[2SLS vs LIML vs Fuller(1)] coef do preço de interesse:\n")
print(tab_LF2, digits = 3)

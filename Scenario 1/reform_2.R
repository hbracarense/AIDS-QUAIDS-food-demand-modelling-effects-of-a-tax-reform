## =========================================================
## Pacotes (apenas checagem; use os que já tem no projeto)
## =========================================================
need <- c("systemfit","sandwich","lmtest","AER")
opts <- sapply(need, requireNamespace, quietly = TRUE)
if (!all(opts)) stop("Faltam pacotes: ", paste(names(opts)[!opts], collapse=", "))

# Opcionais para testes weak-IV/AR/Fuller/Hansen GMM:
opt_pkgs <- c(ivmodel = "ivmodel", gmm = "gmm", ivreg = "ivreg", clubSandwich = "clubSandwich")
opt_has  <- sapply(opt_pkgs, requireNamespace, quietly = TRUE)
df_iv_ok <- df_iv
best_fit <- fit_q_3sls
best_iv_set <- iv_set_full
## =========================================================
## Utilitários básicos já compatíveis com seu ambiente
## (ajustam coeficientes e computam elasticidades)
## =========================================================

# 1) Reparse robusto dos coeficientes do systemfit (usa seus rótulos)
reparse_and_fix_quaids <- function(fit){
  stopifnot(!is.null(fit$fit))
  cf_all <- coef(fit$fit)
  nm_all <- names(cf_all)
  sn <- fit$shareNames; pn <- fit$priceNames; K <- length(pn)
  iomit <- match(fit$omit_share, sn); jdrop <- match(fit$drop_price, pn)
  if (is.na(iomit) || is.na(jdrop)) stop("omit_share/drop_price não encontrados.")
  keep_i <- setdiff(seq_len(K), iomit)
  ln_all <- paste0("ln_", pn)
  keep_ln <- setdiff(ln_all, paste0("ln_", fit$drop_price))
  alpha  <- setNames(rep(NA_real_, K), sn)
  beta   <- setNames(rep(NA_real_, K), sn)
  lambda <- setNames(rep(0, K), sn)
  gamma  <- matrix(NA_real_, K, K, dimnames = list(sn, pn))
  eq_labels <- gsub("[ _]", "", sn)
  
  for (ii in keep_i){
    eqlab <- eq_labels[ii]
    idx_eq <- grep(paste0("^", eqlab, "([_.:])"), nm_all)
    if (!length(idx_eq)) next
    nm_eq <- nm_all[idx_eq]; cf_eq <- cf_all[idx_eq]
    var <- sub(paste0("^", eqlab, "([_.:])"), "", nm_eq)
    a_pos <- which(var == "(Intercept)"); if (length(a_pos) == 1) alpha[ii]  <- unname(cf_eq[a_pos])
    b_pos <- which(var == "z");          if (length(b_pos) == 1) beta[ii]   <- unname(cf_eq[b_pos])
    l_pos <- which(var == "z2");         if (length(l_pos) == 1) lambda[ii] <- unname(cf_eq[l_pos])
    for (lnv in keep_ln){
      j <- match(sub("^ln_", "", lnv), pn)
      g_pos <- which(var == lnv); if (length(g_pos) == 1) gamma[ii, j] <- unname(cf_eq[g_pos])
    }
  }
  # adding-up (linha omitida) e homogeneidade (coluna dropada)
  alpha[iomit]  <- 1 - sum(alpha[keep_i],  na.rm = TRUE)
  beta[iomit]   <- -   sum(beta[ keep_i],  na.rm = TRUE)
  lambda[iomit] <- -   sum(lambda[keep_i], na.rm = TRUE)
  for (j in seq_len(K)) gamma[iomit, j] <- -sum(gamma[keep_i, j], na.rm = TRUE)
  keep_j <- setdiff(seq_len(K), jdrop)
  gamma[, jdrop] <- -rowSums(gamma[, keep_j, drop = FALSE], na.rm = TRUE)
  rs <- rowSums(gamma, na.rm = TRUE)
  if (any(abs(rs) > 1e-10)) gamma <- sweep(gamma, 1, rs / K, "-")
  
  # assegura sem NAs
  alpha[is.na(alpha)] <- 0; beta[is.na(beta)] <- 0; lambda[is.na(lambda)] <- 0
  gamma[is.na(gamma)] <- 0
  
  fit$coef$alpha  <- alpha
  fit$coef$beta   <- beta
  fit$coef$lambda <- lambda
  fit$coef$gamma  <- gamma
  fit
}

# 2) Wrapper de previsão de shares + elasticidades (usa sua função)
elas_at <- function(fit, x, p,
                    dlogp=1e-3,
                    normalize_eval="softmax", normalize_deriv="softmax",
                    enforce_hicks=TRUE, enforce_symmetry=TRUE){
  elas_quaids_manual(
    fit, x=x, p=p, enforce_hicks=enforce_hicks, enforce_symmetry=enforce_symmetry
  )
}

## =========================================================
## (a) Bootstrap por cluster para elasticidades (M, H, η)
## =========================================================

# Função para refit (tenta detectar assinatura c/ shifters)
.refit_quaids_any <- function(df, priceNames, shareNames, x_name, instNames,
                              priceIndex="Ls", estMethod="3SLS",
                              omit_share=1, drop_price=1, maxiter=400,
                              shifters=NULL){
  fmls <- names(formals(fit_quaids_manual_km1))
  args <- list(
    prices     = df[, priceNames,  drop=FALSE],
    shares     = df[, shareNames,  drop=FALSE],
    x          = df[[x_name]],
    priceIndex = priceIndex,
    estMethod  = estMethod,
    omit_share = omit_share,
    drop_price = drop_price,
    instNames  = instNames,
    instData   = df,
    use_z2     = TRUE,
    maxiter    = maxiter
  )
  if ("shifters" %in% fmls) args$shifters <- shifters
  fit <- do.call(fit_quaids_manual_km1, args)
  reparse_and_fix_quaids(fit)
}

# Bootstrap cluster (amostra clusters com reposição, refaz o fit e as elasticidades)
boot_elas_cluster <- function(B = 400, cluster_var,
                              df_all, priceNames, shareNames, x_name, instNames,
                              x_eval = c("median","mean")[1],
                              p_eval = c("sample-mean-log","median")[1],
                              seed = 123,
                              normalize_eval="softmax", normalize_deriv="softmax",
                              store_full = FALSE,
                              shifters = NULL){
  stopifnot(cluster_var %in% names(df_all))
  set.seed(seed)
  cl <- df_all[[cluster_var]]
  K  <- length(priceNames); S <- length(shareNames)
  if (K != S) stop("K!=S (número de preços e shares).")
  
  # ponto empírico
  x_star0 <- switch(x_eval,
                    "median" = median(df_all[[x_name]], na.rm=TRUE),
                    "mean"   = mean(df_all[[x_name]], na.rm=TRUE))
  if (p_eval == "sample-mean-log") {
    p_star0 <- exp(colMeans(df_all[paste0("ln_", priceNames)], na.rm=TRUE))
  } else {
    p_star0 <- apply(df_all[, priceNames, drop=FALSE], 2, median, na.rm=TRUE)
  }
  
  # objetos para acumular
  diagM <- matrix(NA_real_, nrow=B, ncol=K)
  colnames(diagM) <- shareNames
  eta   <- matrix(NA_real_, nrow=B, ncol=K)
  colnames(eta) <- shareNames
  lamaxH <- numeric(B)
  
  # loop bootstrap
  ucl <- unique(cl); ncl <- length(ucl)
  for(b in seq_len(B)){
    cb <- sample(ucl, size=ncl, replace=TRUE)
    idx <- which(cl %in% cb)
    dfb <- df_all[idx, , drop=FALSE]
    
    # refit
    fitb <- try(.refit_quaids_any(dfb, priceNames, shareNames, x_name, instNames,
                                  shifters = shifters), silent=TRUE)
    if (inherits(fitb, "try-error")) next
    
    # ponto de avaliação no bootstrap
    x_star <- if (x_eval == "median") median(dfb[[x_name]], na.rm=TRUE) else mean(dfb[[x_name]], na.rm=TRUE)
    p_star <- if (p_eval == "sample-mean-log") exp(colMeans(fitb$data[paste0("ln_", priceNames)], na.rm=TRUE))
    else apply(dfb[, priceNames, drop=FALSE], 2, median, na.rm=TRUE)
    
    Eb <- try(elas_at(fitb, x=x_star, p=p_star,
                      normalize_eval=normalize_eval, normalize_deriv=normalize_deriv), silent=TRUE)
    if (inherits(Eb, "try-error")) next
    
    diagM[b,] <- suppressWarnings(diag(Eb$marshall))
    eta[b, ]  <- Eb$expenditure
    Hsym <- (Eb$hicks + t(Eb$hicks))/2
    lamaxH[b] <- tryCatch(max(eigen(Hsym, symmetric=TRUE, only.values=TRUE)$values), error=function(e) NA_real_)
  }
  
  # sumários (EPs e ICs percentis)
  summ <- function(x){
    c(mean = mean(x, na.rm=TRUE),
      sd   = stats::sd(x, na.rm=TRUE),
      q025 = stats::quantile(x, 0.025, na.rm=TRUE),
      q975 = stats::quantile(x, 0.975, na.rm=TRUE))
  }
  diagM_tab <- t(apply(diagM, 2, summ))
  eta_tab   <- t(apply(eta,   2, summ))
  H_tab     <- summ(lamaxH)
  
  out <- list(
    boot = if (isTRUE(store_full)) list(diagM=diagM, eta=eta, lamaxH=lamaxH) else NULL,
    diagM = as.data.frame(diagM_tab),
    eta   = as.data.frame(eta_tab),
    hicks_lamax = H_tab
  )
  out
}

## =========================================================
## (b) Negatividade (Hicks) em grade + projeção NSD
## =========================================================

# avalia H em grade de (p, x); p é construído ao redor de p0 multiplicando fatores exp(d)
negativity_grid <- function(fit, x0, p0,
                            dlogp_grid = seq(-0.2, 0.2, length.out = 9),
                            x_mult     = c(0.6, 0.8, 1.0, 1.2, 1.5)){
  K <- length(p0)
  comb <- expand.grid(d = dlogp_grid, xmult = x_mult, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  lamax <- numeric(nrow(comb))
  for (i in seq_len(nrow(comb))){
    pj <- p0 * exp(rep(comb$d[i], K))
    xi <- x0 * comb$xmult[i]
    Ei <- elas_at(fit, x=xi, p=pj)
    Hs <- (Ei$hicks + t(Ei$hicks))/2
    lamax[i] <- max(eigen(Hs, symmetric=TRUE, only.values=TRUE)$values)
  }
  cbind(comb, lamax = lamax)
}

# projeção de uma matriz simétrica para semidefinida-negativa (Higham-like)
project_H_to_NSD <- function(H, eps = 0){
  ev <- eigen((H + t(H))/2, symmetric = TRUE)
  vals <- pmin(ev$values, -abs(eps))
  Hn <- ev$vectors %*% diag(vals, nrow=length(vals)) %*% t(ev$vectors)
  (Hn + t(Hn))/2
}

# reconstrói M via Slutsky após projeção
rebuild_M_from_H <- function(Hproj, w, eta){ Hproj - tcrossprod(eta, as.numeric(w)) }

## =========================================================
## (c) Weak-IV suite (SW-F, AR/CLR/LIML/Fuller quando possível)
## =========================================================
weak_iv_suite_v3 <- function(fit_sys, iv_set, df_all,
                             price_of_interest = NULL,
                             exogs_mode = c("stone","rhs"),
                             include_shifters_in_Z = FALSE) {
  stopifnot(requireNamespace("AER"),
            requireNamespace("sandwich"),
            requireNamespace("lmtest"))
  exogs_mode <- match.arg(exogs_mode)
  
  sn <- fit_sys$shareNames
  pn <- fit_sys$priceNames
  K  <- length(pn)
  
  # Endógenos = todos ln_preco_* exceto o dropado
  endogs_all <- setdiff(paste0("ln_", pn), paste0("ln_", fit_sys$drop_price))
  
  # Equações ativas (ordem por índice)
  i_omit <- match(fit_sys$omit_share, sn)
  keep_i <- setdiff(seq_len(K), i_omit)
  y_names <- sn[keep_i]
  
  # RHS por equação (fallback = só {z,z2})
  rhs_by_eq <- fit_sys$rhs_by_eq
  if (is.null(rhs_by_eq) || length(rhs_by_eq) != length(keep_i)) {
    rhs_by_eq <- replicate(length(keep_i),
                           intersect(c("z","z2"), names(fit_sys$data)),
                           simplify = FALSE)
  }
  
  # Exógenos de cada eq
  if (exogs_mode == "rhs") {
    exogs_by_eq <- lapply(rhs_by_eq, function(rhs) setdiff(rhs, endogs_all))
  } else {
    base_exogs <- intersect(c("z","z2"), names(fit_sys$data))
    if (include_shifters_in_Z) {
      sh_cols <- grep("^SH_", names(fit_sys$data), value = TRUE)
      base_exogs <- unique(c(base_exogs, sh_cols))
    }
    exogs_by_eq <- replicate(length(keep_i), base_exogs, simplify = FALSE)
  }
  
  # Pool de IVs
  iv_pool <- intersect(unique(iv_set), names(df_all))
  if (!length(iv_pool)) stop("Nenhum IV de iv_set encontrado em df_all.")
  
  # Dados combinados (sem duplicar colunas)
  dat0 <- cbind(
    fit_sys$data,
    df_all[, setdiff(iv_pool, names(fit_sys$data)), drop = FALSE]
  )
  
  # ---- helper: F-Wald robusto (HC1) para um bloco de coeficientes do lm grande ----
  F_block <- function(lm_big, vars_to_test) {
    b_all <- coef(lm_big)
    V_all <- sandwich::vcovHC(lm_big, type = "HC1")
    
    # --- ALINHAMENTO CRÍTICO: garante dimensões compatíveis ---
    common <- intersect(names(b_all), colnames(V_all))
    b <- b_all[common]
    V <- V_all[common, common, drop = FALSE]
    
    vars_to_test <- intersect(vars_to_test, common)
    if (!length(vars_to_test)) {
      return(c(F = NA_real_, df1 = NA_integer_, df2 = NA_integer_, p = NA_real_))
    }
    
    q <- length(vars_to_test)
    R <- matrix(0, nrow = q, ncol = length(b))
    colnames(R) <- names(b)
    for (j in seq_len(q)) {
      R[j, vars_to_test[j]] <- 1
    }
    
    Rb <- as.numeric(R %*% b)
    RV <- R %*% V %*% t(R)
    
    # proteções numéricas
    if (anyNA(Rb) || anyNA(RV) || qr(RV)$rank < nrow(RV)) {
      return(c(F = NA_real_, df1 = q, df2 = df.residual(lm_big), p = NA_real_))
    }
    
    stat <- drop(t(Rb) %*% solve(RV, Rb)) / q
    df1  <- q
    df2  <- df.residual(lm_big)
    pval <- 1 - pf(stat, df1, df2)
    c(F = stat, df1 = df1, df2 = df2, p = pval)
  }
  
  
  # ---------- 1) Por equação: Wu–Hausman, Sargan e min-F ----------
  by_eq <- lapply(seq_along(y_names), function(i) {
    y   <- y_names[i]
    exi <- intersect(exogs_by_eq[[i]], colnames(dat0))
    endo_eff <- intersect(endogs_all, colnames(dat0))
    iv_eff   <- setdiff(intersect(iv_pool, colnames(dat0)),
                        union(endo_eff, exi))
    
    # checagens / casos completos
    need <- unique(c(y, endo_eff, exi, iv_eff))
    cc   <- complete.cases(dat0[, need, drop = FALSE])
    dat  <- dat0[cc, , drop = FALSE]
    
    if (length(endo_eff) == 0L || length(iv_eff) == 0L) {
      return(data.frame(
        eq = y, n_cc = nrow(dat), k_endog_eff = length(endo_eff),
        k_excl_eff = length(iv_eff), F_weak_min = NA_real_,
        WuHausman_F = NA_real_, WuHausman_p = NA_real_,
        Sargan_stat = NA_real_, Sargan_p = NA_real_,
        check.names = FALSE
      ))
    }
    
    # Estrutural OLS
    f_ols <- as.formula(paste0(y, " ~ ", paste(c(endo_eff, exi), collapse=" + ")))
    ols   <- lm(f_ols, data = dat)
    
    # 2SLS para resíduos (para Wu-Hausman)
    f_iv  <- as.formula(paste0(
      y, " ~ ", paste(c(endo_eff, exi), collapse = " + "),
      " | ",     paste(c(exi, iv_eff), collapse = " + ")
    ))
    iv <- AER::ivreg(f_iv, data = dat)
    
    # Wu–Hausman residual-inclusion (robusto)
    Vres <- lapply(endo_eff, function(zj) {
      resid(lm(as.formula(paste0(zj, " ~ ", paste(c(exi, iv_eff), collapse = " + "))), data = dat))
    })
    names(Vres) <- endo_eff
    dat_HW <- cbind(dat[, c(y, endo_eff, exi), drop = FALSE],
                    as.data.frame(Vres, check.names = FALSE))
    f_hw   <- as.formula(paste0(
      y, " ~ ",
      paste(c(endo_eff, exi), collapse = " + "), " + ",
      paste(paste0("`", endo_eff, "`"), collapse = " + ")
    ))
    hw <- lm(f_hw, data = dat_HW)
    # Teste conjunto dos resíduos == 0
    res_names <- endo_eff
    wh <- F_block(hw, res_names)
    WuF <- unname(wh["F"]); Wup <- unname(wh["p"])
    
    # Sargan clássico (homoscedástico): J = n * R2(res2sls ~ Z)
    uhat <- resid(iv)
    Z    <- as.matrix(model.matrix(~ 0 + ., data = dat[, c(exi, iv_eff), drop = FALSE]))
    R2   <- summary(lm(uhat ~ Z))$r.squared
    J    <- length(uhat) * R2
    dfJ <- length(iv_eff) - length(endo_eff)
    Jp  <- if (dfJ > 0) 1 - pchisq(J, dfJ) else NA_real_

    
    # min-F fraco por 1º estágios: cada endógeno ~ exi + IVs, testando IVs=0 no modelo grande
    F_list <- sapply(endo_eff, function(zj) {
      # modelo grande com exi+IVs
      m1 <- lm(as.formula(paste0(zj, " ~ ", paste(c(exi, iv_eff), collapse = " + "))), data = dat)
      out <- F_block(m1, intersect(iv_eff, names(coef(m1))))
      unname(out["F"])
    })
    F_min <- suppressWarnings(min(as.numeric(F_list), na.rm = TRUE))
    if (!is.finite(F_min)) F_min <- NA_real_
    
    data.frame(
      eq = y,
      n_cc = nrow(dat),
      k_endog_eff = length(endo_eff),
      k_excl_eff  = length(iv_eff),
      F_weak_min  = F_min,
      WuHausman_F = if (is.finite(WuF))   WuF   else NA_real_,
      WuHausman_p = if (is.finite(Wup))   Wup   else NA_real_,
      Sargan_stat = if (is.finite(J))     J     else NA_real_,
      Sargan_p    = if (is.finite(Jp))    Jp    else NA_real_,
      check.names = FALSE
    )
  })
  tab_eq <- do.call(rbind, by_eq)
  
  # ---------- 2) SW condicional por endógeno ----------
  exi_sw <- intersect(c("z","z2"), names(dat0))
  if (include_shifters_in_Z) {
    exi_sw <- unique(c(exi_sw, grep("^SH_", names(dat0), value = TRUE)))
  }
  sw_rows <- lapply(endogs_all, function(y_end) {
    endo_eff <- intersect(endogs_all, colnames(dat0))
    iv_eff   <- setdiff(intersect(iv_pool, colnames(dat0)),
                        union(endo_eff, exi_sw))
    others   <- setdiff(endo_eff, y_end)
    need <- unique(c(y_end, others, exi_sw, iv_eff))
    cc   <- complete.cases(dat0[, need, drop = FALSE])
    dat  <- dat0[cc, , drop = FALSE]
    
    # modelo grande com (others + exi_sw + IVs); testa IVs == 0
    rhs_big <- unique(c(others, exi_sw, iv_eff))
    m1   <- lm(as.formula(paste0(y_end, " ~ ", paste(rhs_big, collapse = " + "))), data = dat)
    out  <- F_block(m1, intersect(iv_eff, names(coef(m1))))
    
    data.frame(var = y_end, n_cc = nobs(m1),
               F_SW = unname(out["F"]),
               df1  = unname(out["df1"]),
               df2  = unname(out["df2"]),
               p    = unname(out["p"]),
               check.names = FALSE)
  })
  tab_sw <- do.call(rbind, sw_rows)
  if (!is.null(price_of_interest)) {
    tab_sw <- subset(tab_sw, var == price_of_interest)
  }
  
  list(by_equation = tab_eq, sw_conditional = tab_sw)
}

## =========================================================
## (d) Hansen J (robusto) explícito por equação
##     (Na prática, summary(ivreg, vcovHC) já devolve "Sargan" mas é o J robusto.)
## =========================================================
hansenJ_by_equation_v4 <- function(fit, iv_names, data, exogs = c("z","z2")) {
  stopifnot(inherits(fit$fit, "systemfit"))
  
  # Equações estimadas e variáveis
  sn_all <- fit$shareNames
  sn     <- setdiff(sn_all, fit$omit_share)
  pn     <- fit$priceNames
  endogs <- paste0("ln_", setdiff(pn, fit$drop_price))
  
  # Exógenos e IVs efetivos que existem
  exi    <- intersect(exogs, names(fit$data))       # z, z2 vindos do fit$data
  iv_eff <- intersect(iv_names, names(data))        # IVs em df_iv_ok
  
  # Base consolidada: tudo o que pode entrar em Z e X
  base_dat <- cbind(fit$data, data[, iv_eff, drop = FALSE])
  
  out <- vector("list", length(sn))
  
  for (i in seq_along(sn)) {
    y <- sn[i]
    
    # X: regressors estruturais (intercepto + exógenos + endógenos)
    X_names <- c(exi, endogs)
    # Z: instrumentos (intercepto + exógenos + IVs excluídos)
    Z_names <- c(exi, iv_eff)
    
    use_vars <- unique(c(y, X_names, Z_names))
    cc <- stats::complete.cases(base_dat[, use_vars, drop = FALSE])
    dat <- base_dat[cc, , drop = FALSE]
    if (nrow(dat) == 0L) next
    
    # Fórmula IV para estimar resíduos u
    fml_str <- paste0(
      y, " ~ ", paste(X_names, collapse = " + "),
      " | ",  paste(Z_names, collapse = " + ")
    )
    fit_iv <- AER::ivreg(stats::as.formula(fml_str), data = dat)
    
    # Matrizes X e Z (com intercepto)
    X <- model.matrix(stats::reformulate(X_names, intercept = TRUE), data = dat)
    Z <- model.matrix(stats::reformulate(Z_names, intercept = TRUE), data = dat)
    
    # Post-estimation objects
    u  <- stats::resid(fit_iv)
    n  <- length(u)
    kx <- qr(X)$rank
    kz <- qr(Z)$rank
    
    # g = Z' u  (kz x 1)
    g <- crossprod(Z, u)
    
    # S = Z' diag(u^2) Z  (kz x kz) sem materializar diag:
    # t(Z) %*% (Z * u^2) == crossprod(Z, Z * u^2)
    Zu2 <- Z * as.numeric(u^2)
    S   <- crossprod(Z, Zu2)
    
    # Resolver S^{-1} g de forma estável
    S_inv_g <- try(solve(S, g), silent = TRUE)
    if (inherits(S_inv_g, "try-error")) {
      if (!requireNamespace("MASS", quietly = TRUE))
        stop("Instale 'MASS' para usar ginv() quando S é singular.")
      S_inv_g <- MASS::ginv(S) %*% g
    }
    
    # Estatística de J (escala não afeta p-valor)
    J   <- as.numeric(crossprod(g, S_inv_g))
    dfJ <- max(kz - kx, 0)
    if (dfJ <= 0) {
      warning(sprintf("Sem sobre-ID na eq %s (dfJ=%d).", y, dfJ))
      pJ <- NA_real_
    } else {
      pJ <- 1 - pchisq(J, dfJ)
    }
    
    out[[i]] <- data.frame(
      eq = y,
      n_cc = n,
      rankX = kx,
      rankZ = kz,
      dfJ = dfJ,
      HansenJ = J,
      HansenJ_p = pJ,
      row.names = NULL
    )
  }
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

## =========================================================
## (e) Sensibilidade: AIDS vs QUAIDS; Stone vs Ls; com/sem z2
## =========================================================

spec_sensitivity <- function(df, priceNames, shareNames, x_name, instNames,
                             base = list(priceIndex=c("Ls","S"),
                                         estMethod = "3SLS",
                                         z2 = c(TRUE, FALSE),
                                         quaids = c(TRUE, FALSE)),
                             eval_rule = list(x="median", p="sample-mean-log")){
  grid <- expand.grid(priceIndex = base$priceIndex,
                      z2         = base$z2,
                      quaids     = base$quaids,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  res <- list()
  for (i in seq_len(nrow(grid))){
    g <- grid[i,]
    fit <- fit_quaids_manual_km1(
      prices     = df[, priceNames, drop=FALSE],
      shares     = df[, shareNames, drop=FALSE],
      x          = df[[x_name]],
      priceIndex = g$priceIndex,
      estMethod  = "3SLS",
      omit_share = 1, drop_price = 1,
      instNames  = instNames,
      instData   = df,
      use_z2     = isTRUE(g$z2)
    )
    fit <- reparse_and_fix_quaids(fit)
    x_star <- if (eval_rule$x=="median") median(df[[x_name]], na.rm=TRUE) else mean(df[[x_name]], na.rm=TRUE)
    p_star <- if (eval_rule$p=="sample-mean-log") exp(colMeans(fit$data[paste0("ln_", priceNames)], na.rm=TRUE))
    else apply(df[, priceNames, drop=FALSE], 2, median, na.rm=TRUE)
    E <- elas_at(fit, x=x_star, p=p_star)
    res[[i]] <- list(spec=g, diagM = diag(E$marshall), eta = E$expenditure)
  }
  res
}

## =========================================================
## (f) Pesos de amostragem e clusters nos EPs
##     (Pesos entram no re-fit via bootstrap ponderado; EPs via cluster-bootstrap)
## =========================================================

# Se tiver coluna de peso 'peso', use sampling com prob ~ peso/sum(peso)
weighted_cluster_indices <- function(df, cluster_var, weight_var = NULL){
  cl <- df[[cluster_var]]
  ucl <- unique(cl)
  if (!is.null(weight_var) && weight_var %in% names(df)) {
    w <- tapply(df[[weight_var]], cl, mean, na.rm=TRUE)
    w <- w / sum(w, na.rm=TRUE)
  } else {
    w <- rep(1/length(ucl), length(ucl))
  }
  list(ucl = ucl, w = w)
}

## =========================================================
## (g) Zeros de consumo (Shonkwiler–Yen em modo robustez)
## =========================================================

# Para cada share w_i: probit de 1{w_i>0} ~ shifters, gera IMR_i; inclui IMR_i como shifter extra.
# Rodamos sistema com as IMRs adicionadas ao conjunto de shifters do RHS.
shonkwiler_yen_fit <- function(df, priceNames, shareNames, x_name, instNames,
                               shifter_names){
  # 1) Probits e IMRs
  IMR <- list()
  for (si in shareNames){
    yi <- as.numeric(df[[si]] > 0 & is.finite(df[[si]]))
    if (!any(yi==1) || !any(yi==0)) { IMR[[si]] <- rep(0, nrow(df)); next }
    Xp <- as.data.frame(df[, shifter_names, drop=FALSE])
    fml <- as.formula(paste0("yi ~ ", paste(colnames(Xp), collapse=" + ")))
    pr  <- glm(fml, family = binomial(link="probit"), data = data.frame(yi=yi, Xp))
    xb  <- drop(model.matrix(pr) %*% coef(pr))
    phi <- dnorm(xb); Phi <- pnorm(xb)
    imr <- phi / pmax(Phi, 1e-8)       # Mills inversa
    IMR[[si]] <- imr
  }
  IMR <- as.data.frame(IMR); colnames(IMR) <- paste0("IMR_", make.names(colnames(IMR)))
  df2 <- cbind(df, IMR)
  
  # 2) Refit com IMRs como shifters adicionais (se função aceitar 'shifters')
  fmls <- names(formals(fit_quaids_manual_km1))
  extra <- colnames(IMR)
  if ("shifters" %in% fmls) {
    fit <- fit_quaids_manual_km1(
      prices     = df2[, priceNames, drop=FALSE],
      shares     = df2[, shareNames, drop=FALSE],
      x          = df2[[x_name]],
      priceIndex = "Ls",
      estMethod  = "3SLS",
      omit_share = 1, drop_price = 1,
      instNames  = instNames,
      instData   = df2,
      shifters   = extra,
      use_z2     = TRUE
    )
  } else {
    # se sua versão não aceitar 'shifters', injete IMRs manualmente no df usado por systemfit (via rhs_by_eq)
    stop("Sua versão de fit_quaids_manual_km1 não aceita 'shifters'. Use a variante com 'shifters='.")
  }
  reparse_and_fix_quaids(fit)
}

## =========================================================
## (h) Regressores gerados & incerteza — já coberto pelo bootstrap
##     (refit completo a cada draw) — nada extra a fazer aqui.
## =========================================================


## ============================
## ====== EXEMPLO DE USO ======
## ============================

## Pressupõe que você já tem:
##   - df_iv_ok, priceNames, shareNames, best_iv_set
##   - best_fit (estimado na amostra cheia)
##   - funções elas_quaids_manual, make_iv_appendix
##   - variável de cluster (ajuste 'cluster_var' abaixo)

# 0) Ponto empírico (como você já fez)
px_emp <- list(
  x_star = median(df_iv_ok$gasto_total_atualhat, na.rm=TRUE),
  p_star = exp(colMeans(best_fit$data[paste0("ln_", priceNames)], na.rm=TRUE))
)

# (a) Bootstrap por cluster
# defina seu cluster (ex.: "muni_mes" ou "uf_ano" — troque pelo seu)
cluster_var <- "uf"   # <-- SUBSTITUA pelo nome certo na sua base
if (!cluster_var %in% names(df_iv_ok)) {
  message("⚠️ Defina 'cluster_var' com o nome correto do seu cluster em df_iv_ok.")
}

boot_out <- boot_elas_cluster(
  B = 400,
  cluster_var = cluster_var,
  df_all = df_iv_ok,
  priceNames = priceNames,
  shareNames = shareNames,
  x_name = "gasto_total_atualhat",
  instNames = best_iv_set,
  x_eval = "median",
  p_eval = "sample-mean-log",
  seed = 123,
  normalize_eval = "softmax",
  normalize_deriv = "softmax",
  store_full = FALSE
)

cat("\n=== Bootstrap elasticidades (diag Marshall) ===\n"); print(round(boot_out$diagM, 3))
cat("\n=== Bootstrap eta ===\n"); print(round(boot_out$eta, 3))
cat("\n=== Máx autovalor(H) na amostra bootstrap ===\n"); print(round(boot_out$hicks_lamax, 4))

# (b) Negatividade em grade + projeção
neg_scan <- negativity_grid(best_fit, x0 = px_emp$x_star, p0 = px_emp$p_star)
cat("\nMaior autovalor(H) na grade: ", round(max(neg_scan$lamax, na.rm=TRUE), 4), "\n")
# Se positivo: projeta H no ponto e reconstrói M (apenas para relatório de robustez)
E0 <- elas_at(best_fit, x=px_emp$x_star, p=px_emp$p_star)
if (max(neg_scan$lamax, na.rm=TRUE) > 0) {
  Hproj <- project_H_to_NSD(E0$hicks, eps = 0)
  Mproj <- rebuild_M_from_H(Hproj, w=E0$at$w, eta=E0$expenditure)
  cat("Após projeção: λ_max(Hproj) = ",
      round(max(eigen((Hproj+t(Hproj))/2, symmetric=TRUE, only.values=TRUE)$values), 6), "\n")
}

# (c) Weak-IV suite (foco em ln_preco_com_reforma13)
wiv <- weak_iv_suite_v3(
  best_fit, best_iv_set, df_iv_ok,
  price_of_interest   = "ln_preco_com_reforma13",
  exogs_mode          = "stone",   # só {z, z2}
  include_shifters_in_Z = FALSE    # sem SH_* em Z (como no appendix original)
)

cat("\n=== Appendix IV (SW-F, Sargan/Hansen-robusto) ===\n")
print(wiv$by_equation,    digits = 3)
cat("\n=== SW condicional por endógeno ===\n")
print(wiv$sw_conditional, digits = 3)



# (d) Hansen J explícito (robusto) por equação
HJ <- hansenJ_by_equation_v4(
  fit       = best_fit,
  iv_names  = best_iv_set,
  data      = df_iv_ok,
  exogs     = c("z","z2")
)

cat("\n=== Hansen J por equação (robusto/HC) ===\n"); print(HJ, digits=3)

# (e) Sensibilidade (AIDS vs QUAIDS; Stone vs Ls; com/sem z2)
sens <- spec_sensitivity(df_iv_ok, priceNames, shareNames, "gasto_total_atualhat", best_iv_set)
cat("\n=== Sensibilidade (diag(M), eta) por especificação ===\n")
for (i in seq_along(sens)){
  s <- sens[[i]]
  cat(sprintf("\nSpec #%d: index=%s, z2=%s, quaids=%s\n",
              i, s$spec$priceIndex, s$spec$z2, s$spec$quaids))
  print(round(s$diagM, 3))
  print(round(s$eta, 3))
}

# (f) Pesos: se tiver 'peso_amostral' na base, você pode direcionar o bootstrap
# (já fazemos cluster bootstrap; se quiser ponderar a seleção de clusters, use a função 'weighted_cluster_indices').

# (g) Zeros: Shonkwiler–Yen (robustez). Exemplo se seus shifters forem 'SH_*'.
sh_cols <- grep("^SH_", names(df_iv_ok), value = TRUE)
if (length(sh_cols)) {
  shy_fit <- try(shonkwiler_yen_fit(df_iv_ok, priceNames, shareNames, "gasto_total_atualhat",
                                    best_iv_set, shifter_names = sh_cols), silent=TRUE)
  if (!inherits(shy_fit, "try-error")) {
    E_shy <- elas_at(shy_fit, x=px_emp$x_star, p=px_emp$p_star)
    cat("\n=== Robustez Shonkwiler–Yen: diag(M) ===\n"); print(round(diag(E_shy$marshall), 3))
    cat("\n=== Robustez Shonkwiler–Yen: eta ===\n"); print(round(E_shy$expenditure, 3))
  } else {
    message("Shonkwiler–Yen: use a versão de fit_quaids_manual_km1 com argumento 'shifters'.")
  }
}

# (h) Regressor gerado: já coberto pelo (a) — refit em cada draw do bootstrap.

sw_conditional_all <- function(fit, iv_names, data, exogs = c("z","z2")) {
  pn  <- fit$priceNames
  endogs <- paste0("ln_", setdiff(pn, fit$drop_price))
  out <- lapply(endogs, function(v) {
    weak_iv_suite_v3(fit, iv_names, data,
                     price_of_interest = v, exogs_mode = "stone", include_shifters_in_Z = FALSE
    )$sw_conditional
  })
  do.call(rbind, out)
}
# uso:
sw_all <- sw_conditional_all(best_fit, best_iv_set, df_iv_ok)
print(sw_all, digits=3)

# ---------- κ LIML estável ----------
.kappa_liml_stable <- function(y, X, Z,
                               tau_grid = c(0, 1e-12, 1e-10, 1e-8, 1e-6),
                               verbose = TRUE) {
  n <- NROW(X); p <- ncol(X) + 1L
  Y <- cbind(y, X)
  
  # Projeções sem formar n x n
  ZZ <- crossprod(Z); ZZinv <- tryCatch(solve(ZZ), error = function(e) MASS::ginv(ZZ))
  PZy <- Z %*% (ZZinv %*% crossprod(Z, y))
  PZX <- Z %*% (ZZinv %*% crossprod(Z, X))
  
  # M_Z * y, M_Z * X sem formar M_Z:
  MZy <- y - PZy
  MZX <- X - PZX
  
  # A0 = Y' MZ Y ; B0 = Y' MX Y, com MX = I - PX
  A0 <- crossprod(cbind(MZy, MZX))
  PX  <- X %*% solve(crossprod(X), t(X))
  MXy <- y - PX %*% y
  MXX <- X - PX %*% X
  B0 <- crossprod(cbind(MXy, MXX))
  
  A0 <- (A0 + t(A0))/2; B0 <- (B0 + t(B0))/2
  
  # 1) Tenta geigen (sem inversão)
  if (requireNamespace("geigen", quietly = TRUE)) {
    ev <- try(geigen::geigen(B0, A0, symmetric = TRUE, only.values = TRUE)$values,
              silent = TRUE)
    if (!inherits(ev, "try-error")) {
      k_raw <- min(Re(ev))
      if (verbose) cat(sprintf("[κ via geigen] raw=%.6g\n", k_raw))
      return(min(max(k_raw, 1e-10), 1 - 1e-10))
    }
  }
  
  # 2) Fallback: chol(A0 + τI) e autovalor de A^{-1/2} B A^{-1/2}
  for (tau in tau_grid) {
    A <- A0 + diag(tau, p)
    ch <- try(chol(A), silent = TRUE)
    if (inherits(ch, "try-error")) next
    RiT <- forwardsolve(t(ch), diag(p))
    S   <- RiT %*% B0 %*% t(RiT)
    S   <- (S + t(S))/2
    ev  <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    k_raw <- min(Re(ev))
    if (verbose) cat(sprintf("[κ via chol, tau=%.0e] raw=%.6g\n", tau, k_raw))
    return(min(max(k_raw, 1e-10), 1 - 1e-10))
  }
  
  stop("Não consegui estabilizar κ.")
}

# ---------- k-class interno (sem ivreg::kclass) ----------
.kclass_coef <- function(y, X, Z, k, ridge = 0) {
  # W_k = (1-k)I + k PZ  => X'W_kX = (1-k)X'X + k X'PZ X
  ZZ <- crossprod(Z); ZZinv <- tryCatch(solve(ZZ), error = function(e) MASS::ginv(ZZ))
  ZtX <- crossprod(Z, X); Zty <- crossprod(Z, y)
  XPZX <- t(ZtX) %*% ZZinv %*% ZtX
  XPZy <- t(ZtX) %*% ZZinv %*% Zty
  
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  
  Ak <- (1 - k) * XtX + k * XPZX
  bk <- (1 - k) * Xty + k * XPZy
  
  if (ridge > 0) Ak <- Ak + diag(ridge, ncol(Ak))
  beta <- tryCatch(solve(Ak, bk), error = function(e) {
    # último recurso: pequeno ridge adaptativo
    solve(Ak + diag(1e-10, ncol(Ak)), bk)
  })
  drop(beta)
}

qr_prune <- function(M, keep_intercept = TRUE) {
  # remove colunas ~constantes
  sdok <- vapply(as.data.frame(M), function(v) sd(v) > 1e-12, TRUE)
  M <- M[, sdok, drop = FALSE]
  # QR com pivot
  q <- qr(M); r <- q$rank
  idx <- q$pivot[seq_len(r)]
  M2 <- M[, sort(idx), drop = FALSE]
  if (keep_intercept && any(colnames(M2) == "(Intercept)")) {
    # garante intercepto presente
    ii <- which(colnames(M2) == "(Intercept)")[1]
    M2 <- cbind(`(Intercept)` = M2[, ii], M2[, -ii, drop = FALSE])
  }
  M2
}


# ---------- LIML/Fuller por equação ----------
fit_liml_fuller <- function(fit, iv_names, data,
                            exogs = c("z","z2"),
                            fuller_alpha = NULL,
                            verbose = TRUE) {
  eqs    <- setdiff(fit$shareNames, fit$omit_share)
  endogs <- paste0("ln_", setdiff(fit$priceNames, fit$drop_price))
  exi    <- intersect(exogs, names(data))
  iv_eff <- intersect(iv_names, names(data))
  
  out <- vector("list", length(eqs))
  
  for (i in seq_along(eqs)) {
    yname   <- eqs[i]
    X_names <- c(exi, endogs)
    Z_names <- c(exi, iv_eff)
    
    excl <- setdiff(Z_names, X_names)
    if (!length(excl)) {
      warning(sprintf("Eq %s: não há IV excluído. Pulando.", yname))
      next
    }
    
    use <- unique(c(yname, X_names, Z_names))
    dat <- data[complete.cases(data[, use, drop = FALSE]), , drop = FALSE]
    if (nrow(dat) == 0L) next
    
    y <- dat[[yname]]
    X <- model.matrix(stats::reformulate(X_names, intercept = TRUE), dat)
    Z <- model.matrix(stats::reformulate(Z_names, intercept = TRUE), dat)
    X <- qr_prune(X, keep_intercept = TRUE)
    Z <- qr_prune(Z, keep_intercept = TRUE)
    kx <- qr(X)$rank
    kz <- qr(Z)$rank
    if (verbose) cat(sprintf("  (após prune) rank(X)=%d/%d, rank(Z)=%d/%d\n",
                             kx, ncol(X), kz, ncol(Z)))
    
    if (kz <= kx && verbose) cat("  (aviso) rank(Z) <= rank(X)\n")
    
    # κ LIML estável (retorna já clipeado p/ (0,1))
    k_liml <- .kappa_liml_stable(y, X, Z, verbose = verbose)
    if (verbose) cat(sprintf("  κ usado (clip): %.6g\n", k_liml))
    
    # Coeficientes LIML via k-class interno
    b_liml <- .kclass_coef(y, X, Z, k = k_liml)
    
    res <- list(
      eq     = yname,
      kappa  = k_liml,
      coef   = setNames(b_liml, colnames(X)),
      method = "LIML (k-class interno)"
    )
    
    # Fuller(α) opcional
    if (!is.null(fuller_alpha)) {
      kF <- 1 - fuller_alpha / (nrow(X) - kz)
      kF <- min(max(kF, 1e-10), 1 - 1e-10)
      if (verbose) cat(sprintf("  Fuller(%g): k = %.6g\n", fuller_alpha, kF))
      res$fuller <- list(
        alpha = fuller_alpha,
        kappa = kF,
        coef  = setNames(.kclass_coef(y, X, Z, k = kF), colnames(X)),
        method = sprintf("Fuller(%g) (k-class interno)", fuller_alpha)
      )
    }
    
    out[[i]] <- res
  }
  names(out) <- eqs
  out
}
fl  <- fit_liml_fuller(best_fit, best_iv_set, df_iv_ok, exogs = c("z","z2"))
# ou LIML + Fuller(1):
fl1 <- fit_liml_fuller(best_fit, best_iv_set, df_iv_ok, exogs = c("z","z2"),
                       fuller_alpha = 1)

# --- util: parser que não transforma ±Inf em NA ---
bounds2num_robust <- function(x){
  if (is.null(x)) return(c(NA_real_, NA_real_))
  to_val <- function(s){
    if (is.numeric(s)) return(as.numeric(s))
    s <- as.character(s)
    if (grepl("^-?Inf$", s)) return(suppressWarnings(as.numeric(s)))
    suppressWarnings(as.numeric(s))
  }
  v <- vapply(x, to_val, numeric(1))
  if (length(v) != 2) v <- c(NA_real_, NA_real_)
  v
}

# --- CI AR/CLR por grid (robusto/HC) com ivmodel ---
# CI AR/CLR por grid usando ivmodel (objeto), com fallback de centro no 2SLS
ar_clr_ci_grid_ivmodel <- function(Y, D, Z, X = NULL, alpha = 0.05, grid = NULL) {
  stopifnot(is.numeric(Y), is.numeric(D))
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  
  # Centro do grid no 2SLS (se falhar, usa 0)
  if (is.null(grid)) {
    b2 <- tryCatch({
      fit2sls <- AER::ivreg(Y ~ D + X | X + Z)
      unname(coef(fit2sls)["D"])
    }, error = function(e) NA_real_)
    if (!is.finite(b2)) b2 <- 0
    grid <- seq(b2 - 50, b2 + 50, length.out = 4001) 
  }
  
  # Cria o objeto ivmodel uma única vez
  ivm <- ivmodel::ivmodel(Y = Y, D = D, Z = Z, X = X)
  
  # p-values para cada beta0 da grade
  p_ar  <- sapply(grid, function(b0) ivmodel::AR.test(ivm,  beta0 = b0)$p.value)
  p_clr <- sapply(grid, function(b0) ivmodel::CLR(ivm, beta0 = b0)$p.value)
  
  # Conjuntos aceitos (p >= alpha)
  keep_ar  <- which(p_ar  >= alpha)
  keep_clr <- which(p_clr >= alpha)
  
  ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
  ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
  list(ci_ar = ci_ar, ci_clr = ci_clr, grid = grid, p_ar = p_ar, p_clr = p_clr)
}

# --- prepara X/Z "condicionais" corretamente ---
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  Y <- as.numeric(df[[eq_y]])
  stopifnot(poi %in% names(df))
  D <- as.numeric(df[[poi]])
  
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  xnames <- c(intersect(c("z","z2"), names(df)), intersect(ln_others, names(df)))
  znames <- intersect(iv_names, names(df))
  
  use <- unique(c(eq_y, poi, xnames, znames))
  dat <- df[complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  
  list(
    dat = dat,
    Y = as.numeric(dat[[eq_y]]),
    D = as.numeric(dat[[poi]]),
    X = if (length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL,
    Z = if (length(znames)) as.matrix(dat[, znames, drop = FALSE]) else NULL,
    Xnames = xnames,
    Znames = znames
  )
}

library(ivmodel)
library(AER)
library(dplyr)
library(purrr)

poi <- "ln_preco_com_reforma13"
eqs <- setdiff(shareNames, best_fit$omit_share)

res <- map(eqs, function(eq_y){
  pr <- prep_XZ_cond(eq_y, poi, df_iv_ok, best_iv_set,
                     drop_price = best_fit$drop_price,
                     priceNames = priceNames)
  
  # F condicional “correto”: exclusão conjunta de Z no 1º estágio (D ~ X + Z)
  f0 <- lm(pr$D ~ pr$X)           # sem Z
  f1 <- lm(pr$D ~ pr$X + pr$Z)    # com Z
  a  <- anova(f0, f1)
  F_cond <- as.numeric(a$F[2])
  
  # AR/CLR via grid (ivmodel)
  ci <- ar_clr_ci_grid_ivmodel(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05)
  
  # (Opcional) AR robusto via ivpack (hetero = TRUE). Requer fórmula.
  ar_ivpack <- tryCatch({
    # monta fórmulas: y ~ D + X | X + Z
    x_rhs <- if(length(pr$Xnames)) paste(pr$Xnames, collapse = " + ") else "1"
    z_rhs <- if(length(pr$Znames)) paste(pr$Znames, collapse = " + ") else "1"
    frm_y <- as.formula(paste(eq_y, "~", paste(c(poi, pr$Xnames), collapse = " + ")))
    frm_z <- as.formula(paste(eq_y, "~", paste(c(pr$Xnames, pr$Znames), collapse = " + ")))
    # estima 2SLS p/ centralizar beta
    fit2 <- AER::ivreg(frm_y | frm_z, data = pr$dat)
    b2   <- unname(coef(fit2)[poi]); if (!is.finite(b2)) b2 <- 0
    # teste AR hetero com grid curto (ajuste se quiser)
    ivpack::AndersonRubinTest(y = pr$dat[[eq_y]],
                              d = pr$dat[[poi]],
                              z = as.matrix(pr$dat[, pr$Znames, drop = FALSE]),
                              x = if(length(pr$Xnames)) as.matrix(pr$dat[, pr$Xnames, drop = FALSE]) else NULL,
                              beta = b2, hetero = TRUE)$CI
  }, error = function(e) c(NA_real_, NA_real_))
  
  tibble(
    eq = eq_y,
    n  = nrow(pr$dat),
    rankX = qr(cbind(1, pr$X))$rank,
    rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
    k_excl = if(is.null(pr$Z)) 0L else ncol(pr$Z),
    F_cond = F_cond,
    AR_low  = ci$ci_ar[1],  AR_high  = ci$ci_ar[2],
    CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2],
    AR_low_ivpack = ar_ivpack[1], AR_high_ivpack = ar_ivpack[2]
  )
}) %>% bind_rows()
options(tibble.width = Inf)
print(res, digits = 3)

print_ci <- function(lo, hi) {
  fmt <- function(x) ifelse(is.na(x), "NA",
                            ifelse(is.infinite(x), "∞", sprintf("%.3f", x)))
  paste0("[", fmt(lo), ", ", fmt(hi), "]")
}

# depois de calcular 'res' com grid ±50 ou ±100:
res$AR_CI  <- mapply(print_ci, res$AR_low,  res$AR_high)
res$CLR_CI <- mapply(print_ci, res$CLR_low, res$CLR_high)
res_out <- res[, c("eq","F_cond","AR_CI","CLR_CI")]
print(res_out, row.names = FALSE)

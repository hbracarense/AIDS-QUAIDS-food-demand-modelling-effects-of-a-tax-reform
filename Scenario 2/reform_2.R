# -----------------------------
# Pacotes
# -----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(haven)
  library(micEconAids)
  library(systemfit)
  library(AER)        # ivreg (diag 1ª etapa)
  library(zoo)
  library(splines)
  library(lmtest)
})
setDTthreads(0)

`%||%` <- function(x, y) if (!is.null(x)) x else y
.bt <- function(v) paste0("`", gsub("`", "\\`", v), "`")  # protege nomes

# -----------------------------
# Dados
# -----------------------------
path <- 'C:/Users/hbrac/Downloads/'
file <- 'banco_com correcao endogeneidade.dta'
df   <- read_dta(paste0(path, file))

# Conjuntos principais
prices <- as.data.frame(df[, paste0("preco_com_reforma2", 1:6), drop = FALSE])
shares <- as.data.frame(df[, paste0("w_despesahat", 1:6), drop = FALSE])
x      <- df[["gasto_total_atualhat"]]

K <- ncol(prices)
priceNames <- colnames(prices)
shareNames <- colnames(shares)

# Garantir numeric/positivo
for(nm in priceNames) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
for(nm in shareNames) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
df$gasto_total_atualhat <- suppressWarnings(as.numeric(df$gasto_total_atualhat))

ok <- rowSums(sapply(priceNames, function(v) is.finite(df[[v]]) & df[[v]]>0))==K &
  is.finite(df$gasto_total_atualhat) & df$gasto_total_atualhat>0
df_iv <- df[ ok, , drop = FALSE ]

# ln(preços) para RHS
lnp_cols <- paste0("ln_", priceNames)
for(j in seq_len(K)) df_iv[[lnp_cols[j]]] <- log(df_iv[[priceNames[j]]])

# ===============================================================
# IVs de renda: construção robusta
# ===============================================================
cand_income <- c("renda_total","renda_total_atualhat","renda","rendatotal",
                 "rendimento_total","renda_domiciliar")
inc_name <- cand_income[cand_income %in% names(df_iv)][1]
if (is.na(inc_name)) stop("Nao encontrei coluna de renda. Ajuste 'cand_income'.")

df_iv$ln_income <- log(suppressWarnings(as.numeric(df_iv[[inc_name]])))
df_iv <- df_iv[is.finite(df_iv$ln_income), , drop = FALSE]

ln_inc_c <- as.numeric(scale(df_iv$ln_income, TRUE, TRUE))
rng      <- range(df_iv$ln_income, na.rm = TRUE)
ln_inc_u <- (df_iv$ln_income - rng[1]) / (diff(rng)+1e-9)

build_income_IVs <- function(ln_inc_c, ln_inc_u,
                             deg_poly = 8, df_bs = 8, df_ns = 8,
                             fourier_K = 4,
                             hinge_probs = c(.2,.4,.6,.8),
                             add_interactions = TRUE) {
  OP <- poly(ln_inc_c, degree = deg_poly, raw = FALSE) %>% as.data.frame()
  colnames(OP) <- sprintf("iv_op%02d", seq_len(ncol(OP)))
  BS <- as.data.frame(bs(ln_inc_c, df = df_bs))
  NS <- as.data.frame(ns(ln_inc_c, df = df_ns))
  colnames(BS) <- sprintf("iv_bs%02d", seq_len(ncol(BS)))
  colnames(NS) <- sprintf("iv_ns%02d", seq_len(ncol(NS)))
  FT <- NULL
  if (fourier_K > 0) {
    FTlist <- lapply(seq_len(fourier_K), function(k){
      dfk <- data.frame(
        sin(2*pi*k*ln_inc_u),
        cos(2*pi*k*ln_inc_u),
        check.names = FALSE
      )
      names(dfk) <- c(paste0("iv_sin", k), paste0("iv_cos", k))
      dfk
    })
    FT <- do.call(cbind, FTlist)
  }
  qs <- quantile(ln_inc_c, probs = hinge_probs, na.rm = TRUE)
  H  <- sapply(qs, function(qk) pmax(0, ln_inc_c - qk))
  H  <- sweep(H, 2, colMeans(H, na.rm = TRUE), FUN = "-")
  H  <- as.data.frame(H)
  colnames(H) <- sprintf("iv_hg%02d", round(100*hinge_probs))
  INT <- NULL
  if (add_interactions) {
    op2 <- OP[, 1:min(2, ncol(OP)), drop = FALSE]
    h3  <- H[,  1:min(3, ncol(H)),  drop = FALSE]
    INTlist <- lapply(seq_len(ncol(op2)), function(i){
      out <- sweep(as.matrix(h3), 1, op2[[i]], `*`)
      out <- as.data.frame(out)
      colnames(out) <- paste0(colnames(op2)[i], "_x_", colnames(h3))
      out
    })
    INT <- do.call(cbind, INTlist)
  }
  IV <- cbind(OP, BS, NS, if (!is.null(FT)) FT, H, if (!is.null(INT)) INT)
  IV <- as.data.frame(lapply(IV, function(col){
    if (is.list(col)) unlist(col, use.names = FALSE) else as.numeric(col)
  }), check.names = FALSE)
  keep <- vapply(IV, function(v) {
    v <- as.numeric(v); is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12
  }, logical(1))
  IV[, keep, drop = FALSE]
}
IV_income <- build_income_IVs(ln_inc_c, ln_inc_u)

# ===============================================================
# PATCH (1): Tempo, fatores e dummies (robusto, com capital_interior)
# ===============================================================

pick_first <- function(cands) { cands[cands %in% names(df_iv)][1] %||% NA_character_ }
has_col    <- function(nm)    { !is.null(nm) && !is.na(nm) && nm %in% names(df_iv) }

nm_uf   <- pick_first(c("uf","UF","estado","sigla_uf","cod_uf"))
nm_reg  <- pick_first(c("regiao","região","region","macroregiao","macro_regiao"))
nm_area <- pick_first(c("area","area_metropolitana","zona","situacao","situacao_domicilio",
                        "urbano_rural","urbano_rural_flag","urbano_rural_cat"))
nm_cap  <- pick_first(c("capital","is_capital","capital_flag","capital_interior"))
nm_int  <- pick_first(c("interior","is_interior","interior_flag","capital_interior"))
nm_mes  <- pick_first(c("mes","month"))
nm_ano  <- pick_first(c("ano","year"))
nm_data <- pick_first(c("data","date"))

# Derivados de tempo
if (has_col(nm_data)) {
  if (!inherits(df_iv[[nm_data]], "Date")) suppressWarnings(df_iv[[nm_data]] <- as.Date(df_iv[[nm_data]]))
  df_iv$._ym  <- zoo::as.yearmon(df_iv[[nm_data]])
  df_iv$._mes <- as.integer(format(df_iv[[nm_data]], "%m"))
  df_iv$._ano <- as.integer(format(df_iv[[nm_data]], "%Y"))
} else {
  if (has_col(nm_mes)) df_iv$._mes <- suppressWarnings(as.integer(df_iv[[nm_mes]]))
  if (has_col(nm_ano)) df_iv$._ano <- suppressWarnings(as.integer(df_iv[[nm_ano]]))
  if (all(c("._mes","._ano") %in% names(df_iv))) {
    df_iv$._ym <- zoo::as.yearmon(paste(df_iv$._ano, df_iv$._mes), "%Y %m")
  }
}
if ("._mes" %in% names(df_iv)) df_iv$._tri <- ((pmax(1, pmin(12, df_iv$._mes))-1L) %/% 3L) + 1L

# Fatores
if (has_col(nm_uf))   df_iv$f_uf   <- factor(df_iv[[nm_uf]])
if (has_col(nm_reg))  df_iv$f_reg  <- factor(df_iv[[nm_reg]])
if (has_col(nm_area)) df_iv$f_area <- factor(df_iv[[nm_area]])
if (has_col(nm_cap))  df_iv$f_cap  <- factor(df_iv[[nm_cap]])
if (has_col(nm_int))  df_iv$f_int  <- factor(df_iv[[nm_int]])

if ("._mes" %in% names(df_iv)) df_iv$f_mes <- factor(df_iv$._mes)
if ("._ano" %in% names(df_iv)) df_iv$f_ano <- factor(df_iv$._ano)
if ("._ym"  %in% names(df_iv)) df_iv$f_ym  <- factor(df_iv$._ym)
if ("._tri" %in% names(df_iv)) df_iv$f_tri <- factor(df_iv$._tri)

# Interação UF×TRI (se ambos existem)
if (all(c("f_uf","f_tri") %in% names(df_iv))) {
  df_iv$f_uf_tri <- interaction(df_iv$f_uf, df_iv$f_tri, drop = TRUE, lex.order = TRUE)
}

# Dummies helper
make_dummies <- function(x, prefix){
  f <- factor(x, exclude = NULL)
  X <- model.matrix(~ 0 + f)
  colnames(X) <- paste0(prefix, "_", make.names(levels(f), allow_ = FALSE))
  as.data.frame(X, check.names = FALSE)
}

# Dummies (cria apenas se fator existir)
D_uf  <- if ("f_uf"  %in% names(df_iv)) make_dummies(df_iv$f_uf,  "IV_d_uf")  else NULL
D_reg <- if ("f_reg" %in% names(df_iv)) make_dummies(df_iv$f_reg, "IV_d_reg") else NULL
D_mes <- if ("f_mes" %in% names(df_iv)) make_dummies(df_iv$f_mes, "IV_d_mes") else NULL

# capital_interior: gera IV_d_cap_*; evita duplicar "int" se for igual à cap
D_cap <- NULL; D_int <- NULL
if ("f_cap" %in% names(df_iv)) D_cap <- make_dummies(df_iv$f_cap, "IV_d_cap")
if ("f_int" %in% names(df_iv)) {
  same_source <- ("f_cap" %in% names(df_iv)) &&
    identical(as.character(df_iv$f_int), as.character(df_iv$f_cap))
  if (!same_source) D_int <- make_dummies(df_iv$f_int, "IV_d_int")
}

# UF×MÊS (se existirem e com bound de níveis)
D_ufmes <- NULL
if (!is.null(D_uf) && !is.null(D_mes)) {
  f_uf  <- factor(df_iv$f_uf,  exclude = NULL)
  f_mes <- factor(df_iv$f_mes, exclude = NULL)
  ufxm  <- interaction(f_uf, f_mes, drop = TRUE, lex.order = TRUE)
  if (nlevels(ufxm) <= 200) D_ufmes <- make_dummies(ufxm, "IV_d_ufmes")
}

# Interações renda×região/capital (moderadas)
OP_tmp <- as.data.frame(poly(ln_inc_c, degree = 2, raw = FALSE))
colnames(OP_tmp) <- c("IVop1","IVop2")
INT <- NULL
if (!is.null(D_reg)) {
  dR  <- D_reg[, 1:min(3, ncol(D_reg)), drop = FALSE]
  INT <- do.call(cbind, lapply(seq_len(ncol(OP_tmp)), function(i){
    out <- sweep(as.matrix(dR), 1, OP_tmp[[i]], `*`)
    colnames(out) <- paste0(colnames(OP_tmp)[i], "_x_", colnames(dR))
    out
  })) %>% as.data.frame(check.names = FALSE)
}
if (!is.null(D_cap)) {
  INT2 <- sweep(as.matrix(D_cap), 1, OP_tmp[[1]], `*`)
  colnames(INT2) <- paste0(colnames(OP_tmp)[1], "_x_", colnames(D_cap))
  INT <- cbind(INT, as.data.frame(INT2, check.names = FALSE))
}

# Contínuas exógenas adicionais (se existirem)
dem_reg_cont <- names(df_iv)[grepl("^(dem_|reg_)", names(df_iv))]
dem_reg_cont <- dem_reg_cont[!dem_reg_cont %in% c(shareNames, priceNames, "gasto_total_atualhat")]

# Monta matriz de IVs candidatos
IV_cand_df <- do.call(
  cbind,
  Filter(Negate(is.null), list(
    IV_income,
    D_uf, D_mes, D_reg, D_cap, D_int, D_ufmes,
    INT,
    if (length(dem_reg_cont)) df_iv[dem_reg_cont] else NULL
  ))
)
IV_cand_df <- as.data.frame(lapply(IV_cand_df, function(v) as.numeric(v)), check.names = FALSE)
colnames(IV_cand_df) <- make.names(colnames(IV_cand_df), unique = TRUE)
keep_var <- vapply(IV_cand_df, function(v){ vv <- var(v, na.rm = TRUE); is.finite(vv) && vv > 1e-12 }, logical(1))
IV_cand_df <- IV_cand_df[, keep_var, drop = FALSE]

# Cola e atualiza conjunto de candidatos
df_iv <- cbind(df_iv, IV_cand_df)
iv_all_candidates <- colnames(IV_cand_df)

# Posto de Z amplo
Z0  <- model.matrix(stats::reformulate(iv_all_candidates, intercept = FALSE), data = df_iv)
qr0 <- qr(Z0)
iv_base <- colnames(Z0)[qr0$pivot[seq_len(qr0$rank)]]
req_endog <- length(priceNames) - 1L
cat("Endogenas por equacao (K-1):", req_endog, "\n")
cat("Rank(Z) com IVs amplos:", length(iv_base), "\n")

# ===============================================================
# PATCH (2): conjuntos de IVs (core / mid / full) + saneamento
# ===============================================================
head_grep <- function(pat, n = Inf, x = names(df_iv)) head(grep(pat, x, value = TRUE), n)
N_OP <- 3; N_UF <- 5; N_MES <- 6; N_CAP <- 2; N_INT <- 2; N_UFMES <- 20

iv_set_core <- unique(c(
  head_grep("^iv_op",       N_OP),
  head_grep("^IV_d_uf_",    N_UF),
  head_grep("^IV_d_mes_",   N_MES),
  head_grep("^IV_d_cap_",   N_CAP),
  head_grep("^IV_d_int_",   N_INT),
  head_grep("^IV_d_ufmes_", N_UFMES)
))

# Interações geradas (após make.names) que combinam IVop* com IV_d_reg_*
INT_in_df <- grep("^IVop[0-9]+_x_IV_d_reg_", names(df_iv), value = TRUE)

iv_set_mid <- unique(c(
  iv_set_core,
  grep("^IV_d_reg_", names(df_iv), value = TRUE),
  INT_in_df
))

iv_set_full <- unique(intersect(iv_all_candidates %||% character(0), names(df_iv)))

# Diagnóstico opcional:
iv_rank <- function(vars, data = df_iv) {
  if (length(vars) == 0) return(0L)
  out <- try({
    Z <- model.matrix(stats::reformulate(vars, intercept = FALSE), data = data)
    keep <- apply(Z, 2, function(v) is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12)
    Z <- Z[, keep, drop = FALSE]
    qr(Z)$rank
  }, silent = TRUE)
  if (inherits(out, "try-error")) NA_integer_ else out
}
cat(sprintf("IV counts — core: %d | mid: %d | full: %d\n",
            length(iv_set_core), length(iv_set_mid), length(iv_set_full)))
cat(sprintf("IV rank  — core: %s | mid: %s | full: %s\n",
            iv_rank(iv_set_core), iv_rank(iv_set_mid), iv_rank(iv_set_full)))

# ===============================================================
# Benchmarks SUR (opcional)
# ===============================================================

# SL (Stone defasado)
aids_sl <- aidsEst(priceNames = colnames(prices),
                   shareNames = colnames(shares),
                   data = data.frame(prices, shares, x = x),
                   totExpName = "x", method = "LA",
                   priceIndex = "SL", estMethod = "SUR", maxiter = 200)
cat("\n==== Resultado SL-SUR ====\n"); print(summary(aids_sl)); print(elas(aids_sl))

# Ls (Laspeyres simplificado)
aids_ls <- aidsEst(priceNames = colnames(prices),
                   shareNames = colnames(shares),
                   data = data.frame(prices, shares, x = x),
                   totExpName = "x", method = "LA",
                   priceIndex = "Ls", estMethod = "SUR", maxiter = 200)
cat("\n==== Resultado Ls-SUR ====\n"); print(summary(aids_ls)); print(elas(aids_ls))

# ===============================================================
# QUAIDS manual (funções) — com 3SLS tolerante a IVs ausentes
# ===============================================================

fit_quaids_manual <- function(prices, shares, x,
                              priceIndex = c("Ls","S"),
                              estMethod  = c("SUR","3SLS"),
                              baseShares = NULL,
                              instNames  = NULL,
                              maxiter    = 200,
                              base_drop  = 1){
  stopifnot(is.data.frame(prices), is.data.frame(shares), length(x) == nrow(prices))
  stopifnot(nrow(prices) == nrow(shares))
  priceIndex <- match.arg(priceIndex)
  estMethod  <- match.arg(estMethod)
  
  priceNames <- colnames(prices)
  shareNames <- colnames(shares)
  K <- ncol(prices)
  if (!(base_drop %in% seq_len(K))) stop("base_drop deve estar entre 1 e K.")
  base_name <- priceNames[base_drop]
  
  as_num <- function(v) suppressWarnings(as.numeric(v))
  prices[] <- lapply(prices, as_num)
  shares[] <- lapply(shares, as_num)
  x <- as_num(x)
  
  ok <- rowSums(sapply(prices, function(v) is.finite(v) & v > 0)) == K & is.finite(x) & x > 0
  P  <- prices[ok, , drop = FALSE]
  W  <- shares[ok, , drop = FALSE]
  X  <- x[ok]
  
  lnP <- switch(priceIndex,
                "Ls" = {
                  wbar <- if (is.null(baseShares)) colMeans(W, na.rm = TRUE) else as.numeric(baseShares)
                  if (length(wbar) != K) stop("baseShares deve ter K elementos.")
                  wbar <- wbar / sum(wbar)
                  as.numeric(as.matrix(log(P)) %*% wbar)
                },
                "S"  = {
                  w_i <- W / rowSums(W)
                  rowSums(w_i * log(P))
                })
  z  <- log(X) - lnP
  z2 <- z^2
  
  lnPcols <- paste0("ln_", priceNames)
  df <- data.frame(W, setNames(as.data.frame(log(P)), lnPcols), z = z, z2 = z2, check.names = FALSE)
  
  ln_kept_global <- lnPcols[-base_drop]
  eqLabels <- gsub("[ _]", "", shareNames)
  
  build_rhs_eq <- function(lhs) {
    rhs_all <- c(ln_kept_global, "z", "z2")
    fml_all <- as.formula(paste0(lhs, " ~ ", paste(c("1", rhs_all), collapse = " + ")))
    Xall <- model.matrix(fml_all, data = df)
    keep <- colnames(Xall) %in% c("(Intercept)", "z")
    if ("z2" %in% colnames(Xall)) {
      X_try <- Xall[, keep | (colnames(Xall) == "z2"), drop = FALSE]
      r1 <- qr(X_try)$rank
      if (r1 > qr(Xall[, keep, drop = FALSE])$rank) keep[colnames(Xall) == "z2"] <- TRUE
    }
    for (lnv in ln_kept_global) {
      X_try <- Xall[, keep | (colnames(Xall) == lnv), drop = FALSE]
      if (qr(X_try)$rank > qr(Xall[, keep, drop = FALSE])$rank) keep[colnames(Xall) == lnv] <- TRUE
    }
    rhs_vars <- setdiff(colnames(Xall)[keep], "(Intercept)")
    list(rhs = paste(c("1", rhs_vars), collapse = " + "),
         kept = rhs_vars,
         dropped = setdiff(rhs_all, rhs_vars))
  }
  
  rhs_list   <- vector("list", K)
  eqs        <- vector("list", K); names(eqs) <- eqLabels
  kept_by_eq <- vector("list", K)
  drop_by_eq <- vector("list", K)
  
  for (i in seq_len(K)) {
    sel <- build_rhs_eq(shareNames[i])
    rhs_list[[i]] <- sel$rhs
    kept_by_eq[[i]] <- sel$kept
    drop_by_eq[[i]] <- sel$dropped
    eqs[[i]] <- as.formula(paste0(shareNames[i], " ~ ", sel$rhs))
  }
  
  if (estMethod == "SUR") {
    fit <- systemfit::systemfit(eqs, data = df, method = "SUR", maxit = maxiter)
  } else {
    if (is.null(instNames) || length(instNames) == 0) stop("Sem instrumentos (instNames) fornecidos.")
    # Intersecta com df e expande
    inst_in_df <- intersect(instNames, names(df))
    if (!length(inst_in_df)) stop("Sem instrumentos após interseção com dados.")
    Z <- model.matrix(stats::reformulate(inst_in_df, intercept = FALSE), data = df)
    keepZ <- apply(Z, 2, function(v) is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12)
    Z <- Z[, keepZ, drop = FALSE]
    qrZ <- qr(Z)
    Z_ok_names <- colnames(Z)[qrZ$pivot[seq_len(qrZ$rank)]]
    inst_terms <- paste(.bt(Z_ok_names), collapse = " + ")
    inst_fml <- as.formula(paste("~ 0 +", inst_terms))
    fit <- systemfit::systemfit(eqs, data = df, method = "3SLS",
                                inst = inst_fml, maxit = maxiter)
  }
  
  cf <- coef(fit); se <- sqrt(diag(vcov(fit)))
  pick_eq <- function(eqLabel) {
    nm <- names(cf); idx <- grep(paste0("^", eqLabel, "[\\.:]"), nm)
    list(cf = cf[idx], se = se[idx], nm = nm[idx])
  }
  
  alpha <- setNames(numeric(K), shareNames)
  beta  <- setNames(numeric(K), shareNames)
  lambda<- setNames(numeric(K), shareNames)
  gamma <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  
  ln_to_j <- function(lnv) match(sub("^ln_", "", lnv), priceNames)
  for (i in seq_len(K)) {
    pi <- pick_eq(eqLabels[i])
    aidx <- grep("\\(Intercept\\)", pi$nm); if (length(aidx)==1) alpha[i] <- unname(pi$cf[aidx])
    bidx <- grep("([\\.:]|^)z$",  pi$nm);   if (length(bidx)==1) beta[i]   <- unname(pi$cf[bidx])
    lidx <- grep("([\\.:]|^)z2$", pi$nm);   if (length(lidx)==1) lambda[i] <- unname(pi$cf[lidx])
    kept_ln_i <- kept_by_eq[[i]][grepl("^ln_", kept_by_eq[[i]])]
    for (lnv in kept_ln_i) {
      j <- ln_to_j(lnv)
      gidx <- grep(paste0("(^|[\\.:])", lnv, "$"), pi$nm)
      if (length(gidx)==1) gamma[i, j] <- unname(pi$cf[gidx])
    }
    gamma[i, 1] <- - sum(gamma[i, -1], na.rm = TRUE)  # base_drop = 1
  }
  
  obj <- list(
    call         = match.call(),
    method       = "QUAIDS",
    priceIndex   = priceIndex,
    estMethod    = estMethod,
    maxiter      = maxiter,
    priceNames   = priceNames,
    shareNames   = shareNames,
    eqLabels     = eqLabels,
    totExpName   = "x",
    coef         = list(alpha = alpha, beta = beta, lambda = lambda, gamma = gamma),
    fit          = fit,
    data         = df,
    resid        = residuals(fit),
    base_dropped = priceNames[1],
    rhs_by_eq    = rhs_list,
    kept_by_eq   = kept_by_eq,
    dropped_by_eq= drop_by_eq
  )
  class(obj) <- c("aidsEst","quaids_manual","list")
  obj
}

summary.quaids_manual <- function(object, ...){
  fit <- object$fit
  cat("Demand analysis with the QUADRATIC Almost Ideal Demand System (QUAIDS)\n")
  cat("Estimation Method:", object$estMethod, "\n")
  cat("Price Index:", object$priceIndex, "\n")
  if (!is.null(object$base_dropped)) {
    cat("Base price dropped from RHS:", object$base_dropped, "\n\n")
  } else {
    if (!is.null(object$drop_price)) cat("Dropped price from RHS:", object$drop_price, "\n")
    if (!is.null(object$omit_share)) cat("Omitted share equation:", object$omit_share, "\n")
    cat("\n")
  }
  
  co <- coef(fit); se <- sqrt(diag(vcov(fit))); t  <- co/se
  p  <- 2*pt(-abs(t), df = fit$df.residual)
  lab <- names(co)
  relabel <- sub("^([^.:]+)[.:]\\(Intercept\\)$", "alpha \\1", lab)
  relabel <- sub("([.:]|^)z2$", "lambda", sub("([.:]|^)z$", "beta", relabel))
  relabel <- sub("([.:]|^)ln_", "gamma ", relabel)
  relabel <- sub("^([^ ]+)[.:]", "\\1 ", relabel)
  
  tab <- data.frame(Estimate = co, `Std. Error` = se, `t value` = t, `Pr(>|t|)` = p,
                    row.names = relabel, check.names = FALSE)
  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE)
  
  cat("\nR-squared Values of expenditure shares:\n")
  eq_list <- try(fit$eq, silent = TRUE)
  if (inherits(eq_list, "try-error") || length(eq_list) == 0) {
    print(setNames(numeric(0), character(0)))
  } else {
    r2 <- vapply(seq_along(eq_list), function(i){
      ei <- eq_list[[i]]
      y    <- model.response(model.frame(ei))
      yhat <- fitted(ei)
      num <- sum((y - yhat)^2, na.rm = TRUE)
      den <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
      if (!is.finite(den) || den <= 0) NA_real_ else 1 - num/den
    }, numeric(1))
    nm <- vapply(eq_list, function(ei){
      lhs <- tryCatch(as.character(formula(ei))[2], error = function(e) NA_character_)
      if (is.na(lhs)) "eq" else lhs
    }, character(1))
    names(r2) <- nm
    print(r2)
  }
  invisible(tab)
}

# -----------------------------------------------------------------------------
# .inst_qr_names
# Constrói matriz de instrumentos a partir de nomes de colunas em 'data',
# expande fatores via model.matrix, NÃO remove linhas com NA (preenche com 0),
# remove colinearidade por QR com pivoteamento e retorna nomes mantidos.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# .inst_qr_names — compatível (var_tol/qr_tol ou tol) e à prova de 'unused args'
# Retorna APENAS os nomes das colunas independentes de Z, já expandidas.
# -----------------------------------------------------------------------------
.inst_qr_names <- function(instNames,
                           data,
                           ...,
                           var_tol       = NULL,   # tolerância p/ variância ~0
                           qr_tol        = NULL,   # tolerância p/ rank no QR
                           add_intercept = FALSE,  # quase sempre FALSE p/ inst
                           na_fill       = 0,
                           center        = FALSE,
                           scale         = FALSE,
                           reorder_keep  = c("original","pivot"),
                           tol           = NULL) { # compat: usa 'tol' se passado
  reorder_keep <- match.arg(reorder_keep)
  
  if (length(instNames) == 0L)
    stop("instNames vazio.", call. = FALSE)
  
  instNames <- unique(instNames)
  
  # interseção com as colunas existentes
  miss <- setdiff(instNames, colnames(data))
  if (length(miss)) {
    warning("Variáveis ausentes em 'data' serão ignoradas: ",
            paste(miss, collapse = ", "))
    instNames <- intersect(instNames, colnames(data))
  }
  if (length(instNames) == 0L)
    stop("Nenhuma variável de instrumento restante em 'data'.", call. = FALSE)
  
  # defaults de tolerância: prioriza var_tol/qr_tol; cai para 'tol' se houver; senão usa padrões robustos
  if (is.null(var_tol)) var_tol <- if (!is.null(tol)) tol else 1e-12
  if (is.null(qr_tol))  qr_tol  <- if (!is.null(tol)) tol else 1e-10
  
  # fórmula RHS usando reformulate (lida com nomes automaticamente)
  fml <- stats::reformulate(instNames, intercept = isTRUE(add_intercept))
  
  # expande sem dropar linhas com NA
  mm <- stats::model.matrix(fml, data = data, na.action = stats::na.pass)
  storage.mode(mm) <- "double"
  
  # substitui NA/Inf
  bad <- !is.finite(mm)
  if (any(bad)) mm[bad] <- na_fill
  
  # remove quase-constantes
  vv <- apply(mm, 2L, function(v) stats::var(v, na.rm = TRUE))
  keep_var <- is.finite(vv) & (vv > var_tol)
  if (!all(keep_var)) {
    drop0 <- colnames(mm)[!keep_var]
    if (length(drop0))
      message("Colunas quase-constantes removidas: ", paste(drop0, collapse = ", "))
    mm <- mm[, keep_var, drop = FALSE]
  }
  if (ncol(mm) == 0L)
    stop("Todas as colunas foram descartadas por variância ~0.", call. = FALSE)
  
  # padronização opcional
  if (isTRUE(center) || isTRUE(scale)) {
    mm <- scale(mm, center = center, scale = scale)
    mm[!is.finite(mm)] <- 0
  }
  
  # QR com pivoteamento
  qrZ <- qr(mm, tol = qr_tol)
  r   <- qrZ$rank
  if (is.na(r) || r <= 0)
    stop("Rank(Z)=0 após QR; ajuste 'qr_tol' ou revise seus IVs.", call. = FALSE)
  
  piv <- qrZ$pivot
  idx_keep <- piv[seq_len(r)]
  if (reorder_keep == "original") idx_keep <- sort(idx_keep)
  
  kept_names <- colnames(mm)[idx_keep]
  kept_names
}


fit_quaids_manual_km1 <- function(prices, shares, x,
                                  priceIndex = c("S","Ls"),
                                  estMethod  = c("SUR","3SLS"),
                                  baseShares = NULL,
                                  instNames  = NULL,
                                  maxiter    = 200,
                                  omit_share = 1,
                                  drop_price = 1,
                                  instData   = NULL,
                                  use_z2     = TRUE){
  stopifnot(is.data.frame(prices), is.data.frame(shares), length(x) == nrow(prices))
  stopifnot(nrow(prices) == nrow(shares))
  priceIndex <- match.arg(priceIndex)
  estMethod  <- match.arg(estMethod)
  
  priceNames <- colnames(prices)
  shareNames <- colnames(shares)
  K <- ncol(prices)
  if (!(omit_share %in% seq_len(K))) stop("omit_share deve estar entre 1..K")
  if (!(drop_price %in% seq_len(K))) stop("drop_price deve estar entre 1..K")
  
  as_num <- function(v) suppressWarnings(as.numeric(v))
  prices[] <- lapply(prices, as_num)
  shares[] <- lapply(shares, as_num)
  x <- as_num(x)
  
  ok <- rowSums(sapply(prices, function(v) is.finite(v) & v > 0)) == K & is.finite(x) & x > 0
  P <- prices[ok, , drop = FALSE]
  W <- shares[ok, , drop = FALSE]
  X <- x[ok]
  
  lnP <- switch(priceIndex,
                "Ls" = {
                  wbar <- if (is.null(baseShares)) colMeans(W, na.rm = TRUE) else as.numeric(baseShares)
                  if (length(wbar) != K) stop("baseShares deve ter K elementos.")
                  wbar <- wbar / sum(wbar)
                  as.numeric(as.matrix(log(P)) %*% wbar)
                },
                "S"  = {
                  w_i <- W/rowSums(W)
                  rowSums(w_i * log(P))
                })
  z  <- log(X) - lnP
  z2 <- z^2
  
  lnPcols <- paste0("ln_", priceNames)
  df <- data.frame(W, setNames(as.data.frame(log(P)), lnPcols),
                   z = z, z2 = z2, check.names = FALSE)
  
  # Equações (definidas ANTES da estimação)
  keep_sh_idx <- setdiff(seq_len(K), omit_share)
  keep_ln_idx <- setdiff(seq_len(K), drop_price)
  share_km1 <- shareNames[keep_sh_idx]
  ln_km1    <- lnPcols[keep_ln_idx]
  
  rhs_terms <- c(ln_km1, "z", if (use_z2) "z2" else NULL)
  eqLabels <- gsub("[ _]", "", share_km1)
  eqs <- vector("list", length(share_km1)); names(eqs) <- eqLabels
  for (i in seq_along(share_km1)) eqs[[i]] <- stats::reformulate(rhs_terms, response = share_km1[i])
  
  # Estimação
  if (estMethod == "SUR") {
    fit <- systemfit::systemfit(eqs, data = df, method = "SUR", maxit = maxiter)
  } else {
    # ------------- 3SLS robusto a IVs fora da base e colinearidade em Z -------------
    if (is.null(instNames) || length(instNames) == 0) stop("Sem instrumentos (instNames) fornecidos.")
    instNames <- unique(instNames[!is.na(instNames) & nzchar(instNames)])
    
    in_df  <- intersect(instNames, names(df))
    miss   <- setdiff(instNames, in_df)
    
    if (length(miss) > 0) {
      if (is.null(instData)) stop("Para 3SLS, passe instData com as colunas de instrumentos.")
      avail <- intersect(miss, names(instData))
      if (length(avail) == 0L) {
        ex_iv <- head(grep("^(iv_|IV_d_)", names(instData), value = TRUE), 12)
        stop(
          sprintf(
            "Sem instrumentos após interseção com data.\nSolicitados (ex.): %s\nDisponíveis que parecem IVs (ex.): %s",
            paste(head(instNames, 12), collapse = ", "),
            if (length(ex_iv)) paste(ex_iv, collapse = ", ") else "<nenhum padrão ^iv_|^IV_d_ encontrado>"
          )
        )
      }
      df <- cbind(df, as.data.frame(instData[, avail, drop = FALSE]))
      in_df <- c(in_df, avail)
    }
    
    # Afina por QR: só termos linearmente independentes entram na fórmula
    inst_ok <- .inst_qr_names(in_df, data = df, var_tol = 1e-12, qr_tol = 1e-10)
    inst_fml <- stats::reformulate(inst_ok, intercept = FALSE)
    
    # Tenta rodar; se solver reclamar de singularidade, afrouxa o tol e afina de novo
    fit <- try(
      systemfit::systemfit(eqs, data = df, method = "3SLS", inst = inst_fml, maxit = maxiter),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      inst_ok2 <- .inst_qr_names(in_df, data = df, var_tol = 1e-12, qr_tol = 1e-8)
      inst_fml2 <- stats::reformulate(inst_ok2, intercept = FALSE)
      fit <- systemfit::systemfit(eqs, data = df, method = "3SLS", inst = inst_fml2, maxit = maxiter)
    }
  }
  
  # Coeficientes e reconstrução
  cf <- coef(fit); se <- sqrt(diag(vcov(fit)))
  pick_eq <- function(eqLabel){ nm <- names(cf); idx <- grep(paste0("^", eqLabel, "[\\.:]"), nm); list(cf = cf[idx], se = se[idx], nm = nm[idx]) }
  
  alpha <- setNames(numeric(K), shareNames)
  beta  <- setNames(numeric(K), shareNames)
  lambda<- setNames(numeric(K), shareNames)
  gamma <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  
  map_ln_j <- setNames(match(sub("^ln_", "", ln_km1), priceNames), ln_km1)
  for (i in seq_along(share_km1)) {
    eqLab <- eqLabels[i]
    pi <- pick_eq(eqLab)
    aidx <- grep("\\(Intercept\\)", pi$nm); if (length(aidx) == 1) alpha[ keep_sh_idx[i] ] <- unname(pi$cf[aidx])
    bidx <- grep("([\\.:]|^)z$",  pi$nm);   if (length(bidx) == 1) beta[  keep_sh_idx[i] ] <- unname(pi$cf[bidx])
    lidx <- grep("([\\.:]|^)z2$", pi$nm);   if (length(lidx) == 1) lambda[keep_sh_idx[i] ] <- unname(pi$cf[lidx])
    for (lnv in names(map_ln_j)) {
      j <- map_ln_j[[lnv]]
      gidx <- grep(paste0("(^|[\\.:])", lnv, "$"), pi$nm)
      if (length(gidx) == 1) gamma[ keep_sh_idx[i], j ] <- unname(pi$cf[gidx])
    }
  }
  alpha[omit_share]  <- 1 - sum(alpha[keep_sh_idx],  na.rm = TRUE)
  beta[omit_share]   <- -   sum(beta[ keep_sh_idx],  na.rm = TRUE)
  lambda[omit_share] <- -   sum(lambda[keep_sh_idx], na.rm = TRUE)
  for (j in seq_len(K)) gamma[omit_share, j] <- - sum(gamma[keep_sh_idx, j], na.rm = TRUE)
  adj <- sum(gamma[omit_share, ], na.rm = TRUE)
  if (is.finite(adj) && abs(adj) > 1e-10) gamma[omit_share, ] <- gamma[omit_share, ] - adj/K
  
  obj <- list(
    call         = match.call(), method = "QUAIDS",
    priceIndex   = priceIndex,   estMethod = estMethod, maxiter = maxiter,
    priceNames   = priceNames,   shareNames = shareNames,
    omit_share   = shareNames[omit_share], drop_price = priceNames[drop_price],
    coef         = list(alpha = alpha, beta = beta, lambda = lambda, gamma = gamma),
    fit          = fit, data = df, resid = residuals(fit), rhs = paste(rhs_terms, collapse = " + ")
  )
  class(obj) <- c("quaids_manual","aidsEst","list")
  obj
}

# Elasticidades QUAIDS manual
.enforce_hicks_adding <- function(EH, w){ ww <- sum(w^2); if (!is.finite(ww) || ww <= 0) return(EH); adj <- as.numeric(EH %*% w) / ww; EH - adj %o% w }
.symmetrize_slutsky   <- function(EH, w){ if (any(!is.finite(w)) || any(w <= 0)) return(EH); W <- diag(w); M <- W %*% EH; Ms <- (M + t(M))/2; solve(W, Ms) }
.enforce_cournot_marshall <- function(EM, eta, w){ ww <- sum(w^2); if (!is.finite(ww) || ww <= 0) return(EM); target <- -(eta - 1); cur <- as.numeric(EM %*% w); a <- (target - cur) / ww; EM + a %o% w }
.check_identities <- function(EH, EM, eta, w){ add_hicks <- as.numeric(EH %*% w); cour_mar <- cbind(lhs = as.numeric(EM %*% w), rhs = -(eta - 1), diff = as.numeric(EM %*% w) + (eta - 1)); W <- diag(w); symm_err <- max(abs(W %*% EH - t(W %*% EH))); list(adding_hicks = add_hicks, cournot_marshall = cour_mar, symmetry_err = symm_err) }

# ========= PATCH: elas_quaids_manual robusto a NA/Inf nos coeficientes ========
# ====== PATCH: elas_quaids_manual com shares > 0 e soma=1 ======
# ====== PATCH: elas_quaids_manual com passo ajustável e diagnósticos ======
# --- SUBSTITUIR a versão atual ---

# ---------------------------------------------------------------
# Preditor de shares para QUAIDS (compatível com o seu script)
# ---------------------------------------------------------------
predict_shares_quaids <- function(fit_q, x, p = NULL, eps = 1e-9,
                                  normalize = c("clip","softmax","none"),
                                  tau = 1) {
  stopifnot(inherits(fit_q, "quaids_manual"))
  normalize <- match.arg(normalize)
  
  priceNames <- fit_q$priceNames
  shareNames <- fit_q$shareNames
  K <- length(priceNames)
  
  # Coeficientes (robustez a NA/Inf)
  alpha  <- fit_q$coef$alpha;   alpha[!is.finite(alpha)]   <- 0
  beta   <- fit_q$coef$beta;    beta[!is.finite(beta)]     <- 0
  lambda <- fit_q$coef$lambda;  lambda[!is.finite(lambda)] <- 0
  gamma  <- fit_q$coef$gamma;   gamma[!is.finite(gamma)]   <- 0
  
  # Preços-alvo (se p == NULL, usa média dos ln preços da amostra do fit)
  if (is.null(p)) {
    ln_cols <- paste0("ln_", priceNames)
    lv <- colMeans(fit_q$data[, ln_cols, drop = FALSE], na.rm = TRUE)
    lv[!is.finite(lv)] <- median(as.numeric(as.matrix(fit_q$data[, ln_cols, drop = FALSE])), na.rm = TRUE)
    p <- exp(lv)
  }
  p <- pmax(as.numeric(p), eps); names(p) <- priceNames
  lp <- log(p)
  
  # Pesos médios de base (para índice Ls)
  wbar <- colMeans(fit_q$data[, shareNames, drop = FALSE], na.rm = TRUE)
  wbar <- wbar / sum(wbar)
  
  # Índice de preços e índice linear s_i
  if (identical(fit_q$priceIndex, "Ls")) {
    z <- log(x) - sum(wbar * lp)
    s <- as.numeric(alpha + gamma %*% lp + beta * z + lambda * z^2)
  } else if (identical(fit_q$priceIndex, "S")) {
    # passo 1: z com SL (wbar) e pré-share
    z_ls <- log(x) - sum(wbar * lp)
    s1   <- as.numeric(alpha + gamma %*% lp + beta * z_ls + lambda * z_ls^2)
    w1 <- switch(
      normalize,
      "clip"    = { w <- pmax(s1, 1e-6); w/sum(w) },
      "softmax" = { ex <- exp(s1/tau - max(s1/tau)); ex/sum(ex) },
      "none"    = s1
    )
    # passo 2: z com Stone (w1)
    z <- log(x) - sum(w1 * lp)
    s <- as.numeric(alpha + gamma %*% lp + beta * z + lambda * z^2)
  } else {
    stop("Índice de preço desconhecido em fit_q$priceIndex (use 'Ls' ou 'S').")
  }
  
  # Normalização para obter shares
  w <- switch(
    normalize,
    "clip"    = { w <- pmax(s, 1e-6); w/sum(w) },
    "softmax" = { ex <- exp(s/tau - max(s/tau)); ex/sum(ex) },
    "none"    = s
  )
  w[!is.finite(w)] <- 0
  if (normalize == "none") {
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) w <- rep(1/K, K) else w <- w / sw
  }
  names(w) <- shareNames
  w
}


elas_quaids_manual <- function(fit_q, x, p = NULL,
                               dlogp = 1e-3,          # passo inicial em log-preço
                               dlogp_max = 5e-2,      # teto p/ passo adaptativo
                               target_dw = 5e-4,      # alvo de variação de share
                               normalize_eval  = c("clip","softmax","none"),
                               normalize_deriv = c("softmax","clip","none"),
                               tau = 1,
                               enforce_shares = TRUE, # mantém compatibilidade da checagem
                               enforce_hicks = TRUE,
                               enforce_symmetry = FALSE,
                               enforce_cournot = FALSE,
                               return_checks = TRUE,
                               eps = 1e-8, eps_w = 1e-6){
  
  stopifnot(inherits(fit_q, "quaids_manual"))
  priceNames <- fit_q$priceNames; shareNames <- fit_q$shareNames; K <- length(priceNames)
  normalize_eval  <- match.arg(normalize_eval)
  normalize_deriv <- match.arg(normalize_deriv)
  
  if (is.null(p)) {
    ln_cols <- paste0("ln_", priceNames)
    lv <- colMeans(fit_q$data[, ln_cols, drop = FALSE], na.rm = TRUE)
    lv[!is.finite(lv)] <- median(as.numeric(as.matrix(fit_q$data[, ln_cols, drop = FALSE])), na.rm = TRUE)
    p <- exp(lv)
  }
  p <- pmax(as.numeric(p), eps); names(p) <- priceNames
  lp <- log(p)
  
  # shares "para exibir/checar"
  w_eval <- predict_shares_quaids(fit_q, x = x, p = p, normalize = normalize_eval,  tau = tau)
  # shares "para derivar" (regularizados para evitar canto)
  w0     <- predict_shares_quaids(fit_q, x = x, p = p, normalize = normalize_deriv, tau = tau)
  
  if (any(!is.finite(w0)) || any(w0 <= 0)) stop("Shares previstos não-positivos após normalização derivativa.")
  
  # z para o índice escolhido (usa w0 se priceIndex = "S")
  wbar <- colMeans(fit_q$data[, shareNames, drop = FALSE], na.rm = TRUE); wbar <- wbar/sum(wbar)
  z0 <- if (identical(fit_q$priceIndex, "Ls")) sum(log(x)) - sum(wbar * lp) else log(x) - sum(w0 * lp)
  
  alpha  <- fit_q$coef$alpha;   alpha[!is.finite(alpha)]   <- 0
  beta   <- fit_q$coef$beta;    beta[!is.finite(beta)]     <- 0
  lambda <- fit_q$coef$lambda;  lambda[!is.finite(lambda)] <- 0
  gamma  <- fit_q$coef$gamma;   gamma[!is.finite(gamma)]   <- 0
  
  eta <- 1 + (beta + 2*lambda*z0) / w0; names(eta) <- shareNames
  
  q0 <- (w0 * x) / p
  E_M <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  
  # passo adaptativo: aumenta dlogp até mexer no mínimo target_dw
  step0 <- dlogp
  for (j in seq_len(K)) {
    step <- step0
    repeat {
      p_up <- p; p_up[j] <- p_up[j] * exp(step)
      w_up <- predict_shares_quaids(fit_q, x = x, p = p_up, normalize = normalize_deriv, tau = tau)
      if (max(abs(w_up - w0)) >= target_dw || step >= dlogp_max) break
      step <- step * 2
    }
    q_up <- (w_up * x) / p_up
    E_M[, j] <- (log(q_up) - log(q0)) / step
  }
  
  E_H <- E_M + outer(eta, w0)
  
  .enforce_hicks_adding <- function(EH, w){
    ww <- sum(w^2); if (!is.finite(ww) || ww <= 0) return(EH)
    adj <- as.numeric(EH %*% w) / ww; EH - adj %o% w
  }
  .symmetrize_slutsky <- function(EH, w){
    if (any(!is.finite(w)) || any(w <= 0)) return(EH)
    W <- diag(w); M <- W %*% EH; Ms <- (M + t(M))/2; solve(W, Ms)
  }
  .enforce_cournot_marshall <- function(EM, eta, w){
    ww <- sum(w^2); if (!is.finite(ww) || ww <= 0) return(EM)
    target <- -(eta - 1); cur <- as.numeric(EM %*% w); a <- (target - cur) / ww
    EM + a %o% w
  }
  .check_identities <- function(EH, EM, eta, w){
    add_hicks <- as.numeric(EH %*% w)
    cour_mar  <- cbind(lhs = as.numeric(EM %*% w), rhs = -(eta - 1),
                       diff = as.numeric(EM %*% w) + (eta - 1))
    W <- diag(w); symm_err <- max(abs(W %*% EH - t(W %*% EH)))
    S  <- (W %*% EH + t(W %*% EH))/2
    eig_max <- tryCatch(max(eigen(S, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
    list(adding_hicks = add_hicks, cournot_marshall = cour_mar,
         symmetry_err = symm_err, slutsky_max_eig = eig_max)
  }
  
  if (enforce_symmetry) E_H <- .symmetrize_slutsky(E_H, w0)
  if (enforce_hicks)    E_H <- .enforce_hicks_adding(E_H, w0)
  E_M <- E_H - outer(eta, w0)
  if (enforce_cournot) { E_M <- .enforce_cournot_marshall(E_M, eta, w0); E_H <- E_M + outer(eta, w0) }
  
  out <- list(at = list(p = setNames(p, priceNames), x = x,
                        w = w_eval, w_deriv = w0, z = z0),
              expenditure = eta, marshall = E_M, hicks = E_H)
  if (return_checks) out$checks <- .check_identities(E_H, E_M, eta, w0)
  out
}

# ================== fim do PATCH ==================

# ================== fim do PATCH ==================

# ================== fim do PATCH ==================


# ===============================================================
# Fallback escalonado 3SLS (core → mid → full; com/sem z2)
# ===============================================================
safe_fit_quaids_3sls_any <- function(inst_sets, data_Z, ..., use_z2_first = TRUE){
  last_err <- NULL
  req_endog_local <- tryCatch(get("req_endog", inherits = TRUE), error = function(e) NULL)
  if (is.null(req_endog_local)) req_endog_local <- 1L
  
  for (k in seq_along(inst_sets)) {
    S <- inst_sets[[k]]
    S0 <- unique(S)
    S0 <- S0[!is.na(S0) & nzchar(S0)]
    S0 <- intersect(S0, names(data_Z))
    if (length(S0) == 0L) { last_err <- sprintf("Conjunto de IVs #%d ficou vazio após interseção com a base.", k); next }
    
    try_once <- function(use_z2_flag){
      try(
        fit_quaids_manual_km1(instNames = S0, instData = data_Z, use_z2 = use_z2_flag, ...),
        silent = TRUE
      )
    }
    ans <- if (use_z2_first) try_once(TRUE) else try_once(FALSE)
    if (!inherits(ans, "try-error")) return(ans)
    last_err <- ans
    ans2 <- if (use_z2_first) try_once(FALSE) else try_once(TRUE)
    if (!inherits(ans2, "try-error")) return(ans2)
    last_err <- ans2
  }
  stop(last_err)
}

# ===============================================================
# Execução 3SLS (três opções + seleção)
# ===============================================================

Z_chk <- model.matrix(stats::reformulate(iv_all_candidates, intercept = FALSE), data = df_iv)
cat("Rank(Z) final (cheque):", qr(Z_chk)$rank, "\n")

run_quaids_3sls <- function(instNames, label, use_z2 = TRUE,
                            tol_grid = c(1e-10, 1e-8, 1e-6, 1e-4)) {
  cat("\n---- Tentando set:", label, " (|IV|=", length(instNames),
      ", use_z2=", use_z2, ") ----\n", sep = "")
  # Intersecta com df_iv e afina por QR (tol crescente até rodar)
  base_vars <- intersect(unique(instNames), names(df_iv))
  if (length(base_vars) == 0L) {
    cat("Falhou. Erro:\n", "Nenhum instrumento do set existe na base.\n", sep = "")
    return(NULL)
  }
  for (tol in tol_grid) {
    ok_vars <- try(.inst_qr_names(base_vars, data = df_iv, qr_tol = tol), silent = TRUE)
    if (inherits(ok_vars, "try-error")) next
    cat("  • tol=", format(tol, scientific = TRUE),
        " ⇒ IVs efetivos: ", length(ok_vars), "\n", sep = "")
    ans <- try(fit_quaids_manual_km1(
      prices     = df_iv[, priceNames, drop = FALSE],
      shares     = df_iv[, shareNames, drop = FALSE],
      x          = df_iv[["gasto_total_atualhat"]],
      priceIndex = "Ls",
      estMethod  = "3SLS",
      omit_share = 1,
      drop_price = 1,
      instNames  = ok_vars,      # << usa só colunas independentes
      instData   = df_iv,
      use_z2     = use_z2,
      maxiter    = 500
    ), silent = TRUE)
    if (!inherits(ans, "try-error")) return(ans)
    cat("    ↳ ainda falhou (", as.character(ans), ")\n", sep = "")
  }
  cat("Falhou em todos os tolerances.\n")
  NULL
}


# Tenta core/mid/full com z2=TRUE; se falhar, tenta sem z2
fit_core <- run_quaids_3sls(iv_set_core, "core", use_z2 = TRUE)
if (is.null(fit_core)) fit_core <- run_quaids_3sls(iv_set_core, "core", use_z2 = FALSE)

fit_mid  <- run_quaids_3sls(iv_set_mid,  "mid",  use_z2 = TRUE)
if (is.null(fit_mid))  fit_mid  <- run_quaids_3sls(iv_set_mid,  "mid",  use_z2 = FALSE)

fit_full <- run_quaids_3sls(iv_set_full, "full", use_z2 = TRUE)
if (is.null(fit_full)) fit_full <- run_quaids_3sls(iv_set_full, "full", use_z2 = FALSE)

# Seleciona prioridade: mid → full → core
fit_q_3sls <- if (!is.null(fit_mid)) fit_mid else if (!is.null(fit_full)) fit_full else fit_core
if (is.null(fit_q_3sls)) {
  cat("\nNenhum set (core/mid/full) rodou direto — chamando fallback escalonado...\n")
  fit_q_3sls <- safe_fit_quaids_3sls_any(
    inst_sets   = list(iv_set_mid, iv_set_full, iv_set_core),
    data_Z      = df_iv,
    prices      = df_iv[, priceNames, drop = FALSE],
    shares      = df_iv[, shareNames, drop = FALSE],
    x           = df_iv[["gasto_total_atualhat"]],
    priceIndex  = "Ls",
    estMethod   = "3SLS",
    omit_share  = 1,
    drop_price  = 1,
    maxiter     = 500,
    use_z2_first = TRUE
  )
}

cat("\n==== QUAIDS 3SLS selecionado — summary ====\n"); print(summary(fit_q_3sls))

# (Opcional) também imprimir os que rodaram
if (!is.null(fit_core)) { cat("\n==== QUAIDS 3SLS (core) — summary ====\n"); print(summary(fit_core)) }
if (!is.null(fit_mid))  { cat("\n==== QUAIDS 3SLS (mid)  — summary ====\n"); print(summary(fit_mid))  }
if (!is.null(fit_full)) { cat("\n==== QUAIDS 3SLS (full) — summary ====\n"); print(summary(fit_full)) }

# (Opcional) Elasticidades ao ponto médio
# elas_q <- elas_quaids_manual(fit_q_3sls, x = mean(df_iv$gasto_total_atualhat, na.rm=TRUE))
# str(elas_q)

r2_alt <- function(fit){
  eqs <- try(fit$fit$eq, silent = TRUE)
  if (inherits(eqs, "try-error") || length(eqs) == 0) return(setNames(numeric(0), character(0)))
  vapply(eqs, function(ei){
    y    <- model.response(model.frame(ei))
    yhat <- fitted(ei)
    if (var(y, na.rm = TRUE) <= 0) NA_real_ else suppressWarnings(cor(y, yhat, use = "complete.obs")^2)
  }, numeric(1))
}
round(r2_alt(fit_full), 3)

# R² por equação (mesma métrica do summary)
r2_per_eq <- function(fit){
  eqs <- try(fit$fit$eq, silent = TRUE)
  if (inherits(eqs, "try-error") || length(eqs) == 0) return(setNames(numeric(0), character(0)))
  out <- vapply(eqs, function(ei){
    y    <- model.response(model.frame(ei))
    yhat <- fitted(ei)
    den <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    if (!is.finite(den) || den <= 0) NA_real_ else 1 - sum((y - yhat)^2, na.rm = TRUE)/den
  }, numeric(1))
  out
}
mean_r2 <- function(fit) { r <- r2_per_eq(fit); if (!length(r)) NA_real_ else mean(r, na.rm = TRUE) }

cands <- list(core = fit_core, mid = fit_mid, full = fit_full)
scores <- sapply(cands, function(f) if (is.null(f)) NA_real_ else mean_r2(f))
print(round(scores, 3))

best <- names(which.max(scores))
fit_q_3sls <- cands[[best]]
cat("Selecionado:", best, " | média R² =", round(scores[best], 3), "\n")

# summary do escolhido (sem duplicar):
summary(fit_q_3sls)

elas_q <- elas_quaids_manual(
  fit_q_3sls,
  x = mean(df_iv$gasto_total_atualhat, na.rm = TRUE),
  dlogp = 1e-3,              # teste 0.1%; se ficar muito “travado”, tente 5e-3
  enforce_hicks = TRUE,
  enforce_symmetry = TRUE
)

cat("shares previstos: min=", min(elas_q$at$w), " max=", max(elas_q$at$w), " sum=", sum(elas_q$at$w), "\n")
cat("own-price (Marshall) diag:\n"); print(round(diag(elas_q$marshall), 3))
cat("mean cross-price (Marshall) off-diag:\n"); print(round(mean(elas_q$marshall[row(elas_q$marshall)!=col(elas_q$marshall)]), 4))
cat("expenditure elasticities (eta):\n"); print(round(elas_q$expenditure, 3))

# ============================================================
# create_base_rob(): prepara base com ln(preços), índices Ls/S,
#                    z_Ls, z_S, z2_* e metadados úteis.
# ------------------------------------------------------------
# Args:
#   df           : data.frame com shares e preços
#   priceNames   : nomes dos preços (ex.: c("preco1","preco2",...))
#   shareNames   : nomes das shares (ex.: c("w1","w2",...))
#   income_var   : nome da coluna do log da despesa (ex.: "lndespesa")
#   fallback_inc : fallback para log renda, se income_var não existir
#   eps          : piso numérico para shares
# Return: data.frame (base_rob) com:
#   ln_<preco>, lnP_Ls, lnP_S, z_Ls, z_S, z2_Ls, z2_S
#   attr(., "Wbar") com w̄ usado no Ls
# ============================================================
## -------------------------------------------------------------
## Cria o objeto-base usado por fit_km1()  (base_rob)
## -------------------------------------------------------------
# exemplo esquemático — use o seu mapeamento real
inst_map <- list(
  core = c("iv_op01","iv_op02", "iv_hg20","iv_hg40","iv_hg60","iv_hg80"),
  mid  = c("iv_bs01","iv_bs02","iv_ns01","iv_ns02","IV_d_reg_X1","IV_d_cap_X1"),
  full_clean = setdiff(names(df_iv), c(priceNames, shareNames, "z","z2"))
)

make_base_rob <- function(df, priceNames, shareNames, eps = 1e-12) {
  stopifnot(all(shareNames %in% names(df)))
  ## garante ln_<preço>
  for (p in priceNames) {
    ln <- paste0("ln_", p)
    if (!ln %in% names(df) && p %in% names(df)) {
      df[[ln]] <- log(df[[p]])
    }
  }
  lnpriceNames <- paste0("ln_", priceNames)
  
  ## pesos médios (wbar) para índice Ls e para montar z
  wbar <- colMeans(df[, shareNames, drop = FALSE], na.rm = TRUE)
  wbar <- pmax(wbar, eps); wbar <- wbar / sum(wbar)
  
  ## escolhe a base de despesa: lndespesa ou ln_income
  pick_base_ln <- function(d) {
    if ("lndespesa" %in% names(d)) return(d$lndespesa)
    if ("ln_income" %in% names(d)) return(d$ln_income)
    stop("Nem 'lndespesa' nem 'ln_income' estão no data.frame.")
  }
  
  ## calcula z (= ln x - sum(wbar ln p)) e z2
  compute_z <- function(d) {
    for (p in priceNames) {
      ln <- paste0("ln_", p)
      if (!ln %in% names(d) && p %in% names(d)) d[[ln]] <- log(d[[p]])
    }
    LN <- as.matrix(d[, paste0("ln_", priceNames), drop = FALSE])
    lnP_Ls <- rowSums(LN * matrix(as.numeric(wbar), nrow(LN), length(wbar), byrow = TRUE))
    z <- pick_base_ln(d) - lnP_Ls
    d$z  <- z
    d$z2 <- z^2
    d
  }
  
  ## índice de Stone (usa shares observadas da própria linha)
  stone_index <- function(d) {
    LN <- as.matrix(d[, paste0("ln_", priceNames), drop = FALSE])
    W  <- as.matrix(d[, shareNames, drop = FALSE])
    rowSums(LN * W)
  }
  
  list(
    shareNames    = shareNames,
    priceNames    = priceNames,
    lnpriceNames  = lnpriceNames,
    wbar          = wbar,
    compute_z     = compute_z,   # devolve d com z e z2
    stone_index   = stone_index, # ln P^S por linha (opcional em diagnósticos)
    eps           = eps
  )
}

# ============================================================
# fit_km1(): ajusta AIDS/QUAIDS por SUR (fallback 3SLS opcional)
# ------------------------------------------------------------
# Args:
#   model      : "AIDS" | "QUAIDS"
#   index      : "Ls" | "S"   (escolhe z_* e z2_* da base)
#   inst_try   : vetor de chaves do seu 'inst_map' (ex. c("full_clean","mid","core"))
#   base       : data.frame (use create_base_rob() antes)
#   use_z2     : se TRUE, inclui z2_* nas equações (para QUAIDS)
#   try_3sls   : se TRUE, tenta 3SLS com instrumentos e cai para SUR se falhar
# Return:
#   list(fit = <systemfit>, shareNames = ..., rhs = ..., class = "sur_fallback")
#   com attr(., "model") = "AIDS"/"QUAIDS" e attr(., "index") = "Ls"/"S"
# ============================================================
## -------------------------------------------------------------
## Ajusta SUR para AIDS/QUAIDS usando apenas info do `base`
## (sem get0/variáveis globais).
## -------------------------------------------------------------
fit_km1 <- function(model = c("AIDS", "QUAIDS"),
                    index  = c("Ls", "S"),
                    base,
                    inst_try = c("full_clean","mid","core"),  # aceito, mas não usado aqui
                    use_z2 = FALSE,
                    df = NULL) {
  
  stopifnot(is.list(base),
            !is.null(base$priceNames),
            !is.null(base$shareNames),
            !is.null(base$compute_z))
  
  model <- match.arg(model)
  index <- match.arg(index)
  
  # dados de entrada
  dat <- if (!is.null(df)) {
    df
  } else if (exists("df_iv", inherits = TRUE)) {
    get("df_iv", inherits = TRUE)
  } else {
    stop("Forneça os dados via `df=` ou defina `df_iv` no ambiente.")
  }
  
  # garante ln_<preço>
  for (p in base$priceNames) {
    ln <- paste0("ln_", p)
    if (!ln %in% names(dat) && p %in% names(dat)) {
      dat[[ln]] <- log(dat[[p]])
    }
  }
  
  # z e z2 conforme a spec do base
  dat <- base$compute_z(dat)
  
  # monta RHS
  rhs_prices <- paste(paste0("ln_", base$priceNames), collapse = " + ")
  rhs <- paste(rhs_prices, "+ z")
  if (model == "QUAIDS" && isTRUE(use_z2)) {
    rhs <- paste(rhs, "+ z2")
  }
  
  # sistema com k-1 equações (evita soma=1)
  k <- length(base$shareNames)
  eqs <- vector("list", k - 1)
  names(eqs) <- paste0("w", 1:(k - 1))
  for (i in 1:(k - 1)) {
    lhs <- base$shareNames[i]
    eqs[[i]] <- stats::as.formula(paste(lhs, "~", rhs))
  }
  
  fit <- try(systemfit::systemfit(eqs, method = "SUR", data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(fit)
  
  # metadados para o pipeline/predict_shares_generic()
  attr(fit, "model")      <- toupper(model)   # "AIDS" ou "QUAIDS"
  attr(fit, "index")      <- toupper(index)   # "LS" ou "S"
  attr(fit, "shareNames") <- base$shareNames
  attr(fit, "priceNames") <- base$priceNames
  attr(fit, "wbar")       <- base$wbar
  
  fit
}

# 1) montar a base_rob
base_rob <- make_base_rob(df_iv, priceNames, shareNames)

# ===============================================================
# PIPELINE 1–9: AIDS vs QUAIDS, diagnósticos IV, normalização,
# base/reescalas, seleção via holdout, SUR/Hausman, checagens
# microeconômicas, índice Ls x S, higiene de dados.
# v4.1 — CV por fold + Wald(λ) + Slutsky NSD
# ===============================================================
# ============================================================
# Pipeline AIDS vs QUAIDS c/ SUR, CV, holdout e diagnósticos IV
# ============================================================

if (!requireNamespace("systemfit", quietly = TRUE)) {
  stop("Pacote 'systemfit' é necessário. Instale e tente novamente.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
set.seed(123)

# ------------------------------------------------------------
# Fallbacks leves (só ativados se ausentes no ambiente)
# ------------------------------------------------------------

# make_base_rob() mínimo (se você já tem, o seu será usado)
if (!exists("make_base_rob")) {
  make_base_rob <- function(df, priceNames, shareNames) {
    compute_z <- function(dat) {
      # garante ln_<preço>
      for (p in priceNames) {
        ln <- paste0("ln_", p)
        if (!ln %in% names(dat) && p %in% names(dat)) {
          dat[[ln]] <- log(pmax(as.numeric(dat[[p]]), 1e-12))
        }
      }
      lnpriceNames <- paste0("ln_", priceNames)
      base_ln <- if ("lndespesa" %in% names(dat)) dat$lndespesa else dat$ln_income
      if (is.null(base_ln)) stop("Nem 'lndespesa' nem 'ln_income' encontrados p/ montar 'z'.")
      Wbar <- colMeans(dat[, shareNames, drop = FALSE], na.rm = TRUE)
      Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
      LN   <- as.matrix(dat[, lnpriceNames, drop = FALSE])
      lnP  <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
      dat$z  <- as.numeric(base_ln) - lnP
      dat$z2 <- dat$z^2
      dat
    }
    # wbar global da base
    df2  <- compute_z(df)
    wbar <- colMeans(df2[, shareNames, drop = FALSE], na.rm = TRUE)
    wbar <- pmax(wbar, 1e-12); wbar <- wbar / sum(wbar)
    list(priceNames = priceNames, shareNames = shareNames, compute_z = compute_z, wbar = wbar)
  }
}

# inst_map mínimo (se você já tem, o seu será usado)
if (!exists("inst_map")) {
  iv_cols <- grep("^(iv_|IV_)", names(df_iv), value = TRUE)
  inst_map <- list(
    core = iv_cols,
    mid  = iv_cols
  )
}
inst_try_default <- c("mid","core")

# ------------------------------------------------------------
# fit_km1(): SUR para AIDS/QUAIDS usando base_rob
# ------------------------------------------------------------
fit_km1 <- function(model = c("AIDS", "QUAIDS"),
                    index  = c("Ls", "S"),
                    base,
                    inst_try = c("full_clean","mid","core"),  # aceito, mas não usado aqui
                    use_z2 = FALSE,
                    df = NULL) {
  stopifnot(is.list(base),
            !is.null(base$priceNames),
            !is.null(base$shareNames),
            !is.null(base$compute_z))
  model <- match.arg(model)
  index <- match.arg(index)
  
  # dados
  dat <- if (!is.null(df)) df else if (exists("df_iv", inherits = TRUE)) get("df_iv", inherits = TRUE) else stop("Forneça df= ou defina df_iv.")
  
  # garante ln_<preço> e z/z2
  for (p in base$priceNames) {
    ln <- paste0("ln_", p)
    if (!ln %in% names(dat) && p %in% names(dat)) dat[[ln]] <- log(pmax(as.numeric(dat[[p]]), 1e-12))
  }
  dat <- base$compute_z(dat)
  
  # RHS
  rhs_prices <- paste(paste0("ln_", base$priceNames), collapse = " + ")
  rhs <- paste(rhs_prices, "+ z")
  if (model == "QUAIDS" && isTRUE(use_z2)) rhs <- paste(rhs, "+ z2")
  
  # sistema k-1
  k <- length(base$shareNames)
  eqs <- vector("list", k - 1)
  names(eqs) <- paste0("w", 1:(k - 1))
  for (i in 1:(k - 1)) {
    lhs <- base$shareNames[i]
    eqs[[i]] <- stats::as.formula(paste(lhs, "~", rhs))
  }
  
  fit <- try(systemfit::systemfit(eqs, method = "SUR", data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(fit)
  
  # metadados
  attr(fit, "model")      <- toupper(model)
  attr(fit, "index")      <- toupper(index)
  attr(fit, "shareNames") <- base$shareNames
  attr(fit, "priceNames") <- base$priceNames
  attr(fit, "wbar")       <- base$wbar
  fit
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

iv_first_stage_diagnostics <- function(dat, priceNames, inst_map) {
  stopifnot(is.list(inst_map), length(inst_map) > 0)
  out <- list()
  for (tag in names(inst_map)) {
    Z <- inst_map[[tag]]
    if (!length(Z)) next
    res_tag <- list()
    for (p in priceNames) {
      y <- paste0("ln_", p)
      ok <- stats::complete.cases(dat[, c(y, Z), drop = FALSE])
      if (!any(ok)) { res_tag[[y]] <- list(R2 = NA_real_, adjR2 = NA_real_, n_inst = length(Z)); next }
      fit <- stats::lm(stats::reformulate(Z, response = y), data = dat[ok, , drop = FALSE])
      sm  <- summary(fit)
      res_tag[[y]] <- list(R2 = unname(sm$r.squared),
                           adjR2 = unname(sm$adj.r.squared),
                           n_inst = length(Z))
    }
    out[[tag]] <- res_tag
  }
  out
}

.softmax_safe <- function(eta, tau = 1, eps = 1e-6) {
  if (!all(is.finite(eta))) return(rep(NA_real_, length(eta)))
  z  <- (eta - max(eta)) / max(tau, 1e-9)
  ez <- exp(z)
  w  <- ez / sum(ez)
  w  <- pmax(w, eps)
  w / sum(w)
}
wrap_normalize_shares <- function(w, method = c("softmax","clip"), tau = 1, eps = 1e-6) {
  method <- match.arg(method)
  if (any(!is.finite(w))) return(rep(NA_real_, length(w)))
  if (method == "softmax") .softmax_safe(w, tau = tau, eps = eps) else { w <- pmax(w, eps); w / sum(w) }
}

predict_from_systemfit <- function(fit, ln_p, z, z2 = NA_real_,
                                   model = c("AIDS","QUAIDS"),
                                   shareNames, priceNames) {
  model <- match.arg(toupper(model), c("AIDS","QUAIDS"))
  cf <- try(coef(fit), silent = TRUE)
  if (inherits(cf, "try-error")) return(rep(NA_real_, length(shareNames)))
  k  <- length(shareNames); m <- length(priceNames)
  w_hat <- rep(NA_real_, k); names(w_hat) <- shareNames
  present <- grep("^w\\d+_", names(cf), value = TRUE)
  eq_ids  <- unique(sub("^(w\\d+)_.*$", "\\1", present))
  for (eq in eq_ids) {
    i <- as.integer(sub("^w", "", eq))
    if (!is.finite(i) || i < 1 || i > k) next
    alpha <- cf[paste0(eq, "_(Intercept)")]
    betaZ <- cf[paste0(eq, "_z")]
    betaZ2<- if (toupper(model) == "QUAIDS") cf[paste0(eq, "_z2")] else NA_real_
    alpha <- ifelse(is.na(alpha), 0, alpha)
    betaZ <- ifelse(is.na(betaZ), 0, betaZ)
    betaZ2<- ifelse(is.na(betaZ2),0, betaZ2)
    gam <- numeric(m); names(gam) <- priceNames
    for (j in seq_len(m)) {
      nm <- paste0(eq, "_ln_", priceNames[j])
      if (nm %in% names(cf)) gam[j] <- cf[nm]
    }
    w_hat[i] <- alpha + sum(gam * ln_p[priceNames]) + betaZ * z +
      if (toupper(model) == "QUAIDS") betaZ2 * (z2 %||% 0) else 0
  }
  miss <- which(!is.finite(w_hat))
  if (length(miss) == 1L) w_hat[miss] <- 1 - sum(w_hat[is.finite(w_hat)])
  w_hat
}

# Preditor genérico (usa preditores oficiais se existirem; senão, fallback via systemfit)
predict_shares_generic <- function(model_fit, x, p, tau = 1, normalize = "softmax", eps = 1e-6) {
  w_raw <- try({
    mdl <- attr(model_fit, "model")
    if (!is.null(mdl) && toupper(mdl) == "QUAIDS") {
      predict_shares_quaids(model_fit, x = x, p = p, normalize = NULL)
    } else if (!is.null(mdl) && toupper(mdl) == "AIDS") {
      predict_shares_aids(model_fit, x = x, p = p, normalize = NULL)
    } else {
      out <- try(predict_shares_quaids(model_fit, x = x, p = p, normalize = NULL), silent = TRUE)
      if (inherits(out, "try-error")) {
        out <- try(predict_shares_aids(model_fit, x = x, p = p, normalize = NULL), silent = TRUE)
      }
      out
    }
  }, silent = TRUE)
  if (inherits(w_raw, "try-error") || any(!is.finite(w_raw))) {
    mdl  <- (attr(model_fit, "model") %||% "AIDS")
    shN  <- attr(model_fit, "shareNames")  %||% get0("shareNames",  ifnotfound = NULL, inherits = TRUE)
    prN  <- attr(model_fit, "priceNames")  %||% get0("priceNames",  ifnotfound = NULL, inherits = TRUE)
    wbar <- attr(model_fit, "wbar")        %||% get0("wbar",        ifnotfound = NULL, inherits = TRUE)
    if (is.null(prN) || is.null(shN)) return(rep(NA_real_, length(p)))
    if (is.null(wbar)) wbar <- rep(1/length(prN), length(prN))
    wbar <- as.numeric(wbar); wbar <- wbar / sum(wbar)
    ln_p <- stats::setNames(log(p), prN)
    z    <- log(x) - sum(wbar * ln_p[prN])
    z2   <- if (toupper(mdl) == "QUAIDS") z*z else NA_real_
    w_raw <- try(
      predict_from_systemfit(model_fit, ln_p = ln_p, z = z, z2 = z2,
                             model = mdl, shareNames = shN, priceNames = prN),
      silent = TRUE
    )
    if (inherits(w_raw, "try-error")) return(rep(NA_real_, length(p)))
  }
  wrap_normalize_shares(w_raw, method = normalize, tau = tau, eps = eps)
}

.rmse_mean <- function(W_obs, W_hat) {
  ok <- is.finite(W_obs) & is.finite(W_hat)
  if (!any(ok)) return(NA_real_)
  dif <- W_obs[ok] - W_hat[ok]
  sqrt(mean(dif^2))
}
.mae_mean <- function(W_obs, W_hat) {
  ok <- is.finite(W_obs) & is.finite(W_hat)
  if (!any(ok)) return(NA_real_)
  mean(abs(W_obs[ok] - W_hat[ok]))
}

kfold_cv_metrics_safe <- function(dat, shareNames, lnpriceNames,
                                  build_fit, predict_fun,
                                  folds = 5, seed = 123,
                                  tau = 1, normalize = "softmax", eps = 1e-6) {
  set.seed(seed)
  n <- nrow(dat)
  idx <- stats::complete.cases(dat[, c(shareNames, lnpriceNames, "z"), drop = FALSE])
  id <- which(idx)
  if (length(id) < folds) return(list(rmse = NA_real_, mae = NA_real_))
  k_id <- sample(rep_len(1:folds, length.out = length(id)))
  rmses <- maes <- numeric(folds)
  for (k in 1:folds) {
    test_id  <- id[k_id == k]
    train_id <- setdiff(id, test_id)
    fit_k <- try(build_fit(dat[train_id, , drop = TRUE]), silent = TRUE)
    if (inherits(fit_k, "try-error")) { rmses[k] <- NA_real_; maes[k] <- NA_real_; next }
    LN <- as.matrix(dat[test_id, lnpriceNames, drop = FALSE])
    P  <- exp(LN)
    z  <- dat$z[test_id]
    # Stone index com w̄ do treino
    Wbar <- colMeans(dat[train_id, shareNames, drop = FALSE], na.rm = TRUE)
    Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
    base_lnP <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
    X  <- exp(z + base_lnP)
    W_hat <- matrix(NA_real_, nrow = length(test_id), ncol = length(shareNames))
    for (i in seq_along(test_id)) {
      w_i <- try(predict_fun(fit_k, x = X[i], p = P[i, ], tau = tau, normalize = normalize, eps = eps), silent = TRUE)
      if (!inherits(w_i, "try-error") && all(is.finite(w_i))) W_hat[i, ] <- w_i
    }
    W_obs <- as.matrix(dat[test_id, shareNames, drop = FALSE])
    rmses[k] <- .rmse_mean(as.numeric(W_obs), as.numeric(W_hat))
    maes[k]  <- .mae_mean (as.numeric(W_obs), as.numeric(W_hat))
  }
  list(rmse = mean(rmses, na.rm = TRUE), mae = mean(maes, na.rm = TRUE))
}

kfold_cv_by_good <- function(dat, shareNames, lnpriceNames,
                             build_fit, predict_fun,
                             folds = 5, seed = 123,
                             tau = 1, normalize = "softmax", eps = 1e-6) {
  set.seed(seed)
  idx <- stats::complete.cases(dat[, c(shareNames, lnpriceNames, "z"), drop = FALSE])
  id  <- which(idx)
  if (length(id) < folds) {
    return(data.frame(share = shareNames, n = 0, rmse = NA_real_, mae = NA_real_))
  }
  k_id <- sample(rep_len(1:folds, length(out = length(id))))
}

# (fix do bug na linha anterior)
kfold_cv_by_good <- function(dat, shareNames, lnpriceNames,
                             build_fit, predict_fun,
                             folds = 5, seed = 123,
                             tau = 1, normalize = "softmax", eps = 1e-6) {
  set.seed(seed)
  idx <- stats::complete.cases(dat[, c(shareNames, lnpriceNames, "z"), drop = FALSE])
  id  <- which(idx)
  if (length(id) < folds) {
    return(data.frame(share = shareNames, n = 0, rmse = NA_real_, mae = NA_real_))
  }
  k_id <- sample(rep_len(1:folds, length.out = length(id)))
  SSE <- setNames(numeric(length(shareNames)), shareNames)
  SAE <- setNames(numeric(length(shareNames)), shareNames)
  N   <- setNames(integer(length(shareNames)), shareNames)
  for (k in 1:folds) {
    test_id  <- id[k_id == k]; train_id <- setdiff(id, test_id)
    fit_k <- try(build_fit(dat[train_id, , drop = FALSE]), silent = TRUE)
    if (inherits(fit_k, "try-error")) next
    LN  <- as.matrix(dat[test_id, lnpriceNames, drop = FALSE])
    P   <- exp(LN)
    z   <- dat$z[test_id]
    Wbar <- colMeans(dat[train_id, shareNames, drop = FALSE], na.rm = TRUE)
    Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
    base_lnP <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
    X   <- exp(z + base_lnP)
    W_hat <- matrix(NA_real_, nrow = length(test_id), ncol = length(shareNames),
                    dimnames = list(NULL, shareNames))
    for (i in seq_along(test_id)) {
      pr <- try(predict_fun(fit_k, x = X[i], p = P[i, ], tau = tau,
                            normalize = normalize, eps = eps), silent = TRUE)
      if (!inherits(pr, "try-error") && all(is.finite(pr))) W_hat[i, ] <- pr
    }
    W_obs <- as.matrix(dat[test_id, shareNames, drop = FALSE])
    for (j in seq_along(shareNames)) {
      ok <- is.finite(W_obs[, j]) & is.finite(W_hat[, j])
      if (!any(ok)) next
      dif <- W_obs[ok, j] - W_hat[ok, j]
      SSE[j] <- SSE[j] + sum(dif^2)
      SAE[j] <- SAE[j] + sum(abs(dif))
      N[j]   <- N[j]   + sum(ok)
    }
  }
  data.frame(
    share = names(N),
    n     = as.integer(N),
    rmse  = ifelse(N > 0, sqrt(SSE / N), NA_real_),
    mae   = ifelse(N > 0, SAE / N,       NA_real_),
    row.names = NULL
  )
}

cv_tau_grid <- function(dat, shareNames, lnpriceNames, build_fit, predict_fun,
                        tau_grid = seq(0.15, 0.40, by = 0.05),
                        folds = 5, seed = 123, normalize = "softmax", eps = 1e-6) {
  set.seed(seed)
  out <- lapply(tau_grid, function(tau) {
    met <- kfold_cv_metrics_safe(
      dat = dat, shareNames = shareNames, lnpriceNames = lnpriceNames,
      build_fit = build_fit, predict_fun = predict_fun,
      folds = folds, seed = seed, tau = tau, normalize = normalize, eps = eps
    )
    c(tau = tau, rmse = met$rmse, mae = met$mae)
  })
  tab <- as.data.frame(do.call(rbind, out))
  tab <- tab[order(tab$rmse), ]
  list(tau = tab$tau[1], tau_grid = tau_grid, scores = tab$rmse, table = tab)
}

preprocess_min <- function(dat, priceNames, shareNames) {
  for (p in priceNames) {
    ln <- paste0("ln_", p)
    if (ln %in% names(dat)) next
    if (p %in% names(dat)) dat[[ln]] <- log(pmax(as.numeric(dat[[p]]), 1e-12))
  }
  lnpriceNames <- paste0("ln_", priceNames)
  base_ln <- if ("lndespesa" %in% names(dat)) dat$lndespesa else dat$ln_income
  if (is.null(base_ln)) stop("Nem 'lndespesa' nem 'ln_income' encontrados para montar 'z'.")
  Wbar <- colMeans(dat[, shareNames, drop = FALSE], na.rm = TRUE)
  Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
  LN   <- as.matrix(dat[, lnpriceNames, drop = FALSE])
  lnP  <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
  dat$z  <- as.numeric(base_ln) - lnP
  dat$z2 <- dat$z^2
  list(
    data = dat,
    price_sc_names = paste0(priceNames, "_sc"),
    ln_sd = stats::setNames(apply(LN, 2, stats::sd, na.rm = TRUE), priceNames)
  )
}

# ------------------------------------------------------------
# Monta base_rob com sua base
# ------------------------------------------------------------
base_rob <- make_base_rob(df_iv, priceNames, shareNames)

# ------------------------------------------------------------
# STEP 4 — Pré-processamento (antes, para ser usado no Step 2)
# ------------------------------------------------------------
step4_pp <- preprocess_min(df_iv, priceNames, shareNames)
dat_prep <- step4_pp$data
pipeline_results <- list()
pipeline_results$step4 <- list(preprocess = step4_pp, bases = NULL)
lnpriceNames <- paste0("ln_", priceNames)
cat("[Step4] dat_prep criado. Linhas:", nrow(dat_prep), "| Colunas:", ncol(dat_prep), "\n")

# ------------------------------------------------------------
# STEP 1 — Ajustes SUR: QUAIDS (z,z2) e AIDS (z)
# ------------------------------------------------------------
pipeline_results$step1 <- list(
  fit_quaids_full = fit_km1("QUAIDS", "Ls",
                            inst_try = inst_try_default,
                            base = base_rob, use_z2 = TRUE, df = df_iv),
  fit_aids_full   = fit_km1("AIDS",   "Ls",
                            inst_try = inst_try_default,
                            base = base_rob, use_z2 = FALSE, df = df_iv),
  wald_lambda     = NA_real_
)
wald_lambda_test_sur <- function(fit_sur_quaids) {
  cf <- try(stats::coef(fit_sur_quaids), silent = TRUE)
  V  <- try(stats::vcov(fit_sur_quaids), silent = TRUE)
  if (inherits(cf,"try-error") || inherits(V,"try-error")) return(NULL)
  idx <- grepl("_z2$", names(cf))
  if (sum(idx) == 0) return(NULL)
  b  <- cf[idx]; Vb <- V[idx, idx, drop = FALSE]
  W  <- as.numeric(t(b) %*% solve(Vb, b))
  df <- length(b); p  <- 1 - stats::pchisq(W, df)
  list(stat = W, df = df, pval = p)
}
if (!inherits(pipeline_results$step1$fit_quaids_full, "try-error") &&
    !is.null(pipeline_results$step1$fit_quaids_full)) {
  wl <- try(wald_lambda_test_sur(pipeline_results$step1$fit_quaids_full), silent = TRUE)
  pipeline_results$step1$wald_lambda <- if (inherits(wl,"try-error")) NULL else wl
} else {
  pipeline_results$step1$wald_lambda <- NULL
}

# ------------------------------------------------------------
# STEP 2 — Diagnósticos de 1º estágio (IV) + F incremental
# ------------------------------------------------------------
pipeline_results$step2 <- list(
  iv_diagnostics = try(iv_first_stage_diagnostics(dat_prep, priceNames, inst_map), silent = TRUE)
)
# F incremental (IVs | controles)
Z_cols <- grep("^(iv_|IV_)", names(dat_prep), value = TRUE)
C_cols <- intersect(c("z","z2","f_reg","f_cap"), names(dat_prep))
first_stage_F <- function(y, controls, inst, dat) {
  if (length(inst) == 0) return(data.frame(var=y, F=NA_real_, df1=NA_integer_, p=NA_real_))
  rhs_ctrl <- if (length(controls)>0) paste(controls, collapse=" + ") else "1"
  fml_full <- as.formula(paste(y, "~", paste(c(rhs_ctrl, paste(inst, collapse=" + ")), collapse=" + ")))
  fml_rest <- as.formula(paste(y, "~", rhs_ctrl))
  m_full <- lm(fml_full, data = dat); m_rest <- lm(fml_rest, data = dat)
  a <- tryCatch(anova(m_rest, m_full), error = function(e) NULL)
  if (is.null(a) || nrow(a) < 2 || is.na(a$F[2])) return(data.frame(var=y, F=NA_real_, df1=NA_integer_, p=NA_real_))
  data.frame(var=y, F=unname(a$F[2]), df1=unname(a$Df[2]), p=unname(a$`Pr(>F)`[2]))
}
endogs <- intersect(lnpriceNames, names(dat_prep))
if (length(endogs) == 0) endogs <- grep("^ln_preco_com_reforma[0-9]{2}$", names(dat_prep), value = TRUE)
iv_list <- lapply(endogs, function(y) first_stage_F(y, C_cols, Z_cols, dat_prep))
iv_strength <- do.call(rbind, iv_list)
if (is.null(iv_strength)) iv_strength <- data.frame(var=character(), F=numeric(), df1=integer(), p=numeric())
pipeline_results$step2$first_stage_F <- iv_strength
print(iv_strength)
cat("\n[IV] Regra prática: F ≳ 10 (incremental dos IVs dado os controles).\n")

# ------------------------------------------------------------
# STEP 3 — Escolha de τ (softmax) — grade grossa + micro-grid
# ------------------------------------------------------------
build_fit_cv_Ls <- function(train_df) {
  base_local <- make_base_rob(train_df, priceNames, shareNames)
  fit_km1("QUAIDS", "Ls",
          inst_try = inst_try_default,
          base = base_local, use_z2 = TRUE, df = train_df)
}
build_fit_cv_S <- function(train_df) {
  base_local <- make_base_rob(train_df, priceNames, shareNames)
  fit_km1("AIDS", "S",
          inst_try = inst_try_default,
          base = base_local, use_z2 = FALSE, df = train_df)
}

cv_coarse <- cv_tau_grid(
  dat = dat_prep, shareNames = shareNames, lnpriceNames = lnpriceNames,
  build_fit = build_fit_cv_Ls, predict_fun = predict_shares_generic,
  tau_grid = c(0.25, 0.5, 1, 2, 4, 8), folds = 5, seed = 123
)
tau0 <- cv_coarse$tau
tau_grid_fino <- seq(0.17, 0.23, by = 0.01)
compute_cv_Ls_at_tau <- function(dat, tau, K = 5, seed = 123) {
  kfold_cv_metrics_safe(
    dat = dat, shareNames = shareNames, lnpriceNames = lnpriceNames,
    build_fit = build_fit_cv_Ls, predict_fun = predict_shares_generic,
    folds = K, seed = seed, tau = tau, normalize = "softmax", eps = 1e-6
  )[c("rmse","mae")]
}
scores_finos <- lapply(tau_grid_fino, function(t) compute_cv_Ls_at_tau(dat_prep, tau = t))
rmse_finos   <- sapply(scores_finos, `[[`, "rmse")
mae_finos    <- sapply(scores_finos,  `[[`, "mae")
tab_fino <- data.frame(tau = tau_grid_fino, rmse = rmse_finos, mae = mae_finos)
tau_star <- tab_fino$tau[ which.min(tab_fino$rmse) ]
pipeline_results$step3 <- list(
  tau = tau_star,
  tau_grid = tau_grid_fino,
  scores = tab_fino$rmse,
  table = tab_fino[order(tab_fino$rmse), ],
  coarse = cv_coarse
)
cat(sprintf("[Step3] τ* = %.2f (RMSE=%.6f)\n", tau_star, min(tab_fino$rmse)))

# ------------------------------------------------------------
# STEP 5 — CV(Ls) agregada e por bem
# ------------------------------------------------------------
cv_Ls_metrics <- kfold_cv_metrics_safe(
  dat = dat_prep,
  shareNames = shareNames,
  lnpriceNames = lnpriceNames,
  build_fit = build_fit_cv_Ls,
  predict_fun = predict_shares_generic,
  folds = 5, seed = 123,
  tau = pipeline_results$step3$tau %||% 1,
  normalize = "softmax", eps = 1e-6
)
pipeline_results$step5 <- c(cv_Ls_metrics)

by_good_Ls <- kfold_cv_by_good(dat_prep, shareNames, lnpriceNames,
                               build_fit = build_fit_cv_Ls,
                               predict_fun = predict_shares_generic,
                               folds = 5, seed = 123, tau = pipeline_results$step3$tau)
pipeline_results$step5$by_good <- by_good_Ls

# RMSE/MAE ponderado por share médio (na amostra de CV)
share_means_cv <- sapply(shareNames, function(v) mean(dat_prep[[v]], na.rm = TRUE))
wg <- merge(by_good_Ls, data.frame(share = names(share_means_cv), weight = as.numeric(share_means_cv)), by="share", all.x=TRUE)
rmse_weighted <- sqrt( sum( (wg$rmse^2) * wg$weight, na.rm = TRUE ) / sum(wg$weight, na.rm = TRUE) )
mae_weighted  <- sum( wg$mae * wg$weight, na.rm = TRUE ) / sum(wg$weight, na.rm = TRUE)
pipeline_results$step5$rmse_weighted <- rmse_weighted
pipeline_results$step5$mae_weighted  <- mae_weighted
pipeline_results$step5$by_good_weighted <- wg[order(wg$rmse, decreasing = TRUE), ]

# ------------------------------------------------------------
# STEP 5c — Variante FE em w4 (comparativo rápido)
# ------------------------------------------------------------
build_fit_cv_Ls_FEw4 <- function(train_df) {
  base_local <- make_base_rob(train_df, priceNames, shareNames)
  dat <- base_local$compute_z(train_df)
  k <- length(shareNames)
  lnset <- paste0("ln_", priceNames)
  rhs_core <- paste(c(lnset, "z", "z2"), collapse = " + ")
  eqs <- vector("list", k - 1); names(eqs) <- paste0("w", 1:(k - 1))
  for (i in 1:(k - 1)) {
    lhs <- shareNames[i]
    rhs_i <- rhs_core
    if (identical(lhs, "w_despesahat4")) rhs_i <- paste(rhs_i, "+ f_reg + f_cap")
    eqs[[i]] <- stats::as.formula(paste(lhs, "~", rhs_i))
  }
  fit <- systemfit::systemfit(eqs, method = "SUR", data = dat)
  attr(fit, "model")      <- "QUAIDS"
  attr(fit, "index")      <- "LS"
  attr(fit, "shareNames") <- shareNames
  attr(fit, "priceNames") <- priceNames
  attr(fit, "wbar")       <- base_local$wbar
  fit
}
cv_Ls_FEw4 <- kfold_cv_metrics_safe(
  dat = dat_prep, shareNames = shareNames, lnpriceNames = lnpriceNames,
  build_fit = build_fit_cv_Ls_FEw4, predict_fun = predict_shares_generic,
  folds = 5, seed = 123, tau = pipeline_results$step3$tau, normalize = "softmax", eps = 1e-6
)
pipeline_results$step5$cv_FEw4 <- cv_Ls_FEw4$rmse
cat(sprintf("[Step5c] CV(Ls+FE@τ*): %.6f (vs Ls puro: %.6f)\n",
            pipeline_results$step5$cv_FEw4, pipeline_results$step5$rmse))

# ------------------------------------------------------------
# Holdout final (train/test split) — predição linha a linha
# ------------------------------------------------------------
set.seed(2025)
dat_all <- dat_prep
n <- nrow(dat_all)
test_frac <- 0.20
test_ids  <- sample.int(n, size = ceiling(test_frac * n))
train_ids <- setdiff(seq_len(n), test_ids)
dat_train <- dat_all[train_ids, ]
dat_test  <- dat_all[test_ids, ]
tau_holdout <- pipeline_results$step3$tau

sanitize_prices <- function(df, priceNames) {
  out <- df
  for (nm in priceNames) {
    ln_nm  <- paste0("ln_", nm)
    if (!ln_nm %in% names(out)) next
    bad <- !is.finite(out[[ln_nm]])
    if (!any(bad)) next
    alt <- paste0("ln_corr_", nm)
    if (alt %in% names(out)) {
      out[bad, ln_nm] <- out[bad, alt]
    } else if (nm %in% names(out)) {
      v <- suppressWarnings(log(pmax(as.numeric(out[[nm]]), 1e-12)))
      repl <- v
      repl[!is.finite(repl)] <- stats::median(out[[ln_nm]][is.finite(out[[ln_nm]])], na.rm = TRUE)
      out[bad, ln_nm] <- repl[bad]
    } else {
      out[bad, ln_nm] <- stats::median(out[[ln_nm]][is.finite(out[[ln_nm]])], na.rm = TRUE)
    }
  }
  out
}

fit_predict_quaids_Ls <- function(train, test, tau) {
  # treina
  fit <- build_fit_cv_Ls(train)
  # saneia logs/preços e prepara objetos
  test2 <- sanitize_prices(test, priceNames)
  n     <- nrow(test2); k <- length(shareNames)
  LN <- as.matrix(test2[, paste0("ln_", priceNames), drop = FALSE])
  P  <- exp(LN)
  z  <- test2$z
  # Stone index com w̄ do TREINO (consistente com a CV)
  Wbar <- colMeans(train[, shareNames, drop = FALSE], na.rm = TRUE)
  Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
  base_lnP <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
  X  <- exp(z + base_lnP)
  # predições linha a linha
  M  <- matrix(NA_real_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    pr <- try(predict_shares_generic(fit, x = X[i], p = P[i, ], tau = tau, normalize = "softmax", eps = 1e-6), silent = TRUE)
    if (!inherits(pr, "try-error") && all(is.finite(pr))) M[i, ] <- pr
  }
  # fallback: preenche NA com w̄ do treino e renormaliza
  for (j in seq_len(k)) {
    idx <- which(!is.finite(M[, j]))
    if (length(idx)) M[idx, j] <- Wbar[j]
  }
  rs <- rowSums(M)
  bad <- which(!is.finite(rs) | abs(rs - 1) > 1e-6)
  if (length(bad)) M[bad, ] <- sweep(M[bad, , drop = FALSE], 1, rowSums(M[bad, , drop = FALSE]), "/")
  colnames(M) <- shareNames
  M
}

# métricas robustas
rmse_mae <- function(y_true, y_pred) {
  stopifnot(all(colnames(y_true) == colnames(y_pred)))
  yt <- as.matrix(y_true); yp <- as.matrix(y_pred)
  ok <- rowSums(!is.finite(yt)) == 0 & rowSums(!is.finite(yp)) == 0
  yt <- yt[ok, , drop = FALSE]; yp <- yp[ok, , drop = FALSE]
  err  <- yp - yt
  rmse <- sqrt(colMeans(err^2, na.rm = TRUE))
  mae  <- colMeans(abs(err), na.rm = TRUE)
  list(
    agg = c(rmse = sqrt(mean(err^2, na.rm = TRUE)), mae = mean(abs(err), na.rm = TRUE)),
    by_good = data.frame(share = colnames(yt), rmse = rmse, mae = mae, row.names = NULL)
  )
}

# 1) Prever
preds_test <- fit_predict_quaids_Ls(dat_train, dat_test, tau = tau_holdout)

# 2) Verdade-terreno
share_cols <- grep("^w_despesahat[1-9][0-9]*$", names(dat_test), value = TRUE)
y_true <- dat_test[, share_cols, drop = FALSE]
# garante o mesmo ordenamento de colunas
preds_test <- preds_test[, shareNames, drop = FALSE]
colnames(preds_test) <- shareNames
colnames(y_true)     <- shareNames

# 3) Métricas
mm <- rmse_mae(y_true, as.data.frame(preds_test))

# 4) Ponderação por w* do interior (aqui uso o w̄ do treino p/ consistência)
w_star <- colMeans(dat_train[, shareNames, drop = FALSE], na.rm = TRUE)
w_star <- pmax(w_star, 1e-12); w_star <- w_star / sum(w_star)
by_good_weighted <- transform(
  mm$by_good,
  weight = as.numeric(w_star[match(share, names(w_star))])
)
rmse_weighted_holdout <- sqrt(sum(by_good_weighted$weight * (by_good_weighted$rmse^2), na.rm = TRUE))
mae_weighted_holdout  <- sum(by_good_weighted$weight * by_good_weighted$mae, na.rm = TRUE)

pipeline_results$holdout <- list(
  tau             = tau_holdout,
  agg             = mm$agg,
  by_good         = mm$by_good,
  rmse_weighted   = rmse_weighted_holdout,
  mae_weighted    = mae_weighted_holdout,
  by_good_weighted = by_good_weighted,
  n_train         = nrow(dat_train),
  n_test          = nrow(dat_test)
)

cat(sprintf("\n[Holdout @ τ=%.2f] RMSE: %.6f | MAE: %.6f | RMSE_w: %.6f | MAE_w: %.6f\n",
            pipeline_results$holdout$tau,
            pipeline_results$holdout$agg["rmse"],
            pipeline_results$holdout$agg["mae"],
            pipeline_results$holdout$rmse_weighted,
            pipeline_results$holdout$mae_weighted))

# ------------------------------------------------------------
# STEP 6 — SUR “puro” (AIDS-Style: só ln preços; 5 eqns)
# ------------------------------------------------------------
price_ln_set <- intersect(grep("^ln_", names(dat_prep), value = TRUE), paste0("ln_", priceNames))
if (length(price_ln_set) == 0) price_ln_set <- lnpriceNames
eqs <- list()
for (j in 2:length(shareNames)) {
  lhs <- shareNames[j]
  rhs <- paste(price_ln_set, collapse = " + ")
  eqs[[paste0("w", j)]] <- stats::as.formula(paste(lhs, "~", rhs))
}
fit_sur_try <- try(systemfit::systemfit(eqs, method = "SUR", data = dat_prep), silent = TRUE)
pipeline_results$step6 <- list(fit_sur = fit_sur_try)

# ------------------------------------------------------------
# STEP 7 — Ponto interior + elasticidades (SUR → Γ → Slutsky NSD)
# ------------------------------------------------------------
extract_gamma_from_sur <- function(fit_sur, shareNames, priceNames) {
  cf <- try(coef(fit_sur), silent = TRUE)
  if (inherits(cf, "try-error")) return(NULL)
  k <- length(shareNames); m <- length(priceNames)
  G <- matrix(0, nrow = k, ncol = m, dimnames = list(shareNames, priceNames))
  for (i in 2:k) {
    eq_tag <- paste0("w", i, "_")
    for (j in seq_along(priceNames)) {
      nm <- paste0(eq_tag, "ln_", priceNames[j])
      if (nm %in% names(cf)) G[i, j] <- as.numeric(cf[nm])
    }
  }
  if (k >= 2) G[1, ] <- -colSums(G[2:k, , drop = FALSE])
  G
}
project_gamma <- function(G, symmetric = TRUE, enforce_add_up = TRUE, iters = 8) {
  if (is.null(G)) return(NULL)
  H <- G; k <- nrow(H); m <- ncol(H)
  fix_add_up <- function(M) {
    M[1, ] <- -colSums(M[2:k, , drop = FALSE])
    M[, m] <- -rowSums(M[, 1:(m - 1), drop = FALSE])
    M
  }
  if (enforce_add_up) for (t in seq_len(iters)) H <- fix_add_up(H)
  if (symmetric) {
    H <- 0.5 * (H + t(H))
    if (enforce_add_up) for (t in seq_len(iters)) H <- fix_add_up(H)
  }
  H
}
elas_from_G_w <- function(G, w) {
  stopifnot(length(w) == nrow(G))
  w <- pmax(w, 1e-12); w <- w / sum(w)
  eta <- rep(1, length(w))
  k <- length(w)
  E_H <- matrix(0, nrow = k, ncol = k, dimnames = list(rownames(G), colnames(G)))
  for (i in 1:k) for (j in 1:k) E_H[i, j] <- -as.numeric(i == j) + G[i, j] / w[i]
  E_M <- E_H - outer(eta, w)
  adding_hicks <- as.numeric(E_H %*% w)
  symmetry_err <- max(abs(E_H - t(E_H)), na.rm = TRUE)
  W12 <- diag(sqrt(w), nrow = length(w))
  S   <- W12 %*% E_H %*% W12
  ev  <- try(eigen((S + t(S))/2, only.values = TRUE)$values, silent = TRUE)
  slutsky_max_eig <- if (inherits(ev, "try-error")) NA_real_ else max(Re(ev), na.rm = TRUE)
  list(expenditure = eta, hicks = E_H, marshall = E_M,
       checks = list(adding_hicks = adding_hicks, symmetry_err = symmetry_err, slutsky_max_eig = slutsky_max_eig))
}
project_S_slutsky <- function(S, nsdef = TRUE, tol = 1e-12) {
  S <- 0.5 * (S + t(S))
  k <- nrow(S)
  J <- diag(k) - matrix(1 / k, k, k)
  S <- J %*% S %*% J
  S[abs(S) < tol] <- 0
  if (nsdef) {
    ee  <- eigen(S, symmetric = TRUE)
    lam <- pmin(ee$values, 0)
    S   <- ee$vectors %*% (diag(lam, length(lam)) %*% t(ee$vectors))
    S <- 0.5 * (S + t(S)); S <- J %*% S %*% J
  }
  S
}

make_interior_point <- function(dat, lnpriceNames, shareNames) {
  LNbar <- colMeans(dat[, lnpriceNames, drop = FALSE], na.rm = TRUE)
  p  <- exp(LNbar)
  Wbar <- colMeans(dat[, shareNames,   drop = FALSE], na.rm = TRUE)
  Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
  zbar <- mean(dat$z, na.rm = TRUE)
  x  <- exp(zbar + sum(Wbar * LNbar))
  list(p = stats::setNames(p, gsub("^ln_", "", names(LNbar))), x = x, z = zbar, wbar = Wbar, lnPbar = LNbar)
}

project_and_elas <- (function() {
  pt_int <- make_interior_point(dat_prep, lnpriceNames, shareNames)
  G_sur_raw <- if (!inherits(pipeline_results$step6$fit_sur, "try-error")) {
    extract_gamma_from_sur(pipeline_results$step6$fit_sur, shareNames, priceNames)
  } else NULL
  if (is.null(G_sur_raw)) {
    return(list(interior = list(at = list(p = pt_int$p, x = pt_int$x, z = pt_int$z, w = pt_int$wbar)), micro = list(mid = NULL)))
  }
  G_sur <- project_gamma(G_sur_raw, symmetric = TRUE, enforce_add_up = TRUE, iters = 8)
  w  <- pt_int$wbar; names(w) <- shareNames
  k  <- length(w); D  <- diag(w, k)
  el0 <- elas_from_G_w(G_sur, w)
  S_in  <- D %*% el0$hicks
  S_out <- project_S_slutsky(S_in, nsdef = TRUE, tol = 1e-12)
  Hh    <- solve(D, S_out)
  M     <- Hh - outer(el0$expenditure, w)
  add_hicks <- rowSums(Hh)
  sym_err   <- max(abs(D %*% Hh - t(D %*% Hh)))
  Ssym      <- 0.5 * (D %*% Hh + t(D %*% Hh))
  max_eig_S <- max(eigen(Ssym, only.values = TRUE)$values)
  list(
    interior = list(
      at = list(p = pt_int$p, x = pt_int$x, z = pt_int$z, w = w),
      expenditure = el0$expenditure,
      marshall    = `dim<-`(M, dim(M)),
      hicks       = `dim<-`(Hh, dim(Hh)),
      checks      = list(adding_hicks = add_hicks, symmetry_err = sym_err, slutsky_max_eig = max_eig_S),
      S = `dim<-`(D %*% Hh, dim(D %*% Hh))
    ),
    micro = list(mid = NULL)
  )
})()
if (!is.null(project_and_elas$interior$marshall)) {
  dimnames(project_and_elas$interior$marshall) <- list(shareNames, priceNames)
  dimnames(project_and_elas$interior$hicks)    <- list(shareNames, priceNames)
  dimnames(project_and_elas$interior$S)        <- list(shareNames, shareNames)
}
pipeline_results$step7 <- project_and_elas

# ------------------------------------------------------------
# STEP 8 — AIDS (Stone) + CV(S)
# ------------------------------------------------------------
fit_aids_S <- try(fit_km1("AIDS", "S",
                          inst_try = inst_try_default,
                          base = base_rob, use_z2 = FALSE, df = df_iv), silent = TRUE)
cv_S_metrics <- try(kfold_cv_metrics_safe(
  dat = dat_prep,
  shareNames = shareNames, lnpriceNames = lnpriceNames,
  build_fit = build_fit_cv_S,
  predict_fun = predict_shares_generic,
  folds = 5, seed = 123,
  tau = pipeline_results$step3$tau %||% 1,
  normalize = "softmax", eps = 1e-6
), silent = TRUE)
pipeline_results$step8 <- list(
  fit_aids_S   = if (!inherits(fit_aids_S, "try-error")) fit_aids_S else NULL,
  cv_Ls_at_tau = pipeline_results$step5$rmse,
  cv_S_at_tau  = if (!inherits(cv_S_metrics, "try-error")) cv_S_metrics$rmse else NA_real_,
  by_good      = try(kfold_cv_by_good(
    dat = dat_prep, shareNames = shareNames, lnpriceNames = lnpriceNames,
    build_fit = build_fit_cv_S, predict_fun = predict_shares_generic,
    folds = 5, seed = 123, tau = pipeline_results$step3$tau
  ), silent = TRUE)
)
rmse_Ls <- pipeline_results$step5$rmse
rmse_S  <- pipeline_results$step8$cv_S_at_tau
pipeline_results$step8$winner <- if (is.finite(rmse_Ls) && is.finite(rmse_S)) {
  if (rmse_Ls <= rmse_S) "QUAIDS(Ls)" else "AIDS(S)"
} else NA_character_

# -----------------------------------------------------------------------
# STEP 9 — Hygiene
# -----------------------------------------------------------------------
share_zero_frac <- vapply(dat_prep[ , shareNames, drop=FALSE], function(v) mean(v == 0, na.rm=TRUE), 0)
share_one_frac  <- vapply(dat_prep[ , shareNames, drop=FALSE], function(v) mean(v == 1, na.rm=TRUE), 0)
ln_sd <- stats::setNames(vapply(lnpriceNames, function(nm) stats::sd(dat_prep[[nm]], na.rm=TRUE), 0), priceNames)

pipeline_results$step9 <- list(
  hygiene = list(
    ln_price_sd = ln_sd,
    share_zero_frac = share_zero_frac,
    share_one_frac  = share_one_frac
  )
)

# >>> COLE APÓS O Step 8, ANTES dos "prints rápidos", se quiser garantir <<<

if (is.null(pipeline_results$step8$winner)) {
  rmse_Ls <- pipeline_results$step8$cv_Ls_at_tau
  rmse_S  <- pipeline_results$step8$cv_S_at_tau
  pipeline_results$step8$winner <- ifelse(is.finite(rmse_Ls) && rmse_Ls <= rmse_S,
                                          "QUAIDS(Ls)", "AIDS(Stone)")
}

#Calibração--------------------------------------------------------------
# =========================
# Integração oficial da calibração no pipeline
# =========================

# -- util: projeção no simplex (não-neg. e soma=1)
# garante que %||% existe
`%||%` <- get0("%||%", ifnotfound = function(a,b) if (!is.null(a)) a else b, inherits = TRUE)

project_rows_simplex <- get0("project_rows_simplex", ifnotfound = NULL, inherits = TRUE) %||% function(M) {
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M[M < 0] <- 0
  rs <- rowSums(M)
  zero <- !is.finite(rs) | rs <= 0
  if (any(zero)) M[zero, ] <- 1 / ncol(M)
  sweep(M, 1, rowSums(M), "/")
}

# -- OOF preds no treino (reusa seu builder/predictor)
## garante o operador %||%
`%||%` <- get0("%||%", ifnotfound = function(a, b) if (!is.null(a)) a else b, inherits = TRUE)

## projeção simples no simplex (clip+renorm), se você ainda não tiver uma
if (!exists("project_rows_simplex", inherits = TRUE)) {
  project_rows_simplex <- function(M) {
    M <- as.matrix(M)
    M[!is.finite(M)] <- 0
    M[M < 0] <- 0
    rs <- rowSums(M)
    bad <- !is.finite(rs) | rs <= 0
    if (any(bad)) M[bad, ] <- 1 / ncol(M)
    sweep(M, 1, rowSums(M), "/")
  }
}

## OOF para QUAIDS(Ls) com refit por fold
get_oof_preds_Ls <- get0("get_oof_preds_Ls", ifnotfound = NULL, inherits = TRUE) %||% 
  function(dat, K = 5, seed = 123, tau = NULL) {
    stopifnot(exists("shareNames", inherits = TRUE),
              exists("lnpriceNames", inherits = TRUE),
              exists("build_fit_cv_Ls", inherits = TRUE),
              exists("predict_shares_generic", inherits = TRUE))
    
    tau <- tau %||% (get0("tau_holdout", inherits = TRUE) %||% get0("pipeline_results", inherits = TRUE)$step3$tau %||% 1)
    set.seed(seed)
    
    req <- c(shareNames, lnpriceNames, "z")
    miss <- setdiff(req, names(dat))
    if (length(miss)) stop("Colunas ausentes em 'dat': ", paste(miss, collapse = ", "))
    
    n <- nrow(dat); k <- length(shareNames)
    oof <- matrix(NA_real_, n, k, dimnames = list(NULL, shareNames))
    
    ok  <- stats::complete.cases(dat[, req, drop = FALSE])
    ids <- which(ok)
    if (!length(ids)) stop("Nenhuma linha completa para OOF.")
    folds <- sample(rep_len(1:K, length.out = length(ids)))
    
    # fallback para preencher buracos (w* do ponto interior, se existir)
    w_star <- get0("pipeline_results", inherits = TRUE)$step7$interior$at$w
    if (is.null(w_star)) w_star <- rep(1 / k, k)
    w_star <- as.numeric(w_star[shareNames])  # garante ordem
    
    for (f in 1:K) {
      test_id  <- ids[folds == f]
      train_id <- setdiff(ids, test_id)
      if (!length(test_id) || !length(train_id)) next
      
      fit_k <- try(build_fit_cv_Ls(dat[train_id, , drop = FALSE]), silent = TRUE)
      if (inherits(fit_k, "try-error")) next
      
      LN <- as.matrix(dat[test_id, lnpriceNames, drop = FALSE])
      mode(LN) <- "numeric"             # força numeric (labeled/factor → numeric)
      P  <- exp(LN)
      
      Wbar <- colMeans(dat[train_id, shareNames, drop = FALSE], na.rm = TRUE)
      Wbar <- pmax(Wbar, 1e-12); Wbar <- Wbar / sum(Wbar)
      base_lnP <- rowSums(LN * matrix(as.numeric(Wbar), nrow(LN), length(Wbar), byrow = TRUE))
      X <- exp(dat$z[test_id] + base_lnP)
      
      for (i in seq_along(test_id)) {
        wi <- try(
          predict_shares_generic(fit_k, x = X[i], p = P[i, ], tau = tau,
                                 normalize = "softmax", eps = 1e-6),
          silent = TRUE
        )
        if (inherits(wi, "try-error")) next
        wi <- suppressWarnings(as.numeric(wi))
        if (!length(wi)) next
        if (length(wi) == 1L) wi <- rep(wi, k)
        if (length(wi) != k) next
        if (any(!is.finite(wi))) {
          wi[!is.finite(wi)] <- w_star[!is.finite(wi)]
        }
        oof[test_id[i], ] <- wi
      }
    }
    
    # preenche linhas faltantes com w*
    badrow <- which(!stats::complete.cases(oof))
    if (length(badrow)) oof[badrow, ] <- matrix(rep(w_star, length(badrow)), ncol = k, byrow = TRUE)
    
    # projeta no simplex (não-negativo e soma 1 por linha)
    oof <- project_rows_simplex(oof)
    oof
  }

# -- métricas auxiliares
.metrics_simple <- function(pred, yt) {
  rmse <- sqrt(mean((as.matrix(pred) - as.matrix(yt))^2, na.rm=TRUE))
  mae  <- mean(abs(as.matrix(pred) - as.matrix(yt)), na.rm=TRUE)
  c(rmse=rmse, mae=mae)
}
.metrics_full <- function(y_true, y_pred, w_star=NULL) {
  stopifnot(all(colnames(y_true) == colnames(y_pred)))
  yt <- as.matrix(y_true); yp <- as.matrix(y_pred)
  err  <- yp - yt
  agg  <- c(rmse = sqrt(mean(err^2, na.rm = TRUE)),
            mae  = mean(abs(err), na.rm = TRUE))
  by   <- data.frame(
    share = colnames(yt),
    rmse  = sqrt(colMeans(err^2, na.rm=TRUE)),
    mae   = colMeans(abs(err), na.rm=TRUE),
    row.names = NULL
  )
  out <- list(agg=agg, by_good=by)
  if (!is.null(w_star)) {
    ww <- as.numeric(w_star[colnames(yt)])
    out$rmse_weighted <- sqrt(sum(ww * by$rmse^2, na.rm=TRUE))
    out$mae_weighted  <- sum(ww * by$mae,     na.rm=TRUE)
    out$by_good_weighted <- transform(by, weight=ww)
  }
  out
}

# -- treina calibrador no TRAIN (OOF) e escolhe método
fit_share_calibrator <- function(dat_train, K=5, seed=123, tau=tau_holdout, pick=c("auto","shift","affine")) {
  pick <- match.arg(pick)
  oof <- get_oof_preds_Ls(dat_train, K=K, seed=seed, tau=tau)
  y   <- as.matrix(dat_train[, shareNames, drop=FALSE])
  
  # a) SHIFT (remove viés médio)
  cal_shift <- colMeans(y - oof, na.rm=TRUE)
  apply_shift <- function(P) project_rows_simplex(sweep(as.matrix(P), 2, cal_shift, "+"))
  met_shift   <- .metrics_simple(apply_shift(oof), y)
  
  # b) AFIM (y ≈ a + b·pred), por bem
  fit_affine <- function(y, x) {
    do.call(rbind, lapply(seq_along(shareNames), function(j){
      ok <- is.finite(y[,j]) & is.finite(x[,j])
      if (!any(ok)) return(c(a=0,b=1))
      cf <- coef(lm(y[ok,j] ~ x[ok,j])); c(a=unname(cf[1]), b=unname(cf[2]))
    })) |> `rownames<-`(shareNames)
  }
  cal_aff <- fit_affine(y, oof)
  apply_affine <- function(P){
    M <- as.matrix(P)
    A <- matrix(cal_aff[,"a"], nrow=nrow(M), ncol=ncol(M), byrow=TRUE)
    B <- matrix(cal_aff[,"b"], nrow=nrow(M), ncol=ncol(M), byrow=TRUE)
    project_rows_simplex(A + B*M)
  }
  met_aff <- .metrics_simple(apply_affine(oof), y)
  
  # seleção
  choice <- switch(pick,
                   auto   = if (met_aff["rmse"] <= met_shift["rmse"]) "affine" else "shift",
                   shift  = "shift",
                   affine = "affine"
  )
  
  list(
    choice = choice,
    shift  = list(delta = cal_shift, oof_metrics = met_shift),
    affine = list(a_b = cal_aff,     oof_metrics = met_aff),
    apply  = if (choice=="affine") apply_affine else apply_shift
  )
}

# ============== roda e grava no pipeline ==============
calib <- fit_share_calibrator(dat_train, K=5, seed=123, tau=tau_holdout, pick="auto")

preds_test_cal <- calib$apply(preds_test)

# métricas completas (inclui ponderadas por w*)
y_true_test <- dat_test[, shareNames, drop=FALSE]
w_star      <- pipeline_results$step7$interior$at$w
mm_base     <- .metrics_full(y_true_test, preds_test,     w_star)
mm_cal      <- .metrics_full(y_true_test, preds_test_cal, w_star)

# salva
pipeline_results$holdout$calibration <- list(
  method   = calib$choice,
  params   = if (calib$choice=="affine") calib$affine$a_b else calib$shift$delta,
  oof      = list(affine = calib$affine$oof_metrics,
                  shift  = calib$shift$oof_metrics),
  metrics  = list(before = mm_base, after = mm_cal),
  bias_by_good = list(
    before = colMeans(y_true_test - preds_test,     na.rm=TRUE),
    after  = colMeans(y_true_test - preds_test_cal, na.rm=TRUE)
  )
)

# print rápido
cat("\n[Calib] método escolhido:", pipeline_results$holdout$calibration$method, "\n")
cat(sprintf("[Holdout] RMSE base: %.6f  →  cal: %.6f | MAE base: %.6f  →  cal: %.6f\n",
            mm_base$agg["rmse"], mm_cal$agg["rmse"],
            mm_base$agg["mae"],  mm_cal$agg["mae"]))
cat(sprintf("[Holdout] RMSE_w base: %.6f  →  cal: %.6f | MAE_w base: %.6f  →  cal: %.6f\n",
            mm_base$rmse_weighted, mm_cal$rmse_weighted,
            mm_base$mae_weighted,  mm_cal$mae_weighted))

# (opcional) guarda as predições calibradas
pipeline_results$holdout$preds_calibrated <- preds_test_cal


# -----------------------------------------------------------------------
# Prints rápidos
# -----------------------------------------------------------------------
cat("\n[Step1] QUAIDS ok? ", !is.null(pipeline_results$step1$fit_quaids_full),
    " | AIDS ok? ", !is.null(pipeline_results$step1$fit_aids_full), "\n", sep = "")
if (is.list(pipeline_results$step1$wald_lambda)) {
  wl <- pipeline_results$step1$wald_lambda
  cat(sprintf("        Wald(λ): W=%.3f (df=%d), p=%.4g\n", wl$stat, wl$df, wl$pval))
} else {
  cat("        Wald(λ): NA\n")
}
cat("[Step3] tau* = ", pipeline_results$step3$tau %||% NA, "\n", sep = "")
cat(sprintf("[Step5] CV(Ls) — RMSE: %.6f | MAE: %.6f\n", pipeline_results$step5$rmse, pipeline_results$step5$mae))
cat(sprintf("[Step8] CV(Ls@τ*): %.6f | CV(S@τ*): %.6f\n",
            pipeline_results$step8$cv_Ls_at_tau, pipeline_results$step8$cv_S_at_tau))
cat("[Step6] SUR status: ", if (inherits(pipeline_results$step6$fit_sur, "try-error")) "erro" else "ok", "\n", sep = "")
if (!is.null(pipeline_results$step7$interior$w)) {
  cat("[Step7] w(interior) = ", paste(round(pipeline_results$step7$interior$w, 4), collapse = " "), "\n", sep = "")
}

if (!is.null(pipeline_results$step5$by_good)) {
  bg <- pipeline_results$step5$by_good
  cat(sprintf("[Step5b] Melhor/PIOR (Ls) por RMSE: %s / %s\n",
              bg$share[which.min(bg$rmse)], bg$share[which.max(bg$rmse)]))
}
if (!is.null(pipeline_results$step8$by_good)) {
  bgS <- pipeline_results$step8$by_good
  cat(sprintf("[Step8b] Melhor/PIOR (S)  por RMSE: %s / %s\n",
              bgS$share[which.min(bgS$rmse)], bgS$share[which.max(bgS$rmse)]))
}
if (!is.null(pipeline_results$step8$winner)) {
  cat(sprintf("[Resumo] Vencedor (RMSE agregado): %s\n", pipeline_results$step8$winner))
}

cat("\n[Calib] método escolhido:", pipeline_results$holdout$calibration$method, "\n")
cat(sprintf("[Holdout] RMSE base: %.6f  →  cal: %.6f | MAE base: %.6f  →  cal: %.6f\n",
            mm_base$agg["rmse"], mm_cal$agg["rmse"],
            mm_base$agg["mae"],  mm_cal$agg["mae"]))
cat(sprintf("[Holdout] RMSE_w base: %.6f  →  cal: %.6f | MAE_w base: %.6f  →  cal: %.6f\n",
            mm_base$rmse_weighted, mm_cal$rmse_weighted,
            mm_base$mae_weighted,  mm_cal$mae_weighted))

invisible(pipeline_results)

# --- dump opcional (como você faz) ---
print(pipeline_results$step1)
print(pipeline_results$step2)
print(pipeline_results$step3)
print(pipeline_results$step4)
print(pipeline_results$step5)
print(pipeline_results$step6)
print(pipeline_results$step7)
print(pipeline_results$step8)
print(pipeline_results$step9)
print(pipeline_results$step3_refine)
cat(sprintf("\n[Holdout @ τ=%.2f] RMSE: %.6f | MAE: %.6f | RMSE_w: %.6f | MAE_w: %.6f\n",
            pipeline_results$holdout$tau,
            pipeline_results$holdout$agg["rmse"],
            pipeline_results$holdout$agg["mae"],
            pipeline_results$holdout$rmse_weighted,
            pipeline_results$holdout$mae_weighted))
print(pipeline_results$holdout$preds_calibrated)

# ============================================================
# QA do pipeline — Checagens 1 a 7
# ============================================================

# --------- Helpers básicos ----------
`%||%` <- function(x, y) if (is.null(x)) y else x

quiet_try <- function(expr) try(expr, silent = TRUE)

stop_if <- function(cond, msg) { if (isTRUE(cond)) stop(msg, call. = FALSE) }

pkg <- function(p){ if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet=TRUE); suppressPackageStartupMessages(library(p, character.only=TRUE)) }

fmt_num <- function(x) ifelse(is.na(x), "NA", formatC(x, digits=6, format="f"))

# Descobre nomes de colunas de shares "observadas" no dataset
detect_share_cols <- function(dat) {
  s <- grep("^w_", names(dat), value = TRUE)
  # mantém apenas as de despesa (evita colaterais tipo "wbar" etc.)
  grep("^w_despesa", s, value = TRUE)
}

# Pega dataset "base" do pipeline (pré-processado)
get_dat <- function() {
  pipeline_results$step4$preprocess$data %||%
    pipeline_results$step4$preprocess$dat %||%
    stop("Não achei o data.frame pré-processado em pipeline_results$step4$preprocess.")
}

# Pega matriz de predições calibradas do holdout
get_preds_holdout <- function() {
  pr <- pipeline_results$holdout
  stop_if(is.null(pr), "pipeline_results$holdout ausente.")
  stop_if(is.null(pr$preds_calibrated), "preds_calibrated ausente em pipeline_results$holdout.")
  as.matrix(pr$preds_calibrated)
}

# Tenta descobrir os índices (linhas do data) usados no holdout
get_holdout_idx <- function(n_rows_needed, dat_n) {
  pr <- pipeline_results$holdout
  cand <- NULL
  for (nm in c("ids","index","idx","test_idx","test_id","holdout_idx","rows")) {
    if (!is.null(pr[[nm]])) { cand <- pr[[nm]]; break }
  }
  if (!is.null(cand)) return(as.integer(cand))
  warning("Índices do holdout não encontrados no objeto. As análises por subgrupo e gráficos serão puladas (precisam ligar predições aos observados).")
  integer(0)
}

# Métricas simples
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

# ---- 1) Simplex: soma=1 e não-negatividade nas predições calibradas ----
check_simplex <- function(P, tol=1e-8) {
  row_sum_err <- rowSums(P) - 1
  neg_min     <- suppressWarnings(min(P, na.rm=TRUE))
  list(
    n = nrow(P),
    max_abs_sum_err = max(abs(row_sum_err), na.rm = TRUE),
    n_bad_sum = sum(abs(row_sum_err) > tol),
    min_entry = neg_min,
    n_neg = sum(P < -tol, na.rm=TRUE)
  )
}

# ---- 2) Erro por subgrupo (região, capital/interior, quintil) ----
by_group_errors <- function(dat, obs, pred, group_vars, focus_share = NULL) {
  stop_if(nrow(obs) != nrow(pred), "OBS e PRED precisam ter mesmo número de linhas.")
  share_names <- colnames(obs)
  if (!is.null(focus_share)) {
    stop_if(!(focus_share %in% share_names), sprintf("focus_share '%s' não encontrado.", focus_share))
    share_names <- focus_share
    obs <- obs[, focus_share, drop=FALSE]
    pred <- pred[, focus_share, drop=FALSE]
  }
  out <- list()
  for (g in group_vars) {
    stop_if(!(g %in% names(dat)), sprintf("Variável de grupo '%s' não existe no data.", g))
    G <- dat[[g]]
    tab <- do.call(rbind, lapply(share_names, function(sj){
      e <- obs[, sj] - pred[, sj]
      ag <- aggregate(e, by = list(G), FUN = function(z) c(n=length(z), rmse=rmse(z), mae=mae(z), bias=mean(z)))
      # "explode" o vetor retornado
      M <- do.call(data.frame, ag$x |> lapply(as.numeric))
      names(M) <- c("n","rmse","mae","bias")
      data.frame(group_var=g, group=ag[[1]], share=sj, M, row.names=NULL)
    }))
    rownames(tab) <- NULL
    out[[g]] <- tab[order(tab$share, tab$group), ]
  }
  out
}

# ---- 3) Estabilidade do tau* com re-splits ----
assess_tau_stability <- function(B = 10, seeds = NULL) {
  # Requer uma função do seu pipeline que refaça a CV p/ cada tau.
  # Tentamos algumas convenções comuns:
  fn <- NULL
  for (nm in c("select_tau_cv", "cv_quaids_Ls", "recompute_tau_cv", "select_tau")) {
    if (exists(nm, mode = "function")) { fn <- get(nm); break }
  }
  if (is.null(fn)) {
    message("Não encontrei função de seleção de τ via CV no ambiente. Pulando a checagem 3.")
    return(invisible(NULL))
  }
  tau_grid <- pipeline_results$step3$tau_grid %||% seq(0.15, 0.30, by=0.01)
  if (is.null(seeds)) seeds <- sample(1:1e6, B)
  winners <- numeric(B)
  scores  <- matrix(NA_real_, nrow=B, ncol=length(tau_grid), dimnames=list(NULL, paste0("tau_", tau_grid)))
  for (b in seq_len(B)) {
    set.seed(seeds[b])
    res_b <- quiet_try(fn(tau_grid = tau_grid, seed = seeds[b]))
    if (inherits(res_b, "try-error") || is.null(res_b$rmse_by_tau)) next
    # espera-se um named numeric com RMSE por tau
    sc   <- res_b$rmse_by_tau
    scores[b, match(names(sc), colnames(scores))] <- sc
    winners[b] <- as.numeric(names(sc)[which.min(sc)])
  }
  list(tau_grid = tau_grid, winners = winners, table = scores)
}

# ---- 4) Bootstrap de elasticidades (ponto interior) ----
bootstrap_elasticities <- function(B = 200, seed = 123) {
  # Precisa de funções do seu pipeline que (i) refitam o QUAIDS e (ii) calculem elasticidades em 1 ponto.
  for (nm in c("fit_quaids_full", "build_fit_quaids", "fit_quaids")) {
    if (exists(nm, mode="function")) { .fit <- get(nm); break }
  }
  for (nm in c("compute_interior", "elasticities_at", "get_elasticities")) {
    if (exists(nm, mode="function")) { .el  <- get(nm);  break }
  }
  if (!exists(".fit") || !exists(".el")) {
    message("Funções para re-estimar QUAIDS/elasticidades não encontradas. Pulando a checagem 4.")
    return(invisible(NULL))
  }
  dat <- get_dat()
  set.seed(seed)
  HICKS  <- list(); MARSH <- list()
  for (b in seq_len(B)) {
    id <- sample.int(nrow(dat), replace = TRUE)
    fit_b <- quiet_try(.fit(dat[id, , drop=FALSE]))
    if (inherits(fit_b, "try-error")) next
    el_b  <- quiet_try(.el(fit_b))
    if (inherits(el_b, "try-error")) next
    HICKS[[length(HICKS)+1]] <- el_b$hicks
    MARSH[[length(MARSH)+1]] <- el_b$marshall
  }
  # Empilha e tira quantis
  qfun <- function(L, probs=c(0.025, 0.5, 0.975)) {
    if (!length(L)) return(NULL)
    A <- abind::abind(L, along = 3)  # p x p x B
    qs <- apply(A, c(1,2), function(v) quantile(v, probs = probs, na.rm=TRUE))
    list(q = qs, mean = apply(A, c(1,2), mean, na.rm=TRUE), sd = apply(A, c(1,2), sd, na.rm=TRUE))
  }
  if (!requireNamespace("abind", quietly=TRUE)) install.packages("abind", quiet=TRUE)
  list(hicks = qfun(HICKS), marshall = qfun(MARSH))
}

# ---- 5) IV: KP e Hansen J numa eq. estrutural (ex.: w1) ----
iv_tests <- function(dep_share = NULL) {
  # Usa AER::ivreg numa equação (ex.: w_despesahat1 ~ preços + z + z2 | instrumentos + z + z2)
  if (!requireNamespace("AER", quietly=TRUE)) install.packages("AER", quiet=TRUE)
  if (!requireNamespace("sandwich", quietly=TRUE)) install.packages("sandwich", quiet=TRUE)
  if (!requireNamespace("lmtest", quietly=TRUE)) install.packages("lmtest", quiet=TRUE)
  pkg("AER"); pkg("sandwich"); pkg("lmtest")
  
  dat <- get_dat()
  share_cols <- detect_share_cols(dat)
  if (is.null(dep_share)) dep_share <- share_cols[1]
  stop_if(!(dep_share %in% share_cols), sprintf("Share '%s' não encontrado no data.", dep_share))
  
  price_cols <- grep("^ln_preco_com_reforma2[1-6]$", names(dat), value = TRUE)
  stop_if(length(price_cols) == 0, "Colunas ln_preco_com_reforma2# não encontradas.")
  # exógenas "lado direito"
  x_exog <- intersect(c("z", "z2"), names(dat))
  # instrumentos: tudo que começa com iv_/IV_*
  Z_inst <- grep("^(iv_|IV_)", names(dat), value = TRUE)
  stop_if(length(Z_inst) == 0, "Não encontrei colunas de instrumentos (iv_/IV_).")
  
  f_rhs  <- paste(c(price_cols, x_exog), collapse = " + ")
  f_iv   <- paste(c(Z_inst, x_exog),  collapse = " + ")
  f_full <- as.formula(sprintf("%s ~ %s | %s", dep_share, f_rhs, f_iv))
  
  fit_iv <- AER::ivreg(f_full, data = dat)
  # vcov HAC (carteira ampla)
  V <- sandwich::vcovHC(fit_iv, type = "HC1")
  diag_iv <- AER::linearHypothesis(fit_iv, diag(nrow(coef(summary(fit_iv)))), vcov.=V) # força cálculo
  
  s <- summary(fit_iv, vcov = V, diagnostics = TRUE)
  cat("\n[IV - equação estrutural]\n")
  print(s$coefficients)
  
  # Diagnósticos reportados pelo AER::ivreg:
  # (Cragg-Donald, Wu-Hausman, Sargan). Algumas versões mostram Kleibergen-Paap.
  if (!is.null(s$diagnostics)) {
    cat("\n[Diagnósticos IV]\n")
    print(s$diagnostics)
    cat("\nObs.: Se 'Kleibergen-Paap' não aparecer, a tabela mostra Cragg–Donald e Sargan.\n")
  } else {
    cat("\nDiagnósticos não disponíveis nesta versão do AER::ivreg.\n")
  }
  invisible(list(fit = fit_iv, summary = s, formula = f_full))
}

# ---- 6) Gráficos (parity plot + resíduos) para o holdout ----
make_plots <- function() {
  if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2", quiet=TRUE)
  pkg("ggplot2")
  
  dat <- get_dat()
  P   <- get_preds_holdout()
  idx <- get_holdout_idx(nrow(P), nrow(dat))
  if (!length(idx)) { message("Sem índices do holdout → pulando gráficos."); return(invisible(NULL)) }
  
  share_cols <- detect_share_cols(dat)
  stop_if(length(share_cols) != ncol(P), "N° de shares no data ≠ n° de colunas em preds_calibrated.")
  
  OBS <- as.matrix(dat[idx, share_cols, drop=FALSE])
  
  long <- function(M, id, key) {
    data.frame(id = id, share = rep(colnames(M), each=nrow(M)), value = c(M), key = key, row.names = NULL)
  }
  DF <- rbind(
    long(OBS, id = idx, key = "Obs"),
    long(P,   id = idx, key = "Pred")
  )
  DFw <- reshape(DF, timevar = "key", idvar = c("id","share"), direction = "wide")
  names(DFw) <- sub("value\\.", "", names(DFw))
  
  # Parity: Pred vs Obs
  g1 <- ggplot(DFw, aes(x = Obs, y = Pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~ share, scales = "free") +
    labs(title = "Parity plot — holdout", x = "Observado", y = "Previsto")
  
  # Resíduos
  DFw$resid <- DFw$Obs - DFw$Pred
  g2 <- ggplot(DFw, aes(x = resid)) +
    geom_histogram(bins = 40) +
    facet_wrap(~ share, scales = "free") +
    labs(title = "Distribuição dos resíduos — holdout", x = "Residuo (Obs − Pred)", y = "Contagem")
  
  print(g1); print(g2)
  invisible(DFw)
}

# ---- 7) Artefatos: salvar objetos/relatórios ----
save_artifacts <- function(dir = "qa_artifacts") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  P <- get_preds_holdout()
  write.csv(P, file.path(dir, "preds_calibrated_holdout.csv"), row.names = FALSE)
  
  # Exporta elasticidades do interior (se existirem)
  if (!is.null(pipeline_results$step7$interior$hicks)) {
    write.csv(pipeline_results$step7$interior$hicks,  file.path(dir, "elasticidades_hicks_interior.csv"))
  }
  if (!is.null(pipeline_results$step7$interior$marshall)) {
    write.csv(pipeline_results$step7$interior$marshall, file.path(dir, "elasticidades_marshall_interior.csv"))
  }
  
  saveRDS(pipeline_results, file.path(dir, "pipeline_results.rds"))
  message(sprintf("Artefatos salvos em: %s", normalizePath(dir)))
}

# ============================================================
# EXECUÇÃO
# ============================================================

cat("== Checagem 1: simplex em preds_calibrated ==\n")
P <- get_preds_holdout()
simp <- check_simplex(P)
cat(
  "  Linhas:", simp$n,
  " | max|sum-1|:", fmt_num(simp$max_abs_sum_err),
  " | #linhas c/ soma≠1:", simp$n_bad_sum,
  " | min(entry):", fmt_num(simp$min_entry),
  " | #negativos:", simp$n_neg, "\n", sep=" "
)

cat("\n== Checagem 2: erro por subgrupo (holdout) ==\n")
dat <- get_dat()
share_cols <- detect_share_cols(dat)
idx <- get_holdout_idx(nrow(P), nrow(dat))

if (length(idx) && length(share_cols) == ncol(P)) {
  OBS <- as.matrix(dat[idx, share_cols, drop=FALSE])
  groups <- intersect(c("regiao", "capital_interior", "quintil_rendapc"), names(dat))
  if (length(groups)) {
    # foco em w1 para viés (além do resumo geral)
    res_all <- by_group_errors(dat[idx, , drop=FALSE], OBS, P, group_vars = groups, focus_share = NULL)
    res_w1  <- by_group_errors(dat[idx, , drop=FALSE], OBS, P, group_vars = groups, focus_share = share_cols[1])
    
    cat("  (Resumo por grupo — todas as shares)\n")
    for (g in names(res_all)) {
      cat("  •", g, "\n")
      print(res_all[[g]])
    }
    cat("\n  (Foco em ", share_cols[1], " — viés médio por grupo)\n", sep="")
    for (g in names(res_w1)) {
      print(res_w1[[g]][, c("group_var","group","share","n","rmse","mae","bias")])
    }
  } else {
    cat("  Variáveis de grupo não encontradas (esperava regiao/capital_interior/quintil_rendapc). Pulando.\n")
  }
} else {
  cat("  Sem índices de holdout (ou #shares não bate). Pulando.\n")
}

cat("\n== Checagem 3: estabilidade do τ* com re-splits ==\n")
stab <- assess_tau_stability(B = 10)
if (!is.null(stab)) {
  tab <- data.frame(tau = stab$tau_grid,
                    freq = sapply(stab$tau_grid, function(t) sum(stab$winners == t, na.rm=TRUE)))
  tab$freq <- as.integer(tab$freq)
  tab$prop <- tab$freq / sum(tab$freq)
  print(tab[order(-tab$freq), ])
}

cat("\n== Checagem 4: bootstrap de elasticidades (ponto interior) ==\n")
boot <- bootstrap_elasticities(B = 200, seed = 321)
if (!is.null(boot)) {
  cat("  Hicks: matrizes de quantis armazenadas em boot$hicks$q (dim: probs x i x j)\n")
  cat("  Marshall: idem em boot$marshall$q\n")
}

cat("\n== Checagem 5: IV (Kleibergen–Paap / Hansen J) numa eq. estrutural ==\n")
iv_out <- quiet_try(iv_tests(dep_share = NULL))
if (inherits(iv_out, "try-error")) {
  cat("  Falhou a estimação IV (verifique pacotes AER/sandwich/lmtest e especificação).\n")
}

cat("\n== Checagem 6: gráficos (parity + resíduos) ==\n")
quiet_try(make_plots())

cat("\n== Checagem 7: salvar artefatos ==\n")
save_artifacts("qa_artifacts")

cat("\n✅ QA concluído.\n")

# ==== HOLDOUT: use a base própria, se existir, e alinhe linhas ====

get_dat_holdout <- function() {
  pr <- pipeline_results$holdout
  # tente bases comuns de holdout geradas no pipeline
  if (!is.null(pr$preprocess$data)) return(pr$preprocess$data)
  if (!is.null(pr$data))            return(pr$data)
  if (!is.null(pr$X))               return(as.data.frame(pr$X))
  # fallback: usa a base de treino (pior cenário)
  message("⚠️  Holdout sem base própria → usando step4$preprocess$data (pode estar errado).")
  get_dat()
}

# se o número de linhas da base do holdout == nrow(preds),
# não precisamos de índices: assumimos 1:n
get_holdout_idx <- function(n_rows_needed, dat_n) {
  pr  <- pipeline_results$holdout
  datH <- try(get_dat_holdout(), silent = TRUE)
  for (nm in c("ids","index","idx","test_idx","test_id","holdout_idx","rows")) {
    if (!is.null(pr[[nm]])) return(as.integer(pr[[nm]]))
  }
  if (!inherits(datH, "try-error") && nrow(datH) == n_rows_needed) {
    return(seq_len(n_rows_needed))
  }
  warning("Índices do holdout não encontrados e não deu para inferir por tamanho.")
  integer(0)
}

# ==== Checagem 2 e 6: use a base do holdout detectada acima ====

run_check_2_and_6 <- function() {
  cat("\n== Checagem 2 (subgrupos) & 6 (gráficos) — com base do holdout ==\n")
  datH <- get_dat_holdout()
  P    <- get_preds_holdout()
  idx  <- get_holdout_idx(nrow(P), nrow(datH))
  if (!length(idx)) { cat("  Sem como alinhar holdout ↔ observados. Pulando.\n"); return(invisible(NULL)) }
  
  share_cols <- detect_share_cols(datH)
  stop_if(length(share_cols) != ncol(P), "N° de shares da base holdout ≠ ncol(preds_calibrated).")
  OBS <- as.matrix(datH[idx, share_cols, drop = FALSE])
  
  groups <- intersect(c("regiao","capital_interior","quintil_rendapc"), names(datH))
  if (!length(groups)) {
    cat("  Variáveis de grupo ausentes no holdout. Pulando.\n")
  } else {
    res_all <- by_group_errors(datH[idx, , drop=FALSE], OBS, P, group_vars = groups, focus_share = NULL)
    res_w1  <- by_group_errors(datH[idx, , drop=FALSE], OBS, P, group_vars = groups, focus_share = share_cols[1])
    
    cat("  (Resumo por grupo — todas as shares)\n")
    for (g in names(res_all)) { cat("  •", g, "\n"); print(res_all[[g]]) }
    
    cat("\n  (Foco em ", share_cols[1], " — viés médio por grupo)\n", sep="")
    for (g in names(res_w1)) print(res_w1[[g]][, c("group_var","group","share","n","rmse","mae","bias")])
  }
  
  # gráficos (parity + resíduos) com a base holdout
  if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2", quiet=TRUE)
  suppressPackageStartupMessages(library(ggplot2))
  
  long <- function(M, id, key) data.frame(id=id, share=rep(colnames(M), each=nrow(M)), value=c(M), key=key)
  DF  <- rbind(long(OBS, id = idx, key="Obs"), long(P, id = idx, key="Pred"))
  DFw <- reshape(DF, timevar="key", idvar=c("id","share"), direction="wide")
  names(DFw) <- sub("value\\.", "", names(DFw))
  DFw$resid <- DFw$Obs - DFw$Pred
  
  print(ggplot(DFw, aes(Obs, Pred)) + geom_point(alpha=.5) + geom_abline(slope=1, intercept=0, linetype=2) +
          facet_wrap(~share, scales="free") + labs(title="Parity plot — holdout"))
  print(ggplot(DFw, aes(resid)) + geom_histogram(bins=40) + facet_wrap(~share, scales="free") +
          labs(title="Distribuição dos resíduos — holdout", x="Resíduo (Obs − Pred)", y="Contagem"))
  invisible(NULL)
}

# ==== Checagem 5 (IV) — versão robusta (limpa IVs e colinearidade) ====

iv_tests_robust <- function(dep_share = NULL, max_iv = 40) {
  if (!requireNamespace("AER", quietly=TRUE)) install.packages("AER", quiet=TRUE)
  if (!requireNamespace("sandwich", quietly=TRUE)) install.packages("sandwich", quiet=TRUE)
  suppressPackageStartupMessages({ library(AER); library(sandwich) })
  
  dat <- get_dat()  # IV costuma ser feito na base de treino
  share_cols <- detect_share_cols(dat)
  if (is.null(dep_share)) dep_share <- share_cols[1]
  stop_if(!(dep_share %in% share_cols), sprintf("Share '%s' não encontrado no data.", dep_share))
  
  # regressoras endógenas (preços) + exógenas
  price_cols <- grep("^ln_preco_com_reforma2[1-6]$", names(dat), value=TRUE)
  stop_if(length(price_cols)==0, "Colunas ln_preco_com_reforma2# não encontradas.")
  x_exog <- intersect(c("z","z2"), names(dat))
  
  # candidatos a instrumentos (CONTÍNUOS; evita dummies IV_d_* que explodem colinearidade)
  Z_regex <- "^(iv_(op|bs|ns|sin|cos|hg)|iv_op0[12]_x_iv_hg[24680]+)"
  Z_all   <- grep(Z_regex, names(dat), value = TRUE)
  stop_if(length(Z_all)==0, "Não encontrei IVs contínuos (iv_op/bs/ns/sin/cos/hg).")
  
  # remove colunas com NA ou variância ~0
  good <- Z_all[vapply(Z_all, function(cn) {
    v <- dat[[cn]]
    ok <- is.numeric(v) && sum(is.finite(v))>0 && (length(unique(v[is.finite(v)]))>1)
    ok
  }, logical(1))]
  Z <- dat[, good, drop=FALSE]
  
  # completa casos nas variáveis do modelo
  keep_vars <- unique(c(dep_share, price_cols, x_exog, colnames(Z)))
  D <- dat[, keep_vars, drop=FALSE]
  D <- D[complete.cases(D), , drop=FALSE]
  
  # tira colinearidade nos IVs (QR)
  qr_keep <- function(M) { qrM <- qr(as.matrix(M)); as.matrix(M)[, qrM$pivot[seq_len(qrM$rank)], drop=FALSE] }
  Z2 <- qr_keep(D[, colnames(Z), drop=FALSE])
  
  # limita a um conjunto enxuto de IVs por correlação com cada preço (proxy p/ 1º estágio forte)
  cor_score <- function(X, y) {
    s <- suppressWarnings(abs(cor(X, y)))
    if (is.matrix(s)) apply(s, 1, max, na.rm=TRUE) else s
  }
  sc <- Reduce(pmax, lapply(price_cols, function(p) cor_score(Z2, D[[p]])))
  iv_order <- names(sort(sc, decreasing=TRUE))
  Z3 <- Z2[, head(iv_order, min(max_iv, ncol(Z2))), drop=FALSE]
  
  # fórmula IV
  f_rhs  <- paste(c(price_cols, x_exog), collapse=" + ")
  f_iv   <- paste(c(colnames(Z3), x_exog), collapse=" + ")
  f_full <- as.formula(sprintf("%s ~ %s | %s", dep_share, f_rhs, f_iv))
  
  fit_iv <- AER::ivreg(f_full, data = D)
  s <- summary(fit_iv, vcov = sandwich::vcovHC(fit_iv, type="HC1"), diagnostics = TRUE)
  
  cat("\n[IV robusto] fórmula:\n  ", deparse(f_full), "\n")
  print(coef(s))
  if (!is.null(s$diagnostics)) {
    cat("\n[Diagnósticos IV]\n"); print(s$diagnostics)
  } else {
    cat("\n(sem tabela de diagnósticos nesta versão do AER)\n")
  }
  invisible(list(fit=fit_iv, sum=s, used_iv=colnames(Z3)))
}

# ======= EXECUTAR as checagens 2,5,6 de novo com o patch =======

run_check_2_and_6()
iv_out2 <- iv_tests_robust(dep_share = NULL, max_iv = 40)

# --- HOTFIX: registrar explicitamente o holdout para #2 e #6 ---

# use um dos dois caminhos:
# 1) se você tem o vetor de índices do holdout relativo à base de treino:
#    register_holdout(idx = <vetor_de_linhas_no_step4_preprocess_data>)
# 2) se você tem a base observada do holdout (mesmo nº de linhas de preds_calibrated):
#    register_holdout(data = holdout_df)

register_holdout <- function(idx = NULL, data = NULL) {
  stop_if(is.null(idx) && is.null(data),
          "Passe 'idx' (inteiros nas linhas do step4$preprocess$data) ou 'data' (data.frame do holdout).")
  pr <- pipeline_results$holdout %||% list()
  if (!is.null(idx))  pr$index <- as.integer(idx)
  if (!is.null(data)) {
    stop_if(!is.data.frame(data), "'data' deve ser data.frame.")
    pr$data <- data
  }
  pipeline_results$holdout <<- pr
  message("Holdout registrado. Rode run_check_2_and_6() novamente.")
}

# (opcional) se você tiver, por ex., 'test_idx' no ambiente:
# register_holdout(idx = test_idx)

# (ou) se você tiver a base observada do holdout:
# register_holdout(data = df_holdout)  # tem que ter colunas w_despesahat1:6 e mesmo nrow das predições
run_check_2_and_6()
iv_tests_all_shares <- function(max_iv = 40) {
  dat <- get_dat()
  shares <- detect_share_cols(dat)
  out <- lapply(shares, function(sj) {
    cat("\n====================\nShare:", sj, "\n")
    iv_tests_robust(dep_share = sj, max_iv = max_iv)
  })
  invisible(out)
}

# usar:
iv_all <- iv_tests_all_shares(40)

# ------------------------------------------------------------------
# IV robusto com seleção de instrumentos e suporte a drop/keep/add
# ------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

.get_dat <- function(){
  pipeline_results$step4$preprocess$data %||%
    pipeline_results$step4$preprocess$dat %||%
    stop("Não achei o data.frame em step4$preprocess.")
}

.iv_filter <- function(cands, drop=NULL, keep=NULL, add=NULL, regex=TRUE){
  u <- function(v) unique(v[!is.na(v) & v != ""])
  cands <- u(cands)
  if (!is.null(keep)) {
    K <- if (regex) unlist(lapply(keep, \(p) grep(p, cands, value=TRUE))) else intersect(cands, keep)
    cands <- u(K)
  }
  if (!is.null(drop)) {
    D <- if (regex) unique(unlist(lapply(drop, \(p) grep(p, cands, value=TRUE)))) else intersect(cands, drop)
    cands <- setdiff(cands, D)
  }
  if (!is.null(add)) cands <- u(c(cands, add))
  cands
}

# calcula score de relevância: média(|cor(Z,ln_preco_k)|) nos 6 preços
.iv_score <- function(dat, Znames, price_cols){
  sc <- sapply(Znames, function(z){
    vv <- sapply(price_cols, function(p) suppressWarnings(cor(dat[[z]], dat[[p]], use="pairwise.complete.obs")))
    mean(abs(vv), na.rm=TRUE)
  })
  sc[is.na(sc)] <- 0
  sort(sc, decreasing = TRUE)
}

iv_tests_robust <- function(dep_share = NULL,
                            max_iv = 40,
                            drop = NULL,      # ex.: c("^iv_sin","^iv_cos") ou c("iv_sin1","iv_cos3")
                            keep = NULL,      # ex.: c("^iv_op","^iv_hg")
                            add  = NULL,      # nomes extra para sempre incluir
                            always_exog = c("z","z2"),
                            print_formula = TRUE){
  
  if (!requireNamespace("AER", quietly=TRUE)) install.packages("AER", quiet=TRUE)
  if (!requireNamespace("sandwich", quietly=TRUE)) install.packages("sandwich", quiet=TRUE)
  if (!requireNamespace("lmtest", quietly=TRUE)) install.packages("lmtest", quiet=TRUE)
  library(AER); library(sandwich); library(lmtest)
  
  dat <- .get_dat()
  
  # shares e preços
  share_cols  <- grep("^w_despesa", names(dat), value = TRUE)
  if (is.null(dep_share)) dep_share <- share_cols[1]
  stopifnot(dep_share %in% share_cols)
  
  price_cols <- grep("^ln_preco_com_reforma2[1-6]$", names(dat), value = TRUE)
  stopifnot(length(price_cols) == 6)
  
  # instrumentos candidatos
  Z_all <- grep("^(iv_|IV_)", names(dat), value=TRUE)
  # tira colunas sem variância
  Z_all <- Z_all[sapply(Z_all, function(z) sd(dat[[z]], na.rm=TRUE) > 0)]
  # aplica filtros
  Z_all <- .iv_filter(Z_all, drop=drop, keep=keep, add=add, regex=TRUE)
  
  # rank por correlação média com os 6 ln_preco
  sc <- .iv_score(dat, Z_all, price_cols)
  Z_ranked <- names(sc)
  Z_sel <- head(Z_ranked, max_iv)
  
  # exógenas do RHS (sempre entram dos dois lados)
  X_exog <- intersect(always_exog, names(dat))
  
  # fórmulas
  f_rhs  <- paste(c(price_cols, X_exog), collapse = " + ")
  f_iv   <- paste(c(Z_sel,      X_exog), collapse = " + ")
  f_full <- as.formula(sprintf("%s ~ %s | %s", dep_share, f_rhs, f_iv))
  
  if (print_formula) {
    cat("\n[IV robusto] fórmula:\n  ", deparse(f_full), "\n", sep = "")
    cat("  #IV selecionados:", length(Z_sel), " de ", length(Z_all), " candidatos\n", sep = "")
    cat("  (primeiros 10): ", paste(head(Z_sel, 10), collapse=", "), "\n", sep = "")
  }
  
  fit <- ivreg(f_full, data = dat)
  V   <- sandwich::vcovHC(fit, type = "HC1")
  s   <- summary(fit, vcov = V, diagnostics = TRUE)
  
  print(s$coefficients)
  
  if (!is.null(s$diagnostics)) {
    cat("\n[Diagnósticos IV]\n")
    print(s$diagnostics)
  } else {
    cat("\nDiagnósticos IV não disponíveis nesta versão do AER.\n")
  }
  
  invisible(list(fit=fit, vcov=V, summary=s, formula=f_full, iv_selected=Z_sel, iv_scores=sc))
}

scan_iv_count_for_w5 <- function(K = c(15,20,25,30,35,40), drop = NULL, keep = NULL){
  lapply(K, function(k){
    cat("\n--- max_iv =", k, "---\n")
    iv_tests_robust("w_despesahat5", max_iv = k, drop = drop, keep = keep)
  })
}

# (a) igual você tentou: tirar seno/cosseno
scan_iv_count_for_w5(drop = c("^iv_sin", "^iv_cos"))

# (b) manter só famílias op/hg (exemplo)
scan_iv_count_for_w5(keep = c("^iv_op", "^iv_hg"))


plots_in_sample <- function(bins = 20,
                            weight_col = c("peso_final","peso","w"),
                            title_suffix = "") {
  # -------- helpers --------
  `%||%` <- function(x, y) if (is.null(x)) y else x
  rmse   <- function(e, w = NULL) {
    if (is.null(w)) sqrt(mean(e^2, na.rm=TRUE))
    else sqrt(sum(w*e^2, na.rm=TRUE)/sum(w, na.rm=TRUE))
  }
  mae    <- function(e, w = NULL) {
    if (is.null(w)) mean(abs(e), na.rm=TRUE)
    else sum(w*abs(e), na.rm=TRUE)/sum(w, na.rm=TRUE)
  }
  get_dat <- function(){
    pipeline_results$step4$preprocess$data %||%
      pipeline_results$step4$preprocess$dat %||%
      stop("Não encontrei o data.frame pré-processado em step4$preprocess.")
  }
  need_pkg <- function(p){ if (!requireNamespace(p, quietly=TRUE)) install.packages(p, quiet=TRUE); suppressPackageStartupMessages(library(p, character.only=TRUE)) }
  
  need_pkg("ggplot2")
  
  # -------- base --------
  dat <- get_dat()
  
  # pesos (se existir)
  wname <- intersect(weight_col, names(dat))
  wvec  <- if (length(wname)) as.numeric(dat[[wname[1]]]) else NULL
  
  # -------- shares observadas --------
  desp_cols <- paste0("despesa", 1:6)
  stopifnot(all(desp_cols %in% names(dat)))
  OBS_amt <- as.matrix(dat[, desp_cols])
  tot_obs <- rowSums(OBS_amt, na.rm = TRUE)
  OBS     <- sweep(OBS_amt, 1, tot_obs, "/")
  
  # -------- shares preditas --------
  # 1) preferir w_despesahat1:6 se existirem
  hat_share_cols <- paste0("w_despesahat", 1:6)
  if (all(hat_share_cols %in% names(dat))) {
    PRED <- as.matrix(dat[, hat_share_cols])
  } else {
    # 2) fallback: despesa#_hat / soma_hat
    hat_amt_cols <- paste0("despesa", 1:6, "_hat")
    stopifnot(all(hat_amt_cols %in% names(dat)))
    HAT_amt <- as.matrix(dat[, hat_amt_cols])
    tot_hat <- rowSums(HAT_amt, na.rm = TRUE)
    PRED    <- sweep(HAT_amt, 1, tot_hat, "/")
  }
  
  # -------- sanity / alinhamento --------
  if (!all(dim(OBS) == dim(PRED)))
    stop("Dimensões de OBS e PRED diferem.")
  
  # linhas válidas (sem NA e com totais positivos)
  ok <- is.finite(rowSums(OBS)) & is.finite(rowSums(PRED)) &
    is.finite(tot_obs) & (tot_obs > 0)
  OBS  <- OBS[ok, , drop=FALSE]
  PRED <- PRED[ok, , drop=FALSE]
  if (!is.null(wvec)) wvec <- wvec[ok]
  
  colnames(OBS)  <- paste0("w", 1:6)
  colnames(PRED) <- paste0("w", 1:6)
  
  # -------- métricas --------
  shares <- colnames(OBS)
  summ <- do.call(rbind, lapply(shares, function(s){
    e <- OBS[, s] - PRED[, s]
    c(share = s,
      rmse   = rmse(e),
      mae    = mae(e),
      rmse_w = if (is.null(wvec)) NA_real_ else rmse(e, wvec),
      mae_w  = if (is.null(wvec)) NA_real_ else mae(e, wvec))
  }))
  summ <- as.data.frame(summ, stringsAsFactors = FALSE)
  numcols <- c("rmse","mae","rmse_w","mae_w")
  summ[ , numcols] <- lapply(summ[ , numcols], as.numeric)
  print(summ)
  
  # -------- dados longos p/ gráficos --------
  long <- function(M, key){
    data.frame(share = rep(colnames(M), each = nrow(M)),
               value = as.numeric(M),
               key   = key,
               row   = rep(seq_len(nrow(M)), times = ncol(M)),
               stringsAsFactors = FALSE)
  }
  DFo <- long(OBS,  "Obs")
  DFp <- long(PRED, "Pred")
  DF  <- merge(DFo, DFp, by = c("share","row"), suffixes = c("_obs","_pred"))
  DF$resid <- DF$value_obs - DF$value_pred
  if (!is.null(wvec)) {
    DF$w <- rep(wvec, times = length(shares))
  }
  
  # -------- gráficos --------
  g1 <- ggplot(DF, aes(x = value_obs, y = value_pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~ share, scales = "free") +
    labs(title = paste0("Parity plot — in-sample", title_suffix),
         x = "Observado", y = "Previsto")
  
  g2 <- ggplot(DF, aes(x = resid)) +
    geom_histogram(bins = bins) +
    facet_wrap(~ share, scales = "free") +
    labs(title = paste0("Resíduos (Obs − Pred) — in-sample", title_suffix),
         x = "Resíduo", y = "Contagem")
  
  print(g1); print(g2)
  
  invisible(list(summary = summ, data = DF, parity_plot = g1, resid_plot = g2))
}
plots_in_sample()
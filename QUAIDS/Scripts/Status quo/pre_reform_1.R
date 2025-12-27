# ============================
# Pacotes
# ============================
library(tidyverse)
library(data.table)
library(haven)
library(micEconAids)
library(systemfit)
library(AER)        # ivreg
library(zoo)
library(splines)
library(lmtest)     # lrtest
library(openxlsx)
setDTthreads(0)

`%||%` <- function(x, y) if (!is.null(x)) x else y
.pstars <- function(p){
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
.bt <- function(v) paste0("`", gsub("`", "\\\\`", v), "`")  # protege nomes em fórmulas

# ============================
# Dados
# ============================
path <- 'C:/Users/x16610962/Downloads/Reproduction/QUAIDS/Data/'
file <- 'banco_analise_AIDS.dta'
df <- read_dta(paste0(path, file))

prices <- as.data.frame(df[, paste0("preco_por_kg", 1:6), drop = FALSE])
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

# ---- ln(preços) para RHS e para LOO/LAG
lnp_cols <- paste0("ln_", priceNames)
for(j in seq_len(K)) df_iv[[lnp_cols[j]]] <- log(df_iv[[priceNames[j]]])

# ============================
# Construção de candidatos a IV
# ============================

# 1) ln(renda) e suas funções
cand_income <- c("renda_total","renda_total_atualhat","renda",
                 "rendatotal","rendimento_total","renda_domiciliar")
inc_name <- cand_income[cand_income %in% names(df_iv)][1]
if (is.na(inc_name)) stop("Nao encontrei coluna de renda. Ajuste 'cand_income'.")
df_iv$ln_income <- log(suppressWarnings(as.numeric(df_iv[[inc_name]])))
df_iv <- df_iv[is.finite(df_iv$ln_income), , drop = FALSE]

ln_inc_c <- as.numeric(scale(df_iv$ln_income, TRUE, TRUE))
rng <- range(df_iv$ln_income, na.rm = TRUE)
ln_inc_u <- (df_iv$ln_income - rng[1]) / (diff(rng)+1e-9)

build_income_IVs <- function(ln_inc_c, ln_inc_u,
                             deg_poly = 8, df_bs = 8, df_ns = 8,
                             fourier_K = 4,
                             hinge_probs = c(.2,.4,.6,.8),
                             add_interactions = TRUE) {
  # Polinômios ortogonais
  OP <- poly(ln_inc_c, degree = deg_poly, raw = FALSE)
  OP <- as.data.frame(OP, stringsAsFactors = FALSE)
  colnames(OP) <- sprintf("iv_op%02d", seq_len(ncol(OP)))
  
  # B-splines e natural splines
  BS <- as.data.frame(bs(ln_inc_c, df = df_bs), stringsAsFactors = FALSE)
  NS <- as.data.frame(ns(ln_inc_c, df = df_ns), stringsAsFactors = FALSE)
  colnames(BS) <- sprintf("iv_bs%02d", seq_len(ncol(BS)))
  colnames(NS) <- sprintf("iv_ns%02d", seq_len(ncol(NS)))
  
  # Fourier
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
  
  # Hinges (ReLU) centrados
  qs <- quantile(ln_inc_c, probs = hinge_probs, na.rm = TRUE)
  H  <- sapply(qs, function(qk) pmax(0, ln_inc_c - qk))
  H  <- sweep(H, 2, colMeans(H, na.rm = TRUE), FUN = "-")
  H  <- as.data.frame(H, stringsAsFactors = FALSE)
  colnames(H) <- sprintf("iv_hg%02d", round(100*hinge_probs))
  
  # Interações moderadas OP x Hinges
  INT <- NULL
  if (add_interactions) {
    op2 <- OP[, 1:min(2, ncol(OP)), drop = FALSE]
    h3  <- H[,  1:min(3, ncol(H)),  drop = FALSE]
    INTlist <- lapply(seq_len(ncol(op2)), function(i){
      out <- sweep(as.matrix(h3), 1, op2[[i]], `*`)
      out <- as.data.frame(out, stringsAsFactors = FALSE)
      colnames(out) <- paste0(colnames(op2)[i], "_x_", colnames(h3))
      out
    })
    INT <- do.call(cbind, INTlist)
  }
  
  # Empilha tudo e garante numérico
  IV <- cbind(OP, BS, NS,
              if (!is.null(FT)) FT,
              H,
              if (!is.null(INT)) INT)
  
  IV <- as.data.frame(lapply(IV, function(col){
    if (is.list(col)) unlist(col, use.names = FALSE) else as.numeric(col)
  }), check.names = FALSE)
  
  # Remove colunas com variância ~ 0
  keep <- vapply(IV, function(v) {
    v <- as.numeric(v)
    is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12
  }, logical(1))
  IV[, keep, drop = FALSE]
}
IV_income <- build_income_IVs(ln_inc_c, ln_inc_u)

# 2) Variáveis categóricas exógenas
pick_first <- function(cands) { cands[cands %in% names(df_iv)][1] %||% NA_character_ }

nm_uf    <- pick_first(c("uf","UF","estado","sigla_uf","cod_uf"))
nm_reg   <- pick_first(c("regiao","região","region","macroregiao","macro_regiao"))
nm_area  <- pick_first(c("area","area_metropolitana","zona","situacao","situacao_domicilio",
                         "urbano_rural","urbano_rural_flag","urbano_rural_cat"))
nm_cap   <- pick_first(c("capital","is_capital","capital_flag"))
nm_int   <- pick_first(c("interior","is_interior","interior_flag"))
nm_mes   <- pick_first(c("mes","month"))
nm_ano   <- pick_first(c("ano","year"))
nm_data  <- pick_first(c("data","date"))

# Deriva tempo
if (!is.na(nm_data)) {
  if (!inherits(df_iv[[nm_data]], "Date")) suppressWarnings(df_iv[[nm_data]] <- as.Date(df_iv[[nm_data]]))
  df_iv$._ym <- as.yearmon(df_iv[[nm_data]])
  df_iv$._mes <- as.integer(format(df_iv[[nm_data]], "%m"))
  df_iv$._ano <- as.integer(format(df_iv[[nm_data]], "%Y"))
} else {
  if (!is.na(nm_mes)) df_iv$._mes <- suppressWarnings(as.integer(df_iv[[nm_mes]]))
  if (!is.na(nm_ano)) df_iv$._ano <- suppressWarnings(as.integer(df_iv[[nm_ano]]))
  if (!is.null(df_iv$._mes) && !is.null(df_iv$._ano)) {
    df_iv$._ym <- as.yearmon(paste(df_iv$._ano, df_iv$._mes), "%Y %m")
  }
}
if (!is.null(df_iv$._mes)) df_iv$._tri <- ((pmax(1, pmin(12, df_iv$._mes))-1) %/% 3) + 1L

# Fatores prontos para model.matrix (~0+)
if (!is.na(nm_uf))   df_iv$f_uf    <- factor(df_iv[[nm_uf]])
if (!is.na(nm_reg))  df_iv$f_reg   <- factor(df_iv[[nm_reg]])
if (!is.na(nm_area)) df_iv$f_area  <- factor(df_iv[[nm_area]])
if (!is.na(nm_cap))  df_iv$f_cap   <- factor(df_iv[[nm_cap]])
if (!is.na(nm_int))  df_iv$f_int   <- factor(df_iv[[nm_int]])
if (!is.null(df_iv$._mes)) df_iv$f_mes <- factor(df_iv$._mes)
if (!is.null(df_iv$._ano)) df_iv$f_ano <- factor(df_iv$._ano)
if (!is.null(df_iv$._ym))  df_iv$f_ym  <- factor(df_iv$._ym)
if (!is.null(df_iv$._tri)) df_iv$f_tri <- factor(df_iv$._tri)

# Interações moderadas
if (!is.null(df_iv$f_uf) && !is.null(df_iv$f_tri)) {
  df_iv$f_uf_tri <- interaction(df_iv$f_uf, df_iv$f_tri, drop = TRUE)
}

# 3) (Nesta versão não usamos LOO/LAG para simplificar — foco no erro de fórmulas)

# 4) Dummies concretas (tudo já em colunas numéricas válidas)
make_dummies <- function(x, prefix){
  f  <- factor(x, exclude = NULL)
  X  <- model.matrix(~ 0 + f)
  colnames(X) <- paste0(prefix, "_", make.names(levels(f), allow_ = FALSE))
  as.data.frame(X, check.names = FALSE)
}
D_uf     <- if (!is.null(df_iv$f_uf))     make_dummies(df_iv$f_uf,     "IV_d_uf") else NULL
D_mes    <- if (!is.null(df_iv$f_mes))    make_dummies(df_iv$f_mes,    "IV_d_mes") else NULL
D_reg    <- if (!is.null(df_iv$f_reg))    make_dummies(df_iv$f_reg,    "IV_d_reg") else NULL
D_cap    <- if (!is.null(df_iv$f_cap))    make_dummies(df_iv$f_cap,    "IV_d_cap") else NULL
D_int    <- if (!is.null(df_iv$f_int))    make_dummies(df_iv$f_int,    "IV_d_int") else NULL
D_ufmes  <- NULL
if (!is.null(D_uf) && !is.null(D_mes)) {
  f_uf  <- factor(df_iv$f_uf,  exclude = NULL)
  f_mes <- factor(df_iv$f_mes, exclude = NULL)
  ufxm  <- interaction(f_uf, f_mes, drop = TRUE, lex.order = TRUE)
  if (nlevels(ufxm) <= 200) D_ufmes <- make_dummies(ufxm, "IV_d_ufmes")
}

# Interações renda (op1/op2) x região/capital
OP_tmp <- as.data.frame(poly(ln_inc_c, degree = 2, raw = FALSE))
colnames(OP_tmp) <- c("IVop1","IVop2")
INT <- NULL
if (!is.null(D_reg)) {
  dR  <- D_reg[, 1:min(3, ncol(D_reg)), drop = FALSE]
  INT <- do.call(cbind, lapply(seq_len(ncol(OP_tmp)), function(i){
    out <- sweep(as.matrix(dR), 1, OP_tmp[[i]], `*`)
    colnames(out) <- paste0(colnames(OP_tmp)[i], "_x_", colnames(dR))
    out
  }))
  INT <- as.data.frame(INT, check.names = FALSE)
}
if (!is.null(D_cap)) {
  INT2 <- sweep(as.matrix(D_cap), 1, OP_tmp[[1]], `*`)
  colnames(INT2) <- paste0(colnames(OP_tmp)[1], "_x_", colnames(D_cap))
  INT <- cbind(INT, as.data.frame(INT2, check.names = FALSE))
}

# 5) Demográficos/regionais contínuos (exógenos)
dem_reg_cont <- names(df_iv)[grepl("^(dem_|reg_)", names(df_iv))]
dem_reg_cont <- dem_reg_cont[!dem_reg_cont %in% c(shareNames, priceNames, "gasto_total_atualhat")]

# 6) Juntar TUDO em colunas numéricas
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

# ---- SANITIZA nomes dos IVs (ESSENCIAL p/ fórmulas)
colnames(IV_cand_df) <- make.names(colnames(IV_cand_df), unique = TRUE)

# Limpa colunas quase constantes
keep_var <- vapply(IV_cand_df, function(v){ vv <- var(v, na.rm = TRUE); is.finite(vv) && vv > 1e-12 }, logical(1))
IV_cand_df <- IV_cand_df[, keep_var, drop = FALSE]

# Cola no df_iv e define o conjunto de candidatos e a base via QR em matriz
df_iv <- cbind(df_iv, IV_cand_df)
iv_all_candidates <- colnames(IV_cand_df)

Z0  <- as.matrix(IV_cand_df)
qr0 <- qr(Z0)
iv_base <- colnames(Z0)[qr0$pivot[seq_len(qr0$rank)]]

req_endog <- length(priceNames) - 1L
cat("Endogenas por equacao (K-1):", req_endog, "\n")
cat("Rank(Z) com IVs amplos:", length(iv_base), "\n")

# ============================
# Estratégias benchmark (SUR)
# ============================
aids_sl <- aidsEst(priceNames = colnames(prices),
                   shareNames = colnames(shares),
                   data = data.frame(prices, shares, x = x),
                   totExpName = "x",
                   method = "LA",
                   priceIndex = "SL",
                   estMethod = "SUR",
                   maxiter = 200)
cat("\n==== Resultado SL-SUR ====\n")
summary(aids_sl); elas(aids_sl)

aids_ls <- aidsEst(priceNames = colnames(prices),
                   shareNames = colnames(shares),
                   data = data.frame(prices, shares, x = x),
                   totExpName = "x",
                   method = "LA",
                   priceIndex = "Ls",
                   estMethod = "SUR",
                   maxiter = 200)
cat("\n==== Resultado Ls-SUR ====\n")
summary(aids_ls); elas(aids_ls)

# =======================
# QUAIDS MANUAL (K preços) – SUR/3SLS (seleção por QR em cada eq.)
# =======================
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
    # monta Z direto da matriz
    inst_in_df <- intersect(instNames, names(df))
    Z <- as.matrix(df[, inst_in_df, drop = FALSE])
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
    gamma[i, base_drop] <- - sum(gamma[i, setdiff(seq_len(K), base_drop)], na.rm = TRUE)
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
    base_dropped = base_name,
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
  cat("Base price dropped from RHS:", object$base_dropped, "\n\n")
  
  co <- coef(fit); se <- sqrt(diag(vcov(fit)))
  t  <- co/se; p <- 2*pt(-abs(t), df = fit$df.residual)
  lab <- names(co)
  
  relabel <- sub("^([^.:]+)[.:]\\(Intercept\\)$", "alpha \\1", lab)
  relabel <- sub("([.:]|^)z2$", "lambda", sub("([.:]|^)z$", "beta", relabel))
  relabel <- sub("([.:]|^)ln_", "gamma ", relabel)
  relabel <- sub("^([^ ]+)[.:]", "\\1 ", relabel)
  
  tab <- data.frame(Estimate = co, `Std. Error` = se, `t value` = t, `Pr(>|t|)` = p,
                    row.names = relabel, check.names = FALSE)
  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE)
  
  cat("\nR-squared Values of expenditure shares:\n")
  r2 <- sapply(seq_along(object$shareNames), function(i){
    eqLab <- object$eqLabels[i]
    y  <- model.response(model.frame(fit$eq[[eqLab]]))
    yhat <- fitted(fit$eq[[eqLab]])
    1 - sum((y - yhat)^2)/sum( (y - mean(y))^2 )
  })
  names(r2) <- object$shareNames
  print(r2)
  
  cat("\nRHS efetivo por equação:\n")
  for (i in seq_along(object$shareNames)) {
    cat("-", object$shareNames[i], ": ", paste(object$kept_by_eq[[i]], collapse = " + "), "\n", sep = "")
  }
  cat("\nRegressores descartados por equação (se houver):\n")
  for (i in seq_along(object$shareNames)) {
    if (length(object$dropped_by_eq[[i]])>0)
      cat("-", object$shareNames[i], ": ", paste(object$dropped_by_eq[[i]], collapse = " + "), "\n", sep = "")
  }
  invisible(tab)
}

# =======================
# QUAIDS MANUAL (K–1 preços) – SUR/3SLS
# =======================
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
  
  # Se for 3SLS, garanta que todos os instrumentos existam em df (traga de instData se faltar)
  if (estMethod == "3SLS" && !is.null(instNames)) {
    instNames <- instNames[!is.na(instNames) & nzchar(instNames)]
    
    miss <- setdiff(instNames, names(df))
    if (length(miss) > 0) {
      if (is.null(instData)) stop("Para 3SLS, passe instData ...")
      miss <- intersect(miss, names(instData))   # <- só o que existe em instData
      if (!length(miss)) stop("Nenhum instrumento de instNames existe em instData.")
      add <- as.data.frame(instData[, miss, drop = FALSE])  # não precisa [ok, ] se instData já é df_iv filtrado
      df  <- cbind(df, add)
    }
  }
  
  keep_sh_idx <- setdiff(seq_len(K), omit_share)
  keep_ln_idx <- setdiff(seq_len(K), drop_price)
  
  share_km1 <- shareNames[keep_sh_idx]
  ln_km1    <- lnPcols[keep_ln_idx]
  
  # Constrói RHS com reformulate (evita problemas de nomes não-sintáticos)
  rhs_terms <- c(ln_km1, "z", if (use_z2) "z2" else NULL)
  
  eqLabels <- gsub("[ _]", "", share_km1)
  eqs <- vector("list", length(share_km1)); names(eqs) <- eqLabels
  for (i in seq_along(share_km1)) {
    eqs[[i]] <- stats::reformulate(rhs_terms, response = share_km1[i])
  }
  
  if (estMethod == "SUR") {
    fit <- systemfit::systemfit(eqs, data = df, method = "SUR", maxit = maxiter)
  } else {
    if (is.null(instNames) || length(instNames) == 0) stop("Sem instrumentos (instNames) fornecidos.")
    
    # Posto-cheio de Z usando os TERMOS ORIGINAIS dos instrumentos (fatores serão expandidos aqui)
    Z <- model.matrix(stats::reformulate(instNames, intercept = FALSE), data = df)
    keepZ <- apply(Z, 2, function(v) is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12)
    Z <- Z[, keepZ, drop = FALSE]
    qrZ <- qr(Z)
    if (qrZ$rank < 1L) stop("Posto de Z < 1 após filtros de variância.")
    
    # Fórmula de instrumentos deve usar os TERMOS ORIGINAIS (não os nomes expandidos)
    inst_fml <- stats::reformulate(instNames, intercept = FALSE)
    
    fit <- systemfit::systemfit(eqs, data = df, method = "3SLS",
                                inst = inst_fml, maxit = maxiter)
  }
  
  cf <- coef(fit); se <- sqrt(diag(vcov(fit)))
  pick_eq <- function(eqLabel){
    nm <- names(cf); idx <- grep(paste0("^", eqLabel, "[\\.:]"), nm)
    list(cf = cf[idx], se = se[idx], nm = nm[idx])
  }
  
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
  
  # Reconstrução por adding-up/homogeneidade
  alpha[omit_share]  <- 1 - sum(alpha[keep_sh_idx],  na.rm = TRUE)
  beta[omit_share]   <- -   sum(beta[ keep_sh_idx],  na.rm = TRUE)
  lambda[omit_share] <- -   sum(lambda[keep_sh_idx], na.rm = TRUE)
  for (j in seq_len(K)) gamma[omit_share, j] <- - sum(gamma[keep_sh_idx, j], na.rm = TRUE)
  adj <- sum(gamma[omit_share, ], na.rm = TRUE)
  if (is.finite(adj) && abs(adj) > 1e-10) gamma[omit_share, ] <- gamma[omit_share, ] - adj/K
  
  obj <- list(
    call         = match.call(),
    method       = "QUAIDS",
    priceIndex   = priceIndex,
    estMethod    = estMethod,
    maxiter      = maxiter,
    priceNames   = priceNames,
    shareNames   = shareNames,
    omit_share   = shareNames[omit_share],
    drop_price   = priceNames[drop_price],
    coef         = list(alpha = alpha, beta = beta, lambda = lambda, gamma = gamma),
    fit          = fit,
    data         = df,
    resid        = residuals(fit),
    rhs          = paste(rhs_terms, collapse = " + ")
  )
  class(obj) <- c("quaids_manual","aidsEst","list")
  obj
}


summary.quaids_manual <- function(object, ...){
  fit <- object$fit
  cat("Demand analysis with the QUADRATIC Almost Ideal Demand System (QUAIDS)\n")
  cat("Estimation Method:", object$estMethod, "\n")
  cat("Price Index:", object$priceIndex, "\n")
  cat("Omitted share equation:", object$omit_share, "\n")
  cat("Dropped price from RHS:", object$drop_price, "\n\n")
  
  co <- coef(fit); se <- sqrt(diag(vcov(fit)))
  t  <- co/se; p <- 2*pt(-abs(t), df = fit$df.residual)
  lab <- names(co)
  
  relabel <- sub("^([^.:]+)[.:]\\(Intercept\\)$", "alpha \\1", lab)
  relabel <- sub("([.:]|^)z2$", "lambda", sub("([.:]|^)z$", "beta", relabel))
  relabel <- sub("([.:]|^)ln_", "gamma ", relabel)
  relabel <- sub("^([^ ]+)[.:]", "\\1 ", relabel)
  
  tab <- data.frame(Estimate = co, `Std. Error` = se, `t value` = t, `Pr(>|t|)` = p,
                    row.names = relabel, check.names = FALSE)
  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE)
  invisible(tab)
}

# =======================
# Funções auxiliares (reconstrução/elasticidades)
# =======================
rebuild_quaids_coefs <- function(fit_q){
  stopifnot(inherits(fit_q, "quaids_manual"))
  fit        <- fit_q$fit
  priceNames <- fit_q$priceNames
  shareNames <- fit_q$shareNames
  K          <- length(priceNames)
  omit_name  <- fit_q$omit_share %||% fit_q$base_dropped
  omit_idx   <- match(omit_name, shareNames)
  keep_idx   <- setdiff(seq_len(K), omit_idx)
  
  cf <- coef(fit)
  alpha  <- setNames(numeric(K), shareNames)
  beta   <- setNames(numeric(K), shareNames)
  lambda <- setNames(numeric(K), shareNames)
  gamma  <- matrix(0, K, K, dimnames = list(shareNames, priceNames))
  
  lnPcols   <- paste0("ln_", priceNames)
  map_ln_j  <- setNames(match(sub("^ln_", "", lnPcols), priceNames), lnPcols)
  
  for (i in seq_along(keep_idx)) {
    si  <- keep_idx[i]
    lab <- gsub("[ _]", "", shareNames[si])
    nm_eq <- names(cf)[grepl(paste0("^", lab, "[_\\.:]"), names(cf))]
    aidx <- grep("\\(Intercept\\)", nm_eq, value=TRUE); if (length(aidx)==1) alpha[si] <- unname(cf[aidx])
    bidx <- grep("z$",  nm_eq, value=TRUE); if (length(bidx)==1) beta[si]   <- unname(cf[bidx])
    lidx <- grep("z2$", nm_eq, value=TRUE); if (length(lidx)==1) lambda[si] <- unname(cf[lidx])
    ln_here <- nm_eq[grepl("ln_", nm_eq)]
    for (nm in ln_here){
      var <- sub(paste0("^", lab, "[_\\.:]"), "", nm)
      j   <- map_ln_j[[var]]
      if (!is.na(j)) gamma[si, j] <- unname(cf[nm])
    }
  }
  alpha[omit_idx]  <- 1 - sum(alpha[keep_idx],  na.rm=TRUE)
  beta[omit_idx]   <- -   sum(beta[keep_idx],   na.rm=TRUE)
  lambda[omit_idx] <- -   sum(lambda[keep_idx], na.rm=TRUE)
  for (j in seq_len(K)) gamma[omit_idx, j] <- - sum(gamma[keep_idx, j], na.rm=TRUE)
  
  fit_q$coef <- list(alpha=alpha, beta=beta, lambda=lambda, gamma=gamma)
  fit_q
}

fix_gamma_dropped_col <- function(fit_q){
  stopifnot(inherits(fit_q, "quaids_manual"))
  j0 <- match(fit_q$drop_price %||% fit_q$base_dropped, fit_q$priceNames)
  G  <- fit_q$coef$gamma
  for(i in seq_len(nrow(G))) G[i, j0] <- - sum(G[i, -j0], na.rm = TRUE)
  fit_q$coef$gamma <- G
  fit_q
}

sanitize_quaids_coefs <- function(fit_q){
  stopifnot(inherits(fit_q, "quaids_manual"))
  fit        <- fit_q$fit
  priceNames <- fit_q$priceNames
  shareNames <- fit_q$shareNames
  K          <- length(priceNames)
  omit_name  <- fit_q$omit_share
  omit_idx   <- match(omit_name, shareNames)
  keep_idx   <- setdiff(seq_len(K), omit_idx)
  eqLabels   <- gsub("[ _]", "", shareNames[keep_idx])
  
  cf <- coef(fit)
  pick_eq <- function(eqLabel){
    nm  <- names(cf); idx <- grep(paste0("^", eqLabel, "[_\\.:]"), nm)
    list(cf = cf[idx], nm = nm[idx])
  }
  
  alpha  <- setNames(rep(NA_real_, K), shareNames)
  beta   <- setNames(rep(NA_real_, K), shareNames)
  lambda <- setNames(rep(NA_real_, K), shareNames)
  gamma  <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  
  lnPcols <- paste0("ln_", priceNames)
  map_ln_j <- setNames(match(sub("^ln_", "", lnPcols), priceNames), lnPcols)
  
  for (i in seq_along(keep_idx)) {
    si  <- keep_idx[i]
    lab <- eqLabels[i]
    pi  <- pick_eq(lab)
    aidx <- grep("\\(Intercept\\)", pi$nm)
    bidx <- grep("([\\.:_]|^)z$",  pi$nm)
    lidx <- grep("([\\.:_]|^)z2$", pi$nm)
    if (length(aidx)==1) alpha[ si ] <- unname(pi$cf[aidx])
    if (length(bidx)==1) beta[  si ] <- unname(pi$cf[bidx])
    if (length(lidx)==1) lambda[si ] <- unname(pi$cf[lidx])
    
    ln_here <- sub(paste0("^", lab, "[_\\.:]"), "", pi$nm)
    ln_here <- ln_here[grepl("^ln_", ln_here)]
    for (lnv in ln_here) {
      j <- map_ln_j[[lnv]]
      if (!is.na(j)) {
        gidx <- grep(paste0("(^|[\\.:_])", lnv, "$"), pi$nm)
        if (length(gidx)==1) gamma[si, j] <- unname(pi$cf[gidx])
      }
    }
  }
  alpha[omit_idx]  <- 1 - sum(alpha[keep_idx],  na.rm = TRUE)
  beta[omit_idx]   <- -   sum(beta[ keep_idx],  na.rm = TRUE)
  lambda[omit_idx] <- -   sum(lambda[keep_idx], na.rm = TRUE)
  for (j in seq_len(K)) gamma[omit_idx, j] <- - sum(gamma[keep_idx, j], na.rm = TRUE)
  gamma[is.na(gamma)] <- 0
  fit_q$coef <- list(alpha = alpha, beta = beta, lambda = lambda, gamma = gamma)
  fit_q
}

enforce_theory <- function(H, w){
  K <- length(w); W <- diag(w)
  # 1) simetriza em termos de W*H
  A  <- W %*% H
  A  <- 0.5 * (A + t(A))
  Hs <- solve(W, A)  # W^{-1} * A
  # 2) impõe adding-up via correção de posto-1 preservando simetria (em W*H)
  r  <- as.numeric(Hs %*% w)
  denom <- as.numeric(t(w) %*% w)
  C  <- (r %*% t(w)) / denom
  Hc <- Hs - C - solve(W, t(C) %*% W)
  Hc
}

.enforce_hicks_adding <- function(EH, w){
  ww <- sum(w^2)
  if (!is.finite(ww) || ww <= 0) return(EH)
  adj <- as.numeric(EH %*% w) / ww        # vetor (Kx1)
  EH - adj %o% w                          # remove componente na direção de w
}

.symmetrize_slutsky <- function(EH, w){
  # projeta M = diag(w) * EH para o cone simétrico e volta
  if (any(!is.finite(w)) || any(w <= 0)) return(EH)
  W  <- diag(w)
  M  <- W %*% EH
  Ms <- (M + t(M))/2
  solve(W, Ms)                            # W^{-1} * Ms
}

.enforce_cournot_marshall <- function(EM, eta, w){
  ww <- sum(w^2)
  if (!is.finite(ww) || ww <= 0) return(EM)
  target <- -(eta - 1)                    # vetor (Kx1)
  cur    <- as.numeric(EM %*% w)          # (Kx1)
  a      <- (target - cur) / ww           # (Kx1)
  EM + a %o% w
}

.check_identities <- function(EH, EM, eta, w){
  add_hicks <- as.numeric(EH %*% w)                     # ~ 0
  cour_mar  <- cbind(lhs = as.numeric(EM %*% w),        # ~ -(eta-1)
                     rhs = -(eta - 1),
                     diff = as.numeric(EM %*% w) + (eta - 1))
  W  <- diag(w); symm_err <- max(abs(W %*% EH - t(W %*% EH)))
  list(adding_hicks = add_hicks, cournot_marshall = cour_mar,
       symmetry_err = symm_err)
}

# --- ELASTICIDADES com Stone iterativo, diff central e renormalização
elas_quaids_manual <- function(
    fit_q, x, p = NULL, eps = 1e-4, enforce_shares = TRUE,
    enforce_hicks = TRUE,          # NOVO: impõe adding-up em Hicks (default: TRUE)
    enforce_symmetry = FALSE,      # NOVO: projeta simetria de Slutsky (default: FALSE)
    enforce_cournot = FALSE,       # NOVO: força Cournot em Marshall (default: FALSE)
    return_checks = TRUE           # NOVO: devolve diagnóstico das identidades
){
  stopifnot(inherits(fit_q, "quaids_manual"))
  priceNames <- fit_q$priceNames
  shareNames <- fit_q$shareNames
  K <- length(priceNames)
  
  alpha  <- fit_q$coef$alpha
  beta   <- fit_q$coef$beta
  lambda <- fit_q$coef$lambda
  gamma  <- fit_q$coef$gamma
  
  if (is.null(p)) {
    ln_cols <- paste0("ln_", priceNames)
    stopifnot(all(ln_cols %in% names(fit_q$data)))
    p <- exp(colMeans(fit_q$data[, ln_cols, drop = FALSE], na.rm = TRUE))
    names(p) <- priceNames
  } else {
    stopifnot(length(p)==K); if (is.null(names(p))) names(p) <- priceNames
  }
  
  wbar <- colMeans(fit_q$data[, shareNames, drop = FALSE], na.rm = TRUE)
  wbar <- wbar / sum(wbar)
  
  w_hat_fun <- function(p_vec, x_val){
    if (fit_q$priceIndex == "Ls") {
      z <- log(x_val) - sum(wbar * log(p_vec))
      w <- as.numeric(alpha + gamma %*% log(p_vec) + beta * z + lambda * z^2)
    } else if (fit_q$priceIndex == "S") {
      z_ls <- log(x_val) - sum(wbar * log(p_vec))
      w1   <- as.numeric(alpha + gamma %*% log(p_vec) + beta * z_ls + lambda * z_ls^2)
      if (enforce_shares) { w1 <- w1 + (1 - sum(w1))/K }
      z_st <- log(x_val) - sum(w1 * log(p_vec))
      w    <- as.numeric(alpha + gamma %*% log(p_vec) + beta * z_st + lambda * z_st^2)
    } else stop("Índice de preço desconhecido.")
    if (enforce_shares) { w <- w + (1 - sum(w))/K }
    names(w) <- shareNames
    w
  }
  
  w0 <- w_hat_fun(p, x)
  if (any(!is.finite(w0))) stop("Shares previstos ficaram NA/Inf; verifique coeficientes e p/x.")
  
  # z conforme índice
  if (fit_q$priceIndex == "Ls") {
    z0 <- log(x) - sum(wbar * log(p))
  } else {
    z0 <- log(x) - sum(w0 * log(p))
  }
  eta <- 1 + (beta + 2*lambda*z0) / w0
  names(eta) <- shareNames
  
  # Marshall via variação finita em preços
  q0 <- (w0 * x) / p
  E_M <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  for (j in seq_len(K)) {
    p_up <- p; p_up[j] <- p_up[j] * (1 + eps)
    w_up <- w_hat_fun(p_up, x)
    q_up <- (w_up * x) / p_up
    E_M[, j] <- (log(q_up) - log(q0)) / (log(p_up[j]) - log(p[j]))
  }
  E_H <- E_M + outer(eta, w0)   # identidade de Slutsky em elasticidades
  
  # ====== NOVO: correções opcionais ======
  if (enforce_symmetry) {
    E_H <- .symmetrize_slutsky(E_H, w0)        # 1) simetria (opcional)
  }
  if (enforce_hicks) {
    E_H <- .enforce_hicks_adding(E_H, w0)      # 2) adding-up em Hicks (default)
  }
  E_M <- E_H - outer(eta, w0)                  # garante ligação Hicks–Marshall
  
  if (enforce_cournot) {
    # 3) Cournot em Marshall; depois atualiza Hicks a partir de Marshall
    E_M <- .enforce_cournot_marshall(E_M, eta, w0)
    E_H <- E_M + outer(eta, w0)
    # (não reimpomos Hicks após Cournot para evitar conflitos)
  }
  # =======================================
  
  out <- list(at = list(p = setNames(p, priceNames), x = x, w = w0, z = z0),
              expenditure = eta, marshall = E_M, hicks = E_H)
  
  if (return_checks) {
    out$checks <- .check_identities(E_H, E_M, eta, w0)
  }
  out
}
# =======================
# Fallback de estimação 3SLS com conjuntos de IVs escalonados
# =======================
safe_fit_quaids_3sls_any <- function(inst_sets, data_Z, ..., use_z2_first = TRUE){
  last_err <- NULL
  # usa req_endog se existir, senão assume 1
  req_endog_local <- tryCatch(get("req_endog", inherits = TRUE), error = function(e) NULL)
  if (is.null(req_endog_local)) req_endog_local <- 1L
  
  for (S in inst_sets){
    # 1) mantém só nomes válidos que EXISTEM em data_Z
    S0 <- unique(S)
    S0 <- S0[!is.na(S0) & nzchar(S0)]
    S0 <- intersect(S0, names(data_Z))
    if (length(S0) == 0L) { last_err <- "Conjunto de IVs vazio após interseção com data_Z"; next }
    
    # 2) constroi Z com reformulate (expande fatores corretamente)
    Z0 <- try(model.matrix(stats::reformulate(S0, intercept = FALSE), data = data_Z), silent = TRUE)
    if (inherits(Z0, "try-error")) { last_err <- Z0; next }
    
    # 3) filtra variância ~0 e checa posto
    keep <- apply(Z0, 2, function(v) is.finite(var(v, na.rm = TRUE)) && var(v, na.rm = TRUE) > 1e-12)
    Z0 <- Z0[, keep, drop = FALSE]
    if (ncol(Z0) == 0L) { last_err <- "Todas as colunas de Z0 foram filtradas por variância"; next }
    
    rZ <- qr(Z0)$rank
    if (!is.finite(rZ) || rZ < req_endog_local) { last_err <- "rank(Z) < nº de endógenas (K-1)"; next }
    
    # 4) tenta com/sem z2
    if (use_z2_first){
      ans <- try(fit_quaids_manual_km1(instNames = S0, instData = data_Z, use_z2 = TRUE,  ...), silent = TRUE)
      if (!inherits(ans, "try-error")) return(ans)
      last_err <- ans
      ans2 <- try(fit_quaids_manual_km1(instNames = S0, instData = data_Z, use_z2 = FALSE, ...), silent = TRUE)
      if (!inherits(ans2, "try-error")) return(ans2)
      last_err <- ans2
    } else {
      ans <- try(fit_quaids_manual_km1(instNames = S0, instData = data_Z, use_z2 = FALSE, ...), silent = TRUE)
      if (!inherits(ans, "try-error")) return(ans)
      last_err <- ans
      ans2 <- try(fit_quaids_manual_km1(instNames = S0, instData = data_Z, use_z2 = TRUE,  ...), silent = TRUE)
      if (!inherits(ans2, "try-error")) return(ans2)
      last_err <- ans2
    }
  }
  stop(last_err)
}

# Conjuntos escalonados de IVs
head_grep <- function(pat, x, n) head(grep(pat, x, value = TRUE), n)

iv_set_core <- unique(c(
  head_grep("^iv_op",   names(df_iv), 3),   # use df_iv, não IV_income
  head_grep("^IV_d_uf_",  names(df_iv), 5),
  head_grep("^IV_d_mes_", names(df_iv), 6)
))
iv_set_core <- iv_set_core[!is.na(iv_set_core) & nzchar(iv_set_core)]

iv_set_mid <- unique(c(
  iv_set_core,
  grep("^IV_d_reg_", names(df_iv), value = TRUE),
  grep("^IV_d_cap_", names(df_iv), value = TRUE),
  grep("^IV_d_int_", names(df_iv), value = TRUE),
  grep("^IV_d_ufmes_", names(df_iv), value = TRUE),
  if (!is.null(INT)) colnames(INT) else character(0)
))
iv_set_full <- unique(intersect(iv_all_candidates, names(df_iv)))

# =======================
# Estimação 3SLS
# =======================
Z_chk <- model.matrix(stats::reformulate(iv_all_candidates, intercept = FALSE), data = df_iv)
cat("Rank(Z) final (cheque):", qr(Z_chk)$rank, "\n")

# Tentativa direta com base ampla
fit_q_3sls <- try(fit_quaids_manual_km1(
  prices     = df_iv[, priceNames, drop = FALSE],
  shares     = df_iv[, shareNames, drop = FALSE],
  x          = df_iv[["gasto_total_atualhat"]],
  priceIndex = "Ls",
  estMethod  = "3SLS",
  omit_share = 1,
  drop_price = 1,
  instNames  = iv_set_core,
  instData   = df_iv,
  use_z2     = TRUE,
  maxiter    = 500
), silent = TRUE)

if (inherits(fit_q_3sls, "try-error")) {
  fit_q_3sls <- safe_fit_quaids_3sls_any(
    inst_sets   = list(iv_set_full, iv_set_mid, iv_set_core),
    data_Z      = df_iv,  # <-- ADICIONE ISTO
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
summary(fit_q_3sls)

# Reconstrução/ajustes e elasticidades
fit_q <- rebuild_quaids_coefs(fit_q_3sls)
fit_q <- fix_gamma_dropped_col(fit_q)
fit_q <- sanitize_quaids_coefs(fit_q)

# Comparação com Ls-SUR
lrtest(aids_ls$est, fit_q$fit)

x_mean <- mean(df$gasto_total_atualhat, na.rm = TRUE)

# a) Correção padrão (Hicks adding-up; sem Cournot; simetria opcional)
E_fix <- elas_quaids_manual(
  fit_q, x = x_mean,
  enforce_hicks = TRUE,
  enforce_symmetry = TRUE,   # ligue se quiser simetria
  enforce_cournot = FALSE
)
E_fix$checks                         # diagnóstico rápido
diag(E_fix$hicks)                    # próprias compensadas (espera-se < 0)

# b) Se quiser forçar Cournot em Marshall (aceitando leve desvio em Hicks):
E_cour <- elas_quaids_manual(
  fit_q, x = x_mean,
  enforce_hicks = TRUE,      # impõe Hicks antes
  enforce_symmetry = FALSE,
  enforce_cournot = TRUE     # e então Cournot em Marshall
)
E_cour$checks
# =======================
# Bootstrap (opcional)
# =======================
.aggregate_boot_array <- function(arr3){
  est <- apply(arr3, c(1,2), mean, na.rm = TRUE)
  se  <- apply(arr3, c(1,2), sd,   na.rm = TRUE)
  n   <- apply(!is.na(arr3), c(1,2), sum)
  z   <- est / se
  p   <- 2*pnorm(-abs(z))
  list(est=est, se=se, z=z, p=p, n=n)
}

.aggregate_boot_vec <- function(mat2){
  est <- apply(mat2, 1, mean, na.rm = TRUE)
  se  <- apply(mat2, 1, sd,   na.rm = TRUE)
  n   <- apply(!is.na(mat2), 1, sum)
  z   <- est / se
  p   <- 2*pnorm(-abs(z))
  list(est=est, se=se, z=z, p=p, n=n)
}
print_elas_table <- function(est, se, p, digits=3, title=NULL){
  if (!is.null(title)) cat("\n", title, "\n", sep="")
  stars <- .pstars(p)
  fmt   <- function(x) formatC(x, format="f", digits=digits)
  est_s <- fmt(est); se_s  <- paste0("(", fmt(se), ")")
  out   <- matrix(paste0(est_s, " ", se_s, stars), nrow=nrow(est),
                  dimnames = dimnames(est))
  print(noquote(out)); invisible(out)
}

bootstrap_quaids_elasticities <- function(
    prices, shares, x,
    B = 200, seed = 123,
    priceIndex = c("S","Ls")[1],
    estMethod  = c("SUR","3SLS"),
    omit_share = 1,
    drop_price = 1,
    maxiter    = 200,
    verbose    = TRUE,
    instNames  = NULL,
    instData   = NULL){
  set.seed(seed); estMethod <- match.arg(estMethod)
  
  K  <- ncol(prices); sn <- colnames(shares); pn <- colnames(prices)
  E_M_arr <- array(NA_real_, dim = c(K, K, B), dimnames = list(sn, pn, NULL))
  E_H_arr <- array(NA_real_, dim = c(K, K, B), dimnames = list(sn, pn, NULL))
  ETA_mat <- matrix(NA_real_, nrow = K, ncol = B, dimnames = list(sn, NULL))
  
  N <- nrow(prices); idx_all <- seq_len(N)
  
  success <- 0L  # contador de réplicas válidas
  
  for (b in 1:B){
    if (verbose && (b %% max(1, round(B/10))==0)) cat("Bootstrap", b, "de", B, "\n")
    
    idx <- sample(idx_all, size = N, replace = TRUE)
    prices_b <- prices[idx, , drop = FALSE]
    shares_b <- shares[idx, , drop = FALSE]
    x_b      <- x[idx]
    inst_b   <- if (!is.null(instData)) instData[idx, , drop = FALSE] else NULL
    if (is.null(inst_b)) inst_b <- data.frame()  # garantia
    
    # usa o fallback que seleciona IVs e checa posto na PRÓPRIA réplica
    fit_b <- try(
      safe_fit_quaids_3sls_any(
        inst_sets = list(iv_set_full, iv_set_mid, iv_set_core),
        data_Z    = inst_b,
        prices     = prices_b,
        shares     = shares_b,
        x          = x_b,
        priceIndex = priceIndex,
        estMethod  = estMethod,
        omit_share = omit_share,
        drop_price = drop_price,
        maxiter    = maxiter,
        use_z2_first = TRUE
      ),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error")) next
    
    fit_b <- try(rebuild_quaids_coefs(fit_b), silent = TRUE); if (inherits(fit_b, "try-error")) next
    fit_b <- try(fix_gamma_dropped_col(fit_b), silent = TRUE); if (inherits(fit_b, "try-error")) next
    
    x_eval <- mean(x_b, na.rm = TRUE)
    e_b <- try(elas_quaids_manual(fit_b, x = x_eval), silent = TRUE)
    if (inherits(e_b, "try-error")) next
    
    E_M_arr[ , , b] <- e_b$marshall
    E_H_arr[ , , b] <- e_b$hicks
    ETA_mat[ , b]   <- e_b$expenditure
    success <- success + 1L
  }
  
  if (verbose) cat("Réplicas válidas:", success, "de", B, "\n")
  
  
  mar <- .aggregate_boot_array(E_M_arr)
  hix <- .aggregate_boot_array(E_H_arr)
  expd<- .aggregate_boot_vec(ETA_mat)
  
  out <- list(marshall=mar, hicks=hix, expenditure=expd, B=B,
              estMethod=estMethod, instNames=instNames)
  class(out) <- "quaids_boot_elas"
  out
}

# Exemplo (opcional):
# boot_res_3sls <- bootstrap_quaids_elasticities(
#   prices     = df_iv[, priceNames],
#   shares     = df_iv[, shareNames],
#   x          = df_iv$gasto_total_atualhat,
#   B          = 500,
#   seed       = 123,
#   priceIndex = "Ls",
#   estMethod  = "3SLS",
#   omit_share = 1,
#   drop_price = 1,
#   maxiter    = 300,
#   instNames  = iv_set_core,
#   instData   = df_iv
# )
# print(boot_res_3sls$marshall)

 # se ainda não tiver:
 lnp_cols <- paste0("ln_", priceNames)
 
 # pesos médios dos shares (normalizados)
 wbar <- colMeans(df_iv[, shareNames], na.rm = TRUE)
 wbar <- wbar / sum(wbar)
 
 # ln P^Ls, z e z2
 df_iv$z  <- log(as.numeric(df_iv$gasto_total_atualhat)) -
   as.numeric(as.matrix(df_iv[, lnp_cols, drop = FALSE]) %*% wbar)
 df_iv$z2 <- df_iv$z^2
 
 check_iv_strength <- function(data_Z, instNames, drop_price, priceNames,
                               use_z2 = TRUE, shareNames = NULL) {
   stopifnot(is.data.frame(data_Z))
   
   # 1) garante ln_ preços
   lnp_cols <- paste0("ln_", priceNames)
   if (!all(lnp_cols %in% names(data_Z))) {
     if (!all(priceNames %in% names(data_Z))) {
       stop("Nem ln_* nem os preços brutos estão em data_Z.")
     }
     for (j in seq_along(priceNames)) {
       data_Z[[paste0("ln_", priceNames[j])]] <- log(as.numeric(data_Z[[priceNames[j]]]))
     }
   }
   
   # 2) garante z/z2 (índice Ls) se faltar
   if (!("z" %in% names(data_Z))) {
     if (is.null(shareNames)) stop("Passe shareNames para calcular z quando ele faltar.")
     wbar <- colMeans(data_Z[, shareNames, drop = FALSE], na.rm = TRUE)
     wbar <- wbar / sum(wbar)
     lnP  <- as.numeric(as.matrix(data_Z[, lnp_cols, drop = FALSE]) %*% wbar)
     data_Z$z <- log(as.numeric(data_Z$gasto_total_atualhat)) - lnP
   }
   if (use_z2 && !("z2" %in% names(data_Z))) data_Z$z2 <- data_Z$z^2
   
   # 3) define conjuntos
   ln_all  <- paste0("ln_", priceNames)
   ln_keep <- setdiff(ln_all, paste0("ln_", priceNames[drop_price]))
   exogs   <- c("z", if (use_z2) "z2" else NULL)
   
   instNames <- unique(instNames[instNames %in% names(data_Z)])
   if (length(instNames) == 0) stop("Nenhum instrumento de instNames está em data_Z.")
   
   # 4) higieniza nomes
   original_names <- names(data_Z)
   safe_names     <- make.names(original_names, unique = TRUE)
   colnames(data_Z) <- safe_names
   map <- setNames(safe_names, original_names)
   
   ln_keep_s <- unname(map[ln_keep])
   exogs_s   <- unname(map[exogs])
   inst_s    <- unname(map[instNames])
   
   mk_rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"
   
   # 5) loop: first-stage OLS por variável endógena
   rows <- lapply(seq_along(ln_keep_s), function(i){
     y_s <- ln_keep_s[i]; y_o <- ln_keep[i]
     
     cols_need_full <- unique(c(y_s, exogs_s, inst_s))
     dsub <- data_Z[complete.cases(data_Z[, cols_need_full, drop = FALSE]), , drop = FALSE]
     if (nrow(dsub) == 0L) {
       return(data.frame(var=y_o, partial_R2=NA_real_, F_weak=NA_real_, N=0L,
                         df1=NA_integer_, df2=NA_integer_, stringsAsFactors = FALSE))
     }
     
     f_restr <- as.formula(paste(y_s, "~", mk_rhs(exogs_s)))
     f_full  <- as.formula(paste(y_s, "~", mk_rhs(c(exogs_s, inst_s))))
     
     m0 <- lm(f_restr, data = dsub)
     m1 <- lm(f_full,  data = dsub)
     
     # partial R^2 dos IVs condicional aos exogs
     R2_0 <- summary(m0)$r.squared
     R2_1 <- summary(m1)$r.squared
     pr2  <- if (is.finite(R2_0) && R2_0 < 1) (R2_1 - R2_0)/(1 - R2_0) else NA_real_
     
     # F de exclusão conjunta (ANOVA de modelos aninhados)
     a12 <- anova(m0, m1)
     Fst <- if (nrow(a12) == 2) as.numeric(a12$F[2]) else NA_real_
     df1 <- if (nrow(a12) == 2) as.integer(a12$Df[2]) else NA_integer_
     df2 <- if (nrow(a12) == 2) as.integer(a12$Res.Df[2]) else NA_integer_
     
     data.frame(var=y_o, partial_R2=pr2, F_weak=Fst, N=nrow(dsub),
                df1=df1, df2=df2, stringsAsFactors = FALSE)
   })
   
   do.call(rbind, rows)
 }
 # exemplo: força dos IVs no seu data e conjunto parcimonioso
 iv_strength <- check_iv_strength(
   data_Z     = df_iv,
   instNames  = iv_set_core,   # ou outro conjunto (mid/full)
   drop_price = 1,
   priceNames = priceNames,
   use_z2     = TRUE,
   shareNames = shareNames
 )
 print(iv_strength)
 
 library(AER)
 
 overid_per_eq_safe <- function(data_Z, shareNames, priceNames, drop_price,
                                instNames, use_z2 = TRUE) {
   stopifnot(length(drop_price) == 1)
   ln_all  <- paste0("ln_", priceNames)
   ln_drop <- paste0("ln_", priceNames[drop_price])
   ln_keep <- setdiff(ln_all, ln_drop)
   
   # exógenos (entram no RHS e podem entrar também como instrumentos)
   exogs   <- c("z", if (use_z2) "z2" else NULL)
   
   # garante z/z2 (índice Ls) se não existirem
   if (!("z" %in% names(data_Z))) {
     wbar <- colMeans(data_Z[, shareNames, drop = FALSE], na.rm = TRUE)
     wbar <- wbar / sum(wbar)
     lnP  <- as.numeric(as.matrix(data_Z[, ln_all, drop = FALSE]) %*% wbar)
     data_Z$z <- log(as.numeric(data_Z$gasto_total_atualhat)) - lnP
   }
   if (use_z2 && !("z2" %in% names(data_Z))) data_Z$z2 <- data_Z$z^2
   
   # instrumentos realmente disponíveis no data
   inst_in <- intersect(instNames, names(data_Z))
   
   mk_rhs <- function(v) if (length(v) > 0) paste(v, collapse = " + ") else "1"
   len0_to_na <- function(x) if (length(x) == 0 || is.null(x)) NA_real_ else as.numeric(x)
   
   rows <- lapply(shareNames, function(y) {
     # RHS com endógenas (ln_keep) + exógenas (exogs)
     rhs_endog <- mk_rhs(c(ln_keep, exogs))
     # conjunto de instrumentos: exógenas + IVs excluídos (inst_in)
     rhs_inst  <- mk_rhs(c(exogs, inst_in))
     
     # se não há IVs excluídos, ivreg roda como OLS (diagnostics = NULL)
     have_excluded <- length(inst_in) > 0
     
     fml <- stats::as.formula(paste(y, "~", rhs_endog, "|", rhs_inst))
     
     fit_iv <- AER::ivreg(fml, data = data_Z)
     sum_iv <- summary(fit_iv, diagnostics = TRUE)
     
     Sargan <- NA_real_
     p_Sarg <- NA_real_
     df_S   <- NA_integer_
     
     if (have_excluded && !is.null(sum_iv$diagnostics)) {
       # só tenta extrair se existir a linha "Sargan"
       if ("Sargan" %in% rownames(sum_iv$diagnostics)) {
         Sargan <- len0_to_na(unname(sum_iv$diagnostics["Sargan", "statistic"]))
         p_Sarg <- len0_to_na(unname(sum_iv$diagnostics["Sargan", "p-value"]))
         # df pode não existir em algumas versões; calcule se faltar
         if (!is.na(Sargan) && "df" %in% colnames(sum_iv$diagnostics)) {
           df_S <- as.integer(sum_iv$diagnostics["Sargan", "df"])
         } else {
           # #overid = (#inst_excluídos) - (#endógenas)
           df_S <- max(length(inst_in) - length(ln_keep), 0L)
         }
       }
     }
     
     data.frame(
       eq      = y,
       endog   = length(ln_keep),
       inst_ex = length(inst_in),
       overid  = max(length(inst_in) - length(ln_keep), 0L),
       Sargan  = Sargan,
       p_Sarg  = p_Sarg,
       df_S    = df_S,
       stringsAsFactors = FALSE
     )
   })
   
   do.call(rbind, rows)
 }
 
 
 # 1) Diagnóstico de sobre-ID com o set parcimonioso
ov_core <- overid_per_eq_safe(
  data_Z     = df_iv,
  shareNames = fit_q$shareNames,
  priceNames = fit_q$priceNames,
  drop_price = 1,                 # preco_por_kg1 omitido no RHS
  instNames  = iv_set_core,       # seu conjunto parcimonioso
  use_z2     = TRUE
)

 print(ov_core)
 
 inst_core <- intersect(iv_set_core, names(df_iv))
 inst_core <- inst_core[!is.na(inst_core) & nzchar(inst_core)]
 
 fit_q_3sls <- fit_quaids_manual_km1(
   prices     = df_iv[, priceNames, drop = FALSE],
   shares     = df_iv[, shareNames, drop = FALSE],
   x          = df_iv[["gasto_total_atualhat"]],
   priceIndex = "Ls",
   estMethod  = "3SLS",
   omit_share = 1,
   drop_price = 1,
   instNames  = inst_core,
   instData   = df_iv,
   use_z2     = TRUE,
   maxiter    = 500
 )
 summary(fit_q_3sls)
 
 setdiff(inst_core, names(df_iv))          # deve dar character(0)
 length(inst_core)                         # deve ser >= K-1 (aqui, 5)
 qr(model.matrix(~0 + ., df_iv[inst_core]))$rank  # posto de Z
 
 # F e R2 parcial (por ln_preço):
 iv_strength <- check_iv_strength(
   data_Z     = df_iv,
   instNames  = inst_core,
   drop_price = 1,
   priceNames = priceNames,
   use_z2     = TRUE,
   shareNames = shareNames
 )
 iv_strength
 
 fit_q_3sls_nz2 <- fit_quaids_manual_km1(
   prices=df_iv[, priceNames, drop=FALSE],
   shares=df_iv[, shareNames, drop=FALSE],
   x=df_iv[["gasto_total_atualhat"]],
   priceIndex="Ls", estMethod="3SLS",
   omit_share=1, drop_price=1,
   instNames=inst_core, instData=df_iv,
   use_z2=FALSE, maxiter=500
 )
 summary(fit_q_3sls_nz2)
 
 eq_list <- fit_q_3sls_nz2$fit$eq
 idx <- seq_along(eq_list)
 
 # rótulos de fallback (shares exceto o omitido)
 labs <- setdiff(fit_q_3sls_nz2$shareNames, fit_q_3sls_nz2$omit_share)
 if (length(names(eq_list)) == 0L) names(eq_list) <- labs
 
 # (1) Condicionamento (kappa) da matriz de regressão de cada equação
 kappa_eq <- vapply(idx, function(i){
   X <- model.matrix(eq_list[[i]])
   kappa(X)
 }, numeric(1))
 names(kappa_eq) <- names(eq_list)
 kappa_eq

 # (2) R^2 por equação
 r2 <- vapply(idx, function(i){
   mf <- model.frame(eq_list[[i]])
   y  <- model.response(mf)
   yhat <- fitted(eq_list[[i]])
   1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
 }, numeric(1))
 names(r2) <- names(eq_list)
 r2
 
 ln_keep <- setdiff(paste0("ln_", fit_q_3sls_nz2$priceNames),
                    paste0("ln_", fit_q_3sls_nz2$drop_price))
 Xln <- as.matrix(df_iv[, ln_keep, drop = FALSE])
 kappa(scale(Xln, center = TRUE, scale = TRUE))
 
 
 #library(AER); library(sandwich); library(lmtest)
 #f_endog <- reformulate(c(ln_keep, "z"), response = shareNames[2])   # eq. 2 de exemplo
 #f_inst  <- reformulate(c("z", inst_core))
 #m2 <- ivreg(as.formula(paste(deparse(f_endog), "|", deparse(f_inst))), data = df_iv)
 #coeftest(m2, vcov. = vcovHC(m2, type = "HC3"))
 
 eq_list <- fit_q_3sls_nz2$fit$eq
 if (is.null(names(eq_list)) || !length(names(eq_list))) {
   names(eq_list) <- setdiff(fit_q_3sls_nz2$shareNames, fit_q_3sls_nz2$omit_share)
 }
 r2_corr <- vapply(seq_along(eq_list), function(i){
   y <- model.response(model.frame(eq_list[[i]]))
   yhat <- fitted(eq_list[[i]])
   (cor(y, yhat)^2)
 }, numeric(1))
 r2_corr
 
 library(AER); library(sandwich); library(lmtest)
 
 ln_keep <- setdiff(paste0("ln_", priceNames), "ln_preco_por_kg1")
 f_endog <- reformulate(c(ln_keep, "z"), response = shareNames[2])
 f_inst  <- reformulate(c("z", inst_core))
 
 m2 <- ivreg(f_endog, instruments = f_inst, data = df_iv)
 
 # robust HC1 (evita hatvalues)
 coeftest(m2, vcov. = sandwich::vcovHC(m2, type = "HC1"))
 
 
 kappa_eq_std <- vapply(seq_along(eq_list), function(i){
   X <- model.matrix(eq_list[[i]])
   X <- scale(X[, colnames(X)!="(Intercept)"], center = TRUE, scale = TRUE)
   kappa(X)
 }, numeric(1))
 kappa_eq_std
 
 fit_aids_iv <- rebuild_quaids_coefs(fit_q_3sls_nz2)
 fit_aids_iv$coef$lambda[] <- 0  # AIDS
 E <- elas_quaids_manual(
   fit_aids_iv,
   x = mean(df_iv$gasto_total_atualhat, na.rm = TRUE),
   enforce_hicks = TRUE,
   enforce_symmetry = TRUE,
   enforce_cournot = FALSE
 )
 diag(E$hicks)    # próprias compensadas (espera-se < 0)
 E$expenditure    # elasticidades-renda (Engel)
 E$checks
 
 # === EXOG SHIFTERS (recomendado) ===
 exog_shifters <- c("f_reg", "f_area", "p_n_mais65anos", "p_sexofem")
 
 # RHS e instrumentos
 ln_keep    <- setdiff(paste0("ln_", priceNames), "ln_preco_por_kg1")
 rhs_terms  <- c(ln_keep, "z", exog_shifters)
 inst_terms <- c("z", exog_shifters, inst_core)  # exógenas também como IVs
 
 # sistema K-1 (omitindo w_despesahat1)
 keep_sh <- setdiff(shareNames, shareNames[1])
 eqs <- lapply(keep_sh, function(y) reformulate(rhs_terms, response = y))
 names(eqs) <- gsub("[^[:alnum:]]", "", keep_sh)  # rótulos sem '_' p/ systemfit
 
 library(systemfit)
 fit_3sls_exog <- systemfit(
   eqs, data = df_iv,
   method = "3SLS",
   inst   = reformulate(inst_terms, intercept = FALSE),
   maxit  = 500
 )
 summary(fit_3sls_exog)
 
 eqs <- fit_3sls_exog$eq
 kappa_eq <- sapply(eqs, function(f) kappa(scale(model.matrix(f)[,-1], TRUE, TRUE)))
 r2_corr  <- sapply(eqs, function(f){ y <- model.response(model.frame(f)); cor(y, fitted(f))^2 })
 r2_corr
 
 # (A) Garante o fator de capital/interior
 if (!"f_cap" %in% names(df_iv)) {
   stopifnot("capital_interior" %in% names(df_iv))
   df_iv$f_cap <- factor(df_iv$capital_interior)
 }
 
 # (B) Escolhe 1 demográfico extra (p_n_0esc preferido; senão p_n_8esc)
 demo_extra <- if ("p_n_0esc" %in% names(df_iv)) {
   "p_n_0esc"
 } else if ("p_n_8esc" %in% names(df_iv)) {
   "p_n_8esc"
 } else {
   stop("Nem 'p_n_0esc' nem 'p_n_8esc' encontrados no df_iv.")
 }
 
 # (C) Interação única: p_sexofem:f_reg
 # (usa ':' para só a interação, sem as marginais adicionais)
 
 ## Componentes “AIDS-IV” (sem z2)
 ln_keep   <- setdiff(paste0("ln_", priceNames), "ln_preco_por_kg1")
 
 ## RHS = ln_preços K-1 + z + shifters exógenos (inclui os 3 ajustes)
 rhs_terms <- c(
   ln_keep, "z",
   "f_reg", "f_area", "f_cap",              # regionais (inclui capital/interior)
   "p_n_mais65anos", "p_sexofem", demo_extra,  # demográficos (inclui o extra)
   "p_sexofem:f_reg"                        # interação única
 )
 
 ## Instrumentos = z + todos exógenos (incluindo interação) + núcleo forte
 inst_terms <- c(
   "z",
   "f_reg", "f_area", "f_cap",
   "p_n_mais65anos", "p_sexofem", demo_extra,
   "p_sexofem:f_reg",
   inst_core                                     # seus IVs fortes (iv_op..., IV_d_uf_...)
 )
 
 ## K-1 equações (omite w1; base de preço = preco_por_kg1)
 keep_sh <- setdiff(shareNames, shareNames[1])
 eqs <- lapply(keep_sh, function(y) reformulate(rhs_terms, response = y))
 
 ## systemfit não aceita '_' nos rótulos => higieniza nomes das equações
 names(eqs) <- gsub("[^[:alnum:]]", "", keep_sh)
 
 ## Estima 3SLS com os novos shifters e instrumentos
 library(systemfit)
 fit_3sls_exog2 <- systemfit(
   eqs, data = df_iv,
   method = "3SLS",
   inst   = reformulate(inst_terms, intercept = FALSE),
   maxit  = 500
 )
 summary(fit_3sls_exog2)
 
 ## Diagnóstico rápido (condicionamento e R^2_corr por equação)
 eqs_fit  <- fit_3sls_exog2$eq
 kappa_eq <- sapply(eqs_fit, function(m){
   X <- model.matrix(m)
   X <- scale(X[, colnames(X)!="(Intercept)"], TRUE, TRUE)
   kappa(X)
 })
 r2_corr <- setNames(sapply(eqs_fit, function(m){
   y <- model.response(model.frame(m))
   cor(y, fitted(m))^2
 }), names(eqs_fit))
 
 kappa_eq
 r2_corr
 
 ## ========================
 ## 3 ajustes nos shifters
 ## ========================
 
 # (A) Garante o fator de capital/interior
 if (!"f_cap" %in% names(df_iv)) {
   stopifnot("capital_interior" %in% names(df_iv))
   df_iv$f_cap <- factor(df_iv$capital_interior)
 }
 
 # (B) Escolhe 1 demográfico extra (p_n_0esc preferido; senão p_n_8esc)
 demo_extra <- if ("p_n_0esc" %in% names(df_iv)) {
   "p_n_0esc"
 } else if ("p_n_8esc" %in% names(df_iv)) {
   "p_n_8esc"
 } else {
   stop("Nem 'p_n_0esc' nem 'p_n_8esc' encontrados no df_iv.")
 }
 
 # (C) Centraliza p_sexofem e usa apenas a interação com f_reg por padrão
 #     (usa ':' para só a interação, sem marginais adicionais)
 df_iv$p_sexofem_c <- as.numeric(scale(df_iv$p_sexofem, center = TRUE, scale = FALSE))
 
 # Se TRUE, mantém o efeito principal de p_sexofem_c no RHS/IV; se FALSE, usa só a interação
 keep_psf_main <- FALSE
 
 ## Componentes “AIDS-IV” (sem z2)
 ln_keep <- setdiff(paste0("ln_", priceNames), "ln_preco_por_kg1")
 
 ## Base de termos exógenos sempre presentes
 base_exog <- c("f_reg", "f_area", "f_cap", "p_n_mais65anos", demo_extra)
 
 ## Interação única (com p_sexofem centrado)
 inter_term <- "p_sexofem_c:f_reg"
 
 ## RHS = ln_preços K-1 + z + exógenas (+ opcional p_sexofem_c) + interação
 rhs_terms <- c(ln_keep, "z", base_exog, inter_term)
 if (keep_psf_main) rhs_terms <- c(rhs_terms, "p_sexofem_c")
 
 ## Instrumentos = z + todas as exógenas (coerentes com o RHS) + interação + núcleo forte
 inst_terms <- c("z", base_exog, inter_term, inst_core)
 if (keep_psf_main) inst_terms <- c(inst_terms, "p_sexofem_c")
 
 ## K-1 equações (omite w1; base de preço = preco_por_kg1)
 keep_sh <- setdiff(shareNames, shareNames[1])
 eqs <- lapply(keep_sh, function(y) reformulate(rhs_terms, response = y))
 
 ## systemfit não aceita '_' nos rótulos => higieniza nomes das equações
 names(eqs) <- gsub("[^[:alnum:]]", "", keep_sh)
 
 ## Estima 3SLS com os novos shifters e instrumentos
 library(systemfit)
 fit_3sls_exog2 <- systemfit(
   eqs, data = df_iv,
   method = "3SLS",
   inst   = reformulate(inst_terms, intercept = FALSE),
   maxit  = 500
 )
 summary(fit_3sls_exog2)
 
 ## Diagnóstico rápido (condicionamento e R^2_corr por equação)
 eqs_fit  <- fit_3sls_exog2$eq
 kappa_eq <- sapply(eqs_fit, function(m){
   X <- model.matrix(m)
   X <- scale(X[, colnames(X) != "(Intercept)"], TRUE, TRUE)
   kappa(X)
 })
 r2_corr <- setNames(sapply(eqs_fit, function(m){
   y <- model.response(model.frame(m))
   cor(y, fitted(m))^2
 }), names(eqs_fit))
 
 kappa_eq
 r2_corr
 
 library(car)
 
 # endógenas de 1ª etapa (ln-preços K-1)
 ln_keep <- setdiff(paste0("ln_", priceNames), "ln_preco_por_kg1")
 
 # exógenos de controle usados no RHS/IV agora
 exog_fs <- c("z","f_reg","f_area","f_cap","p_n_mais65anos","p_n_0esc","p_sexofem_c:f_reg")
 
 # fórmula de 1ª etapa: ln_preco_k ~ exog_fs + inst_core
 form_fs <- function(y) reformulate(c(exog_fs, inst_core), response = y)
 
 # F conjunto para H0: (todos inst_core) = 0 em cada 1ª etapa
 fs_F <- sapply(ln_keep, function(y) {
   m  <- lm(form_fs(y), data = df_iv)
   lh <- linearHypothesis(m, paste(inst_core, "= 0"))
   unname(lh$F[2])
 })
 fs_F
# sargan_manual <- function(m, inst_rhs, data){
#   # 2SLS com os mesmos instrumentos
#   m_iv <- ivreg::ivreg(formula(m), instruments = inst_rhs, data = data)
#   
#   e   <- residuals(m_iv)
#   Z   <- model.matrix(inst_rhs, data)     # todos os instrumentos (inclui exógenas)
#   aux <- lm(e ~ Z - 1)                    # sem intercepto
#   n   <- nrow(Z)
#   L   <- ncol(Z)                          # nº de instrumentos
#   k   <- length(coef(m_iv))               # nº de regress. estimados (inclui intercepto)
#   J   <- n * summary(aux)$r.squared
#   df  <- L - k                            # sobre-restrições
#   p   <- 1 - pchisq(J, df)
#   c(J_homo = unname(J), df_homo = unname(df), p_Sargan = unname(p))
# }
 

# Sargan_manual_df <- data.frame(
#   eq = vapply(eqs_fit, function(m) as.character(formula(m)[[2]]), ""),
#   t(sapply(eqs_fit, sargan_manual, inst_rhs = inst_rhs, data = df_iv)),
#   row.names = NULL
# )
# Sargan_manual_df
 
 # fitted para w2..w6:
 yhat_mat <- sapply(fit_3sls_exog2$eq, fitted)
 w1_hat   <- 1 - rowSums(yhat_mat)
 
 summary(w1_hat)                       # distribuição do w1 previsto
 range(yhat_mat)                       # checar se há share <0 ou >1
 max_abs_addup <- max(abs(w1_hat + rowSums(yhat_mat) - 1))
 max_abs_addup
 #--------------------------------AQUI
 # ==========================================================
 # QUAIDS (AIDS-IV) com 3SLS e shifters exógenos parsimoniosos
 # - cria f_cap (se preciso)
 # - escolhe demográfico extra (p_n_0esc preferido; senão p_n_8esc)
 # - centraliza p_sexofem -> p_sexofem_c
 # - usa interação única p_sexofem_c:f_reg
 # - permite z2 opcional
 # - recria share omitido (observado e previsto)
 # - devolve diagnósticos úteis
 # ==========================================================
 
 # ==========================================================
 # QUAIDS (AIDS-IV) com 3SLS e shifters exógenos parsimoniosos
 # Versão enxuta (sem testes): sem F 1ª etapa e sem Sargan/Hansen
 # - cria f_cap (se preciso)
 # - escolhe demográfico extra (p_n_0esc preferido; senão p_n_8esc)
 # - centraliza p_sexofem -> p_sexofem_c
 # - usa interação única p_sexofem_c:f_reg
 # - opcionalmente inclui z2
 # - recria share omitido (observado e previsto)
 # - diagnósticos: kappa, R2_corr e adding-up
 # ==========================================================
 
 fit_quaids_manual_km1 <- function(
    prices,                      # data.frame/matrix com preços (colunas = priceNames)
    shares,                      # data.frame/matrix com shares (colunas = shareNames)
    x,                           # gasto total (mantido por compatibilidade)
    priceIndex = "Ls",
    estMethod  = "3SLS",
    omit_share = 1,              # índice do share omitido (w1 por padrão)
    drop_price = 1,              # índice do preço/base omitido no RHS (p1 por padrão)
    instNames,                   # vetor com IVs "core"
    instData,                    # data.frame com todas as variáveis (df_iv)
    use_z2     = FALSE,
    maxiter    = 500,
    # opções:
    include_cap     = TRUE,
    demo_pref       = c("p_n_0esc","p_n_8esc"),
    center_sex      = TRUE,
    add_interaction = TRUE,
    compute_diag    = TRUE       # apenas kappa / R2_corr / adding-up
 ) {
   
   if (!requireNamespace("systemfit", quietly = TRUE)) {
     stop("Pacote 'systemfit' é necessário.")
   }
   
   # ---- Nomes e checagens básicas ----
   priceNames <- colnames(as.data.frame(prices))
   shareNames <- colnames(as.data.frame(shares))
   if (is.null(priceNames) || is.null(shareNames)) {
     stop("As matrizes 'prices' e 'shares' precisam ter nomes de colunas.")
   }
   if (nrow(prices) != nrow(shares) || nrow(prices) != nrow(instData)) {
     stop("As quantidades de linhas de 'prices', 'shares' e 'instData' devem coincidir.")
   }
   if (!all(instNames %in% names(instData))) {
     faltam <- setdiff(instNames, names(instData))
     stop("IVs ausentes em instData: ", paste(faltam, collapse = ", "))
   }
   
   df <- instData
   
   # ---- (A) f_cap (se faltar) ----
   if (include_cap && !"f_cap" %in% names(df)) {
     if (!"capital_interior" %in% names(df)) {
       stop("Nem 'f_cap' nem 'capital_interior' encontrados em instData.")
     }
     df$f_cap <- factor(df$capital_interior)
   }
   
   # Garante fatores existentes
   if (!"f_reg"  %in% names(df)) stop("'f_reg' ausente em instData.")
   if (!"f_area" %in% names(df)) stop("'f_area' ausente em instData.")
   if (!is.factor(df$f_reg))  df$f_reg  <- factor(df$f_reg)
   if (!is.factor(df$f_area)) df$f_area <- factor(df$f_area)
   if (include_cap && !is.factor(df$f_cap)) df$f_cap <- factor(df$f_cap)
   
   # ---- (B) demográfico extra ----
   demo_pref  <- match.arg(demo_pref)
   demo_extra <- if (demo_pref %in% names(df)) {
     demo_pref
   } else if (identical(demo_pref, "p_n_0esc") && "p_n_8esc" %in% names(df)) {
     "p_n_8esc"
   } else if (identical(demo_pref, "p_n_8esc") && "p_n_0esc" %in% names(df)) {
     "p_n_0esc"
   } else {
     stop("Nem 'p_n_0esc' nem 'p_n_8esc' encontrados em instData.")
   }
   
   # ---- (C) centraliza p_sexofem e cria interação única com f_reg ----
   if (!"p_sexofem" %in% names(df)) stop("'p_sexofem' ausente em instData.")
   sex_var <- "p_sexofem"
   if (isTRUE(center_sex)) {
     df$p_sexofem_c <- as.numeric(scale(df$p_sexofem, center = TRUE, scale = FALSE))
     sex_var <- "p_sexofem_c"
   }
   inter_term <- if (isTRUE(add_interaction)) paste0(sex_var, ":f_reg") else NULL
   
   # ---- ln-preços no RHS (K-1) ----
   ln_keep <- setdiff(paste0("ln_", priceNames), paste0("ln_", priceNames[drop_price]))
   if (!all(ln_keep %in% names(df))) {
     faltam <- setdiff(ln_keep, names(df))
     stop("Variáveis de ln-preço ausentes em instData: ", paste(faltam, collapse = ", "))
   }
   
   # ---- RHS e instrumentos (limpando NAs/vazios) ----
   base_exog <- c("z", "f_reg", "f_area", if (include_cap) "f_cap",
                  "p_n_mais65anos", sex_var, demo_extra)
   rhs_terms  <- c(ln_keep, base_exog, inter_term)
   inst_terms <- c("z", "f_reg", "f_area", if (include_cap) "f_cap",
                   "p_n_mais65anos", sex_var, demo_extra, inter_term, instNames)
   if (isTRUE(use_z2) && "z2" %in% names(df)) {
     rhs_terms  <- c(rhs_terms,  "z2")
     inst_terms <- c(inst_terms, "z2")
   }
   rhs_terms  <- unique(rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)])
   inst_terms <- unique(inst_terms[!is.na(inst_terms) & nzchar(inst_terms)])
   
   # ---- Equações K-1 (omite share omitido) ----
   keep_sh <- setdiff(shareNames, shareNames[omit_share])
   if (!all(keep_sh %in% names(df))) {
     faltam <- setdiff(keep_sh, names(df))
     stop("Shares ausentes em instData: ", paste(faltam, collapse = ", "))
   }
   eqs <- lapply(keep_sh, function(y) stats::reformulate(rhs_terms, response = y))
   names(eqs) <- gsub("[^[:alnum:]]", "", keep_sh)  # systemfit não aceita '_' nos rótulos
   
   # ---- Estimação ----
   if (!identical(estMethod, "3SLS")) {
     warning("Estimation method alterado para '3SLS' (único implementado aqui).")
   }
   fit <- systemfit::systemfit(
     eqs, data = df,
     method = "3SLS",
     inst   = stats::reformulate(inst_terms, intercept = FALSE),
     maxit  = maxiter
   )
   
   # ---- Recria share omitido (observado e previsto) ----
   shares_df <- as.data.frame(shares)
   # Observado via adding-up das K-1 observadas
   obs_included <- as.matrix(shares_df[, keep_sh, drop = FALSE])
   w_omit_obs_from_add <- 1 - rowSums(obs_included)
   
   # Previsto via soma dos fitted das K-1
   yhat_mat <- sapply(fit$eq, stats::fitted)
   if (is.list(yhat_mat)) yhat_mat <- do.call(cbind, yhat_mat)
   w_omit_hat <- 1 - rowSums(yhat_mat)
   
   # ---- Diagnósticos leves (opcional) ----
   diag <- NULL
   if (isTRUE(compute_diag)) {
     kappa_eq <- sapply(fit$eq, function(m){
       X <- stats::model.matrix(m)
       X <- scale(X[, colnames(X) != "(Intercept)"], TRUE, TRUE)
       kappa(X)
     })
     r2_corr <- sapply(fit$eq, function(m){
       y <- stats::model.response(stats::model.frame(m))
       stats::cor(y, stats::fitted(m))^2
     })
     add_up_max_abs <- max(abs(w_omit_hat + rowSums(yhat_mat) - 1))
     names(kappa_eq) <- names(fit$eq)
     names(r2_corr)  <- names(fit$eq)
     diag <- list(
       kappa_eq       = kappa_eq,
       r2_corr        = r2_corr,
       add_up_max_abs = add_up_max_abs
     )
   }
   
   # ---- Monta fitted/observed completos (incluindo a omitida) ----
   fitted_included_df <- as.data.frame(yhat_mat, optional = TRUE)
   colnames(fitted_included_df) <- keep_sh
   fitted_full <- fitted_included_df
   fitted_full[[shareNames[omit_share]]] <- w_omit_hat
   fitted_full <- fitted_full[, shareNames, drop = FALSE]
   
   observed_included_df <- as.data.frame(obs_included, optional = TRUE)
   colnames(observed_included_df) <- keep_sh
   observed_full <- observed_included_df
   observed_full[[shareNames[omit_share]]] <- w_omit_obs_from_add
   observed_full <- observed_full[, shareNames, drop = FALSE]
   
   # ---- Retorno ----
   res <- list(
     call        = match.call(),
     method      = "3SLS",
     priceIndex  = priceIndex,
     omit_share  = shareNames[omit_share],
     drop_price  = priceNames[drop_price],
     priceNames  = priceNames,
     shareNames  = shareNames,
     rhs_terms   = rhs_terms,
     inst_terms  = inst_terms,
     fit         = fit,
     eq          = fit$eq,
     fitted_shares   = fitted_full,
     observed_shares = observed_full,
     omitted_recreated = list(
       name     = shareNames[omit_share],
       observed = w_omit_obs_from_add,
       fitted   = w_omit_hat
     ),
     diagnostics = diag
   )
   class(res) <- "quaids_km1_fit"
   return(res)
 }
 
 # -------- summary() enxuto --------
 summary.quaids_km1_fit <- function(object, ...) {
   cat("Demand analysis with the QUADRATIC Almost Ideal Demand System (QUAIDS)\n")
   cat("Estimation Method:", object$method, "\n")
   cat("Price Index:", object$priceIndex, "\n")
   cat("Omitted share equation:", object$omit_share, "\n")
   cat("Dropped price from RHS:", object$drop_price, "\n\n")
   print(summary(object$fit))
   if (!is.null(object$diagnostics)) {
     cat("\nQuick diagnostics (no tests)\n")
     with(object$diagnostics, {
       cat("  kappa (std. X) by equation:\n"); print(kappa_eq)
       cat("  R2 (corr^2) by equation:\n");   print(r2_corr)
       cat("  Adding-up |max abs deviation| (fitted):",
           format(add_up_max_abs, digits = 6), "\n")
     })
   }
   invisible(object)
 }
 
 # -------- summary() para o objeto retornado --------
 summary.quaids_km1_fit <- function(object, ...) {
   cat("Demand analysis with the QUADRATIC Almost Ideal Demand System (QUAIDS)\n")
   cat("Estimation Method:", object$method, "\n")
   cat("Price Index:", object$priceIndex, "\n")
   cat("Omitted share equation:", object$omit_share, "\n")
   cat("Dropped price from RHS:", object$drop_price, "\n\n")
   print(summary(object$fit))
   if (!is.null(object$diagnostics)) {
     cat("\nQuick diagnostics\n")
     diag <- object$diagnostics
     if (!is.null(diag$kappa_eq)) {
       cat("  kappa (std. X) by equation:\n"); print(diag$kappa_eq)
     }
     if (!is.null(diag$r2_corr)) {
       cat("  R2 (corr^2) by equation:\n"); print(diag$r2_corr)
     }
     if (!is.null(diag$fs_F)) {
       cat("  First-stage joint F (IV core) by ln-price:\n"); print(diag$fs_F)
     }
     if (!is.null(diag$add_up_max_abs)) {
       cat("  Adding-up |max abs deviation| (fitted):", format(diag$add_up_max_abs, digits = 5), "\n")
     }
     if (!is.null(diag$overid)) {
       cat("\n  Overidentification (Sargan/Hansen J):\n")
       print(diag$overid, row.names = FALSE)
     }
   }
   invisible(object)
 }
 
 fit_q <- fit_quaids_manual_km1(
   prices     = df_iv[, priceNames, drop = FALSE],
   shares     = df_iv[, shareNames, drop = FALSE],
   x          = df_iv[["gasto_total_atualhat"]],
   priceIndex = "Ls",
   estMethod  = "3SLS",
   omit_share = 1,
   drop_price = 1,
   instNames  = inst_core,
   instData   = df_iv,
   use_z2     = TRUE,    # ou TRUE se quiser usar z2
   maxiter    = 500
 )
 
 summary(fit_q)
 quaids_elasticities <- function(fit,
                                 at = c("means","medians"),
                                 w_source = c("fitted","observed"),
                                 point = NULL) {
   at       <- match.arg(at)
   w_source <- match.arg(w_source)
   
   # --- localizar eqs/nomes/omissões ---
   if (!is.list(fit) || is.null(fit$eq) || is.null(fit$priceNames) || is.null(fit$shareNames))
     stop("Objeto 'fit' precisa vir de fit_quaids_manual_km1() (com $eq, $priceNames, $shareNames).")
   
   priceNames <- fit$priceNames
   shareNames <- fit$shareNames
   K <- length(priceNames)
   
   # índices de omitidos (fornecidos como nomes no seu fit)
   j0 <- match(fit$drop_price,  priceNames)
   i0 <- match(fit$omit_share,  shareNames)
   if (is.na(j0) || is.na(i0)) stop("drop_price/omit_share não reconhecidos dentro de 'fit'.")
   
   eqs <- fit$eq
   if (length(eqs) != K-1) stop("Esperava K-1 equações em fit$eq.")
   
   # --- extrair coeficientes e reconstruir omitidos por adding-up ---
   beta   <- rep(NA_real_, K)           # coef. de z
   lambda <- rep(0, K)                  # coef. de z2 (se existir)
   gamma  <- matrix(0, K, K, dimnames=list(shareNames, priceNames))
   
   ln_rhs <- paste0("ln_", priceNames[-j0])  # ln-preços presentes no RHS
   
   obs_rows <- setdiff(seq_len(K), i0)
   for (ii in seq_along(obs_rows)) {
     i <- obs_rows[ii]
     ci <- coef(eqs[[ii]])
     if (!("z" %in% names(ci))) stop("Coeficiente 'z' ausente na eq. de ", shareNames[i])
     beta[i] <- ci[["z"]]
     if ("z2" %in% names(ci)) lambda[i] <- ci[["z2"]]
     for (nm in ln_rhs) if (nm %in% names(ci)) {
       pj <- sub("^ln_", "", nm)
       gamma[i, pj] <- ci[[nm]]
     }
   }
   # completar coluna de preço omitido e linha de share omitido
   gamma[, j0] <- -rowSums(gamma[, -j0, drop=FALSE])
   gamma[i0, ] <- -colSums(gamma[obs_rows, , drop=FALSE])
   
   # completar beta/lambda para o share omitido (adding-up)
   beta[i0]   <- -sum(beta[obs_rows])
   lambda[i0] <- -sum(lambda[obs_rows])
   
   # --- z (ponto de avaliação) ---
   # se 'point' não foi passado, usamos z do model.frame (não precisamos dos preços brutos)
   z_eval <- {
     zvec <- tryCatch(model.frame(eqs[[1]])[["z"]], error=function(e) NULL)
     if (is.null(zvec)) stop("Não encontrei 'z' no model.frame das equações.")
     if (at == "means") mean(zvec, na.rm=TRUE) else stats::median(zvec, na.rm=TRUE)
   }
   
   # Se o usuário fornecer um ponto custom (prices,x), sobrescrevemos z_eval a partir dele.
   if (!is.null(point)) {
     if (!is.list(point) || is.null(point$prices) || is.null(point$x))
       stop("`point` deve ser list(prices = ..., x = ...).")
     p_vec <- as.numeric(if (!is.null(names(point$prices))) point$prices[priceNames] else point$prices)
     if (length(p_vec) != K) stop("`point$prices` deve ter ", K, " elementos (na ordem de priceNames).")
     
     # pesos para Stone no ponto: usamos shares médios como aproximação
     w_base <- if (w_source == "fitted") colMeans(fit$fitted_shares,  na.rm=TRUE)
     else                      colMeans(fit$observed_shares, na.rm=TRUE)
     ln_p   <- log(p_vec)
     P_st   <- exp(sum(w_base * ln_p))
     z_eval <- log(as.numeric(point$x)[1] / P_st)
   }
   z2_eval <- z_eval^2
   
   # --- shares no ponto (w_eval) ---
   # por padrão, usamos a média colunar dos shares já presentes no objeto
   w_eval <- if (w_source == "fitted") colMeans(fit$fitted_shares,  na.rm=TRUE)
   else                      colMeans(fit$observed_shares, na.rm=TRUE)
   if (any(!is.finite(w_eval)) || length(w_eval) != K) stop("Falha ao obter w_eval.")
   if (any(w_eval <= 0)) stop("Algum w_eval <= 0; não é possível calcular elasticidades (divisão por w).")
   
   # --- elasticidades: η_i, Marshallianas e Hicksianas ---
   num_expend <- beta + 2*lambda*z_eval                     # β_i + 2 λ_i z
   eta <- 1 + num_expend / w_eval                           # η_i = 1 + (...) / w_i
   
   epsM <- matrix(NA_real_, K, K, dimnames=list(shareNames, priceNames))
   for (i in seq_len(K)) {
     for (j in seq_len(K)) {
       epsM[i, j] <- (-as.numeric(i == j)) + (gamma[i, j] - num_expend[i]*w_eval[j]) / w_eval[i]
     }
   }
   epsH <- epsM + tcrossprod(eta, w_eval)                   # ε^H = ε^M + η ⊗ w
   
   list(
     expenditure = eta,
     marshallian = epsM,
     hicksian    = epsH,
     w_eval      = w_eval,
     point       = list(z = z_eval, z2 = z2_eval,
                        prices = if (!is.null(point)) as.numeric(point$prices) else NULL,
                        x = if (!is.null(point)) as.numeric(point$x) else NULL)
   )
 }
 
 
 elas_mean <- quaids_elasticities(fit_q, at = "means", w_source = "fitted")
 elas_mean$expenditure
 elas_mean$marshallian     # ε^M_ij
 elas_mean$hicksian        # ε^H_ij
 
 elas_med_obs  <- quaids_elasticities(fit_q, at = "medians", w_source = "observed")
 elas_med_obs$expenditure
 elas_med_obs$marshallian     # ε^M_ij
 elas_med_obs$hicksian  
 
 px_med <- vapply(df_iv[, fit_q$priceNames, drop=FALSE], median, numeric(1), na.rm=TRUE)
 x_p25  <- as.numeric(quantile(df_iv$gasto_total_atualhat, 0.25, na.rm=TRUE))
 elas_poor <- quaids_elasticities(fit_q, w_source = "fitted",
                                  point = list(prices = px_med, x = x_p25))
 
 elas_poor$expenditure
 elas_poor$marshallian     # ε^M_ij
 elas_poor$hicksian
 
 # ------------------------------------------------------------
 # IC e significância das elasticidades QUAIDS por simulação
 # Requer: MASS (mvrnorm)
 # Depende: quaids_elasticities(fit, at=..., w_source=..., point=...)
 #          e objeto 'fit' vindo de fit_quaids_manual_km1(...)
 # ------------------------------------------------------------
 # ------------------------------------------------------------
 # IC e significância das elasticidades QUAIDS por simulação
 # (agora preservando nomes dos coeficientes ao substituir)
 # ------------------------------------------------------------
 quaids_elasticities_ci <- function(fit,
                                    at = c("medians","means")[1],
                                    w_source = c("observed","fitted")[1],
                                    point = NULL,
                                    R = 2000,
                                    level = 0.95,
                                    seed = NULL) {
   stopifnot(inherits(fit, "quaids_km1_fit"))
   if (!requireNamespace("MASS", quietly = TRUE)) {
     stop("Preciso do pacote MASS para mvrnorm().")
   }
   if (!is.null(seed)) set.seed(seed)
   
   el_hat <- quaids_elasticities(fit, at = at, w_source = w_source, point = point)
   
   .vec_elas <- function(el) {
     c(
       exp = unclass(el$expenditure),
       as.vector(t(el$marshallian)),
       as.vector(t(el$hicksian))
     )
   }
   .unvec_elas <- function(v, fit) {
     K <- length(fit$shareNames)
     nm_g <- fit$shareNames; nm_p <- fit$priceNames
     i1 <- 1:K
     i2 <- (K + 1):(K + K*K)
     i3 <- (K + K*K + 1):(K + 2*K*K)
     list(
       expenditure = stats::setNames(v[i1], nm_g),
       marshallian = `dimnames<-`(matrix(v[i2], nrow=K, byrow=TRUE), list(nm_g, nm_p)),
       hicksian    = `dimnames<-`(matrix(v[i3], nrow=K, byrow=TRUE), list(nm_g, nm_p))
     )
   }
   
   v_hat <- .vec_elas(el_hat)
   
   b_hat <- stats::coef(fit$fit)
   Vb    <- stats::vcov(fit$fit)
   
   # >>>>>>> ALTERADO: preserva nomes ao substituir coeficientes <<<<<<<
   .replace_coefs <- function(fit_obj, b_new) {
     fit2 <- fit_obj
     
     # vetor empilhado do systemfit
     orig_names_vec <- names(fit2$fit$coefficients)
     if (is.null(orig_names_vec)) orig_names_vec <- rep.int(NA_character_, length(fit2$fit$coefficients))
     # substitui mantendo nomes/ordem
     fit2$fit$coefficients[] <- b_new
     names(fit2$fit$coefficients) <- orig_names_vec
     
     # por equação (em fit$fit$eq e no atalho fit$eq)
     off <- 0L
     for (j in seq_along(fit2$fit$eq)) {
       pj <- length(fit2$fit$eq[[j]]$coefficients)
       idx <- (off + 1):(off + pj)
       
       # em fit$fit$eq
       orig_names_eq <- names(fit2$fit$eq[[j]]$coefficients)
       fit2$fit$eq[[j]]$coefficients[] <- b_new[idx]
       names(fit2$fit$eq[[j]]$coefficients) <- orig_names_eq
       
       # em fit$eq (espelho usado por quaids_elasticities)
       if (!is.null(fit2$eq) && length(fit2$eq) == length(fit2$fit$eq)) {
         orig_names_eq2 <- names(fit2$eq[[j]]$coefficients)
         fit2$eq[[j]]$coefficients[] <- b_new[idx]
         names(fit2$eq[[j]]$coefficients) <- orig_names_eq2
       }
       
       off <- off + pj
     }
     fit2
   }
   # >>>>>>> FIM da alteração <<<<<<<
   
   sim_mat <- matrix(NA_real_, nrow = R, ncol = length(v_hat))
   for (r in seq_len(R)) {
     b_draw <- as.numeric(MASS::mvrnorm(1L, mu = b_hat, Sigma = Vb))
     fit_r  <- .replace_coefs(fit, b_draw)
     el_r   <- quaids_elasticities(fit_r, at = at, w_source = w_source, point = point)
     sim_mat[r, ] <- .vec_elas(el_r)
   }
   
   alpha <- 1 - level
   se    <- apply(sim_mat, 2, stats::sd, na.rm = TRUE)
   zval  <- (v_hat - 0) / se
   pval  <- 2 * (1 - stats::pnorm(abs(zval)))
   ql    <- apply(sim_mat, 2, stats::quantile, probs = alpha/2, na.rm = TRUE)
   qu    <- apply(sim_mat, 2, stats::quantile, probs = 1 - alpha/2, na.rm = TRUE)
   
   K <- length(fit$shareNames)
   nm_g <- fit$shareNames; nm_p <- fit$priceNames
   
   pack <- function(est_vec, se_vec, z_vec, p_vec, l_vec, u_vec) {
     exp_df <- data.frame(
       good = nm_g,
       est  = est_vec[1:K],
       se   = se_vec[1:K],
       z    = z_vec[1:K],
       p    = p_vec[1:K],
       lwr  = l_vec[1:K],
       upr  = u_vec[1:K],
       row.names = NULL
     )
     
     m_idx <- (K + 1):(K + K*K)
     M_est <- matrix(est_vec[m_idx], nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     M_se  <- matrix(se_vec[m_idx],  nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     M_z   <- matrix(z_vec[m_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     M_p   <- matrix(p_vec[m_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     M_l   <- matrix(l_vec[m_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     M_u   <- matrix(u_vec[m_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     
     h_idx <- (K + K*K + 1):(K + 2*K*K)
     H_est <- matrix(est_vec[h_idx], nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     H_se  <- matrix(se_vec[h_idx],  nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     H_z   <- matrix(z_vec[h_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     H_p   <- matrix(p_vec[h_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     H_l   <- matrix(l_vec[h_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     H_u   <- matrix(u_vec[h_idx],   nrow=K, byrow=TRUE, dimnames=list(nm_g, nm_p))
     
     list(
       expenditure = exp_df,
       marshallian = list(est = M_est, se = M_se, z = M_z, p = M_p, lwr = M_l, upr = M_u),
       hicksian    = list(est = H_est, se = H_se, z = H_z, p = H_p, lwr = H_l, upr = H_u),
       meta = list(level = level, R = R, at = at, w_source = w_source,
                   point = point, seed = seed)
     )
   }
   
   pack(.vec_elas(el_hat), se, zval, pval, ql, qu)
 }
 
 res_ci <- quaids_elasticities_ci(fit_q,
                                  at = "medians",
                                  w_source = "observed",
                                  R = 2000, level = 0.95, seed = 123)
 
 res_mean_ci <- quaids_elasticities_ci(fit_q,
                                  at = "medians",
                                  w_source = "observed",
                                  R = 2000, level = 0.95, seed = 123)
 
 px_med <- vapply(df_iv[, fit_q$priceNames, drop=FALSE], median, numeric(1), na.rm=TRUE)
 x_p25  <- as.numeric(quantile(df_iv$gasto_total_atualhat, 0.25, na.rm=TRUE))

 res_px_ci <- quaids_elasticities_ci(fit_q,
                                  at = "medians",
                                  w_source = "fitted",
                                  point = list(prices = px_med, x = x_p25),
                                  R = 2000, level = 0.95, seed = 123)
 
 # Expenditure (η_i) com EP, z, p e IC
 res_ci$expenditure
 
 # Marshallianas: matrizes com est, se, z, p, lwr, upr
 res_ci$marshallian$est["w_despesahat1", ]   # linha do bem 1
 res_ci$marshallian$se["w_despesahat1", ]
 res_ci$marshallian$lwr["w_despesahat1", ]
 res_ci$marshallian$upr["w_despesahat1", ]
 
 # Hicksianas idem
 res_ci$marshallian
 
 res_ci$hicksian
 
 res_mean <- quaids_elasticities_ci(fit_q,
                                      at = "means",
                                      w_source = "observed",
                                      R = 2000, level = 0.95, seed = 123)
 
 ## ================= Over-ID por equação (robusto e resiliente) =================
 suppressPackageStartupMessages({
   if (!requireNamespace("MASS", quietly = TRUE)) stop("Instale 'MASS'")
 })
 
 # Constrói matriz de instrumentos a partir de uma lista de termos, filtrando problemas
 .build_Z <- function(terms_chr, data, verbose = TRUE) {
   keep <- list(); dropped <- character(0)
   for (trm in terms_chr) {
     mm <- try(model.matrix(reformulate(trm, intercept = FALSE), data = data), silent = TRUE)
     if (inherits(mm, "try-error") || is.null(dim(mm)) || ncol(mm) == 0L) {
       dropped <- c(dropped, trm); next
     }
     # remove colunas com var ~ 0 ou NA
     ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
     if (!all(ok)) {
       if (verbose) message("Drop colunas zero-var de '", trm, "': ",
                            paste(colnames(mm)[!ok], collapse = ", "))
       mm <- mm[, ok, drop = FALSE]
       if (ncol(mm) == 0L) { dropped <- c(dropped, trm); next }
     }
     keep[[trm]] <- mm
   }
   if (!length(keep)) stop("Nenhum instrumento válido encontrado em 'inst_terms'.")
   Z <- do.call(cbind, keep)
   
   # remove colinearidade em Z via QR
   qz <- qr(Z); r <- qz$rank
   if (r < ncol(Z)) {
     piv <- qz$pivot[seq_len(r)]
     if (verbose) {
       gone <- setdiff(colnames(Z), colnames(Z)[piv])
       message("Drop colinear em Z: ", paste(gone, collapse = ", "))
     }
     Z <- Z[, piv, drop = FALSE]
   }
   Z
 }
 
 # Projeção por QR: P_Z * A
 .PZ <- function(Z, A) qr.fitted(qr(Z), A)
 
 # inversão segura
 .solve_safe <- function(M, b) {
   out <- try(solve(M, b), silent = TRUE)
   if (inherits(out, "try-error") || any(!is.finite(out))) MASS::ginv(M) %*% b else out
 }
 
 # J com matriz de ponderação S (genérica); com regularização
 .J_stat <- function(gbar, S, n, ridge_scale = 1e-8, eig_tol = 1e-10) {
   S <- (S + t(S)) / 2
   dmean <- mean(diag(S)); if (!is.finite(dmean) || dmean <= 0) dmean <- 1
   S_r <- S + diag(ridge_scale * dmean, ncol(S))
   val <- try(as.numeric(n * t(gbar) %*% solve(S_r, gbar)), silent = TRUE)
   if (!inherits(val, "try-error") && is.finite(val)) return(val)
   ee <- eigen(S, symmetric = TRUE)
   pos <- ee$values > (max(ee$values, na.rm = TRUE) * eig_tol)
   if (any(pos)) {
     S_psd <- ee$vectors[, pos, drop = FALSE] %*%
       diag(ee$values[pos], sum(pos)) %*%
       t(ee$vectors[, pos, drop = FALSE])
     S_psd <- S_psd + diag(ridge_scale * dmean, ncol(S_psd))
     val2 <- try(as.numeric(n * t(gbar) %*% solve(S_psd, gbar)), silent = TRUE)
     if (!inherits(val2, "try-error") && is.finite(val2)) return(val2)
   }
   as.numeric(n * t(gbar) %*% MASS::ginv(S) %*% gbar)
 }
 
 # --- Kit mínimo de diagnósticos por equação (sem cluster) --------------------
 # Usa AER::ivreg para evitar conflitos com o pacote 'ivreg'
 
 if (!requireNamespace("AER", quietly = TRUE)) stop("Instale 'AER'.")
 if (!requireNamespace("sandwich", quietly = TRUE)) stop("Instale 'sandwich'.")
 
 # =====================================================================
 # Diagnóstico IV por equação (sem cluster), robusto a falhas pontuais
 # Requer: AER, sandwich
 # =====================================================================
 # ---------------------------------------------------------------
 # Diagnóstico IV por equação (sem cluster) — robusto a nomes faltantes
 # Requer: AER, sandwich
 # ---------------------------------------------------------------
 diag_iv_basico <- function(fit_obj, data, inst_terms, show_msgs = TRUE){
   if (!inherits(fit_obj, "quaids_km1_fit")) stop("fit_obj não é 'quaids_km1_fit'.")
   if (!requireNamespace("AER", quietly=TRUE)) stop("Pacote 'AER' não encontrado.")
   if (!requireNamespace("sandwich", quietly=TRUE)) stop("Pacote 'sandwich' não encontrado.")
   
   eqs <- fit_obj$eq
   if (is.null(eqs) || !length(eqs)) stop("fit_obj$eq vazio.")
   data <- as.data.frame(data)
   inst_rhs <- reformulate(inst_terms)  # ~ z + ...
   
   # helper: extrai valor de tabela de diag; retorna NA_real_ se não existir
   get_val <- function(x, row, col){
     if (is.null(x)) return(NA_real_)
     rn <- rownames(x)
     if (is.null(rn) || !row %in% rn) return(NA_real_)
     as.numeric(x[row, col])
   }
   
   linhas <- lapply(seq_along(eqs), function(i){
     f_yx <- formula(eqs[[i]])
     # nome seguro da equação
     nm <- tryCatch({
       lhs <- try(as.character(f_yx[[2]]), silent=TRUE)
       if (!inherits(lhs, "try-error") && length(lhs)==1 && nzchar(lhs)) lhs else paste0("eq", i)
     }, error=function(e) paste0("eq", i))
     
     out_try <- try({
       m_iv <- AER::ivreg(f_yx, instruments = inst_rhs, data = data)
       
       # summary com diagnósticos (pode faltar)
       sm  <- try(summary(m_iv, diagnostics = TRUE), silent = TRUE)
       diag_tab <- if (!inherits(sm, "try-error")) sm$diagnostics else NULL
       
       # Sargan homoc. e Hansen-J robusto HC3 (podem falhar se eq. exatamente identificada)
       s_h <- try(AER::sargan(m_iv), silent = TRUE)
       s_r <- try(AER::sargan(m_iv, vcov = sandwich::vcovHC(m_iv, type = "HC3")), silent = TRUE)
       
       # monta uma linha, garantindo comprimento 1 em cada coluna
       data.frame(
         eq        = nm,
         F_weak    = get_val(diag_tab, "Weak instruments", "statistic"),
         p_weak    = get_val(diag_tab, "Weak instruments", "p-value"),
         p_wu      = get_val(diag_tab, "Wu-Hausman",      "p-value"),
         J_sargan  = if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$statistic)) else NA_real_,
         df_J      = if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$parameter)) else NA_real_,
         p_sargan  = if (!inherits(s_h, "try-error")) unname(as.numeric(s_h$p.value))   else NA_real_,
         J_hansen  = if (!inherits(s_r, "try-error")) unname(as.numeric(s_r$statistic)) else NA_real_,
         df_Jr     = if (!inherits(s_r, "try-error")) unname(as.numeric(s_r$parameter)) else NA_real_,
         p_hansen  = if (!inherits(s_r, "try-error")) unname(as.numeric(s_r$p.value))   else NA_real_,
         stringsAsFactors = FALSE
       )
     }, silent = TRUE)
     
     if (inherits(out_try, "try-error")){
       if (isTRUE(show_msgs)) message("Falhou em eq '", nm, "': ", attr(out_try, "condition")$message)
       return(NULL)
     }
     out_try
   })
   
   linhas <- Filter(Negate(is.null), linhas)
   
   if (!length(linhas)) {
     # devolve data.frame vazio com colunas definidas (nunca NULL)
     return(data.frame(
       eq=character(0), F_weak=numeric(0), p_weak=numeric(0), p_wu=numeric(0),
       J_sargan=numeric(0), df_J=numeric(0), p_sargan=numeric(0),
       J_hansen=numeric(0), df_Jr=numeric(0), p_hansen=numeric(0)
     ))
   }
   out <- do.call(rbind, linhas)
   rownames(out) <- NULL
   out
 }
 
 options(error = NULL)  # só para garantir que não há handler ativo de debug
 diag_tab <- diag_iv_basico(fit_q, df_iv, inst_terms = fit_q$inst_terms)
 print(diag_tab)
 
 # ==== First-stage joint F por endógeno (IVs excluídos) ====
 fs_block_F <- function(fit_obj, data, robust = TRUE){
   stopifnot(inherits(fit_obj, "quaids_km1_fit"))
   if (!requireNamespace("car", quietly=TRUE)) stop("Precisa do pacote 'car'.")
   if (!requireNamespace("sandwich", quietly=TRUE)) stop("Precisa do pacote 'sandwich'.")
   
   data <- as.data.frame(data)
   
   # termos no RHS das equações (inclui controles + ln_preços)
   rhs_all <- fit_obj$rhs_terms
   # instrumentos passados ao 3SLS
   inst_all <- fit_obj$inst_terms
   
   # endógenos = log-preços no RHS
   endo_ln <- grep("^ln_", rhs_all, value = TRUE)
   
   # controles = o que aparece em RHS e também está em 'inst_terms' (são instrumentos "incluídos")
   controls <- intersect(inst_all, rhs_all)
   # IVs excluídos = instrumentos que NÃO estão no RHS
   iv_excl  <- setdiff(inst_all, controls)
   
   if (length(iv_excl) == 0L) stop("Não há IVs excluídos detectados.")
   
   form_fs <- function(y) reformulate(c(controls, iv_excl), response = y)
   
   get_one <- function(y){
     m  <- stats::lm(form_fs(y), data = data)
     K  <- paste(iv_excl, "= 0")  # conjunto dos IVs excluídos
     lh <- try(
       car::linearHypothesis(
         m, K,
         vcov = if (robust) sandwich::vcovHC(m, type="HC3") else NULL,
         test = "F"
       ),
       silent = TRUE
     )
     if (inherits(lh, "try-error")) {
       return(data.frame(endog = y, F = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_))
     }
     data.frame(
       endog = y,
       F     = unname(lh$F[2]),
       df1   = unname(lh$Df[2]),
       df2   = unname(lh$Res.Df[2]),
       p     = unname(lh$`Pr(>F)`[2]),
       row.names = NULL
     )
   }
   
   do.call(rbind, lapply(endo_ln, get_one))
 }
 
 # ---- usar:
 fs_hc3  <- fs_block_F(fit_q, df_iv, robust = TRUE)   # recomendado para o paper
 fs_homo <- fs_block_F(fit_q, df_iv, robust = FALSE)  # opcional (comparação)
 fs_hc3
 
 # ============================================================
 # Wu-Hausman manual (control function) por equação do QUAIDS
 # - Não usa AER; robusto a heteroscedasticidade (HC3)
 # - Funciona com múltiplos endógenos (todos os ln_preco_*)
 # - Usa os mesmos RHS e instrumentos do seu fit_quaids_manual_km1
 # ============================================================
 
 # ============================================================
 # Wu-Hausman manual (control function) por equação do QUAIDS
 # - Robusto a heteroscedasticidade (HC3)
 # - Múltiplos endógenos (todos ln_preco_*)
 # - Usa os mesmos RHS e IVs do fit_quaids_manual_km1
 # ============================================================
 
 wu_por_eq_manual2 <- function(fit_obj, data, robust = TRUE, verbose = FALSE){
   stopifnot(inherits(fit_obj, "quaids_km1_fit"))
   if (!requireNamespace("sandwich", quietly=TRUE)) stop("Precisa do pacote 'sandwich'.")
   if (!requireNamespace("car", quietly=TRUE))       stop("Precisa do pacote 'car'.")
   
   df <- as.data.frame(data)
   
   # 1) Termos
   rhs_all  <- fit_obj$rhs_terms
   inst_all <- fit_obj$inst_terms
   
   # limpar tokens de intercepto ou lixo textual
   clean_terms <- function(v){
     v <- unique(v)
     v[!(v %in% c("-1", "0", "+", "~", "|"))]
   }
   rhs_all  <- clean_terms(rhs_all)
   inst_all <- clean_terms(inst_all)
   
   endo_ln  <- grep("^ln_", rhs_all, value = TRUE)
   controls <- setdiff(rhs_all, endo_ln)
   
   if (length(endo_ln) == 0L) {
     warning("Nenhum regressor endógeno (^ln_) detectado em rhs_terms.")
   }
   
   # 2) Resíduos de 1ª etapa para TODOS endógenos, de uma vez
   #    (usando a MESMA base df; lm com na.exclude)
   resid_names <- character(0)
   for (x in endo_ln){
     rhs_fs <- unique(c(controls, inst_all))
     if (length(rhs_fs) == 0L) {
       # fallback: só intercepto
       f_fs <- reformulate(termlabels = 1, response = x)
     } else {
       f_fs <- reformulate(termlabels = rhs_fs, response = x)
     }
     mfs <- try(stats::lm(f_fs, data = df, na.action = na.exclude), silent = TRUE)
     rn  <- paste0("resid__", x)
     if (inherits(mfs, "try-error")) {
       if (verbose) message("Falha 1ª etapa para ", x, " (definindo resíduos como NA).")
       df[[rn]] <- NA_real_
     } else {
       r <- stats::residuals(mfs)
       # garantir comprimento N com NA para obs. perdidas
       if (length(r) != nrow(df)) r <- na.exclude(r)
       df[[rn]] <- as.numeric(r)
     }
     resid_names <- c(resid_names, rn)
   }
   
   # 3) LHS reais de cada equação (não confie em names(fit$eq))
   eq_lhs <- vapply(fit_obj$eq, function(m) as.character(formula(m)[[2]]), character(1))
   
   do_one <- function(y){
     # fórmula estrutural aumentada: y ~ rhs_all + resíduos de todos endógenos
     rhs_aug <- unique(c(rhs_all, resid_names))
     if (length(rhs_aug) == 0L) {
       return(data.frame(eq = y, k_endog = length(endo_ln),
                         stat = NA_real_, df = NA_real_, p_wu = NA_real_))
     }
     f_aug <- as.formula(paste(y, "~", paste(rhs_aug, collapse = " + ")))
     
     mols <- try(stats::lm(f_aug, data = df, na.action = na.exclude), silent = TRUE)
     if (inherits(mols, "try-error")) {
       if (verbose) message("Falha OLS aumentado para ", y)
       return(data.frame(eq = y, k_endog = length(endo_ln),
                         stat = NA_real_, df = NA_real_, p_wu = NA_real_))
     }
     
     # manter apenas resíduos que existem na base do modelo e não são todos NA
     present <- intersect(names(stats::model.frame(mols)), resid_names)
     keep_res <- present[vapply(present, function(nm) any(!is.na(df[[nm]])), logical(1))]
     if (length(keep_res) == 0L) {
       if (verbose) message("Sem resíduos válidos para testar em ", y)
       return(data.frame(eq = y, k_endog = length(endo_ln),
                         stat = NA_real_, df = 0, p_wu = NA_real_))
     }
     
     K <- paste(keep_res, "= 0")
     vc <- if (robust) sandwich::vcovHC(mols, type = "HC3") else stats::vcov(mols)
     
     # preferir Qui-quadrado; se falhar, tentar F
     lh <- try(car::linearHypothesis(mols, K, vcov = vc, test = "Chisq"), silent = TRUE)
     if (!inherits(lh, "try-error")) {
       stat <- unname(lh$Chisq[2]);  pval <- unname(lh$`Pr(>Chisq)`[2]); dfH0 <- unname(lh$Df[2])
       return(data.frame(eq = y, k_endog = length(endo_ln),
                         stat = stat, df = dfH0, p_wu = pval, row.names = NULL))
     } else {
       lh2 <- try(car::linearHypothesis(mols, K, vcov = vc, test = "F"), silent = TRUE)
       if (inherits(lh2, "try-error")) {
         if (verbose) message("Falha no linearHypothesis em ", y)
         return(data.frame(eq = y, k_endog = length(endo_ln),
                           stat = NA_real_, df = length(keep_res), p_wu = NA_real_))
       }
       stat <- unname(lh2$F[2]);  pval <- unname(lh2$`Pr(>F)`[2]); dfH0 <- unname(lh2$Df[2])
       return(data.frame(eq = y, k_endog = length(endo_ln),
                         stat = stat, df = dfH0, p_wu = pval, row.names = NULL))
     }
   }
   
   out_list <- lapply(eq_lhs, do_one)
   out <- do.call(rbind, out_list)
   rownames(out) <- NULL
   out
 }
 
 wu_tab <- wu_por_eq_manual2(fit_q, df_iv, robust = TRUE, verbose = TRUE)
 print(wu_tab)
 
 library(sandwich)
 library(lmtest)
 
 # ---------- helper: pega a equação por nome do share ----------
 get_eq_by_share <- function(fit, share_target = "w_despesahat3") {
   eq_list  <- fit$eq
   eq_names <- names(eq_list)
   lhs_vec  <- vapply(eq_list, function(m) as.character(formula(m)[[2]]), character(1))
   # 1) match exato na LHS
   idx <- which(lhs_vec == share_target)
   # 2) variantes sem "_"
   if (length(idx) == 0L) {
     plain <- gsub("_", "", share_target)
     idx <- which(eq_names == plain | eq_names == share_target |
                    grepl(paste0("\\b", plain, "\\b"), eq_names))
   }
   # 3) pelo dígito final (…3)
   if (length(idx) == 0L) {
     k <- sub(".*?(\\d+)\\s*$", "\\1", share_target)
     idx <- which(grepl(paste0(k, "$"), eq_names) | grepl(paste0(k, "$"), lhs_vec))
   }
   # 4) pela posição relativa em shareNames (excluindo a omitida)
   if (length(idx) == 0L && all(c("shareNames", "omit_share") %in% names(fit))) {
     sn    <- fit$shareNames
     omit  <- fit$omit_share
     sn_wo <- setdiff(sn, omit)
     pos   <- match(share_target, sn_wo)
     if (!is.na(pos)) idx <- pos
   }
   if (length(idx) == 0L) stop("Não encontrei a equação de ", share_target, " em fit$eq.")
   eq_list[[ idx[1] ]]
 }
 
 # ---------- helper: constrói fórmula segura ----------
 safe_formula_for <- function(fit, data, share_target = "w_despesahat3") {
   # tenta pegar do objeto fit
   f_try <- try({
     fm <- formula(get_eq_by_share(fit, share_target))
     fm
   }, silent = TRUE)
   if (!inherits(f_try, "try-error")) return(f_try)
   
   # fallback: reconstrói via rhs_terms, mantendo apenas colunas existentes em data
   rhs <- fit$rhs_terms
   rhs <- rhs[rhs %in% names(data)]
   if (length(rhs) == 0L) stop("rhs_terms não encontrado(s) no 'data'.")
   reformulate(rhs, response = share_target)
 }
 
 # ---------- 1) Fórmula da w3 ----------
 f_w3 <- safe_formula_for(fit_q, df_iv, "w_despesahat3")
 cat("Fórmula OLS w3:\n", deparse(f_w3), "\n")
 
 # ---------- 2) OLS ----------
 ols_w3 <- lm(f_w3, data = df_iv)
 cat("\n== Resumo OLS (EPs clássicos) ==\n")
 print(summary(ols_w3))
 
 # ---------- 3) EPs robustos HC3 ----------
 cat("\n== Coeficientes com EPs robustos (HC3) ==\n")
 print(coeftest(ols_w3, vcov = sandwich::vcovHC(ols_w3, type = "HC3")))
 
 # ---------- 4) (Opcional) EPs robustos clusterizados ----------
 if ("uf" %in% names(df_iv) || "f_reg" %in% names(df_iv)) {
   cluster_id <- if ("uf" %in% names(df_iv)) df_iv$uf else df_iv$f_reg
   vc_clust   <- sandwich::vcovCL(ols_w3, cluster = cluster_id, type = "HC3")
   cat("\n== Coeficientes com EPs robustos clusterizados (", 
       if ("uf" %in% names(df_iv)) "uf" else "f_reg", ") ==\n", sep = "")
   print(coeftest(ols_w3, vcov = vc_clust))
 }
 
 # ---------- 5) (Opcional) comparação OLS vs 3SLS na mesma equação ----------
 comp_try <- try({
   eq_w3 <- get_eq_by_share(fit_q, "w_despesahat3")
   comp  <- cbind(OLS = coef(ols_w3), `3SLS` = coef(eq_w3))
   cat("\n== Comparação de coeficientes: OLS vs 3SLS (w3) ==\n")
   print(round(comp, 6))
 }, silent = TRUE)
 
 library(lmtest)
 library(sandwich)
 library(car)
 
 ## --- RAMSEY RESET ---
 
 # 1) RESET clássico (H0: forma funcional correta)
 reset_cl <- resettest(ols_w3, power = 2:3, type = "fitted")
 print(reset_cl)
 
 # 2) RESET robusto (HC3) via comparação restrito vs. irrestrito
 #    Adiciona yhat^2 e yhat^3 (versão "fitted") e usa waldtest com EPs robustos
 df_aux <- within(df_iv, {
   yhat  <- fitted(ols_w3)
   yhat2 <- yhat^2
   yhat3 <- yhat^3
 })
 f_unres <- update(formula(ols_w3), . ~ . + yhat2 + yhat3)
 ols_w3_unres <- lm(f_unres, data = df_aux)
 
 reset_rob <- waldtest(
   ols_w3, ols_w3_unres,
   vcov = function(x) sandwich::vcovHC(x, type = "HC3")
 )
 print(reset_rob)
 
 ## Interpretação breve:
 ## - p alto (ex.: > 0.10) -> não rejeita H0 (sem evidência de má especificação de forma).
 ## - p baixo -> possível má especificação (omissão de termos não lineares/interações, etc.).
 
 
 ## --- VIF (colinearidade) ---
 
 vif_tab <- car::vif(ols_w3)
 print(vif_tab)
 
 # Se houver fatores, o car::vif retorna GVIF. Para comparabilidade:
 if (is.matrix(vif_tab)) {
   # para termos com Df > 1 (fatores), usa-se GVIF^(1/(2*Df))
   adj <- data.frame(
     Term  = rownames(vif_tab),
     GVIF  = vif_tab[,"GVIF"],
     Df    = vif_tab[,"Df"],
     GVIF_adj = vif_tab[,"GVIF"]^(1/(2*vif_tab[,"Df"]))
   )
   rownames(adj) <- NULL
   print(adj)
 }
 
 suppressPackageStartupMessages({
   library(dplyr)
   library(tidyr)
   library(tibble)
 })
 
 # Converte texto com vírgula decimal para numeric
 .as_num <- function(x) {
   if (is.numeric(x)) return(x)
   suppressWarnings(as.numeric(gsub(",", ".", x)))
 }
 
 ## Converte uma matriz (ou df) M para long: (eq, price, nome) sem list-cols
 mat_to_long <- function(M, nome){
   M <- as.matrix(M)
   if (is.null(rownames(M))) rownames(M) <- paste0("row", seq_len(nrow(M)))
   if (is.null(colnames(M))) colnames(M) <- paste0("col", seq_len(ncol(M)))
   data.frame(
     eq    = rep(rownames(M), times = ncol(M)),
     price = rep(colnames(M), each  = nrow(M)),
     tmp   = as.numeric(c(M)),
     stringsAsFactors = FALSE
   ) |>
     `names<-`(c("eq","price", nome))
 }
 
 ## Junta as estatísticas disponíveis (est, se, z, p, lwr, upr)
 to_tidy_elas2 <- function(elas_list){
   stopifnot(is.list(elas_list))
   stats <- intersect(c("est","se","z","p","lwr","upr"), names(elas_list))
   stopifnot(length(stats) > 0)
   dfs <- lapply(stats, function(s) mat_to_long(elas_list[[s]], s))
   out <- Reduce(function(a,b) merge(a,b, by=c("eq","price"), all=TRUE), dfs)
   # Garante tipos simples
   out$eq    <- as.character(out$eq)
   out$price <- as.character(out$price)
   for (s in stats) out[[s]] <- as.numeric(out[[s]])
   rownames(out) <- NULL
   out
 }
 

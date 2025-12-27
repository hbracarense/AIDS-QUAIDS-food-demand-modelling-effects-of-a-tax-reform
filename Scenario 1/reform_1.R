# ===============================================================
# Script consolidado — PREÇOS COM REFORMA (AIDS/QUAIDS + IVs)
# Patches incluídos:
#  (1) dummies f_mes/f_ym/f_tri/f_cap/f_int (c/ capital_interior)
#  (2) conjuntos IV_d_* e IV sets (core/mid/full) + saneamento
#  (3) 3SLS tolerante a IVs ausentes (intersecta/instData) + tentativas core/mid/full
# Data: 2025-09-02 | TZ: America/Sao_Paulo
# ===============================================================

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
path <- 'C:/Users/x16610962/Downloads/Reproduction/QUAIDS/Data/'
file <- 'banco_analise_AIDS.dta'
df   <- read_dta(paste0(path, file))

# Conjuntos principais
prices <- as.data.frame(df[, paste0("preco_com_reforma1", 1:6), drop = FALSE])
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
  eq_names <- names(fit$eq)
  r2 <- sapply(eq_names, function(eqname){
    y    <- model.response(model.frame(fit$eq[[eqname]]))
    yhat <- fitted(fit$eq[[eqname]])
    1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
  })
  names(r2) <- eq_names
  print(r2)
  invisible(tab)
}

.inst_qr_names <- function(vars, data, var_tol = 1e-12, qr_tol = 1e-10){
  if (length(vars) == 0) stop("Sem instrumentos candidatos.")
  Z <- model.matrix(stats::reformulate(vars, intercept = FALSE), data = data)
  # remove quase-constantes
  keep <- apply(Z, 2, function(v){
    vv <- var(v, na.rm = TRUE); is.finite(vv) && vv > var_tol
  })
  Z <- Z[, keep, drop = FALSE]
  if (ncol(Z) == 0) stop("Todas as colunas de Z foram filtradas por variância.")
  q <- qr(Z, tol = qr_tol)
  if (!is.finite(q$rank) || q$rank < 1L) stop("Posto de Z < 1 após QR.")
  colnames(Z)[ q$pivot[seq_len(q$rank)] ]
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

elas_quaids_manual <- function(fit_q, x, p = NULL, eps = 1e-4,
                               enforce_shares = TRUE,
                               enforce_hicks = TRUE,
                               enforce_symmetry = FALSE,
                               enforce_cournot = FALSE,
                               return_checks = TRUE){
  stopifnot(inherits(fit_q, "quaids_manual"))
  priceNames <- fit_q$priceNames; shareNames <- fit_q$shareNames; K <- length(priceNames)
  alpha <- fit_q$coef$alpha; beta <- fit_q$coef$beta; lambda <- fit_q$coef$lambda; gamma <- fit_q$coef$gamma
  
  if (is.null(p)) {
    ln_cols <- paste0("ln_", priceNames)
    stopifnot(all(ln_cols %in% names(fit_q$data)))
    p <- exp(colMeans(fit_q$data[, ln_cols, drop = FALSE], na.rm = TRUE)); names(p) <- priceNames
  } else { stopifnot(length(p)==K); if (is.null(names(p))) names(p) <- priceNames }
  
  wbar <- colMeans(fit_q$data[, shareNames, drop = FALSE], na.rm = TRUE); wbar <- wbar / sum(wbar)
  
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
    names(w) <- shareNames; w
  }
  
  w0 <- w_hat_fun(p, x); if (any(!is.finite(w0))) stop("Shares previstos NA/Inf")
  z0 <- if (fit_q$priceIndex == "Ls") log(x) - sum(wbar * log(p)) else log(x) - sum(w0 * log(p))
  eta <- 1 + (beta + 2*lambda*z0) / w0; names(eta) <- shareNames
  
  q0 <- (w0 * x) / p
  E_M <- matrix(NA_real_, K, K, dimnames = list(shareNames, priceNames))
  for (j in seq_len(K)) {
    p_up <- p; p_up[j] <- p_up[j] * (1 + eps)
    w_up <- w_hat_fun(p_up, x)
    q_up <- (w_up * x) / p_up
    E_M[, j] <- (log(q_up) - log(q0)) / (log(p_up[j]) - log(p[j]))
  }
  E_H <- E_M + outer(eta, w0)
  
  if (enforce_symmetry) E_H <- .symmetrize_slutsky(E_H, w0)
  if (enforce_hicks)    E_H <- .enforce_hicks_adding(E_H, w0)
  E_M <- E_H - outer(eta, w0)
  if (enforce_cournot) { E_M <- .enforce_cournot_marshall(E_M, eta, w0); E_H <- E_M + outer(eta, w0) }
  
  out <- list(at = list(p = setNames(p, priceNames), x = x, w = w0, z = z0),
              expenditure = eta, marshall = E_M, hicks = E_H)
  if (return_checks) out$checks <- .check_identities(E_H, E_M, eta, w0)
  out
}

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

run_quaids_3sls <- function(instNames, label, use_z2 = TRUE) {
  cat("\n---- Tentando set:", label, " (|IV|=", length(instNames), ", use_z2=", use_z2, ") ----\n", sep = "")
  ans <- try(fit_quaids_manual_km1(
    prices     = df_iv[, priceNames, drop = FALSE],
    shares     = df_iv[, shareNames, drop = FALSE],
    x          = df_iv[["gasto_total_atualhat"]],
    priceIndex = "Ls",
    estMethod  = "3SLS",
    omit_share = 1,
    drop_price = 1,
    instNames  = instNames,
    instData   = df_iv,
    use_z2     = use_z2,
    maxiter    = 500
  ), silent = TRUE)
  if (inherits(ans, "try-error")) {
    cat("Falhou. Erro:\n", as.character(ans), "\n", sep = "")
    return(NULL)
  }
  ans
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

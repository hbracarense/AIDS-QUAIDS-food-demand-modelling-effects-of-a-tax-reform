## ================== QUAIDS Delta CIs (plug-and-play) ==================
library(openxlsx)

## Requisitos: numDeriv
if (!requireNamespace("numDeriv", quietly = TRUE))
  stop("Instale: install.packages('numDeriv')")

## =========================================================
## Aux: Just-ID builder (usa objetos/funções já do projeto)
## =========================================================
make_eq_justID2 <- function(fit_base, data, fit_eqspec, eqn,
                            prefer_keep = character(0),
                            avoid = character(0),
                            verbose = TRUE) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  stopifnot(is.list(inst_by_eq), eqn %in% names(inst_by_eq))
  
  ## fórmula da eq-alvo a partir de fit_base
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  f <- formula(fit_base$eq[[ match(eqn, eq_names) ]])
  
  cur <- unique(inst_by_eq[[eqn]])
  if (length(avoid)) cur <- setdiff(cur, avoid)
  
  dfJ_of <- function(S) {
    terms <- unique(c(base_exog, S))
    ivm <- try(AER::ivreg(f, instruments = reformulate(terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) return(NA_real_)
    X <- try(stats::model.matrix(ivm, component = "regressors"),  silent = TRUE)
    Z <- try(stats::model.matrix(ivm, component = "instruments"), silent = TRUE)
    if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NA_real_)
    qr(Z)$rank - qr(X)$rank
  }
  
  cur_dfJ <- dfJ_of(cur)
  if (isTRUE(cur_dfJ == 0)) {
    if (verbose) message("[", eqn, "] já está JUST/UNDER-ID (df_J <= 0). Nada a fazer.")
    return(fit_eqspec)
  }
  
  keep_best <- NULL
  for (k in seq_len(length(cur))) {
    for (drop_set in utils::combn(cur, k, simplify = FALSE)) {
      S <- setdiff(cur, drop_set)
      if (length(prefer_keep) && any(!(prefer_keep %in% S))) next
      dj <- dfJ_of(S)
      if (isTRUE(dj == 0)) { keep_best <- S; break }
    }
    if (!is.null(keep_best)) break
  }
  
  if (is.null(keep_best)) {
    if (verbose) message("[", eqn, "] não encontrei subset JUST-ID. Mantendo conjunto atual.")
    return(fit_eqspec)
  }
  
  inst_by_eq[[eqn]] <- unique(keep_best)
  if (verbose) {
    dropped <- setdiff(cur, keep_best)
    msg <- if (length(dropped)) paste(paste0("-", dropped), collapse = " ") else "(sem drops)"
    message("[", eqn, "] JUST-ID por drop: ", msg)
  }
  refit_quaids_try_eqlist(fit_base, data, inst_by_eq)
}

## --------- Util: detecção de nomes no seu dataset ---------
.detect_shares <- function(dat) {
  s <- grep("^w_", names(dat), value = TRUE)
  if (!length(s)) stop("Não encontrei colunas de shares que comecem com 'w_'.")
  ord <- order(suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", s))))
  s[ord]
}
.detect_lnexp <- function(dat) {
  cand <- intersect(c("lndespesa","ln_despesa","lnexp","ln_expenditure"), names(dat))
  if (length(cand)) return(cand[1])
  pick <- grep("^ln.*(desp|exp)", names(dat), value = TRUE)
  if (!length(pick)) stop("Não encontrei a variável ln(x) (ex.: 'lndespesa').")
  pick[1]
}
.detect_sq <- function(dat) {
  cand <- intersect(c("z2","lndespesa2","lnexp2"), names(dat))
  if (length(cand)) return(cand[1])
  NULL
}

## CORREÇÃO (crítica):
##  - prioriza baseline ln_preco_por_kg# quando existir
##  - caso não exista, cai no regex antigo (mas ordenando por dígito final)
.detect_lnprices <- function(dat) {
  p_base <- grep("^ln_preco_por_kg\\d+$", names(dat), value = TRUE)
  if (length(p_base)) {
    ord <- order(as.integer(gsub(".*?(\\d+)$", "\\1", p_base)))
    return(p_base[ord])
  }
  p <- grep("^ln_.*preco", names(dat), value = TRUE)
  if (!length(p)) stop("Não encontrei colunas de preços em log (ex.: 'ln_preco_por_kg#').")
  ord <- order(suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", p))))
  p[ord]
}

## --------- Empacotamento/Desempacotamento θ ----------
.build_theta_from_systemfit <- function(fit_sys, dat) {
  cf  <- stats::coef(fit_sys)
  Vcf <- stats::vcov(fit_sys)
  if (is.null(names(cf))) stop("Coeficientes sem nomes no systemfit.")
  
  ## IMPORTANTE: os coeficientes do seu systemfit usam eq1_, eq2_, ...
  ## Então usamos labels canônicos por posição, independentemente de names(fit_sys$eq)
  L <- length(fit_sys$eq)
  eq_labels <- paste0("eq", seq_len(L))  # "eq1","eq2",...
  eq_lhs    <- vapply(fit_sys$eq, function(m) as.character(stats::formula(m)[[2]]), character(1)) # w_despesahat2..6
  
  shares <- .detect_shares(dat)  # w_despesahat1..6
  K <- length(shares)
  
  ## share omitida = share no data que não aparece como LHS nas equações estimadas
  share_omit <- setdiff(shares, eq_lhs)
  if (length(share_omit) != 1) {
    stop("Não consegui identificar exatamente 1 equação omitida (shares no data vs. LHS das eqs no systemfit).")
  }
  share_omit <- share_omit[1]
  
  xvar  <- .detect_lnexp(dat)
  z2var <- .detect_sq(dat)
  pvars <- .detect_lnprices(dat)   # idealmente ln_preco_por_kg1..6
  
  ## Detecta preços presentes olhando diretamente para os nomes dos coeficientes:
  ## precisa existir algo como "eqk_<pvar>" para algum k
  present_price_vars <- pvars[vapply(pvars, function(pj) {
    any(grepl(paste0("^eq\\d+_", pj, "$"), names(cf)))
  }, logical(1))]
  
  if (!length(present_price_vars)) {
    stop("Não detectei coeficientes de preços no systemfit com o padrão 'eq#_<pvar>' (ex.: eq1_ln_preco_por_kg2).")
  }
  
  ## θ = [ alpha(eq!=omit) , beta(eq!=omit) , lambda(eq!=omit) , vec(Gamma para eq!=omit e pvars presentes) ]
  theta <- numeric(0); tnames <- character(0); pick_idx <- integer(0)
  
  ## α (intercepto)
  for (e in seq_len(L)) {
    eq <- eq_labels[e]
    nm <- paste0(eq, "_(Intercept)")
    val <- if (nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("alpha|", eq_lhs[e]))
    pick_idx <- c(pick_idx, match(nm, names(cf)))
  }
  
  ## β (ln despesa)
  for (e in seq_len(L)) {
    eq <- eq_labels[e]
    nm <- paste0(eq, "_", xvar)
    val <- if (nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("beta|", eq_lhs[e]))
    pick_idx <- c(pick_idx, match(nm, names(cf)))
  }
  
  ## λ ((ln desp)^2)
  for (e in seq_len(L)) {
    eq <- eq_labels[e]
    nm <- if (!is.null(z2var)) paste0(eq, "_", z2var) else NA_character_
    val <- if (!is.na(nm) && nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("lambda|", eq_lhs[e]))
    pick_idx <- c(pick_idx, if (!is.na(nm)) match(nm, names(cf)) else NA_integer_)
  }
  
  ## Γ (preços)
  for (e in seq_len(L)) {
    eq <- eq_labels[e]
    for (pj in present_price_vars) {
      nm <- paste0(eq, "_", pj)
      val <- if (nm %in% names(cf)) cf[[nm]] else 0
      theta <- c(theta, val)
      tnames <- c(tnames, paste0("gamma|", eq_lhs[e], "|", pj))
      pick_idx <- c(pick_idx, match(nm, names(cf)))
    }
  }
  
  ## V_theta
  p <- length(theta)
  Vth <- matrix(0, p, p, dimnames = list(tnames, tnames))
  have <- which(!is.na(pick_idx))
  if (length(have)) {
    idx <- pick_idx[have]
    Vth[have, have] <- Vcf[idx, idx, drop = FALSE]
  }
  
  list(theta = theta, V = Vth, tnames = tnames,
       eq_labels = eq_labels, eq_lhs = eq_lhs,
       shares = shares, share_omit = share_omit,
       xvar = xvar, z2var = z2var, pvars = pvars,
       present_price_vars = present_price_vars)
}

## Reconstrói α, β, λ, Γ completos (K bens x K preços) a partir de θ
.rebuild_params_quaids <- function(theta, meta) {
  eq_lhs <- meta$eq_lhs           ## bens estimados (LHS): w_despesahat2..6
  shares <- meta$shares           ## todos bens no data: w_despesahat1..6
  present_price_vars <- meta$present_price_vars
  pvars_all <- meta$pvars
  K <- length(shares)
  
  idx_good <- function(nm) as.integer(sub(".*?(\\d+)$", "\\1", nm))
  eq_idx   <- idx_good(eq_lhs)
  omit_idx <- idx_good(meta$share_omit)
  
  L <- length(eq_lhs)
  a_vec <- theta[1:L]
  b_vec <- theta[(L+1):(2*L)]
  l_vec <- theta[(2*L+1):(3*L)]
  rest  <- theta[-(1:(3*L))]
  
  ## Γ observada: L x length(present_price_vars)
  G_obs <- matrix(rest, nrow = L, byrow = TRUE,
                  dimnames = list(eq_lhs, present_price_vars))
  
  ## α, β, λ de tamanho K
  alpha <- beta <- lambda <- rep(NA_real_, K)
  alpha[eq_idx]  <- a_vec
  beta[eq_idx]   <- b_vec
  lambda[eq_idx] <- l_vec
  
  ## adding-up para o bem omitido
  alpha[omit_idx]  <- 1 - sum(alpha[-omit_idx], na.rm = TRUE)
  beta[omit_idx]   <- - sum(beta[-omit_idx],  na.rm = TRUE)
  lambda[omit_idx] <- - sum(lambda[-omit_idx],na.rm = TRUE)
  
  ## Γ completa KxK
  Gamma <- matrix(0, K, K, dimnames = list(shares, pvars_all))
  Gamma[eq_idx, present_price_vars] <- G_obs
  
  ## coluna(s) faltante(s) por homogeneidade
  missing_cols <- setdiff(pvars_all, present_price_vars)
  if (length(missing_cols)) {
    for (i in eq_idx) {
      Gamma[i, missing_cols] <- -sum(Gamma[i, present_price_vars, drop = TRUE])
    }
  }
  
  ## linha omitida por adding-up por coluna
  for (pj in pvars_all) {
    Gamma[omit_idx, pj] <- -sum(Gamma[setdiff(seq_len(K), omit_idx), pj])
  }
  
  ## simetria
  Gamma <- 0.5*(Gamma + t(Gamma))
  
  list(alpha = alpha, beta = beta, lambda = lambda, Gamma = Gamma)
}


## --------- Elasticidades QUAIDS (Banks–Blundell–Lewbel) ----------
.quaids_elasticities_from_params <- function(alpha, beta, lambda, Gamma,
                                             lnx, pvec, wvec,
                                             shares, pvars) {
  K <- length(shares)
  if (length(pvec) != K) stop("Comprimento de pvec diferente de K.")
  if (length(wvec) != K) stop("Comprimento de wvec diferente de K.")
  
  lnp <- log(pvec)
  lna <- sum(alpha * lnp) + 0.5 * as.numeric(t(lnp) %*% Gamma %*% lnp)
  lnb <- sum(beta * lnp); b <- exp(lnb)
  L   <- lnx - lna
  
  mu <- beta + (2 * lambda / b) * L
  dlnA_dlnp <- alpha + as.numeric(Gamma %*% lnp)
  
  L2 <- L^2
  mu_ij <- Gamma - outer(mu, dlnA_dlnp) - outer(lambda / b * L2, beta)
  
  e_x <- mu / pmax(wvec, 1e-12) + 1
  delta <- diag(K)
  e_u <- sweep(mu_ij, 1, pmax(wvec, 1e-12), "/") - delta
  e_c <- e_u + outer(e_x, wvec)
  
  list(marshall = e_u, hicks = e_c, expenditure = e_x,
       lna = lna, b = b, L = L)
}

## --------- Função principal: Delta CI para QUAIDS ----------
delta_quaids_CIs <- function(fit_obj, dat, level = 0.95, method = c("pointwise","bonferroni")) {
  method <- match.arg(method)
  
  ## aceita systemfit direto ou wrapper com $fit
  fit_sys <- if (inherits(fit_obj, "systemfit")) fit_obj else fit_obj$fit
  if (!inherits(fit_sys, "systemfit")) stop("Passe o objeto systemfit (ou um wrapper que tenha $fit).")
  
  ## robustez: nomeia equações se necessário
  if (is.null(names(fit_sys$eq)) || any(names(fit_sys$eq) == "")) {
    lhs <- vapply(fit_sys$eq, function(m) as.character(stats::formula(m)[[2]]), character(1))
    names(fit_sys$eq) <- lhs
  }
  
  meta   <- .build_theta_from_systemfit(fit_sys, dat)
  theta0 <- meta$theta
  Vth    <- meta$V
  
  shares <- meta$shares
  pvars  <- meta$pvars
  K      <- length(shares)
  
  ## ponto de avaliação
  w_star <- sapply(shares, function(s) stats::median(dat[[s]], na.rm = TRUE))
  p_star <- sapply(pvars,  function(v) exp(mean(dat[[v]], na.rm = TRUE)))
  lnx_star <- stats::median(dat[[ meta$xvar ]], na.rm = TRUE)
  
  g_eval <- function(th) {
    par <- .rebuild_params_quaids(th, meta)
    E <- .quaids_elasticities_from_params(
      alpha  = par$alpha,
      beta   = par$beta,
      lambda = par$lambda,
      Gamma  = par$Gamma,
      lnx    = lnx_star,
      pvec   = p_star,
      wvec   = w_star,
      shares = shares,
      pvars  = pvars
    )
    c(as.vector(t(E$hicks)), as.vector(t(E$marshall)), as.numeric(E$expenditure))
  }
  
  g0 <- g_eval(theta0)
  
  ## Jacobiana numérica
  J  <- numDeriv::jacobian(g_eval, theta0)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  lo <- g0 - z*se; hi <- g0 + z*se
  
  rn <- shares; cn <- shares
  k2 <- K*K
  idx_H <- seq_len(k2)
  idx_M <- k2 + seq_len(k2)
  idx_E <- (2*k2) + seq_len(K)
  
  H_hat <- matrix(g0[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_lo  <- matrix(lo[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_hi  <- matrix(hi[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_se  <- matrix(se[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  
  M_hat <- matrix(g0[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_lo  <- matrix(lo[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_hi  <- matrix(hi[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_se  <- matrix(se[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  
  eta_hat <- setNames(g0[idx_E], rn)
  eta_lo  <- setNames(lo[idx_E], rn)
  eta_hi  <- setNames(hi[idx_E], rn)
  eta_se  <- setNames(se[idx_E], rn)
  
  mk_tbl <- function(M, lo, hi, se) {
    out <- do.call(rbind, lapply(seq_len(K), function(i) {
      data.frame(i = rn[i], j = cn, est = M[i,], se = se[i,], lwr = lo[i,], upr = hi[i,], row.names = NULL)
    }))
    rownames(out) <- NULL; out
  }
  
  list(
    marshallian = mk_tbl(M_hat, M_lo, M_hi, M_se),
    hicksian    = mk_tbl(H_hat, H_lo, H_hi, H_se),
    expenditure = data.frame(good = rn, est = eta_hat, se = eta_se, lwr = eta_lo, upr = eta_hi, row.names = NULL),
    level       = level,
    method      = paste0("Delta (QUAIDS) — w* mediana; x*=median(", meta$xvar, "); p*=exp(mean ln p)")
  )
}

## =========================================================
## Exemplo de uso (não falha se objetos não existirem)
## =========================================================
if (exists("fit_q") && exists("df_iv")) {
  
  ## (1) escolher um fit_start já criado no pipeline (ou baseline)
  fit_start <-
    if (exists("fit_sw2")) fit_sw2 else
      if (exists("fit_tryK")) fit_tryK else
        if (exists("fit_bumped")) fit_bumped else
          if (exists("fit_fix")) fit_fix else
            fit_q
  
  ## (2) construir fit_just6 (se as funções do projeto existirem)
  if (exists("refit_quaids_try_eqlist") && exists("get_base_exog")) {
    
    fit_just6 <- make_eq_justID2(
      fit_base     = fit_q,
      data         = df_iv,
      fit_eqspec   = fit_start,
      eqn          = "w_despesahat6",
      prefer_keep  = "IV_d_uf_X12",
      avoid        = c("IV_d_uf_X11"),
      verbose      = TRUE
    )
    
    if (!("ln_preco_por_kg1" %in% names(df_iv))) df_iv$ln_preco_por_kg1 <- 0
    out_delta <- delta_quaids_CIs(fit_just6, df_iv, level = 0.95)
    
    
    cat("\n== ELASTICIDADES — MARSHALL (não compensadas) — IC Delta ==\n")
    print(out_delta$marshallian, row.names = FALSE)
    cat("\n== ELASTICIDADES — HICKS (compensadas) — IC Delta ==\n")
    print(out_delta$hicksian, row.names = FALSE)
    cat("\n== ELASTICIDADES — RENDA/DESPESA — IC Delta ==\n")
    print(out_delta$expenditure, row.names = FALSE)
    
  } else {
    message("Pulando exemplo: faltam get_base_exog() e/ou refit_quaids_try_eqlist() no ambiente.")
  }
}

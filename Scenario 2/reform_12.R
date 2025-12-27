## ============================
## DEPENDÊNCIAS
## ============================
suppressPackageStartupMessages({
  library(AER)
  library(sandwich)
  library(lmtest)
  library(car)
  library(dplyr)
})

## ============================
## HELPERS ROBUSTOS
## ============================
.nzv_drop <- function(M, tol = 1e-12){
  if (is.null(dim(M)) || NCOL(M) == 0) return(M)
  keep <- apply(M, 2, function(x) sd(x, na.rm=TRUE) > tol)
  M[, keep, drop = FALSE]
}

.hansenJ_matrix <- function(ivfit){
  ## J robusto (HC0) diretamente nas matrizes com os MESMOS pesos do 2SLS
  mf <- model.frame(ivfit)
  y  <- model.response(mf)
  X  <- model.matrix(ivfit, component = "regressors")
  Zf <- model.matrix(ivfit, component = "instruments")
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, nrow(X))
  s  <- sqrt(w/mean(w)); yS <- y*s; XS <- X*s; ZS <- Zf*s
  
  ## 2SLS manual
  ZtZ <- crossprod(ZS)
  PZ  <- ZS %*% qr.solve(ZtZ, t(ZS))
  XtPZX <- crossprod(XS, PZ %*% XS)
  XtPZy <- crossprod(XS, PZ %*% yS)
  b2    <- solve(XtPZX, XtPZy)
  e2    <- yS - XS %*% b2
  
  ## separa Z_only = Z \ X para diagnosticar sobre-identificação
  Xn <- colnames(XS)
  Z_only <- ZS[, setdiff(colnames(ZS), Xn), drop = FALSE]
  Z_only <- .nzv_drop(Z_only)
  Z_use  <- cbind(XS, Z_only)   # para momentos completos
  
  n  <- nrow(Z_use)
  m  <- crossprod(Z_use, e2) / n
  S  <- crossprod(Z_use * as.numeric(e2)) / n
  S  <- (S + t(S))/2
  
  ## ridge mínimo se S mal-condicionada
  ridges <- c(0, 1e-8, 1e-6, 1e-4)
  ridge_used <- NA_real_; ok <- FALSE; sol <- NULL
  for (rg in ridges){
    Srg <- S
    if (rg > 0) diag(Srg) <- diag(Srg) + rg * mean(diag(Srg))
    sol <- try(qr.solve(Srg, m), silent = TRUE)
    if (!inherits(sol, "try-error") && all(is.finite(sol))) { ok <- TRUE; ridge_used <- rg; break }
  }
  J   <- if (ok) drop(n * crossprod(m, sol)) else NA_real_
  dfJ <- NCOL(Z_use) - NCOL(XS)
  pJ  <- if (!is.na(J) && dfJ > 0) pchisq(J, dfJ, lower.tail = FALSE) else NA_real_
  
  data.frame(statistic = J, df = dfJ, p.value = pJ, ridge = ridge_used, check.names = FALSE)
}

.wu_hausman_durbin <- function(endog_vars, exog_incl, iv_excl, dat, wcol){
  ## 1º estágios (cada endógeno em Z + demais endógenos)
  Rhat <- lapply(endog_vars, function(v){
    rhs <- c(exog_incl, iv_excl, setdiff(endog_vars, v))
    f1  <- reformulate(rhs, response = v)
    fm1 <- if (is.null(wcol)) lm(f1, dat) else lm(f1, dat, weights = dat[[wcol]])
    resid(fm1)
  })
  names(Rhat) <- endog_vars
  dat_dw <- cbind(dat, setNames(Rhat, paste0(endog_vars, "_res")))
  
  ## 2º estágio “residual-inclusion”
  yvar <- all.vars(attr(terms(reformulate(endog_vars, response="y")), "variables"))[2]  # placeholder
  ## inferir y da notação usual das shares (pega primeira *_res para localizar y na fórmula abaixo)
  ## O y correto será passado de fora; aqui só construímos a fórmula geral:
  f_dw <- reformulate(c(exog_incl, endog_vars, paste0(endog_vars, "_res")),
                      response = attr(dat_dw, "yvar_override"))
  stop("Internal")  # (ver função externa que define yvar corretamente)
}

.reset_wls_2sls <- function(y, endog_vars, exog_incl, iv_excl, dat, wcol){
  ## constroi Xhat com 1º estágios e roda WLS “2º estágio OLS” para RESET
  Xhat <- lapply(endog_vars, function(v){
    rhs <- c(exog_incl, iv_excl, setdiff(endog_vars, v))
    fm1 <- if (is.null(wcol)) lm(reformulate(rhs, v), dat)
    else               lm(reformulate(rhs, v), dat, weights = dat[[wcol]])
    fitted(fm1)
  })
  names(Xhat) <- paste0(endog_vars, "_hat")
  dat2 <- cbind(dat, as.data.frame(Xhat))
  f2   <- reformulate(c(exog_incl, names(Xhat)), response = y)
  wls2 <- if (is.null(wcol)) lm(f2, dat2) else lm(f2, dat2, weights = dat2[[wcol]])
  lmtest::resettest(wls2, power = 2:3, type = "regressor", vcov. = function(x) sandwich::vcovHC(x, type = "HC1"))
}

.sw_F <- function(endog_vars, exog_incl, iv_excl, dat, wcol){
  do.call(rbind, lapply(endog_vars, function(v){
    rhs <- c(exog_incl, iv_excl, setdiff(endog_vars, v))
    fm1 <- if (is.null(wcol)) lm(reformulate(rhs, response=v), dat)
    else               lm(reformulate(rhs, response=v), dat, weights = dat[[wcol]])
    Zpresent <- intersect(iv_excl, colnames(model.matrix(fm1)))
    if (length(Zpresent) == 0)
      return(data.frame(regressor = v, df1 = NA, df2 = NA, `F-statistic` = NA, `p-value` = NA))
    lh <- car::linearHypothesis(fm1, Zpresent, vcov. = sandwich::vcovHC(fm1, type = "HC1"), test = "F")
    data.frame(regressor = v,
               df1 = unname(lh$Df[2]), df2 = unname(lh$Res.Df[2]),
               `F-statistic` = unname(lh$F[2]), `p-value` = unname(lh$`Pr(>F)`[2]),
               check.names = FALSE)
  }))
}

.kp_rkWaldF <- function(ivfit){
  if (requireNamespace("ivpack", quietly = TRUE)) {
    out <- try(ivpack::KP.stat(ivfit), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(data.frame(rkWaldF = unname(out$KP.rk), df = unname(out$df), p.value = unname(out$p.value)))
    }
  }
  ## Fallback informativo (média dos SW-F será inserida depois pelo chamador)
  data.frame(rkWaldF = NA_real_, df = NA_real_, p.value = NA_real_)
}

## ============================
## TESTES PRINCIPAIS
## ============================
# === Requisitos ===
suppressPackageStartupMessages({
  library(AER); library(sandwich); library(lmtest); library(car)
})

run_iv_tests <- function(fit2, dat, endog_vars, exog_incl, iv_excl, wcol=NULL,
                         tol=1e-10, ridge_seq=c(0,1e-8,1e-6,1e-4)) {
  
  nzv_drop <- function(M, thr=1e-12){
    if (is.null(dim(M)) || NCOL(M)==0) return(M)
    keep <- apply(M, 2, function(x) sd(x, na.rm=TRUE) > thr)
    M[, keep, drop=FALSE]
  }
  prune_Ziv_QR <- function(X, Zf){
    Xn  <- colnames(X)
    Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop=FALSE]
    Ziv <- nzv_drop(Ziv)
    if (NCOL(Ziv) > 0){
      Ziv_s <- scale(Ziv, center=TRUE, scale=FALSE)
      qrz   <- qr(Ziv_s, LAPACK=TRUE)
      r     <- qrz$rank
      cols  <- if (r>0) sort(qrz$pivot[seq_len(r)]) else integer(0)
      Zsel  <- if (length(cols)) Ziv[, cols, drop=FALSE] else Ziv[, 0, drop=FALSE]
    } else Zsel <- Ziv
    cbind(X, Zsel)
  }
  safe_2sls <- function(y, X, Z, w){
    s  <- sqrt(w/mean(w)); yS <- y*s; XS <- X*s; ZS <- Z*s
    ZtZ <- crossprod(ZS)
    ok <- FALSE
    for (rg in ridge_seq){
      ZtZ_rg <- ZtZ
      if (rg>0) diag(ZtZ_rg) <- diag(ZtZ_rg) + rg*mean(diag(ZtZ_rg))
      PZ_X <- try(ZS %*% qr.solve(ZtZ_rg, crossprod(ZS, XS)), silent=TRUE)
      PZ_y <- try(ZS %*% qr.solve(ZtZ_rg, crossprod(ZS, yS)),  silent=TRUE)
      if (!inherits(PZ_X,"try-error") && !inherits(PZ_y,"try-error") &&
          all(is.finite(PZ_X)) && all(is.finite(PZ_y))) { ok <- TRUE; break }
    }
    if (!ok) stop("Falha ao inverter Z'Z mesmo com ridge.")
    XtPZ_X <- crossprod(XS, PZ_X); XtPZ_y <- crossprod(XS, PZ_y)
    ok2 <- FALSE
    for (rg2 in ridge_seq){
      XtPZ_X_rg <- XtPZ_X
      if (rg2>0) diag(XtPZ_X_rg) <- diag(XtPZ_X_rg) + rg2*mean(diag(XtPZ_X_rg))
      beta <- try(solve(XtPZ_X_rg, XtPZ_y), silent=TRUE)
      if (!inherits(beta,"try-error") && all(is.finite(beta))) { ok2 <- TRUE; break }
    }
    if (!ok2) stop("Falha ao inverter X'(P_Z)X.")
    e <- as.numeric(yS - XS %*% beta)
    list(beta=beta, e=e, XS=XS, ZS=ZS)
  }
  hansenJ_HC0 <- function(e, ZS, p){
    n <- NROW(ZS); m <- crossprod(ZS, e)/n
    S <- crossprod(ZS*as.numeric(e))/n; S <- (S+t(S))/2
    ok <- FALSE; ridge_used <- NA_real_
    for (rg in ridge_seq){
      Srg <- S; if (rg>0) diag(Srg) <- diag(Srg) + rg*mean(diag(Srg))
      sol <- try(qr.solve(Srg, m), silent=TRUE)
      if (!inherits(sol,"try-error") && all(is.finite(sol))) { ok <- TRUE; ridge_used <- rg; break }
    }
    if (!ok) return(list(J=NA_real_, df=NA_integer_, p=NA_real_, ridge="S_singular"))
    J  <- drop(n * crossprod(m, sol)); df <- NCOL(ZS) - p
    list(J=if (df>0) J else NA_real_, df=if (df>0) df else 0L,
         p=if (df>0) pchisq(J, df, lower.tail=FALSE) else NA_real_, ridge=ridge_used)
  }
  
  # --- objetos do ajuste ---
  vc2 <- vcovHC(fit2, type="HC1")
  mf  <- model.frame(fit2); y <- model.response(mf)
  X   <- model.matrix(fit2, "regressors")
  Zf  <- model.matrix(fit2, "instruments")
  w   <- model.weights(mf); if (is.null(w)) w <- rep(1, nrow(X))
  
  # --- Wald bloco (robusto), tolerante a alias ---
  wald_block <- try(
    linearHypothesis(fit2, paste0(endog_vars, " = 0"),
                     vcov.=vc2, test="F", singular.ok=TRUE),
    silent=TRUE
  )
  
  # --- Hansen-J (HC0) com poda QR ---
  Z   <- prune_Ziv_QR(X, Zf)
  est <- safe_2sls(y, X, Z, w)
  J   <- hansenJ_HC0(est$e, est$ZS, p=ncol(X))
  hansenJ <- data.frame(statistic=J$J, df=J$df, p.value=J$p, ridge=J$ridge)
  
  # --- Wu–Hausman (Durbin-Wu) robusto (permitindo alias) ---
  Rhat <- lapply(endog_vars, function(v){
    rhs <- unique(c(exog_incl, iv_excl, setdiff(endog_vars, v)))
    f1  <- reformulate(rhs, response=v)
    fm  <- if (is.null(wcol)) lm(f1, dat) else lm(f1, dat, weights=dat[[wcol]])
    list(res=resid(fm), fit=fitted(fm))
  })
  names(Rhat) <- endog_vars
  dat_DW <- dat
  for (v in endog_vars) dat_DW[[paste0(v,"_res")]] <- Rhat[[v]]$res
  f_dw  <- reformulate(c(exog_incl, endog_vars, paste0(endog_vars,"_res")),
                       response = all.vars(formula(fit2))[1])
  fm_dw <- if (is.null(wcol)) lm(f_dw, dat_DW) else lm(f_dw, dat_DW, weights=dat_DW[[wcol]])
  # removendo resíduos aliased do teste
  cf_dw <- coef(fm_dw); res_ok <- paste0(endog_vars,"_res")
  res_ok <- res_ok[res_ok %in% names(cf_dw)[!is.na(cf_dw)]]
  dwh <- if (length(res_ok)) {
    linearHypothesis(fm_dw, res_ok, vcov.=vcovHC(fm_dw, type="HC1"), test="F", singular.ok=TRUE)
  } else {
    structure(list(message="todos os resíduos aliased"), class="try-error")
  }
  
  # --- Sanderson–Windmeijer F por endógeno (filtrando IVs aliased) ---
  sw_rows <- lapply(endog_vars, function(v){
    rhs <- unique(c(exog_incl, iv_excl, setdiff(endog_vars, v)))
    f1  <- reformulate(rhs, response=v)
    fm1 <- if (is.null(wcol)) lm(f1, dat) else lm(f1, dat, weights=dat[[wcol]])
    cf1 <- coef(fm1)
    iv_ok <- setdiff(intersect(iv_excl, names(cf1)[!is.na(cf1)]), "(Intercept)")
    if (!length(iv_ok)) {
      return(data.frame(regressor=v, df1=NA, df2=NA, F=NA, p.value=NA,
                        note="IVs todos aliased", check.names=FALSE))
    }
    lh <- linearHypothesis(fm1, iv_ok, vcov.=vcovHC(fm1, type="HC1"), test="F", singular.ok=TRUE)
    data.frame(regressor=v,
               df1=unname(lh$Df[2]), df2=unname(lh$Res.Df[2]),
               F=unname(lh$F[2]), p.value=unname(lh$`Pr(>F)`[2]),
               check.names=FALSE)
  })
  sw <- do.call(rbind, sw_rows)
  
  # --- RESET (com X_hat) ---
  dat_hat <- dat
  for (v in endog_vars) dat_hat[[paste0(v,"_hat")]] <- Rhat[[v]]$fit
  f2sls_ols   <- reformulate(c(exog_incl, paste0(endog_vars,"_hat")),
                             response = all.vars(formula(fit2))[1])
  fm_2sls_ols <- if (is.null(wcol)) lm(f2sls_ols, dat_hat) else lm(f2sls_ols, dat_hat, weights=dat_hat[[wcol]])
  # --- Ramsey RESET (WLS 2SLS stage-2, HC1) ---
  reset <- try({
    dep <- all.vars(formula(fit2))[1]
    mf2 <- model.frame(fit2)
    w   <- model.weights(mf2); if (is.null(w)) w <- rep(1, nrow(mf2))
    
    # termos da equação de 2º estágio (sem o intercepto, que a fórmula inclui)
    base_terms <- unique(c(exog_incl, endog_vars))
    
    # predito linear do 2º estágio (mesmo X e mesmos pesos)
    X2    <- model.matrix(fit2, "regressors")
    yhat2 <- as.numeric(X2 %*% coef(fit2))
    
    # Data para o RESET (mesmo n das observações usadas)
    dat_reset <- cbind(dat, .__yhat__ = yhat2)
    
    f_base <- reformulate(base_terms, response = dep)
    f_aug  <- reformulate(c(base_terms, "I(.__yhat__^2)","I(.__yhat__^3)"), response = dep)
    
    m0 <- if (is.null(wcol)) lm(f_base, dat_reset) else lm(f_base, dat_reset, weights = dat_reset[[wcol]])
    m1 <- if (is.null(wcol)) lm(f_aug,  dat_reset) else lm(f_aug,  dat_reset, weights = dat_reset[[wcol]])
    
    wt <- lmtest::waldtest(m0, m1, vcov = sandwich::vcovHC(m1, type="HC1"), test="F")
    
    data.frame(
      F        = unname(wt$F[2]),
      df1      = unname(wt$Df[2]),
      df2      = unname(wt$Res.Df[2]),
      p.value  = unname(wt$`Pr(>F)`[2]),
      check.names = FALSE
    )
  }, silent = TRUE)
  
  if (inherits(reset, "try-error") || NROW(reset) == 0L) {
    reset <- data.frame(F = NA_real_, df1 = NA_integer_, df2 = NA_integer_, p.value = NA_real_)
  }
  
 
  
  
  # --- Sargan clássico (sem pesos) ---
  fit_unw <- AER::ivreg(formula(fit2), data = dat)
  eu <- resid(fit_unw); Zu <- model.matrix(fit_unw, "instruments"); Xu <- model.matrix(fit_unw, "regressors")
  lm_eZ   <- lm(eu ~ Zu - 1)
  R2      <- summary(lm_eZ)$r.squared
  n_u     <- NROW(Zu); df_sarg <- NCOL(Zu) - NCOL(Xu)
  sargan  <- data.frame(statistic = n_u*R2,
                        df = df_sarg,
                        p.value = if (df_sarg>0) pchisq(n_u*R2, df_sarg, lower.tail=FALSE) else NA_real_)
  
  kp <- try(data.frame(rkWaldF = mean(sw$F, na.rm=TRUE)), silent=TRUE)
  if (inherits(kp, "try-error")) kp <- data.frame(rkWaldF = NA_real_)
  
  list(
    wald_block = wald_block,
    hansenJ    = hansenJ,
    wu_hausman = dwh,
    swF        = sw,
    reset      = reset,
    sargan     = sargan,
    kpF        = kp
  )
}



## ============================
## MONTAGEM DO MODELO (Cenário 2 “melhor”)
## ============================
## Objetos já existentes no seu script:
##   - df_iv          (dados)
##   - priceNames     (endógenos: logs dos preços do Cenário 2)
##   - shareNames     (dependentes: w_despesahat1..6)
##   - iv_set_core / iv_set_mid / iv_set_full  (blocos de IVs)
##
## Exógenas no 2º estágio (iguais às usadas no Cenário 2)
# --- checagens e montagem robusta da especificação ---
exog_alvo   <- c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4))
exog_incl   <- intersect(exog_alvo, names(df_iv))
iv_excl     <- intersect(iv_set_mid, names(df_iv))
price_endog <- intersect(priceNames, names(df_iv))

faltando_exog  <- setdiff(exog_alvo,   exog_incl)
faltando_iv    <- setdiff(iv_set_mid,  iv_excl)
faltando_price <- setdiff(priceNames,  price_endog)

if (length(price_endog) == 0L) stop("Nenhum dos endógenos (priceNames) está no df_iv.")
if (length(faltando_price)) message("Aviso: endógenos ausentes no df_iv: ", paste(faltando_price, collapse=", "))
if (length(faltando_exog))  message("Aviso: exógenas ausentes no df_iv: ",  paste(faltando_exog,   collapse=", "))
if (length(faltando_iv))    message("Aviso: IVs ausentes no df_iv: ",       paste(faltando_iv,     collapse=", "))

dep_best <- shareNames[5]  # equação 5 (modelo “melhor”)

form_best <- as.formula(paste0(
  dep_best, " ~ ", paste(c(price_endog, exog_incl), collapse = " + "),
  " | ",            paste(c(exog_incl,  iv_excl),   collapse = " + ")
))

fit2_best <- AER::ivreg(form_best, data = df_iv, weights = df_iv$peso_final)

# --- roda todos os testes importantes (usa sua função já criada) ---
tests_best <- run_iv_tests(
  fit2       = fit2_best,
  dat        = df_iv,
  endog_vars = price_endog,
  exog_incl  = exog_incl,
  iv_excl    = iv_excl,
  wcol       = "peso_final"
)

# --- tabelas resumidas prontas para reporte ---
wald_tbl   <- data.frame(test="Wald (block of prices)",
                         F=unname(tests_best$wald_block$F[2]),
                         df1=unname(tests_best$wald_block$Df[2]),
                         df2=unname(tests_best$wald_block$Res.Df[2]),
                         p_val=unname(tests_best$wald_block$`Pr(>F)`[2]),
                         check.names=FALSE)

hansen_tbl <- transform(tests_best$hansenJ, test="Hansen-J (HC0, ridge min)")[,c("test","statistic","df","p.value","ridge")]

wu_tbl     <- data.frame(test="Wu–Hausman (residual inclusion, HC1)",
                         F=unname(tests_best$wu_hausman$F[2]),
                         df1=unname(tests_best$wu_hausman$Df[2]),
                         df2=unname(tests_best$wu_hausman$Res.Df[2]),
                         p_val=unname(tests_best$wu_hausman$`Pr(>F)`[2]),
                         check.names=FALSE)

## --- Sanderson–Windmeijer F, por endógeno ---
sw_tbl <- within(tests_best$swF, {
  test      <- paste0("Sanderson–Windmeijer F (", regressor, ")")
  statistic <- F
  `p-value` <- p.value
})
sw_tbl <- sw_tbl[, c("test","df1","df2","statistic","p-value")]

## --- Kleibergen–Paap rk Wald F ---
if (requireNamespace("ivpack", quietly = TRUE)) {
  kp <- try(ivpack::KP.stat(fit2_best), silent = TRUE)
  if (!inherits(kp, "try-error")) {
    kp_tbl <- data.frame(
      test      = "Kleibergen–Paap rk Wald F",
      rkWaldF   = unname(kp$rk),          # estatística rk Wald F
      df_num    = NA_integer_,            # df usuais não são reportados pelo KP
      df_den    = NA_integer_,
      `p-value` = unname(kp$rk.p),
      check.names = FALSE
    )
  } else {
    kp_tbl <- data.frame(
      test    = "Kleibergen–Paap rk Wald F (aprox. média SW)",
      rkWaldF = tests_best$kpF$rkWaldF[1],
      check.names = FALSE
    )
  }
} else {
  kp_tbl <- data.frame(
    test    = "Kleibergen–Paap rk Wald F (aprox. média SW)",
    rkWaldF = tests_best$kpF$rkWaldF[1],
    check.names = FALSE
  )
}

reset_tbl <- within(tests_best$reset, {
  test <- "Ramsey RESET (WLS 2SLS stage-2, HC1)"
})[, c("test","F","df1","df2","p.value")]


sargan_tbl <- transform(tests_best$sargan, test="Sargan (unweighted)")[,c("test","statistic","df","p.value")]

# visualizar rapidamente
print(wald_tbl); print(hansen_tbl); print(wu_tbl); print(sw_tbl); print(kp_tbl); print(reset_tbl); print(sargan_tbl)

run_iv_tests_all <- function(dep_vars = shareNames,
                             df = df_iv,
                             prices = price_endog,
                             exog  = exog_incl,
                             ivs   = iv_excl,
                             wcol  = "peso_final") {
  
  out <- lapply(dep_vars, function(dep){
    form <- as.formula(paste0(
      dep, " ~ ", paste(c(prices, exog), collapse=" + "),
      " | ",       paste(c(exog,   ivs), collapse=" + ")
    ))
    fit  <- AER::ivreg(form, data=df, weights=df[[wcol]])
    tst  <- run_iv_tests(fit2 = fit, dat = df,
                         endog_vars = prices, exog_incl = exog,
                         iv_excl = ivs, wcol = wcol)
    list(dep=dep, tests=tst)
  })
  names(out) <- dep_vars
  
  bind <- function(extract) {
    do.call(rbind, lapply(out, function(o){
      tab <- extract(o$tests)
      if (!is.null(tab) && nrow(tab)) cbind(equation = o$dep, tab) else NULL
    }))
  }
  
  list(
    wald   = bind(function(x)
      data.frame(F = unname(x$wald_block$F[2]),
                 df1 = unname(x$wald_block$Df[2]),
                 df2 = unname(x$wald_block$Res.Df[2]),
                 `p-value` = unname(x$wald_block$`Pr(>F)`[2]),
                 check.names = FALSE)),
    hansen = bind(function(x)
      transform(x$hansenJ, `p-value` = p.value)),
    wu     = bind(function(x)
      data.frame(F = unname(x$wu_hausman$F[2]),
                 df1 = unname(x$wu_hausman$Df[2]),
                 df2 = unname(x$wu_hausman$Res.Df[2]),
                 `p-value` = unname(x$wu_hausman$`Pr(>F)`[2]),
                 check.names = FALSE)),
    sw     = bind(function(x)
      setNames(x$swF[, c("regressor","df1","df2","F","p.value")],
               c("regressor","df1","df2","statistic","p-value"))),
    reset  = bind(function(x) setNames(x$reset, c("F","df1","df2","p-value"))),
    sargan = bind(function(x)
      setNames(x$sargan[, c("statistic","df","p.value")],
               c("statistic","df","p-value"))),
    kp     = bind(function(x) x$kpF)
  )
}

## -------------------------
## Objetos do seu Cenário 2
## -------------------------
iv_excl     <- iv_set_mid          # bloco de IVs do “melhor modelo”
price_endog <- priceNames          # logs dos preços endógenos
exog_incl   <- intersect(c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4)),
                         names(df_iv))   # garante que existem em df_iv

## -------------------------
## Chamada (todas as 6 equações)
## -------------------------
all_tests <- run_iv_tests_all(
  dep_vars = shareNames,   # w_despesahat1..6
  df       = df_iv,
  prices   = price_endog,
  exog     = exog_incl,
  ivs      = iv_excl,
  wcol     = "peso_final"
)

## -------------------------
## Exemplos de acesso/uso
## -------------------------
all_tests$wald     # Wald do bloco de preços por equação
all_tests$hansen   # Hansen-J (HC0, ridge) por equação
all_tests$wu       # Wu–Hausman (residual inclusion) por equação
all_tests$sw       # Sanderson–Windmeijer F por endógeno e equação
all_tests$reset    # Ramsey RESET por equação
all_tests$sargan   # Sargan (sem pesos) por equação
all_tests$kp       # KP rk Wald F (ou aprox), por equação

## (opcional) Filtrar a “equação 5”
subset(all_tests$wald, equation == shareNames[5])
subset(all_tests$hansen, equation == shareNames[5])


## =========================================================
## Weak-IV robust inference (AR) + LIML/Fuller — plug-and-play
## Requer: df_iv, priceNames, shareNames, iv_set_core/mid/full
## Pesos: df_iv$peso_final
## Pacotes: AER (usa 'ivreg' por baixo), sandwich, lmtest, car
## (opcional): 'ivreg' >= 0.5 para k-class LIML/Fuller
## =========================================================

suppressPackageStartupMessages({
  library(AER)
  library(sandwich)
  library(lmtest)
  library(car)
  ## 'ivreg' (novo) é opcional; se existir, tentamos kclass="liml"/Fuller
  have_ivreg <- requireNamespace("ivreg", quietly = TRUE)
})

## -------- 0) Especificação (mesma do cenário 2) ----------------
build_spec_cen2 <- function(dat, iv_set = NULL) {
  stopifnot(all(c("peso_final") %in% names(dat)))
  exog_target  <- c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4))
  exog_incl    <- intersect(exog_target, names(dat))
  price_endog  <- intersect(priceNames, names(dat))
  if (is.null(iv_set)) iv_set <- iv_set_mid
  iv_excl      <- intersect(iv_set, names(dat))
  
  if (!length(price_endog)) stop("Nenhum log de preço de cenário 2 encontrado em df_iv.")
  if (!length(exog_incl))   warning("Nenhuma exógena sazonal/z encontrada em df_iv.")
  if (!length(iv_excl))     warning("Nenhum IV do bloco escolhido encontrado em df_iv.")
  
  list(price_endog = price_endog, exog_incl = exog_incl, iv_excl = iv_excl)
}

## -------- 1) Ajuste 2SLS/LIML/Fuller por equação ---------------
## no topo do script:
have_ivreg_pkg <- requireNamespace("ivreg", quietly = TRUE)
has_kclass_arg <- FALSE
if (have_ivreg_pkg) {
  ## checa se a sua build exporta o argumento 'kclass'
  fmls <- try(formals(ivreg::ivreg), silent = TRUE)
  if (!inherits(fmls, "try-error")) has_kclass_arg <- ("kclass" %in% names(fmls))
}

fit_iv_three <- function(dep, dat, price_endog, exog_incl, iv_excl, wcol="peso_final") {
  rhs_2s  <- paste(c(price_endog, exog_incl), collapse = " + ")
  rhs_iv  <- paste(c(exog_incl, iv_excl),    collapse = " + ")
  fml     <- as.formula(paste0(dep, " ~ ", rhs_2s, " | ", rhs_iv))
  
  fit_2s  <- AER::ivreg(fml, data = dat, weights = dat[[wcol]])
  
  fit_liml <- NULL
  fit_fuller <- NULL
  if (have_ivreg_pkg && has_kclass_arg) {
    fit_liml <- try(ivreg::ivreg(fml, data = dat, weights = dat[[wcol]], kclass = "liml"),
                    silent = TRUE)
    if (inherits(fit_liml, "try-error")) fit_liml <- NULL
    
    fit_fuller <- try(ivreg::ivreg(fml, data = dat, weights = dat[[wcol]],
                                   kclass = "Fuller", alpha = 1),
                      silent = TRUE)
    if (inherits(fit_fuller, "try-error")) fit_fuller <- NULL
  }
  list(fit_2s = fit_2s, fit_liml = fit_liml, fit_fuller = fit_fuller)
}


## -------- 2) Anderson–Rubin (HC1) para o BLOCO de preços -------
## H0: beta_preços = 0 (teste conjunto, robusto a weak-IV)
## Implementação: y* = y − X_preço * 0 = y;
## Regressão WLS de y* em [exógenas + IVs] e teste F de que coeficientes de IVs = 0
## Anderson–Rubin (HC1) para o BLOCO de preços — versão “anti-alias”
ar_block_test <- function(dep, dat, exog_incl, iv_excl, wcol="peso_final") {
  if (!length(iv_excl)) {
    return(data.frame(equation = dep, F = NA_real_, df1 = NA_integer_, df2 = NA_integer_,
                      `p-value` = NA_real_, note = "sem IVs no conjunto", check.names = FALSE))
  }
  f_aux <- as.formula(paste0(dep, " ~ ", paste(c(exog_incl, iv_excl), collapse = " + ")))
  wls   <- lm(f_aux, data = dat, weights = dat[[wcol]])
  
  ## <<< TRUQUE-CHAVE: testar só IVs com coeficiente estimado (não-aliased) >>>
  cf      <- coef(wls)
  iv_used <- intersect(iv_excl, names(cf)[!is.na(cf)])
  iv_used <- setdiff(iv_used, "(Intercept)")
  
  if (!length(iv_used)) {
    return(data.frame(equation = dep, F = NA_real_, df1 = NA_integer_, df2 = NA_integer_,
                      `p-value` = NA_real_, note = "todos os IVs ficaram aliased no WLS", check.names = FALSE))
  }
  
  lh <- car::linearHypothesis(
    wls, iv_used,
    vcov. = sandwich::vcovHC(wls, type = "HC1"),
    test  = "F",
    singular.ok = TRUE
  )
  
  data.frame(
    equation = dep,
    F        = unname(lh$F[2]),
    df1      = unname(lh$Df[2]),
    df2      = unname(lh$Res.Df[2]),
    `p-value`= unname(lh$`Pr(>F)`[2]),
    note     = if (length(setdiff(iv_excl, iv_used))) "IVs aliased removidos do teste" else NA_character_,
    check.names = FALSE
  )
}


## -------- 3) Tabelas de coeficientes (preços) com HC1 ----------
coef_table_prices <- function(fit, price_endog) {
  vc <- sandwich::vcovHC(fit, type = "HC1")
  cf <- lmtest::coeftest(fit, vcov. = vc)
  keep <- rownames(cf) %in% price_endog
  if (!any(keep)) return(NULL)
  out <- data.frame(regressor = rownames(cf)[keep],
                    estimate  = cf[keep, "Estimate"],
                    std.error = cf[keep, "Std. Error"],
                    statistic = cf[keep, "t value"],
                    p.value   = cf[keep, "Pr(>|t|)"],
                    row.names = NULL, check.names = FALSE)
  out
}

## -------- 4) Orquestrador: roda i) AR e ii) LIML/Fuller --------
run_i_and_ii <- function(dat, shareNames, price_endog, exog_incl, iv_excl, wcol = "peso_final") {
  fits <- vector("list", length(shareNames)); names(fits) <- shareNames
  tabs_2s <- tabs_liml <- tabs_full <- vector("list", length(shareNames)); names(tabs_2s) <- names(tabs_liml) <- names(tabs_full) <- shareNames
  ar_list <- vector("list", length(shareNames)); names(ar_list) <- shareNames
  
  for (dep in shareNames) {
    f3 <- fit_iv_three(dep, dat, price_endog, exog_incl, iv_excl, wcol)
    fits[[dep]] <- f3
    
    ## Tabelas de coeficientes (preços) — 2SLS, LIML, Fuller
    tabs_2s[[dep]]   <- coef_table_prices(f3$fit_2s,   price_endog)
    tabs_liml[[dep]] <- if (!is.null(f3$fit_liml))   coef_table_prices(f3$fit_liml, price_endog)   else NULL
    tabs_full[[dep]] <- if (!is.null(f3$fit_fuller)) coef_table_prices(f3$fit_fuller, price_endog) else NULL
    
    ## Anderson–Rubin (HC1) para o bloco de preços
    ar_list[[dep]] <- ar_block_test(dep, dat, exog_incl, iv_excl, wcol)
  }
  
  ## Empilha resultados
  tab_ar   <- do.call(rbind, ar_list)
  tab_2s   <- do.call(rbind, lapply(names(tabs_2s),   function(k) if (is.null(tabs_2s[[k]])) NULL else transform(tabs_2s[[k]],   equation = k)))
  tab_liml <- do.call(rbind, lapply(names(tabs_liml), function(k) if (is.null(tabs_liml[[k]])) NULL else transform(tabs_liml[[k]], equation = k)))
  tab_full <- do.call(rbind, lapply(names(tabs_full), function(k) if (is.null(tabs_full[[k]])) NULL else transform(tabs_full[[k]], equation = k)))
  
  list(
    fits   = fits,
    AR_HC1 = tab_ar[ , c("equation","F","df1","df2","p-value")],
    coef_2SLS   = tab_2s[ , c("equation","regressor","estimate","std.error","statistic","p.value")],
    coef_LIML   = if (!is.null(tab_liml)) tab_liml[ , c("equation","regressor","estimate","std.error","statistic","p.value")] else NULL,
    coef_Fuller = if (!is.null(tab_full)) tab_full[ , c("equation","regressor","estimate","std.error","statistic","p.value")] else NULL
  )
}

## =========================
## >>> CHAMADA (plug-and-play)
## =========================
## Escolha do conjunto de IVs do Cenário 2 (troque se quiser)
spec <- build_spec_cen2(df_iv, iv_set = iv_set_mid)

res_i_ii <- run_i_and_ii(
  dat         = df_iv,
  shareNames  = shareNames,
  price_endog = spec$price_endog,
  exog_incl   = spec$exog_incl,
  iv_excl     = spec$iv_excl,
  wcol        = "peso_final"
)

## Visualização rápida:
print(res_i_ii$AR_HC1)        # Anderson–Rubin (HC1) para o bloco de preços (robusto a weak-IV)
head(res_i_ii$coef_2SLS)      # Coeficientes 2SLS (HC1) para cada preço por equação
head(res_i_ii$coef_LIML)      # (se disponível)
head(res_i_ii$coef_Fuller)    # (se disponível)

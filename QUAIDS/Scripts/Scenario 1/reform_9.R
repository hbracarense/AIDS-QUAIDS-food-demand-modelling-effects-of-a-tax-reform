## =========================================================
##  JUST-IDENTIFIED DEMAND (AIDS/QUAIDS) — TOP-TIER PACK
##  P0: Setup  |  P1: IVs = 1 PC/price  |  P2: 2SLS + CL
##  P3: Força (SW F e R2p)               P4: AR/CLR (poi)
##  P5: Elasticidades                    P6: Bootstrap (opt)
## =========================================================

## ---------- P0) Setup ----------
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("pls","AER","sandwich","lmtest","dplyr","purrr","tibble"))

has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)

stopifnot(exists("df_iv_ok"), exists("priceNames"), exists("shareNames"),
          exists("best_fit"), exists("best_iv_set"), is.character(best_iv_set))

drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
exogs      <- intersect(c("z","z2"), names(df_iv_ok))          # pode ser char(0)
eqs        <- as.character(setdiff(shareNames, omit_share))
endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))    # todos os ln preços endógenos

set.seed(1234)

## ---------- Utils básicos ----------
# residualização robusta (sem apply; evita problemas de classe)
resid_on_vec <- function(v, X) {
  X <- if (is.null(X) || ncol(as.matrix(X))==0) NULL else as.matrix(X)
  v <- as.numeric(v)
  if (is.null(X)) return(v)
  XtX <- crossprod(X)
  v - as.numeric(X %*% solve(XtX, crossprod(X, v)))
}
resid_cols_on <- function(Z, X) {
  Z <- as.matrix(Z)
  X <- if (is.null(X) || ncol(as.matrix(X))==0) NULL else as.matrix(X)
  if (is.null(X)) return(Z)
  XtX <- crossprod(X)
  Z - X %*% solve(XtX, crossprod(X, Z))
}

# F SW condicional (exclusão conjunta de Z no 1º estágio)
first_stage_F <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if(length(X)) X else "1"
  rhs1 <- c(X, Z)
  f0 <- lm(stats::reformulate(rhs0, response=D), data=data)
  f1 <- lm(stats::reformulate(rhs1, response=D), data=data)
  as.numeric(anova(f0, f1)$F[2])
}
partial_R2 <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if(length(X)) X else "1"
  rhs1 <- c(X, Z)
  f0 <- lm(stats::reformulate(rhs0, response=D), data=data)
  f1 <- lm(stats::reformulate(rhs1, response=D), data=data)
  1 - sum(residuals(f1)^2)/sum(residuals(f0)^2)
}

## ---------- P1) IVs just-identified: 1 PC por preço (PLS condicional) ----------
make_pc_for_price <- function(poi, df, iv_pool, priceNames, drop_price, exogs=c("z","z2")){
  stopifnot(poi %in% names(df))
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xn <- c(intersect(exogs, names(df)), intersect(ln_others, names(df)))
  Zn <- intersect(iv_pool, names(df))
  use <- unique(c(poi, Xn, Zn))
  dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
  
  D  <- as.numeric(dat[[poi]])
  X  <- if (length(Xn)) as.matrix(dat[, Xn, drop=FALSE]) else NULL
  Z0 <- as.matrix(dat[, Zn, drop=FALSE])
  
  # residualiza condicional em X
  D_t <- resid_on_vec(D, X)
  Z_t <- resid_cols_on(Z0, X)
  
  # remove quase-constantes
  sd_ok <- apply(Z_t, 2, sd, na.rm=TRUE) > 1e-12
  Z_t   <- Z_t[, sd_ok, drop=FALSE]
  stopifnot(ncol(Z_t) > 0)
  
  # PLS com 1 componente (just-ID)
  fit <- pls::plsr(D_t ~ ., data=as.data.frame(Z_t), ncomp=1, validation="none", scale=TRUE)
  s1  <- drop(pls::scores(fit)[,1])
  nm  <- paste0("pc_", sub("^ln_", "", poi), "_1")
  dat[[nm]] <- s1
  list(data=dat, pc_name=nm)
}

gen_all_pcs <- function(df, priceNames, drop_price, iv_pool, exogs=c("z","z2")){
  endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))
  out_df <- df; pc_names <- character(0)
  for (poi in endogs_all){
    pc <- make_pc_for_price(poi, out_df, iv_pool, priceNames, drop_price, exogs)
    out_df  <- pc$data
    pc_names <- c(pc_names, pc$pc_name)
  }
  list(data=out_df, iv_star=pc_names)
}

# dicionário base (apenas IVs que existem na base)
iv_pool <- intersect(best_iv_set, names(df_iv_ok))
if (!length(iv_pool)) stop("iv_pool vazio: 'best_iv_set' não bate com colunas de df_iv_ok.")

pcs <- gen_all_pcs(df_iv_ok, priceNames, drop_price, iv_pool, exogs)
dfJI    <- pcs$data                   # base com PCs
iv_star <- as.character(pcs$iv_star)  # 1 PC por ln_preço (mesmo comprimento de endogs_all)
stopifnot(length(iv_star) == length(endogs_all))

cat("\n[P1] IVs (just-ID) criados:\n"); print(iv_star)

## ---------- P2) 2SLS eq-a-eq (VCOV cluster robusta) ----------
cluster_var <- c("cluster_id","id_domicilio","id_municipio","id_setor_censitario")
cluster_var <- intersect(cluster_var, names(dfJI))
if (!length(cluster_var)) { dfJI$.__rowid__ <- seq_len(nrow(dfJI)); cluster_var <- ".__rowid__" }
cluster_var <- cluster_var[1]

fit_ivreg_one <- function(dat, y, endogs, exogs, iv_names, cluster){
  rhs_y <- c(endogs, exogs)                 # regressors estruturais
  inst  <- unique(c(exogs, iv_names))       # instrumentos: X + PCs
  fy <- stats::reformulate(rhs_y, response=y)
  fz <- stats::reformulate(inst)
  fit <- AER::ivreg(formula=fy, instruments=fz, data=dat)
  vc  <- tryCatch(sandwich::vcovCL(fit, cluster=dat[[cluster]]),
                  error=function(e) sandwich::vcovHC(fit, type="HC1"))
  list(fit=fit, vcov=vc, f_y=fy, f_z=fz)
}

ivsys <- lapply(eqs, function(y){
  use <- unique(c(y, endogs_all, exogs, iv_star, cluster_var))
  dat <- dfJI[stats::complete.cases(dfJI[, use, drop=FALSE]), , drop=FALSE]
  info <- fit_ivreg_one(dat, y, endogs_all, exogs, iv_star, cluster_var)
  list(eq=y, n=nrow(dat), fit=info$fit, vcov=info$vcov, f_y=info$f_y, f_z=info$f_z)
})
names(ivsys) <- eqs

ivsys_tidy <- function(ivsys){
  do.call(rbind, lapply(ivsys, function(x){
    ct <- lmtest::coeftest(x$fit, vcov.=x$vcov)
    data.frame(eq=x$eq, term=rownames(ct),
               estimate=ct[,1], std.error=ct[,2], z=ct[,3], p.value=ct[,4],
               n=x$n, row.names=NULL)
  }))
}

tab_iv <- ivsys_tidy(ivsys)
cat("\n[P2] 2SLS (just-ID) | VCOV cluster — principais coeficientes:\n")
print(subset(tab_iv, term %in% c("(Intercept)", endogs_all)), digits=3)

## ---------- P3) Força do 1º estágio ----------
# Mapeia cada preço ao seu PC correspondente
map_price_pc <- setNames(iv_star, endogs_all)

sw_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  use <- unique(c(Dv, exogs, iv_star))
  dat <- dfJI[stats::complete.cases(dfJI[, use, drop=FALSE]), , drop=FALSE]
  # F conjunto (todos os PCs) e F do instrumento próprio (PC do Dv)
  F_all <- first_stage_F(Dv, exogs, iv_star, dat)
  F_one <- first_stage_F(Dv, exogs, map_price_pc[Dv], dat)
  tibble::tibble(
    price   = Dv,
    n_cc    = nrow(dat),
    F_SW_all= F_all,
    F_SW_own= F_one,
    R2p_all = partial_R2(Dv, exogs, iv_star, dat),
    df1_all = length(iv_star),
    df1_own = 1L,
    df2     = nrow(dat) - length(c(exogs, iv_star)) - 1
  )
})
cat("\n[P3] Força 1º estágio (SW F e R² parcial):\n"); print(sw_tbl, digits=3)

## ---------- P4) AR/CLR para o preço de interesse ----------
price_of_interest <- "ln_preco_com_reforma13"   # ajuste se quiser

ar_clr_oneeq <- function(eq_y, poi, halfwidth=200){
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  
  # X = exógenos + outros ln_p; Z = TODOS os PCs (instruments para todos endogs)
  Xnames <- c(exogs, ln_others)
  use <- unique(c(eq_y, poi, Xnames, iv_star))
  dat <- dfJI[stats::complete.cases(dfJI[, use, drop=FALSE]), , drop=FALSE]
  
  Y <- as.numeric(dat[[eq_y]])
  D <- as.numeric(dat[[poi]])
  X <- if(length(Xnames)) as.matrix(dat[, Xnames, drop=FALSE]) else NULL
  Z <- as.matrix(dat[, iv_star,  drop=FALSE])
  
  # centro da grade no 2SLS
  b2 <- tryCatch(unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]), error=function(e) 0)
  grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out=4001)
  
  if (has_ivmodel && ncol(as.matrix(D))==1) {
    ivm <- ivmodel::ivmodel(Y=Y, D=D, Z=Z, X=X)
    p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0=b0)$p.value, 0.0)
    p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0=b0)$p.value, 0.0)
    keepA <- which(p_ar  >= 0.05); keepC <- which(p_clr >= 0.05)
    ciA <- if(length(keepA)) c(min(grid[keepA]), max(grid[keepA])) else c(NA,NA)
    ciC <- if(length(keepC)) c(min(grid[keepC]), max(grid[keepC])) else c(NA,NA)
  } else {
    # Fallback robusto (multi-endógenos): AR via sistema empilhado
    p_ar <- vapply(grid, function(b0){
      fit <- lm(cbind(Y - b0*D, X) ~ 0 + Z) # testa exclusão de Z no empilhado
      1 - pf(summary(aov(fit))[[1]][["F value"]][1], ncol(Z), nrow(Z) - ncol(Z))
    }, 0.0)
    keepA <- which(p_ar >= 0.05)
    ciA <- if(length(keepA)) c(min(grid[keepA]), max(grid[keepA])) else c(NA,NA)
    ciC <- ciA
  }
  tibble::tibble(eq=eq_y, AR_low=ciA[1], AR_high=ciA[2], CLR_low=ciC[1], CLR_high=ciC[2])
}

arclr_tbl <- purrr::map_dfr(eqs, ~ar_clr_oneeq(.x, price_of_interest, halfwidth=200))
cat("\n[P4] AR/CLR (weak-IV robust) | coef de ", price_of_interest, ":\n", sep="")
print(arclr_tbl, digits=3)

## ---------- P5) Elasticidades (Γ de 2SLS → objeto do modelo) ----------
# extrai Γ^IV (coeficientes nos ln p_j) de cada equação e impõe adding-up
gamma_from_ivsys <- function(ivsys, priceNames, drop_price, shareNames, omit_share){
  eqs <- setdiff(shareNames, omit_share)
  pj  <- paste0("ln_", setdiff(priceNames, drop_price))
  G   <- matrix(NA_real_, nrow=length(shareNames), ncol=length(pj),
                dimnames=list(shareNames, pj))
  for (eq in eqs) {
    co <- coef(ivsys[[eq]]$fit)
    G[eq, pj] <- unname(co[pj])
  }
  # linha da share omitida pela restrição de soma
  miss <- omit_share
  G[miss, pj] <- -colSums(G[eqs, pj, drop=FALSE], na.rm=TRUE)
  G
}

Giv <- gamma_from_ivsys(ivsys, priceNames, drop_price, shareNames, omit_share)

fitJI <- best_fit
# mapeia nomes (remove "ln_")
cols_fit <- colnames(fitJI$coef$gamma)           # ex.: "preco_com_reforma12", ...
pj_short <- sub("^ln_", "", colnames(Giv))
fitJI$coef$gamma[, match(pj_short, cols_fit)] <- Giv

# ponto de avaliação (x*, p*)
px_emp <- list(
  x_star = median(dfJI$gasto_total_atualhat, na.rm=TRUE),
  p_star = exp(colMeans(dfJI[paste0("ln_", priceNames)], na.rm=TRUE))
)

stopifnot(exists("elas_quaids_manual"))
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


E <- elas_quaids_manual(
  fitJI,
  x = px_emp$x_star, p = px_emp$p_star,
  normalize_eval="softmax", normalize_deriv="softmax",
  dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
)

cat("\n[P5] Elasticidades (just-ID, PCs):\n")
cat("diag(Marshall):\n"); print(round(diag(E$marshall), 3))
cat("\nEta (expenditure):\n"); print(round(E$expenditure, 3))
cat("\nHicks:\n"); print(round(diag(E$hicks), 3))

## ---------- P6) (Opcional) Bootstrap cluster para elasticidades ----------
# Exemplo de esqueleto (desative se não for usar):
# cluster_boot_elas <- function(B=199){
#   set.seed(123)
#   cl <- dfJI[[cluster_var]]
#   cls <- unique(cl)
#   acc_M <- array(NA_real_, dim=c(B, length(priceNames), length(priceNames)))
#   acc_eta <- matrix(NA_real_, nrow=B, ncol=length(priceNames))
#   for (b in 1:B) {
#     pick <- sample(cls, length(cls), replace=TRUE)
#     idx  <- unlist(lapply(pick, function(k) which(cl==k)), use.names=FALSE)
#     d    <- dfJI[idx, , drop=FALSE]
#     # (re-estimar PCs para cada b se quiser total honestidade; aqui reusa PCs)
#     # re-estimar ivsys, gerar Giv_b, substituir em fitJI_b e calcular E_b ...
#     # acc_M[b,,] <- diag(Eb$marshall); acc_eta[b,] <- Eb$expenditure
#   }
#   list(M=acc_M, eta=acc_eta)
# }

cat("\n=== FIM: P0..P5 concluídos (just-ID). KP-rk F não se aplica em just-ID; reporte SW F (own/all), R² parcial, AR/CLR e elasticidades. ===\n")

## ====================== P-SYS: 3SLS + DELTA (sem bootstrap) ======================

## ====================== P-SYS: 3SLS + DELTA (sem bootstrap) ======================

## Requisitos
stopifnot(requireNamespace("systemfit", quietly = TRUE))
stopifnot(requireNamespace("dplyr",     quietly = TRUE))
stopifnot(requireNamespace("tibble",    quietly = TRUE))
stopifnot(requireNamespace("numDeriv",  quietly = TRUE))
stopifnot(exists("best_fit"), exists("priceNames"), exists("shareNames"))
stopifnot(exists("dfJI"), exists("iv_star"))              # base com PCs por preço (just-ID)
stopifnot(exists("elas_quaids_manual"))                   # sua função de elasticidades

## ------------ 0) Utilidades genéricas ------------
replace_gamma_return_E <- function(fit_template, Giv, x_star, p_star,
                                   normalize_eval="softmax", normalize_deriv="softmax"){
  fit2 <- fit_template
  cols_fit <- colnames(fit2$coef$gamma)           # ex: "preco_com_reforma12"
  pj_short <- sub("^ln_", "", colnames(Giv))
  fit2$coef$gamma[, match(pj_short, cols_fit)] <- Giv
  elas_quaids_manual(
    fit2, x = x_star, p = p_star,
    normalize_eval = normalize_eval, normalize_deriv = normalize_deriv,
    dlogp = 1e-3, enforce_hicks = TRUE, enforce_symmetry = TRUE
  )
}

## (gradiente central — não usado, deixo para debug se quiser)
numgrad_central <- function(fun, theta, eps = 1e-6){
  f0 <- fun(theta)
  m  <- length(f0); k <- length(theta)
  J  <- matrix(NA_real_, nrow = m, ncol = k)
  for (j in seq_len(k)) {
    thp <- theta; thm <- theta
    h   <- eps * (abs(theta[j]) + 1)
    thp[j] <- theta[j] + h
    thm[j] <- theta[j] - h
    fp <- fun(thp); fm <- fun(thm)
    J[, j] <- (fp - fm) / (2*h)
  }
  J
}

## ------------ 1) Objetos do pipeline ------------
drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))   # todos ln preços (exceto dropado)
eqs        <- setdiff(shareNames, omit_share)                  # shares estimadas
exogs      <- intersect(c("z","z2"), names(dfJI))              # exógenos no RHS (ajuste se precisar)

## IVs just-ID (precisa 1 PC por preço)
stopifnot(all(iv_star %in% names(dfJI)),
          length(iv_star) == length(endogs_all))

## ------------ 2) Monta e estima o sistema 3SLS just-ID (labels seguros) ------------
eqs_orig <- setdiff(shareNames, omit_share)
eqs_safe <- paste0("eq", seq_along(eqs_orig))                  # labels sem espaço/underscore
orig2safe <- setNames(eqs_safe, eqs_orig)
safe2orig <- setNames(eqs_orig, eqs_safe)

rhs_y <- paste(c(endogs_all, exogs), collapse = " + ")
rhs_z <- paste(c(exogs, iv_star),  collapse = " + ")

eq_list   <- setNames(vector("list", length(eqs_orig)), eqs_safe)
inst_list <- setNames(vector("list", length(eqs_orig)), eqs_safe)

for (i in seq_along(eqs_orig)) {
  y_orig  <- eqs_orig[i]
  y_safe  <- orig2safe[[y_orig]]
  eq_list[[y_safe]]   <- stats::as.formula(paste(y_orig, "~", rhs_y))
  inst_list[[y_safe]] <- stats::as.formula(paste("~", rhs_z))
}

use_all <- unique(c(eqs_orig, endogs_all, exogs, iv_star))
dat_sys <- dfJI[stats::complete.cases(dfJI[, use_all, drop = FALSE]), , drop = FALSE]

fit_sys <- systemfit::systemfit(eq_list, method = "3SLS", inst = inst_list, data = dat_sys)

## ------------ 3) Mapeamento robusto eq-systemfit ↔ shares ------------
## --- PATCH 1: diagnóstico robusto de labels/LHS no systemfit ---
diag_eq_map <- function(fit_sys){
  # 1) Tenta pegar via fit_sys$eq (caminho "bonito")
  eq_list <- tryCatch(fit_sys$eq, error = function(e) NULL)
  if (!is.null(eq_list) && length(eq_list) > 0) {
    labs <- names(eq_list)
    if (is.null(labs) || !any(nzchar(labs))) labs <- paste0("eq", seq_along(eq_list))
    lhs  <- vapply(seq_along(eq_list), function(i){
      f <- try(eq_list[[i]]$formula, silent = TRUE)
      if (inherits(f, "try-error") || is.null(f)) f <- try(eq_list[[i]]$call$formula, silent = TRUE)
      if (inherits(f, "try-error") || is.null(f)) return(NA_character_)
      all.vars(f)[1]
    }, character(1))
    return(data.frame(eq_label = labs,
                      lhs      = lhs,
                      key_label= toupper(gsub("[^A-Za-z0-9]","", labs)),
                      key_lhs  = toupper(gsub("[^A-Za-z0-9]","", lhs)),
                      stringsAsFactors = FALSE))
  }
  
  # 2) Fallback: reconstrói labels a partir de names(coef(fit_sys))
  cn <- names(coef(fit_sys))
  if (length(cn) == 0) {
    return(data.frame(eq_label=character(0), lhs=character(0),
                      key_label=character(0), key_lhs=character(0)))
  }
  # padrão eq<lab>_<var> ou eq<lab>.<var>
  eq_guess <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  data.frame(eq_label = eq_labels,
             lhs       = NA_character_,                  # sem acesso fácil ao LHS aqui
             key_label = toupper(gsub("[^A-Za-z0-9]","", eq_labels)),
             key_lhs   = NA_character_,
             stringsAsFactors = FALSE)
}

.canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))


## --- PATCH 2: mapeamento super-robusto (com fallback por ordem) ---
auto_eq_map2 <- function(fit_sys, shareNames, omit_share){
  diag <- diag_eq_map(fit_sys)
  shares_ok  <- setdiff(shareNames, omit_share)
  key_shares <- .canon(shares_ok)
  
  # Tentativa A/B/C como antes (casa por LHS literal / LHS canônico / LABEL canônico)
  hitA <- if (!all(is.na(diag$lhs))) match(diag$lhs, shares_ok,  nomatch = 0L) else rep(0L, nrow(diag))
  hitB <- if (!all(is.na(diag$key_lhs))) match(diag$key_lhs, key_shares, nomatch = 0L) else rep(0L, nrow(diag))
  hitC <- match(diag$key_label, key_shares, nomatch = 0L)
  
  idx <- unique(c(which(hitA>0), which(hitB>0), which(hitC>0)))
  
  if (length(idx) > 0) {
    eqs_safe <- diag$eq_label[idx]
    eqs_orig <- shares_ok[
      ifelse(idx %in% which(hitA>0), hitA[idx],
             ifelse(idx %in% which(hitB>0), hitB[idx], hitC[idx]))
    ]
    dup <- duplicated(eqs_orig)
    if (any(dup)) { eqs_safe <- eqs_safe[!dup]; eqs_orig <- eqs_orig[!dup] }
    return(list(eqs_safe = eqs_safe, eqs_orig = eqs_orig, diag = diag))
  }
  
  # Fallback FINAL: casa por ORDEM (eqs como aparecem nos coeficientes) ↔ shares_ok
  cn <- names(coef(fit_sys))
  if (length(cn) == 0)
    stop("auto_eq_map2: não foi possível inspecionar coeficientes do systemfit.")
  eq_guess <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  Ke <- length(eq_labels)
  if (Ke != length(shares_ok)) {
    stop(sprintf("auto_eq_map2: fallback por ordem falhou: #eq=%d != #shares=%d.",
                 Ke, length(shares_ok)))
  }
  message("[auto_eq_map2] Usando fallback por ORDEM: ",
          paste(sprintf("%s↔%s", eq_labels, shares_ok), collapse = "; "))
  list(eqs_safe = eq_labels, eqs_orig = shares_ok,
       diag = data.frame(eq_label = eq_labels, lhs = NA_character_,
                         key_label = .canon(eq_labels), key_lhs = NA_character_))
}


## ------------ 4) Extrator θ=vec(Γ) e Var(θ) (robusto a nomes) ------------
stack_gamma_vcov_sys2 <- function(fit_sys, pj, eq_safe, eq_orig) {
  stopifnot(length(eq_safe) == length(eq_orig))
  co <- coef(fit_sys); V <- vcov(fit_sys); cn <- names(co)
  Kp <- length(pj); Ke <- length(eq_safe); K <- Ke * Kp
  theta <- numeric(K); Vth <- matrix(NA_real_, K, K); nm_th <- character(K)
  
  find_idx <- function(eq, var) {
    pats <- c(paste0("^", eq, "_", var, "$"),
              paste0("^", eq, "\\.", var, "$"))
    for (pt in pats) {
      ii <- grep(pt, cn); if (length(ii)==1L) return(ii)
    }
    ii_eq  <- grep(paste0("^", eq, "(_|\\.)"), cn)
    ii_var <- grep(paste0("(_|\\.)", var, "$"), cn)
    ii <- intersect(ii_eq, ii_var)
    if (length(ii)==1L) return(ii)
    key_eq  <- .canon(eq); key_var <- .canon(var); key_cn <- .canon(cn)
    cand <- which(startsWith(key_cn, key_eq) & endsWith(key_cn, key_var))
    if (length(cand)==1L) return(cand)
    NA_integer_
  }
  
  for (i in seq_along(eq_safe)) for (j in seq_along(pj)) {
    pos <- (i - 1L) * Kp + j
    idx <- find_idx(eq_safe[i], pj[j])
    if (is.na(idx))
      stop(sprintf("Coeficiente não encontrado: eq='%s' var='%s'\nVeja names(coef(fit_sys)).",
                   eq_safe[i], pj[j]))
    theta[pos] <- unname(co[idx]); nm_th[pos] <- paste(eq_orig[i], pj[j], sep="::")
  }
  
  for (i in seq_along(eq_safe)) for (k in seq_along(eq_safe)) {
    for (a in seq_along(pj)) for (b in seq_along(pj)) {
      pos1 <- (i - 1L) * Kp + a; pos2 <- (k - 1L) * Kp + b
      if (is.na(Vth[pos1,pos2])) {
        idx1 <- find_idx(eq_safe[i], pj[a]); idx2 <- find_idx(eq_safe[k], pj[b])
        Vth[pos1,pos2] <- V[idx1, idx2]
      }
    }
  }
  names(theta) <- nm_th
  list(theta = theta, V = Vth)
}

## ------------ 5) Delta Method para ICs (sem bootstrap) ------------
## ================= DELTA: ICs para Hicks + Marshall + Eta =================
## ================= DELTA: ICs p/ Hicks + Marshall + Eta (robusto a nomes) ==============
# ==================== Delta CIs v7 (com Marshall FULL) ====================
elasticity_CIs_delta_system_v7 <- function(
    fit_sys, fit_template, dat_sys,
    priceNames, shareNames, drop_price, omit_share,
    level = 0.95, method = c("pointwise","bonferroni")
){
  if (!requireNamespace("numDeriv", quietly = TRUE))
    stop("Instale: install.packages('numDeriv')")
  if (!requireNamespace("tibble", quietly = TRUE))
    stop("Instale: install.packages('tibble')")
  method <- match.arg(method)
  `%||%` <- function(a,b) if(!is.null(a)) a else b
  .canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x %||% ""))
  
  ## --- 1) Mapear equações do systemfit -> suas shares (com fallback) ---
  map <- auto_eq_map2(fit_sys, shareNames, omit_share)  # sua helper
  eqs_safe <- map$eqs_safe
  eqs_orig <- map$eqs_orig
  
  ## --- 2) Empilhar theta = vec(Gamma) e Var(theta) ---
  pj       <- paste0("ln_", setdiff(priceNames, drop_price))
  pj_short <- sub("^ln_", "", pj)
  meta     <- stack_gamma_vcov_sys2(fit_sys, pj, eqs_safe, eqs_orig) # sua helper
  theta    <- meta$theta
  Vth      <- meta$V
  
  ## --- 3) De theta -> Gfull (com adding-up para a share omitida) ---
  theta_to_G <- function(th){
    Ke <- length(eqs_orig); Kp <- length(pj)
    Ghat <- matrix(th, nrow = Ke, byrow = TRUE, dimnames = list(eqs_orig, pj))
    Gfull <- matrix(NA_real_, nrow = length(shareNames), ncol = Kp,
                    dimnames = list(shareNames, pj))
    Gfull[eqs_orig, pj] <- Ghat
    Gfull[omit_share, pj] <- -colSums(Ghat, na.rm = TRUE)
    Gfull
  }
  
  ## --- 4) Ponto de avaliação (x*, p*) ---
  x_star <- if ("gasto_total_atualhat" %in% names(dat_sys))
    median(dat_sys$gasto_total_atualhat, na.rm=TRUE) else
      median(dat_sys[[grep("gasto|exp", names(dat_sys), value=TRUE)[1]]], na.rm=TRUE)
  p_star <- exp(colMeans(dat_sys[paste0("ln_", priceNames)], na.rm=TRUE))
  
  ## --- 5) Avaliação em theta para capturar ordem/nome “nativos” ---
  G0 <- theta_to_G(theta)
  f2 <- fit_template
  cols_fit <- colnames(f2$coef$gamma)
  f2$coef$gamma[, match(pj_short, cols_fit)] <- G0
  E0 <- elas_quaids_manual(
    f2, x=x_star, p=p_star,
    normalize_eval="softmax", normalize_deriv="softmax",
    dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
  )
  
  rn0 <- rownames(E0$hicks); cn0 <- colnames(E0$hicks)
  if (is.null(rn0)) rn0 <- shareNames
  if (is.null(cn0)) cn0 <- shareNames
  k0  <- length(rn0)
  
  ## --- 6) Função g(theta) com vetor na ordem: Hicks FULL, Marshall FULL, Eta ---
  ##     IMPORTANTe: vetorizar por LINHA: as.vector(t(M)) para reconstruir com byrow=TRUE
  g_eval <- function(th){
    G <- theta_to_G(th)
    f3 <- fit_template
    f3$coef$gamma[, match(pj_short, cols_fit)] <- G
    E <- elas_quaids_manual(
      f3, x=x_star, p=p_star,
      normalize_eval="softmax", normalize_deriv="softmax",
      dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
    )
    c(as.vector(t(E$hicks)),
      as.vector(t(E$marshall)),
      as.numeric(E$expenditure))
  }
  
  g0 <- g_eval(theta)
  J  <- numDeriv::jacobian(g_eval, theta)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  lo <- g0 - z*se; hi <- g0 + z*se
  
  ## --- 7) Reconstituir blocos NATIVOS (sem alinhamento) ---
  idx_H    <- seq_len(k0*k0)
  idx_Mall <- (k0*k0) + seq_len(k0*k0)
  idx_eta  <- (k0*k0 + k0*k0) + seq_len(k0)
  
  H_hat0 <- matrix(g0[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_lo0  <- matrix(lo[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_hi0  <- matrix(hi[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_se0  <- matrix(se[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  M_hat0 <- matrix(g0[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_lo0  <- matrix(lo[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_hi0  <- matrix(hi[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_se0  <- matrix(se[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  eta_hat0 <- setNames(g0[idx_eta], rn0)
  eta_lo0  <- setNames(lo[idx_eta], rn0)
  eta_hi0  <- setNames(hi[idx_eta], rn0)
  eta_se0  <- setNames(se[idx_eta], rn0)
  
  ## --- 8) Alinhamento p/ shareNames (se houver match) ---
  key_tar <- .canon(shareNames)
  key_r   <- .canon(rownames(H_hat0))
  key_c   <- .canon(colnames(H_hat0))
  ridx <- match(key_tar, key_r)
  cidx <- match(key_tar, key_c)
  
  all_miss <- all(is.na(ridx)) && all(is.na(cidx))
  k <- length(shareNames)
  
  mkMat <- function() matrix(NA_real_, nrow=k, ncol=k, dimnames=list(shareNames, shareNames))
  H_hat <- mkMat(); H_lo <- mkMat(); H_hi <- mkMat(); H_se <- mkMat()
  M_hat <- mkMat(); M_lo <- mkMat(); M_hi <- mkMat(); M_se <- mkMat()
  
  ok_r <- which(!is.na(ridx)); ok_c <- which(!is.na(cidx))
  if (length(ok_r) && length(ok_c)) {
    H_hat[ok_r, ok_c] <- H_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_lo [ok_r, ok_c] <- H_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_hi [ok_r, ok_c] <- H_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_se [ok_r, ok_c] <- H_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    
    M_hat[ok_r, ok_c] <- M_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_lo [ok_r, ok_c] <- M_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_hi [ok_r, ok_c] <- M_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_se [ok_r, ok_c] <- M_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
  }
  
  map1 <- match(key_tar, key_r)
  Mdiag_hat <- Mdiag_lo <- Mdiag_hi <- Mdiag_se <- rep(NA_real_, k)
  names(Mdiag_hat) <- names(Mdiag_lo) <- names(Mdiag_hi) <- names(Mdiag_se) <- shareNames
  ok <- which(!is.na(map1))
  if (length(ok)) {
    Mdiag_hat[ok] <- diag(M_hat0)[ map1[ok] ]
    Mdiag_lo [ok] <- diag(M_lo0)[ map1[ok] ]
    Mdiag_hi [ok] <- diag(M_hi0)[ map1[ok] ]
    Mdiag_se [ok] <- diag(M_se0)[ map1[ok] ]
  }
  
  eta_hat <- eta_lo <- eta_hi <- eta_se <- rep(NA_real_, k)
  names(eta_hat) <- names(eta_lo) <- names(eta_hi) <- names(eta_se) <- shareNames
  if (length(ok)) {
    eta_hat[ok] <- eta_hat0[ map1[ok] ]
    eta_lo [ok] <- eta_lo0 [ map1[ok] ]
    eta_hi [ok] <- eta_hi0 [ map1[ok] ]
    eta_se [ok] <- eta_se0 [ map1[ok] ]
  }
  
  ## --- 9) Saídas em formato TIDY (alinhado e nativo) ---
  to_tbl <- function(M, lo, hi, se, rows, cols){
    do.call(rbind, lapply(seq_along(rows), function(i)
      data.frame(i=rows[i], j=cols, hat=M[i,], lo=lo[i,], hi=hi[i,], se=se[i,], row.names=NULL)))
  }
  
  hicks_tbl_aligned    <- to_tbl(H_hat, H_lo, H_hi, H_se, shareNames, shareNames)
  hicks_tbl_native     <- to_tbl(H_hat0, H_lo0, H_hi0, H_se0, rn0, cn0)
  marshall_tbl_aligned <- to_tbl(M_hat, M_lo, M_hi, M_se, shareNames, shareNames)
  marshall_tbl_native  <- to_tbl(M_hat0, M_lo0, M_hi0, M_se0, rn0, cn0)
  
  marshall_diag_tbl <- tibble::tibble(
    good = shareNames, M_hat = Mdiag_hat, M_lo = Mdiag_lo, M_hi = Mdiag_hi, M_se = Mdiag_se
  )
  eta_tbl <- tibble::tibble(
    good = shareNames, eta_hat = eta_hat, eta_lo = eta_lo, eta_hi = eta_hi, eta_se = eta_se
  )
  hicks_diag_tbl <- tibble::tibble(
    good = shareNames,
    H_hat = diag(H_hat), H_lo = diag(H_lo), H_hi = diag(H_hi), H_se = diag(H_se)
  )
  
  note <- NULL
  if (all_miss) {
    note <- "Nomes não casaram: use *_native (ordem e rótulos nativos do E)."
    message("[Delta-elasticities v7] ", note)
  }
  
  list(
    ## Hicks
    hicks_full        = hicks_tbl_aligned,
    hicks_full_native = hicks_tbl_native,
    hicks_diag        = hicks_diag_tbl,
    ## Marshall (FULL + DIAG) — agora com ICs completos
    marshall_full        = marshall_tbl_aligned,
    marshall_full_native = marshall_tbl_native,
    marshall_diag        = marshall_diag_tbl,
    ## Expenditure
    eta              = eta_tbl,
    level            = level,
    method           = method,
    mapping_diag     = map$diag,
    note             = note
  )
}
# ================== FIM Delta CIs v7 (com Marshall FULL) ==================

## ======================= FIM v5 =======================

## ========================= FIM DO PATCH =========================

## ------------ 6) Rodar e imprimir ------------
# 6.1 Diagnóstico de como o systemfit rotulou as equações
cat("\n[systemfit] labels e LHS por equação:\n")
print(diag_eq_map(fit_sys), row.names = FALSE)

# 6.2 ICs por Delta (sem bootstrap)
cis_sys <- elasticity_CIs_delta_system_v7(
  fit_sys      = fit_sys,
  fit_template = best_fit,
  dat_sys      = dat_sys,
  priceNames   = priceNames,
  shareNames   = shareNames,
  drop_price   = drop_price,
  omit_share   = omit_share,
  level        = 0.95,
  method       = "bonferroni"  # ou "pointwise"
)

cat("\n[ICs Delta — Marshall]\n"); print(cis_sys$marshall_full_native, digits=3)
cat("\n[ICs Delta — Hicks]\n"); print(cis_sys$hicks_full_native, digits = 3)
cat("\n[ICs Delta — Eta]\n");             print(cis_sys$eta,          digits=3)

# D = endógenos (todos ln preços), X = exógenos, Z = IVs excluídos (PCs)
sw_F_one <- function(j, Dnames, Xnames, Znames, data){
  Dminus <- setdiff(Dnames, j)
  # 1ª etapa para cada D_-j
  fits <- lapply(Dminus, function(d)
    lm(reformulate(c(Xnames, Znames), response=d), data=data))
  Dhat <- as.data.frame(sapply(fits, fitted))
  colnames(Dhat) <- paste0(Dminus, "_hat")
  dat2 <- cbind(data, Dhat)
  
  # Regressão “auxiliar” de D_j em X + Z + D_-j_hat
  f_red  <- lm(reformulate(c(Xnames,        names(Dhat)), response=j), data=dat2)
  f_full <- lm(reformulate(c(Xnames, Znames, names(Dhat)), response=j), data=dat2)
  as.numeric(anova(f_red, f_full)$F[2])  # este é o SW F_j
}

SW_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  use <- unique(c(endogs_all, exogs, iv_star))
  dat <- dfJI[complete.cases(dfJI[, use, drop=FALSE]), , drop=FALSE]
  tibble::tibble(price=Dv, n=nrow(dat),
                 F_SW_j = sw_F_one(Dv, endogs_all, exogs, iv_star, dat))
})
print(SW_tbl, digits=3)
write.csv2(SW_tbl, 'C:/Users/x16610962/Downloads/SW_tbl.csv', row.names = FALSE)
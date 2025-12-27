# =========================================================
# (PREP) Dependências e checagens
# =========================================================
suppressPackageStartupMessages({
  library(AER)
  library(dplyr)
  library(purrr)
  library(tibble)
})

# Se não houver split 80/20 feito no passo (5), cria agora
if (!exists("trn") || !exists("tst")) {
  set.seed(123)
  n   <- nrow(dat)
  ix  <- sample.int(n, size = floor(0.8 * n))
  trn <- dat[ix, ]
  tst <- dat[-ix, ]
}

# Se não houver fits_iv_trn e 3SLS_trn, reestima como no passo (5)
if (!exists("fits_iv_trn")) {
  fits_iv_trn <- map(share_vars, ~ fit_iv_weighted(
    dat = trn, dep = .x, price_vars = price_vars, weights_col = weights_col,
    iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
    z_var = z_var, z2_var = z2_var,
    seasonal_controls = seasonal_controls_full,
    use_spline = FALSE, spline_df = 4
  ))
  names(fits_iv_trn) <- share_vars
}

if (!exists("fit_3sls_trn")) {
  fit_3sls_trn <- fit_3sls_system(
    dat = trn, dep_vec = share_vars, price_vars = price_vars,
    iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
    z_var = z_var, z2_var = z2_var,
    seasonal_controls = seasonal_controls_min,
    use_spline = FALSE, spline_df = 4,
    drop_equation = "w_despesahat6",
    weights_col = weights_col,
    z2_present = FALSE
  )
}

# =========================================================
# (6) HÍBRIDO: normalização + escolha por share (2SLS vs 3SLS_minA)
# =========================================================

# Previsões OOS (se não existirem do passo 5, cria)
if (!exists("pred_iv_tst")) {
  pred_iv_tst <- sapply(share_vars, function(v) {
    as.numeric(predict(fits_iv_trn[[v]]$fit, newdata = tst))
  })
  colnames(pred_iv_tst) <- share_vars
}
if (!exists("pred_3sls_tst")) {
  pred_3sls_tst <- predict_3sls(fit_3sls_trn, tst)
}

w_tst   <- tst[[weights_col]]
obs_tst <- as.data.frame(tst[, share_vars])

rmse_w <- function(obs, pred, w) sqrt(weighted.mean((obs - pred)^2, w))

# Normalizador simplex (garante soma=1, com clip opcional para [0,1])
renorm_simplex <- function(M, clip01 = TRUE) {
  if (clip01) M <- pmin(pmax(M, 0), 1)
  s <- rowSums(M); s[s == 0] <- 1
  M / s
}

# Monta matriz 3SLS (lista -> matriz)
pred_3sls_tst_mat <- do.call(cbind, pred_3sls_tst[share_vars])
colnames(pred_3sls_tst_mat) <- share_vars

# Normaliza
pred_iv_tst_norm   <- renorm_simplex(pred_iv_tst, clip01 = TRUE)
pred_3sls_tst_norm <- renorm_simplex(pred_3sls_tst_mat, clip01 = TRUE)

# RMSE OOS (normalizado)
rmse_tab_norm <- tibble(
  share = share_vars,
  rmse_w_2sls_oos_norm = sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_iv_tst_norm[, v],   w_tst)),
  rmse_w_3sls_oos_norm = sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_3sls_tst_norm[, v], w_tst))
)

# Híbrido: escolhe o melhor por share e renormaliza
winners <- ifelse(rmse_tab_norm$rmse_w_3sls_oos_norm < rmse_tab_norm$rmse_w_2sls_oos_norm,
                  "3SLS_minA", "2SLS")

pred_hybrid <- pred_iv_tst_norm
for (j in seq_along(share_vars)) {
  v <- share_vars[j]
  if (winners[j] == "3SLS_minA") pred_hybrid[, v] <- pred_3sls_tst_norm[, v]
}
pred_hybrid <- renorm_simplex(pred_hybrid, clip01 = TRUE)

rmse_hybrid <- sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_hybrid[, v], w_tst))

cat("\n==== (6) HÍBRIDO (OOS, RMSE ponderado) ====\n")
print(rmse_tab_norm)
cat("\nVencedores por share (2SLS vs 3SLS_minA):\n")
print(tibble(share = share_vars, winner = winners))
cat("\nRMSE do híbrido:\n")
print(tibble(share = share_vars, rmse_w_hybrid_oos = as.numeric(rmse_hybrid)))


# =========================================================
# (7) C-TESTS POR BLOCO DE INSTRUMENTOS (diferença de Hansen-J robusto)
#     Requer utilitários: .get_mats(), .w2sls(), .hansenJ_from_res()
# =========================================================
# --- util: remove colunas de variância ~0 e duplicadas ---
## ===== UTILITÁRIOS ROBUSTOS (substituem os anteriores) =====================

# extrai y, X (regressores de 2ª etapa), Zfull (exógenas + IVs) e pesos usados pelo ivreg
.get_mats <- function(fit){
  mf <- model.frame(fit)                      # usa exatamente as linhas do modelo
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Zf <- model.matrix(fit, component = "instruments")
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, NROW(X))
  X  <- as.matrix(X);  storage.mode(X)  <- "double"
  Zf <- as.matrix(Zf); storage.mode(Zf) <- "double"
  list(y=y, X=X, Zfull=Zf, w=as.numeric(w))
}

# remove colunas com variância ~0 e duplicadas
.drop_nzv <- function(A, tol = 1e-12){
  if (is.null(A) || NCOL(A) == 0) return(A[, 0, drop = FALSE])
  A <- as.matrix(A)
  keep <- apply(A, 2, function(x) sd(x) > tol)
  A <- A[, keep, drop = FALSE]
  A <- A[, !duplicated(colnames(A)), drop = FALSE]
  storage.mode(A) <- "double"
  A
}

# 2SLS ponderado por QR (evita (Z'Z)^{-1}); pré-multiplica por sqrt(w)
.w2sls_qr <- function(y, X, Z, w, ridge = 1e-8){
  stopifnot(length(y) == NROW(X), NROW(X) == NROW(Z), length(w) == NROW(X))
  s  <- sqrt(w / mean(w))
  yS <- y * s; XS <- X * s; ZS <- Z * s
  
  qZ <- qr(ZS)                            # projeção no espaço de Z
  PZ_X <- qr.fitted(qZ, XS)
  PZ_y <- qr.fitted(qZ, yS)
  
  XtPZ_X <- crossprod(XS, PZ_X)
  XtPZ_y <- crossprod(XS, PZ_y)
  
  if (!is.finite(rcond(XtPZ_X)) || rcond(XtPZ_X) < 1e-10){
    diag(XtPZ_X) <- diag(XtPZ_X) + ridge
  }
  
  beta <- solve(XtPZ_X, XtPZ_y)
  e    <- as.numeric(yS - XS %*% beta)
  list(beta = beta, e = e, X = XS, Z = ZS, n = NROW(XS), p = ncol(XS), k = ncol(ZS))
}

# Hansen-J robusto (HC0) a partir de resíduos e ZS
.hansenJ_from_res <- function(e, ZS, p, ridge_seq = c(0, 1e-8, 1e-6, 1e-4)){
  n <- NROW(ZS)
  m <- crossprod(ZS, e) / n
  S <- crossprod(ZS * e) / n; S <- (S + t(S))/2
  note <- NA_character_
  ok <- FALSE; sol <- NULL
  for (rg in ridge_seq){
    Sreg <- S
    if (rg > 0) diag(Sreg) <- diag(Sreg) + rg
    sol <- try(qr.solve(Sreg, m), silent = TRUE)
    if (!inherits(sol, "try-error") && all(is.finite(sol))) { ok <- TRUE; note <- if (rg>0) paste0("ridge=", rg) else NA_character_; break }
  }
  if (!ok) return(tibble::tibble(J = NA_real_, df = NA_real_, p = NA_real_, note = "S_singular"))
  J  <- drop(n * crossprod(m, sol))
  df <- ncol(ZS) - p
  tibble::tibble(
    J = if (df>0) J else NA_real_,
    df = if (df>0) df else 0L,
    p  = if (df>0) pchisq(J, df, lower.tail = FALSE) else NA_real_,
    note = note
  )
}

# util p/ escapar regex do nome do z_var (se necessário)
.esc <- function(x) gsub("([\\W])", "\\\\\\1", x, perl = TRUE)

## ===== C-TEST ROBUSTO (diferença de Hansen-J) ==============================

ctest_block <- function(fit, regex_block, ridge_seq = c(0, 1e-8, 1e-6, 1e-4)) {
  tryCatch({
    mats <- .get_mats(fit)
    X  <- mats$X; Zf <- mats$Zfull; w <- mats$w
    Xn <- colnames(X)
    
    # IVs estritos (os que não estão em X)
    Ziv <- Zf[, setdiff(colnames(Zf), Xn), drop = FALSE]
    Ziv <- .drop_nzv(Ziv)
    
    # constrói Z_full e Z_restr (dropar bloco)
    drop_cols <- if (NCOL(Ziv)) grep(regex_block, colnames(Ziv), value = TRUE) else character(0)
    Z_full  <- cbind(X, Ziv)
    Z_restr <- cbind(X, if (length(drop_cols)) Ziv[, setdiff(colnames(Ziv), drop_cols), drop = FALSE] else Ziv)
    
    storage.mode(Z_full)  <- "double"
    storage.mode(Z_restr) <- "double"
    
    # Se não há IVs além de X:
    if (ncol(Z_full) <= ncol(X))
      return(tibble::tibble(J_full=NA_real_, J_rest=NA_real_, C=NA_real_, df=0L, p=NA_real_,
                            note_full="no_iv", note_rest="no_iv",
                            k_full=ncol(Z_full), k_restr=ncol(Z_restr),
                            dropped=paste(drop_cols, collapse=",")))
    
    est_full <- .w2sls_qr(mats$y, X, Z_full,  w)
    J_full   <- .hansenJ_from_res(est_full$e, est_full$Z, est_full$p, ridge_seq)
    
    # Se Z_restr removeu todos os IVs (fica só X), não há C-test
    if (ncol(Z_restr) <= ncol(X)){
      return(tibble::tibble(J_full=J_full$J, J_rest=NA_real_, C=NA_real_,
                            df=ncol(Z_full)-ncol(X), p=NA_real_,
                            note_full=J_full$note, note_rest="no_iv",
                            k_full=ncol(Z_full), k_restr=ncol(Z_restr),
                            dropped=paste(drop_cols, collapse=",")))
    }
    
    est_rest <- .w2sls_qr(mats$y, X, Z_restr, w)
    J_rest   <- .hansenJ_from_res(est_rest$e, est_rest$Z, est_rest$p, ridge_seq)
    
    df_C <- ncol(Z_full) - ncol(Z_restr)
    C    <- (J_full$J - J_rest$J)
    pC   <- if (is.finite(C) && df_C > 0) pchisq(C, df = df_C, lower.tail = FALSE) else NA_real_
    
    tibble::tibble(
      J_full = J_full$J, J_rest = J_rest$J, C = C, df = df_C, p = pC,
      note_full = J_full$note, note_rest = J_rest$note,
      k_full = ncol(Z_full), k_restr = ncol(Z_restr),
      dropped = paste(drop_cols, collapse = ",")
    )
  }, error = function(e){
    tibble::tibble(J_full=NA_real_, J_rest=NA_real_, C=NA_real_, df=NA_integer_, p=NA_real_,
                   note_full=paste0("err: ", conditionMessage(e)), note_rest="err",
                   k_full=NA_integer_, k_restr=NA_integer_, dropped=NA_character_)
  })
}

## ===== BLOCO: rodar C-tests por share x regex ==============================

# Defina seus blocos; se quiser testar o próprio z_var nos IVs, use .esc(z_var)
blocks_regex <- c(
  "^iv_op",                 # op*
  "^iv_hg",                 # hg*
  "_x_",                    # interações
  "^(iv_sin|iv_cos)",       # sazonais trig
  paste0("^", .esc(z_var), "$")  # z_var, se estiver entre IVs
)

ctests_tbl <- purrr::map_dfr(share_vars, function(dep) {
  fit <- fits_iv[[dep]]$fit
  purrr::map_dfr(blocks_regex, function(rx) ctest_block(fit, rx)) |>
    dplyr::mutate(share = dep, bloco = blocks_regex, .before = 1)
})

print(ctests_tbl)
# Dica de leitura:
# p pequeno (ex.: <0.05) no bloco ⇒ forte evidência de que esse bloco não é exógeno.
# Remover esse bloco e reestimar costuma melhorar J-tests/OOS.


# =========================================================
# (8) LIML e FULLER-k (ponderados) sem pacotes externos
#     Implementação via problema de autovalor generalizado
# =========================================================

# Projeções e matrizes auxiliares
.projP <- function(A, ridge = 0) {
  AtA <- crossprod(A)
  if (!is.finite(rcond(AtA)) || rcond(AtA) < 1e-12) diag(AtA) <- diag(AtA) + ridge * mean(diag(AtA))
  A %*% solve(AtA, t(A))
}
.projM <- function(A, ridge = 0) diag(nrow(A)) - .projP(A, ridge = ridge)

# Estimador LIML/Fuller ponderado
# fuller_a = NULL -> LIML; fuller_a = 1 ou 4 -> Fuller(1)/Fuller(4)
.liml_fuller_w <- function(y, X, Z, w, fuller_a = NULL, ridge = 1e-8) {
  # escala por sqrt(w)
  s  <- sqrt(w / mean(w))
  yS <- y * s; XS <- X * s; ZS <- Z * s
  n  <- NROW(XS)
  
  # Projeções
  MZ <- .projM(ZS, ridge = ridge)
  MX <- .projM(XS, ridge = ridge)
  
  # Monta W = [y, X]
  W  <- cbind(yS, XS)
  
  # A = W' MZ W ; B = W' MX W
  A <- crossprod(W, MZ %*% W); A <- (A + t(A))/2
  B <- crossprod(W, MX %*% W); B <- (B + t(B))/2
  if (!is.finite(rcond(B)) || rcond(B) < 1e-12) diag(B) <- diag(B) + ridge * mean(diag(B))
  
  # autovalor generalizado: eigen(B^{-1} A)
  ev  <- eigen(solve(B, A), symmetric = TRUE, only.values = TRUE)$values
  kappa <- min(Re(ev))  # menor autovalor
  
  # Ajuste Fuller (opcional): kappa* = kappa - a/(n - L), L = ncol(Z)
  if (!is.null(fuller_a)) {
    L <- ncol(ZS)
    kappa <- kappa - fuller_a / (n - L)
  }
  
  S   <- MZ - kappa * MX
  XtS <- crossprod(XS, S)
  XtSX <- XtS %*% XS
  XtSy <- XtS %*% yS
  
  if (!is.finite(rcond(XtSX)) || rcond(XtSX) < 1e-12) diag(XtSX) <- diag(XtSX) + ridge * mean(diag(XtSX))
  
  beta <- drop(solve(XtSX, XtSy))
  e    <- drop(yS - XS %*% beta)
  
  list(beta = beta, e = e, n = n, p = ncol(XS), kappa = kappa)
}

# Wrapper: estima a partir de um ivreg (AER) já montado
fit_liml_from_ivreg <- function(fit, fuller_a = NULL, ridge = 1e-8) {
  mf <- model.frame(fit)
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(fit, component = "regressors")
  Z  <- model.matrix(fit, component = "instruments")
  w  <- model.weights(mf); if (is.null(w)) w <- rep(1, NROW(X))
  
  est <- .liml_fuller_w(y, X, Z, w, fuller_a = fuller_a, ridge = ridge)
  names(est$beta) <- colnames(X)
  est
}

# Predição com coeficientes (usa matriz de regressoras do 'fit')
# ---------- PATCH: predição LIML/Fuller alinhada ao newdata ----------
# Usa os 'terms' e 'contrasts' do X do treino para criar X_new consistente.
predict_liml_from_fit <- function(fit, beta, newdata) {
  X_tr <- model.matrix(fit, component = "regressors") # do treino
  tt   <- attr(X_tr, "terms")
  cc   <- attr(X_tr, "contrasts")
  
  # constrói X_new com a mesma especificação; mantém NAs (vamos tratar no RMSE)
  X_new <- model.matrix(tt, data = newdata, contrasts.arg = cc, na.action = na.pass)
  
  # garante que todas as colunas de beta existam em X_new (insere 0 se faltar)
  miss <- setdiff(names(beta), colnames(X_new))
  if (length(miss)) {
    X_new <- cbind(X_new,
                   matrix(0, nrow(X_new), length(miss),
                          dimnames = list(NULL, miss)))
  }
  # reordena colunas
  X_new <- X_new[, names(beta), drop = FALSE]
  
  drop(X_new %*% beta)
}

predict_liml_from_fit <- function(fit_ivreg, beta, newdata) {
  if (is.null(names(beta)))
    stop("beta precisa ser um vetor nomeado com os mesmos nomes de colunas do X.")
  
  # fórmula só dos regressores (1º RHS do ivreg: y ~ X)
  f_reg <- formula(fit_ivreg, rhs = 1L)
  tt    <- terms(f_reg)                # inclui a resposta
  Xnew  <- model.matrix(delete.response(tt), data = newdata)
  
  # checagem/alinhamento dos nomes
  faltam <- setdiff(colnames(Xnew), names(beta))
  if (length(faltam))
    stop("beta não tem coeficientes para: ", paste(faltam, collapse = ", "))
  
  beta_use <- beta[colnames(Xnew)]     # descarta extras e ordena
  drop(Xnew %*% beta_use)
}

# Helper: estima LIML se disponível; senão, volta em 2SLS
fit_liml_one <- function(fit_ivreg, data) {
  f <- formula(fit_ivreg)  # reaproveita a MESMA fórmula de ivreg
  if (requireNamespace("ivreg", quietly = TRUE) &&
      "method" %in% names(formals(ivreg::ivreg))) {
    fitL <- ivreg::ivreg(f, data = data, method = "LIML")
    list(beta = coef(fitL), method = "LIML")
  } else {
    message("Pacote {ivreg} com method='LIML' não encontrado. Usando 2SLS por ora.")
    list(beta = coef(fit_ivreg), method = "2SLS-fallback")
  }
}

# ---------- LIML robusto (sem (Z'Z)^{-1}) ----------
liml_from_ivreg_fit <- function(fit, data, tol = 1e-10, ridge = 1e-8) {
  mf <- model.frame(fit, data = data)
  y  <- model.response(mf)
  X  <- model.matrix(fit, "regressors")
  Z  <- model.matrix(fit, "instruments")
  n  <- NROW(X)
  
  # pesos (se houver)
  w <- tryCatch(model.weights(mf), error = function(e) NULL)
  if (!is.null(w)) {
    sw <- sqrt(w)
    y  <- as.numeric(sw) * y
    X  <- X * sw
    Z  <- Z * sw
  }
  
  # Projetor via QR (sem (Z'Z)^{-1})
  proj_Q <- function(M, tol = 1e-10) {
    qrm <- qr(M, tol = tol)
    rnk <- qrm$rank
    if (rnk == 0) return(matrix(0, nrow = nrow(M), ncol = nrow(M)))
    Q <- qr.Q(qrm)   # n x r
    Q %*% t(Q)
  }
  
  Pz <- proj_Q(Z, tol);  Mz <- diag(n) - Pz
  Px <- proj_Q(X, tol);  Mx <- diag(n) - Px
  
  Qc <- cbind(y, X)
  A  <- crossprod(Qc, Mz %*% Qc)
  B  <- crossprod(Qc, Mx %*% Qc)
  
  # Simetriza e estabiliza B
  A  <- (A + t(A))/2
  B  <- (B + t(B))/2
  
  # eigen generalizado via B^{-1/2} A B^{-1/2} (evita solve(B, A))
  eigB <- eigen(B, symmetric = TRUE)
  pos  <- eigB$values > max(eigB$values) * tol
  if (!any(pos)) {
    lam <- ridge * mean(diag(B))
    eigB <- eigen(B + lam * diag(ncol(B)), symmetric = TRUE)
    pos  <- eigB$values > max(eigB$values) * tol
  }
  U   <- eigB$vectors[, pos, drop = FALSE]
  d   <- eigB$values [  pos]
  Binv2 <- U %*% diag(1/sqrt(d), nrow = length(d)) %*% t(U)
  
  C   <- Binv2 %*% A %*% Binv2
  k_hat <- min(Re(eigen(C, symmetric = TRUE, only.values = TRUE)$values))
  
  # k-class com k = k_hat
  Wk   <- diag(n) - k_hat * Mz
  XtWX <- crossprod(X, Wk %*% X)
  XtWy <- crossprod(X, Wk %*% y)
  
  # estabiliza XtWX se necessário
  if (qr(XtWX)$rank < ncol(XtWX)) {
    lam <- ridge * mean(diag(XtWX))
    XtWX <- XtWX + lam * diag(ncol(XtWX))
  }
  beta <- solve(XtWX, XtWy)
  
  list(beta = drop(beta),
       k = k_hat,
       terms_reg = terms(fit, "regressors"))
}

# ---------- Predição com os betas ----------
predict_with_betas <- function(fit, beta, newdata) {
  Xnew <- model.matrix(terms(fit, "regressors"), newdata)
  bn   <- names(beta)
  if (!all(bn %in% colnames(Xnew))) {
    stop("Faltam colunas em newdata: ", paste(setdiff(bn, colnames(Xnew)), collapse=", "))
  }
  as.numeric(Xnew[, bn, drop = FALSE] %*% beta)
}

# Mantém o nome que você já usa
predict_liml_from_fit <- function(fit, beta, newdata) {
  predict_with_betas(fit, beta, newdata)
}


# ---- Helpers LIML (k-class) a partir de um ajuste ivreg existente ----

## betas LIML no treino
liml_fits_trn <- lapply(share_vars, function(v) {
  fit_v <- fits_iv_trn[[v]]$fit
  liml_from_ivreg_fit(fit_v, data = trn)
})
names(liml_fits_trn) <- share_vars

## previsões OOS
pred_liml_tst <- sapply(share_vars, function(v) {
  fit_v  <- fits_iv_trn[[v]]$fit
  beta_v <- liml_fits_trn[[v]]$beta
  predict_liml_from_fit(fit_v, beta_v, tst)
})



## betas LIML no treino
liml_fits_trn <- lapply(share_vars, function(v) {
  fit_v <- fits_iv_trn[[v]]$fit
  liml_from_ivreg_fit(fit_v, data = trn)
})
names(liml_fits_trn) <- share_vars

## previsões OOS
pred_liml_tst <- sapply(share_vars, function(v) {
  fit_v  <- fits_iv_trn[[v]]$fit
  beta_v <- liml_fits_trn[[v]]$beta
  predict_liml_from_fit(fit_v, beta_v, tst)
})


# --- util: projetor por QR (não usa (Z'Z)^{-1}) ---
.proj_Q <- function(M, tol = 1e-10){
  qrm <- qr(M, tol = tol)
  r   <- qrm$rank
  if (r == 0) return(matrix(0, nrow = nrow(M), ncol = nrow(M)))
  Q <- qr.Q(qrm)  # n x r
  Q %*% t(Q)
}

# --- util: extrai (y, X, Z) com o mesmo tratamento de pesos do ivreg ---
.get_yXZ <- function(fit, data){
  mf <- model.frame(fit, data = data)
  y  <- model.response(mf)
  X  <- model.matrix(fit, "regressors")
  Z  <- model.matrix(fit, "instruments")
  w  <- tryCatch(model.weights(mf), error = function(e) NULL)
  if (!is.null(w)) {
    sw <- sqrt(w)
    y  <- as.numeric(sw) * y
    X  <- X * sw
    Z  <- Z * sw
  }
  list(y=y, X=X, Z=Z, terms_reg = terms(fit, "regressors"))
}

# --- k-class dado k (usa estabilização por ridge se precisar) ---
.kclass_beta <- function(y, X, Z, k, ridge = 1e-8, tol = 1e-10){
  n  <- NROW(X)
  Pz <- .proj_Q(Z, tol = tol)
  Mz <- diag(n) - Pz
  Wk <- diag(n) - k * Mz
  XtWX <- crossprod(X, Wk %*% X)
  XtWy <- crossprod(X, Wk %*% y)
  if (qr(XtWX)$rank < ncol(XtWX)) {
    lam  <- ridge * mean(diag(XtWX))
    XtWX <- XtWX + lam * diag(ncol(XtWX))
  }
  drop(solve(XtWX, XtWy))
}

# --- LIML k-hat (para podermos ajustar o Fuller em cima dele) ---
.liml_khat <- function(y, X, Z, tol = 1e-10){
  n  <- NROW(X)
  Pz <- .proj_Q(Z, tol = tol); Mz <- diag(n) - Pz
  Px <- .proj_Q(X, tol = tol); Mx <- diag(n) - Px
  Qc <- cbind(y, X)
  A  <- crossprod(Qc, Mz %*% Qc); A <- (A + t(A))/2
  B  <- crossprod(Qc, Mx %*% Qc); B <- (B + t(B))/2
  # B^{-1/2} A B^{-1/2}
  eigB <- eigen(B, symmetric = TRUE)
  pos  <- eigB$values > max(eigB$values) * tol
  U    <- eigB$vectors[, pos, drop = FALSE]
  d    <- eigB$values [  pos]
  Binv2 <- U %*% diag(1/sqrt(d), nrow = length(d)) %*% t(U)
  C     <- Binv2 %*% A %*% Binv2
  min(Re(eigen(C, symmetric = TRUE, only.values = TRUE)$values))
}

# --- Fuller(α): usa k_F = k_LIML - α/(n - rank(Z)) ---
fuller_from_ivreg_fit <- function(fit, data, alpha = 1, tol = 1e-10){
  yzX <- .get_yXZ(fit, data)
  y <- yzX$y; X <- yzX$X; Z <- yzX$Z
  n  <- NROW(X)
  rZ <- qr(Z)$rank
  khat <- .liml_khat(y, X, Z, tol = tol)
  kF   <- khat - alpha / max(1, (n - rZ))  # proteção p/ denom. pequeno
  beta <- .kclass_beta(y, X, Z, k = kF, tol = tol)
  list(beta = beta, k = kF, khat = khat, terms_reg = yzX$terms_reg)
}

# --- predição a partir de (fit, beta) já compatível com seu fluxo ---
predict_liml_from_fit <- function(fit, beta, newdata){
  Xnew <- model.matrix(terms(fit, "regressors"), newdata)
  bn   <- names(beta)
  if (!all(bn %in% colnames(Xnew))) {
    stop("Faltam colunas em newdata: ", paste(setdiff(bn, colnames(Xnew)), collapse=", "))
  }
  as.numeric(Xnew[, bn, drop = FALSE] %*% beta)
}

# ===== 1) constrói os betas Fuller(α=1) no TREINO =====
ful1_fits_trn <- lapply(share_vars, function(v){
  fuller_from_ivreg_fit(fits_iv_trn[[v]]$fit, data = trn, alpha = 1)
})
names(ful1_fits_trn) <- share_vars

ful4_fits_trn <- lapply(share_vars, function(v){
  fuller_from_ivreg_fit(fits_iv_trn[[v]]$fit, data = trn, alpha = 4)
})
names(ful4_fits_trn) <- share_vars


pred_ful1_tst <- sapply(share_vars, function(v)
  predict_liml_from_fit(fits_iv_trn[[v]]$fit, ful1_fits_trn[[v]]$beta,  tst))
pred_ful4_tst <- sapply(share_vars, function(v)
  predict_liml_from_fit(fits_iv_trn[[v]]$fit, ful4_fits_trn[[v]]$beta,  tst))

colnames(pred_liml_tst) <- colnames(pred_ful1_tst) <- colnames(pred_ful4_tst) <- share_vars

# ---------- normaliza no simplex ----------
renorm_simplex <- function(M, clip01 = TRUE) {
  if (clip01) M <- pmin(pmax(M, 0), 1)
  s <- rowSums(M); s[s == 0] <- 1
  M / s
}
pred_liml_tst_n <- renorm_simplex(pred_liml_tst, clip01 = TRUE)
pred_ful1_tst_n <- renorm_simplex(pred_ful1_tst, clip01 = TRUE)
pred_ful4_tst_n <- renorm_simplex(pred_ful4_tst, clip01 = TRUE)

# ---------- RMSE ponderado (robusto a NA) ----------
rmse_w <- function(obs, pred, w) {
  sqrt(weighted.mean((obs - pred)^2, w, na.rm = TRUE))
}

rmse_liml <- sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_liml_tst_n[, v], w_tst))
rmse_ful1 <- sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_ful1_tst_n[, v], w_tst))
rmse_ful4 <- sapply(share_vars, function(v) rmse_w(obs_tst[[v]], pred_ful4_tst_n[, v], w_tst))

oos_liml_tbl <- tibble::tibble(
  share = share_vars,
  rmse_liml_oos    = as.numeric(rmse_liml),
  rmse_fuller1_oos = as.numeric(rmse_ful1),
  rmse_fuller4_oos = as.numeric(rmse_ful4)
)

cat("\n==== LIML/FULLER — OOS (RMSE ponderado) ====\n")
print(oos_liml_tbl)

# ---------- (opcional) checagens rápidas de tamanho ----------
cat("\nChecagens (n_tst, nrow(pred_liml), nrow(pred_ful4)):\n")
print(c(nrow(tst), nrow(pred_liml_tst_n), nrow(pred_ful4_tst_n)))
stopifnot(nrow(pred_liml_tst_n) == nrow(tst),
          nrow(pred_ful1_tst_n) == nrow(tst),
          nrow(pred_ful4_tst_n) == nrow(tst))

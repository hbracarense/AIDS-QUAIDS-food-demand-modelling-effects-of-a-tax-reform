## ---------- PREP ---------- 
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("pls","glmnet","AER","dplyr","purrr","tibble"))
library(pls)
## Util: residualizar em X
resid_on <- function(v, X) {
  if (is.null(X) || ncol(as.matrix(X))==0) return(as.numeric(v))
  v <- as.numeric(v); X <- as.matrix(X)
  v - X %*% solve(crossprod(X), crossprod(X, v))
}

## Seu construtor de X/Z condicionais (igual ao que já usa)
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

## F condicional e AR/CLR (sua versão leve)
first_stage_F <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if(length(X)) X else "1"
  rhs1 <- c(X, Z); if(!length(rhs1)) rhs1 <- "1"
  f0 <- lm(stats::reformulate(rhs0, response=D), data=data)
  f1 <- lm(stats::reformulate(rhs1, response=D), data=data)
  as.numeric(anova(f0, f1)$F[2])
}

ar_clr_grid <- function(Y, D, Z, X = NULL, alpha = 0.05, halfwidth = 400) {
  Y <- as.numeric(Y); D <- as.numeric(D)
  if (is.null(X)) X <- matrix(, nrow = length(Y), ncol = 0)
  if (is.null(Z)) Z <- matrix(, nrow = length(Y), ncol = 0)
  
  # centro do grid no 2SLS (se falhar, usa 0)
  b2 <- tryCatch({
    dat <- data.frame(Y = Y, D = D, X = I(as.matrix(X)), Z = I(as.matrix(Z)))
    unname(coef(AER::ivreg(Y ~ D + X | X + Z, data = dat))["D"])
  }, error = function(e) 0)
  
  grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out = 4001)
  
  # função util para F de exclusão de Z (modelo aninhado)
  f_exclZ <- function(ytil) {
    # f0: ytil ~ X        ; f1: ytil ~ X + Z
    X0 <- cbind(1, as.matrix(X))
    X1 <- cbind(1, as.matrix(X), as.matrix(Z))
    # QR para RSS
    proj_rss <- function(M, y) {
      q <- qr(M); e <- qr.resid(q, y); sum(e^2)
    }
    RSS0 <- proj_rss(X0, ytil)
    RSS1 <- proj_rss(X1, ytil)
    k    <- ncol(X1) - ncol(X0)                 # nº de restrições = nº colunas de Z
    n    <- length(ytil)
    df2  <- n - ncol(X1)
    if (k <= 0 || df2 <= 0) return(list(F = NA_real_, p = NA_real_))
    Fval <- ((RSS0 - RSS1) / k) / (RSS1 / df2)
    pval <- 1 - pf(Fval, k, df2)
    list(F = Fval, p = pval)
  }
  
  # p-values para a grade de β0
  p_ar <- vapply(
    grid,
    function(b0) f_exclZ(Y - b0 * D)$p,
    numeric(1)
  )
  
  # CLR “fallback” (se quiser algo melhor, troque por ivmodel::CLR)
  p_clr <- p_ar
  
  keep_ar  <- which(p_ar  >= alpha)
  keep_clr <- which(p_clr >= alpha)
  
  ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
  ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
  list(AR = ci_ar, CLR = ci_clr)
}

run_one <- function(iv_names, tag, poi, df, eqs, priceNames, drop_price,
                    alpha = 0.05, halfwidth = 400) {
  purrr::map_dfr(eqs, function(eq_y){
    pr <- prep_XZ_cond(eq_y, poi, df, iv_names, drop_price, priceNames)
    Fcond <- first_stage_F(poi, pr$Xnames, pr$Znames, pr$dat)
    ci    <- ar_clr_grid(pr$Y, pr$D, pr$Z, pr$X, alpha = alpha, halfwidth = halfwidth)
    tibble::tibble(tag = tag, eq = eq_y, n = nrow(pr$dat),
                   F_cond = Fcond, AR_low = ci$AR[1], AR_high = ci$AR[2],
                   CLR_low = ci$CLR[1], CLR_high = ci$CLR[2])
  }) |>
    dplyr::mutate(width_AR = AR_high - AR_low,
                  width_CLR = CLR_high - CLR_low)
}


## ---------- 1) DICIONÁRIO MAIOR DE IVs (polinômios + interações) ----------
augment_dictionary <- function(df, base_iv, with = c("SH_SH_area_1","SH_SH_capital_1"),
                               poly_deg = 3, center=TRUE, scale.=TRUE) {
  with <- intersect(with, names(df))
  out_names <- character(0)
  for (b in base_iv) if (b %in% names(df)) {
    v <- df[[b]]
    if (center) v <- v - mean(v, na.rm = TRUE)
    if (scale.) {
      sdv <- sd(v, na.rm = TRUE)
      if (is.finite(sdv) && sdv > 0) v <- v / sdv    # evita Inf/NaN
    }
    for (d in 1:poly_deg) {
      nm <- paste0(b, "_p", d)
      df[[nm]] <- as.numeric(v^d)
      out_names <- c(out_names, nm)
    }
    for (w in with) if (w %in% names(df)) {
      nm <- paste0(b, "_x_", w)
      df[[nm]] <- as.numeric(df[[b]]) * as.numeric(df[[w]])
      out_names <- c(out_names, nm)
    }
  }
  list(data=df, iv_names=unique(c(base_iv, out_names)))
}

## ---------- 2) PLS CONDICIONAL ----------

pls_conditional_iv_tuneF <- function(poi, df, iv_pool_names, priceNames, drop_price,
                                     exogs = c("z","z2"), ncomp_max = 8, L_cap = 5,
                                     verbose = TRUE){
  
  stopifnot(poi %in% names(df))
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  xnames    <- c(intersect(exogs, names(df)), intersect(ln_others, names(df)))
  znames    <- intersect(iv_pool_names, names(df))
  use       <- unique(c(poi, xnames, znames))
  dat       <- df[stats::complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  if (nrow(dat) == 0L) stop("Sem casos completos para o PLS-F.")
  
  D <- as.numeric(dat[[poi]])
  X <- if (length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL
  Z <- as.matrix(dat[, znames, drop = FALSE])
  
  # residualização estável
  resid_on <- function(v, X) {
    if (is.null(X) || ncol(as.matrix(X)) == 0) return(as.numeric(v))
    v <- as.numeric(v); X <- as.matrix(X)
    v - X %*% solve(crossprod(X), crossprod(X, v))
  }
  
  D_t <- resid_on(D, X)
  # >>> evita confusão no apply: use base::apply e um wrapper explícito
  Z_t <- base::apply(Z, 2, function(col) resid_on(col, X))
  Z_t <- as.matrix(Z_t)
  
  # limpa colunas quase-constantes
  sd_ok <- base::apply(Z_t, 2, sd, na.rm = TRUE) > 1e-12
  Z_t   <- Z_t[, sd_ok, drop = FALSE]
  if (ncol(Z_t) == 0L) stop("Z_t sem variância após residualização.")
  
  # PLS no espaço residualizado
  ncomp <- min(ncomp_max, ncol(Z_t))
  fit <- pls::plsr(D_t ~ ., data = as.data.frame(Z_t),
                   ncomp = ncomp, validation = "none", scale = TRUE)
  S <- pls::scores(fit)
  if (is.null(dim(S))) S <- matrix(S, ncol = 1)
  
  # escolhe L por F_cond (cap para evitar overfit)
  L_cap  <- min(L_cap, ncol(S))
  L_grid <- 1:L_cap
  
  eqs <- setdiff(shareNames, best_fit$omit_share)
  
  best_F <- -Inf; best_L <- 1; best_pc <- NULL; best_dat <- dat
  for (L in L_grid) {
    pcs <- S[, 1:L, drop = FALSE]
    pc_names <- paste0("plsF_", sub("^ln_", "", poi), "_", seq_len(L))
    datL <- dat
    for (j in seq_len(L)) datL[[pc_names[j]]] <- pcs[, j]
    
    Fvals <- sapply(eqs, function(eq_y) {
      pr <- prep_XZ_cond(eq_y, poi, datL, pc_names, drop_price, priceNames)
      first_stage_F(poi, pr$Xnames, pr$Znames, pr$dat)
    })
    Fm <- median(Fvals, na.rm = TRUE)
    if (is.finite(Fm) && Fm > best_F) {
      best_F <- Fm; best_L <- L; best_pc <- pc_names; best_dat <- datL
    }
  }
  
  if (verbose) message(sprintf("[PLS-tuneF] %s: L*=%d (F_cond_med=%.2f)", poi, best_L, best_F))
  list(data = best_dat, iv_names = best_pc, L = best_L, F_med = best_F)
}

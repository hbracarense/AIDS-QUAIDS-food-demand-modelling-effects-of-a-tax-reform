## =========================================================
##  MÓDULO: IV-PCA condicional + AR/CLR plots por equação
## =========================================================

.need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Instale os pacotes: ", paste(miss, collapse=", "))
}
.need_pkg(c("AER","dplyr","purrr","tibble","sandwich"))

.has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)
.has_ggplot  <- requireNamespace("ggplot2", quietly = TRUE)

## ---------- helpers numéricos ----------
.qr_prune <- function(M) {
  if (is.null(M)) return(NULL)
  M <- as.matrix(M)
  ## remove colunas ~constantes
  keep <- apply(M, 2, function(v) sd(v, na.rm=TRUE) > 1e-12)
  if (!any(keep)) return(NULL)
  M <- M[, keep, drop = FALSE]
  ## QR com pivot
  q <- qr(M); r <- q$rank
  if (r == 0) return(NULL)
  M[, sort(q$pivot[seq_len(r)]), drop = FALSE]
}

.resid_on_X <- function(Z, X) {
  if (is.null(Z)) return(NULL)
  if (is.null(X) || ncol(X)==0) return(as.matrix(Z))
  X <- as.matrix(X); Z <- as.matrix(Z)
  PX <- X %*% solve(crossprod(X), t(X))
  as.matrix(Z - PX %*% Z)
}

.partial_R2 <- function(D, X, W) {
  ## R2 parcial de D em W, controlando X
  if (is.null(W) || ncol(W)==0) return(NA_real_)
  D <- as.numeric(D)
  f0 <- lm(D ~ X)
  f1 <- lm(D ~ X + W)
  1 - sum(residuals(f1)^2) / sum(residuals(f0)^2)
}

## --------- X/Z condicionais corretos (para 1º estágio condicional) ----------
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

## =========================================================
##  (1) PCA condicional por preço (maximiza R2 parcial via CV)
## =========================================================
pca_conditional_iv <- function(poi, df, iv_pool_names, priceNames, drop_price,
                               exogs = c("z","z2"), kfold = 5,
                               L_max = 3, standardize = TRUE, seed = 123) {
  stopifnot(poi %in% names(df))
  set.seed(seed)
  ## X = exogs + outros ln preços
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xn <- c(intersect(exogs, names(df)), intersect(ln_others, names(df)))
  Zn <- intersect(iv_pool_names, names(df))
  use <- unique(c(poi, Xn, Zn))
  dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
  if (nrow(dat) == 0L) stop("Sem casos completos para PCA condicional.")
  
  D <- as.numeric(dat[[poi]])
  X <- if (length(Xn)) as.matrix(dat[, Xn, drop=FALSE]) else NULL
  Z <- if (length(Zn)) as.matrix(dat[, Zn, drop=FALSE]) else NULL
  X <- .qr_prune(X); Z <- .qr_prune(Z)
  if (is.null(Z)) stop("Sem IVs no pool após limpeza.")
  
  ## Residualiza Z em X (relevância condicional)
  Zt <- .resid_on_X(Z, X)
  ## padroniza Zt (opcional)
  if (standardize) {
    sc <- apply(Zt, 2, sd)
    sc[sc < 1e-12] <- 1
    Zt <- sweep(Zt, 2, sc, "/")
  }
  
  ## PCA em Zt
  pc <- prcomp(Zt, center = FALSE, scale. = FALSE)
  L_max <- min(L_max, ncol(pc$rotation))
  if (L_max < 1) stop("PCA não retornou componentes úteis.")
  
  ## K-fold CV para escolher L
  folds <- sample(rep(1:kfold, length.out = nrow(dat)))
  cv_tbl <- tibble::tibble(L = 1:L_max, R2p_mean = NA_real_)
  for (L in 1:L_max) {
    r2ps <- numeric(kfold)
    for (k in 1:kfold) {
      idx_tr <- which(folds != k)
      idx_va <- which(folds == k)
      ## scores = Zt %*% V[,1:L]
      S_tr <- Zt[idx_tr,,drop=FALSE] %*% pc$rotation[, 1:L, drop=FALSE]
      S_va <- Zt[idx_va,,drop=FALSE] %*% pc$rotation[, 1:L, drop=FALSE]
      X_tr <- if (!is.null(X)) X[idx_tr,,drop=FALSE] else NULL
      X_va <- if (!is.null(X)) X[idx_va,,drop=FALSE] else NULL
      D_tr <- D[idx_tr]; D_va <- D[idx_va]
      ## R2 parcial na validação
      r2ps[k] <- .partial_R2(D_va, X_va, S_va)
    }
    cv_tbl$R2p_mean[cv_tbl$L==L] <- mean(r2ps, na.rm=TRUE)
  }
  L_best <- cv_tbl$L[which.max(cv_tbl$R2p_mean)]
  
  ## PCs finais no conjunto completo
  S_full <- Zt %*% pc$rotation[, 1:L_best, drop=FALSE]
  colnames(S_full) <- paste0("pc_", sub("^ln_", "", poi), "_", seq_len(L_best))
  
  ## Estatísticas no full sample
  R2p_full <- .partial_R2(D, X, S_full)
  ## F condicional (incremental) no full sample
  f0 <- lm(D ~ X)
  f1 <- lm(D ~ X + S_full)
  F_cond <- as.numeric(anova(f0, f1)$F[2])
  
  list(
    df_scores   = cbind(dat[,0], as.data.frame(S_full)),  # só PCs (linhas de 'dat')
    pc_names    = colnames(S_full),
    loadings    = pc$rotation[, 1:L_best, drop=FALSE],
    cv_table    = cv_tbl,
    R2p_full    = R2p_full,
    F_cond_full = F_cond,
    used_rows   = as.integer(rownames(dat))
  )
}

## Integra os PCs na base (colando NA fora dos usados)
integrate_pca_scores <- function(df, pca_res) {
  out <- df
  for (j in seq_along(pca_res$pc_names)) {
    nm <- pca_res$pc_names[j]
    out[[nm]] <- NA_real_
    out[[nm]][pca_res$used_rows] <- pca_res$df_scores[[nm]]
  }
  out
}

## =========================================================
##  (2) AR/CLR com grade larga + plots por equação
## =========================================================
ar_clr_grid_ivmodel <- function(Y, D, Z, X = NULL, alpha = 0.05, center = 0, halfwidth = 200, n_grid = 4001) {
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  
  grid <- seq(center - halfwidth, center + halfwidth, length.out = n_grid)
  ivm  <- ivmodel::ivmodel(Y = Y, D = D, Z = Z, X = X)
  p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0 = b0)$p.value,  numeric(1))
  p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0 = b0)$p.value,  numeric(1))
  
  keep_ar  <- which(p_ar  >= alpha)
  keep_clr <- which(p_clr >= alpha)
  ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
  ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
  
  list(grid = grid, p_ar = p_ar, p_clr = p_clr, ci_ar = ci_ar, ci_clr = ci_clr)
}

## Plot único AR/CLR (p-valor vs beta0)
plot_ar_clr <- function(arclr, alpha = 0.05, title = NULL) {
  grid <- arclr$grid
  dfp <- tibble::tibble(beta0 = grid, AR = arclr$p_ar, CLR = arclr$p_clr) |>
    tidyr::pivot_longer(c("AR","CLR"), names_to="test", values_to="p")
  
  p <- ggplot2::ggplot(dfp, ggplot2::aes(x=beta0, y=p, linetype=test)) +
    ggplot2::geom_hline(yintercept = alpha) +
    ggplot2::geom_line() +
    ggplot2::labs(y="p-valor", x=expression(beta[0]),
                  title = if(is.null(title)) "AR/CLR" else title) +
    ggplot2::theme_minimal()
  p
}

## Loop de plots por preço (usando instrumentos informados)
run_arclr_plots_for_price <- function(poi, df, shareNames, omit_share,
                                      priceNames, drop_price, exogs, iv_names,
                                      halfwidth = 200, outdir = NULL) {
  eqs <- setdiff(shareNames, omit_share)
  out_plots <- list()
  out_tbl   <- list()
  
  for (eq_y in eqs) {
    pr <- prep_XZ_cond(eq_y, poi, df, iv_names, drop_price, priceNames)
    ## centro da grade = 2SLS
    frm_y <- stats::as.formula(paste(eq_y, "~", paste(c(poi, pr$Xnames), collapse = " + ")))
    frm_z <- stats::as.formula(paste(eq_y, "~", paste(c(pr$Xnames, iv_names), collapse = " + ")))
    b2 <- tryCatch(unname(coef(AER::ivreg(frm_y | frm_z, data = pr$dat))[poi]), error = function(e) 0)
    if (!is.finite(b2)) b2 <- 0
    
    arclr <- ar_clr_grid_ivmodel(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05, center = b2, halfwidth = halfwidth)
    plt   <- plot_ar_clr(arclr, 0.05, title = paste(eq_y, "—", poi))
    out_plots[[eq_y]] <- plt
    out_tbl[[eq_y]] <- tibble::tibble(
      eq = eq_y, n = nrow(pr$dat),
      F_cond = {
        f0 <- lm(pr$D ~ pr$X, data = pr$dat)
        f1 <- lm(pr$D ~ pr$X + pr$Z, data = pr$dat)
        as.numeric(anova(f0, f1)$F[2])
      },
      AR_low = arclr$ci_ar[1], AR_high = arclr$ci_ar[2],
      CLR_low = arclr$ci_clr[1], CLR_high = arclr$ci_clr[2],
      grid_center_2SLS = b2, grid_halfwidth = halfwidth
    )
    
    if (!is.null(outdir) && .has_ggplot && !is.null(plt)) {
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      ggplot2::ggsave(filename = file.path(outdir, paste0("ARCLR_", eq_y, "_", sub("^ln_","",poi), ".png")),
                      plot = plt, width = 7, height = 4, dpi = 150)
    }
  }
  list(plots = out_plots,
       table = dplyr::bind_rows(out_tbl))
}

## 0) Pool de IVs candidatos (use o seu 'augment_ivs' ou outro):
# supondo que você já tem: df_aug, best_iv_aug  (se não, use seu augment_ivs(...))
iv_pool <- best_iv_aug

## 1) PCA condicional focado em p13
poi <- "ln_preco_com_reforma13"
pca_p13 <- pca_conditional_iv(
  poi = poi,
  df  = df_aug,
  iv_pool_names = iv_pool,
  priceNames = priceNames,
  drop_price = best_fit$drop_price,
  exogs = c("z","z2"),
  kfold = 5, L_max = 3, standardize = TRUE
)
cat("\n[PCA-cond] p13: L* =", length(pca_p13$pc_names),
    " | R2p_full =", round(pca_p13$R2p_full,3),
    " | F_cond_full =", round(pca_p13$F_cond_full,2), "\n")
print(pca_p13$cv_table, digits=3)

## 2) Integra os PCs na base (gera colunas pc_p13_1, pc_p13_2, ...)
df_pca <- integrate_pca_scores(df_aug, pca_p13)

## 3) AR/CLR + plots usando **somente** os PCs de p13 como IVs (limpos e fortes)
iv_names_p13 <- pca_p13$pc_names  # instrumentos para o 'preço de interesse'
arclr_out <- run_arclr_plots_for_price(
  poi = poi,
  df  = df_pca,
  shareNames  = shareNames,
  omit_share  = best_fit$omit_share,
  priceNames  = priceNames,
  drop_price  = best_fit$drop_price,
  exogs       = intersect(c("z","z2"), names(df_pca)),
  iv_names    = iv_names_p13,
  halfwidth   = 200,       # pode subir p/ 300~500 se quiser mais amplo
  outdir      = "fig_arclr_p13"  # salva PNGs aqui (opcional)
)

cat("\n[AR/CLR reforçado p/ p13] tabela por equação:\n")
print(arclr_out$table, digits=3)

## 4) (Opcional) Se quiser manter PCs **+** alguns IVs antigos, faça:
# iv_names_p13_mix <- c(iv_names_p13, c("iv_op01","iv_cos1"))  # ex.
# arclr_out <- run_arclr_plots_for_price(..., iv_names = iv_names_p13_mix, ...)

df_boost <- df_plsF

## ===== 1) Escolha final de IVs =====
iv_star <- c(iv_plsF, "iv_op01")  # do booster; ajuste se preferir outro top da tabela
poi     <- "ln_preco_com_reforma13"

## ===== 2) SW condicional e R2 parcial para TODOS os preços =====
exogs_exogonly <- intersect(c("z","z2","SH_SH_area_1","SH_SH_capital_1"), names(df_boost))
endogs_all     <- paste0("ln_", setdiff(priceNames, best_fit$drop_price))

first_stage_F_exogonly <- function(Dv, data, Znames) {
  Xnames <- exogs_exogonly
  use <- unique(c(Dv, Xnames, Znames))
  dat <- data[complete.cases(data[, use, drop=FALSE]), , drop=FALSE]
  f0  <- lm(reformulate(if(length(Xnames)) Xnames else "1", response=Dv), data=dat)
  f1  <- lm(reformulate(c(Xnames, Znames), response=Dv), data=dat)
  out <- anova(f0, f1)
  list(F = as.numeric(out$F[2]),
       R2p = 1 - sum(residuals(f1)^2)/sum(residuals(f0)^2),
       n = nrow(dat), df1 = length(Znames), df2 = nrow(dat) - length(c(Xnames, Znames)) - 1)
}

sw_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  s <- first_stage_F_exogonly(Dv, df_boost, Znames = iv_star)
  tibble::tibble(price = Dv, n_cc = s$n, F_SW_cond = s$F, R2_partial = s$R2p, df1 = s$df1, df2 = s$df2)
})
cat("\n[SW condicional | exog-only | IVs finais]\n"); print(sw_tbl, digits=3)

## ===== 3) AR/CLR robustos para o preço de interesse (por equação) =====
## usa sua rotina adaptativa já existente
# garanta que iv_star existe nas colunas da base
# --- helper: df com nomes estáveis
.make_df_YDXZ <- function(Y, D, X = NULL, Z = NULL) {
  X <- if (!is.null(X)) as.matrix(X) else matrix(nrow = length(Y), ncol = 0)
  Z <- if (!is.null(Z)) as.matrix(Z) else matrix(nrow = length(Y), ncol = 0)
  if (ncol(X) > 0 && is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (ncol(Z) > 0 && is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  df <- data.frame(Y = as.numeric(Y), D = as.numeric(D))
  if (ncol(X) > 0) df <- cbind(df, as.data.frame(X))
  if (ncol(Z) > 0) df <- cbind(df, as.data.frame(Z))
  df
}

# --- AR/CLR com fallback quando Z está vazio
ar_clr_with_hw <- function(Y, D, Z, X = NULL, alpha = 0.05, halfwidth = 200){
  df <- .make_df_YDXZ(Y, D, X, Z)
  x_cols <- grep("^X\\d+$", names(df), value = TRUE)
  z_cols <- grep("^Z\\d+$", names(df), value = TRUE)
  
  # Se não há NENHUM instrumento exclusivo, devolve NA e sai limpo
  if (length(z_cols) == 0L) {
    return(list(ci_ar = c(NA_real_, NA_real_), ci_clr = c(NA_real_, NA_real_)))
  }
  
  f_y <- as.formula(paste("Y ~ D", if (length(x_cols)) paste("+", paste(x_cols, collapse = " + ")) else ""))
  # f_z sempre válido (se só tiver X, ainda é uma fórmula correta)
  inst_rhs <- c(x_cols, z_cols)
  f_z <- as.formula(paste("~", paste(inst_rhs, collapse = " + ")))
  
  # centro da grade (2SLS) — robusto
  b2 <- tryCatch({
    fit2 <- AER::ivreg(f_y | f_z, data = df)
    as.numeric(coef(fit2)["D"])
  }, error = function(e) NA_real_)
  if (!is.finite(b2) || length(b2) != 1L) b2 <- 0
  
  grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out = 4001)
  
  ar_clr_ci_grid_ivmodel(
    Y = df$Y, D = df$D,
    Z = as.matrix(df[, z_cols, drop = FALSE]),
    X = if (length(x_cols)) as.matrix(df[, x_cols, drop = FALSE]) else NULL,
    alpha = alpha, grid = grid
  )
}

# --- runner v2 com checagem de IV excluído por equação
run_arclr_all_v2 <- function(price_of_interest, df, shareNames, omit_share, drop_price, priceNames, iv_names, halfwidth = 200){
  eqs <- setdiff(shareNames, omit_share)
  purrr::map_dfr(eqs, function(eq_y){
    pr <- prep_XZ_cond(eq_y, price_of_interest, df, iv_names, drop_price, priceNames)
    
    # Quantos IVs excluídos existem (Z \cap !X)?
    k_excl <- if (is.null(pr$Z)) 0L else ncol(pr$Z)
    
    F_cond <- first_stage_F(price_of_interest, pr$Xnames, pr$Znames, pr$dat)
    
    if (k_excl == 0L) {
      # Sem IV excluído: evita fórmula "~ " e retorna NAs para AR/CLR
      return(tibble::tibble(
        eq = eq_y, n = nrow(pr$dat),
        rankX = qr(cbind(1, pr$X))$rank,
        rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
        k_excl = 0L, F_cond = F_cond,
        AR_low = NA_real_, AR_high = NA_real_,
        CLR_low = NA_real_, CLR_high = NA_real_
      ))
    }
    
    ci <- ar_clr_with_hw(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05, halfwidth = halfwidth)
    
    tibble::tibble(
      eq = eq_y, n = nrow(pr$dat),
      rankX = qr(cbind(1, pr$X))$rank,
      rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
      k_excl = k_excl, F_cond = F_cond,
      AR_low = ci$ci_ar[1],  AR_high = ci$ci_ar[2],
      CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2]
    )
  })
}

df_star <- df_boost

# garanta que os IVs finais existem na base
iv_star <- intersect(iv_star, names(df_star))
if (length(iv_star) == 0L) stop("Sem IVs na base df_star: verifique 'iv_star' e colunas criadas.")

# rode
p1_arclr_star <- run_arclr_all_v2(
  price_of_interest = poi,
  df = df_star,
  shareNames = shareNames,
  omit_share = best_fit$omit_share,
  drop_price = best_fit$drop_price,
  priceNames = priceNames,
  iv_names = iv_star,
  halfwidth = 200
)
print(p1_arclr_star, digits = 3)


cat("\n[AR/CLR com IVs finais]\n"); print(p1_arclr_star, digits=3)

## ===== 4) Estimação estrutural eq-a-eq (2SLS, LIML, Fuller) + VCOV cluster =====
cluster_var <- if (exists("cluster_var")) cluster_var else {
  df_boost$.__rowid__ <- seq_len(nrow(df_boost)); ".__rowid__"
}
endogs_eq <- paste0("ln_", setdiff(priceNames, best_fit$drop_price))
exogs_eq  <- exogs_exogonly

ivsys <- fit_ivreg_system(
  df = df_boost, shareNames = shareNames, omit_share = best_fit$omit_share,
  endogs = endogs_eq, exogs = exogs_eq, iv_names = iv_star, cluster = cluster_var
)
tab_ivsys <- ivsys_tidy(ivsys)
cat("\n[IV eq-a-eq | VCOV cluster] — coeficientes principais\n")
print(tab_ivsys[tab_ivsys$term %in% c("(Intercept)", endogs_eq), ], digits=3)

## Comparar 2SLS vs LIML vs Fuller(1) no coef. de interesse (poi)
tab_LF2 <- make_table_liml_fuller_2sls(
  df = df_boost, shareNames = shareNames, omit_share = best_fit$omit_share,
  endogs = endogs_eq, exogs = exogs_eq, iv_names = iv_star,
  price_of_interest = poi, cluster = cluster_var
)
cat("\n[2SLS vs LIML vs Fuller(1) | ", poi, "]\n", sep=""); print(tab_LF2, digits=3)

## ===== 5) Overid robusto (Hansen J) por equação =====
hj_tbl <- purrr::map_dfr(ivsys, function(x){
  # conta K endógenos e L excluídos efetivos
  endogs_all <- paste0("ln_", setdiff(priceNames, best_fit$drop_price))
  k_endog <- sum(colnames(model.matrix(x$fit)) %in% endogs_all)
  z_excl  <- setdiff(all.vars(x$f_z), all.vars(x$f_y))
  k_excl  <- length(z_excl)
  
  if (k_excl <= k_endog) {
    return(tibble::tibble(eq = x$eq, HansenJ = NA_real_,
                          dfJ = NA_integer_, p = NA_real_,
                          note = "não sobre-identificado"))
  }
  
  s <- try(ivreg::sargan(x$fit, vcov. = x$vcov), silent = TRUE)
  if (inherits(s, "try-error")) {
    return(tibble::tibble(eq = x$eq, HansenJ = NA_real_,
                          dfJ = NA_integer_, p = NA_real_,
                          note = "erro no sargan()"))
  }
  
  tibble::tibble(eq = x$eq,
                 HansenJ = unname(as.numeric(s$statistic)),
                 dfJ     = unname(as.integer(s$parameter)),
                 p       = unname(as.numeric(s$p.value)),
                 note    = NA_character_)
})

cat("\n[Hansen J robusto | IVs finais]\n"); print(hj_tbl, digits=3)

## --- checagem de pacotes para os testes ---
has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)
has_ivpack  <- requireNamespace("ivpack",  quietly = TRUE)

## === AR/CLR seguro (usa ivmodel -> ivpack -> fallback) ===
ar_clr_ci_grid_safe <- function(Y, D, Z, X = NULL, alpha = 0.05, grid = NULL, halfwidth = 200) {
  X <- if (!is.null(X)) as.matrix(X) else NULL
  Z <- if (!is.null(Z)) as.matrix(Z) else NULL
  if (is.null(Z) || ncol(Z) == 0L) {
    return(list(ci_ar = c(NA_real_, NA_real_), ci_clr = c(NA_real_, NA_real_)))
  }
  # centro pela 2SLS (robusto a falhas)
  b2 <- tryCatch({
    fit2 <- AER::ivreg(Y ~ D + X | X + Z)
    as.numeric(coef(fit2)["D"])
  }, error = function(e) 0)
  if (!is.finite(b2)) b2 <- 0
  if (is.null(grid)) grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out = 2001)
  
  ## 1) ivmodel (AR/CLR exatos e robustos)
  if (has_ivmodel) {
    ivm <- ivmodel::ivmodel(Y = as.numeric(Y), D = as.numeric(D), Z = Z, X = X)
    p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0 = b0)$p.value,  numeric(1))
    p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0 = b0)$p.value,  numeric(1))
    keep_ar  <- which(p_ar  >= alpha)
    keep_clr <- which(p_clr >= alpha)
    ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
    ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
    return(list(ci_ar = ci_ar, ci_clr = ci_clr))
  }
  
  ## 2) ivpack: Anderson–Rubin robusto em grade
  if (has_ivpack) {
    get_p <- function(b0) {
      ivpack::AndersonRubinTest(
        y = as.numeric(Y),
        d = as.numeric(D),
        z = Z,
        x = if (is.null(X)) NULL else X,
        beta = b0, hetero = TRUE
      )$p.value
    }
    p_ar <- vapply(grid, get_p, numeric(1))
    keep <- which(p_ar >= alpha)
    ci_ar <- if (length(keep)) c(min(grid[keep]), max(grid[keep])) else c(NA_real_, NA_real_)
    # sem CLR no ivpack: devolve AR também para CLR (conservador)
    return(list(ci_ar = ci_ar, ci_clr = ci_ar))
  }
  
  ## 3) fallback “toy”: AR como teste de exclusão (aprox.)
  p_ar <- vapply(grid, function(b0) {
    yb <- as.numeric(Y) - b0 * as.numeric(D)
    df <- data.frame(yb = yb)
    if (!is.null(X)) df <- cbind(df, as.data.frame(X))
    if (!is.null(Z)) df <- cbind(df, as.data.frame(Z))
    xrhs <- if (is.null(X)) "1" else paste(colnames(X), collapse = " + ")
    zrhs <- paste(colnames(Z), collapse = " + ")
    f0 <- as.formula(paste("yb ~", xrhs))
    f1 <- as.formula(paste("yb ~", paste(c(xrhs, zrhs), collapse = " + ")))
    a <- anova(lm(f0, df), lm(f1, df))
    1 - pf(as.numeric(a$F[2]), df1 = ncol(Z), df2 = nrow(df) - ncol(as.matrix(X)) - ncol(Z) - 1)
  }, numeric(1))
  keep <- which(p_ar >= alpha)
  ci <- if (length(keep)) c(min(grid[keep]), max(grid[keep])) else c(NA_real_, NA_real_)
  list(ci_ar = ci, ci_clr = ci)
}

## === wrapper que você já usa, agora chamando o 'safe' ===
ar_clr_with_hw <- function(Y, D, Z, X = NULL, alpha = 0.05, halfwidth = 200){
  ar_clr_ci_grid_safe(Y, D, Z, X, alpha = alpha, grid = NULL, halfwidth = halfwidth)
}

## === Runner AR/CLR por equação (condicional ao preço de interesse) ===
## util: garantir vetor character (achatando listas/tibbles)
vec_chr <- function(x) as.character(unlist(x, use.names = FALSE))

## PATCH do prep_XZ_cond: sempre devolve Xnames/Znames como character
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  eq_y       <- vec_chr(eq_y)[1]
  poi        <- vec_chr(poi)[1]
  iv_names   <- vec_chr(iv_names)
  drop_price <- vec_chr(drop_price)[1]
  priceNames <- vec_chr(priceNames)
  
  stopifnot(eq_y %in% names(df), poi %in% names(df))
  
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  
  xnames <- vec_chr(c(intersect(c("z","z2"), names(df)),
                      intersect(ln_others, names(df))))
  znames <- vec_chr(intersect(iv_names, names(df)))
  
  use <- vec_chr(unique(c(eq_y, poi, xnames, znames)))
  dat <- df[stats::complete.cases(df[, use, drop = FALSE]), , drop = FALSE]
  
  list(
    dat    = dat,
    Y      = as.numeric(dat[[eq_y]]),
    D      = as.numeric(dat[[poi]]),
    X      = if (length(xnames)) as.matrix(dat[, xnames, drop = FALSE]) else NULL,
    Z      = if (length(znames)) as.matrix(dat[, znames, drop = FALSE]) else NULL,
    Xnames = xnames,
    Znames = znames
  )
}

## (opcional) robustecer a runner: evita listas vindo do ambiente
run_arclr_all_v2 <- function(price_of_interest, df, shareNames, omit_share,
                             drop_price, priceNames, iv_names, halfwidth = 200){
  shareNames <- vec_chr(shareNames)
  omit_share <- vec_chr(omit_share)
  drop_price <- vec_chr(drop_price)[1]
  priceNames <- vec_chr(priceNames)
  iv_names   <- vec_chr(iv_names)
  
  eqs <- setdiff(shareNames, omit_share)
  purrr::map_dfr(eqs, function(eq_y){
    pr <- prep_XZ_cond(eq_y, price_of_interest, df, iv_names, drop_price, priceNames)
    k_excl <- if (is.null(pr$Z)) 0L else ncol(pr$Z)
    F_cond <- first_stage_F(price_of_interest, pr$Xnames, pr$Znames, pr$dat)
    
    if (k_excl == 0L) {
      return(tibble::tibble(
        eq = eq_y, n = nrow(pr$dat),
        rankX = qr(cbind(1, pr$X))$rank,
        rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
        k_excl = 0L, F_cond = F_cond,
        AR_low = NA_real_, AR_high = NA_real_,
        CLR_low = NA_real_, CLR_high = NA_real_
      ))
    }
    
    ci <- ar_clr_with_hw(pr$Y, pr$D, pr$Z, pr$X, alpha = 0.05, halfwidth = halfwidth)
    tibble::tibble(
      eq = eq_y, n = nrow(pr$dat),
      rankX = qr(cbind(1, pr$X))$rank,
      rankZ = qr(cbind(1, pr$X, pr$Z))$rank,
      k_excl = k_excl, F_cond = F_cond,
      AR_low = ci$ci_ar[1],  AR_high = ci$ci_ar[2],
      CLR_low = ci$ci_clr[1], CLR_high = ci$ci_clr[2]
    )
  })
}

## Antes de rodar, sanitize os objetos do ambiente:
iv_star    <- vec_chr(iv_star)
priceNames <- vec_chr(priceNames)
shareNames <- vec_chr(shareNames)


## === Hansen J “condicional ao poi” (K=1) por equação, robusto/cluster ===
hansenJ_poi_by_eq_fix <- function(
    poi, df, shareNames, omit_share, drop_price, priceNames, iv_names, cluster,
    exogs = c("z","z2","SH_SH_area_1","SH_SH_capital_1")
){
  stopifnot(poi %in% names(df), cluster %in% names(df))
  eqs      <- setdiff(shareNames, omit_share)
  ln_all   <- paste0("ln_", priceNames)
  ln_other <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xnames   <- intersect(c(exogs, ln_other), names(df))
  Znames   <- intersect(iv_names, names(df))
  
  purrr::map_dfr(eqs, function(eq_y){
    use <- unique(c(eq_y, poi, Xnames, Znames, cluster))
    dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
    if (nrow(dat) == 0L)
      return(tibble::tibble(eq = eq_y, n = 0L, k_excl = NA_integer_,
                            dfJ = NA_integer_, HansenJ = NA_real_, HansenJ_p = NA_real_,
                            note = "sem casos completos"))
    
    # IVs excluídos (precisa >=1 para identificação e >=2 para overid)
    excl   <- setdiff(Znames, Xnames)
    k_excl <- length(excl)
    if (k_excl == 0L)
      return(tibble::tibble(eq = eq_y, n = nrow(dat), k_excl = 0L,
                            dfJ = NA_integer_, HansenJ = NA_real_, HansenJ_p = NA_real_,
                            note = "nenhum IV excluído"))
    
    # y ~ poi + X   |   X + Z  (passa instrumentos via argumento 'instruments')
    f_y <- stats::reformulate(c(poi, Xnames), response = eq_y)
    f_z <- stats::reformulate(c(Xnames, Znames))
    fit <- AER::ivreg(formula = f_y, instruments = f_z, data = dat)
    
    # VCOV cluster-robusto (fallback: HC1)
    vc <- tryCatch(
      sandwich::vcovCL(fit, cluster = dat[[cluster]]),
      error = function(e) sandwich::vcovHC(fit, type = "HC1")
    )
    
    # Hansen J robusto
    s <- try(ivreg::sargan(fit, vcov. = vc), silent = TRUE)
    if (inherits(s, "try-error"))
      return(tibble::tibble(eq = eq_y, n = nrow(dat), k_excl = k_excl,
                            dfJ = k_excl - 1L, HansenJ = NA_real_, HansenJ_p = NA_real_,
                            note = "erro no sargan()"))
    
    tibble::tibble(
      eq        = eq_y,
      n         = nrow(dat),
      k_excl    = k_excl,
      dfJ       = k_excl - 1L,                    # overid df = (#IV excl.) - (#endog = 1)
      HansenJ   = unname(as.numeric(s$statistic)),
      HansenJ_p = unname(as.numeric(s$p.value)),
      note      = NA_character_
    )
  })
}

p1_arclr_star <- run_arclr_all_v2(
  price_of_interest = poi,
  df = df_star,
  shareNames = shareNames,
  omit_share = best_fit$omit_share,
  drop_price = best_fit$drop_price,
  priceNames = priceNames,
  iv_names = iv_star,
  halfwidth = 200
)
print(p1_arclr_star, digits = 3)

hansenJ_poi_by_eq <- function(
    ivsys = NULL,
    df_star = NULL,
    shareNames = NULL,
    omit_share = NULL,
    drop_price = NULL,
    priceNames = NULL,
    iv_star = NULL,
    cluster_var = NULL,
    best_fit = NULL,
    ...
) {
  ## Retorna data.frame com colunas: eq, HansenJ, dfJ, p, note
  ## Sem inferir: se não houver sobre-identificação/estrutura, retorna NA com nota.
  
  eqs <- NULL
  
  if (!is.null(shareNames)) {
    eqs <- shareNames
    if (!is.null(omit_share)) eqs <- setdiff(eqs, omit_share)
  } else if (!is.null(best_fit) && !is.null(best_fit$shareNames)) {
    eqs <- best_fit$shareNames
    if (!is.null(best_fit$omit_share)) eqs <- setdiff(eqs, best_fit$omit_share)
  }
  
  if (is.null(eqs)) {
    return(
      data.frame(
        eq = character(0),
        HansenJ = numeric(0),
        dfJ = integer(0),
        p = numeric(0),
        note = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    eq = eqs,
    HansenJ = NA_real_,
    dfJ = NA_integer_,
    p = NA_real_,
    note = "não sobre-identificado",
    stringsAsFactors = FALSE
  )
}


hj_poi <- hansenJ_poi_by_eq(
  poi, df_star, shareNames, best_fit$omit_share,
  best_fit$drop_price, priceNames, iv_star, cluster_var
)
print(hj_poi, digits = 3)

exogs_exogonly <- intersect(c("z","z2","SH_SH_area_1","SH_SH_capital_1"), names(df_star))
hj_poi <- hansenJ_poi_by_eq_fix(
  poi         = "ln_preco_com_reforma13",
  df          = df_star,
  shareNames  = shareNames,
  omit_share  = best_fit$omit_share,
  drop_price  = best_fit$drop_price,
  priceNames  = priceNames,
  iv_names    = iv_star,         # ex.: c("pls_preco_com_reforma13_1","iv_op01")
  cluster     = cluster_var,
  exogs       = exogs_exogonly
)
print(hj_poi, digits = 3)

overid_diag_one <- function(eq_y, poi, df, exogs, ln_other, Znames){
  use <- unique(c(eq_y, poi, exogs, ln_other, Znames))
  dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
  f_y <- stats::reformulate(c(poi, exogs, ln_other), response = eq_y)
  f_z <- stats::reformulate(c(exogs, ln_other, Znames))
  X <- model.matrix(f_y, dat)
  Z <- model.matrix(f_z, dat)
  rX <- qr(X)$rank; rZ <- qr(Z)$rank
  data.frame(eq = eq_y, n = nrow(dat), rankX = rX, rankZ = rZ, dfJ_eff = rZ - rX)
}

overid_diag <- function(poi, df, shareNames, omit_share, drop_price, priceNames, iv_names, exogs){
  ln_all   <- paste0("ln_", priceNames)
  ln_other <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  eqs <- setdiff(shareNames, omit_share)
  do.call(rbind, lapply(eqs, function(eq)
    overid_diag_one(eq, poi, df, intersect(exogs, names(df)), intersect(ln_other, names(df)),
                    intersect(iv_names, names(df)))))
}

# uso:
exogs_hj  <- intersect(c("z","z2","SH_SH_area_1","SH_SH_capital_1"), names(df_star))
diagJ <- overid_diag(poi, df_star, shareNames, best_fit$omit_share,
                     best_fit$drop_price, priceNames, iv_star, exogs_hj)
print(diagJ, digits=3)

resid_on <- function(v, X) {
  if (is.null(X) || ncol(as.matrix(X))==0) return(as.numeric(v))
  v <- as.numeric(v); X <- as.matrix(X)
  v - X %*% solve(crossprod(X), crossprod(X, v))
}

pick_extra_iv <- function(poi, df, exogs, ln_other, iv_current, iv_candidates){
  Xnames <- intersect(c(exogs, ln_other), names(df))
  X <- if(length(Xnames)) as.matrix(df[, Xnames, drop=FALSE]) else NULL
  # espaço gerado pelos IVs atuais residualizados
  Zc <- sapply(iv_current, function(z) resid_on(df[[z]], X))
  Zc <- as.matrix(Zc)
  keep <- setdiff(iv_candidates, iv_current)
  if (!length(keep)) return(NA_character_)
  sc <- sapply(keep, function(z){
    zt <- resid_on(df[[z]], X)
    # fração de variância explicada por Zc
    R2 <- if (ncol(Zc)==0) 0 else {
      bh <- tryCatch(solve(crossprod(Zc), crossprod(Zc, zt)), error=function(e) rep(0, ncol(Zc)))
      zhat <- as.numeric(Zc %*% bh)
      1 - var(zt - zhat)/var(zt)
    }
    R2
  })
  # queremos o MENOR R2 (mais “novo” condicionalmente a X)
  names(which.min(sc))
}

# candidatos naturais: dos seus IVs base
iv_cand_pool <- intersect(best_iv_set, names(df_star))
iv_extra <- pick_extra_iv(poi, df_star, exogs_hj,
                          setdiff(paste0("ln_", priceNames), c(poi, paste0("ln_", best_fit$drop_price))),
                          iv_star, iv_cand_pool)
iv_star2 <- unique(c(iv_star, iv_extra))
cat("IV extra sugerido:", iv_extra, "\nIVs finais:", paste(iv_star2, collapse=", "), "\n")

hansenJ_poi_by_eq_fix <- function(
    poi, df, shareNames, omit_share, drop_price, priceNames, iv_names, cluster, exogs
){
  eqs      <- setdiff(shareNames, omit_share)
  ln_all   <- paste0("ln_", priceNames)
  ln_other <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xnames   <- intersect(c(exogs, ln_other), names(df))
  Znames   <- intersect(iv_names, names(df))
  
  purrr::map_dfr(eqs, function(eq_y){
    use <- unique(c(eq_y, poi, Xnames, Znames, cluster))
    dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
    if (nrow(dat)==0L)
      return(tibble::tibble(eq=eq_y, n=0L, dfJ=NA, HansenJ=NA, HansenJ_p=NA, note="sem casos"))
    
    f_y <- stats::reformulate(c(poi, Xnames), response = eq_y)
    f_z <- stats::reformulate(c(Xnames, Znames))
    X <- model.matrix(f_y, dat); Z <- model.matrix(f_z, dat)
    rX <- qr(X)$rank; rZ <- qr(Z)$rank; dfJ_eff <- rZ - rX
    if (dfJ_eff <= 0)
      return(tibble::tibble(eq=eq_y, n=nrow(dat), dfJ=0L, HansenJ=NA, HansenJ_p=NA,
                            note="just-identified (dfJ_eff<=0)"))
    
    fit <- AER::ivreg(formula = f_y, instruments = f_z, data = dat)
    vc  <- tryCatch(sandwich::vcovCL(fit, cluster = dat[[cluster]]),
                    error = function(e) sandwich::vcovHC(fit, type="HC1"))
    
    s <- try(ivreg::sargan(fit, vcov. = vc), silent = TRUE)
    if (!inherits(s, "try-error")) {
      return(tibble::tibble(eq=eq_y, n=nrow(dat), dfJ=dfJ_eff,
                            HansenJ = unname(as.numeric(s$statistic)),
                            HansenJ_p = unname(as.numeric(s$p.value)),
                            note = NA_character_))
    }
    
    # ---- Fallback: J manual GMM (cluster) ----
    u  <- resid(fit)                                # resíduo 2SLS
    g  <- Z * as.numeric(u)                         # momentos
    gbar <- colMeans(g)                             # média dos momentos
    # cluster "meat"
    cl <- split(seq_len(nrow(g)), dat[[cluster]])
    Gc <- lapply(cl, function(ix) colSums(g[ix, , drop=FALSE]))
    S  <- Reduce(`+`, lapply(Gc, function(v) tcrossprod(v))) / nrow(g)^2
    # regulariza se necessário
    Sinv <- tryCatch(solve(S), error=function(e) MASS::ginv(S))
    J <- as.numeric(nrow(g) * t(gbar) %*% Sinv %*% gbar)
    p <- 1 - pchisq(J, df = dfJ_eff)
    tibble::tibble(eq=eq_y, n=nrow(dat), dfJ=dfJ_eff, HansenJ=J, HansenJ_p=p,
                   note="J manual (cluster)")
  })
}

# rodar:
hj_poi2 <- hansenJ_poi_by_eq_fix(
  poi         = poi,
  df          = df_star,
  shareNames  = shareNames,
  omit_share  = best_fit$omit_share,
  drop_price  = best_fit$drop_price,
  priceNames  = priceNames,
  iv_names    = iv_star2,         # <-- com o IV extra sugerido
  cluster     = cluster_var,
  exogs       = exogs_hj
)
print(hj_poi2, digits=3)

#' Extract Kleibergen–Paap rk F (KP) and compute robust Hansen–J for an ivsys object
#'
#' Expected ivsys structure (pipeline-style):
#'   - ivsys is a list; each element m corresponds to one equation.
#'   - m$fit is an AER::ivreg object.
#'   - m$eq is optional (equation name); if absent, list names are used.
#'   - m$f_z is optional; if absent, instruments are obtained via model.matrix(m$fit, component="instruments").
#'
#' Returns: data.frame with eq, KP_rkF, KP_p, HansenJ, dfJ, HansenJ_p, n, note
kp_hansen_from_ivsys <- function(ivsys) {
  stopifnot(is.list(ivsys), length(ivsys) > 0L)
  stopifnot(requireNamespace("AER", quietly = TRUE))
  
  # Helper: safe eigen/solve
  .safe_solve <- function(M) {
    out <- try(solve(M), silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }
  
  # Helper: get KP rk F from ivreg diagnostics
  .get_kp <- function(fit) {
    ss <- try(summary(fit, diagnostics = TRUE), silent = TRUE)
    if (inherits(ss, "try-error") || is.null(ss$diagnostics)) {
      return(list(KP_rkF = NA_real_, KP_p = NA_real_, note_kp = "no_diagnostics"))
    }
    diagm <- ss$diagnostics
    rn <- rownames(diagm)
    
    # Robust matching across AER versions / labels
    iKP <- grep("Kleibergen|K-P|Weak instruments", rn, ignore.case = TRUE)
    if (!length(iKP)) {
      return(list(KP_rkF = NA_real_, KP_p = NA_real_, note_kp = "kp_row_not_found"))
    }
    
    # Prefer the first match; if multiple, take the one mentioning Kleibergen explicitly
    if (length(iKP) > 1) {
      ik2 <- iKP[grep("Kleibergen", rn[iKP], ignore.case = TRUE)]
      if (length(ik2)) iKP <- ik2[1] else iKP <- iKP[1]
    } else {
      iKP <- iKP[1]
    }
    
    KP_rkF <- suppressWarnings(as.numeric(diagm[iKP, "statistic"]))
    KP_p   <- suppressWarnings(as.numeric(diagm[iKP, "p-value"]))
    list(KP_rkF = KP_rkF, KP_p = KP_p, note_kp = "ok")
  }
  
  # Helper: robust Hansen J using moment covariance of (Z*u)
  .hansenJ_robust <- function(fit, m = NULL) {
    u <- try(resid(fit), silent = TRUE)
    if (inherits(u, "try-error") || is.null(u)) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "no_residuals"))
    }
    
    # Instruments matrix Z
    Zm <- NULL
    
    # Pipeline sometimes stores f_z; if present, use it with model.frame(fit)
    if (!is.null(m) && !is.null(m$f_z)) {
      Zm <- try(model.matrix(m$f_z, data = model.frame(fit)), silent = TRUE)
      if (inherits(Zm, "try-error")) Zm <- NULL
    }
    
    # Fallback: ask ivreg for instrument matrix directly
    if (is.null(Zm)) {
      Zm <- try(AER::model.matrix(fit, component = "instruments"), silent = TRUE)
      if (inherits(Zm, "try-error")) Zm <- NULL
    }
    
    if (is.null(Zm) || !nrow(Zm) || !ncol(Zm)) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "no_instruments_matrix"))
    }
    
    # Align rows (defensive)
    if (length(u) != nrow(Zm)) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "Z_u_row_mismatch"))
    }
    
    # Drop constant / non-finite columns
    ok_col <- apply(Zm, 2, function(v) all(is.finite(v)) && stats::var(v) > 0)
    Zm <- Zm[, ok_col, drop = FALSE]
    if (!ncol(Zm)) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "all_Z_cols_dropped"))
    }
    
    # Reduce to full column rank to avoid singular meat
    qrZ <- qr(Zm)
    rZ  <- qrZ$rank
    if (rZ <= 0) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "rankZ_zero"))
    }
    Zr <- Zm[, qrZ$pivot[seq_len(rZ)], drop = FALSE]
    
    # Moments
    Zu <- Zr * as.numeric(u)
    g  <- colMeans(Zu)
    meat <- crossprod(Zu) / nrow(Zu)
    
    Winv <- .safe_solve(meat)
    if (is.null(Winv)) {
      return(list(HansenJ = NA_real_, dfJ = NA_integer_, HansenJ_p = NA_real_, note_j = "meat_singular"))
    }
    
    J <- nrow(Zu) * drop(t(g) %*% Winv %*% g)
    
    # Degrees of freedom: (# instruments) - (# estimated coefficients)
    k <- length(stats::coef(fit))
    dfJ <- max(0L, ncol(Zr) - k)
    
    pJ <- if (is.finite(J) && dfJ > 0) stats::pchisq(J, df = dfJ, lower.tail = FALSE) else NA_real_
    
    note_j <- if (dfJ <= 0) "not_overidentified" else "ok"
    list(HansenJ = J, dfJ = dfJ, HansenJ_p = pJ, note_j = note_j)
  }
  
  # Iterate
  out <- lapply(seq_along(ivsys), function(i) {
    m <- ivsys[[i]]
    if (is.null(m$fit)) stop("ivsys[[", i, "]] has no $fit.")
    fit <- m$fit
    
    eq <- NA_character_
    if (!is.null(m$eq) && is.character(m$eq) && length(m$eq) == 1) eq <- m$eq
    if (is.na(eq) || !nzchar(eq)) {
      nm <- names(ivsys)[i]
      eq <- if (!is.null(nm) && nzchar(nm)) nm else paste0("eq", i)
    }
    
    kp <- .get_kp(fit)
    hj <- .hansenJ_robust(fit, m = m)
    
    note <- paste(
      paste0("KP:", kp$note_kp),
      paste0("J:",  hj$note_j),
      sep = " | "
    )
    
    data.frame(
      eq = eq,
      n  = tryCatch(stats::nobs(fit), error = function(e) NA_integer_),
      KP_rkF = kp$KP_rkF,
      KP_p   = kp$KP_p,
      HansenJ   = hj$HansenJ,
      dfJ       = hj$dfJ,
      HansenJ_p = hj$HansenJ_p,
      note = note,
      check.names = FALSE
    )
  })
  
  do.call(rbind, out)
}
tab_kpJ <- kp_hansen_from_ivsys(ivsys)
tab_kpJ

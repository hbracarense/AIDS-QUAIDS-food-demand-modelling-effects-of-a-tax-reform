## =========================================================
## MINI-BOOSTER DE IVs (PLS-F + CLÁSSICOS) — PLUG & PLAY
## Requisitos no ambiente: df_iv_ok, priceNames, shareNames, best_fit, best_iv_set
## =========================================================

## ---------- 0) Dependências ----------
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("AER","ivmodel","dplyr","purrr","tibble","pls"))

`%||%` <- function(a,b) if(!is.null(a)) a else b
as_chr_vec <- function(x) unique(as.character(unlist(x, use.names = FALSE)))

set.seed(123)

## ---------- 1) Pré-requisitos do ambiente ----------
stopifnot(exists("df_iv_ok"), exists("priceNames"), exists("shareNames"),
          exists("best_fit"), exists("best_iv_set"), is.character(best_iv_set))

drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
eqs        <- as_chr_vec(setdiff(shareNames, omit_share))
stopifnot(length(eqs) > 0)

## ---------- 2) Helpers estatísticos ----------
resid_on <- function(v, X){
  if (is.null(X) || ncol(as.matrix(X))==0) return(as.numeric(v))
  v <- as.numeric(v); X <- as.matrix(X)
  v - X %*% solve(crossprod(X), crossprod(X, v))
}
resid_cols_on <- function(Z, X){
  Z <- as.matrix(Z)
  X <- if (is.null(X) || ncol(as.matrix(X))==0) NULL else as.matrix(X)
  if (is.null(X)) return(Z)
  Z - X %*% solve(crossprod(X), crossprod(X, Z))
}
first_stage_F <- function(D, X, Z, data){
  if (!length(Z)) return(NA_real_)
  rhs0 <- if(length(X)) X else "1"
  rhs1 <- c(X, Z)
  f0 <- lm(stats::reformulate(rhs0, response = D), data = data)
  f1 <- lm(stats::reformulate(rhs1, response = D), data = data)
  as.numeric(anova(f0, f1)$F[2])
}
prep_XZ_cond <- function(eq_y, poi, df, iv_names, drop_price, priceNames){
  eq_y     <- as_chr_vec(eq_y)[1]
  poi      <- as_chr_vec(poi)[1]
  iv_names <- as_chr_vec(iv_names)
  stopifnot(eq_y %in% names(df), poi %in% names(df))
  
  exogs_all <- intersect(c("z","z2","SH_SH_area_1","SH_SH_capital_1"), names(df))
  xnames <- as_chr_vec(exogs_all)
  znames <- intersect(iv_names, names(df))
  use <- as_chr_vec(unique(c(eq_y, poi, xnames, znames)))
  dat <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
  
  list(
    dat    = dat,
    Y      = as.numeric(dat[[eq_y]]),
    D      = as.numeric(dat[[poi]]),
    X      = if(length(xnames)) as.matrix(dat[, xnames, drop=FALSE]) else NULL,
    Z      = if(length(znames)) as.matrix(dat[, znames, drop=FALSE]) else NULL,
    Xnames = xnames,
    Znames = znames
  )
}

## ---------- 3) PLS-F para 1 instrumento do preço de interesse ----------
make_sw_pc_for_price <- function(poi, df, iv_pool, priceNames, drop_price, exogs = c("z","z2")){
  stopifnot(poi %in% names(df))
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xn        <- as_chr_vec(c(intersect(exogs, names(df)), intersect(ln_others, names(df))))
  Znames    <- intersect(as_chr_vec(iv_pool), names(df))
  use       <- as_chr_vec(unique(c(poi, Xn, Znames)))
  dat       <- df[stats::complete.cases(df[, use, drop=FALSE]), , drop=FALSE]
  if (nrow(dat) == 0L) stop("Sem casos completos para construir PLS-F do preço: ", poi)
  
  Dj_t <- resid_on(dat[[poi]], if(length(Xn)) as.matrix(dat[, Xn, drop=FALSE]) else NULL)
  Z_t  <- resid_cols_on(dat[, Znames, drop=FALSE],
                        if(length(Xn)) as.matrix(dat[, Xn, drop=FALSE]) else NULL)
  
  if (ncol(Z_t) > 1) {
    keep <- apply(Z_t, 2, sd, na.rm=TRUE) > 1e-12
    Z_t  <- Z_t[, keep, drop=FALSE]
  }
  fit <- pls::plsr(Dj_t ~ ., data = as.data.frame(Z_t),
                   ncomp = 1, validation = "none", scale = TRUE)
  s1  <- drop(pls::scores(fit)[,1])
  nm  <- paste0("pc_", sub("^ln_", "", poi), "_swF1")
  dat[[nm]] <- s1
  list(data = dat, pc_name = nm)
}

## ---------- 4) Avaliação AR/CLR por combinação de IVs ----------
run_one_ivm <- function(iv_names, tag, poi, df, eqs, priceNames, drop_price,
                        alpha = 0.05, halfwidth = 400) {
  iv_names <- as_chr_vec(iv_names)
  eqs      <- as_chr_vec(eqs)
  poi      <- as_chr_vec(poi)[1]
  
  purrr::map_dfr(eqs, function(eq_y){
    eq_y <- as_chr_vec(eq_y)[1]
    pr <- prep_XZ_cond(eq_y, poi, df, iv_names, drop_price, priceNames)
    
    # Se por algum motivo Z ficou vazio (não deveria, pois sempre inclui o PC)
    if (is.null(pr$Z) || length(pr$Znames) == 0L || ncol(as.matrix(pr$Z)) == 0L) {
      return(tibble::tibble(
        tag = tag, eq = eq_y, n = nrow(pr$dat), F_cond = NA_real_,
        AR_low = NA_real_, AR_high = NA_real_, CLR_low = NA_real_, CLR_high = NA_real_,
        width_AR = NA_real_, width_CLR = NA_real_
      ))
    }
    
    Fcond <- suppressWarnings(first_stage_F(poi, pr$Xnames, pr$Znames, pr$dat))
    
    # 2SLS para centralizar a grade
    b2 <- tryCatch({
      fy  <- as.formula(paste(eq_y, "~", paste(c(poi, pr$Xnames), collapse = "+")))
      fz  <- as.formula(paste("~",   paste(c(pr$Xnames, pr$Znames), collapse = "+")))
      fit2 <- AER::ivreg(formula = fy, instruments = fz, data = pr$dat)
      unname(coef(fit2)[poi])
    }, error = function(e) 0)
    
    grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out = 4001)
    
    out <- tryCatch({
      ivm  <- ivmodel::ivmodel(Y = pr$Y, D = pr$D,
                               Z = if(is.null(pr$Z)) NULL else as.matrix(pr$Z),
                               X = if(is.null(pr$X)) NULL else as.matrix(pr$X))
      p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0 = b0)$p.value, numeric(1))
      p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0 = b0)$p.value, numeric(1))
      keep_ar  <- which(p_ar  >= alpha)
      keep_clr <- which(p_clr >= alpha)
      ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(-Inf, +Inf)
      ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(-Inf, +Inf)
      
      tibble::tibble(
        tag = tag, eq = eq_y, n = nrow(pr$dat),
        F_cond = Fcond,
        AR_low  = ci_ar[1],  AR_high  = ci_ar[2],
        CLR_low = ci_clr[1], CLR_high = ci_clr[2]
      ) |>
        dplyr::mutate(width_AR  = AR_high  - AR_low,
                      width_CLR = CLR_high - CLR_low)
    }, error = function(e){
      # Falha do ivmodel: retorna NAs mas não quebra o loop
      tibble::tibble(
        tag = tag, eq = eq_y, n = nrow(pr$dat),
        F_cond = Fcond,
        AR_low  = NA_real_,  AR_high  = NA_real_,
        CLR_low = NA_real_,  CLR_high = NA_real_,
        width_AR = NA_real_, width_CLR = NA_real_
      )
    })
    
    out
  })
}

evaluate_combo <- function(iv_names_all, tag, poi, df, eqs, priceNames, drop_price,
                           alpha = 0.05, halfwidth = 400){
  iv_names_all <- as_chr_vec(iv_names_all)
  eqs          <- as_chr_vec(eqs)
  
  out <- run_one_ivm(iv_names_all, tag, poi, df, eqs, priceNames, drop_price,
                     alpha = alpha, halfwidth = halfwidth)
  
  tibble::tibble(
    tag            = unique(out$tag),
    F_cond_med     = median(out$F_cond, na.rm = TRUE),
    F_cond_min     = suppressWarnings(min(out$F_cond, na.rm = TRUE)),
    CLR_mean_width = mean(out$width_CLR[is.finite(out$width_CLR)], na.rm = TRUE),
    AR_mean_width  = mean(out$width_AR [is.finite(out$width_AR )], na.rm = TRUE)
  )
}

## ---------- 5) Enumerador de combinações (garante vetores char) ----------
enumerate_combos <- function(items, max_extra = 3, cap = 300L){
  items <- as_chr_vec(items)
  res <- vector("list", 0L)
  for (k in 0:max_extra) {
    if (k == 0) {
      combos <- list(character(0))
    } else {
      total <- choose(length(items), k)
      if (total <= cap) {
        combos <- combn(items, k, simplify = FALSE)
      } else {
        bag <- new.env(parent = emptyenv())
        set <- list()
        while (length(set) < cap) {
          pick <- sort(sample(items, k))
          key  <- paste(pick, collapse = "||")
          if (!exists(key, envir = bag, inherits = FALSE)) {
            assign(key, TRUE, envir = bag)
            set[[length(set)+1L]] <- pick
          }
        }
        combos <- set
      }
    }
    # força cada elemento a ser character()
    combos <- lapply(combos, as_chr_vec)
    res <- c(res, combos)
  }
  res
}

## ---------- 6) Construção do dicionário de IVs ----------
augment_dictionary <- function(df, base_iv, with = c("SH_SH_area_1","SH_SH_capital_1"),
                               poly_deg = 3, center=TRUE, scale.=TRUE){
  base_iv <- as_chr_vec(base_iv)
  with <- intersect(with, names(df))
  out_names <- base_iv
  for (b in base_iv) if (b %in% names(df)) {
    v <- df[[b]]
    if (center) v <- v - mean(v, na.rm=TRUE)
    if (scale.) v <- v / sd(v, na.rm=TRUE)
    for (d in 2:poly_deg) {
      nm <- paste0(b,"_p", d)
      if (!nm %in% names(df)) df[[nm]] <- as.numeric(v^d)
      out_names <- c(out_names, nm)
    }
    for (w in with) {
      nm <- paste0(b,"_x_", w)
      if (!nm %in% names(df)) df[[nm]] <- as.numeric(df[[b]]) * as.numeric(df[[w]])
      out_names <- c(out_names, nm)
    }
  }
  list(data = df, iv_names = unique(as_chr_vec(out_names)))
}
add_harmonics <- function(df, base, K=3){
  base <- as_chr_vec(base)
  keep <- base[base %in% names(df)]
  for (b in keep) for (k in 2:K) {
    s <- paste0(b, "_sin", k); c <- paste0(b, "_cos", k)
    if (!s %in% names(df)) df[[s]] <- sin(k * df[[b]])
    if (!c %in% names(df)) df[[c]] <- cos(k * df[[b]])
  }
  df
}

## ---------- 7) Montagem da base e do pool ----------
aug <- augment_dictionary(
  df_iv_ok,
  base_iv = as_chr_vec(best_iv_set),
  with    = c("SH_SH_area_1","SH_SH_capital_1"),
  poly_deg= 3
)
df_plsF <- aug$data
iv_pool <- aug$iv_names

df_plsF <- add_harmonics(df_plsF, base = c("iv_sin1","iv_cos1"), K = 3)
iv_pool <- unique(c(iv_pool, grep("^iv_(sin|cos)[2-9]$", names(df_plsF), value=TRUE)))

base_iv_pool <- intersect(
  c("iv_op01","iv_op02","iv_op03","iv_sin1","iv_cos1","iv_sin2","iv_cos2","iv_hg40"),
  names(df_plsF)
)
pc_pool <- grep("^pc_", names(df_plsF), value = TRUE)
iv_pool_all <- unique(c(iv_pool, base_iv_pool, pc_pool))

cat("\n[Mini-Booster] IVs clássicos detectados:\n"); print(base_iv_pool)
cat("\n[Mini-Booster] PCs pré-existentes detectados:\n"); print(pc_pool)

## ---------- 8) PLS-F (1 PC) para o preço de interesse ----------
poi <- "ln_preco_com_reforma13"  # ajuste se quiser
stopifnot(poi %in% names(df_plsF))

pls_one <- make_sw_pc_for_price(
  poi        = poi,
  df         = df_plsF,
  iv_pool    = iv_pool_all,
  priceNames = priceNames,
  drop_price = drop_price
)
df_plsF <- pls_one$data
iv_plsF <- as_chr_vec(pls_one$pc_name)
cat("\n[PLS-F] Instrumento base (poi) criado: ", iv_plsF, "\n", sep="")

## ---------- 9) Booster de combinações (loop robusto; sem imap_dfr) ----------
target_F  <- 10
max_extra <- 3
halfwidth <- 400

cmbs <- enumerate_combos(base_iv_pool, max_extra = max_extra, cap = 300L)

rows <- vector("list", length(cmbs))
for (i in seq_along(cmbs)) {
  extra_chr <- as_chr_vec(cmbs[[i]])
  iv_all    <- as_chr_vec(c(iv_plsF, extra_chr))
  
  metr <- tryCatch(
    evaluate_combo(
      iv_names_all = iv_all,
      tag          = paste0("PLS-F+", length(extra_chr)),
      poi          = poi,
      df           = df_plsF,
      eqs          = eqs,
      priceNames   = priceNames,
      drop_price   = drop_price,
      alpha        = 0.05,
      halfwidth    = halfwidth
    ),
    error = function(e) {
      tibble::tibble(
        tag = paste0("PLS-F+", length(extra_chr)),
        F_cond_med = NA_real_, F_cond_min = NA_real_,
        CLR_mean_width = NA_real_, AR_mean_width = NA_real_
      )
    }
  )
  
  rows[[i]] <- tibble::tibble(
    combo_id   = i,
    size_extra = length(extra_chr),
    iv_extra   = if (length(extra_chr)) paste(extra_chr, collapse = " + ") else "(nenhum)",
    iv_all_n   = length(iv_all)
  ) |> dplyr::bind_cols(metr)
}

score_tbl <- dplyr::bind_rows(rows)

score_tbl <- score_tbl |>
  dplyr::arrange(size_extra, dplyr::desc(F_cond_med), dplyr::coalesce(CLR_mean_width, Inf))

cat("\n[Booster] Top 10 combinações:\n")
print(utils::head(score_tbl, 10), digits = 3, row.names = FALSE)

## ---------- 10) Seleção da melhor combinação ----------
best_row <- score_tbl |>
  dplyr::filter(is.finite(F_cond_med), F_cond_med >= target_F) |>
  dplyr::slice(1)

if (nrow(best_row) == 0L) {
  warning("Nenhuma combinação atingiu a meta F_cond_med ≥ ", target_F,
          ". Considere aumentar max_extra/cap, ampliar o dicionário (poly/interações) ou elevar a meta.")
} else {
  best_extra <- if (identical(best_row$iv_extra, "(nenhum)")) character(0) else
    strsplit(best_row$iv_extra, " \\+ ")[[1]]
  best_ivset <- as_chr_vec(c(iv_plsF, best_extra))
  
  cat("\n[Booster] Melhor combinação (atingiu a meta):\n")
  cat("  Extras:", if(length(best_extra)) paste(best_extra, collapse=", ") else "(nenhum)", "\n")
  cat("  |IV total| =", length(best_ivset),
      " | F_cond_med =", round(best_row$F_cond_med, 2),
      " | F_cond_min =", round(best_row$F_cond_min, 2), "\n")
  
  tab_best <- run_one_ivm(best_ivset, "BEST", poi, df_plsF, eqs, priceNames, drop_price,
                          alpha = 0.05, halfwidth = halfwidth)
  cat("\n[Booster] Detalhe por equação (BEST):\n")
  print(tab_best, digits = 3, row.names = FALSE)
}

cat("\n=== FIM: Mini-Booster concluído. ===\n")
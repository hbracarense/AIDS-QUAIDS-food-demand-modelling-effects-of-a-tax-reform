# ============================================================
# Passos 1–5: 2SLS ponderado + sazonais + (z,z2) ou spline(z)
#             + 3SLS (adding-up) com defesas contra singularidade
#             + diagnósticos/plots
# Dataset: dat  (usa nomes fornecidos pelo usuário)
# ============================================================

suppressPackageStartupMessages({
  library(AER)        # ivreg
  library(systemfit)  # 3SLS
  library(splines)    # bs()
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(sandwich)
})

options(na.action = na.exclude)

# ========================
# CONFIGURAÇÃO
# ========================

share_vars  <- paste0("w_despesahat", 1:6)
price_vars  <- paste0("ln_preco_com_reforma2", 1:6)
weights_col <- "peso_final"

seasonal_controls   <- c(paste0("iv_sin", 1:4), paste0("iv_cos", 1:4))
iv_include_patterns <- c("^iv_", "^IV_")
iv_exclude_patterns <- c("^iv_sin", "^iv_cos")

z_var   <- "z"
z2_var  <- "z2"   # já existe em dat
use_spline <- FALSE
spline_df  <- 4

# ============================================================
# HELPERS
# ============================================================

get_cols_by_regex <- function(dat, include_patterns = character(), exclude_patterns = character()) {
  vars <- colnames(dat)
  keep <- if (length(include_patterns)) Reduce(`|`, lapply(include_patterns, grepl, x = vars)) else rep(TRUE, length(vars))
  drop <- if (length(exclude_patterns)) Reduce(`|`, lapply(exclude_patterns, grepl, x = vars)) else rep(FALSE, length(vars))
  vars[keep & !drop]
}

make_rhs <- function(price_vars,
                     z_var = "z",
                     z2_var = "z2",
                     seasonal_controls = character(),
                     use_spline = FALSE,
                     spline_df = 4,
                     centered = FALSE) {
  z_name  <- if (centered) "z_c"  else z_var
  z2_name <- if (centered) "z2_c" else z2_var
  z_terms <- if (isTRUE(use_spline)) sprintf("splines::bs(%s, df=%d)", z_name, spline_df) else c(z_name, z2_name)
  rhs_vec <- unique(c(price_vars, z_terms, seasonal_controls))
  paste(rhs_vec, collapse = " + ")
}

prune_instruments <- function(dat, iv_candidates, tol = 1e-7) {
  stopifnot(length(iv_candidates) > 0)
  f_inst_full <- as.formula(paste("~ -1 +", paste(iv_candidates, collapse = " + ")))
  mm <- model.matrix(f_inst_full, dat)
  if (ncol(mm) == 0) stop("Todos os IVs viraram constantes após model.matrix().")
  keep_var <- apply(mm, 2, function(x) isTRUE(sd(x, na.rm = TRUE) > 0))
  mm <- mm[, keep_var, drop = FALSE]
  if (ncol(mm) == 0) stop("Nenhum IV com variância positiva após limpeza.")
  q <- qr(mm, tol = tol)
  if (q$rank == 0) stop("Rank 0 no bloco de instrumentos.")
  keep_idx <- sort(q$pivot[seq_len(q$rank)])
  iv_keep <- colnames(mm)[keep_idx]
  list(
    iv_keep = iv_keep,
    dropped = setdiff(colnames(mm), iv_keep),
    inst_formula = as.formula(paste("~ -1 +", paste(iv_keep, collapse = " + ")))
  )
}

prune_exog_for_rhs <- function(dat, seasonal_controls, z_names, tol = 1e-7) {
  exog_vec <- unique(c(seasonal_controls, z_names))
  if (!length(exog_vec)) return(character(0))
  f_exog <- as.formula(paste("~ -1 +", paste(exog_vec, collapse = " + ")))
  mm <- model.matrix(f_exog, dat)
  if (ncol(mm) == 0) return(character(0))
  keep_var <- apply(mm, 2, function(x) isTRUE(sd(x, na.rm = TRUE) > 0))
  mm <- mm[, keep_var, drop = FALSE]
  if (ncol(mm) == 0) return(character(0))
  q <- qr(mm, tol = tol)
  keep_idx <- sort(q$pivot[seq_len(q$rank)])
  colnames(mm)[keep_idx]
}

w_rmse <- function(res, w) sqrt(weighted.mean(res^2, w))
w_mae  <- function(res, w) weighted.mean(abs(res), w)

calibration_df <- function(fitted, resid, w, n_bins = 10) {
  tibble(fitted = fitted, resid = resid, w = w) |>
    mutate(bin = ntile(fitted, n_bins)) |>
    group_by(bin) |>
    summarise(
      n           = n(),
      res_medio   = mean(resid),
      res_medio_w = weighted.mean(resid, w),
      .groups = "drop"
    )
}

# ============================================================
# (1–3) 2SLS ponderado por equação
# ============================================================

fit_iv_weighted <- function(dat, dep, price_vars, weights_col,
                            iv_include_patterns, iv_exclude_patterns,
                            z_var = "z", z2_var = "z2",
                            seasonal_controls = character(),
                            use_spline = FALSE, spline_df = 4,
                            iv_qr_tol = 1e-7) {
  
  stopifnot(all(c(dep, price_vars, weights_col, z_var) %in% names(dat)))
  if (!isTRUE(use_spline)) stopifnot(z2_var %in% names(dat))
  
  if (!all(seasonal_controls %in% names(dat))) {
    seasonal_controls <- intersect(seasonal_controls, names(dat))
    if (!length(seasonal_controls)) message("Sem controles sazonais válidos no data.")
  }
  
  rhs_str <- make_rhs(price_vars = price_vars, z_var = z_var, z2_var = z2_var,
                      seasonal_controls = seasonal_controls,
                      use_spline = use_spline, spline_df = spline_df)
  
  iv_external <- get_cols_by_regex(dat, iv_include_patterns, iv_exclude_patterns)
  if (length(iv_external) == 0) stop("Nenhum instrumento candidato após filtros include/exclude.")
  
  exog_simple <- unique(c(seasonal_controls, z_var, if (!isTRUE(use_spline)) z2_var))
  inst_candidates <- unique(c(iv_external, exog_simple))
  
  pr <- prune_instruments(dat, inst_candidates, tol = iv_qr_tol)
  inst_str <- paste(pr$iv_keep, collapse = " + ")
  
  f_full <- as.formula(paste0(dep, " ~ ", rhs_str, " | -1 + ", inst_str))
  
  w <- dat[[weights_col]]
  fit <- AER::ivreg(f_full, data = dat, weights = w)
  
  vc <- sandwich::vcovHC(fit, type = "HC1")
  sm <- summary(fit, vcov. = vc, diagnostics = TRUE)
  
  list(
    fit = fit,
    summary = sm,
    formula = f_full,
    rhs = rhs_str,
    instruments_used = pr$iv_keep,
    instruments_dropped = pr$dropped,
    weights = w
  )
}

# ============================================================
# (4) 3SLS com adding-up (funções internas)
# ============================================================

.fit_3sls_once <- function(dat, dep_vec, price_vars,
                           iv_include_patterns, iv_exclude_patterns,
                           seasonal_controls,
                           weights_col,
                           z_var, use_spline, spline_df,
                           drop_equation,
                           iv_qr_tol_inst = 1e-6,
                           exog_qr_tol = 1e-6) {
  
  keep_deps <- setdiff(dep_vec, drop_equation)
  
  dat_sys <- dat
  center_z <- mean(dat_sys[[z_var]], na.rm = TRUE)
  zc <- dat_sys[[z_var]] - center_z
  dat_sys$z_c  <- zc
  dat_sys$z2_c <- zc^2
  
  seasonal_controls <- intersect(seasonal_controls, names(dat_sys))
  
  z_names <- if (isTRUE(use_spline)) "z_c" else c("z_c", "z2_c")
  exog_keep <- if (!isTRUE(use_spline)) {
    prune_exog_for_rhs(dat_sys, seasonal_controls, z_names, tol = exog_qr_tol)
  } else {
    seasonal_controls
  }
  
  rhs_str <- make_rhs(price_vars,
                      z_var = z_var, z2_var = z2_var,
                      seasonal_controls = exog_keep,
                      use_spline = use_spline, spline_df = spline_df,
                      centered = TRUE)
  
  eq_list <- setNames(
    lapply(keep_deps, function(d) as.formula(sprintf("%s ~ %s", d, rhs_str))),
    paste0("w", seq_along(keep_deps))
  )
  map <- tibble(dep = keep_deps, label = names(eq_list))
  
  iv_external <- get_cols_by_regex(dat_sys, iv_include_patterns, iv_exclude_patterns)
  if (length(iv_external) == 0) stop("Nenhum instrumento candidato após filtros include/exclude.")
  
  exog_simple <- unique(c(exog_keep, "z_c", if (!isTRUE(use_spline)) "z2_c"))
  inst_candidates <- unique(c(iv_external, exog_simple))
  
  pr <- prune_instruments(dat_sys, inst_candidates, tol = iv_qr_tol_inst)
  inst_formula <- pr$inst_formula
  
  message(sprintf("[3SLS] IVs candidatos: %d | após poda (QR): %d | descartados: %d | exogs no RHS: %d",
                  length(inst_candidates), length(pr$iv_keep), length(pr$dropped), length(exog_keep)))
  
  inst_list <- setNames(replicate(length(eq_list), inst_formula, simplify = FALSE),
                        names(eq_list))
  
  weightlist <- NULL
  if (!is.null(weights_col) && weights_col %in% names(dat_sys)) {
    wv <- dat_sys[[weights_col]]
    weightlist <- setNames(replicate(length(eq_list), wv, simplify = FALSE),
                           names(eq_list))
  }
  
  ctrl <- systemfit::systemfit.control(maxit = 200, tol = 1e-7)
  fit3 <- systemfit::systemfit(eq_list,
                               method = "3SLS",
                               inst = inst_list,
                               data = dat_sys,
                               control = ctrl,
                               weights = weightlist)
  
  list(
    fit = fit3,
    rhs = rhs_str,
    instruments_used = pr$iv_keep,
    instruments_dropped = pr$dropped,
    dropped = drop_equation,
    map = map,
    center_z = center_z,
    z_var = z_var,
    use_spline = use_spline,
    spline_df = spline_df
  )
}

fit_3sls_system <- function(dat, dep_vec, price_vars, iv_include_patterns, iv_exclude_patterns,
                            z_var = "z", z2_var = "z2",
                            seasonal_controls = character(),
                            use_spline = FALSE, spline_df = 4,
                            drop_equation = tail(dep_vec, 1),
                            weights_col = NULL) {
  
  out <- try(
    .fit_3sls_once(dat, dep_vec, price_vars,
                   iv_include_patterns, iv_exclude_patterns,
                   seasonal_controls,
                   weights_col,
                   z_var, use_spline, spline_df,
                   drop_equation,
                   iv_qr_tol_inst = 1e-6,
                   exog_qr_tol = 1e-6),
    silent = TRUE
  )
  if (!inherits(out, "try-error")) return(out)
  
  message("⚠️ 3SLS falhou (singular). Tentando reduzir IVs e reforçar poda...")
  
  iv_inc_fb <- c("^iv_op", "^iv_hg")
  out2 <- try(
    .fit_3sls_once(dat, dep_vec, price_vars,
                   iv_inc_fb, iv_exclude_patterns,
                   seasonal_controls,
                   weights_col,
                   z_var, use_spline, spline_df,
                   drop_equation,
                   iv_qr_tol_inst = 1e-4,
                   exog_qr_tol = 1e-5),
    silent = TRUE
  )
  if (!inherits(out2, "try-error")) {
    out2$fallback <- "IVs reduzidos: ^iv_op, ^iv_hg; tolerâncias aumentadas"
    return(out2)
  }
  
  message("⚠️ Ainda singular. Tentando 3º fallback: sem sazonais no RHS.")
  out3 <- try(
    .fit_3sls_once(dat, dep_vec, price_vars,
                   iv_inc_fb, c(iv_exclude_patterns, "^iv_sin", "^iv_cos"),
                   seasonal_controls = character(0),
                   weights_col,
                   z_var, use_spline, spline_df,
                   drop_equation,
                   iv_qr_tol_inst = 1e-4,
                   exog_qr_tol = 1e-5),
    silent = TRUE
  )
  if (!inherits(out3, "try-error")) {
    out3$fallback <- "Sem sazonais no RHS; IVs ^iv_op/^iv_hg; tolerâncias aumentadas"
    return(out3)
  }
  
  stop("3SLS não convergiu: matriz singular mesmo após fallbacks. Tente reduzir IVs/controles.")
}

# ===== PREVISÃO 3SLS: reconstrói fórmula a partir do RHS salvo =====
# ===== PREVISÃO 3SLS (à prova de versão do systemfit) =====
predict_3sls <- function(fit3_obj, newdata, clip01 = FALSE) {
  stopifnot(is.list(fit3_obj), inherits(fit3_obj$fit, "systemfit"))
  nd <- as.data.frame(newdata)
  
  # Recria a centralização usada no ajuste
  nd$z_c  <- nd[[fit3_obj$z_var]] - fit3_obj$center_z
  if (!isTRUE(fit3_obj$use_spline)) nd$z2_c <- nd$z_c^2
  
  k    <- nrow(fit3_obj$map)        # nº de equações estimadas (5)
  cols <- fit3_obj$map$dep          # nomes originais (w_despesahat1..5)
  labs <- fit3_obj$map$label        # rótulos usados no systemfit (w1..w5)
  
  preds <- matrix(NA_real_, nrow = nrow(nd), ncol = k)
  
  for (j in seq_len(k)) {
    # 1) tenta eqn = j (systemfit >= algumas versões)
    pj <- try(stats::predict(fit3_obj$fit, newdata = nd, eqn = j), silent = TRUE)
    
    # 2) tenta equation = j (systemfit mais antigo)
    if (inherits(pj, "try-error")) {
      pj <- try(stats::predict(fit3_obj$fit, newdata = nd, equation = j), silent = TRUE)
    }
    
    # 3) se ainda vier lista ou falhar, pede todas e escolhe j/label
    if (inherits(pj, "try-error") || is.list(pj)) {
      pj_all <- stats::predict(fit3_obj$fit, newdata = nd)
      if (is.list(pj_all)) {
        # prioriza pegar pelo rótulo salvo
        if (!is.null(names(pj_all)) && labs[j] %in% names(pj_all)) {
          pj <- pj_all[[labs[j]]]
        } else {
          pj <- pj_all[[j]]
        }
      } else {
        pj <- pj_all
      }
    }
    
    preds[, j] <- as.numeric(pj)
  }
  
  colnames(preds) <- cols
  out <- as.data.frame(preds)
  
  # Reconstrói a equação dropada pelo adding-up
  out[[fit3_obj$dropped]] <- 1 - rowSums(preds)
  
  # (opcional) clip em [0,1] + renormaliza
  if (isTRUE(clip01)) {
    out[] <- pmax(0, pmin(1, out))
    s <- rowSums(out)
    s[s == 0] <- 1
    out[] <- out / s
  }
  
  tibble::as_tibble(out[, c(cols, fit3_obj$dropped), drop = FALSE])
}

# ============================================================
# (5) Diagnósticos & gráficos
# ============================================================
# sobrescreva a função
# SUBSTITUI a função por uma versão 100% escalar (1 linha por equação)
diagnostics_iv <- function(fit_obj, dep) {
  # usa só as observações efetivamente usadas no ajuste
  res <- stats::residuals(fit_obj)          # comprimento = n_uso
  fit <- stats::fitted(fit_obj)             # idem (não usamos, mas ok p/ checagem)
  mf  <- stats::model.frame(fit_obj)        # só dados usados
  w_used <- stats::model.weights(mf)        # pesos alinhados a res
  
  rmse   <- sqrt(mean(res^2, na.rm = TRUE))
  mae    <- mean(abs(res),    na.rm = TRUE)
  rmse_w <- if (!is.null(w_used)) sqrt(weighted.mean(res^2, w_used, na.rm = TRUE)) else NA_real_
  mae_w  <- if (!is.null(w_used)) weighted.mean(abs(res), w_used, na.rm = TRUE)    else NA_real_
  
  tibble::tibble(
    share  = dep,
    rmse   = rmse,
    mae    = mae,
    rmse_w = rmse_w,
    mae_w  = mae_w
  )
}




plot_resid_hist_grid <- function(resid_list, bins = 20) {
  df <- imap_dfr(resid_list, ~tibble(share = .y, resid = .x))
  ggplot(df, aes(x = resid)) +
    geom_histogram(bins = bins, fill = "grey30") +
    facet_wrap(~ share, scales = "free_y") +
    labs(x = "Resíduo", y = "Contagem", title = "Resíduos (Obs − Pred) — in-sample")
}

plot_calibration <- function(calib_tbl, title = "Calibração (resíduo médio por decis do fitted)") {
  ggplot(calib_tbl, aes(x = bin)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(y = res_medio), linewidth = 0.6) +
    geom_point(aes(y = res_medio), size = 1.8) +
    geom_line(aes(y = res_medio_w), linewidth = 0.6) +
    geom_point(aes(y = res_medio_w), size = 1.8, shape = 21, fill = "white") +
    labs(x = "Decis do fitted", y = "Resíduo médio (simples e ponderado)", title = title)
}

plot_qq <- function(resid, title = "QQ-plot dos resíduos") {
  df <- tibble(sample = resid)
  ggplot(df, aes(sample = sample)) +
    stat_qq() + stat_qq_line() +
    labs(title = title, x = "Quantis teóricos", y = "Quantis amostrais")
}

# ========================
# EXECUÇÃO
# ========================

# (1–3) 2SLS
fits_iv <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls,
  use_spline = use_spline, spline_df = spline_df
))
names(fits_iv) <- share_vars

diag_tbl <- purrr::map2_dfr(fits_iv, share_vars, ~ diagnostics_iv(.x$fit, .y))
print(diag_tbl, n = Inf)

resid_list <- setNames(lapply(fits_iv, \(x) resid(x$fit)), share_vars)
p_hist <- plot_resid_hist_grid(resid_list, bins = 20)
print(p_hist)

for (dep in c("w_despesahat5", "w_despesahat6")) {
  ff <- fits_iv[[dep]]$fit
  
  # vetores alinhados às observações usadas no ajuste
  res    <- as.numeric(stats::residuals(ff))
  fit_v  <- as.numeric(stats::fitted(ff))
  mf     <- stats::model.frame(ff)                 # apenas linhas usadas
  w_used <- stats::model.weights(mf)
  if (is.null(w_used)) w_used <- rep(1, length(res))
  
  stopifnot(length(res) == length(fit_v), length(res) == length(w_used))
  
  calib <- calibration_df(fit_v, res, w_used, n_bins = 10)
  print(plot_calibration(calib, paste0("Calibração — ", dep)))
  print(plot_qq(res,         paste0("QQ-plot — ", dep)))
}


# (4) 3SLS com adding-up
fit_3sls <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls,
  use_spline = use_spline, spline_df = spline_df,
  drop_equation = "w_despesahat6",
  weights_col = weights_col
)
summary(fit_3sls$fit)

# Predições do 3SLS
pred_shares_3sls <- predict_3sls(fit_3sls, dat)
head(pred_shares_3sls)
summary(rowSums(pred_shares_3sls))

# ========================
# DICAS
# ========================
# - Se ainda houver singularidade, reduza manualmente os IVs (ex.: ^iv_op, ^iv_hg).
# - Para spline em z: use_spline <- TRUE (em ambos 2SLS e 3SLS).
# - Salvar gráficos: ggsave("hist_residuos.png", p_hist, width=10, height=6)
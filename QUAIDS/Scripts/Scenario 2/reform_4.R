# ============================================================
# Pipeline 1–6: 2SLS por equação, diagnósticos, 3SLS (adding-up),
#               comparações, testes, validação e reconciliação
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
  library(car)        # linearHypothesis
})

options(na.action = na.exclude)

# ========================
# CONFIGURAÇÃO
# ========================

share_vars  <- paste0("w_despesahat", 1:6)
price_vars  <- paste0("ln_preco_com_reforma2", 1:6)
weights_col <- "peso_final"

seasonal_controls_full <- c(paste0("iv_sin", 1:4), paste0("iv_cos", 1:4))
seasonal_controls_min  <- c("iv_sin1","iv_sin2","iv_cos1","iv_cos2")

# IVs (padrão amplo e versão “mínima”)
iv_include_patterns_full <- c("^iv_", "^IV_")
iv_include_patterns_min  <- c("^iv_op", "^iv_hg")
iv_exclude_patterns      <- c("^iv_sin", "^iv_cos")  # exclui sazonais do bloco externo

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
  z_terms <- if (isTRUE(use_spline)) sprintf("splines::bs(%s, df=%d)", z_name, spline_df) else {
    # permite remover z2 passando z2_var = NA
    if (is.null(z2_var) || is.na(z2_var) || z2_var == "") z_name else c(z_name, z2_name)
  }
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

# (1–3) 2SLS por equação
fit_iv_weighted <- function(dat, dep, price_vars, weights_col,
                            iv_include_patterns, iv_exclude_patterns,
                            z_var = "z", z2_var = "z2",
                            seasonal_controls = character(),
                            use_spline = FALSE, spline_df = 4,
                            allow_no_z2 = FALSE,
                            iv_qr_tol = 1e-7) {
  
  stopifnot(all(c(dep, price_vars, weights_col, z_var) %in% names(dat)))
  if (!isTRUE(use_spline) && !isTRUE(allow_no_z2)) stopifnot(z2_var %in% names(dat))
  
  if (!all(seasonal_controls %in% names(dat))) {
    seasonal_controls <- intersect(seasonal_controls, names(dat))
  }
  
  rhs_str <- make_rhs(price_vars = price_vars, z_var = z_var, z2_var = z2_var,
                      seasonal_controls = seasonal_controls,
                      use_spline = use_spline, spline_df = spline_df)
  
  iv_external <- get_cols_by_regex(dat, iv_include_patterns, iv_exclude_patterns)
  if (length(iv_external) == 0) stop("Nenhum instrumento candidato após filtros include/exclude.")
  
  exog_simple <- unique(c(seasonal_controls, z_var, if (!isTRUE(use_spline) && !isTRUE(allow_no_z2)) z2_var))
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

# (4) 3SLS com adding-up
.fit_3sls_once <- function(dat, dep_vec, price_vars,
                           iv_include_patterns, iv_exclude_patterns,
                           seasonal_controls,
                           weights_col,
                           z_var, use_spline, spline_df,
                           drop_equation,
                           z2_present = TRUE,
                           iv_qr_tol_inst = 1e-6,
                           exog_qr_tol = 1e-6) {
  
  keep_deps <- setdiff(dep_vec, drop_equation)
  
  dat_sys <- dat
  center_z <- mean(dat_sys[[z_var]], na.rm = TRUE)
  zc <- dat_sys[[z_var]] - center_z
  dat_sys$z_c  <- zc
  dat_sys$z2_c <- zc^2
  
  seasonal_controls <- intersect(seasonal_controls, names(dat_sys))
  
  z_names <- if (isTRUE(use_spline)) "z_c" else if (isTRUE(z2_present)) c("z_c", "z2_c") else "z_c"
  exog_keep <- if (!isTRUE(use_spline)) {
    prune_exog_for_rhs(dat_sys, seasonal_controls, z_names, tol = exog_qr_tol)
  } else {
    seasonal_controls
  }
  
  rhs_str <- make_rhs(price_vars,
                      z_var = z_var, z2_var = if (isTRUE(z2_present)) "z2" else NA,
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
  
  exog_simple <- unique(c(exog_keep, "z_c", if (!isTRUE(use_spline) && isTRUE(z2_present)) "z2_c"))
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
                            weights_col = NULL,
                            z2_present = TRUE) {
  
  out <- try(
    .fit_3sls_once(dat, dep_vec, price_vars,
                   iv_include_patterns, iv_exclude_patterns,
                   seasonal_controls,
                   weights_col,
                   z_var, use_spline, spline_df,
                   drop_equation,
                   z2_present = z2_present,
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
                   z2_present = z2_present,
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
                   z2_present = z2_present,
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

# PREVISÃO 3SLS (robusta a versões do systemfit)
predict_3sls <- function(fit3_obj, newdata, clip01 = FALSE) {
  stopifnot(is.list(fit3_obj), inherits(fit3_obj$fit, "systemfit"))
  nd <- as.data.frame(newdata)
  
  nd$z_c  <- nd[[fit3_obj$z_var]] - fit3_obj$center_z
  if (!isTRUE(fit3_obj$use_spline)) nd$z2_c <- nd$z_c^2
  
  k    <- nrow(fit3_obj$map)
  cols <- fit3_obj$map$dep
  labs <- fit3_obj$map$label
  
  preds <- matrix(NA_real_, nrow = nrow(nd), ncol = k)
  for (j in seq_len(k)) {
    pj <- try(stats::predict(fit3_obj$fit, newdata = nd, eqn = j), silent = TRUE)
    if (inherits(pj, "try-error")) {
      pj <- try(stats::predict(fit3_obj$fit, newdata = nd, equation = j), silent = TRUE)
    }
    if (inherits(pj, "try-error") || is.list(pj)) {
      pj_all <- stats::predict(fit3_obj$fit, newdata = nd)
      if (is.list(pj_all)) {
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
  out[[fit3_obj$dropped]] <- 1 - rowSums(preds)
  
  if (isTRUE(clip01)) {
    out[] <- pmax(0, pmin(1, out))
    s <- rowSums(out)
    s[s == 0] <- 1
    out[] <- out / s
  }
  
  tibble::as_tibble(out[, c(cols, fit3_obj$dropped), drop = FALSE])
}

# GRÁFICOS E DIAGNÓSTICOS AUXILIARES
# ---------- helpers seguros para ivreg ----------
get_res_fitted <- function(fit){
  res  <- as.numeric(stats::residuals(fit))
  fitv <- as.numeric(stats::fitted(fit))
  mf   <- stats::model.frame(fit)           # só as linhas usadas
  w    <- stats::model.weights(mf)
  if (is.null(w)) w <- rep(1, length(res))
  stopifnot(length(res) == length(fitv), length(res) == length(w))
  list(res = res, fit = fitv, w = w)
}

w_rmse <- function(res, w) sqrt(weighted.mean(res^2, w))
w_mae  <- function(res, w) weighted.mean(abs(res), w)

diagnostics_iv <- function(fit_obj, dep){
  rf <- get_res_fitted(fit_obj)
  tibble::tibble(
    share  = dep,
    rmse   = sqrt(mean(rf$res^2)),
    mae    = mean(abs(rf$res)),
    rmse_w = w_rmse(rf$res, rf$w),
    mae_w  = w_mae(rf$res, rf$w),
    res    = list(rf$res),   # list-cols para gráficos depois
    fitted = list(rf$fit)
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

# ------------------------------------------------------------
# ======================== EXECUÇÃO 1–6 =======================
# ------------------------------------------------------------

# (1–3) 2SLS por equação + métricas/plots básicos
fits_iv <- map(share_vars, ~ fit_iv_weighted(
  dat = dat, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = use_spline, spline_df = spline_df
))
names(fits_iv) <- share_vars

diag_tbl <- purrr::map2_dfr(fits_iv, names(fits_iv), ~ diagnostics_iv(.x$fit, .y)) |>
  dplyr::select(share, rmse, mae, rmse_w, mae_w)
print(diag_tbl)


resid_list <- setNames(lapply(fits_iv, \(x) resid(x$fit)), share_vars)
print(plot_resid_hist_grid(resid_list, bins = 20))

for (dep in c("w_despesahat5","w_despesahat6")) {
  rf <- get_res_fitted(fits_iv[[dep]]$fit)
  calib <- calibration_df(rf$fit, rf$res, rf$w, n_bins = 10)
  print(plot_calibration(calib, paste0("Calibração — ", dep)))
  print(plot_qq(rf$res, paste0("QQ-plot — ", dep)))
}


# (4) 3SLS com adding-up: dropa a 6ª e reconstrói pela soma
fit_3sls <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = spline_df,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = TRUE
)
print(summary(fit_3sls$fit))

# (4b) Predições 3SLS e comparação com observados (RMSE ponderado)
pred_shares_3sls <- predict_3sls(fit_3sls, dat)  # clip01 = TRUE se quiser forçar [0,1]
obs_mat   <- as.data.frame(dat[, share_vars])
w_vec     <- dat[[weights_col]]

rmse_w_iv <- sapply(share_vars, function(v) {
  res <- resid(fits_iv[[v]]$fit)
  sqrt(weighted.mean(res^2, w_vec))
})

rmse_w_3sls <- sapply(share_vars, function(v) {
  res <- obs_mat[[v]] - pred_shares_3sls[[v]]
  sqrt(weighted.mean(res^2, w_vec))
})

comp_tbl <- tibble(
  share = share_vars,
  rmse_w_2sls = as.numeric(rmse_w_iv),
  rmse_w_3sls = as.numeric(rmse_w_3sls)
)
print(comp_tbl)

# ------------------------------------------------------------
# (1) DIAGNÓSTICOS DE IV (weak/overid) a partir do 2SLS
# ------------------------------------------------------------
extract_iv_diags_explicit <- function(dat, dep_vec,
                                      price_vars  = paste0("ln_preco_com_reforma2", 1:6),
                                      exog        = c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4)),
                                      iv_regex    = "^iv_",
                                      iv_exclude  = "^(iv_sin|iv_cos)") {
  # IVs externos (exclui sazonais)
  ivs_ext <- grep(iv_regex, names(dat), value = TRUE)
  ivs_ext <- ivs_ext[!grepl(iv_exclude, ivs_ext)]
  inst    <- unique(c(exog, ivs_ext))
  
  # checagem rápida: precisamos de pelo menos 1 IV externo
  if (length(setdiff(inst, exog)) == 0) {
    message("⚠️ Nenhum IV externo restante (setdiff(inst, exog) vazio). Sem isso, não há diagnósticos.")
    return(tibble(share=character(), teste=character(),
                  statistic=double(), df1=double(), df2=double(), `p-value`=double()))
  }
  
  purrr::map_dfr(dep_vec, function(dep) {
    rhs1 <- paste(c(exog, price_vars), collapse = " + ")   # exog + endog (preços)
    rhs2 <- paste(inst, collapse = " + ")                  # exog + IVs externos (sem preços!)
    f    <- as.formula(sprintf("%s ~ %s | %s", dep, rhs1, rhs2))
    
    su <- summary(AER::ivreg(f, data = dat), diagnostics = TRUE)  # sem weights para habilitar diag
    dm <- su$diagnostics
    if (is.null(dm)) return(tibble())  # caso extremo
    
    tibble::as_tibble(dm, rownames = "teste") |>
      dplyr::mutate(share = dep, .before = 1)
  })
}

iv_diags <- extract_iv_diags_explicit(dat, dep_vec = paste0("w_despesahat", 1:6))
print(iv_diags)

exog <- c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4))
ivs_ext <- grep("^iv_", names(dat), value=TRUE)
ivs_ext <- ivs_ext[!grepl("^(iv_sin|iv_cos)", ivs_ext)]
setdiff(unique(c(exog, ivs_ext)), exog)

# ---------- 1) Helpers p/ montar fórmulas e checar o que SOBRA ----------
exog_base <- c("z","z2", paste0("iv_sin",1:4), paste0("iv_cos",1:4))

get_iv_externos <- function(dat) {
  ivs <- grep("^iv_", names(dat), value=TRUE)
  # tira sazonais do bloco de IVs externos
  ivs[!grepl("^(iv_sin|iv_cos)", ivs)]
}

fit_iv_diag <- function(dat, dep, price_vars, exog = exog_base) {
  ivs_ext <- get_iv_externos(dat)
  inst    <- unique(c(exog, ivs_ext))              # exog + IVs externos (sem preços!)
  rhs1    <- paste(c(exog, price_vars), collapse = " + ")
  rhs2    <- paste(inst,               collapse = " + ")
  f       <- as.formula(sprintf("%s ~ %s | %s", dep, rhs1, rhs2))
  
  # sem weights de propósito: summary(..., diagnostics=TRUE) não calcula diag com weights
  fit <- AER::ivreg(f, data = dat, x = TRUE, y = TRUE)
  
  # Mostra o que SOBROU após alias/NA:
  Xreg <- try(colnames(model.matrix(fit, "regressors")), silent=TRUE)
  Zinst<- try(colnames(model.matrix(fit, "instruments")), silent=TRUE)
  cat("\n=== ", dep, " ===\n", sep="")
  cat("Regressores usados (primeiro estágio):\n"); print(Xreg)
  cat("Instrumentos usados:\n");                print(Zinst)
  
  su <- summary(fit, diagnostics = TRUE)
  list(fit=fit, diag=su$diagnostics)
}

# ---------- 2) Rodar para todas as 6 equações ----------
share_vars  <- paste0("w_despesahat", 1:6)
price_vars  <- paste0("ln_preco_com_reforma2", 1:6)

diag_list <- lapply(share_vars, function(dep) fit_iv_diag(dat, dep, price_vars))
names(diag_list) <- share_vars

# ---------- 3) Tabela de diagnósticos (quando houver) ----------
iv_diags <- dplyr::bind_rows(lapply(names(diag_list), function(dep){
  dm <- diag_list[[dep]]$diag
  if (is.null(dm)) return(tibble::tibble())
  tibble::as_tibble(dm, rownames="teste") |>
    dplyr::mutate(share = dep, .before = 1)
}))

print(iv_diags)



# Sugestão de leitura rápida:
# dplyr::filter(iv_diags, grepl("Weak|first-stage", teste, ignore.case = TRUE))
# dplyr::filter(iv_diags, grepl("Sargan|Overid", teste, ignore.case = TRUE))

# ------------------------------------------------------------
# (2) SIMPLIFICAÇÃO DO RHS (remover z2 e/ou reduzir sazonais) + (3) 3SLS mínimo
# ------------------------------------------------------------

# Variante A: 3SLS com z apenas (sem z2) + sazonais reduzidos + IVs mínimos
fit_3sls_minA <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = spline_df,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE  # remove z2 do RHS (usa só z_c)
)
print(summary(fit_3sls_minA$fit))

pred_minA <- predict_3sls(fit_3sls_minA, dat)
rmse_w_3sls_minA <- sapply(share_vars, function(v) {
  sqrt(weighted.mean((obs_mat[[v]] - pred_minA[[v]])^2, w_vec))
})
print(tibble(share = share_vars, rmse_w_3sls_minA = as.numeric(rmse_w_3sls_minA)))

# Variante B: 3SLS com spline(z) + sazonais reduzidos + IVs mínimos
fit_3sls_minB <- fit_3sls_system(
  dat = dat, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = TRUE, spline_df = 4,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE
)
print(summary(fit_3sls_minB$fit))

pred_minB <- predict_3sls(fit_3sls_minB, dat)
rmse_w_3sls_minB <- sapply(share_vars, function(v) {
  sqrt(weighted.mean((obs_mat[[v]] - pred_minB[[v]])^2, w_vec))
})
print(tibble(share = share_vars, rmse_w_3sls_minB = as.numeric(rmse_w_3sls_minB)))

# ------------------------------------------------------------
# (4) TESTES CONJUNTOS PARA PREÇOS (2SLS com vcov robusto)
# ------------------------------------------------------------

joint_price_tests <- bind_rows(lapply(share_vars, function(v) {
  fit <- fits_iv[[v]]$fit
  cn  <- names(coef(fit))
  price_cn <- cn[grepl("^ln_preco_com_reforma2", cn)]
  if (length(price_cn) == 0) return(tibble(share = v, F = NA, df = NA, p = NA))
  hyp <- paste0(price_cn, " = 0")
  lh  <- car::linearHypothesis(fit, hyp, white.adjust = "hc1")  # usa vcovHC tipo HC1
  tibble(share = v, F = as.numeric(lh[2, "F"]), df = paste0(lh[2, "Df"], collapse = ","), p = as.numeric(lh[2, "Pr(>F)"]))
}))
print(joint_price_tests)

# ------------------------------------------------------------
# (5) VALIDAÇÃO OUT-OF-SAMPLE (hold-out simples 80/20)
# ------------------------------------------------------------

set.seed(123)
n   <- nrow(dat)
ix  <- sample.int(n, size = floor(0.8 * n))
trn <- dat[ix, ]
tst <- dat[-ix, ]

# Refit 2SLS no treino
fits_iv_trn <- map(share_vars, ~ fit_iv_weighted(
  dat = trn, dep = .x, price_vars = price_vars, weights_col = weights_col,
  iv_include_patterns = iv_include_patterns_full, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_full,
  use_spline = FALSE, spline_df = spline_df
))
names(fits_iv_trn) <- share_vars

# Refit 3SLS mínimo (Variante A) no treino
fit_3sls_trn <- fit_3sls_system(
  dat = trn, dep_vec = share_vars, price_vars = price_vars,
  iv_include_patterns = iv_include_patterns_min, iv_exclude_patterns = iv_exclude_patterns,
  z_var = z_var, z2_var = z2_var,
  seasonal_controls = seasonal_controls_min,
  use_spline = FALSE, spline_df = spline_df,
  drop_equation = "w_despesahat6",
  weights_col = weights_col,
  z2_present = FALSE
)

# Previsões no teste
pred_iv_tst <- sapply(share_vars, function(v) {
  as.numeric(predict(fits_iv_trn[[v]]$fit, newdata = tst))
})
colnames(pred_iv_tst) <- share_vars
pred_3sls_tst <- predict_3sls(fit_3sls_trn, tst)

w_tst <- tst[[weights_col]]
obs_tst <- as.data.frame(tst[, share_vars])

rmse_w_iv_tst <- sapply(share_vars, function(v) sqrt(weighted.mean((obs_tst[[v]] - pred_iv_tst[, v])^2, w_tst)))
rmse_w_3sls_tst <- sapply(share_vars, function(v) sqrt(weighted.mean((obs_tst[[v]] - pred_3sls_tst[[v]])^2, w_tst)))

oos_tbl <- tibble(
  share = share_vars,
  rmse_w_2sls_oos = as.numeric(rmse_w_iv_tst),
  rmse_w_3sls_oos = as.numeric(rmse_w_3sls_tst)
)
print(oos_tbl)

# ------------------------------------------------------------
# (6) CHECAGEM DAS SHARES DO 2SLS (soma=1?) E RECONCILIAÇÃO
# ------------------------------------------------------------

# Previsões 2SLS in-sample
pred_iv_in <- sapply(share_vars, function(v) as.numeric(predict(fits_iv[[v]]$fit)))
colnames(pred_iv_in) <- share_vars
sum_iv <- rowSums(pred_iv_in)

cat("\nResumo da soma das shares previstas (2SLS):\n")
print(summary(sum_iv))

# Histograma da soma
print(
  ggplot(data.frame(sum_iv = sum_iv), aes(x = sum_iv)) +
    geom_histogram(bins = 30, fill = "grey30") +
    geom_vline(xintercept = 1, linetype = 2) +
    labs(title = "Soma das shares previstas (2SLS)", x = "Σ shareŝ", y = "Contagem")
)

# Reconciliação simples: clip [0,1], renormaliza para soma=1
reconcile_01_norm <- function(mat) {
  M <- pmax(0, pmin(1, mat))
  s <- rowSums(M); s[s == 0] <- 1
  sweep(M, 1, s, "/")
}
pred_iv_in_recon <- reconcile_01_norm(pred_iv_in)

# RMSE ponderado após reconciliação
rmse_w_2sls_recon <- sapply(share_vars, function(v) {
  sqrt(weighted.mean((obs_mat[[v]] - pred_iv_in_recon[, v])^2, w_vec))
})
print(tibble(share = share_vars, rmse_w_2sls_in = as.numeric(rmse_w_iv),
             rmse_w_2sls_reconc = as.numeric(rmse_w_2sls_recon)))

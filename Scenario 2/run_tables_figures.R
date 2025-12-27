## ============================================================
## RUN – APPENDIX D (EXPORT TABLES + FIGURES / INSUMOS)
## Versão DEFINITIVA: inclui contrato completo para rodar reform_6.R
## (sem tentativa/erro), usando base_rob$compute_z do reform_1.R.
##
## Requisitos (arquivos do projeto no working directory):
## - reform_1.R ... reform_12.R
##
## Saídas:
## - CSV: write.csv2 em output_tables_D/
## - Figuras PNG: 300 dpi em output_figures_D/
## ============================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

## -----------------------------
## OUTPUT DIRS
## -----------------------------
out_dir_tables  <- "output_tables_D"
out_dir_figures <- "output_figures_D"
if (!dir.exists(out_dir_tables))  dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)

## -----------------------------
## LOG
## -----------------------------
log_file <- file.path(out_dir_tables, "run_appendixD_log.txt")
log_line <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse=""))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

## -----------------------------
## UTIL: write.csv2 seguro
## -----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

wcsv2 <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.vector(x) && !is.list(x)) x <- data.frame(value = as.character(x))
  if (is.list(x) && !is.data.frame(x)) {
    nm <- names(x) %||% seq_along(x)
    x  <- data.frame(name = nm, value = I(x), check.names = FALSE)
  }
  try(write.csv2(x, file = file, row.names = FALSE, na = ""), silent = TRUE)
  invisible(TRUE)
}

## -----------------------------
## UTIL: source só funções (evita top-level que depende de objetos)
## -----------------------------
safe_source_functions <- function(file, envir = .GlobalEnv) {
  exprs <- parse(file = file, keep.source = TRUE)
  
  is_fun_call <- function(x, name) is.call(x) && length(x) >= 1 && identical(as.character(x[[1]])[1], name)
  is_assignment <- function(x) is.call(x) && length(x) >= 3 &&
    (identical(as.character(x[[1]])[1], "<-") || identical(as.character(x[[1]])[1], "="))
  
  for (e in exprs) {
    ## library/require/source
    if (is.call(e) && length(e) >= 2) {
      head <- as.character(e[[1]])[1]
      if (head %in% c("library","require","source")) {
        try(eval(e, envir = envir), silent = TRUE)
        next
      }
    }
    ## name <- function(...)
    if (is_assignment(e)) {
      rhs <- e[[3]]
      if (is_fun_call(rhs, "function")) {
        eval(e, envir = envir)
        next
      }
    }
    ## assign("name", function(...))
    if (is_fun_call(e, "assign") && length(e) >= 3) {
      rhs <- e[[3]]
      if (is_fun_call(rhs, "function")) {
        eval(e, envir = envir)
        next
      }
    }
  }
  invisible(TRUE)
}

## -----------------------------
## UTIL: carregar apenas vars (atribuições) específicas de um script
## -----------------------------
load_contract_vars <- function(file, vars, envir = .GlobalEnv) {
  exprs <- parse(file = file, keep.source = TRUE)
  
  is_assign <- function(x) is.call(x) && length(x) >= 3 &&
    (identical(as.character(x[[1]])[1], "<-") || identical(as.character(x[[1]])[1], "="))
  
  lhs_name <- function(e) {
    lhs <- e[[2]]
    if (is.symbol(lhs)) return(as.character(lhs))
    NA_character_
  }
  
  for (e in exprs) {
    ## library/require/source (para RHS)
    if (is.call(e) && length(e) >= 2) {
      head <- as.character(e[[1]])[1]
      if (head %in% c("library","require","source")) {
        try(eval(e, envir = envir), silent = TRUE)
        next
      }
    }
    if (is_assign(e)) {
      nm <- lhs_name(e)
      if (!is.na(nm) && nm %in% vars) eval(e, envir = envir)
    }
  }
  invisible(TRUE)
}

## ============================================================
## STEP 1) RODAR REFORM_1 (base definitiva)
## ============================================================
log_line("[STEP 1] Sourcing reform_1.R (full run)")
source("reform_1.R", encoding = "UTF-8")

stopifnot(exists("df_iv"))
stopifnot(exists("priceNames"), exists("shareNames"))
stopifnot(exists("fit_q_3sls"))
stopifnot(exists("base_rob"), is.list(base_rob))
stopifnot(!is.null(base_rob$compute_z), is.function(base_rob$compute_z))

fit_q <- fit_q_3sls
stopifnot(!is.null(fit_q$fit))

## ============================================================
## STEP 2) CARREGAR FUNÇÕES DOS DEMAIS SCRIPTS (sem top-level)
## ============================================================
log_line("[STEP 2] Loading functions only from reform_2..reform_12")
for (f in sprintf("reform_%d.R", 2:12)) {
  if (file.exists(f)) {
    log_line("  - functions: ", f)
    safe_source_functions(f)
  } else {
    log_line("  - missing: ", f)
  }
}

## ============================================================
## D.1 (INSUMOS) – IV candidates
## ============================================================
log_line("[D.1] Export IV candidates inputs")
if (exists("IV_cand_df")) {
  wcsv2(as.data.frame(IV_cand_df), file.path(out_dir_tables, "D1_iv_candidates_IV_cand_df.csv"))
} else {
  log_line("[WARN] IV_cand_df não existe no ambiente.")
}
if (exists("iv_all_candidates")) wcsv2(iv_all_candidates, file.path(out_dir_tables, "D1_iv_candidates_iv_all_candidates.csv"))
if (exists("iv_base"))          wcsv2(iv_base,          file.path(out_dir_tables, "D1_iv_candidates_iv_base.csv"))
if (exists("iv_set_core"))      wcsv2(iv_set_core,      file.path(out_dir_tables, "D1_iv_candidates_iv_set_core.csv"))
if (exists("iv_set_mid"))       wcsv2(iv_set_mid,       file.path(out_dir_tables, "D1_iv_candidates_iv_set_mid.csv"))
if (exists("iv_set_full"))      wcsv2(iv_set_full,      file.path(out_dir_tables, "D1_iv_candidates_iv_set_full.csv"))

## ============================================================
## D.2 (INSUMOS) – RMSE/MAE por equação (se existirem)
## ============================================================
log_line("[D.2] Export RMSE/MAE inputs (from reform_1 objects)")
if (exists("pipeline_results")) {
  pr <- pipeline_results
  if (!is.null(pr$step5) && is.list(pr$step5)) {
    wcsv2(as.data.frame(pr$step5), file.path(out_dir_tables, "D2_rmse_mae_pipeline_step5.csv"))
  }
  ok_df <- try(as.data.frame(pr), silent = TRUE)
  if (!inherits(ok_df, "try-error")) {
    wcsv2(ok_df, file.path(out_dir_tables, "D2_rmse_mae_pipeline_results_as_df.csv"))
  }
} else {
  log_line("[WARN] pipeline_results não existe no ambiente (reform_1).")
}

## ============================================================
## Figure D.4 (INSUMO) – Coeficientes QUAIDS 3SLS + coefplot
## ============================================================
log_line("[Fig D.4] Export QUAIDS 3SLS coefficients (inputs + plot)")

cf <- try(coef(fit_q$fit), silent = TRUE)
V  <- try(vcov(fit_q$fit), silent = TRUE)
if (!inherits(cf, "try-error")) {
  df_coef <- data.frame(term = names(cf), estimate = as.numeric(cf), stringsAsFactors = FALSE)
  
  if (!inherits(V, "try-error") && is.matrix(V) && all(dim(V) > 0)) {
    se <- sqrt(pmax(0, diag(V)))
    names(se) <- colnames(V)
    df_coef$std_error <- se[df_coef$term]
    df_coef$conf_low  <- df_coef$estimate - qnorm(0.975) * df_coef$std_error
    df_coef$conf_high <- df_coef$estimate + qnorm(0.975) * df_coef$std_error
  }
  
  wcsv2(df_coef, file.path(out_dir_tables, "D4_quaids3sls_coefficients.csv"))
  
  ord <- if ("std_error" %in% names(df_coef)) {
    df_coef$tval <- df_coef$estimate / df_coef$std_error
    order(abs(df_coef$tval), decreasing = TRUE, na.last = NA)
  } else {
    order(abs(df_coef$estimate), decreasing = TRUE, na.last = NA)
  }
  topN <- min(60, length(ord))
  dfp  <- df_coef[ord[seq_len(topN)], , drop = FALSE]
  dfp  <- dfp[nrow(dfp):1, , drop = FALSE]
  
  png(file.path(out_dir_figures, "Figure_D4_quaids3sls_coefplot.png"),
      width = 2600, height = 1800, res = 300)
  par(mar = c(5, 26, 3, 2))
  y <- seq_len(nrow(dfp))
  plot(dfp$estimate, y, pch = 16, yaxt = "n",
       xlab = "Estimate", ylab = "", main = "QUAIDS 3SLS – Coefficients (Top terms)")
  axis(2, at = y, labels = dfp$term, las = 2, cex.axis = 0.6)
  abline(v = 0, lty = 2)
  if (all(c("conf_low","conf_high") %in% names(dfp))) {
    segments(dfp$conf_low, y, dfp$conf_high, y)
  }
  dev.off()
}

## ============================================================
## Figure D.3 (INSUMOS) – Benchmarking models (DEFINITIVO)
## - NÃO usa shareNames/priceNames aqui
## - Carrega contrato de nomes do reform_4.R
## - Garante z/z2 via base_rob$compute_z (reform_1)
## ============================================================
log_line("[Fig D.3] Deterministic contract + source reform_6.R")

## (A) Carregar funções necessárias do reform_4.R (caso ainda não existam)
need_funs <- c("fit_iv_weighted","fit_3sls_system","predict_3sls")
missing_f <- need_funs[!vapply(need_funs, exists, logical(1))]
if (length(missing_f) > 0) {
  log_line("  - loading missing functions from reform_4.R: ", paste(missing_f, collapse = ", "))
  safe_source_functions("reform_4.R")
}
missing_f2 <- need_funs[!vapply(need_funs, exists, logical(1))]
if (length(missing_f2) > 0) {
  stop(paste0("Faltam funções necessárias de reform_4.R: ", paste(missing_f2, collapse = ", ")))
}

## (B) Carregar variáveis de contrato do reform_4.R (nomes definitivos p/ benchmarking)
need_vars <- c(
  "share_vars","price_vars","weights_col",
  "seasonal_controls_full","seasonal_controls_min",
  "iv_include_patterns_full","iv_include_patterns_min","iv_exclude_patterns",
  "z_var","z2_var"
)
load_contract_vars("reform_4.R", need_vars)
missing_v <- need_vars[!vapply(need_vars, exists, logical(1))]
if (length(missing_v) > 0) {
  stop(paste0("Faltam variáveis de contrato de reform_4.R: ", paste(missing_v, collapse = ", ")))
}

## (C) dat definitivo para reform_6.R
dat <- if (exists("df_iv_ok")) df_iv_ok else df_iv
stopifnot(is.data.frame(dat))

## (D) garantir z/z2 via método oficial (base_rob$compute_z)
if (!all(c(z_var, z2_var) %in% names(dat))) {
  dat <- base_rob$compute_z(dat)
}
if (!all(c(z_var, z2_var) %in% names(dat))) {
  stop("Contrato reform_6: após base_rob$compute_z, dat ainda não contém z/z2.")
}

## (E) validar colunas exigidas por fit_iv_weighted()
need_cols <- unique(c(share_vars, price_vars, weights_col, z_var, z2_var))
miss_cols <- setdiff(need_cols, names(dat))
if (length(miss_cols) > 0) {
  stop(paste0("Contrato reform_6: dat não contém colunas exigidas: ", paste(miss_cols, collapse = ", ")))
}

## (F) rodar reform_6.R (full run)
source("reform_6.R", encoding = "UTF-8")

if (exists("comp_in_tbl")) wcsv2(as.data.frame(comp_in_tbl), file.path(out_dir_tables, "D3_benchmarking_comp_in_tbl.csv"))
if (exists("oos_tbl"))     wcsv2(as.data.frame(oos_tbl),     file.path(out_dir_tables, "D3_benchmarking_oos_tbl.csv"))

## Figura simples RMSE in-sample (se colunas existirem)
if (exists("comp_in_tbl")) {
  dfb <- as.data.frame(comp_in_tbl)
  if (all(c("share","rmse_w_2sls","rmse_w_3sls_minA","rmse_w_3sls_minB") %in% names(dfb))) {
    png(file.path(out_dir_figures, "Figure_D3_benchmarking_rmse_insample.png"),
        width = 2200, height = 1400, res = 300)
    par(mar = c(7, 6, 3, 2))
    x <- seq_len(nrow(dfb))
    ylim <- range(c(dfb$rmse_w_2sls, dfb$rmse_w_3sls_minA, dfb$rmse_w_3sls_minB), na.rm = TRUE)
    plot(x, dfb$rmse_w_2sls, type = "b", pch = 16, xaxt = "n",
         xlab = "", ylab = "Weighted RMSE", main = "Benchmarking (in-sample): RMSE by equation", ylim = ylim)
    axis(1, at = x, labels = dfb$share, las = 2, cex.axis = 0.8)
    lines(x, dfb$rmse_w_3sls_minA, type = "b", pch = 17)
    lines(x, dfb$rmse_w_3sls_minB, type = "b", pch = 15)
    legend("topright", legend = c("2SLS","3SLS minA","3SLS minB"),
           pch = c(16,17,15), lty = 1, bty = "n", cex = 0.9)
    dev.off()
  }
}

## ============================================================
## Tabelas D.3–D.10 (definitivo) – reform_12.R
## ============================================================
log_line("[Tables D.3–D.10] Sourcing reform_12.R (full run)")
source("reform_12.R", encoding = "UTF-8")

if (exists("all_tests")) {
  at <- all_tests
  if (!is.null(at$hansen)) wcsv2(as.data.frame(at$hansen), file.path(out_dir_tables, "D3_hansenJ_HC0_by_equation.csv"))
  if (!is.null(at$wu))     wcsv2(as.data.frame(at$wu),     file.path(out_dir_tables, "D4_wu_hausman_by_equation.csv"))
  if (!is.null(at$reset))  wcsv2(as.data.frame(at$reset),  file.path(out_dir_tables, "D5_ramsey_reset_by_equation.csv"))
  if (!is.null(at$wald))   wcsv2(as.data.frame(at$wald),   file.path(out_dir_tables, "D6_wald_block_prices_by_equation.csv"))
  if (!is.null(at$sargan)) wcsv2(as.data.frame(at$sargan), file.path(out_dir_tables, "D7_sargan_unweighted_by_equation.csv"))
  if (!is.null(at$sw))     wcsv2(as.data.frame(at$sw),     file.path(out_dir_tables, "D8_sanderson_windmeijer_by_equation.csv"))
  if (!is.null(at$kp))     wcsv2(as.data.frame(at$kp),     file.path(out_dir_tables, "D9_kleibergen_paap_rkWaldF_by_equation.csv"))
} else {
  log_line("[WARN] all_tests não existe após reform_12.R.")
}

if (exists("res_i_ii")) {
  rii <- res_i_ii
  if (!is.null(rii$AR_HC1))       wcsv2(as.data.frame(rii$AR_HC1),       file.path(out_dir_tables, "D10_anderson_rubin_HC1_by_equation.csv"))
  if (!is.null(rii$coef_2SLS))    wcsv2(as.data.frame(rii$coef_2SLS),    file.path(out_dir_tables, "D10_coef_prices_2SLS_HC1.csv"))
  if (!is.null(rii$coef_LIML))    wcsv2(as.data.frame(rii$coef_LIML),    file.path(out_dir_tables, "D10_coef_prices_LIML_HC1.csv"))
  if (!is.null(rii$coef_Fuller))  wcsv2(as.data.frame(rii$coef_Fuller),  file.path(out_dir_tables, "D10_coef_prices_Fuller_HC1.csv"))
}

## ============================================================
## Table D.12 (INSUMOS definitivos) – Delta CI elasticities: reform_10.R
## ============================================================
log_line("[Table D.12] Sourcing reform_10.R (full run)")
source("reform_10.R", encoding = "UTF-8")

if (exists("elas_marshall")) wcsv2(as.data.frame(elas_marshall), file.path(out_dir_tables, "D12_delta_elas_marshall.csv"))
if (exists("elas_hicks"))    wcsv2(as.data.frame(elas_hicks),    file.path(out_dir_tables, "D12_delta_elas_hicks.csv"))
if (exists("elas_income"))   wcsv2(as.data.frame(elas_income),   file.path(out_dir_tables, "D12_delta_elas_income_eta.csv"))

## ============================================================
## Table D.11 + Figure D.5 (definitivo) – Bootstrap elasticities: reform_11.R
## ============================================================
log_line("[Table D.11 / Fig D.5] Sourcing reform_11.R (full run)")
source("reform_11.R", encoding = "UTF-8")

if (exists("elas_marshall")) wcsv2(as.data.frame(elas_marshall), file.path(out_dir_tables, "D11_bootstrap_elas_marshall.csv"))
if (exists("elas_hicks"))    wcsv2(as.data.frame(elas_hicks),    file.path(out_dir_tables, "D11_bootstrap_elas_hicks.csv"))
if (exists("elas_income"))   wcsv2(as.data.frame(elas_income),   file.path(out_dir_tables, "D11_bootstrap_elas_income_eta.csv"))

## Heatmap simples (Marshall ponto e_mar) como “Figure D.5”
if (exists("elas_marshall")) {
  
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  library(tidyr)
  
  heat_B5 <- read.csv2('C:/Users/x16610962/Downloads/Reproduction/QUAIDS/Scripts/Scenario 2/output_tables_D/D11_bootstrap_elas_marshall.csv')   # marshallian: colunas i, j, est, se, lwr, upr
  ## ASSUMÇÃO: colunas i e j ou good_i e good_j estão presentes;
  ## se não, adapte conforme a estrutura do seu out_delta.
  names(heat_B5)[names(heat_B5) == 'share'] <- 'i'
  names(heat_B5)[names(heat_B5) == 'price'] <- 'j'
  if (!"good_i" %in% names(heat_B5) && "i" %in% names(heat_B5)) {
    heat_B5$good_i <- as.factor(heat_B5$i)
  }
  if (!"good_j" %in% names(heat_B5) && "j" %in% names(heat_B5)) {
    heat_B5$good_j <- as.factor(heat_B5$j)
  }
  
  
  heat_B5$sig <- ifelse(heat_B5$lo_mar<=0 & 0<=heat_B5$hi_mar, NA, heat_B5$e_mar)
  
  p_B5 <- ggplot(heat_B5, aes(x = good_j, y = good_i, fill = sig)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(    low = "#2166AC",
                             mid = "white",
                             high = "#B2182B",
                             midpoint = -0.44,
                             na.value = "grey80") +
    xlab("j") + ylab("i") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = round(sig, 2)), size = 3) +
    scale_y_discrete(
      labels = c("w exp 1", "w exp 2", "w exp 3", "w exp 4", "w exp 5", "w exp 6")
    ) +
    scale_x_discrete(
      labels = c("log price 1", "log price 2", "log price 3", "log price 4", "log price 5", "log price 6") 
    ) +
    labs(
      x = "Log prices",
      y = "Expenditure shares",
      fill = "Elasticity"
    )
  
  ggsave(file.path("C:/Users/x16610962/Downloads/Figure_D5_heatmap_elasticities.png"), p_B5, dpi = 300, width = 7, height = 5)
}

## ============================================================
## Figures D.1/D.2 (diagnóstico resíduos) – INSUMOS + PNGs (FIX DEFINITIVO)
## - residuals(systemfit) pode retornar lista com elementos não-atômicos
##   (ex.: data.frame/matrix/list aninhada). Este bloco normaliza de forma robusta.
## ============================================================
log_line("[Fig D.1/D.2] Diagnostic residual figures (insumos + png)")

to_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  
  ## já é atômico numérico
  if (is.numeric(x)) return(as.numeric(x))
  
  ## matriz/data.frame
  if (is.matrix(x) || is.data.frame(x)) {
    return(as.numeric(as.matrix(x)))
  }
  
  ## lista (inclusive aninhada): achata e converte
  if (is.list(x)) {
    u <- unlist(x, recursive = TRUE, use.names = FALSE)
    ## se unlist ainda não for numérico, tenta coerção segura
    if (is.numeric(u)) return(as.numeric(u))
    suppressWarnings(return(as.numeric(u)))
  }
  
  ## fallback: coerção direta (character/factor etc.)
  suppressWarnings(as.numeric(x))
}

res_sys <- try(residuals(fit_q$fit), silent = TRUE)

if (inherits(res_sys, "try-error") || is.null(res_sys)) {
  log_line("[WARN] residuals(fit_q$fit) falhou ou retornou NULL.")
} else {
  
  df_res <- NULL
  
  ## Caso padrão: lista por equação
  if (is.list(res_sys) && !is.data.frame(res_sys)) {
    nm_eq <- names(res_sys)
    if (is.null(nm_eq)) nm_eq <- paste0("eq", seq_along(res_sys))
    
    df_res <- do.call(rbind, lapply(seq_along(res_sys), function(i) {
      r <- to_numeric(res_sys[[i]])
      data.frame(
        eq    = nm_eq[i],
        idx   = seq_along(r),
        resid = r,
        stringsAsFactors = FALSE
      )
    }))
    
  } else {
    ## vetor/matriz/data.frame: normaliza tudo
    r <- to_numeric(res_sys)
    df_res <- data.frame(
      eq    = "system",
      idx   = seq_along(r),
      resid = r,
      stringsAsFactors = FALSE
    )
  }
  
  ## Se por algum motivo resid virou tudo NA, registre e pare (evita figuras vazias)
  if (!("resid" %in% names(df_res)) || all(is.na(df_res$resid))) {
    stop("D1/D2: resíduos normalizados resultaram em vetor vazio ou somente NA.")
  }
  
  ## Export insumo (long)
  wcsv2(df_res, file.path(out_dir_tables, "D1D2_insumo_systemfit_residuals_long.csv"))
  
  ## Histogram (global)
  png(file.path(out_dir_figures, "Figure_D1_residuals_hist_systemfit.png"),
      width = 2000, height = 1400, res = 300)
  par(mar = c(5, 5, 3, 2))
  hist(df_res$resid, breaks = 60,
       main = "In-sample residuals (systemfit) – histogram (all eq.)",
       xlab = "Residual")
  dev.off()
  
  ## QQ plot (global)
  png(file.path(out_dir_figures, "Figure_D2_residuals_qq_systemfit.png"),
      width = 2000, height = 1400, res = 300)
  par(mar = c(5, 5, 3, 2))
  qqnorm(df_res$resid, main = "In-sample residuals (systemfit) – QQ plot (all eq.)")
  qqline(df_res$resid, lty = 2)
  dev.off()
  
  ## Histogram por equação (um PNG por eq.)
  eqs <- unique(df_res$eq)
  for (e in eqs) {
    dfe <- df_res[df_res$eq == e, , drop = FALSE]
    fn  <- paste0("Figure_D1_residuals_hist_systemfit_", gsub("[^A-Za-z0-9_]+","_", e), ".png")
    png(file.path(out_dir_figures, fn),
        width = 2000, height = 1400, res = 300)
    par(mar = c(5, 5, 3, 2))
    hist(dfe$resid, breaks = 60,
         main = paste0("Residuals histogram – ", e),
         xlab = "Residual")
    dev.off()
  }
}


log_line("[DONE] Appendix D exports completed.")
cat("\n[OK] Saídas em:\n - ", out_dir_tables, "\n - ", out_dir_figures, "\n\n", sep = "")

## ============================================================
## RUN – APPENDIX C (REFORM PIPELINE)
## Requisitos: reform_1.R ... reform_10.R no working directory
## Saídas: CSV (tabelas) + PNG 300 dpi (figuras)
## ============================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

## -----------------------------
## OUTPUT DIRS
## -----------------------------
out_dir_tables  <- "output_tables_C"
out_dir_figures <- "output_figures_C"
if (!dir.exists(out_dir_tables))  dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)

## -----------------------------
## UTIL: logging simples
## -----------------------------
log_file <- file.path(out_dir_tables, "run_appendixC_log.txt")
log_line <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse=""))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

## -----------------------------
## UTIL: carregar só funções + library/require de um script
## (robusto contra top-level que usa objetos inexistentes: res, best_iv_aug etc.)
## -----------------------------
safe_source_functions <- function(file, envir = .GlobalEnv) {
  exprs <- parse(file = file, keep.source = TRUE)
  
  is_fun_call <- function(x, name) {
    is.call(x) && length(x) >= 1 && identical(as.character(x[[1]])[1], name)
  }
  
  is_assignment <- function(x) {
    is.call(x) && length(x) >= 3 && (identical(as.character(x[[1]])[1], "<-") || identical(as.character(x[[1]])[1], "="))
  }
  
  for (e in exprs) {
    
    ## library(...) / require(...)
    if (is.call(e) && length(e) >= 2) {
      head <- as.character(e[[1]])[1]
      if (identical(head, "library") || identical(head, "require")) {
        try(eval(e, envir = envir), silent = TRUE)
        next
      }
      if (identical(head, "source")) {
        ## manter compatibilidade com scripts que sourceiam helpers
        try(eval(e, envir = envir), silent = TRUE)
        next
      }
    }
    
    ## name <- function(...) { ... }
    if (is_assignment(e)) {
      rhs <- e[[3]]
      if (is_fun_call(rhs, "function")) {
        eval(e, envir = envir)
        next
      }
    }
    
    ## assign("name", function(...) {...})
    if (is_fun_call(e, "assign") && length(e) >= 3) {
      rhs <- e[[3]]
      if (is_fun_call(rhs, "function")) {
        eval(e, envir = envir)
        next
      }
    }
    
    ## qualquer outra coisa: IGNORA (não executa top-level)
  }
  
  invisible(TRUE)
}


## -----------------------------
## UTIL: chamada adaptativa (assinaturas variáveis)
## -----------------------------
call_with_formals <- function(fun, ...) {
  args <- list(...)
  fml  <- names(formals(fun))
  do.call(fun, args[names(args) %in% fml])
}

## -----------------------------
## UTIL: escolher cluster_var
## -----------------------------
pick_first_existing <- function(candidates, nms) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[1] else NA_character_
}

## ============================================================
## 1) RODAR REFORM_1.R (pipeline base: dados + fit principal)
## ============================================================
log_line("[STEP 1] Sourcing reform_1.R (full run)")
source("reform_1.R", encoding = "UTF-8")

## checagens mínimas (sem inferir)
stopifnot(exists("df_iv"))
stopifnot(exists("priceNames"), exists("shareNames"))
stopifnot(exists("fit_q_3sls"))

fit_q <- fit_q_3sls
stopifnot(!is.null(fit_q$fit), inherits(fit_q$fit, "systemfit"))

## IV set: usar o que existe no ambiente (não inventar)
iv_set <- NULL
if (exists("iv_set_mid"))  iv_set <- unique(as.character(iv_set_mid))
if (is.null(iv_set) && exists("iv_set_full")) iv_set <- unique(as.character(iv_set_full))
if (is.null(iv_set) && exists("iv_set_core")) iv_set <- unique(as.character(iv_set_core))
stopifnot(!is.null(iv_set), length(iv_set) > 0)

cluster_var <- pick_first_existing(
  c("psu", "cluster_id", "id_municipio", "id_setor_censitario", "id_domicilio"),
  names(df_iv)
)
if (is.na(cluster_var)) {
  df_iv$.__rowid__ <- seq_len(nrow(df_iv))
  cluster_var <- ".__rowid__"
  log_line("[INFO] Nenhum cluster encontrado; usando rowid como cluster (conservador).")
} else {
  log_line("[INFO] cluster_var = ", cluster_var)
}

## ============================================================
## 2) CARREGAR FUNÇÕES DOS DEMAIS SCRIPTS (sem executar top-level)
## ============================================================
log_line("[STEP 2] Loading functions only from reform_2..reform_10")

for (f in sprintf("reform_%d.R", 2:10)) {
  if (file.exists(f)) {
    log_line("  - functions: ", f)
    safe_source_functions(f)
  } else {
    log_line("  - missing: ", f)
  }
}

## ============================================================
## 3) FIGURE C.1 – Coefficients CI 95% (systemfit)
## ============================================================
log_line("[STEP 3] Figure C.1 (coefficients CI 95%)")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

cf <- coef(fit_q$fit)
V  <- tryCatch(vcov(fit_q$fit), error = function(e) NULL)

ln_price_prefixes <- paste0("ln_", priceNames)
df_coef <- tibble::tibble(term = names(cf), estimate = as.numeric(cf)) |>
  dplyr::filter(grepl(paste(ln_price_prefixes, collapse="|"), term))

## separa eq/var (systemfit pode usar eq.var ou eq_var)
df_coef <- df_coef |>
  tidyr::separate(term, into = c("eq","var"), sep = "[\\._]", extra = "merge", fill = "right")

if (!is.null(V)) {
  se <- sqrt(pmax(0, diag(V)))
  names(se) <- names(cf)
  
  ## tenta casar nomes do vcov com (eq.var) e (eq_var)
  df_coef <- df_coef |>
    dplyr::mutate(k1 = paste0(eq, ".", var),
                  k2 = paste0(eq, "_", var),
                  se = dplyr::coalesce(se[k1], se[k2])) |>
    dplyr::mutate(ci_low = estimate - 1.96*se,
                  ci_high= estimate + 1.96*se)
} else {
  df_coef <- df_coef |>
    dplyr::mutate(se = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  log_line("[WARN] vcov(systemfit) indisponível; CIs NA.")
}

p_C1 <- ggplot(df_coef, aes(x = var, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  coord_flip() +
  facet_wrap(~ eq, scales = "free_y") +
  labs(x = NULL, y = "Estimate")

ggsave(file.path(out_dir_figures, "Figure_C1_coefficients_CI95.png"),
       p_C1, width = 11, height = 8, dpi = 300)

write.csv2(df_coef, file.path(out_dir_tables, "Figure_C1_coefficients_CI95_data.csv"),
          row.names = FALSE)

## ============================================================
## 4) TABLES C.1–C.3 via weak_iv_suite_v3 (reform_2.R)
## ============================================================
log_line("[STEP 4] Tables C.1–C.3 via weak_iv_suite_v3")

if (!exists("weak_iv_suite_v3")) {
  write.csv2(data.frame(note="NOT_AVAILABLE: weak_iv_suite_v3 not found in reform_2.R"),
            file.path(out_dir_tables, "Table_C1_to_C3_NOT_AVAILABLE.csv"),
            row.names=FALSE)
  log_line("[WARN] weak_iv_suite_v3 não existe; C.1–C.3 não exportadas.")
} else {
  
  wiv <- tryCatch(
    weak_iv_suite_v3(
      fit_sys = fit_q,
      iv_set  = iv_set,
      df_all  = df_iv,
      price_of_interest = NULL,
      exogs_mode = "stone",
      include_shifters_in_Z = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(wiv, "error")) {
    ## erro típico do F_block (R %*% V mismatch). Não inventar: exporta diagnóstico e aborta esse bloco.
    write.csv2(data.frame(error=conditionMessage(wiv),
                         note="FAIL: weak_iv_suite_v3 erro. Verifique F_block em reform_2.R (mismatch coef/vcov)"),
              file.path(out_dir_tables, "Table_C1_to_C3_FAIL.csv"),
              row.names=FALSE)
    log_line("[FAIL] weak_iv_suite_v3: ", conditionMessage(wiv))
  } else {
    ## C.1
    tab_C1 <- wiv$sw_conditional |>
      dplyr::rename(Regressor = var,
                    df1 = df1,
                    df2 = df2,
                    F_statistic = F_SW,
                    p_value = p)
    write.csv2(tab_C1, file.path(out_dir_tables, "Table_C1_instrument_strength_Ftest.csv"),
              row.names = FALSE)
    
    ## C.2
    tab_C2 <- wiv$by_equation |>
      dplyr::select(eq, WuHausman_F, WuHausman_p, n_cc, k_endog_eff, k_excl_eff) |>
      dplyr::rename(Equation = eq,
                    WuHausman_F = WuHausman_F,
                    WuHausman_p = WuHausman_p,
                    n = n_cc,
                    k_endog = k_endog_eff,
                    k_excl  = k_excl_eff)
    write.csv2(tab_C2, file.path(out_dir_tables, "Table_C2_endogeneity_WuHausman.csv"),
              row.names = FALSE)
    
    ## C.3
    tab_C3 <- wiv$by_equation |>
      dplyr::select(eq, Sargan_stat, Sargan_p, n_cc, k_endog_eff, k_excl_eff) |>
      dplyr::rename(Equation = eq,
                    Sargan_stat = Sargan_stat,
                    Sargan_p = Sargan_p,
                    n = n_cc,
                    k_endog = k_endog_eff,
                    k_excl  = k_excl_eff)
    write.csv2(tab_C3, file.path(out_dir_tables, "Table_C3_overid_Sargan.csv"),
              row.names = FALSE)
    
    log_line("[OK] Exported C.1–C.3")
  }
}

## ============================================================
## 5) TABLE C.4 – Ramsey RESET (OLS benchmark)
## (Se não houver função explícita no reform_*: exporta com note UNVERIFIED)
## ============================================================
log_line("[STEP 5] Table C.4 RESET (OLS benchmark)")

if (!requireNamespace("lmtest", quietly = TRUE)) {
  write.csv2(data.frame(note="NOT_AVAILABLE: package lmtest not installed"),
            file.path(out_dir_tables, "Table_C4_RESET_NOT_AVAILABLE.csv"),
            row.names=FALSE)
  log_line("[WARN] lmtest ausente; C.4 RESET não exportada.")
} else {
  
  eqs <- setdiff(shareNames, fit_q$omit_share)
  endogs_all <- paste0("ln_", setdiff(fit_q$priceNames, fit_q$drop_price))
  rhs_base   <- intersect(c("z","z2"), names(df_iv))
  rhs_terms  <- c(endogs_all, rhs_base)
  
  reset_by_eq <- function(df, y_names, rhs_terms, powers = 2:3) {
    out <- lapply(y_names, function(y) {
      cols <- c(y, rhs_terms); cols <- cols[cols %in% names(df)]
      dat  <- df[complete.cases(df[, cols, drop=FALSE]), , drop=FALSE]
      f <- as.formula(paste0(y, " ~ ", paste(rhs_terms, collapse=" + ")))
      fit <- lm(f, data = dat)
      rt  <- tryCatch(lmtest::resettest(fit, power=powers, type="fitted"),
                      error=function(e) NULL)
      data.frame(
        Equation   = y,
        n          = nobs(fit),
        RESET_stat = if (!is.null(rt)) unname(rt$statistic) else NA_real_,
        df1        = if (!is.null(rt)) unname(rt$parameter[1]) else NA_real_,
        df2        = if (!is.null(rt)) unname(rt$parameter[2]) else NA_real_,
        p_value    = if (!is.null(rt)) unname(rt$p.value) else NA_real_,
        note       = "UNVERIFIED_BY_SCRIPTS: RESET not implemented explicitly in reform_1..10; computed via lmtest::resettest(type='fitted')",
        check.names = FALSE
      )
    })
    do.call(rbind, out)
  }
  
  tab_C4 <- reset_by_eq(df_iv, eqs, rhs_terms)
  write.csv2(tab_C4, file.path(out_dir_tables, "Table_C4_RESET_UNVERIFIED.csv"),
            row.names = FALSE)
}



## ============================================================
## 6) TABLE C.5 – AR/CLR testing
## (usa a função existente nos reform_*; se não existir, exporta placeholder)
## ============================================================
has_ivmodel      <- requireNamespace("ivmodel", quietly = TRUE)
has_AER          <- requireNamespace("AER", quietly = TRUE)
has_sandwich     <- requireNamespace("sandwich", quietly = TRUE)
has_lmtest       <- requireNamespace("lmtest", quietly = TRUE)
has_clubSandwich <- requireNamespace("clubSandwich", quietly = TRUE)

## Se o runner usa funções via ::, não precisa attach.
## Se usa sem ::, é seguro carregar quando disponível:
if (has_ivmodel)      suppressPackageStartupMessages(library(ivmodel))
if (has_AER)          suppressPackageStartupMessages(library(AER))
if (has_sandwich)     suppressPackageStartupMessages(library(sandwich))
if (has_lmtest)       suppressPackageStartupMessages(library(lmtest))
if (has_clubSandwich) suppressPackageStartupMessages(library(clubSandwich))


log_line("[STEP 6] Table C.5 AR/CLR")

## escolher POI existente (não inventar)
poi_candidates <- c("ln_preco_com_reforma13", "ln_preco_com_reforma12", "ln_preco_com_reforma11",
                    "ln_preco_com_reforma14", "ln_preco_com_reforma15", "ln_preco_com_reforma16")
poi_val <- poi_candidates[poi_candidates %in% names(df_iv)][1]

if (is.na(poi_val)) {
  write.csv2(data.frame(note="FAIL: no POI candidate found in df_iv"),
            file.path(out_dir_tables, "Table_C5_ARCLR_FAIL.csv"),
            row.names = FALSE)
  log_line("[FAIL] C.5: POI não encontrado no df_iv.")
} else {
  
  ## escolher runner disponível
  fun <- NULL
  fun_name <- NA_character_
  if (exists("run_arclr_all_v2")) { fun <- get("run_arclr_all_v2"); fun_name <- "run_arclr_all_v2" }
  else if (exists("run_arclr_all")) { fun <- get("run_arclr_all"); fun_name <- "run_arclr_all" }
  
  if (is.null(fun)) {
    write.csv2(data.frame(note="NOT_AVAILABLE: no AR/CLR runner found (run_arclr_all_v2 or run_arclr_all)"),
              file.path(out_dir_tables, "Table_C5_ARCLR_NOT_AVAILABLE.csv"),
              row.names = FALSE)
    log_line("[WARN] C.5: nenhuma função runner AR/CLR encontrada.")
  } else {
    
    fml <- names(formals(fun))
    
    ## mapear POI para o nome certo (sem suposição)
    poi_arg_name <- intersect(
      c("price_of_interest", "poi", "target", "endog", "endogenous", "price", "var", "price_var"),
      fml
    )[1]
    
    if (is.na(poi_arg_name)) {
      write.csv2(data.frame(
        note = paste0("FAIL: runner ", fun_name, " não possui argumento reconhecível para POI. Formals: ",
                      paste(fml, collapse=", "))
      ), file.path(out_dir_tables, "Table_C5_ARCLR_FAIL.csv"), row.names = FALSE)
      log_line("[FAIL] C.5: runner sem argumento reconhecível para POI: ", fun_name)
    } else {
      
      ## preparar argumentos candidatos (serão filtrados por formals)
      args <- list(
        df = df_iv,
        data = df_iv,
        shareNames = fit_q$shareNames,
        omit_share = fit_q$omit_share,
        drop_price = fit_q$drop_price,
        priceNames = fit_q$priceNames,
        iv_names   = iv_set,
        iv_set     = iv_set,
        halfwidth  = 200
      )
      args[[poi_arg_name]] <- poi_val
      
      ## filtrar pelos argumentos realmente aceitos
      args <- args[names(args) %in% fml]
      
      restore_arclr <- NULL
      
      if (exists("ar_clr_ci_grid_safe", mode = "function")) {
        
        ar_clr_ci_grid_safe_orig <- get("ar_clr_ci_grid_safe", mode = "function")
        
        ar_clr_ci_grid_safe_fixed <- function(Y, D, Z, X = NULL, alpha = 0.05,
                                              grid = NULL, halfwidth = 200) {
          # grid default: centrado em 0 (mesma lógica típica do projeto)
          if (is.null(grid)) {
            grid <- seq(-halfwidth, halfwidth, length.out = 801)
          }
          
          has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)
          has_ivpack  <- requireNamespace("ivpack",  quietly = TRUE)
          
          ## 1) ivmodel (preferencial)
          if (has_ivmodel) {
            ivm <- ivmodel::ivmodel(Y = as.numeric(Y), D = as.numeric(D), Z = Z, X = X)
            
            p_ar  <- vapply(grid, function(b0) {
              tryCatch(ivmodel::AR.test(ivm, beta0 = b0)$p.value, error = function(e) NA_real_)
            }, numeric(1))
            
            p_clr <- vapply(grid, function(b0) {
              tryCatch(ivmodel::CLR(ivm, beta0 = b0)$p.value, error = function(e) NA_real_)
            }, numeric(1))
            
            keep_ar  <- which(is.finite(p_ar)  & (p_ar  >= alpha))
            keep_clr <- which(is.finite(p_clr) & (p_clr >= alpha))
            
            ci_ar  <- if (length(keep_ar))  c(min(grid[keep_ar]),  max(grid[keep_ar]))  else c(NA_real_, NA_real_)
            ci_clr <- if (length(keep_clr)) c(min(grid[keep_clr]), max(grid[keep_clr])) else c(NA_real_, NA_real_)
            
            return(list(ci_ar = ci_ar, ci_clr = ci_clr))
          }
          
          ## 2) ivpack (se existir)
          if (has_ivpack) {
            get_p <- function(b0) {
              tryCatch(
                ivpack::AndersonRubinTest(
                  y = as.numeric(Y),
                  d = as.numeric(D),
                  z = Z,
                  x = if (is.null(X)) NULL else X,
                  beta = b0,
                  hetero = TRUE
                )$p.value,
                error = function(e) NA_real_
              )
            }
            p_ar <- vapply(grid, get_p, numeric(1))
            keep <- which(is.finite(p_ar) & (p_ar >= alpha))
            ci_ar <- if (length(keep)) c(min(grid[keep]), max(grid[keep])) else c(NA_real_, NA_real_)
            # sem CLR no ivpack: devolve AR também para CLR (conservador)
            return(list(ci_ar = ci_ar, ci_clr = ci_ar))
          }
          
          ## 3) fallback toy (ANOVA), NA-safe
          p_ar <- vapply(grid, function(b0) {
            yb <- as.numeric(Y) - b0 * as.numeric(D)
            df0 <- data.frame(yb = yb)
            
            if (!is.null(X)) df0 <- cbind(df0, as.data.frame(X))
            if (!is.null(Z)) df0 <- cbind(df0, as.data.frame(Z))
            
            xrhs <- if (is.null(X)) "1" else paste(colnames(X), collapse = " + ")
            zrhs <- paste(colnames(Z), collapse = " + ")
            
            f0 <- stats::as.formula(paste("yb ~", xrhs))
            f1 <- stats::as.formula(paste("yb ~", paste(c(xrhs, zrhs), collapse = " + ")))
            
            a <- tryCatch(anova(stats::lm(f0, df0), stats::lm(f1, df0)),
                          error = function(e) NULL)
            if (is.null(a)) return(NA_real_)
            
            Fv <- as.numeric(a$F[2])
            if (!is.finite(Fv)) return(NA_real_)
            
            df1 <- ncol(Z)
            df2 <- nrow(df0) - (if (is.null(X)) 1 else ncol(as.matrix(X))) - ncol(Z) - 1
            if (!is.finite(df2) || df2 <= 0) return(NA_real_)
            
            1 - stats::pf(Fv, df1 = df1, df2 = df2)
          }, numeric(1))
          
          keep <- which(is.finite(p_ar) & (p_ar >= alpha))
          ci <- if (length(keep)) c(min(grid[keep]), max(grid[keep])) else c(NA_real_, NA_real_)
          list(ci_ar = ci, ci_clr = ci)
        }
        
        assign("ar_clr_ci_grid_safe", ar_clr_ci_grid_safe_fixed, envir = .GlobalEnv)
        restore_arclr <- function() assign("ar_clr_ci_grid_safe", ar_clr_ci_grid_safe_orig, envir = .GlobalEnv)
        
        log_line("[INFO] C.5: applied temporary NA-safe override for ar_clr_ci_grid_safe() in run_tables_figures.R")
      }
      
      tab_C5 <- tryCatch(
        do.call(fun, args),
        error = function(e) e
      )
      
      if (!is.null(restore_arclr)) {
        restore_arclr()
        log_line("[INFO] C.5: restored original ar_clr_ci_grid_safe() after runner call")
      }
      
      if (inherits(tab_C5, "error")) {
        write.csv2(data.frame(error = conditionMessage(tab_C5),
                             note = paste0("SOURCE: ", fun_name, " | POI_ARG: ", poi_arg_name)),
                  file.path(out_dir_tables, "Table_C5_ARCLR_FAIL.csv"),
                  row.names = FALSE)
        log_line("[FAIL] C.5: ", conditionMessage(tab_C5))
      } else {
        ## acrescenta rastreabilidade
        tab_C5$note <- paste0("SOURCE: ", fun_name, " | POI_ARG: ", poi_arg_name, " | POI: ", poi_val)
        write.csv2(tab_C5, file.path(out_dir_tables, "Table_C5_ARCLR_by_equation.csv"),
                  row.names = FALSE)
        log_line("[OK] Exported C.5 via ", fun_name, " (", poi_arg_name, "=", poi_val, ")")
      }
    }
  }
}


## ============================================================
## 7) TABLE C.6 – 2SLS IV-system with VCOV clusters
## ============================================================
log_line("[STEP 7] Table C.6 IV system")

if (!(exists("fit_ivreg_system") && exists("ivsys_tidy"))) {
  write.csv2(data.frame(note="NOT_AVAILABLE: fit_ivreg_system/ivsys_tidy not found in reform_1..10"),
            file.path(out_dir_tables, "Table_C6_NOT_AVAILABLE.csv"),
            row.names=FALSE)
  log_line("[WARN] C.6 indisponível: funções não encontradas.")
} else {
  endogs_all <- paste0("ln_", setdiff(fit_q$priceNames, fit_q$drop_price))
  exogs_eq <- intersect(c("z","z2"), names(df_iv))
  
  ivsys <- tryCatch(
    fit_ivreg_system(df=df_iv, shareNames=fit_q$shareNames, omit_share=fit_q$omit_share,
                     endogs=endogs_all, exogs=exogs_eq, iv_names=iv_set,
                     cluster=cluster_var),
    error=function(e) e
  )
  
  if (inherits(ivsys, "error")) {
    write.csv2(data.frame(error=conditionMessage(ivsys),
                         note="FAIL: fit_ivreg_system erro"),
              file.path(out_dir_tables, "Table_C6_FAIL.csv"),
              row.names=FALSE)
    log_line("[FAIL] C.6: ", conditionMessage(ivsys))
  } else {
    tab_C6 <- ivsys_tidy(ivsys)
    write.csv2(tab_C6, file.path(out_dir_tables, "Table_C6_IV_system_cluster.csv"),
              row.names=FALSE)
    log_line("[OK] Exported C.6")
  }
}

## ============================================================
## 8) TABLE C.7 – KP rk F and robust Hansen-J
## ============================================================
kp_hansen_table_from_ivsys <- function(ivsys) {
  do.call(rbind, lapply(ivsys, function(m) {
    
    ## KP rk F (igual)
    ss  <- try(summary(m$fit, diagnostics = TRUE), silent = TRUE)
    kpF <- NA_real_
    kpP <- NA_real_
    if (!inherits(ss, "try-error") && !is.null(ss$diagnostics)) {
      rn  <- rownames(ss$diagnostics)
      iKP <- grep("Kleibergen|Weak instruments", rn, ignore.case = TRUE)
      if (length(iKP)) {
        kpF <- unname(ss$diagnostics[iKP[1], "statistic"])
        kpP <- unname(ss$diagnostics[iKP[1], "p-value"])
      }
    }
    
    ## Hansen J robusto via momentos (com QR-prune para evitar singularidade)
    u  <- resid(m$fit)
    Zm <- model.matrix(m$f_z, data = model.frame(m$fit))
    
    ## remover colunas constantes/NA e reduzir a posto completo
    ok_col <- apply(Zm, 2, function(v) all(is.finite(v)) && var(v) > 0)
    Zm <- Zm[, ok_col, drop = FALSE]
    
    qrZ <- qr(Zm)
    rZ  <- qrZ$rank
    if (rZ == 0) {
      J <- NA_real_; dfJ <- NA_integer_; pJ <- NA_real_
    } else {
      ZmR <- Zm[, qrZ$pivot[seq_len(rZ)], drop = FALSE]
      
      g    <- colMeans(ZmR * u)
      meat <- crossprod(ZmR * u) / nrow(ZmR)
      
      Winv <- try(solve(meat), silent = TRUE)
      J    <- if (inherits(Winv, "try-error")) NA_real_ else nrow(ZmR) * drop(t(g) %*% Winv %*% g)
      
      ## graus de liberdade devem usar o posto efetivo de Z, não ncol bruto
      dfJ  <- max(0L, ncol(ZmR) - length(coef(m$fit)))
      pJ   <- if (is.finite(J) && dfJ > 0) pchisq(J, df = dfJ, lower.tail = FALSE) else NA_real_
    }
    
    data.frame(eq = m$eq, KP_rkF = kpF, KP_p = kpP, HansenJ = J, dfJ = dfJ, HansenJ_p = pJ)
  }))
}

log_line("[STEP 8] Table C.7 KP/Hansen")

if (!(exists("fit_ivreg_system") && exists("kp_hansen_table_from_ivsys"))) {
  write.csv2(data.frame(note="NOT_AVAILABLE: fit_ivreg_system or kp_hansen_table_from_ivsys missing"),
            file.path(out_dir_tables, "Table_C7_NOT_AVAILABLE.csv"),
            row.names = FALSE)
  log_line("[WARN] C.7 indisponível: funções não encontradas.")
} else {
  
  ## mesmos insumos do fit_ivreg_system em reform_3.R
  endogs_all <- paste0("ln_", setdiff(fit_q$priceNames, fit_q$drop_price))
  exogs_eq   <- intersect(c("z","z2"), names(df_iv))
  
  ivsys <- tryCatch(
    fit_ivreg_system(
      df         = df_iv,
      shareNames = fit_q$shareNames,
      omit_share = fit_q$omit_share,
      endogs     = endogs_all,
      exogs      = exogs_eq,
      iv_names   = iv_set,
      cluster    = cluster_var
    ),
    error = function(e) e
  )
  
  if (inherits(ivsys, "error")) {
    write.csv2(data.frame(error=conditionMessage(ivsys), note="FAIL: fit_ivreg_system error"),
              file.path(out_dir_tables, "Table_C7_FAIL.csv"),
              row.names = FALSE)
    log_line("[FAIL] C.7: fit_ivreg_system: ", conditionMessage(ivsys))
  } else {
    
    tab_C7 <- tryCatch(
      kp_hansen_table_from_ivsys(ivsys),
      error = function(e) e
    )
    
    if (inherits(tab_C7, "error")) {
      write.csv2(data.frame(error=conditionMessage(tab_C7),
                           note="FAIL: kp_hansen_table_from_ivsys error"),
                file.path(out_dir_tables, "Table_C7_FAIL.csv"),
                row.names = FALSE)
      log_line("[FAIL] C.7: ", conditionMessage(tab_C7))
    } else {
      tab_C7$note <- "SOURCE: reform_3.R top-level block (re-implemented in run_* without executing top-level)"
      write.csv2(tab_C7,
                file.path(out_dir_tables, "Table_C7_KP_HansenJ.csv"),
                row.names = FALSE)
      log_line("[OK] Exported C.7")
    }
  }
}


## ============================================================
## 9) TABLE C.8 – leave-one-IV-out & split-sample
## ============================================================
log_line("[STEP 9] Table C.8 LOO / split-sample")

if (!(exists("loo_iv") && exists("split_sample_compare"))) {
  write.csv2(data.frame(note="NOT_AVAILABLE: loo_iv/split_sample_compare not found in reform_1..10"),
            file.path(out_dir_tables, "Table_C8_NOT_AVAILABLE.csv"),
            row.names=FALSE)
  log_line("[WARN] C.8 indisponível: funções não encontradas.")
} else {
  
  endogs_all <- paste0("ln_", setdiff(fit_q$priceNames, fit_q$drop_price))
  exogs_eq <- intersect(c("z","z2"), names(df_iv))
  
  target <- poi_candidates[poi_candidates %in% names(df_iv)][1]
  if (is.na(target)) {
    write.csv2(data.frame(note="FAIL: no target POI found for LOO"),
              file.path(out_dir_tables, "Table_C8_FAIL.csv"),
              row.names=FALSE)
    log_line("[FAIL] C.8: POI target não encontrado.")
  } else {
    
    tab_C8_loo <- tryCatch(
      loo_iv(df=df_iv, shareNames=fit_q$shareNames, omit_share=fit_q$omit_share,
             endogs=endogs_all, exogs=exogs_eq, iv_names=iv_set,
             cluster=cluster_var, target=target),
      error=function(e) e
    )
    if (inherits(tab_C8_loo, "error")) {
      write.csv2(data.frame(error=conditionMessage(tab_C8_loo), note="FAIL: loo_iv error"),
                file.path(out_dir_tables, "Table_C8_leave_one_IV_out_FAIL.csv"),
                row.names=FALSE)
      log_line("[FAIL] C.8 LOO: ", conditionMessage(tab_C8_loo))
    } else {
      write.csv2(tab_C8_loo, file.path(out_dir_tables, "Table_C8_leave_one_IV_out.csv"),
                row.names=FALSE)
    }
    
    ## --- garante compatibilidade com split_sample_compare() do reform_3.R ---
    price_of_interest <- target
    poi <- target
    
    tab_C8_split <- tryCatch(
      {
        fml <- names(formals(split_sample_compare))
        args <- list(
          df         = df_iv,
          shareNames = fit_q$shareNames,
          omit_share = fit_q$omit_share,
          endogs     = endogs_all,
          exogs      = intersect(c("z","z2"), names(df_iv)),
          iv_names   = iv_set,
          cluster    = cluster_var
        )
        ## se a função aceitar explicitamente o POI, passa; se não, o global acima cobre
        if ("price_of_interest" %in% fml) args$price_of_interest <- target
        if ("poi" %in% fml) args$poi <- target
        if ("target" %in% fml) args$target <- target
        
        args <- args[names(args) %in% fml]
        do.call(split_sample_compare, args)
      },
      error = function(e) e
    )
    
    if (inherits(tab_C8_split, "error")) {
      write.csv2(
        data.frame(error = conditionMessage(tab_C8_split), note = "FAIL: split_sample_compare error"),
        file.path(out_dir_tables, "Table_C8_split_sample_compare_FAIL.csv"),
        row.names = FALSE
      )
    } else {
      write.csv2(
        tab_C8_split,
        file.path(out_dir_tables, "Table_C8_split_sample_compare.csv"),
        row.names = FALSE
      )
    }
    
    
    log_line("[OK] Exported C.8 (as available)")
  }
}

## ============================================================
## 10) TABLE C.9 – Alternate normalization
## Sem rotina inequívoca em reform_1..10 => exporta placeholder com note
## ============================================================
log_line("[STEP 10] Table C.9 alternate normalization (from reform_3.R P6.1)")

## pré-condições (sem inferir)
stopifnot(exists("fit_ivreg_system"))
stopifnot(exists("df_iv"), exists("fit_q"), exists("iv_set"))
stopifnot(!is.null(fit_q$priceNames), !is.null(fit_q$shareNames))

## objetos esperados no bloco P6.1 do reform_3.R
df_aug      <- df_iv
priceNames  <- fit_q$priceNames
shareNames  <- fit_q$shareNames
drop_price  <- fit_q$drop_price
omit_share  <- fit_q$omit_share
exogs       <- intersect(c("z","z2"), names(df_aug))
best_iv_aug <- iv_set

## preço de interesse: usar o que já foi definido no pipeline (C.8), senão escolher candidato existente (sem inventar)
if (!exists("price_of_interest") || is.na(price_of_interest) || !(price_of_interest %in% names(df_aug))) {
  poi_candidates <- c("ln_preco_com_reforma13","ln_preco_com_reforma12","ln_preco_com_reforma11",
                      "ln_preco_com_reforma14","ln_preco_com_reforma15","ln_preco_com_reforma16")
  price_of_interest <- poi_candidates[poi_candidates %in% names(df_aug)][1]
}
if (is.na(price_of_interest)) {
  write.csv2(data.frame(note="FAIL: price_of_interest not found and no candidate exists in df_iv"),
            file.path(out_dir_tables, "Table_C9_alternate_normalization_FAIL.csv"),
            row.names = FALSE)
  log_line("[FAIL] C.9: price_of_interest não encontrado.")
} else {
  
  ## exatamente como no reform_3.R: dp_alt/os_alt pegam o 1º elemento do setdiff
  dp_candidates <- setdiff(priceNames, drop_price)
  os_candidates <- setdiff(shareNames, omit_share)
  
  if (!length(dp_candidates) || !length(os_candidates)) {
    write.csv2(data.frame(note="FAIL: cannot build dp_alt/os_alt from setdiff(priceNames,drop_price) or setdiff(shareNames,omit_share)"),
              file.path(out_dir_tables, "Table_C9_alternate_normalization_FAIL.csv"),
              row.names = FALSE)
    log_line("[FAIL] C.9: dp_alt/os_alt não disponíveis.")
  } else {
    
    dp_alt <- dp_candidates[1]
    os_alt <- os_candidates[1]
    
    ## replica alt_norm_results() do reform_3.R
    end_alt <- paste0("ln_", setdiff(priceNames, dp_alt))
    eqs_alt <- setdiff(shareNames, os_alt)
    
    ivsys_alt <- tryCatch(
      fit_ivreg_system(df_aug, shareNames, os_alt, end_alt, exogs, best_iv_aug, cluster_var),
      error = function(e) e
    )
    
    if (inherits(ivsys_alt, "error")) {
      write.csv2(data.frame(error=conditionMessage(ivsys_alt),
                           note="FAIL: fit_ivreg_system under alternate normalization"),
                file.path(out_dir_tables, "Table_C9_alternate_normalization_FAIL.csv"),
                row.names = FALSE)
      log_line("[FAIL] C.9: fit_ivreg_system: ", conditionMessage(ivsys_alt))
    } else {
      
      tab_C9 <- tibble::tibble(
        Equation = names(ivsys_alt),
        Price_of_interest_coefficient = purrr::map_dbl(ivsys_alt, function(m) {
          cc <- coef(m$fit)
          if (price_of_interest %in% names(cc)) unname(cc[price_of_interest]) else NA_real_
        }),
        dp_alt = dp_alt,
        os_alt = os_alt,
        price_of_interest = price_of_interest,
        note = "SOURCE: reform_3.R P6.1 alt_norm_results (swap drop_price / omit_share via first setdiff element)"
      )
      
      write.csv2(tab_C9,
                file.path(out_dir_tables, "Table_C9_alternate_normalization.csv"),
                row.names = FALSE)
      
      log_line("[OK] Exported C.9 (dp_alt=", dp_alt, ", os_alt=", os_alt, ", poi=", price_of_interest, ")")
    }
  }
}

## ============================================================
## 11) TABLE C.10 – Price suites (reform_4.R)
## ============================================================

log_line("[STEP 11] Table C.10 price suites")

if (!exists("run_iv_suite_all_prices")) {
  
  write.csv2(
    data.frame(note="NOT_AVAILABLE: run_iv_suite_all_prices not found (reform_4.R)"),
    file.path(out_dir_tables, "Table_C10_NOT_AVAILABLE.csv"),
    row.names = FALSE
  )
  log_line("[WARN] C.10 indisponível: run_iv_suite_all_prices não encontrada.")
  
} else {
  
  ## Chamada conforme assinatura real no reform_4.R:
  ## run_iv_suite_all_prices(df, priceNames, shareNames, omit_share, drop_price, base_iv_set, ...)
  
  result_all <- tryCatch(
    call_with_formals(
      run_iv_suite_all_prices,
      df         = df_iv,
      priceNames = fit_q$priceNames,
      shareNames = fit_q$shareNames,
      omit_share = fit_q$omit_share,
      drop_price = fit_q$drop_price,
      base_iv_set = iv_set,
      
      ## opcionais (só entram se estiverem no formals)
      exogs = c("z","z2"),
      alpha = 0.05,
      halfwidth = 50,
      shifters_for_iv = c("SH_SH_area_1","SH_SH_capital_1")
    ),
    error = function(e) e
  )
  
  if (inherits(result_all, "error")) {
    
    write.csv2(
      data.frame(error = conditionMessage(result_all),
                 note  = "FAIL: run_iv_suite_all_prices error"),
      file.path(out_dir_tables, "Table_C10_FAIL.csv"),
      row.names = FALSE
    )
    log_line("[FAIL] C.10: ", conditionMessage(result_all))
    
  } else {
    
    ## Verificação mínima: o objeto deve ter first_stage e ser tabular
    if (!("first_stage" %in% names(result_all))) {
      
      write.csv2(
        data.frame(note="FAIL: run_iv_suite_all_prices returned object without $first_stage"),
        file.path(out_dir_tables, "Table_C10_FAIL.csv"),
        row.names = FALSE
      )
      log_line("[FAIL] C.10: retorno sem first_stage.")
      
    } else {
      
      tab_C10 <- result_all$first_stage
      
      if (!inherits(tab_C10, c("data.frame","tbl_df","tbl"))) {
        write.csv2(
          data.frame(note="FAIL: result_all$first_stage is not a data.frame/tibble",
                     class = paste(class(tab_C10), collapse=",")),
          file.path(out_dir_tables, "Table_C10_FAIL.csv"),
          row.names = FALSE
        )
        log_line("[FAIL] C.10: first_stage não é tabular.")
      } else {
        
        ## Exporta a tabela C.10 (como no docx)
        tab_C10$note <- "SOURCE: reform_4.R::run_iv_suite_all_prices() -> $first_stage"
        write.csv2(
          tab_C10,
          file.path(out_dir_tables, "Table_C10_price_suites.csv"),
          row.names = FALSE
        )
        log_line("[OK] Exported C.10 (first_stage conditional by price).")
        
        ## (opcional) exporta insumos auxiliares do mesmo run (não são a C.10 do docx)
        if ("arclr_byprice_eq" %in% names(result_all) &&
            inherits(result_all$arclr_byprice_eq, c("data.frame","tbl_df","tbl"))) {
          tmp <- result_all$arclr_byprice_eq
          tmp$note <- "AUX: reform_4.R::run_iv_suite_all_prices() -> $arclr_byprice_eq"
          write.csv2(tmp, file.path(out_dir_tables, "C10_aux_arclr_byprice_eq.csv"), row.names = FALSE)
        }
        
        if ("summary_by_price" %in% names(result_all) &&
            inherits(result_all$summary_by_price, c("data.frame","tbl_df","tbl"))) {
          tmp <- result_all$summary_by_price
          tmp$note <- "AUX: reform_4.R::run_iv_suite_all_prices() -> $summary_by_price"
          write.csv2(tmp, file.path(out_dir_tables, "C10_aux_summary_by_price.csv"), row.names = FALSE)
        }
      }
    }
  }
}

## ============================================================
## [PATCH] TABLE C.5 – Sanderson–Windmeijer test (por preço endógeno)
## Fonte: reform_4.R::run_iv_suite_all_prices() -> result_all$first_stage
## Objetivo: exportar Table_C05_Sanderson_Windmeijer.csv em out_dir_tables
## ============================================================

## Requer:
## - out_dir_tables definido
## - result_all já calculado (saída de run_iv_suite_all_prices)
## - dplyr carregado (ou use dplyr:: prefixado, como abaixo)

tab_C5 <- NULL
if (exists("result_all") && is.list(result_all) && "first_stage" %in% names(result_all)) {
  tab_C5 <- result_all$first_stage
}

# Helper: retorna o 1º nome de coluna existente dentre candidatas
pick_col <- function(df, candidates) {
  nm <- candidates[candidates %in% names(df)]
  if (length(nm) == 0) return(NULL)
  nm[1]
}

if (!inherits(tab_C5, c("data.frame", "tbl_df", "tbl"))) {
  
  write.csv2(
    data.frame(
      note  = "FAIL: result_all$first_stage is not a data.frame/tibble",
      class = if (is.null(tab_C5)) "NULL" else paste(class(tab_C5), collapse = ",")
    ),
    file.path(out_dir_tables, "Table_C05_Sanderson_Windmeijer_FAIL.csv"),
    row.names = FALSE
  )
  log_line("[FAIL] C.5: result_all$first_stage não é tabular.")
  
} else {
  
  ## Nomes possíveis (robusto a variações)
  nm_F  <- pick_col(tab_C5, c("F_SW_cond", "F_SW", "F_SW_j"))
  nm_R2 <- pick_col(tab_C5, c("R2_partial", "R_Partial2", "R2p_j", "R_Partial_sq"))
  
  ## Checagens mínimas de colunas-base
  base_ok <- all(c("price", "n_cc", "df1", "df2") %in% names(tab_C5))
  
  if (!base_ok || is.null(nm_F) || is.null(nm_R2)) {
    
    write.csv2(
      data.frame(
        note = "FAIL: required columns not found to build Table C.5",
        has_price = "price" %in% names(tab_C5),
        has_n_cc  = "n_cc"  %in% names(tab_C5),
        has_df1   = "df1"   %in% names(tab_C5),
        has_df2   = "df2"   %in% names(tab_C5),
        picked_F  = if (is.null(nm_F))  "NONE" else nm_F,
        picked_R2 = if (is.null(nm_R2)) "NONE" else nm_R2,
        available_cols = paste(names(tab_C5), collapse = " | ")
      ),
      file.path(out_dir_tables, "Table_C05_Sanderson_Windmeijer_FAIL.csv"),
      row.names = FALSE
    )
    log_line("[FAIL] C.5: colunas necessárias não encontradas em first_stage.")
    
  } else {
    
    tab_C5_out <- dplyr::transmute(
      tab_C5,
      Variable      = .data$price,
      n_CC          = .data$n_cc,
      F_SW          = .data[[nm_F]],
      `R_Partial^2` = .data[[nm_R2]],
      df1           = .data$df1,
      df2           = .data$df2
    )
    
    tab_C5_out$note <- paste0(
      "SOURCE: reform_4.R::run_iv_suite_all_prices() -> $first_stage. ",
      "F col=", nm_F, "; R2 col=", nm_R2, "."
    )
    
    write.csv2(
      tab_C5_out,
      file.path(out_dir_tables, "Table_C05_Sanderson_Windmeijer.csv"),
      row.names = FALSE
    )
    log_line("[OK] Exported C.5 (Sanderson–Windmeijer test by endogenous price).")
  }
}


## ============================================================
## 12) TABLE C.11 – IV boosting
## Só exporta se o projeto tiver objetos/funções inequívocos já no ambiente.
## ============================================================

log_line("[STEP 12] Table C.11 IV boosting")

if (!(exists("df_iv") && exists("fit_q") && exists("iv_set"))) {
  
  write.csv2(
    data.frame(note="FAIL: missing df_iv / fit_q / iv_set in pipeline environment"),
    file.path(out_dir_tables, "Table_C11_IV_boosting_FAIL.csv"),
    row.names = FALSE
  )
  log_line("[FAIL] C.11: faltam df_iv/fit_q/iv_set.")
  
} else {
  
  ## Ambiente isolado, mas com acesso ao search path padrão
  env_boost <- new.env(parent = globalenv())
  
  ## garantir funções básicas
  env_boost$sd <- stats::sd
  
  ## 1) carregar dependências de FUNÇÕES (sem top-level) no mesmo env
  for (f in c("reform_2.R", "reform_3.R", "reform_4.R", "reform_5.R", "reform_6.R")) {
    if (file.exists(f)) {
      try(safe_source_functions(f, envir = env_boost), silent = TRUE)
    }
  }
  
  ## 2) injetar pré-requisitos declarados no reform_7.R / reform_8.R
  env_boost$df_iv_ok    <- df_iv
  env_boost$priceNames  <- fit_q$priceNames
  env_boost$shareNames  <- fit_q$shareNames
  env_boost$best_fit    <- fit_q
  env_boost$best_iv_set <- iv_set
  
  ## ------------------------------------------------------------
  ## FIX: hansenJ_poi_by_eq compatível com a chamada do reform_8.R
  ##      (aceita args nomeados e ignora os não usados via ...)
  ## ------------------------------------------------------------
  env_boost$hansenJ_poi_by_eq <- function(
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
      return(data.frame(eq=character(0), HansenJ=numeric(0), dfJ=integer(0), p=numeric(0), note=character(0)))
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
  
  ## 3) executar os scripts top-level do booster no env isolado
  c11_obj <- tryCatch({
    
    source("reform_7.R", local = env_boost)
    source("reform_8.R", local = env_boost)
    
    if (!exists("sw_tbl", envir = env_boost, inherits = FALSE)) {
      stop("reform_8.R did not create object 'sw_tbl'.")
    }
    
    tab <- get("sw_tbl", envir = env_boost, inherits = FALSE)
    
    if (!inherits(tab, c("data.frame","tbl_df","tbl"))) {
      stop("Object 'sw_tbl' exists but is not a data.frame/tibble.")
    }
    
    tab
    
  }, error = function(e) e)
  
  if (inherits(c11_obj, "error")) {
    
    write.csv2(
      data.frame(error = conditionMessage(c11_obj),
                 note  = "FAIL: executing reform_7.R + reform_8.R to produce sw_tbl"),
      file.path(out_dir_tables, "Table_C11_IV_boosting_FAIL.csv"),
      row.names = FALSE
    )
    log_line("[FAIL] C.11: ", conditionMessage(c11_obj))
    
  } else {
    
    tab_C11 <- as.data.frame(c11_obj)
    
    ## Padronizar nomes para o docx (sem inventar conteúdo)
    if ("price" %in% names(tab_C11)) {
      names(tab_C11)[names(tab_C11) == "price"] <- "Variable"
    }
    
    tab_C11$note <- "SOURCE: reform_8.R::sw_tbl (after loading function-deps into env_boost; then running reform_7.R+reform_8.R top-level)"
    
    write.csv2(
      tab_C11,
      file.path(out_dir_tables, "Table_C11_IV_boosting.csv"),
      row.names = FALSE
    )
    
    log_line("[OK] Exported C.11 (IV boosting) from reform_8.R sw_tbl.")
  }
}

## ============================================================
## 14) TABLE C.12 – Bootstrapped elasticities (FULL matrices)
##     + ADD: abs-cap + winsorization (configurable, no inference)
##     Policy:
##       - abs-cap is applied at the REPLICA level (before storing) by setting
##         entries with |x| > abs_cap_C12 to NA (so they are excluded by summ_vec()).
##       - winsorization is applied at the SUMMARIZATION level inside summ_vec()
##         (optional; if winsor_p_C12 is NULL, it is disabled).
##
##     Exports:
##       - FULL matrices (mean/q025/q975 in wide + long) + eta
##       - diagnostics include cap counts and effective n_ok per cell
## ============================================================

log_line("[STEP 14] Table C.12 bootstrap elasticities FULL matrices (Marshall/Hicks + eta) with abs-cap/winsor options")
abs_cap_C12 <- 10
winsor_p_C12 <- 0.02

## ----------------------------
## 0) Pre-conditions (no inference)
## ----------------------------
need_objs <- c("df_iv", "fit_q", "iv_set", "cluster_var", "out_dir_tables")
miss <- need_objs[!vapply(need_objs, exists, logical(1))]
if (length(miss)) {
  write.csv2(
    data.frame(note=paste0("FAIL: missing objects: ", paste(miss, collapse=", "))),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: missing objects: ", paste(miss, collapse=", "))
}

need_funs <- c(".refit_quaids_any", "elas_quaids_manual")
miss_f <- need_funs[!vapply(need_funs, exists, logical(1))]
if (length(miss_f)) {
  write.csv2(
    data.frame(note=paste0("FAIL: missing functions: ", paste(miss_f, collapse=", "))),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: missing functions: ", paste(miss_f, collapse=", "))
}

if (!(cluster_var %in% names(df_iv))) {
  write.csv2(
    data.frame(note=paste0("FAIL: cluster_var not in df_iv: ", cluster_var)),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: cluster_var not present in df_iv.")
}

best_fit <- fit_q

if (is.null(best_fit$priceNames) || !length(best_fit$priceNames)) {
  write.csv2(
    data.frame(note="FAIL: fit_q$priceNames is missing/empty"),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: fit_q$priceNames missing/empty.")
}

priceNames0 <- best_fit$priceNames
K0 <- length(priceNames0)

if (is.null(best_fit$shareNames) || !length(best_fit$shareNames)) {
  write.csv2(
    data.frame(note="FAIL: fit_q$shareNames is missing/empty"),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: fit_q$shareNames missing/empty.")
}
shareNames0 <- best_fit$shareNames

## ----------------------------
## 1) Determine x_name (strict, non-ambiguous)
## ----------------------------
x_name <- NULL
for (nm in c("x_name","xName","expName","expendName","xvar","x_var")) {
  if (!is.null(best_fit[[nm]]) && is.character(best_fit[[nm]]) && length(best_fit[[nm]]) == 1) {
    if (best_fit[[nm]] %in% names(df_iv)) {
      x_name <- best_fit[[nm]]
      break
    }
  }
}
if (is.null(x_name)) {
  explicit <- c("gasto_total_atualhat", "gasto_total_hat", "gasto_total_atual", "gasto_total")
  explicit <- explicit[explicit %in% names(df_iv)]
  if (length(explicit) == 1) x_name <- explicit[1]
}
if (is.null(x_name)) {
  candidates <- unique(c(
    grep("gasto_total|despesa_total|total_gasto|total_despesa", names(df_iv), value = TRUE),
    intersect(c("gasto_total_atualhat","gasto_total_hat","gasto_total_atual","gasto_total"), names(df_iv))
  ))
  write.csv2(
    data.frame(
      note="FAIL: could not determine x_name from fit_q or df_iv (ambiguous or missing).",
      candidate_vars=paste(candidates, collapse=" | ")
    ),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: x_name not found. See Table_C12_FAIL.csv for candidate_vars.")
}
log_line("[C12] Using x_name = ", x_name)

## ----------------------------
## 2) Bootstrap + robustness settings (pipeline vars if exist)
## ----------------------------
B_boot    <- if (exists("B_C12")) get("B_C12") else 2000
seed_boot <- if (exists("seed_C12")) get("seed_C12") else 123

normalize_eval  <- if (exists("normalize_eval"))  get("normalize_eval")  else "softmax"
normalize_deriv <- if (exists("normalize_deriv")) get("normalize_deriv") else "softmax"

x_eval <- if (exists("x_eval")) get("x_eval") else "median"
p_eval <- if (exists("p_eval")) get("p_eval") else "sample-mean-log"

shifters <- if (exists("shifters")) get("shifters") else NULL

## ---- NEW: abs-cap and winsorization controls ----
## abs_cap_C12: numeric; if finite, any |elas| > abs_cap_C12 is set to NA BEFORE storing
abs_cap_C12 <- if (exists("abs_cap_C12")) get("abs_cap_C12") else Inf
## winsor_p_C12: numeric in (0,0.5); if not NULL, winsorize inside summ_vec at [p, 1-p]
winsor_p_C12 <- if (exists("winsor_p_C12")) get("winsor_p_C12") else NULL

log_line("[C12] Robustness: abs_cap_C12=", abs_cap_C12,
         " winsor_p_C12=", if (is.null(winsor_p_C12)) "NULL" else winsor_p_C12)

## ----------------------------
## 3) Signature-aware caller for elas_quaids_manual (no guessing)
## ----------------------------
call_with_formals_local <- function(fun, ...) {
  dots <- list(...)
  fmls <- names(formals(fun))
  keep <- intersect(names(dots), fmls)
  do.call(fun, dots[keep])
}

elas_at_safe <- function(fit, x, p,
                         dlogp = 1e-3,
                         normalize_eval = "softmax",
                         normalize_deriv = "softmax") {
  call_with_formals_local(
    elas_quaids_manual,
    fit_q = fit,
    fit   = fit,
    x = x,
    p = p,
    eps = dlogp,
    dlogp = dlogp,
    target_dw = dlogp,
    normalize_eval  = normalize_eval,
    normalize_deriv = normalize_deriv
  )
}

## ----------------------------
## 4) p_star builder: MUST match priceNames (K) in order
## ----------------------------
get_priceNames_from_fit <- function(fit_obj, fallback_priceNames) {
  if (!is.null(fit_obj$priceNames) && length(fit_obj$priceNames)) return(fit_obj$priceNames)
  fallback_priceNames
}

get_p_star <- function(best_fit_for_eval, dfb, priceNamesK, p_eval = "sample-mean-log") {
  ln_cols <- paste0("ln_", priceNamesK)
  
  if (!is.null(best_fit_for_eval$data) && all(ln_cols %in% names(best_fit_for_eval$data))) {
    p <- exp(colMeans(best_fit_for_eval$data[, ln_cols, drop = FALSE], na.rm = TRUE))
    names(p) <- priceNamesK
    return(p[priceNamesK])
  }
  if (all(ln_cols %in% names(dfb))) {
    p <- exp(colMeans(dfb[, ln_cols, drop = FALSE], na.rm = TRUE))
    names(p) <- priceNamesK
    return(p[priceNamesK])
  }
  if (all(priceNamesK %in% names(dfb))) {
    if (identical(p_eval, "sample-mean-log")) {
      p <- exp(colMeans(log(dfb[, priceNamesK, drop = FALSE]), na.rm = TRUE))
    } else {
      p <- apply(dfb[, priceNamesK, drop = FALSE], 2, median, na.rm = TRUE)
    }
    names(p) <- priceNamesK
    return(p[priceNamesK])
  }
  
  stop("C12: cannot construct p_star for priceNamesK (missing ln_<price> in best_fit$data and dfb, and levels missing in dfb).")
}

## ----------------------------
## 5) Robustness helpers: abs-cap (replica level) + winsor (summarization)
## ----------------------------
cap_mat <- function(A, cap) {
  if (is.null(cap) || !is.finite(cap)) return(A)
  A2 <- A
  A2[is.finite(A2) & abs(A2) > cap] <- NA_real_
  A2
}

winsor_vec <- function(y, p) {
  if (is.null(p)) return(y)
  if (!is.numeric(p) || length(p) != 1 || !(p > 0 && p < 0.5)) return(y)
  lo <- as.numeric(stats::quantile(y, p,     na.rm = TRUE, type = 7))
  hi <- as.numeric(stats::quantile(y, 1 - p, na.rm = TRUE, type = 7))
  pmin(pmax(y, lo), hi)
}

## ----------------------------
## 6) Summarizers (for B x K x K arrays)  [WITH winsor option]
## ----------------------------
summ_vec <- function(x, winsor_p = NULL){
  ok <- is.finite(x)
  if (!any(ok)) return(c(mean=NA_real_, sd=NA_real_, q025=NA_real_, q975=NA_real_, n_ok=0))
  
  y <- x[ok]
  y <- winsor_vec(y, winsor_p)
  
  c(
    mean = mean(y),
    sd   = stats::sd(y),
    q025 = as.numeric(stats::quantile(y, 0.025, na.rm = TRUE, type = 7)),
    q975 = as.numeric(stats::quantile(y, 0.975, na.rm = TRUE, type = 7)),
    n_ok = length(y)
  )
}

summ_array_long <- function(A, what, winsor_p = NULL){
  dn2 <- dimnames(A)[[2]]; dn3 <- dimnames(A)[[3]]
  out <- vector("list", length(dn2) * length(dn3))
  k <- 1L
  for (i in seq_along(dn2)) for (j in seq_along(dn3)) {
    s <- summ_vec(A[, i, j], winsor_p = winsor_p)
    out[[k]] <- data.frame(
      i = dn2[i], j = dn3[j],
      estimate = s["mean"],
      sd = s["sd"],
      q025 = s["q025"],
      q975 = s["q975"],
      n_ok = s["n_ok"],
      type = what,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, out)
}

summ_array_wide <- function(A, winsor_p = NULL){
  dn2 <- dimnames(A)[[2]]; dn3 <- dimnames(A)[[3]]
  K <- length(dn2)
  M_mean <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_q025 <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_q975 <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_nok  <- matrix(0L,       K, K, dimnames=list(dn2, dn3))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    s <- summ_vec(A[,i,j], winsor_p = winsor_p)
    M_mean[i,j] <- s["mean"]
    M_q025[i,j] <- s["q025"]
    M_q975[i,j] <- s["q975"]
    M_nok[i,j]  <- as.integer(s["n_ok"])
  }
  list(mean = M_mean, q025 = M_q025, q975 = M_q975, n_ok = M_nok)
}

## ----------------------------
## 7) Bootstrap core: refit + elasticities; store FULL matrices on priceNames0 (K0)
##    + abs-cap applied BEFORE storing
## ----------------------------
boot_C12_full <- function(B, seed,
                          df_all, cluster_var,
                          best_fit_for_eval,
                          priceNames0, shareNames0,
                          x_name, instNames,
                          x_eval, p_eval,
                          normalize_eval, normalize_deriv,
                          shifters = NULL,
                          abs_cap = Inf) {
  
  set.seed(seed)
  
  cl  <- df_all[[cluster_var]]
  ucl <- unique(cl)
  ncl <- length(ucl)
  
  K0 <- length(priceNames0)
  
  acc_M   <- array(NA_real_, dim = c(B, K0, K0), dimnames = list(NULL, priceNames0, priceNames0))
  acc_H   <- array(NA_real_, dim = c(B, K0, K0), dimnames = list(NULL, priceNames0, priceNames0))
  acc_eta <- matrix(NA_real_, nrow = B, ncol = K0, dimnames = list(NULL, priceNames0))
  
  n_refit_ok  <- 0L
  n_elas_ok   <- 0L
  n_store_M   <- 0L
  n_store_H   <- 0L
  n_store_eta <- 0L
  
  ## NEW: cap counters
  cap_count_M <- 0L
  cap_count_H <- 0L
  
  err_log <- list()
  
  for (b in seq_len(B)) {
    
    cb  <- sample(ucl, size = ncl, replace = TRUE)
    idx <- which(cl %in% cb)
    dfb <- df_all[idx, , drop = FALSE]
    
    fitb <- try(
      .refit_quaids_any(
        dfb,
        priceNames = priceNames0,
        shareNames = shareNames0,
        x_name     = x_name,
        instNames  = instNames,
        shifters   = shifters
      ),
      silent = TRUE
    )
    if (inherits(fitb, "try-error")) {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="refit", error=as.character(fitb), stringsAsFactors=FALSE)
      next
    }
    n_refit_ok <- n_refit_ok + 1L
    
    priceNamesK <- get_priceNames_from_fit(fitb, priceNames0)
    
    ## hard requirement: same system dimension as the target table
    if (!setequal(priceNamesK, priceNames0) || length(priceNamesK) != length(priceNames0)) {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(
        draw=b, stage="priceNames_mismatch",
        error=paste0("priceNamesK != priceNames0. priceNamesK=", paste(priceNamesK, collapse=","),
                     " | priceNames0=", paste(priceNames0, collapse=",")),
        stringsAsFactors=FALSE
      )
      next
    }
    priceNamesK <- priceNames0  ## enforce order
    
    x_star <- if (identical(x_eval, "median")) median(dfb[[x_name]], na.rm = TRUE) else mean(dfb[[x_name]], na.rm = TRUE)
    
    p_star <- try(get_p_star(best_fit_for_eval, dfb, priceNamesK, p_eval = p_eval), silent = TRUE)
    if (inherits(p_star, "try-error")) {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="p_star", error=as.character(p_star), stringsAsFactors=FALSE)
      next
    }
    if (length(p_star) != length(priceNamesK)) {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(
        draw=b, stage="p_len",
        error=paste0("length(p_star)=", length(p_star), " != K(priceNamesK)=", length(priceNamesK)),
        stringsAsFactors=FALSE
      )
      next
    }
    
    Eb <- try(
      elas_at_safe(
        fitb,
        x = x_star,
        p = p_star,
        normalize_eval  = normalize_eval,
        normalize_deriv = normalize_deriv
      ),
      silent = TRUE
    )
    if (inherits(Eb, "try-error")) {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="elas", error=as.character(Eb), stringsAsFactors=FALSE)
      next
    }
    n_elas_ok <- n_elas_ok + 1L
    
    ## ---- Store FULL KxK (by priceNames0), with abs-cap applied BEFORE storing ----
    if (is.matrix(Eb$marshall)) {
      M_raw <- Eb$marshall
      if (!is.null(rownames(M_raw)) && !is.null(colnames(M_raw)) &&
          all(priceNames0 %in% rownames(M_raw)) && all(priceNames0 %in% colnames(M_raw))) {
        M <- M_raw[priceNames0, priceNames0, drop = FALSE]
      } else if (all(dim(M_raw) == c(K0, K0))) {
        M <- M_raw
        dimnames(M) <- list(priceNames0, priceNames0)
      } else {
        M <- NULL
        if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="store_M_dim", error=paste0("dim(marshall)=", paste(dim(M_raw), collapse="x")), stringsAsFactors=FALSE)
      }
      
      if (!is.null(M)) {
        if (is.finite(abs_cap)) {
          cap_count_M <- cap_count_M + sum(is.finite(M) & abs(M) > abs_cap)
          M <- cap_mat(M, abs_cap)
        }
        acc_M[b,,] <- M
        n_store_M <- n_store_M + 1L
      }
    } else {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="marshall_not_matrix", error=class(Eb$marshall)[1], stringsAsFactors=FALSE)
    }
    
    if (is.matrix(Eb$hicks)) {
      H_raw <- Eb$hicks
      if (!is.null(rownames(H_raw)) && !is.null(colnames(H_raw)) &&
          all(priceNames0 %in% rownames(H_raw)) && all(priceNames0 %in% colnames(H_raw))) {
        H <- H_raw[priceNames0, priceNames0, drop = FALSE]
      } else if (all(dim(H_raw) == c(K0, K0))) {
        H <- H_raw
        dimnames(H) <- list(priceNames0, priceNames0)
      } else {
        H <- NULL
        if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="store_H_dim", error=paste0("dim(hicks)=", paste(dim(H_raw), collapse="x")), stringsAsFactors=FALSE)
      }
      
      if (!is.null(H)) {
        if (is.finite(abs_cap)) {
          cap_count_H <- cap_count_H + sum(is.finite(H) & abs(H) > abs_cap)
          H <- cap_mat(H, abs_cap)
        }
        acc_H[b,,] <- H
        n_store_H <- n_store_H + 1L
      }
    } else {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="hicks_not_matrix", error=class(Eb$hicks)[1], stringsAsFactors=FALSE)
    }
    
    ## Expenditure elasticities (eta): store length K0 if available
    if (!is.null(Eb$expenditure)) {
      et <- Eb$expenditure
      if (!is.null(names(et)) && all(priceNames0 %in% names(et))) {
        acc_eta[b,] <- et[priceNames0]
        n_store_eta <- n_store_eta + 1L
      } else if (length(et) == K0) {
        acc_eta[b,] <- as.numeric(et)
        n_store_eta <- n_store_eta + 1L
      } else {
        if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="eta_len", error=paste0("length(expenditure)=", length(et)), stringsAsFactors=FALSE)
      }
    } else {
      if (length(err_log) < 100) err_log[[length(err_log)+1]] <- data.frame(draw=b, stage="eta_missing", error="Eb$expenditure is NULL", stringsAsFactors=FALSE)
    }
  }
  
  err_df <- if (length(err_log)) do.call(rbind, err_log) else data.frame()
  
  diag <- data.frame(
    B = B,
    refit_ok  = n_refit_ok,
    elas_ok   = n_elas_ok,
    stored_M  = n_store_M,
    stored_H  = n_store_H,
    stored_eta = n_store_eta,
    abs_cap = abs_cap,
    cap_count_M = cap_count_M,
    cap_count_H = cap_count_H,
    stringsAsFactors = FALSE
  )
  
  list(M = acc_M, H = acc_H, eta = acc_eta, diag = diag, err_sample = err_df)
}

## ----------------------------
## 8) Run bootstrap
## ----------------------------
boot_full <- tryCatch(
  boot_C12_full(
    B = B_boot,
    seed = seed_boot,
    df_all = df_iv,
    cluster_var = cluster_var,
    best_fit_for_eval = fit_q,
    priceNames0 = priceNames0,
    shareNames0 = shareNames0,
    x_name = x_name,
    instNames = iv_set,
    x_eval = x_eval,
    p_eval = p_eval,
    normalize_eval = normalize_eval,
    normalize_deriv = normalize_deriv,
    shifters = shifters,
    abs_cap = abs_cap_C12
  ),
  error = function(e) e
)

if (inherits(boot_full, "error")) {
  write.csv2(
    data.frame(error = conditionMessage(boot_full), note = "FAIL: boot_C12_full()"),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  log_line("[FAIL] C.12: ", conditionMessage(boot_full))
  stop("C12 failed.")
}

write.csv2(boot_full$diag, file.path(out_dir_tables, "Table_C12_boot_diagnostics.csv"), row.names = FALSE)

has_err_sample <- !is.null(boot_full$err_sample) &&
  inherits(boot_full$err_sample, "data.frame") &&
  nrow(boot_full$err_sample) > 0

if (has_err_sample) {
  write.csv2(
    boot_full$err_sample,
    file.path(out_dir_tables, "Table_C12_boot_error_sample.csv"),
    row.names = FALSE
  )
}

log_line("[C12] diag: refit_ok=", boot_full$diag$refit_ok,
         " elas_ok=", boot_full$diag$elas_ok,
         " stored_M=", boot_full$diag$stored_M,
         " stored_H=", boot_full$diag$stored_H,
         " stored_eta=", boot_full$diag$stored_eta,
         " cap_count_M=", boot_full$diag$cap_count_M,
         " cap_count_H=", boot_full$diag$cap_count_H,
         " (B=", boot_full$diag$B, ")")

if (boot_full$diag$stored_M == 0 && boot_full$diag$stored_H == 0) {
  write.csv2(
    data.frame(
      note="FAIL: bootstrap completed but stored zero elasticity matrices (see error sample).",
      refit_ok=boot_full$diag$refit_ok,
      elas_ok=boot_full$diag$elas_ok,
      stored_M=boot_full$diag$stored_M,
      stored_H=boot_full$diag$stored_H
    ),
    file.path(out_dir_tables, "Table_C12_FAIL.csv"),
    row.names = FALSE
  )
  stop("C12: no stored elasticities; see Table_C12_boot_error_sample.csv")
}

## ----------------------------
## 9) Export FULL matrices (LONG + WIDE) + eta
##    NOTE: winsorization is applied ONLY here (summarization stage) via winsor_p_C12
## ----------------------------
tab_M_long <- summ_array_long(boot_full$M, "Marshall", winsor_p = winsor_p_C12)
tab_H_long <- summ_array_long(boot_full$H, "Hicks",    winsor_p = winsor_p_C12)

write.csv2(tab_M_long, file.path(out_dir_tables, "Table_C12_boot_marshall_FULL_long.csv"), row.names = FALSE)
write.csv2(tab_H_long, file.path(out_dir_tables, "Table_C12_boot_hicks_FULL_long.csv"),    row.names = FALSE)

M_w <- summ_array_wide(boot_full$M, winsor_p = winsor_p_C12)
H_w <- summ_array_wide(boot_full$H, winsor_p = winsor_p_C12)

df_M_mean <- data.frame(Price = rownames(M_w$mean), M_w$mean, check.names = FALSE)
df_M_q025 <- data.frame(Price = rownames(M_w$q025), M_w$q025, check.names = FALSE)
df_M_q975 <- data.frame(Price = rownames(M_w$q975), M_w$q975, check.names = FALSE)
df_M_nok  <- data.frame(Price = rownames(M_w$n_ok), M_w$n_ok,  check.names = FALSE)

df_H_mean <- data.frame(Price = rownames(H_w$mean), H_w$mean, check.names = FALSE)
df_H_q025 <- data.frame(Price = rownames(H_w$q025), H_w$q025, check.names = FALSE)
df_H_q975 <- data.frame(Price = rownames(H_w$q975), H_w$q975, check.names = FALSE)
df_H_nok  <- data.frame(Price = rownames(H_w$n_ok), H_w$n_ok,  check.names = FALSE)

write.csv2(df_M_mean, file.path(out_dir_tables, "Table_C12_boot_marshall_FULL_mean_wide.csv"), row.names = FALSE)
write.csv2(df_M_q025, file.path(out_dir_tables, "Table_C12_boot_marshall_FULL_q025_wide.csv"), row.names = FALSE)
write.csv2(df_M_q975, file.path(out_dir_tables, "Table_C12_boot_marshall_FULL_q975_wide.csv"), row.names = FALSE)
write.csv2(df_M_nok,  file.path(out_dir_tables, "Table_C12_boot_marshall_FULL_nok_wide.csv"),  row.names = FALSE)

write.csv2(df_H_mean, file.path(out_dir_tables, "Table_C12_boot_hicks_FULL_mean_wide.csv"), row.names = FALSE)
write.csv2(df_H_q025, file.path(out_dir_tables, "Table_C12_boot_hicks_FULL_q025_wide.csv"), row.names = FALSE)
write.csv2(df_H_q975, file.path(out_dir_tables, "Table_C12_boot_hicks_FULL_q975_wide.csv"), row.names = FALSE)
write.csv2(df_H_nok,  file.path(out_dir_tables, "Table_C12_boot_hicks_FULL_nok_wide.csv"),  row.names = FALSE)

## eta summary (K0) -- optionally winsorize eta too (same winsor_p_C12)
eta_summ <- function(x) summ_vec(x, winsor_p = winsor_p_C12)
eta_tab <- t(apply(boot_full$eta, 2, eta_summ))
eta_df <- data.frame(
  Variable = rownames(eta_tab),
  estimate = eta_tab[, "mean"],
  sd       = eta_tab[, "sd"],
  q025     = eta_tab[, "q025"],
  q975     = eta_tab[, "q975"],
  n_ok     = eta_tab[, "n_ok"],
  stringsAsFactors = FALSE
)
write.csv2(eta_df, file.path(out_dir_tables, "Table_C12_boot_eta_FULL.csv"), row.names = FALSE)

## record robustness params (so tables are self-describing)
robust_note <- data.frame(
  abs_cap_C12 = abs_cap_C12,
  winsor_p_C12 = if (is.null(winsor_p_C12)) NA_real_ else winsor_p_C12,
  note = "abs_cap applied pre-store; winsor applied in summarization only",
  stringsAsFactors = FALSE
)
write.csv2(robust_note, file.path(out_dir_tables, "Table_C12_boot_robustness_settings.csv"), row.names = FALSE)

log_line("[OK] Exported C.12 FULL matrices (Marshall/Hicks KxK) + eta (K) with abs-cap/winsor options.")




## --- Figura C.2 – Heatmap de elasticidades (ASSUNÇÃO) ---
## ASSUNÇÃO: a figura usa as elasticidades marshallianas “pontuais”
##           (paramétricas ou bootstrap). Aqui uso as paramétricas
##           (out_delta$marshallian) por serem diretamente associadas
##           ao “modelo final”.

if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)

heat_B5 <- tab_M_long   # marshallian: colunas i, j, est, se, lwr, upr
## ASSUMÇÃO: colunas i e j ou good_i e good_j estão presentes;
## se não, adapte conforme a estrutura do seu out_delta.

if (!"good_i" %in% names(heat_B5) && "i" %in% names(heat_B5)) {
  heat_B5$good_i <- as.factor(heat_B5$i)
}
if (!"good_j" %in% names(heat_B5) && "j" %in% names(heat_B5)) {
  heat_B5$good_j <- as.factor(heat_B5$j)
}

heat_B5$sig <- ifelse(heat_B5$q025<=0 & 0<=heat_B5$q975, NA, heat_B5$estimate)

p_B5 <- ggplot(heat_B5, aes(x = good_j, y = good_i, fill = sig)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(    low = "#2166AC",
                           mid = "white",
                           high = "#B2182B",
                           midpoint = 0.095,
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

ggsave(file.path(out_dir_figures, "Figure_B2_heatmap_elasticities.png"), p_B5, dpi = 300, width = 7, height = 5)

## ============================================================
## Table C.12 — WILD BOOTSTRAP (cluster) — FULL PATCHED VERSION
## Adds:
##   - abs cap (replica-level, pre-store): |elas| > abs_cap_C12 => NA
##   - winsorization (summarization-level): winsor at [p, 1-p] inside summarizers
##
## Requires (from your pipeline):
##   df_iv, fit_q, iv_set, cluster_var, out_dir_tables
##   functions: .refit_quaids_any(), auto_eq_map2()
##   and either elas_at_safe() OR elas_quaids_manual() available
## ============================================================

need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("stats","utils"))

## ---------------------------
## Helpers (ALL embedded)
## ---------------------------

.wild_weights_cluster <- function(cluster, type = c("rademacher","mammen")) {
  type <- match.arg(type)
  g <- as.factor(cluster)
  G <- levels(g)
  
  if (type == "rademacher") {
    vG <- sample(c(-1, 1), length(G), replace = TRUE)
  } else {
    ## Mammen 2-point
    s5 <- sqrt(5)
    a  <- (1 - s5) / 2
    b  <- (1 + s5) / 2
    pA <- (s5 + 1) / (2 * s5)
    vG <- ifelse(runif(length(G)) < pA, a, b)
  }
  vG[match(as.character(g), G)]
}

.enforce_simplex <- function(W, eps = 1e-10) {
  W2 <- as.matrix(W)
  W2[!is.finite(W2)] <- NA_real_
  W2[is.na(W2)] <- 0
  W2 <- pmax(W2, eps)
  
  rs <- rowSums(W2)
  rs[rs == 0] <- 1
  W2 <- W2 / rs
  
  W2 <- pmin(pmax(W2, eps), 1 - eps)
  rs <- rowSums(W2)
  rs[rs == 0] <- 1
  W2 / rs
}

.resolve_x_name <- function(fit_q, df_iv){
  x_name <- NULL
  for (nm in c("x_name","xName","expName","expendName","xvar","x_var")) {
    if (!is.null(fit_q[[nm]]) && is.character(fit_q[[nm]]) && length(fit_q[[nm]]) == 1) {
      if (fit_q[[nm]] %in% names(df_iv)) { x_name <- fit_q[[nm]]; break }
    }
  }
  if (is.null(x_name)) {
    explicit <- c("gasto_total_atualhat","gasto_total_hat","gasto_total_atual","gasto_total","x")
    explicit <- explicit[explicit %in% names(df_iv)]
    if (length(explicit)) x_name <- explicit[1]
  }
  if (is.null(x_name)) stop("wild-C12: não consegui resolver x_name.")
  x_name
}

.get_p_star_fallback <- function(df, priceNames){
  ln_prices <- paste0("ln_", priceNames)
  if (all(ln_prices %in% names(df))) {
    return(exp(colMeans(df[, ln_prices, drop=FALSE], na.rm = TRUE)))
  }
  if (all(priceNames %in% names(df))) {
    return(apply(df[, priceNames, drop=FALSE], 2, median, na.rm = TRUE))
  }
  stop("wild-C12: não encontrei preços/ln_preços para construir p_star.")
}

.get_fitted_shares <- function(fit_q, df, shareNames){
  fit_sys <- fit_q$fit
  if (is.null(fit_sys) || !inherits(fit_sys, "systemfit")) stop("wild-C12: fit_q$fit não é systemfit.")
  if (!exists("auto_eq_map2", mode = "function")) stop("wild-C12: falta auto_eq_map2() no ambiente.")
  
  omit_share <- if (!is.null(fit_q$omit_share)) as.character(fit_q$omit_share) else character(0)
  
  map <- auto_eq_map2(fit_sys, shareNames = shareNames, omit_share = omit_share)
  eq_safe <- map$eqs_safe
  eq_orig <- map$eqs_orig
  
  fh <- try(fitted(fit_sys), silent = TRUE)
  if (inherits(fh, "try-error")) stop("wild-C12: fitted(systemfit) falhou.")
  
  if (is.list(fh)) {
    if (is.null(names(fh))) stop("wild-C12: fitted retornou lista sem names.")
    if (!all(eq_safe %in% names(fh))) stop("wild-C12: fitted(list) não contém todas as eqs_safe.")
    M <- do.call(cbind, fh[eq_safe])
  } else {
    if (is.data.frame(fh)) fh <- as.matrix(fh)
    if (!is.matrix(fh)) stop("wild-C12: fitted retornou objeto não-matriz e não-lista.")
    M <- fh
    if (is.null(colnames(M))) {
      if (ncol(M) != length(eq_safe)) stop("wild-C12: fitted sem colnames e ncol != #eqs_safe.")
      colnames(M) <- eq_safe
    }
    if (!all(eq_safe %in% colnames(M))) stop("wild-C12: fitted(matriz) não contém colnames compatíveis com eqs_safe.")
    M <- M[, eq_safe, drop = FALSE]
  }
  
  if (nrow(M) != nrow(df)) stop("wild-C12: desalinhamento nrow(fitted) != nrow(df).")
  
  colnames(M) <- eq_orig
  
  Wh <- matrix(NA_real_, nrow = nrow(df), ncol = length(shareNames), dimnames = list(NULL, shareNames))
  Wh[, eq_orig] <- M
  
  if (length(omit_share) == 1 && omit_share %in% shareNames) {
    Wh[, omit_share] <- 1 - rowSums(Wh[, eq_orig, drop=FALSE])
  }
  
  if (anyNA(Wh)) stop("wild-C12: NA em fitted shares após mapeamento/adding-up.")
  Wh
}

## signature-aware wrapper for elas_quaids_manual (if elas_at_safe not preloaded)
.call_with_formals_local <- function(fun, ...) {
  dots <- list(...)
  fmls <- names(formals(fun))
  keep <- intersect(names(dots), fmls)
  do.call(fun, dots[keep])
}

.estimate_elasticities <- function(fit, x_star, p_star,
                                   normalize_eval = TRUE, normalize_deriv = TRUE) {
  if (exists("elas_at_safe", mode = "function")) {
    return(elas_at_safe(
      fit, x = x_star, p = p_star,
      normalize_eval = normalize_eval,
      normalize_deriv = normalize_deriv
    ))
  }
  if (exists("elas_quaids_manual", mode = "function")) {
    Eb <- .call_with_formals_local(
      elas_quaids_manual,
      fit_q = fit,
      fit   = fit,
      x = x_star,
      p = p_star,
      dlogp = 1e-3,
      eps   = 1e-3,
      target_dw = 1e-3,
      normalize_eval  = if (isTRUE(normalize_eval)) "softmax" else "none",
      normalize_deriv = if (isTRUE(normalize_deriv)) "softmax" else "none",
      enforce_hicks = TRUE,
      enforce_symmetry = TRUE
    )
    return(Eb)
  }
  stop("wild-C12: não encontrei elas_at_safe() nem elas_quaids_manual().")
}

## abs-cap at replica-level (pre-store)
.cap_mat <- function(A, cap) {
  if (is.null(cap) || !is.finite(cap)) return(A)
  A2 <- A
  A2[is.finite(A2) & abs(A2) > cap] <- NA_real_
  A2
}

## winsor at summarization-level
.winsor_vec <- function(y, p) {
  if (is.null(p)) return(y)
  if (!is.numeric(p) || length(p) != 1 || !(p > 0 && p < 0.5)) return(y)
  lo <- as.numeric(stats::quantile(y, p,     na.rm = TRUE, type = 7))
  hi <- as.numeric(stats::quantile(y, 1 - p, na.rm = TRUE, type = 7))
  pmin(pmax(y, lo), hi)
}

.summ_vec <- function(x, winsor_p = NULL){
  ok <- is.finite(x)
  if (!any(ok)) return(c(mean=NA_real_, sd=NA_real_, q025=NA_real_, q975=NA_real_, n_ok=0))
  y <- x[ok]
  y <- .winsor_vec(y, winsor_p)
  c(
    mean = mean(y),
    sd   = stats::sd(y),
    q025 = as.numeric(stats::quantile(y, 0.025, na.rm=TRUE, type=7)),
    q975 = as.numeric(stats::quantile(y, 0.975, na.rm=TRUE, type=7)),
    n_ok = length(y)
  )
}

.summ_array_wide <- function(A, winsor_p = NULL){
  dn2 <- dimnames(A)[[2]]; dn3 <- dimnames(A)[[3]]
  K <- length(dn2)
  M_mean <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_q025 <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_q975 <- matrix(NA_real_, K, K, dimnames=list(dn2, dn3))
  M_nok  <- matrix(0L,       K, K, dimnames=list(dn2, dn3))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    s <- .summ_vec(A[,i,j], winsor_p = winsor_p)
    M_mean[i,j] <- s["mean"]
    M_q025[i,j] <- s["q025"]
    M_q975[i,j] <- s["q975"]
    M_nok[i,j]  <- as.integer(s["n_ok"])
  }
  list(mean = M_mean, q025 = M_q025, q975 = M_q975, n_ok = M_nok)
}

.write_wide <- function(M, out_csv){
  df <- as.data.frame(M, check.names = FALSE)
  df$good_i <- rownames(M)
  df <- df[, c("good_i", setdiff(names(df), "good_i")), drop = FALSE]
  utils::write.csv2(df, out_csv, row.names = FALSE)
  invisible(TRUE)
}

.combine_ci_long <- function(meanM, loM, hiM, out_csv){
  tab <- expand.grid(
    good_i = rownames(meanM),
    good_j = colnames(meanM),
    stringsAsFactors = FALSE
  )
  tab$estimate <- as.numeric(meanM[cbind(tab$good_i, tab$good_j)])
  tab$q025     <- as.numeric(loM  [cbind(tab$good_i, tab$good_j)])
  tab$q975     <- as.numeric(hiM  [cbind(tab$good_i, tab$good_j)])
  utils::write.csv2(tab, out_csv, row.names = FALSE)
  invisible(tab)
}

## ============================================================
## Main function
## ============================================================
wild_bootstrap_C12 <- function(df_iv, fit_q, iv_set, cluster_var,
                               B = 2000, seed = 20251221,
                               wild_type = c("rademacher","mammen"),
                               out_dir_tables = ".",
                               normalize_eval = TRUE,
                               normalize_deriv = TRUE,
                               eps_share = 1e-10,
                               ## NEW:
                               abs_cap_C12 = Inf,
                               winsor_p_C12 = NULL,
                               export_wide_triplet = TRUE) {
  
  wild_type <- match.arg(wild_type)
  set.seed(seed)
  
  stopifnot(is.data.frame(df_iv))
  if (!(cluster_var %in% names(df_iv))) stop("wild-C12: cluster_var não está em df_iv.")
  if (is.null(fit_q$fit) || !inherits(fit_q$fit, "systemfit")) stop("wild-C12: fit_q$fit não é systemfit.")
  if (is.null(fit_q$priceNames) || !length(fit_q$priceNames)) stop("wild-C12: fit_q$priceNames vazio.")
  if (is.null(fit_q$shareNames) || !length(fit_q$shareNames)) stop("wild-C12: fit_q$shareNames vazio.")
  if (!exists(".refit_quaids_any", mode = "function")) stop("wild-C12: falta .refit_quaids_any() no ambiente.")
  if (!exists("auto_eq_map2", mode = "function")) stop("wild-C12: falta auto_eq_map2() no ambiente.")
  if (is.null(iv_set) || !length(iv_set)) stop("wild-C12: iv_set vazio.")
  dir.create(out_dir_tables, showWarnings = FALSE, recursive = TRUE)
  
  priceNames0 <- as.character(fit_q$priceNames)
  shareNames0 <- as.character(fit_q$shareNames)
  drop_price  <- if (!is.null(fit_q$drop_price)) fit_q$drop_price else character(0)
  omit_share  <- if (!is.null(fit_q$omit_share)) fit_q$omit_share else character(0)
  endogs_all  <- paste0("ln_", setdiff(priceNames0, drop_price))
  exogs       <- intersect(c("z","z2"), names(df_iv))
  shifters    <- if (!is.null(fit_q$shifters)) intersect(as.character(fit_q$shifters), names(df_iv)) else character(0)
  
  x_name <- .resolve_x_name(fit_q, df_iv)
  
  ## align cc with fit (reduce nrow mismatches)
  eqs_est <- setdiff(shareNames0, omit_share)
  use_all <- unique(c(eqs_est, endogs_all, exogs, iv_set, shifters,
                      paste0("ln_", priceNames0), priceNames0, x_name, cluster_var))
  use_all <- use_all[use_all %in% names(df_iv)]
  df_cc <- df_iv[stats::complete.cases(df_iv[, use_all, drop = FALSE]), , drop = FALSE]
  if (nrow(df_cc) == 0) stop("wild-C12: complete.cases resultou em 0 linhas (use_all).")
  
  ## eval point
  x_star <- stats::median(df_cc[[x_name]], na.rm = TRUE)
  p_star <- if (exists("get_p_star", mode = "function")) {
    ps <- try(get_p_star(fit_q, df_cc, priceNames0), silent = TRUE)
    if (inherits(ps, "try-error")) .get_p_star_fallback(df_cc, priceNames0) else ps
  } else {
    .get_p_star_fallback(df_cc, priceNames0)
  }
  
  ## fitted shares + residuals
  W_hat <- .get_fitted_shares(fit_q, df_cc, shareNames0)
  W_obs <- as.matrix(df_cc[, shareNames0, drop = FALSE])
  U_hat <- W_obs - W_hat
  
  ## storage
  K0 <- length(priceNames0)
  acc_M <- array(NA_real_, dim = c(B, K0, K0), dimnames = list(NULL, priceNames0, priceNames0))
  acc_H <- array(NA_real_, dim = c(B, K0, K0), dimnames = list(NULL, priceNames0, priceNames0))
  
  ok_refit <- logical(B)
  ok_elas  <- logical(B)
  ok_store <- logical(B)
  
  ## NEW: cap counters
  cap_count_M <- 0L
  cap_count_H <- 0L
  
  err_log <- data.frame(draw = integer(0), stage = character(0), error = character(0), stringsAsFactors = FALSE)
  
  cl <- df_cc[[cluster_var]]
  
  for (b in seq_len(B)) {
    
    v <- .wild_weights_cluster(cl, type = wild_type)
    
    Wb <- W_hat + (U_hat * v)
    Wb <- .enforce_simplex(Wb, eps = eps_share)
    
    dfb <- df_cc
    for (j in seq_along(shareNames0)) dfb[[shareNames0[j]]] <- Wb[, j]
    
    ## refit (pipeline contract)
    fitb <- try(
      .refit_quaids_any(
        dfb,
        priceNames = priceNames0,
        shareNames = shareNames0,
        x_name     = x_name,
        instNames  = iv_set,
        shifters   = shifters
      ),
      silent = TRUE
    )
    if (inherits(fitb, "try-error")) {
      if (nrow(err_log) < 200) err_log <- rbind(err_log, data.frame(draw=b, stage="refit", error=as.character(fitb)))
      next
    }
    ok_refit[b] <- TRUE
    
    Eb <- try(.estimate_elasticities(fitb, x_star, p_star, normalize_eval, normalize_deriv), silent = TRUE)
    if (inherits(Eb, "try-error")) {
      if (nrow(err_log) < 200) err_log <- rbind(err_log, data.frame(draw=b, stage="elas", error=as.character(Eb)))
      next
    }
    ok_elas[b] <- TRUE
    
    if (is.null(Eb$marshall) || is.null(Eb$hicks) || !is.matrix(Eb$marshall) || !is.matrix(Eb$hicks)) {
      if (nrow(err_log) < 200) err_log <- rbind(err_log, data.frame(draw=b, stage="elas_shape", error="Eb$marshall/Eb$hicks não são matrizes"))
      next
    }
    
    ## KxK full, aligned by dimnames if present
    if (any(dim(Eb$marshall) != c(K0, K0)) || any(dim(Eb$hicks) != c(K0, K0))) {
      
      if (!is.null(rownames(Eb$marshall)) && !is.null(colnames(Eb$marshall)) &&
          all(priceNames0 %in% rownames(Eb$marshall)) && all(priceNames0 %in% colnames(Eb$marshall))) {
        M <- Eb$marshall[priceNames0, priceNames0, drop = FALSE]
      } else {
        if (nrow(err_log) < 200) err_log <- rbind(err_log, data.frame(draw=b, stage="store_M_dim", error="dim(M) != K0xK0 e sem dimnames compatíveis"))
        next
      }
      
      if (!is.null(rownames(Eb$hicks)) && !is.null(colnames(Eb$hicks)) &&
          all(priceNames0 %in% rownames(Eb$hicks)) && all(priceNames0 %in% colnames(Eb$hicks))) {
        H <- Eb$hicks[priceNames0, priceNames0, drop = FALSE]
      } else {
        if (nrow(err_log) < 200) err_log <- rbind(err_log, data.frame(draw=b, stage="store_H_dim", error="dim(H) != K0xK0 e sem dimnames compatíveis"))
        next
      }
      
    } else {
      M <- Eb$marshall
      H <- Eb$hicks
      dimnames(M) <- list(priceNames0, priceNames0)
      dimnames(H) <- list(priceNames0, priceNames0)
    }
    
    ## NEW: abs-cap pre-store (counts how many entries capped)
    if (is.finite(abs_cap_C12)) {
      cap_count_M <- cap_count_M + sum(is.finite(M) & abs(M) > abs_cap_C12)
      cap_count_H <- cap_count_H + sum(is.finite(H) & abs(H) > abs_cap_C12)
      M <- .cap_mat(M, abs_cap_C12)
      H <- .cap_mat(H, abs_cap_C12)
    }
    
    acc_M[b,,] <- M
    acc_H[b,,] <- H
    ok_store[b] <- TRUE
  }
  
  ## diagnostics
  diag_tbl <- data.frame(
    B = B,
    seed = seed,
    wild_type = wild_type,
    n_cc = nrow(df_cc),
    n_refit_ok = sum(ok_refit),
    n_elas_ok  = sum(ok_elas),
    n_store_ok = sum(ok_store),
    share_store_ok = mean(ok_store),
    x_name = x_name,
    x_star = x_star,
    abs_cap_C12 = abs_cap_C12,
    winsor_p_C12 = if (is.null(winsor_p_C12)) NA_real_ else winsor_p_C12,
    cap_count_M = cap_count_M,
    cap_count_H = cap_count_H,
    stringsAsFactors = FALSE
  )
  utils::write.csv2(diag_tbl, file.path(out_dir_tables, "Table_C12_wild_diagnostics.csv"), row.names = FALSE)
  if (nrow(err_log)) utils::write.csv2(err_log, file.path(out_dir_tables, "Table_C12_wild_errors_sample.csv"), row.names = FALSE)
  
  if (!any(ok_store)) stop("wild-C12: nenhuma réplica válida foi armazenada; veja Table_C12_wild_errors_sample.csv")
  
  ## summarize (winsor happens here)
  acc_M_ok <- acc_M[ok_store,,, drop = FALSE]
  acc_H_ok <- acc_H[ok_store,,, drop = FALSE]
  
  meanM <- apply(acc_M_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["mean"])
  loM   <- apply(acc_M_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["q025"])
  hiM   <- apply(acc_M_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["q975"])
  nokM  <- apply(acc_M_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["n_ok"])
  
  meanH <- apply(acc_H_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["mean"])
  loH   <- apply(acc_H_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["q025"])
  hiH   <- apply(acc_H_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["q975"])
  nokH  <- apply(acc_H_ok, c(2,3), function(v) .summ_vec(v, winsor_p = winsor_p_C12)["n_ok"])
  
  ## coerce to matrices with dimnames
  meanM <- matrix(as.numeric(meanM), K0, K0, dimnames=list(priceNames0, priceNames0))
  loM   <- matrix(as.numeric(loM),   K0, K0, dimnames=list(priceNames0, priceNames0))
  hiM   <- matrix(as.numeric(hiM),   K0, K0, dimnames=list(priceNames0, priceNames0))
  nokM  <- matrix(as.integer(nokM),  K0, K0, dimnames=list(priceNames0, priceNames0))
  
  meanH <- matrix(as.numeric(meanH), K0, K0, dimnames=list(priceNames0, priceNames0))
  loH   <- matrix(as.numeric(loH),   K0, K0, dimnames=list(priceNames0, priceNames0))
  hiH   <- matrix(as.numeric(hiH),   K0, K0, dimnames=list(priceNames0, priceNames0))
  nokH  <- matrix(as.integer(nokH),  K0, K0, dimnames=list(priceNames0, priceNames0))
  
  ## export: ONE TABLE (estimate/q025/q975) in long
  .combine_ci_long(meanM, loM, hiM, file.path(out_dir_tables, "Table_C12_wild_marshall_CI_long.csv"))
  .combine_ci_long(meanH, loH, hiH, file.path(out_dir_tables, "Table_C12_wild_hicks_CI_long.csv"))
  
  ## export: wide (optional, compatible with previous outputs)
  if (isTRUE(export_wide_triplet)) {
    .write_wide(meanM, file.path(out_dir_tables, "Table_C12_wild_marshall_mean_wide.csv"))
    .write_wide(loM,   file.path(out_dir_tables, "Table_C12_wild_marshall_q025_wide.csv"))
    .write_wide(hiM,   file.path(out_dir_tables, "Table_C12_wild_marshall_q975_wide.csv"))
    .write_wide(nokM,  file.path(out_dir_tables, "Table_C12_wild_marshall_nok_wide.csv"))
    
    .write_wide(meanH, file.path(out_dir_tables, "Table_C12_wild_hicks_mean_wide.csv"))
    .write_wide(loH,   file.path(out_dir_tables, "Table_C12_wild_hicks_q025_wide.csv"))
    .write_wide(hiH,   file.path(out_dir_tables, "Table_C12_wild_hicks_q975_wide.csv"))
    .write_wide(nokH,  file.path(out_dir_tables, "Table_C12_wild_hicks_nok_wide.csv"))
  }
  
  ## export robustness settings
  utils::write.csv2(
    data.frame(
      abs_cap_C12 = abs_cap_C12,
      winsor_p_C12 = if (is.null(winsor_p_C12)) NA_real_ else winsor_p_C12,
      note = "abs-cap applied pre-store; winsor applied in summarization only",
      stringsAsFactors = FALSE
    ),
    file.path(out_dir_tables, "Table_C12_wild_robustness_settings.csv"),
    row.names = FALSE
  )
  
  invisible(list(
    diagnostics = diag_tbl,
    errors = err_log,
    marshall = list(mean = meanM, q025 = loM, q975 = hiM, n_ok = nokM),
    hicks    = list(mean = meanH, q025 = loH, q975 = hiH, n_ok = nokH)
  ))
}

## ============================================================
## Example call
## ============================================================
res_wild <- wild_bootstrap_C12(
   df_iv = df_iv, fit_q = fit_q, iv_set = iv_set, cluster_var = cluster_var,
   B = 400, seed = 20251221,
   wild_type = "rademacher",
   out_dir_tables = out_dir_tables,
   abs_cap_C12 = 10,          # e.g., 10 (or Inf to disable)
   winsor_p_C12 = 0.01        # e.g., 0.01 (or NULL to disable)
)
## ============================================================
## ============================================================
## 15) TABLE C.13 – Parametric (Delta) CIs for elasticities + ETA (delta)
##     PIPELINE "OFICIAL" (reform_9.R / reform_10.R style), plug-and-play
##
##     IMPORTANTE (seu pedido):
##       - NÃO mexe em Hicks e Marshall (continuam Delta/Jacobian como estava)
##       - Corrige APENAS ETA: IC agora via SIMULAÇÃO paramétrica (percentil),
##         usando theta~N(theta_hat, Vth) para capturar não-linearidade/clamps
##
##     Outputs (out_dir_tables):
##       - Table_C13_delta_marshall_FULL.csv
##       - Table_C13_delta_hicks_FULL.csv
##       - Table_C13_delta_eta.csv
##       - Table_C13_delta_mapping_diag.csv
## ============================================================

options(scipen = 999)

`%||%` <- function(a,b) if(!is.null(a)) a else b

.need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
.need_pkg(c("glmnet","pls","systemfit","numDeriv","tibble","MASS"))

## ----------------------------
## 0) Aliases (compatibilidade com reform_9/10)
## ----------------------------
if (!exists("df_iv_ok", inherits=TRUE)) {
  if (exists("df_iv", inherits=TRUE)) {
    df_iv_ok <- get("df_iv", inherits=TRUE)
  } else {
    stop("C13: não encontrei df_iv_ok nem df_iv no ambiente.")
  }
}

if (!exists("best_fit", inherits=TRUE)) {
  if (exists("fit_q", inherits=TRUE)) {
    best_fit <- get("fit_q", inherits=TRUE)
  } else {
    stop("C13: não encontrei best_fit nem fit_q no ambiente.")
  }
}

if (!exists("best_iv_set", inherits=TRUE)) {
  if (exists("iv_set", inherits=TRUE)) {
    best_iv_set <- get("iv_set", inherits=TRUE)
  } else {
    stop("C13: não encontrei best_iv_set nem iv_set no ambiente.")
  }
}

if (!exists("out_dir_tables", inherits=TRUE)) {
  stop("C13: não encontrei out_dir_tables no ambiente.")
}

if (!exists("elas_quaids_manual", mode="function")) {
  stop("C13: elas_quaids_manual() não está no ambiente. Rode os steps do projeto que a definem antes.")
}

## ----------------------------
## 0b) priceNames/shareNames SEM "ternário"
## ----------------------------
priceNames <- NULL
if (exists("priceNames", inherits=TRUE)) {
  tmp <- get("priceNames", inherits=TRUE)
  if (is.character(tmp) && length(tmp) > 0) priceNames <- tmp
}
if (is.null(priceNames)) priceNames <- best_fit$priceNames

shareNames <- NULL
if (exists("shareNames", inherits=TRUE)) {
  tmp <- get("shareNames", inherits=TRUE)
  if (is.character(tmp) && length(tmp) > 0) shareNames <- tmp
}
if (is.null(shareNames)) shareNames <- best_fit$shareNames

drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share

if (is.null(priceNames) || !length(priceNames)) stop("C13: priceNames vazio.")
if (is.null(shareNames) || !length(shareNames)) stop("C13: shareNames vazio.")
if (length(drop_price) != 1L) stop("C13: drop_price deve ter comprimento 1.")
if (length(omit_share) != 1L) stop("C13: omit_share deve ter comprimento 1.")

## exogs (como no pipeline: usa apenas se existirem)
exogs <- intersect(c("z","z2"), names(df_iv_ok))

endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))

## ----------------------------
## 1) Dicionário expandido de IVs (reform_10 style)
## ----------------------------
augment_dictionary <- function(df, base_iv, with=c("SH_SH_area_1","SH_SH_capital_1"),
                               poly_deg=3, center=TRUE, scale.=TRUE){
  with <- intersect(with, names(df))
  out_names <- character(0)
  for (b in base_iv) if (b %in% names(df)) {
    v <- df[[b]]
    if (center) v <- v - mean(v, na.rm=TRUE)
    if (scale.) {
      s <- stats::sd(v, na.rm=TRUE)
      if (is.finite(s) && s > 0) v <- v / s
    }
    for (d in 1:poly_deg) {
      nm <- paste0(b,"_p",d); df[[nm]] <- as.numeric(v^d); out_names <- c(out_names, nm)
    }
    for (w in with) {
      nm <- paste0(b,"_x_",w); df[[nm]] <- as.numeric(df[[b]])*as.numeric(df[[w]])
      out_names <- c(out_names, nm)
    }
  }
  list(data=df, iv_names=unique(c(base_iv, out_names)))
}

base_iv <- intersect(best_iv_set, names(df_iv_ok))
if (!length(base_iv)) stop("C13: iv_pool vazio: best_iv_set não bate com colunas de df_iv_ok.")

aug <- augment_dictionary(
  df_iv_ok,
  base_iv = base_iv,
  with    = c("SH_SH_area_1","SH_SH_capital_1"),
  poly_deg= 3
)
df_big  <- aug$data
iv_pool <- aug$iv_names

## ----------------------------
## 2) SW-PCs por preço (reform_10 style)
## ----------------------------
resid_on <- function(v, X){
  if (is.null(X)) return(as.numeric(v))
  X <- as.matrix(X)
  if (!ncol(X)) return(as.numeric(v))
  v <- as.numeric(v)
  XtX <- crossprod(X)
  eps <- 1e-10
  beta <- tryCatch(solve(XtX, crossprod(X, v)), error=function(e){
    solve(XtX + diag(eps, ncol(X)), crossprod(X, v))
  })
  as.numeric(v - X %*% beta)
}
resid_cols_on <- function(Z, X){
  Z <- as.matrix(Z)
  if (is.null(X)) return(Z)
  X <- as.matrix(X)
  if (!ncol(X)) return(Z)
  XtX <- crossprod(X)
  eps <- 1e-10
  B <- tryCatch(solve(XtX, crossprod(X, Z)), error=function(e){
    solve(XtX + diag(eps, ncol(X)), crossprod(X, Z))
  })
  Z - X %*% B
}

make_sw_pc_for_price <- function(poi, df, iv_pool, endogs_all, exogs){
  stopifnot(poi %in% names(df))
  Znames <- intersect(iv_pool, names(df))
  Dminus <- setdiff(endogs_all, poi)
  Xn     <- intersect(exogs, names(df))
  use    <- unique(c(endogs_all, Xn, Znames))
  dat    <- df[stats::complete.cases(df[,use,drop=FALSE]), , drop=FALSE]
  
  Dhat <- NULL
  if (length(Dminus)) {
    XZ <- as.matrix(dat[, c(Xn, Znames), drop=FALSE])
    Dhat <- sapply(Dminus, function(d){
      y <- as.numeric(dat[[d]])
      fit <- glmnet::cv.glmnet(XZ, y, alpha=0, intercept=TRUE, standardize=TRUE)
      as.numeric(predict(fit, XZ, s="lambda.min"))
    })
    Dhat <- as.data.frame(Dhat)
    colnames(Dhat) <- paste0(Dminus,"_hat")
  }
  
  Xc <- if (!is.null(Dhat)) cbind(dat[, Xn, drop=FALSE], Dhat) else dat[, Xn, drop=FALSE]
  Dj_t <- resid_on(dat[[poi]], Xc)
  Z_t  <- resid_cols_on(dat[, Znames, drop=FALSE], Xc)
  
  fit <- pls::plsr(Dj_t ~ ., data=as.data.frame(Z_t), ncomp=1, validation="none", scale=TRUE)
  s1  <- drop(pls::scores(fit)[,1])
  nm  <- paste0("pc_", sub("^ln_", "", poi), "_sw1")
  dat[[nm]] <- s1
  
  list(data=dat, pc_name=nm)
}

gen_all_pcs_sw <- function(df, priceNames, drop_price, iv_pool, exogs){
  endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))
  out_df <- df; pc_names <- character(0)
  for (poi in endogs_all){
    pc <- make_sw_pc_for_price(poi, out_df, iv_pool, endogs_all, exogs)
    out_df   <- pc$data
    pc_names <- c(pc_names, pc$pc_name)
  }
  list(data=out_df, iv_star=pc_names)
}

pcs_sw <- gen_all_pcs_sw(df_big, priceNames, drop_price, iv_pool, exogs)
dfJI    <- pcs_sw$data
iv_star <- as.character(pcs_sw$iv_star)

stopifnot(length(iv_star) == length(endogs_all))
if (!all(iv_star %in% names(dfJI))) stop("C13: PCs não encontrados em dfJI (inconsistência).")

## ----------------------------
## 3) Systemfit 3SLS (reform_10 style)
## ----------------------------
eqs_orig  <- setdiff(shareNames, omit_share)
eqs_safe  <- paste0("eq", seq_along(eqs_orig))
orig2safe <- setNames(eqs_safe, eqs_orig)

rhs_y <- paste(c(endogs_all, exogs), collapse = " + ")
rhs_z <- paste(c(exogs, iv_star),  collapse = " + ")

eq_list   <- setNames(vector("list", length(eqs_orig)), eqs_safe)
inst_list <- setNames(vector("list", length(eqs_orig)), eqs_safe)

for (i in seq_along(eqs_orig)) {
  y_orig <- eqs_orig[i]
  y_safe <- orig2safe[[y_orig]]
  eq_list[[y_safe]]   <- stats::as.formula(paste(y_orig, "~", rhs_y))
  inst_list[[y_safe]] <- stats::as.formula(paste("~", rhs_z))
}

use_all <- unique(c(eqs_orig, endogs_all, exogs, iv_star))
dfsw <- dfJI[stats::complete.cases(dfJI[, use_all, drop=FALSE]), , drop=FALSE]
if (nrow(dfsw) == 0) stop("C13: complete.cases gerou 0 linhas em dfsw.")

fit_sys <- systemfit::systemfit(eq_list, method = "3SLS", inst = inst_list, data = dfsw)

## ----------------------------
## 4) Mapeamento/stack para Gamma
## ----------------------------
.canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x %||% ""))

diag_eq_map <- function(fit_sys){
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
    return(data.frame(eq_label=labs, lhs=lhs,
                      key_label=.canon(labs), key_lhs=.canon(lhs),
                      stringsAsFactors=FALSE))
  }
  cn <- names(coef(fit_sys))
  if (!length(cn)) return(data.frame(eq_label=character(0), lhs=character(0),
                                     key_label=character(0), key_lhs=character(0)))
  eq_guess  <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  data.frame(eq_label=eq_labels, lhs=NA_character_,
             key_label=.canon(eq_labels), key_lhs=NA_character_,
             stringsAsFactors=FALSE)
}

auto_eq_map2 <- function(fit_sys, shareNames, omit_share){
  diag <- diag_eq_map(fit_sys)
  shares_ok  <- setdiff(shareNames, omit_share)
  key_shares <- .canon(shares_ok)
  
  hitA <- if (!all(is.na(diag$lhs))) match(diag$lhs, shares_ok,  nomatch = 0L) else rep(0L, nrow(diag))
  hitB <- if (!all(is.na(diag$key_lhs))) match(diag$key_lhs, key_shares, nomatch = 0L) else rep(0L, nrow(diag))
  hitC <- match(diag$key_label, key_shares, nomatch = 0L)
  
  idx <- unique(c(which(hitA>0), which(hitB>0), which(hitC>0)))
  if (length(idx) > 0) {
    eqs_safe2 <- diag$eq_label[idx]
    eqs_orig2 <- shares_ok[
      ifelse(idx %in% which(hitA>0), hitA[idx],
             ifelse(idx %in% which(hitB>0), hitB[idx], hitC[idx]))
    ]
    dup <- duplicated(eqs_orig2)
    if (any(dup)) { eqs_safe2 <- eqs_safe2[!dup]; eqs_orig2 <- eqs_orig2[!dup] }
    return(list(eqs_safe=eqs_safe2, eqs_orig=eqs_orig2, diag=diag))
  }
  
  cn <- names(coef(fit_sys))
  if (!length(cn)) stop("auto_eq_map2: não foi possível inspecionar coeficientes do systemfit.")
  eq_guess  <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  Ke <- length(eq_labels)
  if (Ke != length(shares_ok)) stop(sprintf("auto_eq_map2: fallback por ordem falhou: #eq=%d != #shares=%d", Ke, length(shares_ok)))
  message("[auto_eq_map2] Usando fallback por ORDEM: ",
          paste(sprintf("%s↔%s", eq_labels, shares_ok), collapse="; "))
  list(eqs_safe=eq_labels, eqs_orig=shares_ok,
       diag=data.frame(eq_label=eq_labels, lhs=NA_character_,
                       key_label=.canon(eq_labels), key_lhs=NA_character_))
}

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
    if (is.na(idx)) stop(sprintf("Coeficiente não encontrado: eq='%s' var='%s'", eq_safe[i], pj[j]))
    theta[pos] <- unname(co[idx]); nm_th[pos] <- paste(eq_orig[i], pj[j], sep="::")
  }
  
  for (i in seq_along(eq_safe)) for (k in seq_along(eq_safe)) {
    for (a in seq_along(pj)) for (b in seq_along(pj)) {
      pos1 <- (i - 1L) * Kp + a; pos2 <- (k - 1L) * Kp + b
      idx1 <- find_idx(eq_safe[i], pj[a]); idx2 <- find_idx(eq_safe[k], pj[b])
      Vth[pos1,pos2] <- V[idx1, idx2]
    }
  }
  
  names(theta) <- nm_th
  list(theta=theta, V=Vth)
}

## ----------------------------
## 5) Elasticidades estáveis (Delta) — clamps/caps
## ----------------------------
softmax <- function(v){
  v <- as.numeric(v)
  m <- max(v)
  ev <- exp(v - m)
  ev / sum(ev)
}
.abs_cap_vec <- function(x, cap) pmax(pmin(x, cap), -cap)
.abs_cap_mat <- function(M, cap){ M[] <- .abs_cap_vec(as.numeric(M), cap); M }

elas_quaids_manual_delta_stable <- function(fit_q, x, p,
                                            eps = 1e-3,
                                            w_min = 1e-6,
                                            abs_cap_eta = 20,
                                            abs_cap_elas = 10){
  stopifnot(inherits(fit_q, "quaids_manual"))
  priceNames <- fit_q$priceNames
  shareNames <- fit_q$shareNames
  K <- length(priceNames)
  
  if (is.null(names(p))) names(p) <- priceNames
  p <- as.numeric(p[priceNames])
  if (length(p) != K) stop("C13 stable: p com dimensão errada.")
  
  alpha  <- fit_q$coef$alpha
  beta   <- fit_q$coef$beta
  lambda <- fit_q$coef$lambda
  gamma  <- fit_q$coef$gamma
  
  wbar <- colMeans(fit_q$data[, shareNames, drop=FALSE], na.rm=TRUE)
  wbar[!is.finite(wbar)] <- 0
  if (sum(wbar) <= 0) wbar <- rep(1/length(wbar), length(wbar))
  wbar <- wbar / sum(wbar)
  
  w_hat_fun <- function(p_vec, x_val){
    lp <- log(p_vec)
    if (fit_q$priceIndex == "Ls") {
      z <- log(x_val) - sum(wbar * lp)
      w_raw <- as.numeric(alpha + gamma %*% lp + beta * z + lambda * z^2)
      w <- softmax(w_raw)
    } else if (fit_q$priceIndex == "S") {
      z_ls  <- log(x_val) - sum(wbar * lp)
      w1raw <- as.numeric(alpha + gamma %*% lp + beta * z_ls + lambda * z_ls^2)
      w1    <- softmax(w1raw)
      z_st  <- log(x_val) - sum(w1 * lp)
      w_raw <- as.numeric(alpha + gamma %*% lp + beta * z_st + lambda * z_st^2)
      w <- softmax(w_raw)
    } else stop("C13 stable: priceIndex desconhecido.")
    
    w <- pmax(w, w_min)
    w <- w / sum(w)
    names(w) <- shareNames
    w
  }
  
  w0 <- w_hat_fun(p, x)
  if (any(!is.finite(w0))) stop("C13 stable: w0 não-finito.")
  
  z0 <- if (fit_q$priceIndex == "Ls") log(x) - sum(wbar * log(p)) else log(x) - sum(w0 * log(p))
  
  eta <- 1 + (beta + 2*lambda*z0) / pmax(w0, w_min)
  eta <- .abs_cap_vec(eta, abs_cap_eta)
  names(eta) <- shareNames
  
  q0 <- (w0 * x) / p
  q_floor <- 1e-12
  
  E_M <- matrix(NA_real_, K, K, dimnames=list(shareNames, priceNames))
  for (j in seq_len(K)) {
    p_up <- p
    p_up[j] <- p_up[j] * (1 + eps)
    w_up <- w_hat_fun(p_up, x)
    q_up <- (w_up * x) / p_up
    
    den_logp <- (log(p_up[j]) - log(p[j]))
    
    if (any(!is.finite(q0)) || any(!is.finite(q_up)) ||
        any(q0 <= q_floor) || any(q_up <= q_floor)) {
      den_lvlp <- (p_up[j] - p[j]) / p[j]
      E_M[, j] <- ((q_up - q0) / pmax(q0, q_floor)) / den_lvlp
    } else {
      E_M[, j] <- (log(q_up) - log(q0)) / den_logp
    }
  }
  
  E_M <- .abs_cap_mat(E_M, abs_cap_elas)
  E_H <- E_M + outer(eta, w0)
  E_H <- .abs_cap_mat(E_H, abs_cap_elas)
  
  if (any(!is.finite(E_M)) || any(!is.finite(E_H)) || any(!is.finite(eta))) {
    stop("C13 stable: apareceu NA/Inf mesmo após clamps/caps.")
  }
  
  list(expenditure=eta, marshall=E_M, hicks=E_H,
       at=list(x=x, p=setNames(p, priceNames), w=w0, z=z0))
}

## ----------------------------
## 6) Delta CIs (Hicks FULL, Marshall FULL) + ETA (SIMULAÇÃO)
## ----------------------------
elasticity_CIs_delta_system_C13 <- function(
    fit_sys, fit_template, dfsw,
    priceNames, shareNames, drop_price, omit_share,
    level = 0.95, method = c("pointwise","bonferroni"),
    eps = 1e-3, w_min = 1e-6, abs_cap_eta = 20, abs_cap_elas = 10,
    B_eta = 2000, seed_eta = 123
){
  method <- match.arg(method)
  
  map <- auto_eq_map2(fit_sys, shareNames, omit_share)
  eqs_safe2 <- map$eqs_safe
  eqs_orig2 <- map$eqs_orig
  
  pj       <- paste0("ln_", setdiff(priceNames, drop_price))
  pj_short <- sub("^ln_", "", pj)
  meta     <- stack_gamma_vcov_sys2(fit_sys, pj, eqs_safe2, eqs_orig2)
  theta    <- meta$theta
  Vth      <- meta$V
  
  theta_to_G <- function(th){
    Ke <- length(eqs_orig2); Kp <- length(pj)
    Ghat <- matrix(th, nrow=Ke, byrow=TRUE, dimnames=list(eqs_orig2, pj))
    Gfull <- matrix(NA_real_, nrow=length(shareNames), ncol=Kp, dimnames=list(shareNames, pj))
    Gfull[eqs_orig2, pj] <- Ghat
    Gfull[omit_share, pj] <- -colSums(Ghat, na.rm=TRUE)
    Gfull
  }
  
  x_star <- if ("gasto_total_atualhat" %in% names(dfsw))
    median(dfsw$gasto_total_atualhat, na.rm=TRUE) else
      median(dfsw[[grep("gasto|exp", names(dfsw), value=TRUE)[1]]], na.rm=TRUE)
  
  p_star <- exp(colMeans(dfsw[paste0("ln_", priceNames)], na.rm=TRUE))
  p_star <- setNames(as.numeric(p_star), priceNames)
  
  message(sprintf("[C13] eval-point source=dfsw | x_star=%.6g | min(p_star)=%.6g",
                  x_star, min(p_star, na.rm=TRUE)))
  
  cols_fit <- colnames(fit_template$coef$gamma)
  
  apply_G_to_fit <- function(fit_obj, Gfull){
    f3 <- fit_obj
    f3$coef$gamma[, match(pj_short, cols_fit)] <- Gfull
    jdrop_ix <- match(drop_price, colnames(f3$coef$gamma))
    f3$coef$gamma[!is.finite(f3$coef$gamma)] <- 0
    if (is.finite(jdrop_ix) && !is.na(jdrop_ix)) {
      keep_ix <- setdiff(seq_len(ncol(f3$coef$gamma)), jdrop_ix)
      f3$coef$gamma[, jdrop_ix] <- -rowSums(f3$coef$gamma[, keep_ix, drop=FALSE])
    }
    f3
  }
  
  ## dimensão/nomes nativos
  G0 <- theta_to_G(theta)
  f2 <- apply_G_to_fit(fit_template, G0)
  
  E0 <- elas_quaids_manual_delta_stable(
    f2, x=x_star, p=p_star, eps=eps, w_min=w_min,
    abs_cap_eta=abs_cap_eta, abs_cap_elas=abs_cap_elas
  )
  
  rn0 <- rownames(E0$hicks); cn0 <- colnames(E0$hicks)
  if (is.null(rn0)) rn0 <- shareNames
  if (is.null(cn0)) cn0 <- priceNames
  k0 <- length(rn0)
  
  ## g(theta) (mantém como estava: vec(H), vec(M), eta)
  g_eval <- function(th){
    f3 <- apply_G_to_fit(fit_template, theta_to_G(th))
    E  <- elas_quaids_manual_delta_stable(
      f3, x=x_star, p=p_star, eps=eps, w_min=w_min,
      abs_cap_eta=abs_cap_eta, abs_cap_elas=abs_cap_elas
    )
    c(as.vector(t(E$hicks)),
      as.vector(t(E$marshall)),
      as.numeric(E$expenditure))
  }
  
  ## --- Delta para TODOS os elementos (mantém Hicks/Marshall como você tinha) ---
  g0 <- g_eval(theta)
  J  <- numDeriv::jacobian(g_eval, theta)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  
  lo <- g0 - z*se
  hi <- g0 + z*se
  
  idx_H    <- seq_len(k0*k0)
  idx_Mall <- (k0*k0) + seq_len(k0*k0)
  idx_eta  <- (2*k0*k0) + seq_len(k0)
  
  ## Hicks/Marshall (DELTA) — NÃO ALTERAR
  H_hat0 <- matrix(g0[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_lo0  <- matrix(lo[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_hi0  <- matrix(hi[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_se0  <- matrix(se[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  M_hat0 <- matrix(g0[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_lo0  <- matrix(lo[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_hi0  <- matrix(hi[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_se0  <- matrix(se[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  ## ETA (SIMULAÇÃO paramétrica) — ÚNICA MUDANÇA RELEVANTE
  set.seed(seed_eta)
  theta_draws <- MASS::mvrnorm(n = B_eta, mu = theta, Sigma = Vth)
  
  eta_draws <- matrix(NA_real_, nrow = B_eta, ncol = k0)
  colnames(eta_draws) <- rn0
  for (b in seq_len(B_eta)) {
    gb <- g_eval(theta_draws[b, ])
    eta_draws[b, ] <- gb[idx_eta]
  }
  
  if (method == "bonferroni") {
    a2 <- alpha / k0
    qlo <- a2/2
    qhi <- 1 - a2/2
  } else {
    qlo <- alpha/2
    qhi <- 1 - alpha/2
  }
  
  eta_hat0 <- as.numeric(g0[idx_eta]); names(eta_hat0) <- rn0
  eta_lo0  <- as.numeric(apply(eta_draws, 2, stats::quantile, probs=qlo, na.rm=TRUE)); names(eta_lo0) <- rn0
  eta_hi0  <- as.numeric(apply(eta_draws, 2, stats::quantile, probs=qhi, na.rm=TRUE)); names(eta_hi0) <- rn0
  eta_se0  <- as.numeric(apply(eta_draws, 2, stats::sd, na.rm=TRUE)); names(eta_se0) <- rn0
  
  to_tbl <- function(M, lo, hi, se, rows, cols){
    do.call(rbind, lapply(seq_along(rows), function(i)
      data.frame(i=rows[i], j=cols, hat=M[i,], lo=lo[i,], hi=hi[i,], se=se[i,], row.names=NULL)))
  }
  
  list(
    hicks_full_native    = to_tbl(H_hat0, H_lo0, H_hi0, H_se0, rn0, cn0),
    marshall_full_native = to_tbl(M_hat0, M_lo0, M_hi0, M_se0, rn0, cn0),
    eta_native = tibble::tibble(good=rn0, eta_hat=eta_hat0, eta_lo=eta_lo0, eta_hi=eta_hi0, eta_se=eta_se0),
    meta = list(level=level, method=method, eps=eps, w_min=w_min,
                abs_cap_eta=abs_cap_eta, abs_cap_elas=abs_cap_elas,
                B_eta=B_eta, seed_eta=seed_eta,
                mapping_diag=map$diag)
  )
}

## ----------------------------
## 7) Run + exports
## ----------------------------
cis_sys <- elasticity_CIs_delta_system_C13(
  fit_sys      = fit_sys,
  fit_template = best_fit,
  dfsw         = dfsw,
  priceNames   = priceNames,
  shareNames   = shareNames,
  drop_price   = drop_price,
  omit_share   = omit_share,
  level        = 0.95,
  method       = "bonferroni",
  eps          = 1e-3,
  w_min        = 1e-6,
  abs_cap_eta  = 20,
  abs_cap_elas = 10,
  B_eta        = 2000,
  seed_eta     = 123
)

utils::write.csv2(cis_sys$marshall_full_native, file.path(out_dir_tables, "Table_C13_delta_marshall_FULL.csv"), row.names = FALSE)
utils::write.csv2(cis_sys$hicks_full_native,    file.path(out_dir_tables, "Table_C13_delta_hicks_FULL.csv"),    row.names = FALSE)
utils::write.csv2(cis_sys$eta_native,           file.path(out_dir_tables, "Table_C13_delta_eta.csv"),           row.names = FALSE)
utils::write.csv2(as.data.frame(cis_sys$meta$mapping_diag), file.path(out_dir_tables, "Table_C13_delta_mapping_diag.csv"), row.names = FALSE)

cat("\n[C13] OK: exports escritos em out_dir_tables:\n",
    " - Table_C13_delta_marshall_FULL.csv\n",
    " - Table_C13_delta_hicks_FULL.csv\n",
    " - Table_C13_delta_eta.csv (ETA via simulação)\n",
    " - Table_C13_delta_mapping_diag.csv\n", sep="")

## --- Figura C.2 – Heatmap de elasticidades (ASSUNÇÃO) ---
## ASSUNÇÃO: a figura usa as elasticidades marshallianas “pontuais”
##           (paramétricas ou bootstrap). Aqui uso as paramétricas
##           (out_delta$marshallian) por serem diretamente associadas
##           ao “modelo final”.

if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)

heat_B5 <- cis_sys$marshall_full_native # marshallian: colunas i, j, est, se, lwr, upr
## ASSUMÇÃO: colunas i e j ou good_i e good_j estão presentes;
## se não, adapte conforme a estrutura do seu out_delta.

if (!"good_i" %in% names(heat_B5) && "i" %in% names(heat_B5)) {
  heat_B5$good_i <- as.factor(heat_B5$i)
}
if (!"good_j" %in% names(heat_B5) && "j" %in% names(heat_B5)) {
  heat_B5$good_j <- as.factor(heat_B5$j)
}

heat_B5$sig <- ifelse(heat_B5$q025<=0 & 0<=heat_B5$q975, NA, heat_B5$estimate)

p_B5 <- ggplot(heat_B5, aes(x = good_j, y = good_i, fill = sig)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(    low = "#2166AC",
                           mid = "white",
                           high = "#B2182B",
                           midpoint = -0.82,
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

ggsave(file.path(out_dir_figures, "Figure_C2_heatmap_elasticities.png"), p_B5, dpi = 300, width = 7, height = 5)

log_line("=== DONE Appendix C export ===")
log_line("Tables:  ", normalizePath(out_dir_tables))
log_line("Figures: ", normalizePath(out_dir_figures))





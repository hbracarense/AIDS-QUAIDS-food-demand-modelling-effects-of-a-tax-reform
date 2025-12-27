## =========================================================
##  TOP-TIER: J robusto por eq + C-tests por blocos + poda
## =========================================================
suppressPackageStartupMessages({
  library(AER); library(sandwich); library(dplyr); library(tidyr); library(openxlsx)
})

## ---------- helpers ----------
.mm <- function(formula_or_terms, data, add_intercept = TRUE) {
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && stats::sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}
.solve_psd <- function(S) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric = TRUE)
  pos <- ee$values > max(ee$values, 0, na.rm=TRUE) * 1e-10
  if (!any(pos)) return(MASS::ginv(S))
  ee$vectors[, pos, drop=FALSE] %*% diag(1/ee$values[pos], sum(pos)) %*% t(ee$vectors[, pos, drop=FALSE])
}
J_hansen <- function(uhat, Z, cluster = NULL) {
  # remove intercept de Z (se existir) para alinhar com testes padrão
  Z <- Z[, setdiff(colnames(Z), "(Intercept)"), drop = FALSE]
  n  <- nrow(Z)
  g  <- Z * as.numeric(uhat)              # N x L, sem reciclagem
  if (!is.null(cluster)) {
    # cluster-robusto "tipo Arellano": soma g por cluster
    cg <- rowsum(g, group = cluster, reorder = FALSE)
    S  <- crossprod(cg) / n
    gbar <- colMeans(g)
  } else {
    S  <- crossprod(g) / n
    gbar <- colMeans(g)
  }
  Winv <- .solve_psd(S)
  as.numeric(n * t(gbar) %*% Winv %*% gbar)
}
rank_dfJ <- function(X, Z) {
  rX <- qr(X)$rank
  rZ <- qr(Z)$rank
  c(rX = rX, rZ = rZ, dfJ = max(rZ - rX, 0))
}

## ---------- 1) J por equação (Sargan e Hansen robusto/cluster) ----------
jtests_by_eq_fixed <- function(fit_obj, data, cluster_var = NULL) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  inst_terms <- unique(fit_obj$inst_terms)
  eqs <- fit_obj$eq
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  out <- lapply(seq_along(eqs), function(i){
    f_yx <- formula(eqs[[i]])
    y    <- model.frame(f_yx, data = data)[[1]]
    X    <- .mm(delete.response(terms(f_yx)), data, add_intercept = TRUE)
    Z    <- .mm(inst_terms, data, add_intercept = TRUE)
    
    rk   <- rank_dfJ(X, Z)
    # 2SLS para obter resíduos (alinha com J padrão)
    ivm  <- AER::ivreg(f_yx, instruments = reformulate(inst_terms), data = data)
    uh   <- residuals(ivm)
    
    # Sargan homocedástico (se sobre-ID)
    sarg <- try(AER::sargan(ivm), silent = TRUE)
    J_s  <- if (!inherits(sarg,"try-error")) unname(as.numeric(sarg$statistic)) else NA_real_
    p_s  <- if (!inherits(sarg,"try-error")) unname(as.numeric(sarg$p.value))   else NA_real_
    
    # Hansen-J robusto (HC / cluster-robusto)
    J_h  <- if (rk["dfJ"] > 0) J_hansen(uh, Z, cluster = cl) else NA_real_
    p_h  <- if (rk["dfJ"] > 0) pchisq(J_h, df = rk["dfJ"], lower.tail = FALSE) else NA_real_
    
    data.frame(eq = as.character(f_yx[[2]]),
               rank_X = rk["rX"], rank_Z = rk["rZ"], df_J = rk["dfJ"],
               J_sargan = J_s, p_sargan = p_s,
               J_hansen = J_h, p_hansen = p_h,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

## ---------- 2) C-tests (diferença de Hansen) por blocos de IVs ----------
make_default_blocks <- function(fit_obj) {
  iv <- unique(fit_obj$inst_terms)
  rhs <- unique(fit_obj$rhs_terms)
  iv_excl <- setdiff(iv, intersect(iv, rhs))      # somente IVs excluídos
  # agrupamentos comuns no seu setup:
  blk <- list(
    z       = intersect(iv_excl, iv[grepl("^z$", iv)]),
    z2      = intersect(iv_excl, iv[grepl("z2|\\^2", iv)]),
    region  = intersect(iv_excl, iv[grepl("^f_reg", iv)]),
    area    = intersect(iv_excl, iv[grepl("^f_area", iv)]),
    capital = intersect(iv_excl, iv[grepl("^f_cap", iv)]),
    demo    = intersect(iv_excl, iv[grepl("p_n_mais65anos|p_sexofem_c|p_n_0esc", iv)]),
    inter   = intersect(iv_excl, iv[grepl(":", iv)])
  )
  # remove blocos vazios
  blk[ lengths(blk) > 0 ]
}

library(dplyr)

# ---- DROP-IN FIX ----
c_tests_blocks <- function(fit_obj, data,
                           blocks = make_default_blocks(fit_obj),
                           cluster_var = NULL) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  inst_all <- unique(fit_obj$inst_terms)
  eqs <- fit_obj$eq
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  rows <- vector("list", 0L)
  
  for (i in seq_along(eqs)) {
    f_yx   <- formula(eqs[[i]])
    eqname <- as.character(f_yx[[2]])
    
    X     <- .mm(delete.response(terms(f_yx)), data, add_intercept = TRUE)
    Zfull <- .mm(inst_all, data, add_intercept = TRUE)
    
    ivm <- AER::ivreg(f_yx, instruments = reformulate(inst_all), data = data)
    u   <- residuals(ivm)
    
    rkF <- rank_dfJ(X, Zfull)
    JF  <- if (rkF["dfJ"] > 0) J_hansen(u, Zfull, cluster = cl) else NA_real_
    
    for (bn in names(blocks)) {
      inst_restr <- setdiff(inst_all, blocks[[bn]])
      Zr  <- .mm(inst_restr, data, add_intercept = TRUE)
      rkR <- rank_dfJ(X, Zr)
      
      Jr <- if (rkR["dfJ"] > 0) {
        ivm_r <- AER::ivreg(f_yx, instruments = reformulate(inst_restr), data = data)
        J_hansen(residuals(ivm_r), Zr, cluster = cl)
      } else NA_real_
      
      Cstat <- if (is.finite(JF) && is.finite(Jr)) JF - Jr else NA_real_
      ddf   <- if (rkF["dfJ"] > rkR["dfJ"]) rkF["dfJ"] - rkR["dfJ"] else NA_real_
      pC    <- if (is.finite(Cstat) && is.finite(ddf) && ddf > 0) pchisq(Cstat, df = ddf, lower.tail = FALSE) else NA_real_
      
      rows[[length(rows) + 1L]] <- data.frame(
        eq = eqname,
        block = bn,
        k_removed = length(blocks[[bn]]),
        J_full = JF, dfJ_full = rkF["dfJ"],
        J_restr = Jr, dfJ_restr = rkR["dfJ"],
        C_stat = Cstat, df_C = ddf, p_C = pC,
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- if (length(rows)) dplyr::bind_rows(rows) else
    tibble::tibble(eq = character(), block = character(), k_removed = integer(),
                   J_full = numeric(), dfJ_full = integer(),
                   J_restr = numeric(), dfJ_restr = integer(),
                   C_stat = numeric(), df_C = integer(), p_C = numeric())
  out
}

# ---- roda novamente e ordena de forma à prova de falhas ----
blocks <- make_default_blocks(fit_q)

ctab <- c_tests_blocks(
  fit_q, df_iv, blocks = blocks,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"
)

# arrange seguro (se não tiver 'eq', ordena só por p_C)
ctab <- if ("eq" %in% names(ctab)) {
  dplyr::arrange(ctab, .data$eq, .data$p_C)
} else {
  dplyr::arrange(ctab, .data$p_C)
}



## ---------- 3) Poda sugerida a partir dos C-tests ----------
suggest_prune <- function(ctab, alpha = 0.05) {
  bad <- ctab %>%
    filter(!is.na(p_C), p_C < alpha) %>%
    count(block, sort = TRUE)
  bad$block
}
new_inst_after_prune <- function(fit_obj, blocks_to_drop, blocks = make_default_blocks(fit_obj)) {
  inst_all <- unique(fit_obj$inst_terms)
  drop_terms <- unique(unlist(blocks[blocks_to_drop]))
  setdiff(inst_all, drop_terms)
}

## ---------- 4) Export prático (xlsx) ----------
export_top_tier_diags <- function(path_xlsx, jtab, ctab, fs_tab, pR2_tab) {
  wb <- createWorkbook()
  addWorksheet(wb, "J_tests")
  writeData(wb, "J_tests", jtab)
  addWorksheet(wb, "C_tests_blocks")
  writeData(wb, "C_tests_blocks", ctab)
  addWorksheet(wb, "FirstStage_F")
  writeData(wb, "FirstStage_F", fs_tab)
  addWorksheet(wb, "Partial_R2")
  writeData(wb, "Partial_R2", pR2_tab)
  saveWorkbook(wb, path_xlsx, overwrite = TRUE)
  normalizePath(path_xlsx)
}

## =======================
## COMO USAR (plug-and-play)
## =======================

## 1) J por equação (robusto). Se tiver 'psu' use cluster_var="psu"; caso contrário "uf".
jtab <- jtests_by_eq_fixed(fit_q, df_iv, cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf")
print(jtab)

## 2) C-tests por blocos (diferença de Hansen)
blocks <- make_default_blocks(fit_q)
ctab   <- c_tests_blocks(fit_q, df_iv, blocks = blocks,
                         cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf")
ctab <- ctab %>% arrange(eq, p_C)
print(ctab)

## 3) Sugerir poda e produzir nova lista de IVs (para você refitar o sistema)
to_drop <- suggest_prune(ctab, alpha = 0.05)
cat("\nBlocos sugeridos para poda (p_C < 0.05):", paste(to_drop, collapse=", "), "\n")
inst_new <- new_inst_after_prune(fit_q, to_drop, blocks)
cat("\nInst terms (NOVOS) após poda:\n", paste(inst_new, collapse = " + "), "\n")

## 4) (Opcional) Rode novamente seu estimador do sistema com 'inst_new'
##    Ex.: fit_q2 <- fit_quaids_manual_km1(..., inst_terms = inst_new)
##    Depois repita J/C-tests para mostrar que não rejeita Hansen-J.

## 5) Exportar tudo
## (use suas tabelas já calculadas para força de IVs)
fs_tab <- fs_block_F(fit_q, df_iv, robust = TRUE)
pR2_tab <- partial_R2_exclIV(fit_q, df_iv)
#xlsx_path <- export_top_tier_diags("diagnosticos_top_tier.xlsx", jtab, ctab, fs_tab, pR2_tab)
#cat("\nArquivo salvo em:", xlsx_path, "\n")

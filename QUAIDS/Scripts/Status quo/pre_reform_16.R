# ============================================================
# QU-AIDS — JUST-ID na eq.6 (UF) + Bootstrap por UF (robusto)
# Saídas: Tabelas Marshall/Hicks/Renda com SE & IC (percentis)
# Requer: fit_q, df_iv, get_base_exog, refit_quaids_try_eqlist, quaids_elasticities
# ============================================================
library(openxlsx)
suppressWarnings({
  if (!requireNamespace("AER", quietly = TRUE)) stop("Precisa do pacote 'AER'.")
})

# -------------------------------
# 0) Utils
# -------------------------------
.nice_intersect <- function(x, y) { z <- intersect(x, y); z[!duplicated(z)] }

# -------------------------------
# 1) Força JUST-ID na equação alvo derrubando o mínimo de IVs
# -------------------------------
make_eq_justID2 <- function(fit_base, data, fit_eqspec, eqn,
                            prefer_keep = character(0),
                            avoid = character(0),
                            verbose = TRUE) {
  data <- as.data.frame(data)
  base_exog <- get_base_exog(fit_base)
  
  inst_by_eq <- attr(fit_eqspec, "inst_excl_by_eq")
  stopifnot(is.list(inst_by_eq), eqn %in% names(inst_by_eq))
  
  # fórmula da eq-alvo a partir de fit_base
  eq_names <- vapply(fit_base$eq, function(m) as.character(formula(m)[[2]]), character(1))
  f <- formula(fit_base$eq[[ match(eqn, eq_names) ]])
  
  cur <- unique(inst_by_eq[[eqn]])
  if (length(avoid)) cur <- setdiff(cur, avoid)
  
  dfJ_of <- function(S) {
    terms <- unique(c(base_exog, S))
    ivm <- try(AER::ivreg(f, instruments = reformulate(terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) return(NA_real_)
    X <- try(stats::model.matrix(ivm, component = "regressors"),  silent = TRUE)
    Z <- try(stats::model.matrix(ivm, component = "instruments"), silent = TRUE)
    if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NA_real_)
    qr(Z)$rank - qr(X)$rank
  }
  
  cur_dfJ <- dfJ_of(cur)
  if (isTRUE(cur_dfJ == 0)) {
    if (verbose) message("[", eqn, "] já está JUST/UNDER-ID (df_J <= 0). Nada a fazer.")
    return(fit_eqspec)
  }
  
  keep_best <- NULL
  for (k in seq_len(length(cur))) {
    for (drop_set in utils::combn(cur, k, simplify = FALSE)) {
      S <- setdiff(cur, drop_set)
      if (length(prefer_keep) && any(!(prefer_keep %in% S))) next
      dj <- dfJ_of(S)
      if (isTRUE(dj == 0)) { keep_best <- S; break }
    }
    if (!is.null(keep_best)) break
  }
  
  if (is.null(keep_best)) {
    if (verbose) message("[", eqn, "] não encontrei subset JUST-ID. Mantendo conjunto atual.")
    return(fit_eqspec)
  }
  
  inst_by_eq[[eqn]] <- unique(keep_best)
  if (verbose) {
    dropped <- setdiff(cur, keep_best)
    msg <- if (length(dropped)) paste(paste0("-", dropped), collapse = " ") else "(sem drops)"
    message("[", eqn, "] JUST-ID por drop: ", msg)
  }
  refit_quaids_try_eqlist(fit_base, data, inst_by_eq)
}

# -------------------------------
# 2) Bootstrap clusterizado por UF (Poisson+1), com filtros
#    - descarta réplicas com |elasticidade| > abs_cap (default 20)
#    - winsoriza colunas das réplicas (opcional)
# -------------------------------
.pack_all_elas <- function(fobj, at = "medians", w_source = "observed") {
  el <- quaids_elasticities(fobj, at = at, w_source = w_source)
  c(as.vector(as.matrix(el$marshallian)),
    as.vector(as.matrix(el$hicksian)),
    as.numeric(el$expenditure))
}

.fix_levels_like <- function(d, template) {
  fcs <- intersect(names(d), names(template))
  for (nm in fcs) if (is.factor(template[[nm]])) {
    d[[nm]] <- factor(d[[nm]], levels = levels(template[[nm]]))
  }
  d
}

.winsorize_cols <- function(M, p = 0.01) {
  if (p <= 0) return(M)
  for (j in seq_len(ncol(M))) {
    col <- M[, j]
    if (all(is.finite(col))) {
      lo <- stats::quantile(col, p,     type = 7, na.rm = TRUE)
      hi <- stats::quantile(col, 1 - p, type = 7, na.rm = TRUE)
      M[, j] <- pmin(pmax(col, lo), hi)
    }
  }
  M
}

boot_quaids_elasticities_byUF <- function(
    fit_base, data, fit_eqspec,
    cluster_var = "uf",
    R = 1000, seed = 42, verbose = TRUE,
    at = "medians", w_source = "observed",
    abs_cap = 20,         # filtro de plausibilidade |E| <= abs_cap
    winsor_p = 0.01       # winsorização por coluna nas réplicas válidas
){
  set.seed(seed)
  data <- as.data.frame(data)
  if (!(cluster_var %in% names(data)))
    stop("Cluster '", cluster_var, "' ausente em 'data'.")
  
  cl  <- factor(data[[cluster_var]])
  G   <- nlevels(cl)
  idx <- split(seq_len(nrow(data)), cl)
  
  base_vec <- .pack_all_elas(fit_eqspec, at = at, w_source = w_source)
  q <- length(base_vec)
  
  inst_by_eq_ref <- attr(fit_eqspec, "inst_excl_by_eq")
  draws <- matrix(NA_real_, nrow = R, ncol = q)
  fail_reason <- character(R)
  
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = R, style = 3)
  
  for (b in seq_len(R)) {
    # Poisson+1 por cluster (garante >=1 por grupo)
    cnt_g <- stats::rpois(G, lambda = 1) + 1L
    names(cnt_g) <- levels(cl)
    
    idx_rep <- unlist(
      Map(function(ii, c) rep(ii, each = c), idx, as.list(cnt_g[names(idx)])),
      use.names = FALSE
    )
    
    d_b <- data[idx_rep, , drop = FALSE]
    d_b <- .fix_levels_like(d_b, data)
    
    f_b <- try(refit_quaids_try_eqlist(fit_base, d_b, inst_by_eq_ref), silent = TRUE)
    if (inherits(f_b, "try-error")) { fail_reason[b] <- "refit"; next }
    
    vb <- try(.pack_all_elas(f_b, at = at, w_source = w_source), silent = TRUE)
    if (inherits(vb, "try-error") || length(vb) != q || any(!is.finite(vb))) {
      fail_reason[b] <- "pack_all"; next
    }
    
    # filtro de plausibilidade
    if (max(abs(vb)) > abs_cap) {
      fail_reason[b] <- "abs_cap"; next
    }
    
    draws[b, ] <- vb
    if (verbose) utils::setTxtProgressBar(pb, b)
  }
  if (verbose) close(pb)
  
  keep <- rowSums(is.finite(draws)) == q
  if (verbose) {
    n_ok <- sum(keep); n_tot <- R
    message(sprintf("Réplicas válidas: %d/%d (%.1f%%)", n_ok, n_tot, 100*n_ok/n_tot))
    bad <- table(fail_reason[!keep])
    if (length(bad)) { message("Falhas por motivo:"); print(sort(bad, decreasing = TRUE)) }
  }
  
  draws <- draws[keep, , drop = FALSE]
  if (winsor_p > 0 && nrow(draws) > 0) draws <- .winsorize_cols(draws, p = winsor_p)
  
  list(
    draws = draws,
    base  = base_vec,
    # estatísticas baseadas nas réplicas
    se_sd  = apply(draws, 2, stats::sd,  na.rm = TRUE),
    se_mad = 1.4826 * apply(draws, 2, stats::mad, na.rm = TRUE),
    q025   = apply(draws, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    q975   = apply(draws, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}

# -------------------------------
# 3) Formatação das tabelas finais
#    - se_method: "ciwidth" (padrão), "sd" ou "mad"
#    - IC: percentil (q025, q975)
# -------------------------------
format_boot_tables <- function(B, fit_obj, se_method = c("ciwidth","sd","mad")) {
  se_method <- match.arg(se_method)
  stopifnot(is.list(B), !is.null(B$base), !is.null(B$q025), !is.null(B$q975))
  
  EL <- quaids_elasticities(fit_obj, at = "medians", w_source = "observed")
  goods <- rownames(as.matrix(EL$marshallian))
  K <- length(goods)
  lenM <- K*K; lenH <- K*K; lenE <- K
  
  z975 <- stats::qnorm(0.975)
  se_from_ci <- (B$q975 - B$q025) / (2*z975)
  
  pick_se <- switch(se_method,
                    ciwidth = se_from_ci,
                    sd      = B$se_sd,
                    mad     = B$se_mad
  )
  
  mk <- function(off, len, lab) data.frame(
    lab,
    est = B$base[(off+1):(off+len)],
    se  = pick_se[(off+1):(off+len)],
    lwr = B$q025[(off+1):(off+len)],
    upr = B$q975[(off+1):(off+len)],
    row.names = NULL, check.names = FALSE
  )
  
  lab_M <- expand.grid(i = goods, j = goods, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  lab_H <- lab_M
  lab_E <- data.frame(good = goods, stringsAsFactors = FALSE)
  
  list(
    marshallian = mk(0,         lenM, lab_M),
    hicksian    = mk(lenM,      lenH, lab_H),
    expenditure = mk(lenM+lenH, lenE, lab_E)
  )
}

# -------------------------------
# 4) EXECUÇÃO
# -------------------------------
stopifnot(exists("fit_q"), exists("df_iv"))

# 4.1) Ponto de partida para eq-specific
fit_start <-
  if (exists("fit_sw2")) fit_sw2 else
    if (exists("fit_tryK")) fit_tryK else
      if (exists("fit_bumped")) fit_bumped else
        if (exists("fit_fix")) fit_fix else
          fit_q

# 4.2) Força JUST-ID na eq.6 (UF), preserva e evita conforme pedido
fit_just6 <- make_eq_justID2(
  fit_base     = fit_q,
  data         = df_iv,
  fit_eqspec   = fit_start,
  eqn          = "w_despesahat6",
  prefer_keep  = "IV_d_uf_X12",
  avoid        = c("IV_d_uf_X11"),
  verbose      = TRUE
)

cat("\n== RELATÓRIO (UF) APÓS JUST-ID ==\n")
if (exists("post_refit_report_auto", mode="function")) {
  invisible(post_refit_report_auto(fit_just6, df_iv, cluster_var = "uf", R = 1000, level = 0.95))
} else if (exists("post_refit_report", mode="function")) {
  invisible(post_refit_report(fit_just6, df_iv, cluster_var = "uf", R = 1000, level = 0.95))
}

# 4.3) Bootstrap clusterizado por UF (Poisson+1; sem pesos) + filtros robustos
cat("\n== Rodando bootstrap por UF (Poisson+1; sem pesos) ==\n")
B <- boot_quaids_elasticities_byUF(
  fit_base   = fit_q,
  data       = df_iv,
  fit_eqspec = fit_just6,
  cluster_var= "uf",
  R          = 500,      # ajuste conforme CPU
  seed       = 42,
  verbose    = TRUE,
  abs_cap    = 10,        # descarta réplicas com |E| > 20
  winsor_p   = 0.02       # winsoriza caudas de 1%
)

# 4.4) Tabelas finais (SE pelo método do IC por padrão)
elas_bs <- format_boot_tables(B, fit_just6, se_method = "ciwidth")

cat("\n== ELASTICIDADES — MARSHALL (não compensadas) — IC Bootstrap 95% (UF) ==\n")
print(elas_bs$marshallian, row.names = FALSE)

cat("\n== ELASTICIDADES — HICKS (compensadas) — IC Bootstrap 95% (UF) ==\n")
print(elas_bs$hicksian, row.names = FALSE)

cat("\n== ELASTICIDADES — RENDA/DESPESA — IC Bootstrap 95% (UF) ==\n")
print(elas_bs$expenditure, row.names = FALSE)

# (Opcional) salvar CSV
# write.csv(elas_bs$marshallian, "elas_marshall_bsUF.csv", row.names = FALSE)
# write.csv(elas_bs$hicksian,    "elas_hicks_bsUF.csv",    row.names = FALSE)
# write.csv(elas_bs$expenditure, "elas_renda_bsUF.csv",    row.names = FALSE)
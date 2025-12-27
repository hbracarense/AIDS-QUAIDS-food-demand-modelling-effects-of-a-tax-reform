# ==========================================================
# Bootstrap QU(A)IDS — sumarização robusta + ICs (BC/basic/percentile com fallback)
# Espera: B_sq$draws = matriz (R x 78), B_sq$base = vetor length 78
# Saídas: results_status_quo/
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(readr); library(openxlsx)
})
options(dplyr.summarise.inform = FALSE)
B_sq <- B
# ----- Pré-requisitos -----
stopifnot(exists("B_sq"),
          is.matrix(B_sq$draws),
          length(B_sq$base) == 78)

# ----- Parâmetros -----
alpha      <- 0.05
ci_type    <- "basic"   # "bc", "basic" ou "percentile"
winsor_p   <- 0.05   # winsorização bilateral
abs_cap    <- 6     # cap absoluto |x| <= abs_cap
trim_frac  <- 0.02   # trimming extra (0 = off)

dir_out <- "results_status_quo"
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

# Rótulos (fallback)
if (!exists("goods") || length(goods) != 6) goods <- paste0("good", 1:6)

# ==========================
# Helpers
# ==========================

.apply_winsor_cap <- function(x, p = winsor_p, cap = abs_cap){
  x <- x[is.finite(x)]
  if (!is.finite(cap)) cap <- Inf
  if (is.finite(cap)) x <- pmax(pmin(x, cap), -cap)
  if (is.finite(p) && p > 0 && length(x) >= 5){
    qs <- stats::quantile(x, c(p, 1 - p), na.rm = TRUE, type = 7)
    x  <- pmin(pmax(x, qs[1]), qs[2])
  }
  x
}

.after_trim <- function(x, trim = trim_frac){
  x <- x[is.finite(x)]
  if (is.finite(trim) && trim > 0 && length(x) >= 10){
    qs <- stats::quantile(x, c(trim, 1 - trim), na.rm = TRUE, type = 7)
    x  <- x[x >= qs[1] & x <= qs[2]]
  }
  x
}

.diag_mean <- function(v36){
  stopifnot(length(v36) == 36)
  idx <- 1 + (0:5)*7
  mean(v36[idx], na.rm = TRUE)
}

.detect_block_order <- function(base_vec78){
  stopifnot(length(base_vec78) == 78)
  segs <- list(A = base_vec78[1:36], B = base_vec78[37:72], C = base_vec78[73:78])
  
  # bloco de 6 = renda (se houver dúvida, escolhe o mais próximo de 1)
  lens   <- sapply(segs, length)
  cand6  <- names(segs)[lens == 6]
  if (length(cand6) == 1) {
    where6 <- cand6
  } else {
    closeness <- sapply(cand6, function(k) mean(abs(segs[[k]] - 1), na.rm = TRUE))
    where6 <- names(which.min(closeness))
  }
  
  rest36 <- setdiff(names(segs), where6)
  diag_means <- sapply(rest36, function(k) .diag_mean(segs[[k]]))
  marshall_key <- names(which.min(diag_means))
  hicks_key    <- setdiff(rest36, marshall_key)
  
  list(mh = marshall_key, hx = hicks_key, yl = where6)
}

get_draws_from_matrix <- function(B){
  X <- B$draws
  stopifnot(is.matrix(X), length(B$base) == 78)
  ord <- .detect_block_order(B$base)
  pick <- function(tag) switch(tag, A = 1:36, B = 37:72, C = 73:78)
  idx_mh <- pick(ord$mh); idx_hx <- pick(ord$hx); idx_yl <- pick(ord$yl)
  
  # Marshall (36)
  mat_mh <- X[, idx_mh, drop = FALSE]
  colnames(mat_mh) <- paste0("mh_", seq_len(ncol(mat_mh)))
  mh <- tibble::as_tibble(mat_mh, .name_repair = "unique") |>
    dplyr::mutate(.b = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.b, names_to = "col", values_to = "val") |>
    dplyr::mutate(
      col = as.integer(sub("mh_", "", col)),
      i = ((col - 1L) %/% 6L) + 1L,
      j = ((col - 1L) %% 6L) + 1L
    ) |>
    dplyr::select(.b, i, j, val)
  
  # Hicks (36)
  mat_hx <- X[, idx_hx, drop = FALSE]
  colnames(mat_hx) <- paste0("hx_", seq_len(ncol(mat_hx)))
  hx <- tibble::as_tibble(mat_hx, .name_repair = "unique") |>
    dplyr::mutate(.b = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.b, names_to = "col", values_to = "val") |>
    dplyr::mutate(
      col = as.integer(sub("hx_", "", col)),
      i = ((col - 1L) %/% 6L) + 1L,
      j = ((col - 1L) %% 6L) + 1L
    ) |>
    dplyr::select(.b, i, j, val)
  
  # Income (6)
  mat_yl <- X[, idx_yl, drop = FALSE]
  colnames(mat_yl) <- paste0("yl_", seq_len(ncol(mat_yl)))
  yl <- tibble::as_tibble(mat_yl, .name_repair = "unique") |>
    dplyr::mutate(.b = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.b, names_to = "col", values_to = "val") |>
    dplyr::mutate(i = as.integer(sub("yl_", "", col))) |>
    dplyr::select(.b, i, val)
  
  base_mh <- B$base[idx_mh]; base_hx <- B$base[idx_hx]; base_yl <- B$base[idx_yl]
  message(sprintf("Block mapping: MH=%s, HX=%s, YL=%s", ord$mh, ord$hx, ord$yl))
  list(mh = mh, hx = hx, yl = yl,
       base_mh = base_mh, base_hx = base_hx, base_yl = base_yl)
}

# ---- ICs (versões "safe" com base R) ----
.safe_quantile_pair <- function(v, a = alpha){
  v <- v[is.finite(v)]
  if (length(v) < 2) return(c(NA_real_, NA_real_))
  q <- stats::quantile(v, c(a/2, 1 - a/2), na.rm = TRUE, type = 7)
  c(q[1], q[2])
}

.safe_basic_pair <- function(v, theta_hat, a = alpha){
  v <- v[is.finite(v)]
  if (length(v) < 2 || !is.finite(theta_hat)) return(c(NA_real_, NA_real_))
  q <- stats::quantile(v, c(a/2, 1 - a/2), na.rm = TRUE, type = 7)
  c(2*theta_hat - q[2], 2*theta_hat - q[1])
}

.safe_bc_pair <- function(v, theta_hat, a = alpha){
  v <- v[is.finite(v)]
  n <- length(v); if (n < 10 || !is.finite(theta_hat)) return(c(NA_real_, NA_real_))
  eps <- 1/(n + 1)
  p_lt <- mean(v < theta_hat, na.rm = TRUE); p_lt <- min(max(p_lt, eps), 1 - eps)
  z0  <- qnorm(p_lt); zal <- qnorm(a/2); zau <- qnorm(1 - a/2)
  a1  <- pnorm(2*z0 + zal); a2 <- pnorm(2*z0 + zau)
  q   <- stats::quantile(v, c(a1, a2), na.rm = TRUE, type = 7)
  c(q[1], q[2])
}

.choose_ci_with_fallback <- function(pref_pair, pct_pair){
  bad <- any(!is.finite(pref_pair)) || (pref_pair[2] < pref_pair[1])
  if (bad) list(lo = pct_pair[1], hi = pct_pair[2], method = "percentile(fallback)")
  else     list(lo = pref_pair[1], hi = pref_pair[2], method = "preferred")
}

summarise_block <- function(df_long, base_vec, a = alpha, ci = ci_type){
  stopifnot(all(c(".b","i","j","val") %in% names(df_long)), length(base_vec) == 36)
  to_pos <- function(i,j) (i - 1L) * 6L + j
  
  out <- df_long |>
    dplyr::group_by(i, j) |>
    dplyr::summarise(
      draws   = list(.after_trim(.apply_winsor_cap(val))),
      n_kept  = length(draws[[1]]),
      est     = stats::median(draws[[1]], na.rm = TRUE),
      se_boot = stats::sd(draws[[1]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(pos = to_pos(i, j),
                  theta_hat = base_vec[pos])
  
  # Percentil
  q_pct <- lapply(out$draws, .safe_quantile_pair, a = a)
  out$lo_pct <- vapply(q_pct, function(z) z[1], numeric(1))
  out$hi_pct <- vapply(q_pct, function(z) z[2], numeric(1))
  
  # Basic
  q_bas <- Map(function(v, th) .safe_basic_pair(v, th, a = a), out$draws, out$theta_hat)
  out$lo_basic <- vapply(q_bas, function(z) z[1], numeric(1))
  out$hi_basic <- vapply(q_bas, function(z) z[2], numeric(1))
  
  # BC
  q_bc <- Map(function(v, th) .safe_bc_pair(v, th, a = a), out$draws, out$theta_hat)
  out$lo_bc <- vapply(q_bc, function(z) z[1], numeric(1))
  out$hi_bc <- vapply(q_bc, function(z) z[2], numeric(1))
  
  # Escolha + fallback (sem SIMPLIFY aqui!)
  chosen <- Map(function(pref, pct){
    .choose_ci_with_fallback(pref_pair = pref, pct_pair = pct)
  },
  pref = switch(ci,
                "bc"         = Map(function(a,b) c(a,b), out$lo_bc,   out$hi_bc),
                "basic"      = Map(function(a,b) c(a,b), out$lo_basic,out$hi_basic),
                "percentile" = Map(function(a,b) c(a,b), out$lo_pct,  out$hi_pct)),
  pct  = Map(function(a,b) c(a,b), out$lo_pct, out$hi_pct))
  
  out$lo      <- vapply(chosen, function(x) x$lo,     numeric(1))
  out$hi      <- vapply(chosen, function(x) x$hi,     numeric(1))
  out$ci_used <- vapply(chosen, function(x) x$method, character(1))
  
  out
}

summarise_income <- function(df_long, base_vec, a = alpha, ci = ci_type){
  stopifnot(all(c(".b","i","val") %in% names(df_long)), length(base_vec) == 6)
  
  out <- df_long |>
    dplyr::group_by(i) |>
    dplyr::summarise(
      draws   = list(.after_trim(.apply_winsor_cap(val))),
      n_kept  = length(draws[[1]]),
      eta     = stats::median(draws[[1]], na.rm = TRUE),
      se_boot = stats::sd(draws[[1]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(i) |>
    dplyr::mutate(theta_hat = base_vec)
  
  # Percentil
  q_pct <- lapply(out$draws, .safe_quantile_pair, a = a)
  out$lo_pct <- vapply(q_pct, function(z) z[1], numeric(1))
  out$hi_pct <- vapply(q_pct, function(z) z[2], numeric(1))
  
  # Basic
  q_bas <- Map(function(v, th) .safe_basic_pair(v, th, a = a), out$draws, out$theta_hat)
  out$lo_basic <- vapply(q_bas, function(z) z[1], numeric(1))
  out$hi_basic <- vapply(q_bas, function(z) z[2], numeric(1))
  
  # BC
  q_bc <- Map(function(v, th) .safe_bc_pair(v, th, a = a), out$draws, out$theta_hat)
  out$lo_bc <- vapply(q_bc, function(z) z[1], numeric(1))
  out$hi_bc <- vapply(q_bc, function(z) z[2], numeric(1))
  
  # Escolha + fallback (sem SIMPLIFY aqui!)
  chosen <- Map(function(pref, pct){
    .choose_ci_with_fallback(pref_pair = pref, pct_pair = pct)
  },
  pref = switch(ci,
                "bc"         = Map(function(a,b) c(a,b), out$lo_bc,   out$hi_bc),
                "basic"      = Map(function(a,b) c(a,b), out$lo_basic,out$hi_basic),
                "percentile" = Map(function(a,b) c(a,b), out$lo_pct,  out$hi_pct)),
  pct  = Map(function(a,b) c(a,b), out$lo_pct, out$hi_pct))
  
  out$lo      <- vapply(chosen, function(x) x$lo,     numeric(1))
  out$hi      <- vapply(chosen, function(x) x$hi,     numeric(1))
  out$ci_used <- vapply(chosen, function(x) x$method, character(1))
  
  out
}

# ============================
# 1) Draws longos e sumarização
# ============================
DW <- get_draws_from_matrix(B_sq)

mh_sum <- summarise_block(DW$mh, base_vec = DW$base_mh, a = alpha, ci = ci_type)
hx_sum <- summarise_block(DW$hx, base_vec = DW$base_hx, a = alpha, ci = ci_type)
yl_sum <- summarise_income(DW$yl, base_vec = DW$base_yl, a = alpha, ci = ci_type)

# ============================
# 2) Tabelas finais e rótulos
# ============================
marshall_full <- mh_sum |>
  dplyr::mutate(
    good_i = factor(i, levels = 1:6, labels = goods),
    good_j = factor(j, levels = 1:6, labels = goods),
    z = dplyr::if_else(is.finite(se_boot) & se_boot > 0, est / se_boot, NA_real_)
  ) |>
  dplyr::select(i, j, est, se_boot, lo, hi,
                lo_pct, hi_pct, lo_bc, hi_bc, lo_basic, hi_basic,
                ci_used, n_kept, good_i, good_j, z)

hicks_full <- hx_sum |>
  dplyr::mutate(
    good_i = factor(i, levels = 1:6, labels = goods),
    good_j = factor(j, levels = 1:6, labels = goods),
    z = dplyr::if_else(is.finite(se_boot) & se_boot > 0, est / se_boot, NA_real_)
  ) |>
  dplyr::select(i, j, est, se_boot, lo, hi,
                lo_pct, hi_pct, lo_bc, hi_bc, lo_basic, hi_basic,
                ci_used, n_kept, good_i, good_j, z)

eta_out <- yl_sum |>
  dplyr::mutate(good = factor(i, levels = 1:6, labels = goods)) |>
  dplyr::select(i, good, eta, se_boot, lo, hi,
                lo_pct, hi_pct, lo_bc, hi_bc, lo_basic, hi_basic,
                ci_used, n_kept, theta_hat)

marshall_diag <- marshall_full |> dplyr::filter(i == j)
hicks_diag    <- hicks_full    |> dplyr::filter(i == j)

# ============================
# 3) Heatmap (Marshall z-score)
# ============================
p_heat <- ggplot(marshall_full, aes(x = good_j, y = good_i, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(name = "z-score", low = "#2166AC", mid = "white",
                       high = "#B2182B", midpoint = 0) +
  labs(title = "Standardized Marshallian Price Elasticities (status quo)",
       x = "Price of j", y = "Quantity of i") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(dir_out, "heatmap_status_quo_marshall_z.png"),
       p_heat, width = 7.5, height = 6.2, dpi = 300)

# ============================
# 4) Exporta CSVs
# ============================
#readr::write_csv(marshall_full, file.path(dir_out, "marshall_full_status_quo.csv"))
#readr::write_csv(hicks_full,    file.path(dir_out, "hicks_full_status_quo.csv"))
#readr::write_csv(marshall_diag, file.path(dir_out, "marshall_diag_status_quo.csv"))
#readr::write_csv(hicks_diag,    file.path(dir_out, "hicks_diag_status_quo.csv"))
#readr::write_csv(eta_out,       file.path(dir_out, "eta_income_status_quo.csv"))

# ============================
# 5) Prints
# ============================
message("\nMarshall (diag, bootstrap median point, CI pref = ", ci_type, "):")
print(marshall_diag)

message("\nHicks (diag, bootstrap median point, CI pref = ", ci_type, "):")
print(hicks_diag)

message("\nIncome elasticities (bootstrap median point, CI pref = ", ci_type, "):")
print(eta_out)
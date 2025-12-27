## ====== Pacotes ======
suppressPackageStartupMessages({
  library(AER)
  library(sandwich)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

## ====== Helpers ======

.align_cluster <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)
  mf  <- stats::model.frame(model)
  rid <- rownames(mf)
  if (!is.null(names(cluster))) {
    return(cluster[rid])
  } else {
    idx <- suppressWarnings(as.integer(rid))
    if (all(is.finite(idx))) return(cluster[idx])
  }
  cluster
}

.get_XZ_u <- function(ivm) {
  X <- try(model.matrix(ivm, "regressors"),  silent = TRUE)
  Z <- try(model.matrix(ivm, "instruments"), silent = TRUE)
  if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NULL)
  list(X = X, Z = Z, u = residuals(ivm))
}

.j_overid_manual <- function(ivm, cluster = NULL, ridge = 1e-8) {
  mats <- .get_XZ_u(ivm); if (is.null(mats)) return(c(J = NA_real_, df = NA_real_, p = NA_real_))
  X <- mats$X; Z <- mats$Z; u <- as.numeric(mats$u)
  n <- NROW(Z); kz <- ncol(Z); kx <- ncol(X); df <- kz - kx
  if (!is.finite(df) || df <= 0) return(c(J = NA_real_, df = df, p = NA_real_))
  cl <- .align_cluster(ivm, cluster)
  Zu <- Z * u
  if (!is.null(cl)) {
    G  <- rowsum(Zu, group = as.factor(cl))
    S  <- crossprod(as.matrix(G)) / n
  } else {
    S  <- crossprod(Zu) / n
  }
  S <- (S + t(S))/2
  d <- mean(diag(S)); if (!is.finite(d) || d <= 0) d <- 1
  S_r <- S + diag(ridge * d, ncol(S))
  gbar <- colMeans(Zu)
  J    <- as.numeric(n * crossprod(gbar, solve(S_r, gbar)))
  p    <- stats::pchisq(J, df = df, lower.tail = FALSE)
  c(J = J, df = df, p = p)
}

.clean_terms <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(x)]
  setdiff(x, c("1","-1","0","+","~","|"))
}

.build_inst_terms <- function(fit_obj) {
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  unique(c(exogs, inst))
}

.get_excluded_IVs <- function(fit_obj) {
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  setdiff(inst, exogs)
}

.fit_iv_safe <- function(f_yx, inst_terms, data) {
  try(ivreg(f_yx, instruments = reformulate(inst_terms), data = data), silent = TRUE)
}

## ====== 1) J por equação (robusto/cluster) ======
jtests_by_eq_manual <- function(fit_obj, data, cluster_var = NULL, verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_terms <- .build_inst_terms(fit_obj)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  rows <- lapply(seq_along(fit_obj$eq), function(i){
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    ivm  <- .fit_iv_safe(f_yx, inst_terms, data)
    if (inherits(ivm, "try-error")) {
      if (verbose) message("ivreg falhou na eq ", eqn, ": ", attr(ivm, "condition")$message)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    mats <- .get_XZ_u(ivm); if (is.null(mats)) {
      if (verbose) message("Matrizes X/Z indisponíveis para ", eqn)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    J_h <- try(AER::sargan(ivm), silent = TRUE)
    J_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$statistic)) else NA_real_
    p_s <- if (!inherits(J_h, "try-error")) unname(as.numeric(J_h$p.value))   else NA_real_
    Jrob <- .j_overid_manual(ivm, cluster = cl)
    data.frame(
      eq = eqn, rank_X = rX, rank_Z = rZ, df_J = dfJ,
      J_sargan = J_s, p_sargan = p_s,
      J_hansen = unname(Jrob["J"]), p_hansen = unname(Jrob["p"]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

## ====== 2) Blocos de IVs ======
make_blocks_explicit <- function(fit_obj) {
  excl <- .get_excluded_IVs(fit_obj)
  if (!length(excl)) return(list())
  base <- gsub("\\d+$", "", excl)
  base <- gsub("\\s+", "", base)
  sp   <- split(excl, base)
  names(sp) <- make.unique(names(sp))
  sp
}

## ====== 3) C-tests ADAPTATIVOS (diferença de Hansen) ======
# Remove subconjuntos de tamanho k <= max_k_drop e só calcula C quando o modelo restrito fica sobre-ID.
c_tests_blocks_adaptive <- function(fit_obj, data, blocks,
                                    cluster_var = NULL,
                                    max_k_drop = 2,  # evita explosão combinatória
                                    verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_base <- .build_inst_terms(fit_obj)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  out <- list()
  
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    
    iv_full <- .fit_iv_safe(f_yx, inst_base, data)
    if (inherits(iv_full, "try-error")) {
      if (verbose) message("ivreg (full) falhou em ", eqn)
      next
    }
    Jf <- .j_overid_manual(iv_full, cluster = cl)
    full_mats <- .get_XZ_u(iv_full)
    if (is.null(full_mats)) next
    df_full <- ncol(full_mats$Z) - ncol(full_mats$X)
    if (!is.finite(df_full) || df_full <= 0) {
      if (verbose) message("Eq ", eqn, ": full não sobre-ID (df_J=", df_full, "). Pulando C-tests.")
      next
    }
    k_cap <- max(1, min(max_k_drop, df_full - 1))
    
    for (bk in names(blocks)) {
      in_block <- intersect(inst_base, blocks[[bk]])
      if (!length(in_block)) next
      # gera todas as combinações de 1..k_cap instrumentos a remover dentro do bloco
      for (k in 1:k_cap) {
        if (length(in_block) < k) next
        cmb <- utils::combn(in_block, k, simplify = FALSE)
        for (drop_set in cmb) {
          inst_restr <- setdiff(inst_base, drop_set)
          iv_r <- .fit_iv_safe(f_yx, inst_restr, data)
          if (inherits(iv_r, "try-error")) next
          mats_r <- .get_XZ_u(iv_r); if (is.null(mats_r)) next
          df_restr <- ncol(mats_r$Z) - ncol(mats_r$X)
          if (!is.finite(df_restr) || df_restr <= 0) {
            # restrito não sobre-ID: não dá para calcular J_restr (e logo C)
            next
          }
          Jr <- .j_overid_manual(iv_r, cluster = cl)
          C  <- as.numeric(Jr["J"] - Jf["J"])
          dfC <- length(drop_set)
          pC <- stats::pchisq(C, df = dfC, lower.tail = FALSE)
          out[[length(out) + 1L]] <- data.frame(
            eq = eqn,
            block = bk,
            drop_set = paste(drop_set, collapse = " + "),
            k_removed = dfC,
            J_full = unname(Jf["J"]), dfJ_full = unname(Jf["df"]),
            J_restr = unname(Jr["J"]), dfJ_restr = unname(Jr["df"]),
            C_stat = C, df_C = dfC, p_C = pC,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(out)) {
    return(data.frame(
      eq=character(), block=character(), drop_set=character(), k_removed=integer(),
      J_full=double(), dfJ_full=integer(), J_restr=double(), dfJ_restr=integer(),
      C_stat=double(), df_C=integer(), p_C=double()
    ))
  }
  dplyr::bind_rows(out)
}

## ====== 4) Poda e nova lista de IVs ======
suggest_prune <- function(ctab, alpha = 0.05) {
  if (!NROW(ctab)) return(character(0))
  # retorna os conjuntos específicos (drop_set) com evidência forte
  bad <- subset(ctab, is.finite(p_C) & p_C < alpha)
  unique(bad$drop_set)
}

new_inst_after_prune <- function(fit_obj, to_drop_sets, blocks) {
  inst_base <- .build_inst_terms(fit_obj)
  rhs_all   <- .clean_terms(fit_obj$rhs_terms)
  endo_ln   <- grep("^ln_", rhs_all, value = TRUE)
  exogs     <- setdiff(rhs_all, endo_ln)
  drop_vec  <- character(0)
  if (length(to_drop_sets)) {
    drop_vec <- unique(unlist(strsplit(to_drop_sets, "\\s\\+\\s")))
  }
  excl_ivs  <- setdiff(inst_base, exogs)
  excl_new  <- setdiff(excl_ivs, drop_vec)
  unique(c(exogs, excl_new))
}

## ====== Execução/Exemplo ======

rhs_all  <- fit_q$rhs_terms
inst_all <- fit_q$inst_terms
cat("Endógenos no RHS (esperado ficarem fora de Z): ",
    paste(intersect(rhs_all, grep("^ln_", rhs_all, value=TRUE)), collapse=", "), "\n", sep="")

## 1) J por equação (cluster se existir)
jtab <- jtests_by_eq_manual(
  fit_q, df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"
)
print(jtab)

## 2) C-tests adaptativos (remove 1 e 2 IVs por bloco, pois df_J=3)
blocks <- make_blocks_explicit(fit_q)
message("Blocos detectados: ",
        if (length(blocks)) paste0(names(blocks), " (", vapply(blocks, length, integer(1)), " IVs)", collapse=" | ")
        else "nenhum (sem IVs excluídos)")

ctab <- c_tests_blocks_adaptive(
  fit_q, df_iv, blocks = blocks,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf",
  max_k_drop = 2
)

if (nrow(ctab)) {
  ctab <- dplyr::arrange(ctab, eq, p_C)
  print(ctab, digits = 4)
  cat("\nTOP hits por equação (menor p_C):\n")
  print(ctab %>% group_by(eq) %>% slice_min(order_by = p_C, n = 3, with_ties = FALSE) %>% ungroup())
} else {
  message("Nenhum C-test executável (provável df_J pequeno ou blocos vazios).")
}

## 3) Sugerir poda e nova lista de IVs
to_drop_sets <- suggest_prune(ctab, alpha = 0.05)
cat("\nSubconjuntos a considerar poda (p_C < 0.05): ",
    if (length(to_drop_sets)) paste(to_drop_sets, collapse=" || ") else "(nenhum)", "\n", sep="")

inst_new <- new_inst_after_prune(fit_q, to_drop_sets, blocks)
cat("\nInst terms (NOVOS) após poda:\n", paste(inst_new, collapse = " + "), "\n")

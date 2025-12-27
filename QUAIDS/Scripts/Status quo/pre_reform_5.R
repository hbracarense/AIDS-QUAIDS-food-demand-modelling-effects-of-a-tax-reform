## ========= Requisitos =========
suppressPackageStartupMessages({
  library(AER); library(sandwich); library(dplyr); library(openxlsx)
})

## ========= Helpers robustos =========
.mm <- function(terms, data, add_intercept = TRUE){
  if (length(terms) == 0L) return(matrix(1, nrow(data), 1))
  mm <- model.matrix(reformulate(terms, intercept = add_intercept), data = data)
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}
rank_dfJ <- function(X, Z) {
  rX <- qr(X)$rank; rZ <- qr(Z)$rank
  c(rank_X = rX, rank_Z = rZ, dfJ = rZ - rX)
}
J_hansen <- function(ivmod, cluster = NULL){
  vc <- if (is.null(cluster)) sandwich::vcovHC(ivmod, type = "HC3")
  else sandwich::vcovCL(ivmod, cluster = cluster, type = "HC3")
  s  <- try(AER::sargan(ivmod, vcov = vc), silent = TRUE)
  if (inherits(s, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  c(stat = unname(as.numeric(s$statistic)),
    df   = unname(as.numeric(s$parameter)),
    p    = unname(as.numeric(s$p.value)))
}

get_iv_sets <- function(fit){
  inst_all <- unique(fit$inst_terms)
  rhs_all  <- unique(fit$rhs_terms)
  controls <- intersect(inst_all, rhs_all)     # exógenos incluídos no RHS (instrumentos "incluídos")
  iv_excl  <- setdiff(inst_all, controls)      # IVs excluídos (geram sobre-ID)
  list(inst_all = inst_all, controls = controls, iv_excl = iv_excl)
}

## ========= Blocos explícitos (ajuste se quiser) =========
make_blocks_explicit <- function(fit){
  S <- get_iv_sets(fit)
  iv <- S$iv_excl
  blocks <- list(
    op_shift = grep("^iv_op",  iv, value = TRUE),      # ex.: iv_op01:iv_op03
    uf_dums  = grep("^IV_d_uf_", iv, value = TRUE)     # ex.: IV_d_uf_X11:...X15
    # você pode adicionar outros grupos aqui, p.ex.:
    # other    = setdiff(iv, c(grep("^iv_op", iv, value=TRUE),
    #                          grep("^IV_d_uf_", iv, value=TRUE)))
  )
  # remove blocos vazios
  blocks[vapply(blocks, length, FUN.VALUE = 1L) > 0]
}

## ========= C-tests por blocos (diferença de Hansen) =========
c_tests_blocks2 <- function(fit, data, blocks, cluster_var = NULL){
  stopifnot(inherits(fit, "quaids_km1_fit"))
  sets <- get_iv_sets(fit)
  inst_all <- sets$inst_all
  eqs <- fit$eq
  cl  <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  out <- vector("list", 0L)
  
  for (i in seq_along(eqs)) {
    f_yx   <- formula(eqs[[i]])
    eqname <- as.character(f_yx[[2]])
    
    X     <- .mm(delete.response(terms(f_yx)), data, add_intercept = TRUE)
    Zfull <- .mm(inst_all, data, add_intercept = TRUE)
    rkF   <- rank_dfJ(X, Zfull)
    
    # precisa estar sobre-identificado no conjunto cheio
    if (rkF["dfJ"] <= 0) {
      out[[length(out) + 1L]] <- data.frame(
        eq = eqname, block = NA_character_, k_removed = NA_integer_,
        J_full = NA_real_, dfJ_full = rkF["dfJ"],
        J_restr = NA_real_, dfJ_restr = NA_real_,
        C_stat = NA_real_, df_C = NA_real_, p_C = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    iv_full <- AER::ivreg(f_yx, instruments = reformulate(inst_all), data = data)
    JF      <- J_hansen(iv_full, cluster = cl)  # c(stat, df, p)
    
    for (bn in names(blocks)) {
      inst_restr <- setdiff(inst_all, blocks[[bn]])
      Zr  <- .mm(inst_restr, data, add_intercept = TRUE)
      rkR <- rank_dfJ(X, Zr)
      
      Jr <- c(stat = NA_real_, df = NA_real_, p = NA_real_)
      if (rkR["dfJ"] > 0) {
        iv_r <- AER::ivreg(f_yx, instruments = reformulate(inst_restr), data = data)
        Jr   <- J_hansen(iv_r, cluster = cl)
      }
      
      Cstat <- if (is.finite(JF["stat"]) && is.finite(Jr["stat"])) JF["stat"] - Jr["stat"] else NA_real_
      ddf   <- if (is.finite(JF["df"])  && is.finite(Jr["df"]))  JF["df"]  - Jr["df"]  else NA_real_
      pC    <- if (is.finite(Cstat) && is.finite(ddf) && ddf > 0) pchisq(Cstat, df = ddf, lower.tail = FALSE) else NA_real_
      
      out[[length(out) + 1L]] <- data.frame(
        eq = eqname,
        block = bn,
        k_removed = length(blocks[[bn]]),
        J_full = unname(JF["stat"]), dfJ_full = unname(JF["df"]),
        J_restr = unname(Jr["stat"]), dfJ_restr = unname(Jr["df"]),
        C_stat = Cstat, df_C = ddf, p_C = pC,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(out)
}

## ========= Sugerir poda e produzir novos inst_terms =========
suggest_prune <- function(ctab, alpha = 0.05){
  if (!nrow(ctab)) return(character(0))
  bad <- subset(ctab, is.finite(p_C) & p_C < alpha)
  if (!nrow(bad)) return(character(0))
  # ordena blocos pelos p_C mais baixos (agregando sobre equações)
  ord <- bad %>%
    dplyr::group_by(block) %>%
    dplyr::summarise(min_p = min(p_C, na.rm = TRUE), n_eq = dplyr::n()) %>%
    dplyr::arrange(min_p)
  ord$block
}

new_inst_after_prune <- function(fit, to_drop_blocks, blocks){
  inst_all <- unique(fit$inst_terms)
  drop_iv  <- unique(unlist(blocks[to_drop_blocks], use.names = FALSE))
  setdiff(inst_all, drop_iv)
}

## ========= COMO USAR =========

# 1) J por equação (robusto com cluster se existir)
suppressPackageStartupMessages({
  library(AER); library(sandwich); library(dplyr)
})

## --- helpers robustos ---
.mm <- function(formula_or_terms, data, add_intercept = TRUE){
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    if (length(formula_or_terms) == 0L) {
      return(matrix(1, nrow(data), 1,
                    dimnames = list(NULL, "(Intercept)")))  # fallback
    }
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}

.clean_terms <- function(v, data){
  v <- unique(v)
  v <- v[!is.na(v) & nzchar(v)]
  v <- setdiff(v, c("1","-1","0","+","~","|"))
  v[v %in% names(data)]
}

.J_hansen <- function(ivmod, cluster = NULL){
  if (!is.null(cluster)) {
    mf  <- stats::model.frame(ivmod)
    rid <- rownames(mf)
    if (!is.null(names(cluster))) {
      cluster <- cluster[rid]
    } else {
      idx <- suppressWarnings(as.integer(rid))
      if (all(is.finite(idx))) cluster <- cluster[idx]
    }
  }
  
  if (!is.null(cluster)) {
    # Requer: library(clubSandwich)
    vc <- clubSandwich::vcovCR(ivmod, cluster = cluster, type = "CR2")
  } else {
    vc <- sandwich::vcovHC(ivmod, type = "HC1")
  }
  
  s <- try(AER::sargan(ivmod, vcov = vc), silent = TRUE)
  if (inherits(s, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  c(stat = unname(as.numeric(s$statistic)),
    df   = unname(as.numeric(s$parameter)),
    p    = unname(as.numeric(s$p.value)))
}

## --- J por equação (robusto com cluster se existir) ---
suppressPackageStartupMessages({ library(AER); library(sandwich); library(dplyr) })

.mm <- function(formula_or_terms, data, add_intercept = TRUE){
  if (inherits(formula_or_terms, "formula")) {
    mm <- model.matrix(formula_or_terms, data = data)
  } else {
    # usa fórmula, deixando model.matrix expandir fatores/interações
    mm <- model.matrix(reformulate(formula_or_terms, intercept = add_intercept), data = data)
  }
  ok <- apply(mm, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  mm <- mm[, ok, drop = FALSE]
  q  <- qr(mm); if (q$rank < ncol(mm)) mm <- mm[, q$pivot[seq_len(q$rank)], drop = FALSE]
  mm
}

.J_hansen <- function(ivmod, cluster = NULL){
  vc <- if (is.null(cluster)) sandwich::vcovHC(ivmod, type = "HC3")
  else sandwich::vcovCL(ivmod, cluster = cluster, type = "HC3")
  s  <- try(AER::sargan(ivmod, vcov = vc), silent = TRUE)
  if (inherits(s, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  c(stat = unname(as.numeric(s$statistic)),
    df   = unname(as.numeric(s$parameter)),
    p    = unname(as.numeric(s$p.value)))
}

## ---------- 1) vcov robusto p/ Hansen-J que funciona com ivreg ----------

# alinha o cluster ao modelo (evita comprimentos diferentes com NA handling)
.align_cluster <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)
  mf  <- stats::model.frame(model)
  rid <- rownames(mf)
  if (!is.null(names(cluster))) {
    return(cluster[rid])                    # por nomes
  } else {
    idx <- suppressWarnings(as.integer(rid))
    if (all(is.finite(idx))) return(cluster[idx])
  }
  cluster
}

# Preferir CR2 (clubSandwich) com cluster; senão HC1/CL sem 'type="HC3"'
.J_hansen <- function(ivmod, cluster = NULL){
  cl_aligned <- .align_cluster(ivmod, cluster)
  
  vc <- try({
    if (!is.null(cl_aligned) && requireNamespace("clubSandwich", quietly = TRUE)) {
      clubSandwich::vcovCR(ivmod, cluster = cl_aligned, type = "CR2")
    } else if (!is.null(cl_aligned)) {
      # NADA de type="HC3" aqui: para ivreg o 'type' de vcovCL tem outro significado
      sandwich::vcovCL(ivmod, cluster = cl_aligned)
    } else {
      sandwich::vcovHC(ivmod, type = "HC1")
    }
  }, silent = TRUE)
  
  if (inherits(vc, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  s <- try(AER::sargan(ivmod, vcov = vc), silent = TRUE)
  if (inherits(s, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  
  c(stat = unname(as.numeric(s$statistic)),
    df   = unname(as.numeric(s$parameter)),
    p    = unname(as.numeric(s$p.value)))
}

## ---------- 2) J por equação (robusto/cluster), garantindo Z completo ----------

## ====== Pacotes ======
suppressPackageStartupMessages({
  library(AER)
  library(sandwich)
  library(dplyr)
})

## ====== Helpers ======

# Alinha o vetor de cluster ao modelo (linhas após na.omit etc.)
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

# Matriz de instrumentos e regressores do ivreg (alinhadas ao modelo)
.get_XZ_u <- function(ivm) {
  X <- try(model.matrix(ivm, "regressors"),  silent = TRUE)
  Z <- try(model.matrix(ivm, "instruments"), silent = TRUE)
  if (inherits(X, "try-error") || inherits(Z, "try-error")) return(NULL)
  list(X = X, Z = Z, u = residuals(ivm))
}

# J robusto (Arellano). Se cluster=NULL, vira HC0.
# S regularizado levemente para evitar singularidade.
.j_overid_manual <- function(ivm, cluster = NULL, ridge = 1e-8) {
  mats <- .get_XZ_u(ivm); if (is.null(mats)) return(c(J = NA_real_, df = NA_real_, p = NA_real_))
  X <- mats$X; Z <- mats$Z; u <- as.numeric(mats$u)
  n <- NROW(Z); kz <- ncol(Z); kx <- ncol(X); df <- kz - kx
  if (!is.finite(df) || df <= 0) return(c(J = NA_real_, df = df, p = NA_real_))
  
  # alinhar cluster às observações usadas no modelo
  cl <- .align_cluster(ivm, cluster)
  
  Zu <- Z * u  # n x kz
  if (!is.null(cl)) {
    G  <- rowsum(Zu, group = as.factor(cl))   # soma por cluster
    S  <- crossprod(as.matrix(G)) / n         # (kz x kz)
  } else {
    S  <- crossprod(Zu) / n                   # HC0
  }
  S <- (S + t(S)) / 2
  d <- mean(diag(S)); if (!is.finite(d) || d <= 0) d <- 1
  S_r <- S + diag(ridge * d, ncol(S))
  
  gbar <- colMeans(Zu)                        # (1/n) Z'u
  J    <- as.numeric(n * crossprod(gbar, solve(S_r, gbar)))
  p    <- stats::pchisq(J, df = df, lower.tail = FALSE)
  c(J = J, df = df, p = p)
}

# Limpa termos textuais triviais
.clean_terms <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(x)]
  setdiff(x, c("1","-1","0","+","~","|"))
}

# Constrói Z final = exógenos do RHS + IVs do fit (sem lixo de símbolos)
.build_inst_terms <- function(fit_obj) {
  rhs_all <- fit_obj$rhs_terms
  rhs_all <- .clean_terms(rhs_all)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  unique(c(exogs, inst))  # garante exógenos do RHS em Z
}

# Retorna apenas os IVs excluídos (em termos de rótulos de fórmula)
.get_excluded_IVs <- function(fit_obj) {
  rhs_all <- .clean_terms(fit_obj$rhs_terms)
  endo_ln <- grep("^ln_", rhs_all, value = TRUE)
  exogs   <- setdiff(rhs_all, endo_ln)
  inst    <- .clean_terms(fit_obj$inst_terms)
  setdiff(inst, exogs)
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
    
    # estima com Z completo (exógenos + IVs excluídos)
    ivm <- try(ivreg(f_yx, instruments = reformulate(inst_terms), data = data), silent = TRUE)
    if (inherits(ivm, "try-error")) {
      if (verbose) message("ivreg falhou na eq ", eqn, ": ", attr(ivm, "condition")$message)
      return(data.frame(eq = eqn, rank_X = NA_integer_, rank_Z = NA_integer_, df_J = NA_integer_,
                        J_sargan = NA_real_, p_sargan = NA_real_,
                        J_hansen = NA_real_, p_hansen = NA_real_))
    }
    
    mats <- .get_XZ_u(ivm); if (is.null(mats)) {
      if (verbose) message("Matrizes X/Z não disponíveis para ", eqn)
      return(data.frame(eq = eqn, rank_X = NA, rank_Z = NA, df_J = NA,
                        J_sargan = NA, p_sargan = NA, J_hansen = NA, p_hansen = NA))
    }
    rX <- qr(mats$X)$rank; rZ <- qr(mats$Z)$rank; dfJ <- rZ - rX
    
    # J homocedástico (Sargan) e robusto/cluster (Hansen)
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

## ====== 2) C-tests por blocos (diferença de Hansen) ======

# Cria blocos de IVs excluídos por "raiz" do nome (remove dígitos finais),
# mantendo apenas termos que NÃO estão no RHS como exógenos.
make_blocks_explicit <- function(fit_obj) {
  excl <- .get_excluded_IVs(fit_obj)
  if (!length(excl)) return(list())
  base <- gsub("\\d+$", "", excl)      # remove sufixos numéricos
  base <- gsub("\\s+", "", base)       # tira espaços acidentais
  sp <- split(excl, base)
  # mantém nomes limpos
  names(sp) <- make.unique(names(sp))
  sp
}

# 'blocks' deve ser uma lista nomeada: ex. list(op = c("iv_op01","iv_op02"), duf = c("IV_d_uf_X11", ...))
c_tests_blocks_manual <- function(fit_obj, data, blocks, cluster_var = NULL, verbose = TRUE) {
  stopifnot(inherits(fit_obj, "quaids_km1_fit"))
  data <- as.data.frame(data)
  inst_base <- .build_inst_terms(fit_obj)
  cl <- if (!is.null(cluster_var) && cluster_var %in% names(data)) data[[cluster_var]] else NULL
  
  out <- list()
  for (i in seq_along(fit_obj$eq)) {
    f_yx <- formula(fit_obj$eq[[i]])
    eqn  <- as.character(f_yx[[2]])
    
    # Full model
    iv_full <- try(ivreg(f_yx, instruments = reformulate(inst_base), data = data), silent = TRUE)
    if (inherits(iv_full, "try-error")) {
      if (verbose) message("ivreg (full) falhou em ", eqn)
      next
    }
    Jf <- .j_overid_manual(iv_full, cluster = cl)
    kz_full <- length(colnames(model.matrix(iv_full, "instruments")))
    kx_full <- length(colnames(model.matrix(iv_full, "regressors")))
    df_full <- as.integer(kz_full - kx_full)
    if (!is.finite(df_full) || df_full <= 0) {
      if (verbose) message("Eq ", eqn, " não sobre-ID no full; pulando C-tests.")
      next
    }
    
    for (bk in names(blocks)) {
      drop_set <- intersect(inst_base, blocks[[bk]])
      if (!length(drop_set)) next
      inst_restr <- setdiff(inst_base, drop_set)
      
      iv_r <- try(ivreg(f_yx, instruments = reformulate(inst_restr), data = data), silent = TRUE)
      if (inherits(iv_r, "try-error")) {
        if (verbose) message("ivreg (restr) falhou em ", eqn, " bloco ", bk)
        next
      }
      Jr <- .j_overid_manual(iv_r, cluster = cl)
      
      # Diferença de Hansen: C = J_restr - J_full, df = |drop_set|
      C   <- as.numeric(Jr["J"] - Jf["J"])
      dfC <- length(drop_set)
      pC  <- stats::pchisq(C, df = dfC, lower.tail = FALSE)
      
      out[[length(out) + 1L]] <- data.frame(
        eq = eqn, block = bk, k_removed = dfC,
        J_full = unname(Jf["J"]), dfJ_full = unname(Jf["df"]),
        J_restr = unname(Jr["J"]), dfJ_restr = unname(Jr["df"]),
        C_stat = C, df_C = dfC, p_C = pC, stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) {
    return(data.frame(
      eq=character(), block=character(), k_removed=integer(),
      J_full=double(), dfJ_full=integer(), J_restr=double(), dfJ_restr=integer(),
      C_stat=double(), df_C=integer(), p_C=double()
    ))
  }
  dplyr::bind_rows(out)
}

# Sugestão de poda: blocos com p_C < alpha
suggest_prune <- function(ctab, alpha = 0.05) {
  if (!NROW(ctab)) return(character(0))
  unique(ctab$block[is.finite(ctab$p_C) & ctab$p_C < alpha])
}

# Nova lista de instrumentos após remover blocos
new_inst_after_prune <- function(fit_obj, to_drop_blocks, blocks) {
  inst_base <- .build_inst_terms(fit_obj)
  rhs_all   <- .clean_terms(fit_obj$rhs_terms)
  endo_ln   <- grep("^ln_", rhs_all, value = TRUE)
  exogs     <- setdiff(rhs_all, endo_ln)      # sempre ficam
  drop_set  <- unique(unlist(blocks[to_drop_blocks], use.names = FALSE))
  excl_ivs  <- setdiff(inst_base, exogs)
  excl_new  <- setdiff(excl_ivs, drop_set)
  unique(c(exogs, excl_new))
}

## ====== Execução/Exemplo de uso ======

rhs_all  <- fit_q$rhs_terms
inst_all <- fit_q$inst_terms
# Idealmente vazio (exógenos do RHS devem estar em Z). Se aparecer algo, .build_inst_terms resolve.
print(setdiff(rhs_all, inst_all))

## 1) J por equação (robusto com cluster se existir)
jtab <- jtests_by_eq_manual(
  fit_q, df_iv,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"
)
print(jtab)

## 2) C-tests por blocos (diferença de Hansen)
blocks <- make_blocks_explicit(fit_q)
message("Blocos detectados: ",
        if (length(blocks)) paste0(names(blocks), " (", vapply(blocks, length, integer(1)), " IVs)", collapse=" | ")
        else "nenhum (sem IVs excluídos)")

ctab <- c_tests_blocks_manual(
  fit_q, df_iv, blocks = blocks,
  cluster_var = if ("psu" %in% names(df_iv)) "psu" else "uf"
)

# Ordena e inspeciona
if (nrow(ctab)) {
  ctab <- dplyr::arrange(ctab, eq, p_C)
  print(ctab, digits = 4)
} else {
  message("ATENÇÃO: nenhum bloco testável (provavelmente sem IVs excluídos ou sem sobre-ID).")
}

## 3) Sugestão de poda e nova lista de IVs
to_drop <- suggest_prune(ctab, alpha = 0.05)
cat("\nBlocos sugeridos para poda (p_C < 0.05): ",
    if (length(to_drop)) paste(to_drop, collapse=", ") else "(nenhum)", "\n", sep="")

inst_new <- new_inst_after_prune(fit_q, to_drop, blocks)
cat("\nInst terms (NOVOS) após poda:\n", paste(inst_new, collapse = " + "), "\n")

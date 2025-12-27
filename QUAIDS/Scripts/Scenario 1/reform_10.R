## ========= 0) Pré-requisitos =========
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("glmnet","pls","AER","sandwich","lmtest","dplyr","purrr","tibble","systemfit","numDeriv"))

stopifnot(exists("df_iv_ok"), exists("priceNames"), exists("shareNames"),
          exists("best_fit"), exists("best_iv_set"))
drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))
exogs      <- intersect(c("z","z2"), names(df_iv_ok))
eqs        <- setdiff(shareNames, omit_share)

## ========= 1) Dicionário expandido de IVs =========
augment_dictionary <- function(df, base_iv, with=c("SH_SH_area_1","SH_SH_capital_1"),
                               poly_deg=3, center=TRUE, scale.=TRUE){
  with <- intersect(with, names(df))
  out_names <- character(0)
  for (b in base_iv) if (b %in% names(df)) {
    v <- df[[b]]
    if (center) v <- v - mean(v, na.rm=TRUE)
    if (scale.) v <- v / sd(v, na.rm=TRUE)
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

aug <- augment_dictionary(df_iv_ok, base_iv=intersect(best_iv_set, names(df_iv_ok)),
                          with=c("SH_SH_area_1","SH_SH_capital_1"), poly_deg=3)
df_big  <- aug$data
iv_pool <- aug$iv_names

## ========= 2) Helpers de residualização =========
resid_on <- function(v, X){
  if (is.null(X) || ncol(as.matrix(X))==0) return(as.numeric(v))
  X <- as.matrix(X); v <- as.numeric(v)
  v - X %*% solve(crossprod(X), crossprod(X, v))
}
resid_cols_on <- function(Z, X){
  Z <- as.matrix(Z)
  if (is.null(X) || ncol(as.matrix(X))==0) return(Z)
  X <- as.matrix(X); Z - X %*% solve(crossprod(X), crossprod(X, Z))
}

## ========= 3) SW-PLS por preço (1 PC que mira o SW_j) =========
make_sw_pc_for_price <- function(poi, df, iv_pool, endogs_all, exogs){
  stopifnot(poi %in% names(df))
  Znames <- intersect(iv_pool, names(df))
  Dminus <- setdiff(endogs_all, poi)
  Xn     <- intersect(exogs, names(df))
  use    <- unique(c(endogs_all, Xn, Znames))
  dat    <- df[stats::complete.cases(df[,use,drop=FALSE]), , drop=FALSE]
  
  ## 1) Prever D_-j via ridge usando X + Z (todas as informações IV)
  Dhat <- NULL
  if (length(Dminus)) {
    XZ <- as.matrix(dat[, c(Xn, Znames), drop=FALSE])
    Dhat <- sapply(Dminus, function(d){
      y <- as.numeric(dat[[d]])
      fit <- glmnet::cv.glmnet(XZ, y, alpha=0, intercept=TRUE, standardize=TRUE)
      as.numeric(predict(fit, XZ, s="lambda.min"))
    })
    Dhat <- as.data.frame(Dhat); colnames(Dhat) <- paste0(Dminus,"_hat")
  }
  
  ## 2) Residualizar D_j e Z em [X, D_-j_hat] (mira o SW_j)
  Xc <- if (!is.null(Dhat)) cbind(dat[, Xn, drop=FALSE], Dhat) else dat[, Xn, drop=FALSE]
  Dj_t <- resid_on(dat[[poi]], Xc)
  Z_t  <- resid_cols_on(dat[, Znames, drop=FALSE], Xc)
  
  ## 3) 1 componente PLS
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
    out_df  <- pc$data
    pc_names <- c(pc_names, pc$pc_name)
  }
  list(data=out_df, iv_star=pc_names)
}

pcs_sw <- gen_all_pcs_sw(df_big, priceNames, drop_price, iv_pool, exogs)
dfSW    <- pcs_sw$data
iv_star <- pcs_sw$iv_star  # 1 PC por preço (just-ID)

## ========= 4) SW_j correto com os novos PCs =========
sw_F_one <- function(j, Dnames, Xnames, Znames, data){
  Dminus <- setdiff(Dnames, j)
  # 1ª etapas p/ D_-j
  Dhat <- NULL
  if (length(Dminus)){
    fits <- lapply(Dminus, function(d) lm(reformulate(c(Xnames, Znames), response=d), data=data))
    Dhat <- as.data.frame(sapply(fits, fitted))
    colnames(Dhat) <- paste0(Dminus, "_hat")
  }
  dat2 <- if (!is.null(Dhat)) cbind(data, Dhat) else data
  
  f_red  <- lm(reformulate(c(Xnames, if(!is.null(Dhat)) names(Dhat) else NULL), response=j), data=dat2)
  f_full <- lm(reformulate(c(Xnames, Znames, if(!is.null(Dhat)) names(Dhat) else NULL), response=j), data=dat2)
  as.numeric(anova(f_red, f_full)$F[2])
}

SW_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  use <- unique(c(endogs_all, exogs, iv_star))
  dat <- dfSW[complete.cases(dfSW[,use,drop=FALSE]), , drop=FALSE]
  tibble::tibble(price=Dv, n=nrow(dat),
                 F_SW_j = sw_F_one(Dv, endogs_all, exogs, iv_star, dat))
})
cat("\n[SW_j | PCs otimizados p/ SW]\n"); print(SW_tbl, digits=3)

## ========= 5) Reestimar 2SLS/3SLS + elasticidades/ICs =========
# 2SLS eq-a-eq (cluster robusto)
cluster_var <- intersect(c("cluster_id","id_domicilio","id_municipio","id_setor_censitario"), names(dfSW))
if (!length(cluster_var)) { dfSW$.__rowid__ <- seq_len(nrow(dfSW)); cluster_var <- ".__rowid__" }
cluster_var <- cluster_var[1]

fit_ivreg_one <- function(dat, y, endogs, exogs, iv_names, cluster){
  fy <- stats::reformulate(c(endogs, exogs), response=y)
  fz <- stats::reformulate(unique(c(exogs, iv_names)))
  fit <- AER::ivreg(fy, instruments=fz, data=dat)
  vc  <- tryCatch(sandwich::vcovCL(fit, cluster=dat[[cluster]]),
                  error=function(e) sandwich::vcovHC(fit, type="HC1"))
  list(fit=fit, vcov=vc)
}

ivsys <- lapply(eqs, function(y){
  use <- unique(c(y, endogs_all, exogs, iv_star, cluster_var))
  dat <- dfSW[complete.cases(dfSW[,use,drop=FALSE]), , drop=FALSE]
  info <- fit_ivreg_one(dat, y, endogs_all, exogs, iv_star, cluster_var)
  list(eq=y, n=nrow(dat), fit=info$fit, vcov=info$vcov)
}); names(ivsys) <- eqs

gamma_from_ivsys <- function(ivsys, priceNames, drop_price, shareNames, omit_share){
  eqs <- setdiff(shareNames, omit_share)
  pj  <- paste0("ln_", setdiff(priceNames, drop_price))
  G   <- matrix(NA_real_, nrow=length(shareNames), ncol=length(pj),
                dimnames=list(shareNames, pj))
  for (eq in eqs) {
    co <- coef(ivsys[[eq]]$fit); G[eq, pj] <- unname(co[pj])
  }
  G[omit_share, pj] <- -colSums(G[eqs, pj, drop=FALSE], na.rm=TRUE)
  G
}
Giv <- gamma_from_ivsys(ivsys, priceNames, drop_price, shareNames, omit_share)

fitJI <- best_fit
cols_fit <- colnames(fitJI$coef$gamma)
pj_short <- sub("^ln_", "", colnames(Giv))
fitJI$coef$gamma[, match(pj_short, cols_fit)] <- Giv

x_star <- if ("gasto_total_atualhat" %in% names(dfSW))
  median(dfSW$gasto_total_atualhat, na.rm=TRUE) else
    median(dfSW[[grep("gasto|exp", names(dfSW), value=TRUE)[1]]], na.rm=TRUE)
p_star <- exp(colMeans(dfSW[paste0("ln_", priceNames)], na.rm=TRUE))

E <- elas_quaids_manual(
  fitJI, x=x_star, p=p_star,
  normalize_eval="softmax", normalize_deriv="softmax",
  dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
)
cat("\n[Marshall com SW-PCs]:\n"); print(round(E$marshall,3))
cat("\n[Hicks com SW-PCs]:\n"); print(round(E$hicks,3))
cat("\n[Eta]:\n"); print(round(E$expenditure,3))

## (Opcional) Se quiser manter os CIs por Delta (3SLS), reaproveite seu bloco v7
## bastando trocar dat_sys <- dfSW e usar os mesmos iv_star.

## ------------ 0) Utilidades genéricas ------------
replace_gamma_return_E <- function(fit_template, Giv, x_star, p_star,
                                   normalize_eval="softmax", normalize_deriv="softmax"){
  fit2 <- fit_template
  cols_fit <- colnames(fit2$coef$gamma)           # ex: "preco_com_reforma12"
  pj_short <- sub("^ln_", "", colnames(Giv))
  fit2$coef$gamma[, match(pj_short, cols_fit)] <- Giv
  elas_quaids_manual(
    fit2, x = x_star, p = p_star,
    normalize_eval = normalize_eval, normalize_deriv = normalize_deriv,
    dlogp = 1e-3, enforce_hicks = TRUE, enforce_symmetry = TRUE
  )
}

## (gradiente central — não usado, deixo para debug se quiser)
numgrad_central <- function(fun, theta, eps = 1e-6){
  f0 <- fun(theta)
  m  <- length(f0); k <- length(theta)
  J  <- matrix(NA_real_, nrow = m, ncol = k)
  for (j in seq_len(k)) {
    thp <- theta; thm <- theta
    h   <- eps * (abs(theta[j]) + 1)
    thp[j] <- theta[j] + h
    thm[j] <- theta[j] - h
    fp <- fun(thp); fm <- fun(thm)
    J[, j] <- (fp - fm) / (2*h)
  }
  J
}

## ------------ 1) Objetos do pipeline ------------
# simples e direto, se já rodou o bloco SW-PCs:
dfJI    <- dfSW
iv_star <- pcs_sw$iv_star
stopifnot(all(iv_star %in% names(dfJI)),
          length(iv_star) == length(endogs_all))


drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))   # todos ln preços (exceto dropado)
eqs        <- setdiff(shareNames, omit_share)                  # shares estimadas
exogs      <- intersect(c("z","z2"), names(dfJI))              # exógenos no RHS (ajuste se precisar)

## IVs just-ID (precisa 1 PC por preço)
stopifnot(all(iv_star %in% names(dfJI)),
          length(iv_star) == length(endogs_all))

## ------------ 2) Monta e estima o sistema 3SLS just-ID (labels seguros) ------------
eqs_orig <- setdiff(shareNames, omit_share)
eqs_safe <- paste0("eq", seq_along(eqs_orig))                  # labels sem espaço/underscore
orig2safe <- setNames(eqs_safe, eqs_orig)
safe2orig <- setNames(eqs_orig, eqs_safe)

rhs_y <- paste(c(endogs_all, exogs), collapse = " + ")
rhs_z <- paste(c(exogs, iv_star),  collapse = " + ")

eq_list   <- setNames(vector("list", length(eqs_orig)), eqs_safe)
inst_list <- setNames(vector("list", length(eqs_orig)), eqs_safe)

for (i in seq_along(eqs_orig)) {
  y_orig  <- eqs_orig[i]
  y_safe  <- orig2safe[[y_orig]]
  eq_list[[y_safe]]   <- stats::as.formula(paste(y_orig, "~", rhs_y))
  inst_list[[y_safe]] <- stats::as.formula(paste("~", rhs_z))
}

use_all <- unique(c(eqs_orig, endogs_all, exogs, iv_star))
dfsw <- dfJI[stats::complete.cases(dfJI[, use_all, drop = FALSE]), , drop = FALSE]

fit_sys <- systemfit::systemfit(eq_list, method = "3SLS", inst = inst_list, data = dfsw)

## ------------ 3) Mapeamento robusto eq-systemfit ↔ shares ------------
## --- PATCH 1: diagnóstico robusto de labels/LHS no systemfit ---
diag_eq_map <- function(fit_sys){
  # 1) Tenta pegar via fit_sys$eq (caminho "bonito")
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
    return(data.frame(eq_label = labs,
                      lhs      = lhs,
                      key_label= toupper(gsub("[^A-Za-z0-9]","", labs)),
                      key_lhs  = toupper(gsub("[^A-Za-z0-9]","", lhs)),
                      stringsAsFactors = FALSE))
  }
  
  # 2) Fallback: reconstrói labels a partir de names(coef(fit_sys))
  cn <- names(coef(fit_sys))
  if (length(cn) == 0) {
    return(data.frame(eq_label=character(0), lhs=character(0),
                      key_label=character(0), key_lhs=character(0)))
  }
  # padrão eq<lab>_<var> ou eq<lab>.<var>
  eq_guess <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  data.frame(eq_label = eq_labels,
             lhs       = NA_character_,                  # sem acesso fácil ao LHS aqui
             key_label = toupper(gsub("[^A-Za-z0-9]","", eq_labels)),
             key_lhs   = NA_character_,
             stringsAsFactors = FALSE)
}

.canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))


## --- PATCH 2: mapeamento super-robusto (com fallback por ordem) ---
auto_eq_map2 <- function(fit_sys, shareNames, omit_share){
  diag <- diag_eq_map(fit_sys)
  shares_ok  <- setdiff(shareNames, omit_share)
  key_shares <- .canon(shares_ok)
  
  # Tentativa A/B/C como antes (casa por LHS literal / LHS canônico / LABEL canônico)
  hitA <- if (!all(is.na(diag$lhs))) match(diag$lhs, shares_ok,  nomatch = 0L) else rep(0L, nrow(diag))
  hitB <- if (!all(is.na(diag$key_lhs))) match(diag$key_lhs, key_shares, nomatch = 0L) else rep(0L, nrow(diag))
  hitC <- match(diag$key_label, key_shares, nomatch = 0L)
  
  idx <- unique(c(which(hitA>0), which(hitB>0), which(hitC>0)))
  
  if (length(idx) > 0) {
    eqs_safe <- diag$eq_label[idx]
    eqs_orig <- shares_ok[
      ifelse(idx %in% which(hitA>0), hitA[idx],
             ifelse(idx %in% which(hitB>0), hitB[idx], hitC[idx]))
    ]
    dup <- duplicated(eqs_orig)
    if (any(dup)) { eqs_safe <- eqs_safe[!dup]; eqs_orig <- eqs_orig[!dup] }
    return(list(eqs_safe = eqs_safe, eqs_orig = eqs_orig, diag = diag))
  }
  
  # Fallback FINAL: casa por ORDEM (eqs como aparecem nos coeficientes) ↔ shares_ok
  cn <- names(coef(fit_sys))
  if (length(cn) == 0)
    stop("auto_eq_map2: não foi possível inspecionar coeficientes do systemfit.")
  eq_guess <- sub("(_|\\.).*$", "", cn)
  eq_labels <- unique(eq_guess)
  Ke <- length(eq_labels)
  if (Ke != length(shares_ok)) {
    stop(sprintf("auto_eq_map2: fallback por ordem falhou: #eq=%d != #shares=%d.",
                 Ke, length(shares_ok)))
  }
  message("[auto_eq_map2] Usando fallback por ORDEM: ",
          paste(sprintf("%s↔%s", eq_labels, shares_ok), collapse = "; "))
  list(eqs_safe = eq_labels, eqs_orig = shares_ok,
       diag = data.frame(eq_label = eq_labels, lhs = NA_character_,
                         key_label = .canon(eq_labels), key_lhs = NA_character_))
}


## ------------ 4) Extrator θ=vec(Γ) e Var(θ) (robusto a nomes) ------------
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
    if (is.na(idx))
      stop(sprintf("Coeficiente não encontrado: eq='%s' var='%s'\nVeja names(coef(fit_sys)).",
                   eq_safe[i], pj[j]))
    theta[pos] <- unname(co[idx]); nm_th[pos] <- paste(eq_orig[i], pj[j], sep="::")
  }
  
  for (i in seq_along(eq_safe)) for (k in seq_along(eq_safe)) {
    for (a in seq_along(pj)) for (b in seq_along(pj)) {
      pos1 <- (i - 1L) * Kp + a; pos2 <- (k - 1L) * Kp + b
      if (is.na(Vth[pos1,pos2])) {
        idx1 <- find_idx(eq_safe[i], pj[a]); idx2 <- find_idx(eq_safe[k], pj[b])
        Vth[pos1,pos2] <- V[idx1, idx2]
      }
    }
  }
  names(theta) <- nm_th
  list(theta = theta, V = Vth)
}

## ------------ 5) Delta Method para ICs (sem bootstrap) ------------
## ================= DELTA: ICs para Hicks + Marshall + Eta =================
## ================= DELTA: ICs p/ Hicks + Marshall + Eta (robusto a nomes) ==============
# ==================== Delta CIs v7 (com Marshall FULL) ====================
elasticity_CIs_delta_system_v7 <- function(
    fit_sys, fit_template, dfsw,
    priceNames, shareNames, drop_price, omit_share,
    level = 0.95, method = c("pointwise","bonferroni")
){
  if (!requireNamespace("numDeriv", quietly = TRUE))
    stop("Instale: install.packages('numDeriv')")
  if (!requireNamespace("tibble", quietly = TRUE))
    stop("Instale: install.packages('tibble')")
  method <- match.arg(method)
  `%||%` <- function(a,b) if(!is.null(a)) a else b
  .canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x %||% ""))
  
  ## --- 1) Mapear equações do systemfit -> suas shares (com fallback) ---
  map <- auto_eq_map2(fit_sys, shareNames, omit_share)  # sua helper
  eqs_safe <- map$eqs_safe
  eqs_orig <- map$eqs_orig
  
  ## --- 2) Empilhar theta = vec(Gamma) e Var(theta) ---
  pj       <- paste0("ln_", setdiff(priceNames, drop_price))
  pj_short <- sub("^ln_", "", pj)
  meta     <- stack_gamma_vcov_sys2(fit_sys, pj, eqs_safe, eqs_orig) # sua helper
  theta    <- meta$theta
  Vth      <- meta$V
  
  ## --- 3) De theta -> Gfull (com adding-up para a share omitida) ---
  theta_to_G <- function(th){
    Ke <- length(eqs_orig); Kp <- length(pj)
    Ghat <- matrix(th, nrow = Ke, byrow = TRUE, dimnames = list(eqs_orig, pj))
    Gfull <- matrix(NA_real_, nrow = length(shareNames), ncol = Kp,
                    dimnames = list(shareNames, pj))
    Gfull[eqs_orig, pj] <- Ghat
    Gfull[omit_share, pj] <- -colSums(Ghat, na.rm = TRUE)
    Gfull
  }
  
  ## --- 4) Ponto de avaliação (x*, p*) ---
  x_star <- if ("gasto_total_atualhat" %in% names(dfsw))
    median(dfsw$gasto_total_atualhat, na.rm=TRUE) else
      median(dfsw[[grep("gasto|exp", names(dfsw), value=TRUE)[1]]], na.rm=TRUE)
  p_star <- exp(colMeans(dfsw[paste0("ln_", priceNames)], na.rm=TRUE))
  
  ## --- 5) Avaliação em theta para capturar ordem/nome “nativos” ---
  G0 <- theta_to_G(theta)
  f2 <- fit_template
  cols_fit <- colnames(f2$coef$gamma)
  f2$coef$gamma[, match(pj_short, cols_fit)] <- G0
  E0 <- elas_quaids_manual(
    f2, x=x_star, p=p_star,
    normalize_eval="softmax", normalize_deriv="softmax",
    dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
  )
  
  rn0 <- rownames(E0$hicks); cn0 <- colnames(E0$hicks)
  if (is.null(rn0)) rn0 <- shareNames
  if (is.null(cn0)) cn0 <- shareNames
  k0  <- length(rn0)
  
  ## --- 6) Função g(theta) com vetor na ordem: Hicks FULL, Marshall FULL, Eta ---
  ##     IMPORTANTe: vetorizar por LINHA: as.vector(t(M)) para reconstruir com byrow=TRUE
  g_eval <- function(th){
    G <- theta_to_G(th)
    f3 <- fit_template
    f3$coef$gamma[, match(pj_short, cols_fit)] <- G
    E <- elas_quaids_manual(
      f3, x=x_star, p=p_star,
      normalize_eval="softmax", normalize_deriv="softmax",
      dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
    )
    c(as.vector(t(E$hicks)),
      as.vector(t(E$marshall)),
      as.numeric(E$expenditure))
  }
  
  g0 <- g_eval(theta)
  J  <- numDeriv::jacobian(g_eval, theta)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  lo <- g0 - z*se; hi <- g0 + z*se
  
  ## --- 7) Reconstituir blocos NATIVOS (sem alinhamento) ---
  idx_H    <- seq_len(k0*k0)
  idx_Mall <- (k0*k0) + seq_len(k0*k0)
  idx_eta  <- (k0*k0 + k0*k0) + seq_len(k0)
  
  H_hat0 <- matrix(g0[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_lo0  <- matrix(lo[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_hi0  <- matrix(hi[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_se0  <- matrix(se[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  M_hat0 <- matrix(g0[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_lo0  <- matrix(lo[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_hi0  <- matrix(hi[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_se0  <- matrix(se[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  eta_hat0 <- setNames(g0[idx_eta], rn0)
  eta_lo0  <- setNames(lo[idx_eta], rn0)
  eta_hi0  <- setNames(hi[idx_eta], rn0)
  eta_se0  <- setNames(se[idx_eta], rn0)
  
  ## --- 8) Alinhamento p/ shareNames (se houver match) ---
  key_tar <- .canon(shareNames)
  key_r   <- .canon(rownames(H_hat0))
  key_c   <- .canon(colnames(H_hat0))
  ridx <- match(key_tar, key_r)
  cidx <- match(key_tar, key_c)
  
  all_miss <- all(is.na(ridx)) && all(is.na(cidx))
  k <- length(shareNames)
  
  mkMat <- function() matrix(NA_real_, nrow=k, ncol=k, dimnames=list(shareNames, shareNames))
  H_hat <- mkMat(); H_lo <- mkMat(); H_hi <- mkMat(); H_se <- mkMat()
  M_hat <- mkMat(); M_lo <- mkMat(); M_hi <- mkMat(); M_se <- mkMat()
  
  ok_r <- which(!is.na(ridx)); ok_c <- which(!is.na(cidx))
  if (length(ok_r) && length(ok_c)) {
    H_hat[ok_r, ok_c] <- H_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_lo [ok_r, ok_c] <- H_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_hi [ok_r, ok_c] <- H_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_se [ok_r, ok_c] <- H_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    
    M_hat[ok_r, ok_c] <- M_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_lo [ok_r, ok_c] <- M_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_hi [ok_r, ok_c] <- M_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_se [ok_r, ok_c] <- M_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
  }
  
  map1 <- match(key_tar, key_r)
  Mdiag_hat <- Mdiag_lo <- Mdiag_hi <- Mdiag_se <- rep(NA_real_, k)
  names(Mdiag_hat) <- names(Mdiag_lo) <- names(Mdiag_hi) <- names(Mdiag_se) <- shareNames
  ok <- which(!is.na(map1))
  if (length(ok)) {
    Mdiag_hat[ok] <- diag(M_hat0)[ map1[ok] ]
    Mdiag_lo [ok] <- diag(M_lo0)[ map1[ok] ]
    Mdiag_hi [ok] <- diag(M_hi0)[ map1[ok] ]
    Mdiag_se [ok] <- diag(M_se0)[ map1[ok] ]
  }
  
  eta_hat <- eta_lo <- eta_hi <- eta_se <- rep(NA_real_, k)
  names(eta_hat) <- names(eta_lo) <- names(eta_hi) <- names(eta_se) <- shareNames
  if (length(ok)) {
    eta_hat[ok] <- eta_hat0[ map1[ok] ]
    eta_lo [ok] <- eta_lo0 [ map1[ok] ]
    eta_hi [ok] <- eta_hi0 [ map1[ok] ]
    eta_se [ok] <- eta_se0 [ map1[ok] ]
  }
  
  ## --- 9) Saídas em formato TIDY (alinhado e nativo) ---
  to_tbl <- function(M, lo, hi, se, rows, cols){
    do.call(rbind, lapply(seq_along(rows), function(i)
      data.frame(i=rows[i], j=cols, hat=M[i,], lo=lo[i,], hi=hi[i,], se=se[i,], row.names=NULL)))
  }
  
  hicks_tbl_aligned    <- to_tbl(H_hat, H_lo, H_hi, H_se, shareNames, shareNames)
  hicks_tbl_native     <- to_tbl(H_hat0, H_lo0, H_hi0, H_se0, rn0, cn0)
  marshall_tbl_aligned <- to_tbl(M_hat, M_lo, M_hi, M_se, shareNames, shareNames)
  marshall_tbl_native  <- to_tbl(M_hat0, M_lo0, M_hi0, M_se0, rn0, cn0)
  
  marshall_diag_tbl <- tibble::tibble(
    good = shareNames, M_hat = Mdiag_hat, M_lo = Mdiag_lo, M_hi = Mdiag_hi, M_se = Mdiag_se
  )
  eta_tbl <- tibble::tibble(
    good = shareNames, eta_hat = eta_hat, eta_lo = eta_lo, eta_hi = eta_hi, eta_se = eta_se
  )
  hicks_diag_tbl <- tibble::tibble(
    good = shareNames,
    H_hat = diag(H_hat), H_lo = diag(H_lo), H_hi = diag(H_hi), H_se = diag(H_se)
  )
  
  note <- NULL
  if (all_miss) {
    note <- "Nomes não casaram: use *_native (ordem e rótulos nativos do E)."
    message("[Delta-elasticities v7] ", note)
  }
  
  list(
    ## Hicks
    hicks_full        = hicks_tbl_aligned,
    hicks_full_native = hicks_tbl_native,
    hicks_diag        = hicks_diag_tbl,
    ## Marshall (FULL + DIAG) — agora com ICs completos
    marshall_full        = marshall_tbl_aligned,
    marshall_full_native = marshall_tbl_native,
    marshall_diag        = marshall_diag_tbl,
    ## Expenditure
    eta              = eta_tbl,
    level            = level,
    method           = method,
    mapping_diag     = map$diag,
    note             = note
  )
}
# ================== FIM Delta CIs v7 (com Marshall FULL) ==================

## ======================= FIM v5 =======================

## ========================= FIM DO PATCH =========================

## ------------ 6) Rodar e imprimir ------------
# 6.1 Diagnóstico de como o systemfit rotulou as equações
cat("\n[systemfit] labels e LHS por equação:\n")
print(diag_eq_map(fit_sys), row.names = FALSE)

# 6.2 ICs por Delta (sem bootstrap)
cis_sys <- elasticity_CIs_delta_system_v7(
  fit_sys      = fit_sys,
  fit_template = best_fit,
  dfsw      = dfsw,
  priceNames   = priceNames,
  shareNames   = shareNames,
  drop_price   = drop_price,
  omit_share   = omit_share,
  level        = 0.95,
  method       = "bonferroni"  # ou "pointwise"
)

cat("\n[ICs Delta — Marshall]\n"); print(cis_sys$marshall_full_native, digits=3)
cat("\n[ICs Delta — Hicks]\n"); print(cis_sys$hicks_full_native, digits = 3)
cat("\n[ICs Delta — Eta]\n");             print(cis_sys$eta,          digits=3)

# D = endógenos (todos ln preços), X = exógenos, Z = IVs excluídos (PCs)
sw_F_one <- function(j, Dnames, Xnames, Znames, data){
  Dminus <- setdiff(Dnames, j)
  # 1ª etapa para cada D_-j
  fits <- lapply(Dminus, function(d)
    lm(reformulate(c(Xnames, Znames), response=d), data=data))
  Dhat <- as.data.frame(sapply(fits, fitted))
  colnames(Dhat) <- paste0(Dminus, "_hat")
  dat2 <- cbind(data, Dhat)
  
  # Regressão “auxiliar” de D_j em X + Z + D_-j_hat
  f_red  <- lm(reformulate(c(Xnames,        names(Dhat)), response=j), data=dat2)
  f_full <- lm(reformulate(c(Xnames, Znames, names(Dhat)), response=j), data=dat2)
  as.numeric(anova(f_red, f_full)$F[2])  # este é o SW F_j
}

SW_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  use <- unique(c(endogs_all, exogs, iv_star))
  dat <- dfJI[complete.cases(dfJI[, use, drop=FALSE]), , drop=FALSE]
  tibble::tibble(price=Dv, n=nrow(dat),
                 F_SW_j = sw_F_one(Dv, endogs_all, exogs, iv_star, dat))
})
print(SW_tbl, digits=3)

## =========================
## POLIMENTO TOP TIER (3x)
## =========================
need_pkg <- function(pkgs){
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("Instale: ", paste(miss, collapse=", "))
}
need_pkg(c("AER","sandwich","tibble","dplyr"))
has_ivmodel <- requireNamespace("ivmodel", quietly = TRUE)

stopifnot(exists("dfSW"), exists("iv_star"),
          exists("priceNames"), exists("shareNames"),
          exists("best_fit"))
drop_price <- best_fit$drop_price
omit_share <- best_fit$omit_share
endogs_all <- paste0("ln_", setdiff(priceNames, drop_price))
exogs      <- intersect(c("z","z2"), names(dfSW))
eqs        <- setdiff(shareNames, omit_share)

## ---------------------------------------------------------
## (1) AR/CLR com SW-PCs (weak-IV robust) p/ um 'poi'
## ---------------------------------------------------------
ar_clr_oneeq_sw <- function(eq_y, poi, halfwidth=200){
  ln_all    <- paste0("ln_", priceNames)
  ln_others <- setdiff(ln_all, c(poi, paste0("ln_", drop_price)))
  Xnames <- c(exogs, ln_others)
  
  use <- unique(c(eq_y, poi, Xnames, iv_star))
  dat <- dfSW[stats::complete.cases(dfSW[, use, drop=FALSE]), , drop=FALSE]
  
  Y <- as.numeric(dat[[eq_y]])
  D <- as.numeric(dat[[poi]])
  X <- if(length(Xnames)) as.matrix(dat[, Xnames, drop=FALSE]) else NULL
  Z <- as.matrix(dat[, iv_star, drop=FALSE])
  
  # centro no 2SLS
  b2 <- tryCatch(unname(coef(AER::ivreg(Y ~ D + X | X + Z))["D"]), error=function(e) 0)
  grid <- seq(b2 - halfwidth, b2 + halfwidth, length.out=4001)
  
  if (has_ivmodel && is.numeric(D) && length(D)==nrow(dat)) {
    ivm <- ivmodel::ivmodel(Y=Y, D=D, Z=Z, X=X)
    p_ar  <- vapply(grid, function(b0) ivmodel::AR.test(ivm, beta0=b0)$p.value, 0.0)
    p_clr <- vapply(grid, function(b0) ivmodel::CLR(ivm,     beta0=b0)$p.value, 0.0)
    keepA <- which(p_ar  >= 0.05); keepC <- which(p_clr >= 0.05)
    ciA <- if(length(keepA)) c(min(grid[keepA]), max(grid[keepA])) else c(NA_real_, NA_real_)
    ciC <- if(length(keepC)) c(min(grid[keepC]), max(grid[keepC])) else c(NA_real_, NA_real_)
  } else {
    # fallback robusto (multi-endógenos): exclusão de Z no empilhado
    p_ar <- vapply(grid, function(b0){
      fit <- lm(cbind(Y - b0*D, X) ~ 0 + Z)
      1 - pf(summary(aov(fit))[[1]][["F value"]][1], ncol(Z), nrow(Z) - ncol(Z))
    }, 0.0)
    keepA <- which(p_ar >= 0.05)
    ciA <- if(length(keepA)) c(min(grid[keepA]), max(grid[keepA])) else c(NA_real_, NA_real_)
    ciC <- ciA
  }
  tibble::tibble(eq=eq_y, AR_low=ciA[1], AR_high=ciA[2], CLR_low=ciC[1], CLR_high=ciC[2])
}

# rode para o seu preço de interesse (mantém 'poi' se já existir)
price_of_interest <- if (exists("poi")) poi else "ln_preco_com_reforma13"
arclr_sw_tbl <- purrr::map_dfr(eqs, ~ar_clr_oneeq_sw(.x, price_of_interest, halfwidth=200))
cat("\n[AR/CLR | SW-PCs | coef de ", price_of_interest, "]\n", sep=""); print(arclr_sw_tbl, digits=3)

## -----------------------------------------------------------------
## (2) Mini-tabela: SW F_j + R² parcial_j (definição Stock–Yogo j)
##     Regressão de D_j em [X, D_-j_hat] vs [X, D_-j_hat, Z]
## -----------------------------------------------------------------
sw_F_R2_one <- function(j, Dnames, Xnames, Znames, data){
  Dminus <- setdiff(Dnames, j)
  
  # 1ª etapas p/ D_-j com X+Z
  Dhat <- NULL
  if (length(Dminus)){
    fits <- lapply(Dminus, function(d) lm(reformulate(c(Xnames, Znames), response=d), data=data))
    Dhat <- as.data.frame(sapply(fits, fitted))
    colnames(Dhat) <- paste0(Dminus, "_hat")
  }
  
  dat2 <- if (!is.null(Dhat)) cbind(data, Dhat) else data
  Xred <- c(Xnames, if(!is.null(Dhat)) names(Dhat) else NULL)
  
  f_red  <- lm(reformulate(Xred,                 response=j), data=dat2)
  f_full <- lm(reformulate(c(Xred, Znames),      response=j), data=dat2)
  
  an  <- anova(f_red, f_full)
  F_j <- as.numeric(an$F[2])
  
  # R² parcial dos Z: 1 - SSR_full/SSR_red
  ssr_red  <- sum(residuals(f_red )^2)
  ssr_full <- sum(residuals(f_full)^2)
  R2p_j <- 1 - ssr_full/ssr_red
  
  tibble::tibble(F_SW_j = F_j, R2p_j = R2p_j,
                 df1 = length(Znames),
                 df2 = nrow(dat2) - length(Xred) - length(Znames) - 1)
}

SWmini_tbl <- purrr::map_dfr(endogs_all, function(Dv){
  use <- unique(c(endogs_all, exogs, iv_star))
  dat <- dfSW[complete.cases(dfSW[,use,drop=FALSE]), , drop=FALSE]
  out <- sw_F_R2_one(Dv, endogs_all, exogs, iv_star, dat)
  tibble::tibble(price=Dv, n=nrow(dat), !!!out)
})
cat("\n[SW mini | F_SW_j e R²_parcial_j (SW-PCs)]\n"); print(SWmini_tbl, digits=3)

## ----------------------------------------------------------------------
## (3) Congelar labels do systemfit (sem '_' ou ' ') e exibir mapeamento
## ----------------------------------------------------------------------
need_pkg(c("systemfit"))
eqs_orig <- setdiff(shareNames, omit_share)

# labels seguros: só letras/números, prefixados por 'EQ'
sanitize <- function(x){
  s <- gsub("[^A-Za-z0-9]", "", x)
  s <- ifelse(nzchar(s), s, paste0("X", seq_along(s))) # garante não vazio
  paste0("EQ", s)
}
eqs_safe <- sanitize(eqs_orig)
orig2safe <- setNames(eqs_safe, eqs_orig)
safe2orig <- setNames(eqs_orig, eqs_safe)

rhs_y <- paste(c(endogs_all, exogs), collapse = " + ")
rhs_z <- paste(c(exogs, iv_star),  collapse = " + ")

eq_list   <- setNames(vector("list", length(eqs_orig)), eqs_safe)
inst_list <- setNames(vector("list", length(eqs_orig)), eqs_safe)
for (i in seq_along(eqs_orig)) {
  y_orig <- eqs_orig[i]; y_safe <- orig2safe[[y_orig]]
  eq_list[[y_safe]]   <- stats::as.formula(paste(y_orig, "~", rhs_y))
  inst_list[[y_safe]] <- stats::as.formula(paste("~", rhs_z))
}

use_all <- unique(c(eqs_orig, endogs_all, exogs, iv_star))
dat_sys <- dfSW[stats::complete.cases(dfSW[, use_all, drop=FALSE]), , drop=FALSE]

fit_sys <- systemfit::systemfit(eq_list, method="3SLS", inst=inst_list, data=dat_sys)

cat("\n[systemfit] Mapeamento de labels congelados:\n")
map_df <- tibble::tibble(eq_label=names(eq_list), lhs_original=eqs_orig)
print(map_df, n=Inf)

# (opcional) se quiser rodar o Delta-CIs com estes labels, reaproveite seu bloco v7
# trocando 'dat_sys' (ou 'dfsw') para 'dat_sys' acima e mantendo iv_star.

## ========= Consistência R2p vs F_SW_j =========
## (usa a identidade R2p = (q*F) / (q*F + df2))
SWmini_tbl <- SWmini_tbl |>
  dplyr::mutate(
    R2p_from_F = (df1 * F_SW_j) / (df1 * F_SW_j + df2),
    gap = R2p_j - R2p_from_F
  )

cat("\n[Checagem: R2p_j (direto) vs R2p_from_F (identidade)]\n")
print(dplyr::select(SWmini_tbl, price, F_SW_j, df1, df2, R2p_j, R2p_from_F, gap), digits = 4)

## Se preferir padronizar pelo valor teórico da identidade:
SWmini_tbl <- SWmini_tbl |>
  dplyr::mutate(R2p_j = R2p_from_F) |>
  dplyr::select(-R2p_from_F, -gap)

cat("\n[SW mini | versão final (R2p ajustado pela identidade)]\n")
print(SWmini_tbl, digits = 4)

elasticity_CIs_delta_system_v7 <- function(
    fit_sys, fit_template, dat_sys,
    priceNames, shareNames, drop_price, omit_share,
    level = 0.95, method = c("pointwise","bonferroni")
){
  if (!requireNamespace("numDeriv", quietly = TRUE))
    stop("Instale: install.packages('numDeriv')")
  if (!requireNamespace("tibble", quietly = TRUE))
    stop("Instale: install.packages('tibble')")
  method <- match.arg(method)
  `%||%` <- function(a,b) if(!is.null(a)) a else b
  .canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x %||% ""))
  
  ## --- 1) Mapear equações do systemfit -> suas shares (com fallback) ---
  map <- auto_eq_map2(fit_sys, shareNames, omit_share)  # sua helper
  eqs_safe <- map$eqs_safe
  eqs_orig <- map$eqs_orig
  
  ## --- 2) Empilhar theta = vec(Gamma) e Var(theta) ---
  pj       <- paste0("ln_", setdiff(priceNames, drop_price))
  pj_short <- sub("^ln_", "", pj)
  meta     <- stack_gamma_vcov_sys2(fit_sys, pj, eqs_safe, eqs_orig) # sua helper
  theta    <- meta$theta
  Vth      <- meta$V
  
  ## --- 3) De theta -> Gfull (com adding-up para a share omitida) ---
  theta_to_G <- function(th){
    Ke <- length(eqs_orig); Kp <- length(pj)
    Ghat <- matrix(th, nrow = Ke, byrow = TRUE, dimnames = list(eqs_orig, pj))
    Gfull <- matrix(NA_real_, nrow = length(shareNames), ncol = Kp,
                    dimnames = list(shareNames, pj))
    Gfull[eqs_orig, pj] <- Ghat
    Gfull[omit_share, pj] <- -colSums(Ghat, na.rm = TRUE)
    Gfull
  }
  
  ## --- 4) Ponto de avaliação (x*, p*) ---
  x_star <- if ("gasto_total_atualhat" %in% names(dat_sys))
    median(dat_sys$gasto_total_atualhat, na.rm=TRUE) else
      median(dat_sys[[grep("gasto|exp", names(dat_sys), value=TRUE)[1]]], na.rm=TRUE)
  p_star <- exp(colMeans(dat_sys[paste0("ln_", priceNames)], na.rm=TRUE))
  
  ## --- 5) Avaliação em theta para capturar ordem/nome “nativos” ---
  G0 <- theta_to_G(theta)
  f2 <- fit_template
  cols_fit <- colnames(f2$coef$gamma)
  f2$coef$gamma[, match(pj_short, cols_fit)] <- G0
  E0 <- elas_quaids_manual(
    f2, x=x_star, p=p_star,
    normalize_eval="softmax", normalize_deriv="softmax",
    dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
  )
  
  rn0 <- rownames(E0$hicks); cn0 <- colnames(E0$hicks)
  if (is.null(rn0)) rn0 <- shareNames
  if (is.null(cn0)) cn0 <- shareNames
  k0  <- length(rn0)
  
  ## --- 6) Função g(theta) com vetor na ordem: Hicks FULL, Marshall FULL, Eta ---
  ##     IMPORTANTe: vetorizar por LINHA: as.vector(t(M)) para reconstruir com byrow=TRUE
  g_eval <- function(th){
    G <- theta_to_G(th)
    f3 <- fit_template
    f3$coef$gamma[, match(pj_short, cols_fit)] <- G
    E <- elas_quaids_manual(
      f3, x=x_star, p=p_star,
      normalize_eval="softmax", normalize_deriv="softmax",
      dlogp=1e-3, enforce_hicks=TRUE, enforce_symmetry=TRUE
    )
    c(as.vector(t(E$hicks)),
      as.vector(t(E$marshall)),
      as.numeric(E$expenditure))
  }
  
  g0 <- g_eval(theta)
  J  <- numDeriv::jacobian(g_eval, theta)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  lo <- g0 - z*se; hi <- g0 + z*se
  
  ## --- 7) Reconstituir blocos NATIVOS (sem alinhamento) ---
  idx_H    <- seq_len(k0*k0)
  idx_Mall <- (k0*k0) + seq_len(k0*k0)
  idx_eta  <- (k0*k0 + k0*k0) + seq_len(k0)
  
  H_hat0 <- matrix(g0[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_lo0  <- matrix(lo[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_hi0  <- matrix(hi[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  H_se0  <- matrix(se[idx_H],    nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  M_hat0 <- matrix(g0[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_lo0  <- matrix(lo[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_hi0  <- matrix(hi[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  M_se0  <- matrix(se[idx_Mall], nrow=k0, byrow=TRUE, dimnames=list(rn0, cn0))
  
  eta_hat0 <- setNames(g0[idx_eta], rn0)
  eta_lo0  <- setNames(lo[idx_eta], rn0)
  eta_hi0  <- setNames(hi[idx_eta], rn0)
  eta_se0  <- setNames(se[idx_eta], rn0)
  
  ## --- 8) Alinhamento p/ shareNames (se houver match) ---
  key_tar <- .canon(shareNames)
  key_r   <- .canon(rownames(H_hat0))
  key_c   <- .canon(colnames(H_hat0))
  ridx <- match(key_tar, key_r)
  cidx <- match(key_tar, key_c)
  
  all_miss <- all(is.na(ridx)) && all(is.na(cidx))
  k <- length(shareNames)
  
  mkMat <- function() matrix(NA_real_, nrow=k, ncol=k, dimnames=list(shareNames, shareNames))
  H_hat <- mkMat(); H_lo <- mkMat(); H_hi <- mkMat(); H_se <- mkMat()
  M_hat <- mkMat(); M_lo <- mkMat(); M_hi <- mkMat(); M_se <- mkMat()
  
  ok_r <- which(!is.na(ridx)); ok_c <- which(!is.na(cidx))
  if (length(ok_r) && length(ok_c)) {
    H_hat[ok_r, ok_c] <- H_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_lo [ok_r, ok_c] <- H_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_hi [ok_r, ok_c] <- H_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    H_se [ok_r, ok_c] <- H_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    
    M_hat[ok_r, ok_c] <- M_hat0[ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_lo [ok_r, ok_c] <- M_lo0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_hi [ok_r, ok_c] <- M_hi0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
    M_se [ok_r, ok_c] <- M_se0 [ridx[ok_r], cidx[ok_c], drop=FALSE]
  }
  
  map1 <- match(key_tar, key_r)
  Mdiag_hat <- Mdiag_lo <- Mdiag_hi <- Mdiag_se <- rep(NA_real_, k)
  names(Mdiag_hat) <- names(Mdiag_lo) <- names(Mdiag_hi) <- names(Mdiag_se) <- shareNames
  ok <- which(!is.na(map1))
  if (length(ok)) {
    Mdiag_hat[ok] <- diag(M_hat0)[ map1[ok] ]
    Mdiag_lo [ok] <- diag(M_lo0)[ map1[ok] ]
    Mdiag_hi [ok] <- diag(M_hi0)[ map1[ok] ]
    Mdiag_se [ok] <- diag(M_se0)[ map1[ok] ]
  }
  
  eta_hat <- eta_lo <- eta_hi <- eta_se <- rep(NA_real_, k)
  names(eta_hat) <- names(eta_lo) <- names(eta_hi) <- names(eta_se) <- shareNames
  if (length(ok)) {
    eta_hat[ok] <- eta_hat0[ map1[ok] ]
    eta_lo [ok] <- eta_lo0 [ map1[ok] ]
    eta_hi [ok] <- eta_hi0 [ map1[ok] ]
    eta_se [ok] <- eta_se0 [ map1[ok] ]
  }
  
  ## --- 9) Saídas em formato TIDY (alinhado e nativo) ---
  to_tbl <- function(M, lo, hi, se, rows, cols){
    do.call(rbind, lapply(seq_along(rows), function(i)
      data.frame(i=rows[i], j=cols, hat=M[i,], lo=lo[i,], hi=hi[i,], se=se[i,], row.names=NULL)))
  }
  
  hicks_tbl_aligned    <- to_tbl(H_hat, H_lo, H_hi, H_se, shareNames, shareNames)
  hicks_tbl_native     <- to_tbl(H_hat0, H_lo0, H_hi0, H_se0, rn0, cn0)
  marshall_tbl_aligned <- to_tbl(M_hat, M_lo, M_hi, M_se, shareNames, shareNames)
  marshall_tbl_native  <- to_tbl(M_hat0, M_lo0, M_hi0, M_se0, rn0, cn0)
  
  marshall_diag_tbl <- tibble::tibble(
    good = shareNames, M_hat = Mdiag_hat, M_lo = Mdiag_lo, M_hi = Mdiag_hi, M_se = Mdiag_se
  )
  eta_tbl <- tibble::tibble(
    good = shareNames, eta_hat = eta_hat, eta_lo = eta_lo, eta_hi = eta_hi, eta_se = eta_se
  )
  hicks_diag_tbl <- tibble::tibble(
    good = shareNames,
    H_hat = diag(H_hat), H_lo = diag(H_lo), H_hi = diag(H_hi), H_se = diag(H_se)
  )
  
  note <- NULL
  if (all_miss) {
    note <- "Nomes não casaram: use *_native (ordem e rótulos nativos do E)."
    message("[Delta-elasticities v7] ", note)
  }
  
  list(
    ## Hicks
    hicks_full        = hicks_tbl_aligned,
    hicks_full_native = hicks_tbl_native,
    hicks_diag        = hicks_diag_tbl,
    ## Marshall (FULL + DIAG) — agora com ICs completos
    marshall_full        = marshall_tbl_aligned,
    marshall_full_native = marshall_tbl_native,
    marshall_diag        = marshall_diag_tbl,
    ## Expenditure
    eta              = eta_tbl,
    level            = level,
    method           = method,
    mapping_diag     = map$diag,
    note             = note
  )
}

cis_sys <- elasticity_CIs_delta_system_v7(
  fit_sys      = fit_sys,
  fit_template = best_fit,
  dat_sys      = dfSW[stats::complete.cases(dfSW[, use_all, drop=FALSE]), , drop=FALSE],
  priceNames   = priceNames,
  shareNames   = shareNames,
  drop_price   = drop_price,
  omit_share   = omit_share,
  level        = 0.95,
  method       = "bonferroni"  # ou "pointwise"
)

cat("\n[ICs Delta — Marshall]\n"); print(cis_sys$marshall_full_native, digits=3)
cat("\n[ICs Delta — Hicks]\n"); print(cis_sys$hicks_full_native, digits = 3)
cat("\n[ICs Delta — Eta]\n");             print(cis_sys$eta,          digits=3)

write.csv2(cis_sys$marshall_full_native, 'C:/Users/x16610962/Downloads/Table_C12_param_marshall_FULL_long.csv', row.names = FALSE)
write.csv2(cis_sys$hicks_full_native, 'C:/Users/x16610962/Downloads/Table_C12_param_hicks_FULL_long.csv', row.names = FALSE)
write.csv2(cis_sys$eta, 'C:/Users/x16610962/Downloads/Table_C12_param_eta.csv', row.names = FALSE)

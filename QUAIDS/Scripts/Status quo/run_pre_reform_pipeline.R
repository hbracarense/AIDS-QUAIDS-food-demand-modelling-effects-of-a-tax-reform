## =========================================================
## PIPELINE ÚNICO – PRÉ-REFORMA (ANEXO B)
## Requer: pre_reform_1.R ... pre_reform_17.R no working dir
## =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------
## 0. Caminhos e opções de saída
## ---------------------------------------------------------
out_dir_tables  <- "output_tables_B"
out_dir_figures <- "output_figures_B"
if (!dir.exists(out_dir_tables))  dir.create(out_dir_tables)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures)

## ---------------------------------------------------------
## 1. Estimação principal e diagnósticos “baseline”
##    (pre_reform_1.R)
## ---------------------------------------------------------
source("pre_reform_1.R", encoding = "UTF-8")
source("pre_reform_2.R", encoding = "UTF-8")
source("pre_reform_3.R", encoding = "UTF-8")
source("pre_reform_9.R",  encoding = "UTF-8")  # .j_by_eq_with_excl, fs_block_F (cluster), partial_R2_exclIV, wu_hausman_cf etc.
source("pre_reform_10.R", encoding = "UTF-8")  # utilidades adicionais de J por eq, Wu-Hausman, partial R2
source("pre_reform_11.R", encoding = "UTF-8")
source("pre_reform_12.R", encoding = "UTF-8")
source("pre_reform_13.R", encoding = "UTF-8")# post_refit_report_auto e relatórios eq-specific
source("pre_reform_14.R", encoding = "UTF-8")
source("pre_reform_15.R", encoding = "UTF-8")  # delta_quaids_CIs()
source("pre_reform_16.R", encoding = "UTF-8")  # boot_quaids_elasticities_byUF(), make_eq_justID2(), format_boot_tables()
source("pre_reform_17.R", encoding = "UTF-8") 
## Este script:
## - Lê o banco .dta
## - Constrói df_iv
## - Define fit_quaids_manual_km1()
## - Estima o modelo principal: fit_q (quaids_km1_fit)
## - Define diversas funções de diagnóstico usadas abaixo

stopifnot(exists("df_iv"), exists("fit_q"))

## =========================================================
## TABELA B.1 – First-stage strength (joint F-tests)
## =========================================================
## Fonte: pre_reform_1.R, função fs_block_F() e uso em:
##   fs_hc3  <- fs_block_F(fit_q, df_iv, robust = TRUE)   # recomendado para o paper
##   fs_homo <- fs_block_F(fit_q, df_iv, robust = FALSE)

fs_B1 <- fs_block_F(fit_q, df_iv)
## Ajuste de nomes para bater com o layout do anexo
tab_B1 <- fs_B1
colnames(tab_B1) <- c("Endogenous", "F", "df1", "df2", "p_value")
write.csv2(tab_B1, file.path(out_dir_tables, "B1_first_stage_baseline.csv"), row.names = FALSE)

## =========================================================
## TABELA B.2 – Exogeneity test (control function)
## =========================================================
## Fonte: pre_reform_1.R, função wu_por_eq_manual2()
## (estatística Chi2, df, p-valor por equação)

tab_B2 <- wu_por_eq_manual2(fit_q, df_iv, robust = TRUE, verbose = FALSE)
## Organização conforme texto do anexo
tab_B2_export <- tab_B2
write.csv2(tab_B2_export, file.path(out_dir_tables, "B2_exogeneity_control_function.csv"), row.names = FALSE)

## =========================================================
## TABELA B.3 – Over-identification tests per equation
## =========================================================
## Fonte: pre_reform_1.R, função diag_iv_basico()
## Esta função produz, entre outros, df_J por equação.

diag_baseline <- diag_iv_basico(fit_q, df_iv, inst_terms = fit_q$inst_terms)

## Se df_J == 0, equação é exatamente identificada (sem sobre-ID).
tab_B3 <- data.frame(
  Equation          = diag_baseline$eq,
  Over_ID_available = ifelse(diag_baseline$df_J > 0, "Yes", "No (exactly identified)"),
  J                 = ifelse(diag_baseline$df_J > 0, diag_baseline$J_hansen, NA_real_),
  df                = ifelse(diag_baseline$df_J > 0, diag_baseline$df_J, NA_integer_),
  p_value           = ifelse(diag_baseline$df_J > 0, diag_baseline$p_hansen, NA_real_)
)
write.csv2(tab_B3, file.path(out_dir_tables, "B3_over_id_per_equation.csv"), row.names = FALSE)

## =========================================================
## EXTRA – RESET (forma funcional): clássico e robusto (HC3)
## =========================================================

# checagens mínimas (objetivas)
stopifnot(exists("df_iv"), exists("ols_w3"))
if (!requireNamespace("lmtest", quietly = TRUE)) stop("Falta pacote: lmtest")
if (!requireNamespace("sandwich", quietly = TRUE)) stop("Falta pacote: sandwich")

# 1) RESET clássico (tipo fitted; power 2:3) — igual ao pre_reform_1.R
reset_cl2 <- lmtest::resettest(ols_w3, power = 2:3, type = "fitted")

# Extrai estatística e p-valor (htest)
reset_cl_F <- as.numeric(reset_cl2$statistic[1])
reset_cl_p <- as.numeric(reset_cl2$p.value)

# 2) RESET robusto (HC3): restrito vs irrestrito com yhat^2 e yhat^3
df_aux2 <- within(df_iv, {
  yhat  <- fitted(ols_w3)
  yhat2 <- yhat^2
  yhat3 <- yhat^3
})
f_unres2 <- update(stats::formula(ols_w3), . ~ . + yhat2 + yhat3)
ols_w3_unres2 <- stats::lm(f_unres2, data = df_aux2)

reset_rob2 <- lmtest::waldtest(
  ols_w3, ols_w3_unres2,
  vcov = function(x) sandwich::vcovHC(x, type = "HC3")
)

# Extrai F e p-valor do waldtest (última linha é a comparação relevante)
reset_rob_F <- as.numeric(reset_rob2[nrow(reset_rob2), "F"])
reset_rob_p <- as.numeric(reset_rob2[nrow(reset_rob2), "Pr(>F)"])

# Monta tabela (formato simples e auditável)
tab_RESET <- data.frame(
  test = c("RESET_classical_fitted_power2_3", "RESET_robust_HC3_restricted_vs_unrestricted"),
  statistic_F = c(reset_cl_F, reset_rob_F),
  p_value = c(reset_cl_p, reset_rob_p),
  row.names = NULL
)

write.csv2(tab_RESET, file.path(out_dir_tables, "B_RESET_functional_form.csv"), row.names = FALSE)


## =========================================================
## TABELA B.4 – GVIF e GVIF adjusted by predictor
## =========================================================
## Fonte: pre_reform_1.R, trecho com car::vif() e cálculo de GVIF_adj.
## Aqui reproduzimos o cálculo em forma programática.

if (!requireNamespace("car", quietly = TRUE)) stop("Precisa do pacote 'car' para GVIF.")

# ASSUNÇÃO: Figura baseada na regressão de uma equação representativa (ex.: eq2).
f_eq2 <- formula(fit_q$eq[[2]])
ols_eq2 <- lm(f_eq2, data = df_iv)
v <- car::vif(ols_eq2)
## car::vif retorna GVIF (ou GVIF^(1/(2*Df))) dependendo da classe;
## no anexo, há GVIF e GVIF_adj. Aqui seguimos o padrão usual de GVIF_adj:

make_gvif_table <- function(v) {
  if (is.matrix(v)) {
    df  <- v[, "Df"]
    g   <- v[, 1]
  } else {
    df  <- rep(1, length(v))
    g   <- as.numeric(v)
  }
  gvif_adj <- g^(1 / (2 * df))
  data.frame(
    Predictor = names(g),
    GVIF      = as.numeric(g),
    Df        = as.integer(df),
    GVIF_adj  = as.numeric(gvif_adj),
    row.names = NULL
  )
}

tab_B4 <- make_gvif_table(v)
write.csv2(tab_B4, file.path(out_dir_tables, "B4_GVIF_table.csv"), row.names = FALSE)

## =========================================================
## TABELA B.5 – Coefficients (3SLS) – equation w_exp2
## =========================================================
## =========================================================
## TABELA B.5 – Coefficients (3SLS) – equation w_exp2
## =========================================================
sf <- fit_q$fit

# LHS das equações
lhs_vec <- vapply(sf$eq, function(m) as.character(stats::formula(m)[[2]]), character(1))

# índice da equação-alvo (w_despesahat2)
idx2 <- which(lhs_vec == "w_despesahat2")
if (length(idx2) != 1) stop("Não achei w_despesahat2 no LHS: ", paste(lhs_vec, collapse=", "))

# prefixo dos coeficientes no summary(sf): remove "_" do LHS para bater com "wdespesahat2_..."
coef_prefix <- paste0(gsub("_", "", lhs_vec[idx2]), "_")   # "wdespesahat2_"

# matriz completa de coeficientes do systemfit
C <- summary(sf)$coefficients
rn <- rownames(C)

# filtrar apenas a equação alvo
sel <- startsWith(rn, coef_prefix)
if (!any(sel)) stop("Não encontrei coeficientes com prefixo ", coef_prefix, ". Ex.: ", paste(head(rn, 10), collapse=", "))

C2 <- C[sel, , drop = FALSE]

tab_B5 <- data.frame(
  Equation   = lhs_vec[idx2],
  Variable   = sub(coef_prefix, "", rownames(C2)),
  Estimate   = C2[, "Estimate"],
  Std.error  = C2[, "Std. Error"],
  t          = C2[, "t value"],
  p_value    = C2[, "Pr(>|t|)"],
  row.names  = NULL
)

write.csv2(tab_B5, file.path(out_dir_tables, "B5_coefficients_eq_w_despesahat2.csv"), row.names = FALSE)

## =========================================================
## TABELAS B.6 e B.7 – Instrument relevance & Endogeneity
## (baseline specification)
## =========================================================
## B.6 – ASSUNÇÃO: usa check_iv_strength(...) do pre_reform_1.R
##       no conjunto "core" (iv_set_core), conforme exemplo.
## B.7 – ASSUNÇÃO: usa versão clusterizada Wu-Hausman (wu_hausman_cf)
##       do baseline (cluster por psu/uf se disponível).

## B.6:
iv_strength_B6 <- check_iv_strength(
  data_Z     = df_iv,
  instNames  = iv_set_core,   # definido em pre_reform_1.R
  drop_price = 1,
  priceNames = fit_q$priceNames,
  use_z2     = TRUE,
  shareNames = fit_q$shareNames
)
write.csv2(iv_strength_B6, file.path(out_dir_tables, "B6_instrument_relevance_baseline.csv"), row.names = FALSE)

## B.7:
tab_B7 <- wu_hausman_cf(fit_q, df_iv, cluster_var = cluster_used)
write.csv2(tab_B7, file.path(out_dir_tables, "B7_endogeneity_baseline_cluster.csv"), row.names = FALSE)

## =========================================================
## TAB B.8 – Hansen-J p-values: baseline vs specific
## =========================================================
## Fonte: pre_reform_1.R, função diag_iv_basico()
## Baseline = fit_q
## Specific/final = fit_final (top_out$fit_final)
##
## Saída:
##  - B3_hansenJ_pvalues_baseline_vs_specific.csv  (tabela larga)
##  - B3_hansenJ_pvalues_long.csv                  (insumo para gráfico)
## =========================================================


stopifnot(exists("fit_q"), exists("fit_final"), exists("df_iv"))
stopifnot(exists("jtests_by_eq_manual"))

# cluster_var: usa se existir e se estiver no df
cluster_var <- NULL
if (exists("cluster_used") && !is.null(cluster_used) && cluster_used %in% names(df_iv)) {
  cluster_var <- cluster_used
}

# 1) Baseline (fit_q)
j_base <- jtests_by_eq_manual(
  fit_obj     = fit_q,
  data        = df_iv,
  cluster_var = cluster_var,
  verbose     = FALSE
)

# 2) Specific/final (fit_final)
j_spec <- jtests_by_eq_manual(
  fit_obj     = fit_final,
  data        = df_iv,
  cluster_var = cluster_var,
  verbose     = FALSE
)

# Checagem: colunas necessárias
need_cols <- c("eq","df_J","p_hansen","J_hansen")
if (!all(need_cols %in% names(j_base))) stop("j_base sem colunas esperadas: ", paste(setdiff(need_cols, names(j_base)), collapse=", "))
if (!all(need_cols %in% names(j_spec))) stop("j_spec sem colunas esperadas: ", paste(setdiff(need_cols, names(j_spec)), collapse=", "))

# Alinha por equação (sem renomear)
if (!setequal(j_base$eq, j_spec$eq)) {
  stop("Equações diferem entre baseline e specific.\n",
       "Baseline: ", paste(j_base$eq, collapse=", "), "\n",
       "Specific: ", paste(j_spec$eq, collapse=", "))
}
j_spec <- j_spec[match(j_base$eq, j_spec$eq), ]

# Tabela larga (para apêndice e auditoria)
tab_B3 <- data.frame(
  Equation          = j_base$eq,
  dfJ_baseline      = j_base$df_J,
  hansenJ_p_baseline= j_base$p_hansen,
  dfJ_specific      = j_spec$df_J,
  hansenJ_p_specific= j_spec$p_hansen,
  stringsAsFactors = FALSE
)

# Exporta
if (exists("out_dir_tables")) {
  write.csv2(tab_B3, file.path(out_dir_tables, "B8_hansenJ_pvalues_baseline_vs_specific.csv"), row.names = FALSE)
  
  # Versão “dfJ=1” no specific (igual ao título do seu apêndice)
  tab_B3_dfJ1 <- subset(tab_B3, dfJ_specific == 1)
  write.csv2(tab_B3_dfJ1, file.path(out_dir_tables, "B8_hansenJ_pvalues_specific_dfJ1_only.csv"), row.names = FALSE)
}

## =========================================================
## TABELA B.9 (ROUND 2) – Hansen J p-values: PSU vs UF clustering
## =========================================================
## Fonte: pre_reform_13.R
## - sweep2$fit é criado no ROUND 2 (eqspec_sweep_all)
## - compare_cluster_specs_auto() calcula jtests_overid_auto() por cluster
## Saída desejada:
## Equation | dfJ | p (PSU) | p (UF)

stopifnot(exists("sweep2"), exists("df_iv"))
stopifnot(exists("compare_cluster_specs_auto"))

# calcula J por equação sob diferentes clusters (psu/uf)
j_cl <- compare_cluster_specs_auto(
  fit_obj      = sweep2$fit,
  data         = df_iv,
  cluster_vars = c("psu", "uf")
)

# checagens mínimas de colunas esperadas do jtests_overid_auto
need_cols <- c("eq", "df_J", "p_hansen", "cluster")
miss <- setdiff(need_cols, names(j_cl))
if (length(miss)) {
  stop("compare_cluster_specs_auto não retornou colunas esperadas: faltando ",
       paste(miss, collapse = ", "),
       "\nCols disponíveis: ", paste(names(j_cl), collapse = ", "))
}

# monta tabela wide: p(PSU) e p(UF)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

tab_B4_round2_cluster <- j_cl %>%
  filter(cluster %in% c("psu", "uf")) %>%
  mutate(cluster = toupper(cluster)) %>%
  tidyr::pivot_wider(
    names_from  = cluster,
    values_from = p_hansen,
    names_prefix = "p_"
  ) %>%
  rename(`p (PSU)` = p_PSU, `p (UF)` = p_UF)


write.csv2(
  tab_B4_round2_cluster,
  file.path(out_dir_tables, "B9_hansenJ_psu_vs_uf_round2.csv"),
  row.names = FALSE
)


## ---------------------------------------------------------
## 2. Ferramentas TOP TIER e modelo “final specification”
## ---------------------------------------------------------
## Aqui entram os scripts com rotinas eq-specific, swaps, C-tests etc.

 # run_top_tier(), eqspec_sweep_all, c_tests_eqspec_blocks, iv_moment_diag_auto etc.

## Conjunto inicial de IVs excluídos para o TOP TIER
## (explicitamente especificado no pre_reform_14.R)
excl_now <- c("IV_d_uf_X11","IV_d_uf_X12","IV_d_uf_X13","IV_d_uf_X15",
              "iv_op01","iv_op02","iv_op03")

top_out <- run_top_tier(
  fit_q  = fit_q,
  df_iv  = df_iv,
  excl_now = excl_now,
  must_keep_global   = c("iv_op01","iv_op02","iv_op03"),
  cluster_preference = c("psu","uf", NA),
  alpha_flag         = 0.05,
  sweep_dfJ_max      = 1,
  do_round2_relax    = TRUE,
  round2_dfJ_max     = 2,
  do_ctests          = TRUE,
  enforce_just_id    = "if_reject"
)

fit_final   <- top_out$fit_final      # quaids_km1_fit final
tables_top  <- top_out$tables
cluster_used <- top_out$cluster_used

## =========================================================
## TABELAS B.10 e B.11 – Final specification
## =========================================================
## B.8 – Robust first-stage strength in final specification:
##       usa fs_block_F clusterizado no modelo final.
## B.9 – Endogeneity test in final specification:
##       usa wu_hausman_cf clusterizado no modelo final.

tab_B8 <- tables_top$F_block
write.csv2(tab_B8, file.path(out_dir_tables, "B10_first_stage_final_spec.csv"), row.names = FALSE)

tab_B9 <- tables_top$wu
write.csv2(tab_B9, file.path(out_dir_tables, "B11_endogeneity_final_spec.csv"), row.names = FALSE)

## =========================================================
## ELASTICIDADES – Tabelas B.12, B.13, B.14 e Figura B.1
## =========================================================
##  - B.12 (parametric simulation): delta_quaids_CIs() em pre_reform_15.R
##  - B.11 e B.13 (bootstrap por UF): boot_quaids_elasticities_byUF()
##    e format_boot_tables() em pre_reform_16.R, com sumarização em
##    pre_reform_17.R.
##  - Figura B.5 (heatmap de elasticidades): construída a partir
##    das tabelas de elasticidades.

 # sumarização dos 78 parâmetros em objetos tidy

## --- B.12: elasticidades via delta / simulação paramétrica ---
out_delta <- delta_quaids_CIs(fit_final, df_iv, level = 0.95)  # fit_final é quaids_km1_fit (wrapper com $fit)

tab_B12_M <- out_delta$marshallian
tab_B12_H <- out_delta$hicksian
tab_B12_E <- out_delta$expenditure

write.csv2(tab_B12_M, file.path(out_dir_tables, "B13_elasticities_parametric_marshallian.csv"), row.names = FALSE)
write.csv2(tab_B12_H, file.path(out_dir_tables, "B13_elasticities_parametric_hicksian.csv"),   row.names = FALSE)
write.csv2(tab_B12_E, file.path(out_dir_tables, "B14_elasticities_parametric_expenditure.csv"),row.names = FALSE)

## --- B.11 e B.13: bootstrap por UF (just-ID eq.6) ---
## Fonte: pre_reform_16.R (make_eq_justID2, boot_quaids_elasticities_byUF)
## ATENÇÃO: pre_reform_16.R já possui uma EXECUÇÃO no final,
## que constrói:
##   fit_just6  (quaids_km1_fit eq-specific just-ID na eq.6)
##   B         (lista com draws, base, etc.)
##   elas_bs   (tabelas já formatadas resumidas)
##
## Aqui apenas reaproveitamos esses objetos.

stopifnot(exists("fit_just6"), exists("B"), exists("elas_bs"))

## B.11 – conjunto completo (Marshall, Hicks, renda) com IC bootstrap:
write.csv2(elas_bs$marshallian, file.path(out_dir_tables, "B12_elasticities_bootstrap_marshallian.csv"), row.names = FALSE)
write.csv2(elas_bs$hicksian,    file.path(out_dir_tables, "B12_elasticities_bootstrap_hicksian.csv"),    row.names = FALSE)
write.csv2(elas_bs$expenditure, file.path(out_dir_tables, "B14_elasticities_bootstrap_expenditure.csv"), row.names = FALSE)

## B.13 – efeitos do bootstrap sobre ICs:
## ASSUNÇÃO: B.13 é derivada da comparação entre B.11 (bootstrap) e
##           B.12 (paramétrico), focando nos own-prices, como descrito
##           no texto. Esta agregação está **parcialmente** implementada
##           em pre_reform_17.R, que sumariza os 78 parâmetros.
##           Aqui, deixo o output de pre_reform_17.R como fonte direta
##           para montar B.13 (por você), pois o script já produz:
##           - marshall_full, hicks_full, yl_sum, etc.

## Objetos de pre_reform_17:
##   marshall_full, hicks_full, yl_sum, ...
## Você pode filtrá-los (own-prices) e resumir diferenças de largura
## de IC para montar B.13. Para não extrapolar além do código fornecido,
## não imponho aqui um formato adicional — apenas salvo como está:



## =========================================================
## TABELA B.14 – Significancy share of log-prices (CI 95%) by equation
## =========================================================
## Fonte: coeficientes e vcov dos systemfit (3SLS) dos 4 blocos (Spec 1..4)
## Critério de significância (IC 95%): lwr>0 OU upr<0  (IC não cruza 0)
## =========================================================

stopifnot(exists("out_dir_tables"))

## 1) Localizar a lista de fits dos 4 blocos (Spec 1..4)
##    IMPORTANTÍSSIMO: este patch NÃO inventa nomes.
##    Ele tenta achar um objeto-lista já existente; caso não ache, aborta com mensagem clara.

.pick_block_fits <- function() {
  # candidatos comuns (se existir no seu run_pipeline)
  cand_names <- c("fit_blocks", "fits_blocks", "blocks_fits", "block_fits", "fits_by_block")
  for (nm in cand_names) {
    if (exists(nm, inherits = TRUE)) {
      obj <- get(nm, inherits = TRUE)
      if (is.list(obj) && length(obj) >= 4) return(obj)
    }
  }
  NULL
}

block_fits <- .pick_block_fits()
if (is.null(block_fits)) {
  stop(
    "Não encontrei a lista com os 4 fits dos blocos (Spec 1..4). ",
    "Crie no pipeline um objeto-lista (ex.: fit_blocks[[1:4]]) ",
    "no mesmo ponto em que a Tabela B.1 é montada, e rode novamente."
  )
}

## Normaliza: queremos exatamente 4 specs, em ordem 1..4
block_fits <- block_fits[1:4]

## 2) Extrair systemfit de cada spec (aceita wrapper com $fit)
.get_sys <- function(fit_obj) {
  if (inherits(fit_obj, "systemfit")) return(fit_obj)
  if (is.list(fit_obj) && !is.null(fit_obj$fit) && inherits(fit_obj$fit, "systemfit")) return(fit_obj$fit)
  stop("Um dos blocos não é systemfit nem wrapper com $fit.")
}

sys_list <- lapply(block_fits, .get_sys)

## 3) Função para obter CI 95% (normal) para um termo em uma equação
##    Observação: usa vcov(systemfit) do próprio objeto (não força robust/cluster).
.get_ci_term <- function(sf, eq_lhs, term, level = 0.95) {
  # nomes de coef no seu projeto: prefixo = lhs sem "_" (ex.: w_despesahat2 -> wdespesahat2)
  px <- gsub("_", "", eq_lhs)
  nm <- paste0(px, "_", term)
  
  cf <- stats::coef(sf)
  V  <- stats::vcov(sf)
  
  if (!(nm %in% names(cf))) return(NULL)
  if (!(nm %in% rownames(V))) return(NULL)
  
  est <- unname(cf[[nm]])
  se  <- sqrt(V[nm, nm])
  z   <- qnorm(1 - (1 - level)/2)
  lwr <- est - z * se
  upr <- est + z * se
  
  list(est = est, se = se, lwr = lwr, upr = upr)
}

## 4) Montar base longa: spec_id x equação x preço
price_terms <- paste0("ln_preco_por_kg", 2:6)

rows <- list()
for (s in seq_along(sys_list)) {
  sf <- sys_list[[s]]
  
  # equações pelo LHS (forma segura, independente de nomes em sf$eq)
  lhs_vec <- vapply(sf$eq, function(m) as.character(stats::formula(m)[[2]]), character(1))
  
  for (eq_lhs in lhs_vec) {
    for (pt in price_terms) {
      ci <- .get_ci_term(sf, eq_lhs, pt, level = 0.95)
      if (is.null(ci)) next
      
      sig <- as.integer(ci$lwr > 0 | ci$upr < 0)
      
      rows[[length(rows) + 1L]] <- data.frame(
        spec_id    = s,
        equation   = eq_lhs,
        variable   = pt,
        lwr        = ci$lwr,
        upr        = ci$upr,
        significant = sig,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_long <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (!nrow(df_long)) stop("Não consegui extrair nenhum coeficiente de ln_preco_por_kg2..6 dos 4 blocos.")

## 5) Agregações “By specification” e “Overall”
suppressPackageStartupMessages({
  library(dplyr)
})

tab_by_spec <- df_long %>%
  group_by(spec_id, variable) %>%
  summarise(
    Significant = sum(significant, na.rm = TRUE),
    Total = n(),
    `Share significant` = Significant / Total,
    .groups = "drop"
  ) %>%
  arrange(spec_id, variable)

tab_overall <- df_long %>%
  group_by(variable) %>%
  summarise(
    Significant = sum(significant, na.rm = TRUE),
    Total = n(),
    `Share significant` = Significant / Total,
    .groups = "drop"
  ) %>%
  arrange(variable)

## 6) Exportar
write.csv2(tab_overall, file.path(out_dir_tables, "B14_signif_share_logprices_overall.csv"), row.names = FALSE)
write.csv2(tab_by_spec, file.path(out_dir_tables, "B14_signif_share_logprices_by_spec.csv"), row.names = FALSE)

## =========================================================
## FIGURE B.1 – Distribution of t-statistics for log-price coefficients
## =========================================================

sf <- fit_q$fit   # objeto systemfit (baseline), vindo de pre_reform_1.R

get_t_ln <- function(mod) {
  sm <- summary(mod)
  coefs <- sm$coefficients   # matriz: Estimate, Std. Error, t value, Pr(>|t|)
  rn <- rownames(coefs)
  idx <- grepl("^ln_", rn)
  if (!any(idx)) return(numeric(0))
  coefs[idx, "t value"]
}

t_list <- lapply(sf$eq, get_t_ln)
t_all  <- unlist(t_list)

df_B7 <- data.frame(t_value = t_all)

p_B7 <- ggplot(df_B7, aes(x = t_value)) +
  geom_histogram(bins = 30, aes(y = ..density..)) +
  geom_density() +
  xlab("t-statistics for log-price coefficients") +
  ylab("Density") 

ggsave(
  file.path(out_dir_figures, "Figure_B1_tstats_log_prices.png"),
  p_B7, dpi = 300, width = 7, height = 5
)

## --- Figura B.2 – Heatmap de elasticidades (ASSUNÇÃO) ---
## ASSUNÇÃO: a figura usa as elasticidades marshallianas “pontuais”
##           (paramétricas ou bootstrap). Aqui uso as paramétricas
##           (out_delta$marshallian) por serem diretamente associadas
##           ao “modelo final”.

if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)

heat_B5 <- elas_bs$marshallian   # marshallian: colunas i, j, est, se, lwr, upr
## ASSUMÇÃO: colunas i e j ou good_i e good_j estão presentes;
## se não, adapte conforme a estrutura do seu out_delta.

if (!"good_i" %in% names(heat_B5) && "i" %in% names(heat_B5)) {
  heat_B5$good_i <- as.factor(heat_B5$i)
}
if (!"good_j" %in% names(heat_B5) && "j" %in% names(heat_B5)) {
  heat_B5$good_j <- as.factor(heat_B5$j)
}

heat_B5$sig <- ifelse(heat_B5$lwr<=0 & 0<=heat_B5$upr, NA, heat_B5$est)

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

ggsave(file.path(out_dir_figures, "Figure_B2_heatmap_elasticities.png"), p_B5, dpi = 300, width = 7, height = 5)
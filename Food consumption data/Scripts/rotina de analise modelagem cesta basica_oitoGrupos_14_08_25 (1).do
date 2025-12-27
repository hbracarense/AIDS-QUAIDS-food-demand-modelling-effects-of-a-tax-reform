*construção da rotina para preços de alimentos para modelagem de cesta basica

use  "C:\Users\GEIPS1\OneDrive\INQUERITOS\POF\POF 2017_2018\Bancos_stata\caderneta coletiva.dta", clear 

	*criar identificação do domicílio
	gen double vident= ((v04*100)+ v05) 
	la var vident "identificacao domicilio"
	format vident %16.0f

	*criar identificação da unidade de consumo
	gen double ident= (((v04*100)+ v05) *10+v06)
	la var ident "identificacao da unidade de consumo"
	format ident %16.0f

gen cod_1718=int(v09/100)

save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\Banco_Caderneta_TCaldeira_2025.dta", replace

************************************************************
use "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\Banco_Caderneta_TCaldeira_2025.dta", clear
* 1. Selecionar os  alimentos com alíquotas conhecidas
**************juntar com dados de aliquotas
merge m:1 cod_1718 using "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_alimentos_aliquotas.dta"

drop if _merge ==1
drop if _merge ==2
	drop _merge
	
*****juntar com composição
	
merge m:1 cod_1718 using "C:\Users\GEIPS1\OneDrive\Doutorado\Banco de dados POF-IPCA\banco_completo_aliemntos_POF17_18_composicao.dta"

drop if _merge ==1
drop if _merge ==2
	drop _merge
	

merge m:1 vident using "C:\Users\GEIPS1\OneDrive\INQUERITOS\POF\POF 2017_2018\Bancos_stata\uf_area_pond.dta"

drop if _merge ==1
drop if _merge ==2
	drop _merge
	
************************************************************************
*v13 valor da despesa de um domicilio 
*v22 quantidade adiquirida por kg
summarize v13  v22
drop if v22 <= 0
drop if v13 <= 0
*gasto diario
*No modelo AIDS, a "despesa" (total expenditure) representa o total gasto monetário com os bens incluídos na análise, ou seja, o valor total gasto com os grupos alimentares que você está modelando.
gen despesa_diaria = v13 / 7
gen despesa = v13

*Quantidade diária em kg: 
gen quantidade_diaria_kg = v22 / 7

*Preço por kg: 
gen preco_por_kg = v13 / v22

*total de calorias
gen qgrama_dia = (v22*1000)/7
gen qliqg_dia = qgrama_dia * fc
gen calreal = (energia_kcal_media*qliqg_dia)/100


*por IC 95%
*usar despesa ou despesa_diaria
bysort cod_1718 : summarize despesa , detail

bysort cod_1718: egen media = mean(despesa)
bysort cod_1718: egen dp = sd(despesa)

gen lower_bound = media - 1.96 * dp
gen upper_bound = media + 1.96 * dp

gen dentro_95ci = (despesa >= lower_bound & despesa <= upper_bound)

list cod_1718 despesa if dentro_95ci == 0  // Verificar os outliers
drop if dentro_95ci == 0  // Remover os outliers

drop media dp lower_bound upper_bound dentro_95ci


*************************************

***********************
******teste: usar as médias de preços de subgrupos de alimentos ao inves dos alimentos em si para os calculos de retirada dos impostos


*subgrupos: 11 Leguminosas; 12 Arroz; 13 Farinha de trigo; 14 Féculas; 15 macarrão; 16 outras farinhas; 17 outros cereais; 18 pão frances; 19 reizes e tuberculos; 20 conserva legumes e verduras; 21 FLV secas/desidratadas; 22 hortaliças folhosas; 23 hortaliças frutosas; 24 doce a base de frutas; 25 frutas; 27 nozes e sementes; 28 aves; 29 boi; 30 frutos do mar; 31 ovos; 32 peixes; 33 porco; 34 viceras; 35 leite; 36 queijos processados; 37 açucar; 38 gordura animal; 39 gordura e oleo vegetal; 40 margarina; 41 outros açucares; 42 sal; 43 chás; 44 carnes salgadas; 45 peixes secos/ salgados; 48 condimentos; 49 bebidas lacteas; 50 biscoitos; 51 carnes reconstituidas (AUP); 52 cerais matinais; 53 chocolates, balas; 54 frutas AUP; 55 macarrão instantaneo; 56 pao doce, bolos;  57 produtos lacetos; 58 queijos AUP; 59 refeição pronta; 60 refrigerante; 61 salgadinhos; 62 sopas e misturas AUP; 63 sorvetes; 64 sucos artificiais

*cesta basica: 1 feijões (leguminosas); 2 cereais; 3 raizes; 4 legumes_verduras; 5 frutas; 6 oleaginosas; 7 carnes_ovos; 8 leites_queijos; 9 ACUCAR_SAL_OLEOS; 10 CAFE_CHA_MATE; 12 outros G2 (processados); 13 outrosG3 (ingredientes); 14 AUP; 



save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_completo_para_analises.dta", replace

use "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_completo_para_analises.dta", clear
****2. Calcular o preço sem impostos

***calculo de preços sem impostos- 
* Impostos totais (soma simples — ICMS por dentro, os demais por fora)

gen aqli_icms   = ICMS_Médio / 100
gen aliq_pis    = PIS_PASEP / 100
gen ali_cofins = COFINS / 100
gen ali_ipi    = IPI / 100

/*calcula a fração do preço final que corresponde ao total de impostos indiretos (ICMS, PIS, Cofins, IPI). Ou seja, ela responde à pergunta:

"Qual a fração do preço pago que é só imposto?"*/

gen imposto_total = (1 - ((1 / (1 + aqli_icms)) * (1 - (aliq_pis + ali_cofins + ali_ipi))))

* Preço líquido (sem impostos)
gen preco_liquido = preco_por_kg * (1 - imposto_total)

* Alternativamente, passo a passo:
gen preco_sem_icms = preco_por_kg / (1 + aqli_icms)
gen preco_liquido_v2 = preco_sem_icms / (1 + (aliq_pis + ali_cofins + ali_ipi))

mean preco_por_kg if cod_1718==82035
mean preco_liquido if cod_1718==82035
mean preco_liquido_v2 if cod_1718 == 82035

*tese se preco_liquido e preco_liquido_v2 batem
gen dif = preco_liquido - preco_liquido_v2
sum dif

*usar preco_liquido_v2
********calculo de preços com impostos reforma tributaria
*modelo 1
* Defina as alíquotas (em decimal) modelo 1
gen aliq_cbs1 = CBS_modelo1 / 100 // CBS
gen aliq_ibs1 = IBS_modelo1 / 100 // IBS 
gen aliq_is  = imposto_seletivo  / 100 // imposto seletivo (para refrigerantes)

   * 2. Calcular o fator de multiplicação com os novos tributos
gen carga_tributaria_nova1 = 1 + aliq_cbs1 + aliq_ibs1 + aliq_is

* 3. Calcular o preço final com a nova carga tributária
gen preco_com_reforma1 = preco_liquido_v2 * carga_tributaria_nova1


*calcular a diferença de preço entre o atual e com reforma, adicione:
gen dif_preco_reforma1 = preco_com_reforma1 - preco_liquido_v2
gen perc_variacao_preco1 = 100 * (dif_preco_reforma1 / preco_liquido_v2)
mean perc_variacao_preco1

*******************************
*modelo 2

gen IS_modelo2=  imposto_seletivo
replace IS_modelo2 =25 if cesta_basicaID==14

* Defina as alíquotas (em decimal) modelo 2
gen aliq_cbs2 = CBS_modelo2 / 100 // CBS
gen aliq_ibs2 = IBS_modelo2 / 100 // IBS 
gen aliq_is2 = IS_modelo2 / 100 //  IS 

   * 2. Calcular o fator de multiplicação com os novos tributos_ mais se reforma colocasse o IS para todos AUP
gen carga_tributaria_nova2 = 1 + aliq_cbs2 + aliq_ibs2 + aliq_is2

* 3. Calcular o preço final com a nova carga tributária
gen preco_com_reforma2 = preco_liquido_v2 * carga_tributaria_nova2

*calcular a diferença de preço entre o atual e com reforma, adicione:
gen dif_preco_reforma2 = preco_com_reforma2 - preco_liquido_v2
gen perc_variacao_preco2 = 100 * (dif_preco_reforma2 / preco_liquido_v2)
mean perc_variacao_preco2



*usar preco_com_reforma1 preco_com_reforma2
save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_calculo_aliquotas_21_07_25.dta", replace

use "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_calculo_aliquotas_21_07_25.dta", clear
*agrupar grupos de alimentos de acordo o banco
*cesta basica: 1 feijões (leguminosas); 2 cereais; 3 raizes; 4 legumes_verduras; 5 frutas; 6 oleaginosas; 7 carnes_ovos; 8 leites_queijos; 9 ACUCAR_SAL_OLEOS; 10 CAFE_CHA_MATE; 12 outros G2 (processados); 13 outrosG3 (ingredientes); 14 AUP; 
mean preco_por_kg  [pw=peso_final], over(cesta_basicaID)
sum preco_por_kg, detail

mean preco_com_reforma1   [pw=peso_final], over(cesta_basicaID)
mean preco_com_reforma2  [pw=peso_final], over(cesta_basicaID)
*******************************************************************
/* Criando a variável class_nova_subgrupo com a nova categorização */
drop subgrupos
rename class_nova_subgrupo subgrupos
gen class_nova_subgrupo = .

/* Grupo 1: Básicos (Cereais, raízes, tubérculos) pão frances, cafe, cha*/
replace class_nova_subgrupo = 1 if subgrupos >= 11 & subgrupos <= 17 | subgrupos == 19 | subgrupos == 52| subgrupos == 18  | subgrupos == 43

/* Grupo 2: Frutas/Verduras/Leguminosas */
replace class_nova_subgrupo = 2 if subgrupos == 20 | subgrupos == 21 | subgrupos >= 22 & subgrupos <= 25 | subgrupos == 27

/* Grupo 3: Carnes/Peixes/Ovos */
replace class_nova_subgrupo = 3 if subgrupos >= 28 & subgrupos <= 34 | subgrupos == 44 | subgrupos == 45

/* Grupo 4: Laticínios */
replace class_nova_subgrupo = 4 if subgrupos == 35 | subgrupos == 36 | subgrupos == 46

/* Grupo 5: Óleos/Gorduras/Açúcares/Sal */
replace class_nova_subgrupo = 5 if subgrupos >= 37 & subgrupos <= 42 | subgrupos == 48

/* Grupo 6: Alimentos Ultraprocessados (AUP) gerais */
replace class_nova_subgrupo = 6 if subgrupos == 51 | subgrupos == 54 | subgrupos == 55 | subgrupos == 57 | subgrupos == 58 | subgrupos == 59 | subgrupos == 61 | subgrupos == 62 | subgrupos == 63

/* Grupo 7: Bebidas Adoçadas */
replace class_nova_subgrupo = 7 if subgrupos == 49 | subgrupos == 60 | subgrupos == 64

/* Grupo 8: Panificados e Confeitados */
replace class_nova_subgrupo = 8 if  subgrupos == 50 | subgrupos == 53 | subgrupos == 56

tab class_nova_subgrupo
 drop if class_nova_subgrupo==.
 
collapse (sum)  despesa preco_por_kg  preco_com_reforma1  preco_com_reforma2 qgrama_dia  calreal , by(class_nova_subgrupo vident)


compress

reshape wide despesa preco_por_kg  preco_com_reforma1  preco_com_reforma2 qgrama_dia  calreal , i(vident) j(class_nova_subgrupo)

**********************************************************************

save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_reshape_vident.dta", replace
use "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_reshape_vident.dta", clear
*Juntando com o banco de variáveis sociodemográficas
joinby vident using "C:\Users\GEIPS1\OneDrive\Doutorado\Otimização linear\Analises Thaís\bancos\renda_regiao_area_ponderacao.dta", _merge(_merge)  unmatched(both)


*estrato uf area regiao nfam vident renda_total v49 v50 v04 v05 id_capital       

joinby vident using "C:\Users\GEIPS1\OneDrive\INQUERITOS\POF\Análises Bruna\caract_sociodemográficas.dta", _merge(_merge2)  unmatched(both)

mvencode _all, mv(0) override
drop _merge _merge2
*area 1 – Urbano; 2 – Rural 


rename  v49 peso
rename  v50 peso_final /* fexpdom  */

drop regiao
gen regiao=.
*norte
replace regiao =1  if  (uf == 11 | uf ==  12  | uf == 13  | uf == 14  | uf == 15  | uf == 16  | uf == 17)
*Nordeste
replace regiao =2  if (uf == 21 | uf ==22 | uf ==23 | uf ==24 | uf ==25 | uf ==26 | uf ==27 | uf ==28 | uf ==29)
*Sudeste
replace regiao =3  if (uf == 31 | uf == 32 | uf == 33 | uf == 35) 
*Sul
replace regiao =4  if (uf == 41 | uf == 42 | uf == 43) 
*Centro-oeste
replace regiao =5  if (uf == 50 | uf == 51  | uf == 52 | uf ==53)


*****************
*imputação de dados de preços faltantes

*variavel capital e região metropolitana ou interior
gen capital_interior = .

* Capital
replace capital_interior = 1 if inlist(estrato, 1101, 1102, 1201, 1301, 1302, 1303, 1304, 1305, 1306, 1401, 1402, 1501, 1502, 1503, 1601, 1602, 1701, ///
                                      2101, 2102, 2103, 2201, 2202, 2203, 2301, 2302, 2303, 2304, 2305, 2306, 2401, 2402, 2501, 2502, 2503, 2601, 2602, 2603, ///
                                      2701, 2702, 2703, 2801, 2802, 2901, 2902, 2903, 2904, 2905, 2906, ///
                                      3101, 3102, 3103, 3104, 3105, 3106, 3201, 3202, 3301, 3302, 3303, 3304, 3305, 3306, 3307, 3308, 3309, 3501, 3502, 3503, 3504, 3505, 3506, 3507, 3508, 3509, ///
                                      4101, 4102, 4103, 4104, 4105, 4201, 4202, 4301, 4302, 4303, 4304, 4305, 4306, ///
                                      5001, 5002, 5003, 5101, 5102, 5201, 5202, 5203, 5301, 5302, 5303, 5304, 5305, 5306)

* Resto da Região Metropolitana (RM)
replace capital_interior = 2 if inlist(estrato, 1307, 1504, 1505, 1603, 2104, 2307, 2308, 2309, 2403, 2504, 2505, 2604, 2605, 2606, 2704, 2803, 2907, 2908, 2909, ///
                                      3107, 3108, 3109, 3203, 3204, 3205, 3310, 3311, 3312, 3313, 3314, 3315, 3316, 3317, 3318, 3510, 3511, 3512, 3513, 3514, 3515, ///
                                      4106, 4107, 4108, 4203, 4204, 4307, 4308, 4309, 5103, 5204, 5205, 5206)

* Interior (Resto da UF e Rural)
replace capital_interior = 3 if missing(capital_interior)
recode capital_interior 2=1
recode capital_interior 3=2
tab capital_interior


********olhar se está rodando por dois grupos guia ou por grupos cesta basica. selecionar qual vai usar
collapse (sum) despesa1 -calreal8 nfam n_mais65anos n_5anos  n_8esc n_0esc n_12esc sexofem renda_total (mean) anos_estudo_dom uf peso peso_final area regiao  capital_interior [pw=peso_final], by(estrato)

*todos apresentam 1 estrato sem preço 
****************************************************

gen renda_pc_est = renda_total / nfam
mean renda_pc_est [pw=peso_final]
xtile quintil_rendapc = renda_pc_est, n(5) // 
xtile tercil_rendapc = renda_pc_est, n(3) // 
xtile duo_rendapc = renda_pc_est, n(2) // 
mean renda_pc_est [pw=peso_final], over(quintil_rendapc)
*Transformando as variáveis em percentual

foreach x of varlist n_mais65anos-sexofem {
	gen p_`x' =(`x'*100)/nfam
	
} 

mean p_n_mais65anos [pw=peso_final]
mean p_n_5anos[pw=peso_final]

mean p_n_0esc [pw=peso_final]
mean p_n_8esc [pw=peso_final]
mean p_n_12esc [pw=peso_final]

mean p_sexofem[pw=peso_final]

prop regiao [pw=peso_final]

prop capital_interior [pw=peso_final]
  
**************************************************************
*teste de imputação
*imputar as medias de preços
egen media_preco3 = mean(despesa3), by(quintil_renda regiao area capital_interior)

*imputar valores ausentes
gen preco3_imp = despesa3
replace preco3_imp = media_preco3 if preco3_imp == 0 |preco3_imp ==.


egen media_preco3_nacional = mean(despesa3) if despesa3 > 0 | preco3_imp ==.
replace media_preco3 = media_preco3_nacional if missing(media_preco3)

replace despesa3 = media_preco3 if despesa3 == 0 |missing(despesa3) 
*verificar os valores imputados

list despesa3 quintil_renda regiao area capital_interior ///
    if despesa3 <= 0 | missing(despesa3)

	*media para os valores validos
	egen media_preco3_ = mean(despesa3) if despesa3 > 0, ///
    by(quintil_renda regiao area capital_interior)
	*Substitua todos os valores ausentes ou ≤ 0 pela média imputada
	gen despesa3imp = despesa3
replace despesa3imp = media_preco3_ if missing(despesa3) | despesa3 <= 0
*flag para imputação
gen flag_preco3_imp = (missing(despesa3) | despesa3<= 0)

replace despesa3 = despesa3imp if despesa3 == 0 |missing(despesa3) 
*verificar os valores imputados

list despesa3 quintil_renda regiao area capital_interior ///
    if despesa3 <= 0 | missing(despesa3)
	edit despesa3
	**********************************************
*imputar as medias de preços
egen media_preco4 = mean(despesa4), by(quintil_renda regiao area capital_interior)

*imputar valores ausentes
gen preco4_imp = despesa4
replace preco4_imp = media_preco4 if preco4_imp == 0 |preco4_imp ==.


egen media_preco4_nacional = mean(despesa4) if despesa4 > 0 | preco4_imp ==.
replace media_preco4 = media_preco4_nacional if missing(media_preco4)

replace despesa4 = media_preco4 if despesa4 == 0 |missing(despesa4) 
*verificar os valores imputados

list despesa4 quintil_renda regiao area capital_interior ///
    if despesa4 <= 0 | missing(despesa4)

	*media para os valores validos
	egen media_preco4_ = mean(despesa4) if despesa4 > 0, ///
    by(quintil_renda regiao area capital_interior)
	*Substitua todos os valores ausentes ou ≤ 0 pela média imputada
	gen despesa4imp = despesa4
replace despesa4imp = media_preco4_ if missing(despesa4) | despesa4 <= 0
*flag para imputação
gen flag_preco4_imp = (missing(despesa4) | despesa4<= 0)

replace despesa4 = despesa4imp if despesa4 == 0 | despesa4 ==.
*verificar os valores imputados

list despesa4 quintil_renda regiao area capital_interior ///
    if despesa4 <= 0 | missing(despesa4)
	edit despesa4
	*************************************************************
	
	*imputar as medias de preços

egen media_preco6 = mean(despesa6), by(quintil_renda regiao area capital_interior)

*imputar valores ausentes
gen preco6_imp = despesa6
replace preco6_imp = media_preco6 if preco6_imp == 0 |preco6_imp ==.


egen media_preco6_nacional = mean(despesa6) if despesa6 > 0 | preco6_imp ==.
replace media_preco6 = media_preco6_nacional if missing(media_preco6)

replace despesa6 = media_preco6 if despesa6 == 0 |missing(despesa6) 
*verificar os valores imputados

list despesa6 quintil_renda regiao area capital_interior ///
    if despesa6 <= 0 | missing(despesa6)

	*media para os valores validos
	egen media_preco6_ = mean(despesa6) if despesa6 > 0, ///
    by(quintil_renda regiao area capital_interior)
	*Substitua todos os valores ausentes ou ≤ 0 pela média imputada
	gen despesa6imp = despesa6
replace despesa6imp = media_preco6_ if missing(despesa6) | despesa6 <= 0
*flag para imputação
gen flag_preco6_imp = (missing(despesa6) | despesa6<= 0)

replace despesa6 = despesa6imp if despesa6 == 0 |missing(despesa6) 
*verificar os valores imputados
list despesa6 quintil_renda regiao area capital_interior ///
    if despesa6 <= 0 | missing(despesa6)
	edit despesa6
	
	
*imputar as medias de preços

egen media_preco7 = mean(despesa7), by(quintil_renda regiao area capital_interior)

*imputar valores ausentes
gen preco7_imp = despesa7
replace preco7_imp = media_preco7 if preco7_imp == 0 |preco7_imp ==.


egen media_preco7_nacional = mean(despesa7) if despesa7 > 0 | preco7_imp ==.
replace media_preco7 = media_preco7_nacional if missing(media_preco7)

replace despesa7 = media_preco7 if despesa7 == 0 |missing(despesa7) 
*verificar os valores imputados

list despesa7 quintil_renda regiao area capital_interior ///
    if despesa7 <= 0 | missing(despesa7)

	*media para os valores validos
	egen media_preco7_ = mean(despesa7) if despesa7 > 0, ///
    by(quintil_renda regiao area capital_interior)
	*Substitua todos os valores ausentes ou ≤ 0 pela média imputada
	gen despesa7imp = despesa7
replace despesa7imp = media_preco7_ if missing(despesa7) | despesa7 <= 0
*flag para imputação
gen flag_preco7_imp = (missing(despesa7) | despesa7<= 0)

replace despesa7 = despesa7imp if despesa7 == 0 |missing(despesa7) 
*verificar os valores imputados
list despesa7 quintil_renda regiao area capital_interior ///
    if despesa7 <= 0 | missing(despesa7)
	edit despesa7
*/
*************************************************************

*Contrução do preço de mercado- preço atual
summ preco_por_kg1 preco_por_kg2 preco_por_kg3 preco_por_kg4 ///
                     preco_por_kg5 preco_por_kg6 preco_por_kg7 preco_por_kg8, detail
					 
					 
foreach var of varlist preco_por_kg1 preco_por_kg2 preco_por_kg3 preco_por_kg4 preco_por_kg5 preco_por_kg6 preco_por_kg7 preco_por_kg8 {
    
    * 1. Tratar zeros ou negativos
    replace `var' = . if `var' <= 0

    * 2. Criar log do preço
    gen ln_`var' = log(`var')

    * 3. Calcular média do log do preço
    quietly sum ln_`var', meanonly
    local media_log = r(mean)

    * 4. Rodar regressão de ln(preço) nos controles
    regress ln_`var' area nfam regiao n_mais65anos n_5anos ///
        sexofem uf anos_estudo_dom

    * 5. Obter resíduos
    predict resid_`var', resid

    * 6. Preço corrigido no log
    gen ln_corr_`var' = `media_log' + resid_`var'

    * 7. Voltar para nível
    gen preco_corr_`var' = exp(ln_corr_`var')
}
	
	
	
* Agora você terá preco_corrigido_preco_por_kg1 ... preco_corrigido_preco_por_kg8
summ preco_corr_preco_por_kg1 preco_corr_preco_por_kg2 preco_corr_preco_por_kg3 preco_corr_preco_por_kg4 preco_corr_preco_por_kg5 ///
preco_corr_preco_por_kg6 preco_corr_preco_por_kg7 preco_corr_preco_por_kg8, detail
****************************************
*****************************************
*Contrução do preço de mercado modelo 1
summ preco_com_reforma11 preco_com_reforma12 preco_com_reforma13 preco_com_reforma14 ///
                     preco_com_reforma15 preco_com_reforma16 preco_com_reforma17 preco_com_reforma18, detail
					 
					 
foreach var of varlist preco_com_reforma11 preco_com_reforma12 preco_com_reforma13 preco_com_reforma14  preco_com_reforma15 preco_com_reforma16 preco_com_reforma17 preco_com_reforma18 {
    
    * 1. Tratar zeros ou negativos
    replace `var' = . if `var' <= 0

    * 2. Criar log do preço
    gen ln_`var' = log(`var')

    * 3. Calcular média do log do preço
    quietly sum ln_`var', meanonly
    local media_log = r(mean)

    * 4. Rodar regressão de ln(preço) nos controles
    regress ln_`var' area nfam regiao n_mais65anos n_5anos ///
        sexofem uf anos_estudo_dom

    * 5. Obter resíduos
    predict resid_`var', resid

    * 6. Preço corrigido no log
    gen ln_corr_`var' = `media_log' + resid_`var'

    * 7. Voltar para nível
    gen preco_corr_`var' = exp(ln_corr_`var')
}
	
	
	
* Agora você terá preco_corrigido modelo 1
summ preco_corr_preco_com_reforma11 preco_corr_preco_com_reforma12 preco_corr_preco_com_reforma13 preco_corr_preco_com_reforma14 preco_corr_preco_com_reforma15 preco_corr_preco_com_reforma16 preco_corr_preco_com_reforma17 preco_corr_preco_com_reforma18, detail
****************************************
*****************************************
*Contrução do preço de mercado modelo 2
summ preco_com_reforma21 preco_com_reforma22 preco_com_reforma23 preco_com_reforma24 ///
                     preco_com_reforma25 preco_com_reforma26 preco_com_reforma27 preco_com_reforma28, detail
					 
					 
foreach var of varlist preco_com_reforma21 preco_com_reforma22 preco_com_reforma23 preco_com_reforma24 ///
                     preco_com_reforma25 preco_com_reforma26 preco_com_reforma27 preco_com_reforma28 {
    
    * 1. Tratar zeros ou negativos
    replace `var' = . if `var' <= 0

    * 2. Criar log do preço
    gen ln_`var' = log(`var')

    * 3. Calcular média do log do preço
    quietly sum ln_`var', meanonly
    local media_log = r(mean)

    * 4. Rodar regressão de ln(preço) nos controles
    regress ln_`var' area nfam regiao n_mais65anos n_5anos ///
        sexofem uf anos_estudo_dom

    * 5. Obter resíduos
    predict resid_`var', resid

    * 6. Preço corrigido no log
    gen ln_corr_`var' = `media_log' + resid_`var'

    * 7. Voltar para nível
    gen preco_corr_`var' = exp(ln_corr_`var')
}
	
	
	
* Agora você terá preco_corrigido modelo 2
summ preco_corr_preco_com_reforma21 preco_corr_preco_com_reforma22 preco_corr_preco_com_reforma23 preco_corr_preco_com_reforma24 preco_corr_preco_com_reforma25 preco_corr_preco_com_reforma26 preco_corr_preco_com_reforma27 preco_corr_preco_com_reforma28, detail


*A. Majumder, R. Ray, K. Sinha Calculating rural-urban food price differentials from unit values in household expenditure ///
	*surveys: a comparison with existing methods and a new procedure- uso do log do preço

*correção da endogeneidade de preços

/*preços não são fixos nem determinados externamente, mas são influenciados pelas próprias decisões de consumo das famílias ou por fatores não observados que afetam simultaneamente consumo e preço. 
a quantidade demandada afeta o preço (lei da oferta e demanda), ou porque há características não observadas (como preferências ou choques regionais) que influenciam simultaneamente o preço pago e a quantidade consumida.
Cox, T. L., & Wohlgenant, M. K. (1986). Prices and quality effects in cross-sectional demand analysis. American Journal of Agricultural Economics, 68(4), 908–919. */
**

****************************
****intermediario
gen lndespesa1 = log(despesa1)
reg lndespesa1 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos     sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat1
gen despesa1_hat = exp(ln_hat1)

gen lndespesa2 = log(despesa2)
reg lndespesa2 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos   sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat2
gen despesa2_hat = exp(ln_hat2)

gen lndespesa3 = log(despesa3)
reg lndespesa3 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos     sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat3
gen despesa3_hat = exp(ln_hat3)

gen lndespesa4 = log(despesa4)
reg lndespesa4 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos  sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat4
gen despesa4_hat = exp(ln_hat4)


gen lndespesa5 = log(despesa5)
reg lndespesa5 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos  sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat5
gen despesa5_hat = exp(ln_hat5)


gen lndespesa6 = log(despesa6)
reg lndespesa6 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos  sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat6
gen despesa6_hat = exp(ln_hat6)


gen lndespesa7 = log(despesa7)
reg lndespesa7 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos  sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat7
gen despesa7_hat = exp(ln_hat7)


gen lndespesa8 = log(despesa8)
reg lndespesa8 quintil_rendapc  area  nfam  regiao n_mais65anos n_5anos  sexofem p_n_0esc p_n_8esc  p_n_12esc
predict ln_hat8
gen despesa8_hat = exp(ln_hat8)

summ despesa1_hat despesa2_hat despesa3_hat despesa4_hat despesa5_hat despesa6_hat despesa7_hat despesa8_hat, detail

*usar o resultado no AIDS 

*****************************************************

save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_reshape_vident_pond.dta", replace


 *3. Agregar por domicílio e calcular participação no gasto
*************Avaliar o resultado com e sem impostos
use "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_reshape_vident_pond.dta",clear

* Participação atual no gasto 
*Correção da endogeneidade do dispêndio total
/*
o gasto total com alimentos é geralmente incluído como variável explicativa — porque a demanda por um bem depende da renda (ou do total disponível para gastar em alimentos).
gasto total:

Pode ser correlacionado com o erro do modelo de demanda, pois choca não observados (como preferências ou hábitos) influenciam simultaneamente o quanto a pessoa gasta e a maneira como ela distribui esse gasto entre os grupos alimentares;

Pode não ser exógeno, especialmente se a variável de interesse (como imposto ou preço) também afeta o gasto total.
so viola os pressupostos da regressão, gerando viés nas elasticidades estimadas.
Blundell & Robin (1999): A solução
Blundell e Robin propuseram uma regressão aumentada ("augmented regression") que utiliza uma variável instrumental para o gasto total
Blundell, R., & Robin, J.-M. (1999). Estimation in large and disaggregated demand systems: an estimator for conditionally linear systems. Journal of Applied Econometrics, 14(3), 209–232.
*/

egen gasto_total_atualhat = rowtotal(despesa1_hat despesa2_hat despesa3_hat despesa4_hat despesa5_hat despesa6_hat despesa7_hat despesa8_hat)
gen lndespesa = log(gasto_total_atualhat)


* 2. Calcular participação orçamentária (w_i)
*usando a correção por endogeneidade 
foreach var in 1 2 3 4 5 6 7 8  {
    gen w_despesahat`var' = despesa`var'_hat / gasto_total_atualhat
}
replace w_despesahat1 = despesa1_hat / gasto_total_atualhat
replace w_despesahat2 = despesa2_hat / gasto_total_atualhat
replace w_despesahat3 = despesa3_hat / gasto_total_atualhat
replace w_despesahat4 = despesa4_hat / gasto_total_atualhat
replace w_despesahat5 = despesa5_hat / gasto_total_atualhat
replace w_despesahat6 = despesa6_hat / gasto_total_atualhat
replace w_despesahat7 = despesa7_hat / gasto_total_atualhat
replace w_despesahat8 = despesa8_hat / gasto_total_atualhat
* 3. Checar se todas as w_ estão entre 0 e 1

foreach var in 1 2 3 4 5 6 7 8 {
    summ w_despesahat`var'
}


save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_com correcao endogeneidade.dta", replace
*o que vai ter na função AIDS: 

*participação do gasto com o bem 𝑖 i no total gasto (share)=

*𝑝𝑗  = preço do bem 𝑗j: preco_com_reforma1
*X = gasto total do consumidor: gasto_total_reforma1
*P = índice de preço (geralmente aproximado por um índice de Stone): lnm

* 4. o que eu preciso ter: 

*variavel de gasto total em nivel monetário : gasto_total_reforma1
*participações dos gastos (w_*):w_preco_com_reforma11 w_preco_com_reforma12 ... w_preco_com_reforma114

*****estou usando esse modelo nas analises	
demandsys quaids w_despesahat1 w_despesahat2 w_despesahat3 w_despesahat4 w_despesahat5 w_despesahat6 w_despesahat7 w_despesahat8, prices(preco_corr_preco_por_kg1 preco_corr_preco_por_kg2 preco_corr_preco_por_kg3 preco_corr_preco_por_kg4 preco_corr_preco_por_kg5 preco_corr_preco_por_kg6 preco_corr_preco_por_kg7 preco_corr_preco_por_kg8 ) expenditures(gasto_total_atualhat) 
	
		
*Mede a variação percentual na demanda (participação no gasto, Marshallianas) de um grupo de alimentos quando o preço desse ou de outro grupo muda em 1%.
*São chamadas "não compensadas" porque não controlam a variação na utilidade (bem-estar) — consideram efeito de substituição e efeito renda juntos.
estat elasticities, uncompensated	

estat elasticities if quintil_renda ==1, uncompensated
estat elasticities if quintil_renda ==5, uncompensated
	
estat elasticities , expenditure
estat elasticities if quintil_renda ==1, expenditure	
estat elasticities if quintil_renda ==5, expenditure	
	
*elasticidades-preço compensadas (Hicksianas) — essas medidas são mais úteis para analisar efeitos de substituição puros, ou seja, como o consumidor ajusta a demanda mantendo o mesmo nível de bem-estar (desconsidera o efeito renda).
*Mostram apenas o efeito substituição.
estat elasticities, compensated	
estat elasticities if quintil_renda ==1, compensated	
estat elasticities if quintil_renda ==5, compensated	

	save "C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_analise_AIDS.dta", replace
use	"C:\Users\GEIPS1\OneDrive\POS DOC\PROJETO IDEC\manuscrito 02- modelagem de taxas\Construção da analise\bancos\banco_analise_AIDS.dta", clear

*Calcular a Mudança na Participação de Gasto

*variavel de preço de mercado modelo 1: (as aliquotas aplicadas são as estimadas pela reforma):preco_corr_preco_com_reforma11 preco_corr_preco_com_reforma12 preco_corr_preco_com_reforma13 preco_corr_preco_com_reforma14 preco_corr_preco_com_reforma15 preco_corr_preco_com_reforma16 preco_corr_preco_com_reforma17 preco_corr_preco_com_reforma18


*variavel de preço de mercado modelo 2: (as aliquotas aplicadas são as estimadas pela reforma, mas adicona-se o IS para todos os alimentos ultraprocessados):preco_corr_preco_com_reforma21 preco_corr_preco_com_reforma22 preco_corr_preco_com_reforma23 preco_corr_preco_com_reforma24 preco_corr_preco_com_reforma25 preco_corr_preco_com_reforma26 preco_corr_preco_com_reforma27 preco_corr_preco_com_reforma28
# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")


# Table 1 - description of the subjects ----
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1)) |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Table 2 - effect of covariates on ALS occurence ----
# Estimated risk of ALS occurrence attributed to subject characteristics in the Danish Diet, Cancer and Health study cohort (logistic regression models, sample size: 498). 
table_2 <- results_POPs_ALS_occurrence$main$covar |> 
  as_flex_table()|>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")


# Figure 1 - boxplot expo ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
figure_1 <- results_descriptive$danish$POPs_group_boxplot_danish_by_als

# Figure 2 - forest plot quartiles ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (conditional logistic regression models, sample size: 498).
figure_2 <- results_POPs_ALS_occurrence$main$plot_quart
figure_2_bis_bis <- results_POPs_ALS_occurrence$main$plot_quart_bis_bis

# Figure 3 - gam models ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (generalized additive mixed models, sample size: 498).
figure_3_a <- wrap_plots(
  results_POPs_ALS_occurrence$main$plot_base_gam$PCB_DL,
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$PCB_DL,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$PCB_DL,
  
  results_POPs_ALS_occurrence$main$plot_base_gam$PCB_NDL, 
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$PCB_NDL,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$PCB_NDL, 
  
  results_POPs_ALS_occurrence$main$plot_base_gam$PCB_4, 
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$PCB_4, 
  wrap_elements(grid::nullGrob()),
  
  results_POPs_ALS_occurrence$main$plot_base_gam$OCP_HCB,
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$OCP_HCB,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$OCP_HCB, 
  ncol = 3)

figure_3_b <- wrap_plots(
  results_POPs_ALS_occurrence$main$plot_base_gam$ΣDDT,
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$ΣDDT,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$ΣDDT, 
  
  results_POPs_ALS_occurrence$main$plot_base_gam$OCP_β_HCH, 
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$OCP_β_HCH,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$OCP_β_HCH,
  
  results_POPs_ALS_occurrence$main$plot_base_gam$Σchlordane,
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$Σchlordane,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$Σchlordane,
  
  results_POPs_ALS_occurrence$main$plot_base_gam$ΣPBDE,
  results_POPs_ALS_occurrence$main$plot_adjusted_gam$ΣPBDE,
  results_POPs_ALS_occurrence$main$plot_copollutant_gam$ΣPBDE,
  ncol = 3, nrow = 4)

# Figure 4 - metanalysis ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort, the Finnish Mobile Clinic Health Examination Survey (FMC), its Follow-up Stage (FMCF), and the Mini-Finland Health Survey (MFH) (metanalysis, sample size: 788).
figure_4 <- results_POPs_ALS_occurrence$metanalysis$plot_metanalysis_quart


# Table S1 - exposure distribution ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_S1 <- 
  results_descriptive$danish$POPs_table_danish |>
  select(-Zero.count, -"% > LOQ", -"LOQ") |>
  rename(Exposures = variable) |>
  flextable()  |>
  add_footer_lines(
    'pg/ml.
    POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    PeCB, α-HCH and γ-HCH had more than 95% of the observations below the limit of quantification (LOQ) and were therefore excluded from our analyses. 
    HCB and β-HCH were not included in any group and were therefore studied alone.')|>
  theme_vanilla() |>  
  bold(part = "header") |>  
  bold(j = "Exposures", part = "body") |> 
  align(align = "center", part = "all") |>  
  align(j = "Exposures", align = "left", part = "all") |>
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# Table S2 - quartiles results ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (conditional logistic regressions, sample size: 498).
extra_rows <- results_POPs_ALS_occurrence$main$results_quart |>
  distinct(variable) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95%CI_base" = '-', "p.value_base" = '', "p.value_heterogeneity_base" = '', "p.value_trend_base" = '',
    "OR_adjusted" = '-', "95%CI_adjusted" = '-', "p.value_adjusted" = '', , "p.value_heterogeneity_adjusted" = '', "p.value_trend_adjusted" = '',
    "OR_copollutant" = '-', "95%CI_copollutant" ='-', "p.value_copollutant" = '', , "p.value_heterogeneity_copollutant" = '', "p.value_trend_copollutant" = '')

table_S2 <- 
  results_POPs_ALS_occurrence$main$main_results |>
  select(-lower_CI, -upper_CI, - p.value_raw) |>
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value, p.value_heterogeneity, p.value_trend)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart"), 
         contains("copollutant_quart")) 
colnames(table_S2) <- gsub('_quart', '', colnames(table_S2))

table_S2 <- table_S2 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  select(-p.value_heterogeneity_copollutant, -p.value_trend_copollutant) |>
  mutate(
    variable = gsub("OCP_", "", variable), 
    variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = fct_recode(variable, 
                          "Most prevalent PCBs" = "PCB_4",
                          "Dioxin-like PCBs" = "PCB_DL",
                          "Non-dioxin-like PCBs" = "PCB_NDL",
                          "β-HCH" = "β_HCH")) |>
  arrange(variable, quartiles) |>
  group_by(variable) |>
  mutate(p.value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p.value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p.value_trend_base = ifelse(quartiles == 'Quartile 1', p.value_trend_base[quartiles == 'Quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p.value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p.value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename('OR' = 'OR_base', '95% CI' = '95%CI_base', 'p-value' = 'p.value_base', 
         'Heterogeneity test' = 'p.value_heterogeneity_base', 'Trend test' = 'p.value_trend_base',
         'OR ' = 'OR_adjusted', '95% CI ' = '95%CI_adjusted', 'p-value ' = 'p.value_adjusted', 
         'Heterogeneity test ' = 'p.value_heterogeneity_adjusted', 'Trend test ' = 'p.value_trend_adjusted',
         ' OR ' = 'OR_copollutant', ' 95% CI ' = '95%CI_copollutant', ' p-value ' = 'p.value_copollutant')
table_S2 <- table_S2 |> flextable() |>
  add_footer_lines(
  "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as continuous variables.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.
  5Heterogeneity tests in outcome value across POP quartiles, matched on age and sex.
  6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, matched on age and sex.
  7Heterogeneity tests in outcome value across POP quartiles, matched on age and sex, and adjusted for smoking, BMI, serum total cholesterol, marital status, and education.
  8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, matched on age and sex, and adjusted for smoking, BMI, serum total cholesterol, marital status, and education.") |>
  add_header(
    "variable" = "Exposures", "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", "Heterogeneity test" = 'Base Model', "Trend test" = 'Base Model',
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", "Heterogeneity test " = 'Adjusted Model', "Trend test " = 'Adjusted Model',
    " OR " = "Co-pollutant model", " 95% CI " = "Co-pollutant model", " p-value " = "Co-pollutant model") |>
  merge_h(part = "header") |>
  merge_v(j = "variable") |>
  theme_vanilla() |>
  bold(j = "variable", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "variable", align = "left", part = "all") |> 
  merge_at(j = "variable", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


table_S2_bis_bis <- 
  results_POPs_ALS_occurrence$main$main_results |>
  select(-lower_CI, -upper_CI, - p.value_raw) |>
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart_bis_bis')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value, p.value_heterogeneity, p.value_trend)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart"), 
         contains("copollutant_quart_bis_bis")) 
colnames(table_S2_bis_bis) <- gsub('_quart', '', colnames(table_S2_bis_bis))
colnames(table_S2_bis_bis) <- gsub('_bis_bis', '', colnames(table_S2_bis_bis))

table_S2_bis_bis <- table_S2_bis_bis |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  # select(-p.value_heterogeneity_copollutant, -p.value_trend_copollutant) |>
  mutate(
    variable = gsub("OCP_", "", variable), 
    variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = fct_recode(variable, 
                          "Most prevalent PCBs" = "PCB_4",
                          "Dioxin-like PCBs" = "PCB_DL",
                          "Non-dioxin-like PCBs" = "PCB_NDL",
                          "β-HCH" = "β_HCH")) |>
  arrange(variable, quartiles) |>
  group_by(variable) |>
  mutate(p.value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p.value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p.value_trend_base = ifelse(quartiles == 'Quartile 1', p.value_trend_base[quartiles == 'Quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p.value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p.value_trend_adjusted[quartiles == 'Quartile 2'], ''), 
         p.value_heterogeneity_copollutant = ifelse(quartiles == 'Quartile 1', p.value_heterogeneity_copollutant[quartiles == 'Quartile 2'], ''), 
         p.value_trend_copollutant = ifelse(quartiles == 'Quartile 1', p.value_trend_copollutant[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename('OR' = 'OR_base', '95% CI' = '95%CI_base', 'p-value' = 'p.value_base', 
         'Heterogeneity test' = 'p.value_heterogeneity_base', 'Trend test' = 'p.value_trend_base',
         'OR ' = 'OR_adjusted', '95% CI ' = '95%CI_adjusted', 'p-value ' = 'p.value_adjusted', 
         'Heterogeneity test ' = 'p.value_heterogeneity_adjusted', 'Trend test ' = 'p.value_trend_adjusted',
         ' OR ' = 'OR_copollutant', ' 95% CI ' = '95%CI_copollutant', ' p-value ' = 'p.value_copollutant',
         ' Heterogeneity test ' = 'p.value_heterogeneity_copollutant', ' Trend test ' = 'p.value_trend_copollutant')
table_S2_bis_bis <- table_S2_bis_bis |> flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where a GAM is used in which the chemical of interest is categorized into quartiles, while the other pollutants are included as smooth terms.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.
  5Heterogeneity tests in outcome value across POP quartiles, matched on age and sex.
  6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, matched on age and sex.
  7Heterogeneity tests in outcome value across POP quartiles, matched on age and sex, and adjusted for smoking, BMI, serum total cholesterol, marital status, and education.
  8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, matched on age and sex, and adjusted for smoking, BMI, serum total cholesterol, marital status, and education.") |>
  add_header(
    "variable" = "Exposures", "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", "Heterogeneity test" = 'Base Model', "Trend test" = 'Base Model',
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", "Heterogeneity test " = 'Adjusted Model', "Trend test " = 'Adjusted Model',
    " OR " = "Co-pollutant model", " 95% CI " = "Co-pollutant model", " p-value " = "Co-pollutant model", " Heterogeneity test " = 'Co-pollutant model', " Trend test " = 'Co-pollutant model') |>
  merge_h(part = "header") |>
  merge_v(j = "variable") |>
  theme_vanilla() |>
  bold(j = "variable", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "variable", align = "left", part = "all") |> 
  merge_at(j = "variable", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(extra_rows)

# Table S3 - sensitivity analyse pollutant not summed (quartiles) ----
extra_rows <- results_POPs_ALS_occurrence$sensitivity_not_summed$sensitivity_results_not_summed_quart |>
  distinct(variable) |> 
  mutate(
    quartiles = " Quartile 1",
    "OR_base_not_summed" = '-', "95%CI_base_not_summed" = '-', "p.value_base_not_summed" = '',
    "OR_adjusted_not_summed" = '-', "95%CI_adjusted_not_summed" = '-', "p.value_adjusted_not_summed" = '')

table_S3 <- 
  results_POPs_ALS_occurrence$sensitivity_not_summed$sensitivity_results_not_summed_quart |>
  select(-lower_CI, -upper_CI, - p.value_raw) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart")) 
colnames(table_S3) <- gsub('_quart', '', colnames(table_S3))

table_S3 <- table_S3 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  mutate(
    quartiles = fct_recode(quartiles, 
      "Quartile 2" = "2",
      "Quartile 3" = "3",
      "Quartile 4" = "4"), 
    variable = gsub("OCP-", "", variable), 
    variable = gsub("_", "-", variable), 
    variable = fct_recode(variable, "Oxychlordane" = "oxychlordane", "Transnonachlor"= "transnonachlor"),
    variable = fct_relevel(variable, 
                           "PCB-118", "PCB-156", "PCB-28", "PCB-52", "PCB-74", "PCB-99",
                           "PCB-101", "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183",
                           "PCB-187", "HCB", "p,p'-DDE", "p,p'-DDT", "β-HCH", "Transnonachlor",
                           "Oxychlordane", "PBDE-47", "PBDE-99", "PBDE-153")) |>
  arrange(variable, quartiles) |>
  rename('OR' = 'OR_base_not_summed', '95% CI' = '95%CI_base_not_summed', 'p-value' = 'p.value_base_not_summed', 
         'OR ' = 'OR_adjusted_not_summed', '95% CI ' = '95%CI_adjusted_not_summed', 'p-value ' = 'p.value_adjusted_not_summed') |> 
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as continuous variables.
  2Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "variable" = "Exposures", "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
  merge_h(part = "header") |>
  merge_v(j = "variable") |>
  theme_vanilla() |>
  bold(j = "variable", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "variable", align = "left", part = "all") |> 
  merge_at(j = "variable", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(extra_rows)

# Table S4 - metanalysis ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort, the Finnish Mobile Clinic Health Examination Survey (FMC), its Follow-up Stage (FMCF), and the Mini-Finland Health Survey (MFH) (metanalysis, sample size: 788).
extra_rows <- results_POPs_ALS_occurrence$metanalysis$metanalysis_quart |>
  distinct(explanatory) |> 
  mutate(
    term = "Quartile 1",
      "OR_base" = '-', "95%CI_base" = '-', "p-value_base" = '',
      "OR_adjusted" = '-', "95%CI_adjusted" = '-', "p-value_adjusted" = '',
    "OR_copollutant" = '-', "95%CI_copollutant" = '-', "p-value_copollutant" = '')
    

table_S4 <- results_POPs_ALS_occurrence$metanalysis$metanalysis_quart |>
  select(model, explanatory, term, OR, `95%CI`, "p-value") |>
  pivot_wider(
    names_from = model,  
    values_from = c("OR", "95%CI", "p-value")) |>
  bind_rows(extra_rows) |>
  select(Exposures = explanatory, Quartiles = term,
         contains("base"), 
         contains("adjusted"), 
         contains("copollutant")) |>
  mutate(
    Exposures = gsub("OCP_", "", Exposures), 
    Exposures = fct_relevel(Exposures, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane"),
    Exposures = fct_recode(Exposures, 
                          "Most prevalent PCBs" = "PCB_4",
                          "Dioxin-like PCBs" = "PCB_DL",
                          "Non-dioxin-like PCBs" = "PCB_NDL",
                          "β-HCH" = "β_HCH")) |>
  arrange(Exposures, Quartiles) 

table_S4 <- table_S4 |> 
  rename('OR' = 'OR_base', '95% CI' = '95%CI_base', 'p-value' = 'p-value_base', 
         'OR ' = 'OR_adjusted', '95% CI ' = '95%CI_adjusted', 'p-value ' = 'p-value_adjusted', 
         ' OR ' = 'OR_copollutant', ' 95% CI ' = '95%CI_copollutant', ' p-value ' = 'p-value_copollutant') |>
  flextable() |>
  add_footer_lines(
  "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2Base and adjusted models are conditional logistic regressions matched for sex and age, and for municipality and serum thawing history (Finnish cohorts). Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC cohort for which education data was not available). Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as smooth terms in GAMs.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.") |>
  add_header(
    "Exposures" = "Exposures", "Quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", 
    " OR " = "Co-pollutant model", " 95% CI " = "Co-pollutant model", " p-value " = "Co-pollutant model") |>
  merge_h(part = "header") |>
  merge_v(j = "Exposures") |>
  theme_vanilla() |>
  bold(j = "Exposures", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "Exposures", align = "left", part = "all") |> 
  merge_at(j = "Exposures", part = "header") |>
  merge_at(j = "Quartiles", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(extra_rows)

# Figure S1 - DAG ----
# Directed acyclic graph of the relationship between pre-disease POP serum concentrations and ALS occurrence.

# Figure S2 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/figure_S2.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(results_descriptive$danish$POPs_heatmap_danish, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         number.digits = 1,
         tl.cex = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()


# Figure S3 - outliers sensitivity analysis ---- 
# Sensitivity analyses - Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations after removing extreme values in the Danish Diet, Cancer and Health study cohort (generalized additive mixed models, sample size: 498).
figure_S3_a <- wrap_plots(
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$PCB_DL,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$PCB_DL,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$PCB_DL,
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$PCB_NDL,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$PCB_NDL,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$PCB_NDL,
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$PCB_4,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$PCB_4,
  wrap_elements(grid::nullGrob()),
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$OCP_HCB,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$OCP_HCB,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$OCP_HCB,
  ncol = 3)

figure_S3_b <- wrap_plots(
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$ΣDDT,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$ΣDDT,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$ΣDDT,
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$OCP_β_HCH,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$OCP_β_HCH,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$OCP_β_HCH,
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$Σchlordane,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$Σchlordane,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$Σchlordane,
  
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_base_gam_outlier$ΣPBDE,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_adjusted_gam_outlier$ΣPBDE,
  results_POPs_ALS_occurrence$sensitivity_outliers$plot_copollutant_gam_outlier$ΣPBDE,
  ncol = 3)

# Figure S4 - not summed sensitivity analysis quartiles ----
figure_S4 <- results_POPs_ALS_occurrence$sensitivity_not_summed$plot_quart_sensi_not_summed 

# Figure S5 - not summed sensitivity analysis gam ----
figure_S5 <- wrap_plots(results_POPs_ALS_occurrence$sensitivity_not_summed$plot_adjusted_gam_not_summed)

# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_1.docx")
table_2 <- read_docx() |> body_add_flextable(table_2) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_2.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S1.docx")
table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S2.docx")
table_S2_bis_bis <- read_docx() |> body_add_flextable(table_S2_bis_bis)
print(table_S2_bis_bis, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S2_bis_bis.docx")
table_S3 <- read_docx() |> body_add_flextable(table_S3)
print(table_S3, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S3.docx")
table_S4 <- read_docx() |> body_add_flextable(table_S4)
print(table_S4, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S4.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_1.tiff",
  figure_1,
  height = 4,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_2.tiff",
  figure_2,
  height = 8,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_2_bis_bis.tiff",
  figure_2_bis_bis,
  height = 8,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_3_a.tiff",
  figure_3_a,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_3_b.tiff",
  figure_3_b,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_4.tiff",
  figure_4,
  height = 7,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_S3_a.tiff",
  figure_S3_a,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_S3_b.tiff",
  figure_S3_b,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_S4.tiff",
  figure_S4,
  height = 15,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_S5.tiff",
  figure_S5,
  height = 15,
  width = 9,
  units = "in")


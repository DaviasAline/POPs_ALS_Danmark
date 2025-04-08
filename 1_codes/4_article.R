# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2_analyses.R")


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
  # font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Table 2 - effect of covariates on ALS occurence ----
# Estimated risk of ALS occurrence attributed to subject characteristics in the Danish Diet, Cancer and Health study cohort (logistic regression models, sample size: 498). 
table_2 <- covar |> 
  as_flex_table()|>
  # font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")


# Figure 1 - boxplot expo ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
figure_1 <- descrip_expo_group_by_als

# Figure 2 - forest plot quartiles ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (conditional logistic regression models, sample size: 498).
figure_2 <- plot_quart

# Figure 3 - gamm models ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (generalized additive mixed models, sample size: 498).
figure_3_a <- wrap_plots(
  plot_base_gamm$PCB_DL,
  plot_adjusted_gamm$PCB_DL,
  plot_copollutant_gamm$PCB_DL,
  
  plot_base_gamm$PCB_NDL,
  plot_adjusted_gamm$PCB_NDL,
  plot_copollutant_gamm$PCB_NDL,
  
  plot_base_gamm$PCB_4,
  plot_adjusted_gamm$PCB_4,
  wrap_elements(grid::nullGrob()),
  
  plot_base_gamm$HCB,
  plot_adjusted_gamm$HCB,
  plot_copollutant_gamm$HCB, 
  ncol = 3)

figure_3_b <- wrap_plots(
  plot_base_gamm$ΣDDT,
  plot_adjusted_gamm$ΣDDT,
  plot_copollutant_gamm$ΣDDT,
  
  plot_base_gamm$β_HCH,
  plot_adjusted_gamm$β_HCH,
  plot_copollutant_gamm$β_HCH,
  
  plot_base_gamm$Σchlordane,
  plot_adjusted_gamm$Σchlordane,
  plot_copollutant_gamm$Σchlordane,
  
  plot_base_gamm$ΣPBDE,
  plot_adjusted_gamm$ΣPBDE,
  plot_copollutant_gamm$ΣPBDE,
  ncol = 3)


# Table S1 - exposure distribution ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_S1 <- descrip_num(data = bdd_danish, 
                        vars = c(POPs, "PCB_DL", "PCB_NDL", "ΣDDT", 
                                 "Σchlordane", "ΣPBDE" )) |>
  select(-Zero.count) |>
  rename(Exposures = variable) |>
  mutate(
    Exposures = fct_recode(Exposures, 
                           "Dioxin-like PCBs" = 'PCB_DL',
                           "Non-dioxin-like PCBs" = 'PCB_NDL', 
                           'p,p’-DDT'  = 'pp_DDT', 
                           'p,p’-DDE' = "pp_DDE"), 
    Exposures = gsub("_", "-", Exposures),
    Exposures = fct_relevel(Exposures, 
                            "Dioxin-like PCBs", "PCB-118", "PCB-156", "Non-dioxin-like PCBs", "PCB-28", "PCB-52",
                            "PCB-74", "PCB-99", "PCB-101", "PCB-138", "PCB-153", "PCB-170",
                            "PCB-180", "PCB-183", "PCB-187", "HCB", "PeCB", "ΣDDT", 'p,p’-DDE',
                            'p,p’-DDT', "α-HCH", "β-HCH", "γ-HCH", "Σchlordane", "Transnonachlor",
                            "Oxychlordane", "ΣPBDE", "BDE-47", "BDE-99", "BDE-153")) |>
  arrange(Exposures) |>
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
extra_rows <- results_quart %>%
  distinct(variable) %>% 
  mutate(
    quartiles = "1",
    "OR_base" = '-', "95%CI_base" = '-', "p.value_base" = '', "p.value_heterogeneity_base" = '', "p.value_trend_base" = '',
    "OR_adjusted" = '-', "95%CI_adjusted" = '-', "p.value_adjusted" = '', , "p.value_heterogeneity_adjusted" = '', "p.value_trend_adjusted" = '',
    "OR_copollutant" = '-', "95%CI_copollutant" ='-', "p.value_copollutant" = '', , "p.value_heterogeneity_copollutant" = '', "p.value_trend_copollutant" = '')

table_S2 <- 
  main_results |>
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
    variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = fct_recode(variable, 
                          "Most prevalent PCBs" = "PCB_4",
                          "Dioxin-like PCBs" = "PCB_DL",
                          "Non-dioxin-like PCBs" = "PCB_NDL",
                          "β-HCH" = "β_HCH")) |>
  arrange(variable, quartiles) |>
  group_by(variable) %>%
  mutate(p.value_heterogeneity_base = ifelse(quartiles == '1', p.value_heterogeneity_base[quartiles == '2'], ''), 
         p.value_trend_base = ifelse(quartiles == '1', p.value_trend_base[quartiles == '2'], ''),
         p.value_heterogeneity_adjusted = ifelse(quartiles == '1', p.value_heterogeneity_adjusted[quartiles == '2'], ''), 
         p.value_trend_adjusted = ifelse(quartiles == '1', p.value_trend_adjusted[quartiles == '2'], '')) %>%
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
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(extra_rows)

# Table S3 - sensitivity analyse pollutant not summed (quartiles) ----
extra_rows <- sensitivity_results_not_summed_quart %>%
  distinct(variable) %>% 
  mutate(
    quartiles = "1",
    "OR_base_not_summed" = '-', "95%CI_base_not_summed" = '-', "p.value_base_not_summed" = '',
    "OR_adjusted_not_summed" = '-', "95%CI_adjusted_not_summed" = '-', "p.value_adjusted_not_summed" = '')

Table_S3 <- 
  sensitivity_results_not_summed_quart |>
  select(-lower_CI, -upper_CI, - p.value_raw) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart")) 
colnames(Table_S3) <- gsub('_quart', '', colnames(Table_S3))

Table_S3 <- Table_S3 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  mutate(
    variable = gsub("_", "-", variable), 
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
    "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as continuous variables.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.") |>
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

# Figure S1 - DAG ----
# Directed acyclic graph of the relationship between pre-disease POP serum concentrations and ALS occurrence.

# Figure S2 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
figure_S2 <- bdd_danish %>% 
  select(all_of(POPs), all_of(POPs_group)) |>
  rename(
    "Dioxin-like PCBs" = PCB_DL,
    "Non-dioxin-like PCBs" = PCB_NDL, 
    "p,p’-DDT"  = pp_DDT, 
    "p,p’-DDE" = pp_DDE) |>
  rename_with(~ gsub("_", "-", .x)) |>  
  rename_with(~ gsub("BDE", "PBDE", .x)) |>  
  select(
    "Dioxin-like PCBs", "PCB-118", "PCB-156", "Non-dioxin-like PCBs", "PCB-28", "PCB-52",
    "PCB-74", "PCB-99", "PCB-101", "PCB-138", "PCB-153", "PCB-170",
    "PCB-180", "PCB-183", "PCB-187", "HCB", "ΣDDT", "p,p’-DDE",
    "p,p’-DDT",  "β-HCH",  "Σchlordane", "Transnonachlor",
    "Oxychlordane", "ΣPBDE" = "ΣPPBDE", "PBDE-47", "PBDE-99", "PBDE-153")

figure_S2 <- cor(figure_S2, 
                 use = "pairwise.complete.obs", 
                 method = "pearson")

plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/figure_S2.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(figure_S2, 
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
  plot_base_gamm_outlier$PCB_DL,
  plot_adjusted_gamm_outlier$PCB_D,
  plot_copollutant_gamm_outlier$PCB_DL,
  
  plot_base_gamm_outlier$PCB_NDL,
  plot_adjusted_gamm_outlier$PCB_NDL,
  plot_copollutant_gamm_outlier$PCB_NDL,
  
  plot_base_gamm_outlier$PCB_4,
  plot_adjusted_gamm_outlier$PCB_4,
  wrap_elements(grid::nullGrob()),
  
  plot_base_gamm_outlier$HCB,
  plot_adjusted_gamm_outlier$HCB,
  plot_copollutant_gamm_outlier$HCB,
  ncol = 3)

figure_S3_b <- wrap_plots(
  plot_base_gamm_outlier$ΣDDT,
  plot_adjusted_gamm_outlier$ΣDDT,
  plot_copollutant_gamm_outlier$ΣDDT,
  
  plot_base_gamm_outlier$β_HCH,
  plot_adjusted_gamm_outlier$β_HCH,
  plot_copollutant_gamm_outlier$β_HCH,
  
  plot_base_gamm_outlier$Σchlordane,
  plot_adjusted_gamm_outlier$Σchlordane,
  plot_copollutant_gamm_outlier$Σchlordane,
  
  plot_base_gamm_outlier$ΣPBDE,
  plot_adjusted_gamm_outlier$ΣPBDE,
  plot_copollutant_gamm_outlier$ΣPBDE,
  ncol = 3)

# Figure S4 - not summed sensitivity analysis quartiles ----
figure_S4 <- plot_quart_sensi_not_summed 

# Figure S5 - not summed sensitivity analysis gamm ----
figure_S5 <- wrap_plots(plot_adjusted_gamm_not_summed)

# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_1.docx")
table_2 <- read_docx() |> body_add_flextable(table_2) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_2.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S1.docx")
table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_S2.docx")

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
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_S2.tiff",
  figure_S2,
  height = 10,
  width = 10,
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


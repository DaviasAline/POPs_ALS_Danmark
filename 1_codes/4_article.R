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
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Table 2 - effect of covariates on ALS occurence ----
# Estimated risk of ALS occurrence attributed to subject characteristics in the Danish Diet, Cancer and Health study cohort (logistic regression models, sample size: 498). 
table_2 <- covar |> 
  as_flex_table()|>
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")


# Figure 1 - boxplot expo ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
figure_1 <- descrip_expo_group_by_als

# Figure 2 - gamm models ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (generalized additive mixed models, sample size: 498).
figure_2_a <- wrap_plots(
  plot_base_gamm$PCB_DL + ggtitle('Base model'),
  plot_adjusted_gamm$PCB_DL + ggtitle('Adjusted model'),
  plot_copollutant_gamm$PCB_DL + ggtitle('Copollutant model'),
  
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

figure_2_b <- wrap_plots(
  plot_base_gamm$ΣDDT + ggtitle('Base model'),
  plot_adjusted_gamm$ΣDDT + ggtitle('Adjusted model'),
  plot_copollutant_gamm$ΣDDT + ggtitle('Copollutant model'),
  
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

# Figure 3 - forest plot quartiles ----
# Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (conditional logistic regression models, sample size: 498).
figure_3 <- plot_quart


# Table S1 - exposure distribution ----
# Distribution of pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_s1 <- descrip_num(data = bdd_danish, 
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
    "OR_base" = '-', "95%CI_base" = '-', "p.value_base" = '', 
    "OR_adjusted" = '-', "95%CI_adjusted" = '-', "p.value_adjusted" = '',
    "OR_copollutant" = '-', "95%CI_copollutant" ='-', "p.value_copollutant" = '')

table_s2 <- results_quart |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  mutate(
    variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = fct_recode(variable, 
                          "Most prevalent PCBs" = "PCB_4",
                          "Dioxin-like PCBs" = "PCB_DL",
                          "Non-dioxin-like PCBs" = "PCB_NDL",
                          "β-HCH" = "β_HCH")) |>
  arrange(variable, quartiles) |>
  rename('OR' = 'OR_base', '95% CI' = '95%CI_base', 'p-value' = 'p.value_base',
         'OR ' = 'OR_adjusted', '95% CI ' = '95%CI_adjusted', 'p-value ' = 'p.value_adjusted',
         ' OR ' = 'OR_copollutant', ' 95% CI ' = '95%CI_copollutant', ' p-value ' = 'p.value_copollutant') |>
  flextable() |>
  add_footer_lines(
  "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as continuous variables.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.") |>
  add_header(
    "variable" = "Exposures", "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model",
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

# Figure S1 - DAG ----
# Directed acyclic graph of the relationship between pre-disease POP serum concentrations and ALS occurrence.

# Figure S2 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP serum concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
figure_S2 <- heatmap_POPs

# Figure S3 - outliers sensitivity analysis ---- 
# Sensitivity analyses - Estimated risk of ALS occurrence attributed to pre-disease POP serum concentrations after removing extreme values in the Danish Diet, Cancer and Health study cohort (generalized additive mixed models, sample size: 498).
figure_S3_a <- wrap_plots(
  plot_base_gamm_outlier$PCB_DL + ggtitle('Base model without outliers'),
  plot_adjusted_gamm_outlier$PCB_DL + ggtitle('Adjusted model without outliers'),
  plot_copollutant_gamm_outlier$PCB_DL + ggtitle('Copollutant model without outliers'),
  
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
  plot_base_gamm_outlier$ΣDDT + ggtitle('Base model without outliers'),
  plot_adjusted_gamm_outlier$ΣDDT + ggtitle('Adjusted model without outliers'),
  plot_copollutant_gamm_outlier$ΣDDT + ggtitle('Copollutant model without outliers'),
  
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

# Figure S4 - lipid profile sensitivity analysis ----

# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_1.docx")
table_2 <- read_docx() |> body_add_flextable(table_2) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_2.docx")

table_s1 <- read_docx() |> body_add_flextable(table_s1)
print(table_s1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_s1.docx")
table_s2 <- read_docx() |> body_add_flextable(table_s2)
print(table_s2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_s2.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_1.tiff",
  figure_1,
  height = 4,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_2_a.tiff",
  figure_2_a,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_2_b.tiff",
  figure_2_b,
  height = 12,
  width = 9,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_3.tiff",
  figure_3,
  height = 8,
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


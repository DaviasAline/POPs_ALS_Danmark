# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2_analyses.R")


# Figure 1 (boxplot expo) ----
# Distribution of the POP exposures in the Danish Diet, Cancer and Health study cohort (sample size = 498).
figure_1 <- descrip_expo_group

# Figure 2 (gamm models) ----
figure_2 <- wrap_plots(
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
  ncol = 3
)

# Figure 3 (forest plot quartiles) ----
# Associations between POPs exposure (quartiles) and ALS occurence in conditional logistic regressions 
figure_3 <- plot_quart

# Table 1 ----
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size = 498).
table_1 <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no') |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() 

table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_1.docx")

# Table 2 - effect of covariates on ALS occurence ----
table_2 <- covar |> as_flex_table()
table_2 <- read_docx() |> body_add_flextable(table_2) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/table_2.docx")

# Supplementary figure 1 (dag) ----
# dag 

# Supplementary figure 2 - outliers sensitivity analysis ---- 
figure_S2 <- wrap_plots(
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
  ncol = 3
)

# Supplementary figure 3 - lipid profile sensitivity analysis ----

# Supplementary table 1 - exposure distribution ----
table_s1 <- descrip_num(data = bdd_danish, vars = POPs_group) |>
  select(-Zero.count) |>
  mutate(
    variable = fct_recode(variable, 
                          "Dioxin_like PCBs" = 'PCB_DL',
                          "Non-dioxin-like PCBs" = 'PCB_NDL', 
                          "Most prevalent PCBs" = 'PCB_4', 
                          "β-HCH" = 'β_HCH')) |>
  rename(Exposure = variable)

# Supplementary table 2 - quartiles results ----
table_s2 <- results_quart |>
  mutate(
    variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = fct_recode(variable, 
        "Most prevalent PCBs" = "PCB_4",
        "Dioxin-like PCBs" = "PCB_DL",
        "Non-dioxin-like PCBs" = "PCB_NDL",
        "β-HCH" = "β_HCH")) |>
  arrange(variable)

# Export ----
write_xlsx(
  list("Table 1" = table_1, 
       "Table 2" = table_2, 
       "Table S1" = table_s1, 
       "Table S2" = table_s2), 
  path = "~/Documents/POP_ALS_2025_02_03/2_output/article.xlsx"
)




ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_gamm_main (figure 2).tiff",
  figure_2,
  height = 24,
  width = 12,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/figure_gamm_outlier (figure s2).tiff",
  figure_S2,
  height = 24,
  width = 12,
  units = "in")



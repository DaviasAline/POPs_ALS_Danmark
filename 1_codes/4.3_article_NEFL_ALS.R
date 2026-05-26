# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.6_analyses_proteomic_ALS_occurrence.R")


# Table 1 - Subject characteristics description ---- 
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- 
  bdd_danish |>
  filter(match != 159) |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"),
    follow_up = follow_up/12, 
    follow_up_death = follow_up_death/12) |>
  select(
    als, birth_year, baseline_age, diagnosis_age, follow_up,follow_up_death, 
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, fS_Kol, 
    proteomic_neuro_explo_NEFL) |>
  tbl_summary(by = als, 
              missing = 'no', 
              statistic = list(all_continuous() ~ "{mean} ({sd})", 
                               follow_up ~ "{median} ({p25}, {p75})", 
                               follow_up_death ~ "{median} ({p25}, {p75})", 
                               alcohol ~ "{median} ({p25}, {p75})", 
                               all_categorical() ~ "{n} ({p})"),
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 1, 
                            follow_up_death ~ 1, 
                            bmi ~ 1, 
                            fS_Kol ~ 1, 
                            proteomic_neuro_explo_NEFL ~ 1), 
              label = list(proteomic_neuro_explo_NEFL ~ "Neurofilament light polypeptide (NEFL)")) |>
  bold_labels() |>
  # add_p(include = -c('diagnosis_age', 'follow_up')) |>
  add_overall() |>
  add_n(statistic = "{N_miss} ({p_miss})", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# Figure 1 - Descriptive figure of NEFL distribution (density and boxplots) ----
figure_1 <- 
  results_descriptive$danish$figure_NEFL_2

# Figure 2 - Base and adjusted logistic regressions (ALS risk) ----
figure_2 <- 
  results_proteomic_ALS_occurrence$Nfl_results$NfL_sd_ALS_figure_sensi_follow_up_adj

# Figure 3 - Additional analysis (NfL level over follow-up time) ----
figure_3 <- 
  results_proteomic_ALS_occurrence$additional_analysis_2$figure_NEFL_over_time_sensi_1

# Figure 4 - Additional analysis (AUC) ----
figure_4 <- 
  results_proteomic_ALS_occurrence$additional_analysis_4$additional_analysis_4_figure_unadjusted_pattern +
  labs(title = "") + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE))


# Table S1 ----
table_S1 <-
  bind_rows(
    bdd_danish |> filter(als == 1) |> filter(!match == 159) |> descrip_num(vars = "proteomic_neuro_explo_NEFL") |> mutate(variable = fct_recode(variable, "Cases" = "proteomic_neuro_explo_NEFL")), 
    bdd_danish |> filter(als == 0) |> filter(!match == 159) |> descrip_num(vars = "proteomic_neuro_explo_NEFL") |> mutate(variable = fct_recode(variable, "Controls" = "proteomic_neuro_explo_NEFL"))) |>
  select("ALS status" = variable, N, everything(), -Zero.count) |>
  flextable()  |>
  theme_vanilla() |>  
  bold(part = "header") |>  
  bold(j = "ALS status", part = "body") |> 
  align(align = "center", part = "all") |>  
  align(j = "ALS status", align = "left", part = "all") |>
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# Table S2 ----
table_S2 <- 
  results_proteomic_ALS_occurrence$NfL_results$NfL_sd_ALS_table_sensi_follow_up_sex_base_adj

# Table S3 ----
table_S3 <- 
  results_descriptive$danish$table_other_diag_1

# Figure S1 ----
figure_S1 <- 
  results_descriptive$danish$figure_other_diag


# Export ----
## tables ----
table_1 <- read_docx() |> body_add_flextable(table_1)                           # covariate descriptive table 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/table_1.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)                         # NfL descriptive table 
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)                         # Conditional logistic regressions (als risk)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/table_S2.docx")

table_S3 <- read_docx() |> body_add_flextable(table_S3)                         # Conditional logistic regressions (als risk)
print(table_S3, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/table_S3.docx")


## figures ----
ggsave(                                                                         # NfL descriptive figure 
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/figure_1.tiff",
  figure_1,
  height = 5,
  width = 10,
  units = "in")

ggsave(                                                                         # Forest plots of logistic regressions results (ALS risk)
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/figure_2.tiff",
  figure_2,
  height = 5,
  width = 8,
  units = "in")

ggsave(                                                                         # LOESS curve of NfL ratios depending on time to diagnosis
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/figure_3.tiff",
  figure_3,
  height = 5.5,
  width = 10,
  units = "in")

ggsave(                                                                         # AUC curve
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/figure_4.tiff",
  figure_4,
  height = 9,
  width = 12,
  units = "in")

ggsave(                                                                         # AUC curve
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NfL_ALS/figure_S1.tiff",
  figure_S1,
  height = 6,
  width = 9,
  units = "in")


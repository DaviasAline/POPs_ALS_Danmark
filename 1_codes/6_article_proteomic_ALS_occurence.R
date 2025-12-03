# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.7_analyses_proteomic_ALS_occurrence.R")


# Table 1 - Subject characteristics description ---- 
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- bdd |>
  filter(study == "Danish") |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"),
    follow_up = follow_up/12) |>
  select(
    als, birth_year, baseline_age, diagnosis_age, follow_up,  
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, fS_Kol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 0, 
                            bmi ~ 1, 
                            fS_Kol ~ 1)) |>
  bold_labels() |>
  add_p(include = -c('diagnosis_age', 'follow_up')) |>
  add_overall() |>
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# Table 2 - Main results ----
table_2 <- results_proteomic_ALS_occurrence$sensi_1$proteomic_sd_ALS_table_sensi_1


# Figure 1 - Main results (volcano plot) ----
figure_1 <- wrap_plots(
  list(results_proteomic_ALS_occurrence$sensi_1$proteomic_sd_ALS_base_figure_sensi_1 + theme(legend.position = "none") + labs(title = "Base models"), 
       results_proteomic_ALS_occurrence$sensi_1$proteomic_sd_ALS_adjusted_figure_sensi_1 + theme(legend.position = "bottom") + labs(title = "Adjusted models")), 
  ncol = 2)


# Figure 2 - Sensitivity analysis (volcano plot filtered to > 5 years of follow-up) ----
figure_2 <- wrap_plots(
  list(results_proteomic_ALS_occurrence$sensi_1_2$proteomic_sd_ALS_base_figure_sensi_1 + theme(legend.position = "none") + labs(title = "Base models"), 
       results_proteomic_ALS_occurrence$sensi_1_2$proteomic_sd_ALS_adjusted_figure_sensi_1  + theme(legend.position = "bottom") + labs(title = "Adjusted models")), 
  ncol = 2)


# Figure 3 - Additional analysis (NEFL level over follow-up time) ----
figure_3 <- results_proteomic_ALS_occurrence$additional_analysis_2$figure_NEFL_over_time_sensi_1


# Supplementary table 1 - Sensitivity analysis (volcano plot filtered to > 5 years of follow-up) ----
table_S1 <- results_proteomic_ALS_occurrence$sensi_1_2$proteomic_sd_ALS_table_sensi_1_2

# Supplementary figure 1 - heatmap proteomic neuro explo ----
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S1.tiff", 
     units = "mm", 
     width = 200, 
     height = 200, 
     res = 300)
corrplot(results_descriptive$danish$proteomic_heatmap_danish_neuro_explo, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black',
         order = "hclust",         # regroupement des variables selon niveau de cor
         tl.srt = 45, 
         # addCoef.col = "black",
         # number.cex = 0.8,
         # number.digits = 1,
         tl.cex = 0.4,             # taille police des variables
         col = rev(COL2(diverging = "RdYlBu")))


# Supplementary figure 2 - heatmap proteomic immune response  ----
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S2.tiff", 
     units = "mm", 
     width = 200, 
     height = 200, 
     res = 300)
corrplot(results_descriptive$danish$proteomic_heatmap_danish_immun_res, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black',
         order = "hclust",         # regroupement des variables selon niveau de cor
         tl.srt = 45, 
         # addCoef.col = "black",
         # number.cex = 0.8,
         # number.digits = 1,
         tl.cex = 0.4,             # taille police des variables
         col = rev(COL2(diverging = "RdYlBu")))

# Supplementary figure 3 - heatmap proteomic metabolism ----
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S3.tiff", 
     units = "mm", 
     width = 200, 
     height = 200, 
     res = 300)
corrplot(results_descriptive$danish$proteomic_heatmap_danish_metabolism, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black',
         order = "hclust",         # regroupement des variables selon niveau de cor
         tl.srt = 45, 
         # addCoef.col = "black",
         # number.cex = 0.8,
         # number.digits = 1,
         tl.cex = 0.4,             # taille police des variables
         col = rev(COL2(diverging = "RdYlBu")))


# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/table_1.docx")
table_2 <- read_docx() |> body_add_flextable(table_2) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/table_2.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/table_S1.docx")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_1.tiff",
  figure_1,
  height = 8,
  width = 15,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_2.tiff",
  figure_2,
  height = 8,
  width = 15,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_3.tiff",
  figure_3,
  height = 4,
  width = 8,
  units = "in")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S1.tiff",
  figure_S1,
  height = 10,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S2.tiff",
  figure_S2,
  height = 10,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_proteomics_ALS_occurence/figure_S3.tiff",
  figure_S3,
  height = 10,
  width = 10,
  units = "in")


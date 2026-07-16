# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

results_proteomic_ALS_occurrence <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.1_results_proteomic_ALS_occurrence.rds")
results_proteomic_ALS_occurrence_SL <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.2_results_proteomic_ALS_occurrence_SL.rds")
results_proteomic_ALS_occurrence_tidymodels <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.3_results_proteomic_ALS_occurrence_tidymodels.rds")
results_proteomic_ALS_occurrence_y_to_als <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.4_results_proteomic_ALS_occurrence_y_to_als.rds") 


# Table 1 - Subject characteristics description ---- 
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- bdd_danish |>
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


# Figure 1 - Main results (volcano plot) ----
figure_1 <- wrap_plots(
  list(results_proteomic_ALS_occurrence$sensi_1$proteomic_sd_ALS_adjusted_figure_sensi_1 +  # All + removing NfL outlier 
         theme(legend.position = "none") + 
         labs(title = "All cases and controls (N=495)", x = "") +
         scale_x_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5)) +
         scale_y_continuous(limits = c(0, 4.6), breaks = seq(0, 4.6, by = 1)), 
       results_proteomic_ALS_occurrence$sensi_2$proteomic_sd_ALS_adjusted_figure_sensi_2 +  # < 5 years of follow-up
         theme(legend.position = "none") + 
         labs(title = "Years to ALS < 5 years (N=51)", x = "", y = "") +
         scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
         scale_y_continuous(limits = c(0, 4.6), breaks = seq(0, 4.6, by = 1)), 
       results_proteomic_ALS_occurrence$sensi_1_3_4$proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 + # > 5 years follow- up + filtering follow-up <= 50% + removing NfL outlier
         theme(legend.position = "none") + 
         labs(title = "Years to ALS between 5 and 14.6 years (N=225)") +
         scale_x_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5)) +
         scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)),  
       results_proteomic_ALS_occurrence$sensi_1_3_5$proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 +  # > 5 years follow- up + filtering follow-up > 50% + removing NfL outlier
         theme(legend.position = "bottom") + 
         labs(title = "Years to ALS > 14.6 years (N=219)", y = "") +
         scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
         scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1))), 
  ncol = 2, widths = c(2.7, 4), heights = c(4.6, 3.2))


# Supplementary table S1 - Main results ----
table_S1 <- 
  results_proteomic_ALS_occurrence$main$main_results |> 
  filter(term == "Continuous") |>
  filter(model == "adjusted") |>
  filter(analysis %in% c("sensi_1", "sensi_2", "sensi_1_3_4", "sensi_1_3_5", "sensi_1_7_male", "sensi_1_7_female")) |>
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(analysis, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value")) |>
  select(
    protein_group, explanatory,
    OR_sensi_1, `95% CI_sensi_1`, p_value_sensi_1,
    OR_sensi_2, `95% CI_sensi_2`, p_value_sensi_2,
    OR_sensi_1_3_4, `95% CI_sensi_1_3_4`, p_value_sensi_1_3_4,
    OR_sensi_1_3_5, `95% CI_sensi_1_3_5`, p_value_sensi_1_3_5,
    OR_sensi_1_7_male, `95% CI_sensi_1_7_male`, p_value_sensi_1_7_male,
    OR_sensi_1_7_female, `95% CI_sensi_1_7_female`, p_value_sensi_1_7_female) |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", 
         "OR " = "OR_sensi_2", "95% CI " = "95% CI_sensi_2", "p-value " = "p_value_sensi_2", 
         " OR" = "OR_sensi_1_3_4", " 95% CI" = "95% CI_sensi_1_3_4", " p-value" = "p_value_sensi_1_3_4",
         " OR " = "OR_sensi_1_3_5", " 95% CI " = "95% CI_sensi_1_3_5", " p-value " = "p_value_sensi_1_3_5", 
         "OR  " = "OR_sensi_1_7_male", "95% CI  " = "95% CI_sensi_1_7_male", "p-value  " = "p_value_sensi_1_7_male",
         " OR  " = "OR_sensi_1_7_female", " 95% CI  " = "95% CI_sensi_1_7_female", " p-value  " = "p_value_sensi_1_7_female") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex (except the sex-stratified models), and adjusted for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Main analyses\nAll cases and controls (N=495)", "95% CI" = "Main analyses\nAll cases and controls (N=495)", "p-value" = "Main analyses\nAll cases and controls (N=495)", 
    "OR " = "Sensitivity analyses\nFiltered to years to ALS < 5 years (N=51)", "95% CI " = "Sensitivity analyses\nFiltered to years to ALS < 5 years (N=51)", "p-value " = "Sensitivity analyses\nFiltered to years to ALS < 5 years (N=51)", 
    " OR" = "Sensitivity analyses\nFiltered to years to ALS between 5 and 14.6 years (N=225)", " 95% CI" = "Sensitivity analyses\nFiltered to years to ALS between 5 and 14.6 years (N=225)", " p-value" = "Sensitivity analyses\nFiltered to years to ALS between 5 and 14.6 years (N=225)", 
    " OR " = "Sensitivity analyses\nFiltered to years to ALS > 14.6 years (N=219)", " 95% CI " = "Sensitivity analyses\nFiltered to years to ALS > 14.6 years (N=219)", " p-value " = "Sensitivity analyses\nFiltered to years to ALS > 14.6 years (N=219)", 
    "OR  " = "Sensitivity analyses\nFiltered to males (N=303)", "95% CI  " = "Sensitivity analyses\nFiltered to males (N=303)", "p-value  " = "Sensitivity analyses\nFiltered to males (N=303)", 
    " OR  " = "Sensitivity analyses\nFiltered to females (N=192)", " 95% CI  " = "Sensitivity analyses\nFiltered to females (N=192)", " p-value  " = "Sensitivity analyses\nFiltered to females (N=192)") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")



# Supplementary figure 1 - heatmap proteomic neuro explo ----
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S1.tiff", 
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
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S2.tiff", 
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
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S3.tiff", 
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


# Supplementary figure 4 - Sensitivity analysis (vocanoplot filtered on sex) ----
figure_S4 <-
  results_proteomic_ALS_occurrence$main$main_results |>
  filter(model == "adjusted" &                                                  # select only adjusted results
           term != "Continuous" &                                               # select only continuous results
           analysis %in% c("main", "sensi_1_7_female", "sensi_1_7_male")) |>     
  mutate(
    analysis = fct_recode(analysis, 
                          "All case and controls (N=495)" = "main", 
                          "Females (N=192)" = "sensi_1_7_female", 
                          "Males (N=303)" = "sensi_1_7_male"), 
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(
    data = ~filter(.x, p_value_raw < 0.05),       
    aes(label = explanatory), 
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "grey20", 
    color = "black") +
  scale_color_manual(
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  facet_grid(cols = vars(analysis))

# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_1.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_1.tiff",
  figure_1,
  height = 10,
  width = 13, 
  units = "in")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_S2.docx")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S1.tiff",
  figure_S1,
  height = 10,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S2.tiff",
  figure_S2,
  height = 10,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S3.tiff",
  figure_S3,
  height = 10,
  width = 10,
  units = "in")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S4.tiff",
  figure_S4,
  height = 8,
  width = 15,
  units = "in")



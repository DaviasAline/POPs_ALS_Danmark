# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")


# Table 1 - description of the subjects ----
# Description of the subject characteristics of the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263)
table_1 <- results_descriptive$danish$covar_danish_cases |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Figure 1 - description of the POP levels (boxplots) ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
figure_1 <- results_descriptive$danish$POPs_group_boxplot_danish_by_death

# Figure 2 - forest plot expo - ALS survival (danish EPIC cohort) ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the Danish Diet, Cancer and Health cohort (cox models by exposure quartiles; n = 166).
figure_2 <- 
  results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(term == "Continuous") |>
  filter(analysis %in% c("main", "sensi_1")) |>
  filter(model %in% c("adjusted", "ERS")) |>
  mutate(explanatory = factor(explanatory, levels = c(POPs_group_labels, "Environmental risk score" = "ERS_score_from_elastic_net_sensi_1")),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!c(POPs_group_labels, "Environmental risk score" = "ERS_score_from_elastic_net_sensi_1"))) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR_raw, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  scale_color_manual(values = c("p-value≤0.05" = "red", "p-value>0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()


# Table S1 - description of the POP levels (table) ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
table_S1 <- results_descriptive$danish$POPs_table_danish_by_death |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  |>
  add_footer_lines("POPs were summed as follows (pg/ml): most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE; Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.")

# Table S2 - POPs - ALS survival among the Danish cohort ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the Danish Diet, Cancer and Health cohort (cox models by exposure quartiles; n = 166).
table_S2 <- 
  results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(term == "Continuous") |>
  filter(!model == "copollutant") |> 
  filter(analysis %in% c("main", "sensi_1")) |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |> 
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |> 
  mutate(explanatory = fct_recode(explanatory, !!!c(POPs_group_labels, "Environmental risk score" = "ERS_score_from_elastic_net_sensi_1"))) |> 
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.
    2All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking, BMI and marital status. 
    3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")



# Table S3 - POPs - ALS survival among the Danish cohort (mixture model) ----
# Association between pre-diagnostic POP mixture and survival among ALS cases from the Danish Diet, Cancer and Health cohort (elastic net model; n = 166).
table_S3 <- results_POPs_ALS_survival$sensi2$POPs_sd_ALS_table_sensi2_ERS_danish

# Table S4 - Sensitivity analysis - POPs - ALS survival among the Danish cohort ----
# Association between pre-diagnostic POP groups and survival among ALS cases from the Danish Diet, Cancer and Health cohort (Cox regression modelsl; n = 166).
table_S4 <- results_POPs_ALS_survival$sensi4$POPs_sd_ALS_table_danish_sensi_4


# Figure S1 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP plasma concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/figure_S1a.tiff", units = "mm", width = 300, height = 270, res = 300)
corrplot(results_descriptive$danish$POPs_heatmap_danish_cases, 
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
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/figure_S1b.tiff", units = "mm", width = 130, height = 120, res = 300)
corrplot(results_descriptive$danish$POPs_heatmap_danish_cases_group, 
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

# Figure S2 - Sensitivity analysis - POPs - ALS survival among the Danish cohort ----
# Association between pre-diagnostic POP groups and survival among ALS cases from the Danish Diet, Cancer and Health cohort (Cox regression modelsl; n = 166).
figure_S2 <- results_POPs_ALS_survival$sensi4$POPs_sd_ALS_figure_danish_sensi_4


# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/table_1.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/table_S2.docx")

table_S3 <- read_docx() |> body_add_flextable(table_S3)
print(table_S3, target = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/table_S3.docx")

table_S4 <- read_docx() |> body_add_flextable(table_S4)
print(table_S4, target = "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/table_S4.docx")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/figure_1.tiff",
  figure_1,
  height = 4,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/figure_2.tiff",
  figure_2,
  height = 5,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/2.Article_POPs_ALS_survival/figure_S2.tiff",
  figure_S2,
  height = 5,
  width = 8,
  units = "in")


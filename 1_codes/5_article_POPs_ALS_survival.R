# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")


# Table 1 - description of the subjects ----
# Description of the subject characteristics of the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263)
table_1 <- results_descriptive$comp$covar_comp_cases_2cat |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Figure 1 - description of the POP levels (boxplots) ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
figure_1 <- results_descriptive$comp$POPs_boxplot_comp_cases_2cat

# Figure 2 - forest plot expo - ALS survival (danish EPIC cohort) ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the Danish Diet, Cancer and Health cohort (cox models by exposure quartiles; n = 166).
figure_2 <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(!term == "Continuous") |>
  filter(!model == "copollutant") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

# Figure 3 - POPs - ALS survival among the Danish cohort (mixture model) ----
# Association between pre-diagnostic POP mixture and survival among ALS cases from the Danish Diet, Cancer and Health cohort (ridge model; n = 166).
figure_3 <- results_POPs_ALS_survival$sensi5$POPs_quart_ALS_figure_sensi5_danish

# Table S1 - description of the POP levels (table) ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
table_S1 <- results_descriptive$comp$POPs_table_comp_cases_2cat |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  |>
  add_footer_lines("pg/ml.
                   POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣDDT corresponds to p,p’-DDT and p,p’-DDE; Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.")

# Table S2 - POPs - ALS survival among the Danish cohort ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the Danish Diet, Cancer and Health cohort (cox models by exposure quartiles; n = 166).
quartile1_rows <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '')

table_S2 <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend") |>
  mutate(across(everything(), as.character))

table_S2 <- 
  bind_rows(quartile1_rows, table_S2) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  group_by(explanatory) |>
  mutate(p.value_heterogeneity_base = ifelse(term == 'quartile 1', p.value_heterogeneity_base[term == 'quartile 2'], ''), 
         p.value_trend_base = ifelse(term == 'quartile 1', p.value_trend_base[term == 'quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(term == 'quartile 1', p.value_heterogeneity_adjusted[term == 'quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(term == 'quartile 1', p.value_trend_adjusted[term == 'quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Heterogeneity test" = "p.value_heterogeneity_base", "Trend test" = "p.value_trend_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted",  "Heterogeneity test " = "p.value_heterogeneity_adjusted", "Trend test " = "p.value_trend_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
  "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.
  2All models are adjusted for sex and age at diagnosis. Adjusted models further account for smoking, BMI and marital status.
  3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.
  5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex and age at diagnosis.
  6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex and age at diagnosis.
  7Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
  8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", 
    term = "Quartiles",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",  "Trend test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model",  "Heterogeneity test " = "Adjusted model",  "Trend test " = "Adjusted model") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  merge_v(j = "term") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  align(j = "term", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "term", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(quartile1_rows)
 
# Table S3 - POPs - ALS survival among the Finnish cohorts ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the Finnish Health Surveys (cox models by exposure quartiles, n = 97).
quartile1_rows <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '')

table_S3 <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
  filter(study_design == "raw unadjusted on cohort") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend") |>
  mutate(across(everything(), as.character))

table_S3 <- 
  bind_rows(quartile1_rows, table_S3) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  group_by(explanatory) |>
  mutate(p.value_heterogeneity_base = ifelse(term == 'quartile 1', p.value_heterogeneity_base[term == 'quartile 2'], ''), 
         p.value_trend_base = ifelse(term == 'quartile 1', p.value_trend_base[term == 'quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(term == 'quartile 1', p.value_heterogeneity_adjusted[term == 'quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(term == 'quartile 1', p.value_trend_adjusted[term == 'quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Heterogeneity test" = "p.value_heterogeneity_base", "Trend test" = "p.value_trend_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted",  "Heterogeneity test " = "p.value_heterogeneity_adjusted", "Trend test " = "p.value_trend_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels_finnish), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2All models are adjusted for sex and age at diagnosis. Adjusted models further account for smoking, BMI and marital status.
    3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.
    5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex and age at diagnosis.
    6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex and age at diagnosis.
    7Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", 
    term = "Quartiles",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",  "Trend test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model",  "Heterogeneity test " = "Adjusted model",  "Trend test " = "Adjusted model") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  merge_v(j = "term") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  align(j = "term", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "term", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(quartile1_rows)

# Table S4 - POPs - ALS survival among the Danish cohort (mixture model) ----
# Association between pre-diagnostic POP mixture and survival among ALS cases from the Danish Diet, Cancer and Health cohort (ridge model; n = 166).
table_S4 <- results_POPs_ALS_survival$sensi5$POPs_quart_ALS_table_sensi5_danish

# Figure S1 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP plasma concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S1a.tiff", units = "mm", width = 300, height = 270, res = 300)
corrplot(results_descriptive$comp$POPs_heatmap_cases, 
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
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S1b.tiff", units = "mm", width = 130, height = 120, res = 300)
corrplot(results_descriptive$comp$POPs_heatmap_cases_group, 
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

# Figure S2 - forest plot expo - ALS survival (finnish cohorts) ----
# Association between pre-diagnostic POP concentrations and survival among ALS cases from the FMC and the FMCF Finnish cohorts (total sample size=86).
figure_S2 <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
  filter(!model == "copollutant") |>
  filter(study_design == "raw unadjusted on cohort") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

# Figure S3 - Sensitivity analysis - distribution of baseline age and follow-up among the cohorts ----
figure_S3 <-
  results_POPs_ALS_survival$sensi1$plot_justif_2cat

# Figure S4 - Sensitivity analysis - effect of follow up duration and baseline age among the Finnish cohorts ----
figure_S4 <- 
  results_POPs_ALS_survival$sensi2$results_sensi2 |>
  filter(box_adj %in% c("not adjusted", "in_orange", "in_purple", "in_green")) |>
  mutate(
    box_adj =  fct_relevel(box_adj, "not adjusted", "in_orange", "in_purple", "in_green"),
    box_adj = fct_recode(
      box_adj,
      "Main analysis\n(n=97)" = "not adjusted",
      "Stratified to baseline\nage > 40 years &\nfollow-up < 350\nmonths (n=41)" = "in_orange",
      "Stratified to follow-up\n< 300 months\n(n=42)" = "in_purple",
      "Stratified to baseline\nage > 45 years\n(n=38)" = "in_green"),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid( ~ box_adj)


# Figure S5 - Sensitivity analysis - POPs - ALS survival among the Danish cohort (mixture model) ----
# Association between pre-diagnostic POP mixture and survival among ALS cases from the Danish Diet, Cancer and Health cohort (elastic net model; n = 166).
figure_S5 <- results_POPs_ALS_survival$sensi6$POPs_quart_ALS_figure_sensi6_danish

# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_1.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_S2.docx")

table_S3 <- read_docx() |> body_add_flextable(table_S3)
print(table_S3, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_S3.docx")

table_S4 <- read_docx() |> body_add_flextable(table_S4)
print(table_S4, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_S4.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_1.tiff",
  figure_1,
  height = 4,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_2.tiff",
  figure_2,
  height = 8,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_3.tiff",
  figure_3,
  height = 3,
  width = 4,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S2.tiff",
  figure_S2,
  height = 8,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S3.tiff",
  figure_S3,
  height = 8, 
  width = 12,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S4.tiff",
  figure_S4,
  height = 5,
  width = 7,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S5.tiff",
  figure_S5,
  height = 6,
  width = 5,
  units = "in")

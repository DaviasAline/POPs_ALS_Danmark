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
                               all_categorical() ~ "{n} ({p})"),
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 0, 
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
figure_1 <-   bdd_danish |>                                                     # densityplot all prot
  filter(match != 159) |>
  select(proteomic_neuro_explo_NEFL, als) |>
  mutate(als = factor(as.character(als)), 
         als = fct_recode(als, "Controls (n=330)" = "0", "Cases (n=165)" = "1")) |>
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_density(alpha = 0.4) +
  scale_fill_hue(direction = -1) + 
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Density", fill = "ALS") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        title = element_text(size = 14)) +
  
  
  bdd_danish |>                                                                 # boxplot all prot
  filter(match != 159) |>
  select(proteomic_neuro_explo_NEFL, als) |>
  mutate(als = factor(as.character(als)), 
         als = fct_recode(als, "Controls" = "0", "Cases" = "1")) |>
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_hue(direction = -1) +
  scale_x_continuous(limits = c(1, 5)) +
  labs(x = "Neurofilament light polypeptide (NPX)", fill = "ALS") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.y = element_blank()) +
  
  plot_layout(heights = c(5, 1), ncol = 1)


# Figure 2 - Base and adjusted logistic regressions (ALS risk) ----
figure_2 <- 
  results_proteomic_ALS_occurrence$main$main_results |>
  filter(analysis %in% c("sensi_1", 
                         "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", 
                         "sensi_1_7_female", "sensi_1_7_male"), 
         term == "Continuous", 
         explanatory == "NEFL", 
         model == "adjusted") |> 
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         analysis = fct_recode(analysis, 
             "Main analysis\n(n=495)" = "sensi_1",
             "Follow-up < 5 years\n (n=51)" = "sensi_2",
             #"Filtered to\nfollow-up > 5 years\n (n=444)" = "sensi_1_3",
             "Follow-up\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
             "Follow-up > 14.6 years\n (n=219)" = "sensi_1_3_5", 
             "Females (n=192)" = "sensi_1_7_female", 
             "Males (n=303)" = "sensi_1_7_male"), 
         analysis = fct_relevel(analysis, 
                                "Main analysis\n(n=495)", 
                                "Follow-up < 5 years\n (n=51)", 
                                #"Filtered to\nfollow-up > 5 years\n (n=444)", 
                                "Follow-up\nbetween 5 and 14.6 years\n(n=225)",
                                "Follow-up > 14.6 years\n (n=219)", 
                                "Females (n=192)",
                                "Males (n=303)")) |> 
  ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis),  switch = "y") +   
  labs( y = "Odd Ratios") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()



# Figure 3 - Additional analysis (NfL level over follow-up time) ----
figure_3 <- results_proteomic_ALS_occurrence$additional_analysis_2$figure_NEFL_over_time_sensi_1

# Figure 4 - Additional analysis (AUC) ----
figure_4 <- results_proteomic_ALS_occurrence$additional_analysis_4$additional_analysis_4_figure_unadjusted
figure_4 <- figure_4 + theme(legend.direction = "vertical") + labs(title = "")

# Figure 5 - Base and adjusted Cox regressions (ALS survival) ----
figure_5 <- results_proteomic_ALS_survival$main$main_results |> 
  filter(term == "Continuous", 
         explanatory == "NEFL", 
         analysis %in% c("main", #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", "sensi_1_7_female", "sensi_1_7_male"), 
         model == "adjusted") |> 
  mutate(
         #signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         # model = fct_recode(model,
         #                    "Adjusted models" = "adjusted",
         #                    "Base models" = "base"),
         # model = fct_relevel(model, "Base models", "Adjusted models"),
         analysis = fct_recode(analysis,
                               "Main analyis\n(n=165)" = "main",
                               #"Follow-up > 5 years\n (n=148)" = "sensi_1_3",
                               "Follow-up\nbetween 5 and 14.6 years\n(n=75)" = "sensi_1_3_4",
                               "Follow-up > 14.6 years\n (n=73)" = "sensi_1_3_5",
                               "Females (n=64)" = "sensi_1_7_female",
                               "Males (n=101)" = "sensi_1_7_male"),
         analysis = fct_relevel(analysis,
                                "Main analyis\n(n=165)",
                                #"Follow-up > 5 years\n (n=148)",
                                "Follow-up\nbetween 5 and 14.6 years\n(n=75)",
                                "Follow-up > 14.6 years\n (n=73)", 
                                "Females (n=64)", 
                                "Males (n=101)")) |> 
  ggplot(aes(x = explanatory, y = HR_raw, ymin = lower_CI, ymax = upper_CI)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(
    rows = dplyr::vars(analysis), 
    #cols = dplyr::vars(model), 
    switch = "y") +                         
  labs( y = "Hazard Ratios") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()



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
  results_proteomic_ALS_occurrence$main$main_results |>
  filter(analysis %in% c("sensi_1", "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", "sensi_1_7_female", "sensi_1_7_male"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(analysis = fct_recode(analysis, 
                               "Main analysis (n=495)" = "sensi_1",
                               "Follow-up < 5 years (n=51)" = "sensi_2", 
                               #"Follow-up > 5 years (n=444)" = "sensi_1_3",
                               "Follow-up between 5 and 14.6 years (n=225)" = "sensi_1_3_4",
                               "Follow-up > 14.6 years (n=219)" = "sensi_1_3_5", 
                               "Females (n=192)" = "sensi_1_7_female", 
                               "Males (n=303)" = "sensi_1_7_male"), 
         analysis = fct_relevel(analysis, 
                                "Main analysis (n=495)", 
                                "Follow-up < 5 years (n=51)", 
                                #"Follow-up > 5 years (n=444)", 
                                "Follow-up between 5 and 14.6 years (n=225)",
                                "Follow-up > 14.6 years (n=219)", 
                                "Females (n=192)", 
                                "Males (n=303)")) |> 
  select(analysis, model, explanatory, OR, "95% CI", p_value) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 'OR' = 'OR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'OR ' = 'OR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and birth year. Adjusted models further account for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "analysis" = "Analyses",
    "explanatory" = "Explanatory variable", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  merge_h(part = "header") |>
  theme_vanilla() |>
  bold(j = "analysis", part = "body") |>
  bold(j = "explanatory", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "analysis", align = "left", part = "all") |> 
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "analysis", part = "header") |>
  merge_at(j = "explanatory", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# Table S3 ----
table_S3 <- 
  results_proteomic_ALS_survival$main$main_results |>
  filter(term == "Continuous", 
         explanatory == "NEFL", 
         analysis %in% c("main", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", "sensi_1_7_female", "sensi_1_7_male")) |> 
  select(analysis, model, explanatory, HR, "95% CI", p_value) |> 
  mutate(
    analysis = fct_recode(analysis,
                          "Main analyis (n=165)" = "main",
                          "Follow-up > 5 years (n=148)" = "sensi_1_3",
                          "Follow-up between 5 and 14.6 years (n=75)" = "sensi_1_3_4",
                          "Follow-up > 14.6 years (n=73)" = "sensi_1_3_5",
                          "Females (n=64)" = "sensi_1_7_female",
                          "Males (n=101)" = "sensi_1_7_male"),
    analysis = fct_relevel(analysis,
                           "Main analyis (n=165)",
                           "Follow-up > 5 years (n=148)",
                           "Follow-up between 5 and 14.6 years (n=75)",
                           "Follow-up > 14.6 years (n=73)", 
                           "Females (n=64)", 
                           "Males (n=101)")) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(HR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 
         'HR' = 'HR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'HR ' = 'HR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and birth year. Adjusted models further account for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Explanatory variable",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  merge_h(part = "header") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


# Export ----
## tables ----
table_1 <- read_docx() |> body_add_flextable(table_1)                           # covariate descriptive table 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_1.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)                         # NEFL descriptive table 
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)                         # Conditional logistic regressions (als risk)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_S2.docx")

table_S3 <- read_docx() |> body_add_flextable(table_S3)                         # Cox regressions (als survival)
print(table_S3, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_S3.docx")


## figures ----
ggsave(                                                                         # NEFL descriptive figure 
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_1.tiff",
  figure_1,
  height = 5,
  width = 10,
  units = "in")

ggsave(                                                                         # Forest plots of logistic regressions results (ALS risk)
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_2.tiff",
  figure_2,
  height = 4,
  width = 7,
  units = "in")

ggsave(                                                                         # LOESS curve of NEFL ratios depending on time to diagnosis
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_3.tiff",
  figure_3,
  height = 5,
  width = 10,
  units = "in")

ggsave(                                                                         # AUC curve
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_4.tiff",
  figure_4,
  height = 10,
  width = 7,
  units = "in")

ggsave(                                                                         # Forest plot of Cox regression results (ALS survival)
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_5.tiff",
  figure_5,
  height = 4,
  width = 7,
  units = "in")




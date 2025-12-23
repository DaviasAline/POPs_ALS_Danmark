# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.7_analyses_proteomic_ALS_occurrence.R")


# Table 1 - Subject characteristics description ---- 
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- 
  bdd_danish |>
  filter(match != 159) |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"),
    follow_up = follow_up/12) |>
  select(
    als, birth_year, baseline_age, diagnosis_age, follow_up,  
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, fS_Kol, 
    proteomic_neuro_explo_NEFL) |>
  tbl_summary(by = als, 
              missing = 'no', 
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
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")




# Figure 1 - Descriptive figure of NEFL distribution (density and boxplots) ----
figure_1 <-   bdd_danish |>                                                                 # densityplot all prot
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


# Figure 2 - Base and adjusted logistic regressions (als risk) ----

figure_2 <- results_proteomic_ALS_occurrence$main$main_results |>
  filter(analysis %in% c("sensi_1", "sensi_2", "sensi_1_3", "sensi_1_3_4", "sensi_1_3_5"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted models" = "adjusted",
                            "Base models" = "base"), 
         model = fct_relevel(model, "Base models", "Adjusted models"), 
         analysis = fct_recode(analysis, 
             "Main analyis\n(n=495)" = "sensi_1",
             "Filtered to\nfollow-up < 5 years\n (n=51)" = "sensi_2", 
             "Filtered to\nfollow-up > 5 years\n (n=447)" = "sensi_1_3",
             "Filtered to follow-up\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
             "Filtered to\nfollow-up > 14.6 years\n (n=227)" = "sensi_1_3_5"), 
         analysis = fct_relevel(analysis, 
                                "Main analyis\n(n=495)", 
                                "Filtered to\nfollow-up < 5 years\n (n=51)",
                                "Filtered to\nfollow-up > 5 years\n (n=447)", 
                                "Filtered to follow-up\nbetween 5 and 14.6 years\n(n=225)",
                                "Filtered to\nfollow-up > 14.6 years\n (n=227)")) |> 
  ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI, color = signif)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis), cols = dplyr::vars(model), switch = "y") +                         # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs( y = "Odd ratios (ORs)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()
figure_2



# Figure 3 - Additional analysis (NEFL level over follow-up time) ----
figure_4 <- results_proteomic_ALS_occurrence$additional_analysis_2$figure_NEFL_over_time_sensi_1

# Figure 4 - Additional analysis (AUC) ----


# Figure 5 - Base and adjusted GAMs (ALS survival) ----


# Export ----
## tables ----
table_1 <- read_docx() |> body_add_flextable(table_1)                           # covariate descriptive table 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_1.docx")

table_2 <- read_docx() |> body_add_flextable(table_2)                           # logistic regressions results (ALS risk) 
print(table_2, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_2.docx")

table_3 <- read_docx() |> body_add_flextable(table_3)                           # Cox regressions results (ALS survival) 
print(table_3, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_3.docx")


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
  height = 5,
  width = 10,
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
  height = 4,
  width = 8,
  units = "in")

ggsave(                                                                         # Forest plot of Cox regression results (ALS survival)
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_5.tiff",
  figure_5,
  height = 5,
  width = 10,
  units = "in")



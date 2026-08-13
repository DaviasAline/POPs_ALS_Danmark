# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

results_proteomic_ALS_occurrence <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.1_results_proteomic_ALS_occurrence.rds")
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
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 0, 
                            bmi ~ 1)) |>
  bold_labels() |>
  add_p(include = -c('diagnosis_age', 'follow_up')) |>
  add_overall() |>
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


# Figure 1 -  Main analysis (volcanoplot overall and time-to-diag stratified)----
figure_1 <- results_proteomic_ALS_occurrence$main$main_results |> 
  filter(term == "Continuous" &
         model == "adjusted" &
         analysis %in% c("sensi_1", "sensi_2", "sensi_1_3_4", "sensi_1_3_5")) |> 
  mutate(
    Significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05"), 
    analysis = fct_recode(analysis, 
                          "<b>All cases and controls</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=495)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>" = "sensi_1",
                          "<b>Years to ALS diagnostic < 5 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=51)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>" = "sensi_2",
                          "<b>Years to ALS diagnostic between 5 and 14.6 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=225)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>" = "sensi_1_3_4",
                          "<b>Years to ALS diagnostic > 14.6 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=219)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>" = "sensi_1_3_5"), 
    analysis = fct_relevel(analysis, 
                           "<b>All cases and controls</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=495)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>", 
                           "<b>Years to ALS diagnostic < 5 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=51)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>", 
                           "<b>Years to ALS diagnostic between 5 and 14.6 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=225)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>",  
                           "<b>Years to ALS diagnostic > 14.6 years</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=219)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>")) |>
  ggplot(aes(x = OR_raw, y = -log10(p_value_raw), color = Significance)) +
  geom_point(alpha = 0.8, size = 2, show.legend = FALSE) +
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
  facet_wrap(vars(analysis), ncol = 2L, nrow = 3L) +
  theme(
    strip.text = element_markdown(hjust = 0, size = 12, lineheight = 1.1),
    strip.background = element_blank(), 
    panel.spacing.y = unit(1.5, "lines"))



# Figure 2 - Prediction - Pre-disease protein levels and risk of developping ALS ----
## Figure 2 a ----
figure_2a <- 
  autoplot(results_proteomic_ALS_occurrence_tidymodels$test_1$results_t1, metric = "roc_auc") +
  geom_text(aes(label = round(mean, 3)), hjust = -0.1, size = 3) +
  labs(title = "a) Algorithm predictive performance comparison", 
       subtitle = "10-fold cross-validation - mean AUC", 
       y = "Mean AUC", 
       x = NULL, 
       color = "Algorithms") +
  scale_color_discrete(
    labels = c(
      "boost_tree" = "Gradient Boosted Trees",
      "logistic_reg" = "Penalized Logistic Regression",
      "mars" = "Multivariate Adaptive Regression Splines",
      "rand_forest" = "Random Forest",
      "svm_linear" = "Linear Support Vector Machine")) +
  guides(
    shape = "none",  
    color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(0, 6, by = 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10), 
    axis.text.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    legend.position = "bottom",          
    legend.direction = "vertical", 
    legend.key.spacing.x = unit(0, "cm"))

## Figure 2 b ----
figure_2b <- results_proteomic_ALS_occurrence_tidymodels$test_1$roc_curve_data_t1 +  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + # plot
  labs(
    title = "b) ROC AUC - Penalized logistic regression",
    subtitle = paste0("10-fold cross-validation - mean AUC: ", 
                      round(results_proteomic_ALS_occurrence_tidymodels$test_1$best_auc_value, 3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10), 
    legend.position = "none")

## Figure 2 c ----
figure_2c <- results_proteomic_ALS_occurrence_tidymodels$test_1$glmnet_t1_result$coefs_protein |> 
  slice_head(n = 50) |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_metabolism_|proteomic_immun_res_"), 
         Feature = fct_reorder(Feature, abs_coef)) |>
  ggplot(aes(x = coefficient, y = Feature, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("↑ ALS risk" = "firebrick", "↓ ALS risk" = "steelblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "c) Top 50 of selected proteins by penalized logistic regression", 
       subtitle = "Total number of selected proteins: 89",
       x = "Coefficients", 
       y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10), 
    legend.position = "bottom")

figure_2 <- (figure_2a / figure_2b) | figure_2c 
figure_2 <- figure_2 + plot_layout(widths = c(4,6)) 
rm(figure_2a, figure_2b, figure_2c)

# Figure 3 - Prediction - Pre-disease protein levels with follow-up interaction hard coded in glmnet ----
results_proteomic_ALS_occurrence_tidymodels$test_3


## Figure 3 a ----
figure_3a <- results_proteomic_ALS_occurrence_tidymodels$test_3$roc_curve_data_t3 +  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + # plot
  labs(
    title = "a) ROC AUC - Penalized logistic regression",
    subtitle = paste0("10-fold cross-validation - mean AUC: ", 
                      round(results_proteomic_ALS_occurrence_tidymodels$test_3$best_auc_value, 3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10), 
    legend.position = "none")

## Figure 3 b ----
figure_3b <- results_proteomic_ALS_occurrence_tidymodels$test_3$glmnet_t3_result$coefs_protein |> 
  slice_head(n = 50) |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_"), 
         Feature = str_replace(Feature, "_x_follow_up_no_na_y", " x follow-up")) |>
  ggplot(aes(x = coefficient, y = Feature, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("↑ ALS risk" = "firebrick", "↓ ALS risk" = "steelblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "b) Top 50 of selected proteins by penalized logistic regression", 
       subtitle = "Total number of selected proteins: ?",
       x = "Coefficients", 
       y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10), 
    legend.position = "bottom")

figure_3 <- figure_3a + figure_3b
rm(figure_3a, figure_3b)


# Supplementary table S1 - Main results ----
table_S1 <- 
  results_proteomic_ALS_occurrence$main$main_results |> 
  filter(term == "Continuous") |>
  filter(model == "adjusted") |>
  filter(analysis %in% c("sensi_1", "sensi_2", "sensi_1_3_4", "sensi_1_3_5")) |>
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
    OR_sensi_1_3_5, `95% CI_sensi_1_3_5`, p_value_sensi_1_3_5) |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", 
         "OR " = "OR_sensi_2", "95% CI " = "95% CI_sensi_2", "p-value " = "p_value_sensi_2", 
         " OR" = "OR_sensi_1_3_4", " 95% CI" = "95% CI_sensi_1_3_4", " p-value" = "p_value_sensi_1_3_4",
         " OR " = "OR_sensi_1_3_5", " 95% CI " = "95% CI_sensi_1_3_5", " p-value " = "p_value_sensi_1_3_5") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex, and adjusted for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Main analyses\nAll cases and controls (N=495)", 
    "95% CI" = "Main analyses\nAll cases and controls (N=495)", 
    "p-value" = "Main analyses\nAll cases and controls (N=495)", 
    "OR " = "Filtered to years to ALS diagnosis < 5 years (N=51)", 
    "95% CI " = "Filtered to years to ALS diagnosis < 5 years (N=51)", 
    "p-value " = "Filtered to years to ALS diagnosis < 5 years (N=51)", 
    " OR" = "Filtered to years to ALS diagnosis between 5 and 14.6 years (N=225)", 
    " 95% CI" = "Filtered to years to ALS diagnosis between 5 and 14.6 years (N=225)", 
    " p-value" = "Filtered to years to ALS diagnosis between 5 and 14.6 years (N=225)", 
    " OR " = "Filtered to years to ALS diagnosis > 14.6 years (N=219)", 
    " 95% CI " = "Filtered to years to ALS diagnosis > 14.6 years (N=219)", 
    " p-value " = "Filtered to years to ALS diagnosis > 14.6 years (N=219)") |>
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



# Supplementary table S2 - Sensitivity results (fixed 5-year time-to-diagnosis windows and sex stratified) ----
table_S2 <- 
  bind_rows(results_proteomic_ALS_occurrence$main$main_results |> 
               filter(term == "Continuous") |>
               filter(model == "adjusted") |>
               filter(analysis %in% c("sensi_1", "sensi_1_7_male", "sensi_1_7_female")), 
            results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results |> 
              filter(term == "Continuous") |>
              filter(model == "adjusted") |>
              filter(analysis %in% c("main_y_5_10ans", "main_y_10_15ans", "main_y_15_20ans"))) |>
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(analysis, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value")) |>
  select(
    protein_group, explanatory,
    OR_sensi_1, `95% CI_sensi_1`, p_value_sensi_1,
    OR_main_y_5_10ans, `95% CI_main_y_5_10ans`, p_value_main_y_5_10ans,
    OR_main_y_10_15ans, `95% CI_main_y_10_15ans`, p_value_main_y_10_15ans,
    OR_main_y_15_20ans, `95% CI_main_y_15_20ans`, p_value_main_y_15_20ans,
    OR_sensi_1_7_male, `95% CI_sensi_1_7_male`, p_value_sensi_1_7_male,
    OR_sensi_1_7_female, `95% CI_sensi_1_7_female`, p_value_sensi_1_7_female) |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", 
         "OR " = "OR_main_y_5_10ans", "95% CI " = "95% CI_main_y_5_10ans", "p-value " = "p_value_main_y_5_10ans", 
         " OR" = "OR_main_y_10_15ans", " 95% CI" = "95% CI_main_y_10_15ans", " p-value" = "p_value_main_y_10_15ans",
         " OR " = "OR_main_y_15_20ans", " 95% CI " = "95% CI_main_y_15_20ans", " p-value " = "p_value_main_y_15_20ans", 
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
    "OR" = "Main analyses\nAll cases and controls (N=495)", 
    "95% CI" = "Main analyses\nAll cases and controls (N=495)", 
    "p-value" = "Main analyses\nAll cases and controls (N=495)", 
    
    "OR " = "Sensitivity analyses\nCases diagnosed within 5 to 10 years (27 cases - 468 controls)", 
    "95% CI " = "Sensitivity analyses\nCases diagnosed within 5 to 10 years (27 cases - 468 controls)", 
    "p-value " = "Sensitivity analyses\nCases diagnosed within 5 to 10 years (27 cases - 468 controls)", 
    " OR" = "Sensitivity analyses\nCases diagnosed within 10 to 15 years (54 cases - 441 controls)", 
    " 95% CI" = "Sensitivity analyses\nCases diagnosed within 10 to 15 years (54 cases - 441 controls)", 
    " p-value" = "Sensitivity analyses\nCases diagnosed within 10 to 15 years (54 cases - 441 controls)", 
    " OR " = "Sensitivity analyses\nCases diagnosed within 15 to 20 years (57 cases - 438 controls)", 
    " 95% CI " = "Sensitivity analyses\nCases diagnosed within 15 to 20 years (57 cases - 438 controls)", 
    " p-value " = "Sensitivity analyses\nCases diagnosed within 15 to 20 years (57 cases - 438 controls)", 
    "OR  " = "Sensitivity analyses\nFiltered to males (N=303)", 
    "95% CI  " = "Sensitivity analyses\nFiltered to males (N=303)", 
    "p-value  " = "Sensitivity analyses\nFiltered to males (N=303)", 
    " OR  " = "Sensitivity analyses\nFiltered to females (N=192)", 
    " 95% CI  " = "Sensitivity analyses\nFiltered to females (N=192)", 
    " p-value  " = "Sensitivity analyses\nFiltered to females (N=192)") |>
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



# Supplementary figure 1 - Sensitivity analysis (volcanoplot fixed 5-year time-to-diagnosis windows) ----
figure_S1 <- 
  bind_rows(results_proteomic_ALS_occurrence$main$main_results |> filter(term == "Continuous") |> filter(model == "adjusted") |> filter(analysis == "sensi_1"), 
            results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results |> filter(term == "Continuous") |> filter(model == "adjusted") |> filter(analysis == "main_y_5_10ans"), 
            results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results |> filter(term == "Continuous") |> filter(model == "adjusted") |> filter(analysis == "main_y_10_15ans"), 
            results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results |> filter(term == "Continuous") |> filter(model == "adjusted") |> filter(analysis == "main_y_15_20ans")) |> 
  mutate(
    Significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05"), 
    analysis = fct_recode(analysis, 
                          "<b>All cases and controls</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=495)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>" = "sensi_1",
                          "<b>Cases diagnosed within 5 to 10 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (27 cases - 468 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>" = "main_y_5_10ans",
                          "<b>Cases diagnosed within 10 to 15 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (54 cases - 441 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>" = "main_y_10_15ans",
                          "<b>Cases diagnosed within 15 to 20 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (57 cases - 438 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>" = "main_y_15_20ans"), 
    analysis = fct_relevel(analysis, 
                           "<b>All cases and controls</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=495)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>", 
                           "<b>Cases diagnosed within 5 to 10 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (27 cases - 468 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>", 
                           "<b>Cases diagnosed within 10 to 15 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (54 cases - 441 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>",  
                           "<b>Cases diagnosed within 15 to 20 years</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (57 cases - 438 controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>")) |>
  ggplot(aes(x = OR_raw, y = -log10(p_value_raw), color = Significance)) +
  geom_point(alpha = 0.8, size = 2, show.legend = FALSE) +
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
  facet_wrap(vars(analysis), ncol = 2L, nrow = 3L) +
  theme(
    strip.text = element_markdown(hjust = 0, size = 12, lineheight = 1.1),
    strip.background = element_blank(), 
    panel.spacing.y = unit(1.5, "lines"))


# Supplementary figure 2 - Sensitivity analysis (vocanoplot filtered on sex) ----
figure_S2 <-
  results_proteomic_ALS_occurrence$main$main_results |>
  filter(model == "adjusted" &                                                  # select only adjusted results
           term == "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1", "sensi_1_7_female", "sensi_1_7_male")) |>  
  mutate(
    Significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05"), 
    analysis = fct_recode(analysis, 
                          "<b>All case and controls (N=495)</b>" = "sensi_1",
                          "<b>Females (N=192)</b>" = "sensi_1_7_female",
                          "<b>Males (N=303)</b>" = "sensi_1_7_male"), 
    analysis = fct_relevel(analysis, 
                           "<b>All case and controls (N=495)</b>", 
                           "<b>Females (N=192)</b>",  
                           "<b>Males (N=303)</b>")) |>
  ggplot(aes(x = OR_raw, y = -log10(p_value_raw), color = Significance)) +
  geom_point(alpha = 0.8, size = 2, show.legend = FALSE) +
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
  facet_wrap(vars(analysis), ncol = 3L, nrow = 1L) +
  theme(
    strip.text = element_markdown(hjust = 0, size = 12, lineheight = 1.1),
    strip.background = element_blank(), 
    panel.spacing.y = unit(1.5, "lines"))


# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_1.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_1.tiff",
  figure_1,
  height = 10,
  width = 13, 
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_2.tiff",
  figure_2,
  height = 8,
  width = 10, 
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_3.tiff",
  figure_3,
  height = 8,
  width = 10, 
  units = "in")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_S1.docx")

table_S2 <- read_docx() |> body_add_flextable(table_S2)
print(table_S2, target = "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/table_S2.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S1.tiff",
  figure_S1,
  height = 10,
  width = 13, 
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/4.Article_proteomics_ALS_occurence/figure_S2.tiff",
  figure_S2,
  height = 7,
  width = 13,
  units = "in")


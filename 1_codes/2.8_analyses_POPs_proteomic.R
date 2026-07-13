# Aline Davias
# July 08, 2026
# Analysis of proteomic levels depending on POPs exposure

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# Préparation des données ----
bdd_danish <- bdd_danish |>
  mutate(across(all_of(c(POPs_group, POPs_included)),                           # creation of cohort and case specific standardized exposures
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  

# Préparation des combinaisons de variables ----
analysis_grid <- expand_grid(
  exposure = POPs_group_sd,
  outcome = proteomic)

# Automatisation des modèles (Régressions Linéaires) ----
main_results <- analysis_grid |>
  mutate(
    # base model main             
    base_model = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + as.factor(als) + baseline_age + sex"))
      lm(formule, data = bdd_danish) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "main", 
               model = "base")}),
    
    # adjusted model main
    adjusted_model = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, "  + as.factor(als) + baseline_age + sex + smoking_2cat_i + bmi"))
      lm(formule, data = bdd_danish) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "main", 
               model = "adjusted")}), 
    
    # base model sensi sex female 
    base_model_sensi_sex_f = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + as.factor(als) + baseline_age"))
      lm(formule, data = bdd_danish |> filter(sex == "Female")) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_sex_f", 
               model = "base")}),
    
    # adjusted model sensi sex female 
    adjusted_model_sensi_sex_f = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, "  + as.factor(als) + baseline_age + smoking_2cat_i + bmi"))
      lm(formule, data = bdd_danish |> filter(sex == "Female")) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_sex_f", 
               model = "adjusted")}), 
    
    # base model sensi sex male 
    base_model_sensi_sex_m = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + as.factor(als) + baseline_age"))
      lm(formule, data = bdd_danish |> filter(sex == "Male")) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_sex_m", 
               model = "base")}),
    
    # adjusted model sensi sex male 
    adjusted_model_sensi_sex_m = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, "  + as.factor(als) + baseline_age + smoking_2cat_i + bmi"))
      lm(formule, data = bdd_danish |> filter(sex == "Male")) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_sex_m", 
               model = "adjusted")}), 
    
    # base model sensi als cases 
    base_model_sensi_als_cases = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + baseline_age + sex"))
      lm(formule, data = bdd_danish |> filter(als == 1)) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_als_cases", 
               model = "base")}),
    
    # adjusted model sensi als cases 
    adjusted_model_sensi_als_cases = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + baseline_age + sex + smoking_2cat_i + bmi"))
      lm(formule, data = bdd_danish |> filter(als == 1)) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_als_cases", 
               model = "adjusted")}), 
    
    
    # base model sensi als controls  
    base_model_sensi_als_controls = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + baseline_age + sex"))
      lm(formule, data = bdd_danish |> filter(als == 0)) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_als_controls", 
               model = "base")}),
    
    # adjusted model sensi als controls 
    adjusted_model_sensi_als_controls = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + baseline_age + sex + smoking_2cat_i + bmi"))
      lm(formule, data = bdd_danish |> filter(als == 0)) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "sensi_als_controls", 
               model = "adjusted")}))
rm(analysis_grid)

# Mise en forme et calcul du False Discovery Rate (FDR) ----
main_results <- main_results |>
  pivot_longer(
    cols = c(base_model, adjusted_model, 
             base_model_sensi_sex_f, adjusted_model_sensi_sex_f, 
             base_model_sensi_sex_m, adjusted_model_sensi_sex_m, 
             base_model_sensi_als_cases, adjusted_model_sensi_als_cases, 
             base_model_sensi_als_controls, adjusted_model_sensi_als_controls), 
    names_to = NULL, 
    values_to = "model_summary") |>
  unnest(model_summary) |>
  # Calcul de la q-value (FDR) par type de modèle
  mutate(fdr_group = case_when(
    analysis == "main" ~ "main",
    analysis %in% c("sensi_sex_f", "sensi_sex_m") ~ "sensi_sex",
    analysis %in% c("sensi_als_cases", "sensi_als_controls") ~ "sensi_als")) |>
  group_by(model, fdr_group) |>
  mutate(q.value_raw = p.adjust(p.value, method = "fdr")) |>
  ungroup() |>
  select(-fdr_group) |> 
  mutate(
    estimate_raw = estimate, 
    estimate = sprintf("%.2f", estimate),
    conf.low = sprintf("%.2f", conf.low),
    conf.high = sprintf("%.2f", conf.high),
    p.value_raw = p.value, 
    "p-value" = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")), 
    "FDR-corrected p-value" = ifelse(q.value_raw < 0.01, "<0.01", number(q.value_raw, accuracy = 0.01, decimal.mark = ".")), 
    "95% CI" = paste(conf.low, ", ", conf.high, sep = ''), 
    is_q_significant = if_else(q.value_raw < 0.05, "FDR-corrected p-value < 0.05", "FDR-corrected p-value ≥ 0.05"),
    is_p_significant = if_else(p.value_raw < 0.05, "p-value < 0.05", "p-value ≥ 0.05"),
    outcome_group = case_when(str_detect(outcome, 'proteomic_immun_res') ~ "Immune response", 
                              str_detect(outcome, 'proteomic_metabolism') ~ "Metabolism", 
                              str_detect(outcome, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
    outcome = str_replace(outcome, "proteomic_metabolism_|proteomic_neuro_explo_|proteomic_immun_res_", ""), 
    exposure_group = fct_recode(exposure, 
                                "OCP" = "OCP_HCB_sd",
                                "OCP" = "OCP_β_HCH_sd",
                                "PCB" = "PCB_4_sd",
                                "PCB" = "PCB_DL_sd",
                                "PCB" = "PCB_NDL_sd",
                                "OCP" = "Σchlordane_sd",
                                "OCP" = "ΣDDT_sd",
                                "PBDE" = "ΣPBDE_sd"), 
    exposure = str_replace(exposure, "_sd", ""), 
    exposure = fct_recode(exposure, 
                          "HCB" = "OCP_HCB",
                          "β-HCH" = "OCP_β_HCH",
                          "PCB-4" = "PCB_4",
                          "PCB-DL" = "PCB_DL",
                          "PCB-NDL" = "PCB_NDL")) |>
  select(
    analysis, 
    model,
    exposure_group, exposure,
    outcome_group, outcome,
    Beta = estimate, estimate_raw,      
    conf.low, conf.high, "95% CI", 
    "p-value", p.value_raw, 
    "FDR-corrected p-value", q.value_raw, 
    is_p_significant, is_q_significant)

# Visualisation ----
# main_results |> 
#   filter(q.value_raw <0.05) |> 
#   filter(model == "adjusted") |> 
#   View()
# main_results |> 
#   filter(outcome == "NEFL") |> 
#   filter(model == "adjusted") |> 
#   filter(p.value_raw <0.05) |> 
#   View()


POPs_sd_proteomic_table <-                                                       # select both base and adjusted results
  main_results |>
  filter(analysis == "main" & 
           model %in% c("base", "adjusted")) |>            
  group_by(exposure, outcome) |>                                                      # select explanatory vars significant                
  filter(any(q.value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, exposure, outcome_group, outcome, Beta, "95% CI", "p-value", "FDR-corrected p-value") |>
  pivot_wider(names_from = "model", values_from = c("Beta", "95% CI", "p-value", "FDR-corrected p-value")) |>
  select(exposure, outcome_group, outcome, contains("base"), contains("adjusted")) |>
  rename("Beta" = "Beta_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "FDR-corrected p-value" = "FDR-corrected p-value_base", 
         "Beta " = "Beta_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "FDR-corrected p-value " = "FDR-corrected p-value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome_group" = "Protein groups", 
    "outcome" = "Proteins", 
    "Beta" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "FDR-corrected p-value" = "Base model", 
    "Beta " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "FDR-corrected p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = 1:3) |>
  bold(j = 1:3, part = "body") |>
  align(j = 1:3, align = "left", part = "all") |> 
  merge_at(j = "outcome", part = "header") |>
  merge_at(j = "exposure", part = "header") |>
  merge_at(j = "outcome_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

POPs_sd_NfL_table <-                                                       # select both base and adjusted results
  main_results |>
  filter(analysis == "main" & 
           outcome == "NEFL" &
           model %in% c("base", "adjusted")) |>    
  mutate(outcome = fct_recode(outcome, "Neurofilament light polypeptide" = "NEFL")) |>        
  select(model, exposure, outcome, Beta, "95% CI", "p-value", "FDR-corrected p-value") |>
  pivot_wider(names_from = "model", values_from = c("Beta", "95% CI", "p-value", "FDR-corrected p-value")) |>
  select(exposure, outcome, contains("base"), contains("adjusted")) |>
  rename("Beta" = "Beta_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "FDR-corrected p-value" = "FDR-corrected p-value_base", 
         "Beta " = "Beta_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "FDR-corrected p-value " = "FDR-corrected p-value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome" = "Protein", 
    "Beta" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "FDR-corrected p-value" = "Base model", 
    "Beta " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "FDR-corrected p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "outcome") |>
  bold(j = "outcome", part = "body") |>
  align(j = "outcome", align = "left", part = "all") |> 
  merge_at(j = "outcome", part = "header") |>
  merge_v(j = "exposure") |>
  bold(j = "exposure", part = "body") |>
  align(j = "exposure", align = "left", part = "all") |> 
  merge_at(j = "exposure", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


POPs_sd_proteomic_figure <- main_results |>
  filter(analysis == "main" & model == "adjusted") |>
  ggplot(aes(x = estimate_raw, y = -log10(p.value_raw))) +
  geom_point(aes(color = is_q_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", alpha = 0.3) +
  geom_text_repel(
    aes(label = outcome),
    size = 2.5,
    max.overlaps = 20,
    segment.color = 'grey50') +
  facet_wrap(~ exposure, scales = "fixed", ncol = 3) +
  scale_color_manual(values = c("FDR-corrected p-value < 0.05" = "firebrick3", 
                                "FDR-corrected p-value ≥ 0.05" = "black")) +
  labs(x = "Beta", y = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank())


POPs_sd_proteomic_table_sensi_sex <-                                                       # select adjusted results
  main_results |>
  filter(analysis %in% c("main", "sensi_sex_f", "sensi_sex_m") & 
           model == "adjusted") |>            
  group_by(exposure, outcome) |>                                                      # select explanatory vars significant                
  filter(any(q.value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(analysis, exposure, outcome_group, outcome, Beta, "95% CI", "p-value", "FDR-corrected p-value") |>
  pivot_wider(names_from = "analysis", values_from = c("Beta", "95% CI", "p-value", "FDR-corrected p-value")) |>
  select(exposure, outcome_group, outcome, contains("main"), contains("sensi_sex_m"), contains("sensi_sex_f")) |>
  rename("Beta" = "Beta_main", "95% CI" = "95% CI_main", "p-value" = "p-value_main", "FDR-corrected p-value" = "FDR-corrected p-value_main", 
         " Beta " = "Beta_sensi_sex_m", " 95% CI " = "95% CI_sensi_sex_m", " p-value " = "p-value_sensi_sex_m", " FDR-corrected p-value " = "FDR-corrected p-value_sensi_sex_m", 
         "Beta " = "Beta_sensi_sex_f", "95% CI " = "95% CI_sensi_sex_f", "p-value " = "p-value_sensi_sex_f", "FDR-corrected p-value " = "FDR-corrected p-value_sensi_sex_f") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome_group" = "Protein groups", 
    "outcome" = "Proteins", 
    "Beta" = "Main analysis (N = 498)", "95% CI" = "Main analysis (N = 498)", "p-value" = "Main analysis (N = 498)", "FDR-corrected p-value" = "Main analysis (N = 498)",
    " Beta " = "Sensitivity analysis\nMales (N = 303)", " 95% CI " = "Sensitivity analysis\nMales (N = 303)", " p-value " = "Sensitivity analysis\nMales (N = 303)", " FDR-corrected p-value " = "Sensitivity analysis\nMales (N = 303)", 
    "Beta " = "Sensitivity analysis\nFemales (N = 195)", "95% CI " = "Sensitivity analysis\nFemales (N = 195)", "p-value " = "Sensitivity analysis\nFemales (N = 195)", "FDR-corrected p-value " = "Sensitivity analysis\nFemales (N = 195)") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = 1:3) |>
  bold(j = 1:3, part = "body") |>
  align(j = 1:3, align = "left", part = "all") |> 
  merge_at(j = "outcome", part = "header") |>
  merge_at(j = "exposure", part = "header") |>
  merge_at(j = "outcome_group", part = "header") |>
  bold(i = ~ as.numeric(gsub("<", "", `FDR-corrected p-value`)) < 0.05, j = "FDR-corrected p-value", part = "body") |>
  bold(i = ~ as.numeric(gsub("<", "", `FDR-corrected p-value `)) < 0.05, j = "FDR-corrected p-value ", part = "body") |>
  bold(i = ~ as.numeric(gsub("<", "", ` FDR-corrected p-value `)) < 0.05, j = " FDR-corrected p-value ", part = "body") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


POPs_sd_proteomic_figure_sensi_sex <- 
  main_results |>
  filter(model == "adjusted" & 
           analysis %in% c("main",  "sensi_sex_m", "sensi_sex_f")) |>
  mutate(analysis = fct_recode(analysis,
                               "Main analysis\n(N = 498)" = "main",
                               "Sensitivity analysis\nMales (N = 303)" = "sensi_sex_m", 
                               "Sensitivity analysis\nFemales (N = 195)" = "sensi_sex_f"), 
         analysis = fct_relevel(analysis, 
                                 "Main analysis\n(N = 498)",
                                 "Sensitivity analysis\nMales (N = 303)",
                                 "Sensitivity analysis\nFemales (N = 195)")) |>
  
  ggplot(aes(x = estimate_raw, y = -log10(p.value_raw))) +
  geom_point(aes(color = is_q_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", alpha = 0.3) +
  geom_text_repel(
    aes(label = outcome),
    size = 2.5,
    max.overlaps = 20,
    segment.color = 'grey50') +
  facet_grid(exposure ~ analysis, scales = "fixed", switch = "y") +
  scale_color_manual(values = c("FDR-corrected p-value < 0.05" = "firebrick3", 
                                "FDR-corrected p-value ≥ 0.05" = "black")) +
  labs(x = "Beta", y = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    strip.text = element_text(face = "bold"), 
    strip.text.y.left = element_text(angle = 0, hjust = 1), # angle = 0 rend le texte horizontal
    strip.placement = "outside")


POPs_sd_proteomic_table_sensi_als <-                                                       # select adjusted results
  main_results |>
  filter(analysis %in% c("main", "sensi_als_cases", "sensi_als_controls") & 
           model == "adjusted") |>            
  group_by(exposure, outcome) |>                                                      # select explanatory vars significant                
  filter(any(q.value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(analysis, exposure, outcome_group, outcome, Beta, "95% CI", "p-value", "FDR-corrected p-value") |>
  pivot_wider(names_from = "analysis", values_from = c("Beta", "95% CI", "p-value", "FDR-corrected p-value")) |>
  select(exposure, outcome_group, outcome, contains("main"), contains("sensi_als_cases"), contains("sensi_als_controls")) |>
  rename("Beta" = "Beta_main", "95% CI" = "95% CI_main", "p-value" = "p-value_main", "FDR-corrected p-value" = "FDR-corrected p-value_main", 
         "Beta " = "Beta_sensi_als_cases", "95% CI " = "95% CI_sensi_als_cases", "p-value " = "p-value_sensi_als_cases", "FDR-corrected p-value " = "FDR-corrected p-value_sensi_als_cases", 
         " Beta " = "Beta_sensi_als_controls", " 95% CI " = "95% CI_sensi_als_controls", " p-value " = "p-value_sensi_als_controls", " FDR-corrected p-value " = "FDR-corrected p-value_sensi_als_controls") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome_group" = "Protein groups", 
    "outcome" = "Proteins", 
    "Beta" = "Main analysis\n(N = 498)", "95% CI" = "Main analysis\n(N = 498)", "p-value" = "Main analysis\n(N = 498)", "FDR-corrected p-value" = "Main analysis\n(N = 498)", 
    "Beta " = "Sensitivity analysis\nCases (N = 166)", "95% CI " = "Sensitivity analysis\nCases (N = 166)", "p-value " = "Sensitivity analysis\nCases (N = 166)", "FDR-corrected p-value " = "Sensitivity analysis\nCases (N = 166)",
    " Beta " = "Sensitivity analysis\nControls (N = 332)", " 95% CI " = "Sensitivity analysis\nControls (N = 332)", " p-value " = "Sensitivity analysis\nControls (N = 332)", " FDR-corrected p-value " = "Sensitivity analysis\nControls (N = 332)") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = 1:3) |>
  bold(j = 1:3, part = "body") |>
  align(j = 1:3, align = "left", part = "all") |> 
  merge_at(j = "outcome", part = "header") |>
  merge_at(j = "exposure", part = "header") |>
  merge_at(j = "outcome_group", part = "header") |>
  bold(i = ~ as.numeric(gsub("<", "", `FDR-corrected p-value`)) < 0.05, j = "FDR-corrected p-value", part = "body") |>
  bold(i = ~ as.numeric(gsub("<", "", `FDR-corrected p-value `)) < 0.05, j = "FDR-corrected p-value ", part = "body") |>
  bold(i = ~ as.numeric(gsub("<", "", ` FDR-corrected p-value `)) < 0.05, j = " FDR-corrected p-value ", part = "body") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


POPs_sd_proteomic_figure_sensi_als <- 
  main_results |>
  filter(model == "adjusted" & 
           analysis %in% c("main", "sensi_als_cases", "sensi_als_controls")) |>
  mutate(
    analysis = fct_recode(analysis,
                          "Main analysis\n(N = 498)" = "main",
                          "Sensitivity analysis\nCases (N = 166)" = "sensi_als_cases",    
                          "Sensitivity analysis\nControls (N = 332)" = "sensi_als_controls"), 
    analysis = fct_relevel(analysis, 
                           "Main analysis\n(N = 498)",
                           "Sensitivity analysis\nCases (N = 166)",
                           "Sensitivity analysis\nControls (N = 332)")) |>
  ggplot(aes(x = estimate_raw, y = -log10(p.value_raw))) +
  geom_point(aes(color = is_q_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", alpha = 0.3) +
  geom_text_repel(
    aes(label = outcome),
    size = 2.5,
    max.overlaps = 20,
    segment.color = 'grey50') +
  facet_grid(exposure ~ analysis, scales = "fixed", switch = "y") +
  scale_color_manual(values = c("FDR-corrected p-value < 0.05" = "firebrick3", 
                                "FDR-corrected p-value ≥ 0.05" = "black")) +
  labs(x = "Beta", y = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y.left = element_text(angle = 0, hjust = 1), 
    strip.placement = "outside")


results_POPs_proteomic <- list(
  main = list(main_results = main_results, 
              POPs_sd_proteomic_table = POPs_sd_proteomic_table, 
              POPs_sd_proteomic_figure = POPs_sd_proteomic_figure, 
              POPs_sd_NfL_table = POPs_sd_NfL_table), 
  sensi_sex = list(POPs_sd_proteomic_table_sensi_sex = POPs_sd_proteomic_table_sensi_sex, 
                   POPs_sd_proteomic_figure_sensi_sex = POPs_sd_proteomic_figure_sensi_sex),
  sensi_als = list(POPs_sd_proteomic_table_sensi_als = POPs_sd_proteomic_table_sensi_als, 
                   POPs_sd_proteomic_figure_sensi_als = POPs_sd_proteomic_figure_sensi_als))

saveRDS(results_POPs_proteomic, file = "~/Documents/POP_ALS_2025_02_03/2_output/2.8_results_POPs_proteomic.rds")

rm(main_results, 
   POPs_sd_proteomic_table, 
   POPs_sd_proteomic_figure, 
   POPs_sd_NfL_table, 
   POPs_sd_proteomic_table_sensi_sex, 
   POPs_sd_proteomic_figure_sensi_sex, 
   POPs_sd_proteomic_table_sensi_als, 
   POPs_sd_proteomic_figure_sensi_als)

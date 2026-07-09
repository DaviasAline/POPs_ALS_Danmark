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
    # base model
    base_model = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + as.factor(match)"))
      lm(formule, data = bdd_danish) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "main", 
               model = "base")}),
    
    # adjusted model
    adjusted_model = map2(exposure, outcome, function(exp, out) {
      formule <- as.formula(paste0(out, " ~ ", exp, " + smoking_2cat_i + bmi + as.factor(match)"))
      lm(formule, data = bdd_danish) |> 
        tidy(conf.int = TRUE) |> 
        filter(term == exp) |> 
        mutate(analysis = "main", 
               model = "adjusted")}))
rm(analysis_grid)

# Mise en forme et calcul du False Discovery Rate (FDR) ----
main_results <- main_results |>
  pivot_longer(
    cols = c(base_model, adjusted_model), 
    names_to = NULL, 
    values_to = "model_summary") |>
  unnest(model_summary) |>
  # Calcul de la q-value (FDR) par type de modèle
  group_by(model) |>
  mutate(q.value_raw = p.adjust(p.value, method = "fdr")) |>
  ungroup() |>
  mutate(
    estimate_raw = estimate, 
    estimate = sprintf("%.2f", estimate),
    conf.low = sprintf("%.2f", conf.low),
    conf.high = sprintf("%.2f", conf.high),
    p.value_raw = p.value, 
    "p-value" = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")), 
    "q-value" = ifelse(q.value_raw < 0.01, "<0.01", number(q.value_raw, accuracy = 0.01, decimal.mark = ".")), 
    "95% CI" = paste(conf.low, ", ", conf.high, sep = ''), 
    is_q_significant = if_else(q.value_raw < 0.05, "FDR < 0.05", "FDR ≥ 0.05"),
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
    "q-value", q.value_raw, 
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
  select(model, exposure, outcome_group, outcome, Beta, "95% CI", "p-value", "q-value") |>
  pivot_wider(names_from = "model", values_from = c("Beta", "95% CI", "p-value", "q-value")) |>
  select(exposure, outcome_group, outcome, contains("base"), contains("adjusted")) |>
  rename("Beta" = "Beta_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "q-value" = "q-value_base", 
         "Beta " = "Beta_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "q-value " = "q-value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome_group" = "Protein groups", 
    "outcome" = "Proteins", 
    "Beta" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "q-value" = "Base model", 
    "Beta " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "q-value " = "Adjusted model") |>
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
  select(model, exposure, outcome, Beta, "95% CI", "p-value", "q-value") |>
  pivot_wider(names_from = "model", values_from = c("Beta", "95% CI", "p-value", "q-value")) |>
  select(exposure, outcome, contains("base"), contains("adjusted")) |>
  rename("Beta" = "Beta_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "q-value" = "q-value_base", 
         "Beta " = "Beta_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "q-value " = "q-value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated increase of plasma protein level (NPX) associated with a one standard deviation increase in pre-disease POP.
    3CI: Confidence interval.") |>
  add_header(
    "exposure" = "Pre-disease exposures", 
    "outcome" = "Protein", 
    "Beta" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "q-value" = "Base model", 
    "Beta " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "q-value " = "Adjusted model") |>
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
  filter(model == "adjusted") |>
  ggplot(aes(x = estimate_raw, y = -log10(p.value_raw))) +
  geom_point(aes(color = is_q_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", alpha = 0.3) +
  geom_text_repel(
    aes(label = outcome),
    size = 2.5,
    max.overlaps = 20,
    segment.color = 'grey50') +
  facet_wrap(~ exposure, scales = "fixed", ncol = 4) +
  scale_color_manual(values = c("FDR < 0.05" = "firebrick3", 
                                "FDR ≥ 0.05" = "black")) +
  labs(x = "Beta", y = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank())


results_POPs_proteomic <- list(
  main = list(main_results = main_results, 
              POPs_sd_proteomic_table = POPs_sd_proteomic_table, 
              POPs_sd_proteomic_figure = POPs_sd_proteomic_figure, 
              POPs_sd_NfL_table = POPs_sd_NfL_table))

saveRDS(results_POPs_proteomic, file = "~/Documents/POP_ALS_2025_02_03/2_output/2.8_results_POPs_proteomic.rds")

rm(main_results, 
   POPs_sd_proteomic_table, 
   POPs_sd_proteomic_figure, 
   POPs_sd_NfL_table)

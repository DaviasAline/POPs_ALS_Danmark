# Chargement des packages requis ----
library(tidyverse)
library(broom)
library(ggrepel) # Pour éviter que les noms des protéines se chevauchent

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
results <- analysis_grid |>
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

# Mise en forme et calcul du False Discovery Rate (FDR) ----
main_results <- results |>
  pivot_longer(
    cols = c(base_model, adjusted_model), 
    names_to = NULL, 
    values_to = "model_summary") |>
  unnest(model_summary) |>
  select(
    analysis, 
    model,
    exposure,
    outcome,
    Beta = estimate,      
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value) |>
  # Calcul de la q-value (FDR) par type de modèle
  group_by(model) |>
  mutate(q_value = p.adjust(p_value, method = "fdr")) |>
  ungroup()

# Visualisation ----
main_results |> 
  filter(q_value <0.05) |> 
  filter(model == "adjusted") |> 
  View()
main_results |> 
  filter(outcome == "proteomic_neuro_explo_NEFL") |> 
  filter(model == "adjusted") |> 
  filter(p_value <0.05) |> 
  View()

plot <- main_results |>
  filter(model == "adjusted") |>
  mutate(
    is_significant = if_else(q_value < 0.05, "FDR < 0.05", "FDR ≥ 0.05"),
    label_proteine = if_else(q_value < 0.05, outcome, ""), 
    label_proteine = str_replace(label_proteine, "proteomic_metabolism_|proteomic_neuro_explo_|proteomic_immun_res_", ""), 
    exposure = str_replace(exposure, "_sd", ""), 
    exposure = fct_recode(exposure, 
        "HCB" = "OCP_HCB",
        "β-HCH" = "OCP_β_HCH",
        "PCB-4" = "PCB_4",
        "PCB-DL" = "PCB_DL",
        "PCB-NDL" = "PCB_NDL")) |>
  ggplot(aes(x = Beta, y = -log10(p_value))) +
  geom_point(aes(color = is_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", alpha = 0.3) +
  geom_text_repel(
    aes(label = label_proteine),
    size = 2.5,
    max.overlaps = 20,
    segment.color = 'grey50') +
  facet_wrap(~ exposure, scales = "fixed", ncol = 4) +
  scale_color_manual(values = c("FDR < 0.05" = "firebrick3", 
                                "FDR ≥ 0.05" = "black")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank())


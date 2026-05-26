# Analyses Years to ALS 
# Aline Davias 
# 20/05/2026
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
library(xgboost)

# 1. Cox analyses among cases ----
## 1.1 Proteome wide associations ----
fit_and_plot_cox_gam <- function(
    data,
    vars = proteomic,
    time = "follow_up",
    status = "als", 
    covariates = c("sex", "baseline_age"),  # default
    vars_labels = NULL,
    y_limits = c(-12, 12),
    title = "Base model"                     # default
) {
  
  # 1. Gestion des labels
  if (is.null(vars_labels)) {
    vars_labels <- set_names(vars, vars)
  }
  
  # Sous-fonction pour traiter une protéine
  fit_one <- function(var) {
    
    outcome <- cbind(data[[time]], data[[status]])
    formula_str <- paste0(
      "outcome ~ s(", var, ") + ",
      paste(covariates, collapse = " + "))
    
    model <- mgcv::gam(as.formula(formula_str),
                       data = data,
                       family = cox.ph())
    
    smry <- summary(model)
    edf_val <- smry$s.table[1, "edf"]
    pval_raw <- smry$s.table[1, "p-value"]
    
    edf_label <- format(edf_val, nsmall = 1, digits = 1)
    pval_label <- case_when(
      pval_raw < 0.01 ~ "< 0.01",
      pval_raw > 0.99 ~ "> 0.99",
      TRUE ~ format(pval_raw, nsmall = 2, digits = 2))
    
    pdf(NULL)
    plot_data_raw <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE)
    dev.off()
    
    if (length(plot_data_raw) == 0) return(NULL)
    
    plot_data <- plot_data_raw[[1]]
    smooth_df <- data.frame(
      x = plot_data$x,
      fit = plot_data$fit,
      se = plot_data$se)
    
    p1 <- ggplot(smooth_df, aes(x = x, y = fit)) +
      geom_line(linewidth = 1.2, color = "steelblue") +
      geom_ribbon(aes(ymin = fit - 2 * se, 
                      ymax = fit + 2 * se),
                  alpha = 0.2, fill = "steelblue") +
      labs(
        title = title,
        y = "LogHR (smooth estimate)",
        x = NULL) +
      annotate(
        "text",
        x = -Inf, y = Inf,
        hjust = -0.1, vjust = 1.5,
        label = paste0("EDF: ", edf_label, "\np-value: ", pval_label),
        size = 4.2) +
      scale_y_continuous(limits = y_limits) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
    
    p2 <- ggplot(data, aes(x = .data[[var]], y = 1)) +
      geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
      labs(x = vars_labels[[var]]) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
    
    list(
      model = model,
      pval_raw = pval_raw,
      plot = p1 / p2 + plot_layout(heights = c(10, 1)))
  }
  
  results <- map(vars, fit_one)
  names(results) <- vars_labels[vars]
  
  return(results)
}

bdd_danish_cases <- bdd_danish |> filter(als == 1)

surv_obj <- Surv(time = bdd_danish_cases$follow_up, event = bdd_danish_cases$als)

## Cox model (sd) ----
### Base ----
model1_cox_sd <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age"))  # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_danish_cases) |> summary()      # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "main", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Adjusted ----
model2_cox_sd <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age + smoking_2cat_i + bmi"))  # set the formulas     
  
  model_summary <- coxph(formula_danish, data = bdd_danish_cases) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "main", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

## Cox model (quart) ----
### Base ----
model1_cox_quart <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age"))  # set the formulas   
  
  model_summary <- coxph(formula_danish, data = bdd_danish_cases) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "main", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Adjusted ----
model2_cox_quart <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age + smoking_2cat_i + bmi"))  # set the formulas    
  
  model_summary <- coxph(formula_danish, data = bdd_danish_cases) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "main", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Heterogeneity tests ----
#### base 

heterogeneity_base <- map_dfr(
  proteomic_quart,
  function(var) {
    
    surv_obj <- Surv(time = bdd_danish_cases$follow_up,                       # set the outcomes
                     event = bdd_danish_cases$als)
    
    formula_raw <- as.formula("surv_obj ~ sex + baseline_age")
    test_1 <- coxph(formula_raw, data = bdd_danish_cases) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ sex + baseline_age"))  
    test_2 <- coxph(formula, data = bdd_danish_cases)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "main",
      p_value_heterogeneity = p_lr)
  })

#### adjusted 
heterogeneity_adjusted <- map_dfr(
  proteomic_quart,
  function(var) {
    
    surv_obj <- Surv(time = bdd_danish_cases$follow_up,                       # set the outcomes
                     event = bdd_danish_cases$als)
    
    formula_raw <- as.formula("surv_obj ~ sex + baseline_age + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = bdd_danish_cases) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ sex + baseline_age + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = bdd_danish_cases)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "main",
      p_value_heterogeneity = p_lr)
  })


heterogeneity_tests <- 
  bind_rows(heterogeneity_base, 
            heterogeneity_adjusted) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base 
trend_base <- data.frame(explanatory = character(),
                         model = factor(), 
                         p_value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age"))  
  model <- coxph(formula, data = bdd_danish_cases) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base <- rbind(trend_base, 
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "main", 
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted 
trend_adjusted <- data.frame(explanatory = character(),
                             model = factor(), 
                             p_value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ sex + baseline_age + smoking_2cat_i + bmi"))  
  model <- coxph(formula, data = bdd_danish_cases) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "main", 
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base, heterogeneity_adjusted,
   trend_base, trend_adjusted)


## Cox model (GAMs) ----
vars_labels <- 
  set_names(str_replace(proteomic, 
                        "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                        ""), 
            proteomic)


### Base ----
cox_gam_results_base <- 
  fit_and_plot_cox_gam(data = bdd_danish_cases, 
                       vars_labels = vars_labels, 
                       covariates = c("sex", "baseline_age"), 
                       title = "Base model")

### Adjusted ----
cox_gam_results_adjusted <- 
  fit_and_plot_cox_gam(data = bdd_danish_cases, 
                       vars_labels = vars_labels, 
                       covariates = c("sex", "baseline_age", "smoking_2cat_i", "bmi"), 
                       title = "Adjusted model")
rm(vars_labels, surv_obj)


## Assemblage ----
main_results <-        
  bind_rows(
    model1_cox_sd, model2_cox_sd,
    model1_cox_quart, model2_cox_quart) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se), 
    term = case_when(
      str_detect(term, "_sd") ~ "Continuous", 
      str_detect(term, "Q2") ~ "Quartile 2",
      str_detect(term, "Q3") ~ "Quartile 3",
      str_detect(term, "Q4") ~ "Quartile 4"), 
    explanatory = gsub("_sd", "", explanatory), 
    explanatory = gsub("_quart", "", explanatory), 
    HR_raw = HR, 
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    p_value_raw = p_value, 
    p_value_shape = ifelse(p_value_raw<0.05, "p_value<0.05", "p_value≥0.05"), 
    p_value = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
    p_value = ifelse(p_value == "1.00", ">0.99", p_value)) |>
  select(model, explanatory, analysis, term,  HR, HR_raw, `95% CI`, p_value, p_value_raw, p_value_shape, lower_CI, upper_CI)

main_results <- main_results |>
  left_join(heterogeneity_tests, by = c("explanatory", "analysis", "model")) |>
  mutate(p_value_heterogeneity = ifelse(term == "Continuous", NA, p_value_heterogeneity),
         p_value_heterogeneity = ifelse(p_value_heterogeneity < 0.01, "<0.01", number(p_value_heterogeneity, accuracy = 0.01, decimal.mark = ".")),
         p_value_heterogeneity = ifelse(p_value_heterogeneity == "1.00", ">0.99", p_value_heterogeneity))

main_results <- main_results |> 
  left_join(trend_tests, by = c("explanatory", "analysis", "model")) |>
  mutate(p_value_trend = ifelse(term == "Continuous", NA, p_value_trend),
         p_value_trend = ifelse(p_value_trend < 0.01, "<0.01", number(p_value_trend, accuracy = 0.01, decimal.mark = ".")),
         p_value_trend = ifelse(p_value_trend == "1.00", ">0.99", p_value_trend), 
         protein_group = case_when(str_detect(explanatory, 'proteomic_immun_res') ~ "Immune response", 
                                   str_detect(explanatory, 'proteomic_metabolism') ~ "Metabolism", 
                                   str_detect(explanatory, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
         explanatory = str_replace(explanatory, 'proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_', ""))

results_proteomic_ALS_occurrence_y_to_als <- list()
results_proteomic_ALS_occurrence_y_to_als$main_results_cox <- 
  list(main_results = main_results, 
       cox_gam_results_base = cox_gam_results_base, 
       cox_gam_results_adjusted = cox_gam_results_adjusted)

rm(model1_cox_sd, model2_cox_sd, 
   model1_cox_quart, model2_cox_quart, 
   heterogeneity_tests, trend_tests, 
   cox_gam_results_base, cox_gam_results_adjusted, 
   main_results)


## Tables and figures ----

### Table proteomic (sd) - als survival (main) ----
results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_sd_ALS_table_cox <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "main") |>            # select continuous results
  group_by(explanatory) |>                                                      # select explanatory vars significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, explanatory, protein_group, term, HR, "95% CI", "p_value") |>
  arrange(protein_group, explanatory) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS baseline associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "HR" = "Base models", "95% CI" = "Base models", "p-value" = "Base models", 
    "HR " = "Adjusted models", "95% CI " = "Adjusted models", "p-value " = "Adjusted models") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  merge_v(j = "protein_group") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic (quart) - als survival (main) ----
extra_rows <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_quart_ALS_table_cox <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_quart_ALS_table_cox <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_quart_ALS_table_cox |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p_value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p_value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "HR" = "Base models", "95% CI" = "Base models", "p-value" = "Base models", 
    "Heterogeneity test" = "Base models",  "Trend test" = "Base models",
    "HR " = "Adjusted models", "95% CI " = "Adjusted models", "p-value " = "Adjusted models", 
    "Heterogeneity test " = "Adjusted models",  "Trend test " = "Adjusted models") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

rm(extra_rows)


### Figure proteomic - als survival - base sd (main) ----
results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_sd_ALS_base_figure_cox <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "main") |>
  mutate(
    log2HR = log2(HR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & HR > 1 ~ "HR>1 & p-value<0.05",
      p_value_raw < 0.05 & HR < 1 ~ "HR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = log2HR, y = neg_log10_p, color = significance)) +
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
    values = c("HR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "HR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base cox models",
    x = "log2HR",
    y = "-log10(p-value)", 
    color = "")

### Figure proteomic - als survival - adjusted sd (main) ----
results_proteomic_ALS_occurrence_y_to_als$main_results_cox$proteomic_sd_ALS_adjusted_figure_cox <- 
  results_proteomic_ALS_occurrence_y_to_als$main_results_cox$main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "main") |>
  mutate(
    log2HR = log2(HR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & HR > 1 ~ "HR>1 & p-value<0.05",
      p_value_raw < 0.05 & HR < 1 ~ "HR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = log2HR, y = neg_log10_p, color = significance)) +
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
    values = c("HR<1 & p-value<0.05" = "blue", "p-value>0.05" = "grey70", "HR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted cox models",
    x = "log2HR",
    y = "-log10(p-value)", 
    color = "")



### Figure proteomic - als survival - base gam (main) ----
results_proteomic_ALS_occurrence_y_to_als$main_results_cox$plot_base_cox_gam_danish <- list(
  signif = map(keep(results_proteomic_ALS_occurrence_y_to_als$main_results_cox$cox_gam_results_base, ~ .x$pval_raw <= 0.05), "plot"),
  not_signif = map(keep(results_proteomic_ALS_occurrence_y_to_als$main_results_cox$cox_gam_results_base, ~ .x$pval_raw > 0.05), "plot"))


### Figure proteomic - als survival - adjusted gam (main) ----
results_proteomic_ALS_occurrence_y_to_als$main_results_cox$plot_adjusted_cox_gam_danish <- list(
  signif = map(keep(results_proteomic_ALS_occurrence_y_to_als$main_results_cox$cox_gam_results_adjusted, ~ .x$pval_raw <= 0.05), "plot"),
  not_signif = map(keep(results_proteomic_ALS_occurrence_y_to_als$main_results_cox$cox_gam_results_adjusted, ~ .x$pval_raw > 0.05), "plot"))

## 1.2 XGBoost -----

# 2. Analyse logistique à horizon fixe ----
bdd_danish <- bdd_danish |> 
  group_by(match) |>
  mutate(follow_up_no_na = follow_up) |>          # attention follow_up est en mois 
  fill(follow_up_no_na, .direction = "downup") |> 
  ungroup() |>
  mutate(
    y_10ans = ifelse(als == 1 & follow_up_no_na <= 120, 1, 0),      # attention follow_up_no_na est en mois 
    y_20ans = ifelse(als == 1 & follow_up_no_na <= 240, 1, 0))       # attention follow_up_no_na est en mois 

bdd_danish |> select(als, follow_up, follow_up_bis, follow_up_no_na, y_10ans, y_20ans) |> 
  tbl_summary(by = als)

## 2.1 Proteome wide associations (logistic regressions) ----
### Horizon 10 ans ----
#### model 1 sd ----
model1_sd_y_10ans <- data.frame(explanatory = character(),
                        term = integer(),
                        OR = numeric(),
                        lower_CI = numeric(),
                        upper_CI = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ birth_year + sex"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_y_10ans <- rbind(model1_sd_y_10ans, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model1_sd_y_10ans <- model1_sd_y_10ans |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main_y_10ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### model 1 quartiles ----
model1_quart_y_10ans <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ birth_year + sex"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_y_10ans <- rbind(model1_quart_y_10ans, data.frame(explanatory = var,
                                                 term = term, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 p_value = p_value))
}

model1_quart_y_10ans <- model1_quart_y_10ans |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main_y_10ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_y_10ans <- data.frame(explanatory = character(),
                                 model = factor(),
                                 p_value_heterogeneity = numeric(), 
                                 stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- glm(y_10ans ~ birth_year + sex, family = binomial, data = bdd_danish)
  
  formula <- as.formula(paste("y_10ans ~", var, "+ birth_year + sex"))
  test_2 <- glm(formula, family = binomial, data = bdd_danish)
  
  anova <- broom::tidy(anova(test_1, test_2, test = "LR"))
  p_value_heterogeneity <- anova$p.value[2]
  
  heterogeneity_base_y_10ans <- rbind(heterogeneity_base_y_10ans, 
                              data.frame(explanatory = var,
                                         model = "base", 
                                         analysis = "main_y_10ans",
                                         p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base_y_10ans <- data.frame(explanatory = character(),
                         model = factor(), 
                         p_value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ birth_year + sex"))
  test <- 
    glm(formula, family = binomial, data = bdd_danish) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_y_10ans <- rbind(trend_base_y_10ans, 
                      data.frame(explanatory = var,
                                 model = "base", 
                                 analysis = "main_y_10ans",
                                 p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)

#### model 2 sd ----
model2_sd_y_10ans <- data.frame(explanatory = character(),
                        term = integer(),
                        OR = numeric(),
                        lower_CI = numeric(),
                        upper_CI = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_sd_y_10ans <- rbind(model2_sd_y_10ans, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model2_sd_y_10ans <- model2_sd_y_10ans |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main_y_10ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

#### model 2 quartiles ----
model2_quart_y_10ans <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_y_10ans <- rbind(model2_quart_y_10ans, data.frame(explanatory = var,
                                                 term = term, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 p_value = p_value))
}

model2_quart_y_10ans <- model2_quart_y_10ans |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main_y_10ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_y_10ans <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_heterogeneity = numeric(), 
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- glm(y_10ans ~ sex + birth_year + smoking_2cat_i + bmi, family = binomial, data = bdd_danish)
  
  formula <- as.formula(paste("y_10ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  test_2 <- glm(formula, family = binomial, data = bdd_danish)
  
  anova <- broom::tidy(anova(test_1, test_2, test = "LR"))
  p_value_heterogeneity <- anova$p.value[2]
  
  heterogeneity_adjusted_y_10ans <- rbind(heterogeneity_adjusted_y_10ans, 
                                  data.frame(explanatory = var,
                                             model = "adjusted", 
                                             analysis = "main_y_10ans",
                                             p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted_y_10ans <- data.frame(explanatory = character(),
                             model = factor(), 
                             p_value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("y_10ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  test <- glm(formula, family = binomial, data = bdd_danish) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_y_10ans <- rbind(trend_adjusted_y_10ans, 
                          data.frame(explanatory = var,
                                     model = "adjusted", 
                                     analysis = "main_y_10ans",
                                     p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)


### Horizon 20 ans ----
#### model 1 sd ----
model1_sd_y_20ans <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ birth_year + sex"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_y_20ans <- rbind(model1_sd_y_20ans, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model1_sd_y_20ans <- model1_sd_y_20ans |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main_y_20ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### model 1 quartiles ----
model1_quart_y_20ans <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ birth_year + sex"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_y_20ans <- rbind(model1_quart_y_20ans, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_quart_y_20ans <- model1_quart_y_20ans |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main_y_20ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_y_20ans <- data.frame(explanatory = character(),
                                         model = factor(),
                                         p_value_heterogeneity = numeric(), 
                                         stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- glm(y_20ans ~ birth_year + sex, family = binomial, data = bdd_danish)
  
  formula <- as.formula(paste("y_20ans ~", var, "+ birth_year + sex"))
  test_2 <- glm(formula, family = binomial, data = bdd_danish)
  
  anova <- broom::tidy(anova(test_1, test_2, test = "LR"))
  p_value_heterogeneity <- anova$p.value[2]
  
  heterogeneity_base_y_20ans <- rbind(heterogeneity_base_y_20ans, 
                              data.frame(explanatory = var,
                                         model = "base", 
                                         analysis = "main_y_20ans",
                                         p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base_y_20ans <- data.frame(explanatory = character(),
                                 model = factor(), 
                                 p_value_trend = numeric(), 
                                 stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ birth_year + sex"))
  test <- 
    glm(formula, family = binomial, data = bdd_danish) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_y_20ans <- rbind(trend_base_y_20ans, 
                              data.frame(explanatory = var,
                                         model = "base", 
                                         analysis = "main_y_20ans",
                                         p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)

#### model 2 sd ----
model2_sd_y_20ans <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_sd_y_20ans <- rbind(model2_sd_y_20ans, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model2_sd_y_20ans <- model2_sd_y_20ans |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main_y_20ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

#### model 2 quartiles ----
model2_quart_y_20ans <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_y_20ans <- rbind(model2_quart_y_20ans, data.frame(explanatory = var,
                                                 term = term, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 p_value = p_value))
}

model2_quart_y_20ans <- model2_quart_y_20ans |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main_y_20ans") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_y_20ans <- data.frame(explanatory = character(),
                                             model = factor(), 
                                             p_value_heterogeneity = numeric(), 
                                             stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- glm(y_20ans ~ sex + birth_year + smoking_2cat_i + bmi, family = binomial, data = bdd_danish)
  
  formula <- as.formula(paste("y_20ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  test_2 <- glm(formula, family = binomial, data = bdd_danish)
  
  anova <- broom::tidy(anova(test_1, test_2, test = "LR"))
  p_value_heterogeneity <- anova$p.value[2]
  
  heterogeneity_adjusted_y_20ans <- rbind(heterogeneity_adjusted_y_20ans, 
                                          data.frame(explanatory = var,
                                                     model = "adjusted", 
                                                     analysis = "main_y_20ans",
                                                     p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted_y_20ans <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("y_20ans ~", var, "+ sex + birth_year + smoking_2cat_i + bmi"))
  test <- glm(formula, family = binomial, data = bdd_danish) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_y_20ans <- rbind(trend_adjusted_y_20ans, 
                                  data.frame(explanatory = var,
                                             model = "adjusted", 
                                             analysis = "main_y_20ans",
                                             p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)


### Assemblage ----
main_results_y_10_20 <- 
  bind_rows(model1_sd_y_10ans,                                            # main analyses
            model2_sd_y_10ans, 
            model1_quart_y_10ans, 
            model2_quart_y_10ans, 
            
            model1_sd_y_20ans,                                            # main analyses
            model2_sd_y_20ans, 
            model1_quart_y_20ans, 
            model2_quart_y_20ans) |>
  
  mutate(explanatory = if_else(str_detect(explanatory, "NEFL"), "proteomic_neuro_explo_NEFL", explanatory), 
         
         explanatory = gsub("_quart_med", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory), 
         explanatory = gsub("_sd", "", explanatory), 
         
         OR_raw = OR, 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p_value_raw = p_value, 
         p_value = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
         "95% CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  group_by(model) |>                               
  mutate(
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main_y_10ans" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      NA_real_), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main_y_20ans" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main_y_10ans" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main_y_20ans" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr)) |>
  ungroup() |>
  arrange(explanatory) |>
  select(analysis, 
         model,
         explanatory, 
         term,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p_value"), 
         lower_CI, upper_CI) 

heterogeneity_tests_y_10ans <-                                                  # y_10ans
  bind_rows(heterogeneity_base_y_10ans, heterogeneity_adjusted_y_10ans)  |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_y_10ans <-      
  bind_rows(trend_base_y_10ans, trend_adjusted_y_10ans) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_y_20ans <-                                                  # y_20ans
  bind_rows(heterogeneity_base_y_20ans, heterogeneity_adjusted_y_20ans) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_y_20ans <- 
  bind_rows(trend_base_y_20ans, trend_adjusted_y_20ans) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))


heterogeneity_tests_y_10_20 <- bind_rows(heterogeneity_tests_y_10ans, heterogeneity_tests_y_20ans)
trend_tests_y_10_20 <- bind_rows(trend_tests_y_10ans, trend_tests_y_20ans)


main_results_y_10_20 <- left_join(main_results_y_10_20, heterogeneity_tests_y_10_20, by = c("explanatory", "model", "analysis"))
main_results_y_10_20 <- left_join(main_results_y_10_20, trend_tests_y_10_20, by = c("explanatory", "model", "analysis"))

main_results_y_10_20 <- main_results_y_10_20 |>
  mutate(
    p_value_heterogeneity = ifelse(term == "Continuous", NA, p_value_heterogeneity), 
    p_value_trend = ifelse(term == "Continuous", NA, p_value_trend), 
    p_value_heterogeneity = ifelse(p_value_heterogeneity < 0.01, "<0.01", number(p_value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p_value_trend = ifelse(p_value_trend < 0.01, "<0.01", number(p_value_trend, accuracy = 0.01, decimal.mark = ".")), 
    protein_group = case_when(str_detect(explanatory, 'proteomic_immun_res') ~ "Immune response", 
                              str_detect(explanatory, 'proteomic_metabolism') ~ "Metabolism", 
                              str_detect(explanatory, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
    explanatory = str_replace(explanatory, 'proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_', ""))

results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide <-  
                list(main_results_y_10_20 = main_results_y_10_20)


rm(model1_sd_y_10ans, model1_quart_y_10ans, model2_sd_y_10ans, model2_quart_y_10ans, 
   model1_sd_y_20ans, model1_quart_y_20ans, model2_sd_y_20ans, model2_quart_y_20ans, 
   
   heterogeneity_base_y_10ans, heterogeneity_adjusted_y_10ans, 
   trend_base_y_10ans, trend_adjusted_y_10ans, 
   heterogeneity_base_y_20ans, heterogeneity_adjusted_y_20ans, 
   trend_base_y_20ans, trend_adjusted_y_20ans, 
   heterogeneity_tests_y_10ans, heterogeneity_tests_y_20ans, 
   trend_tests_y_10ans, trend_tests_y_20ans, 
   heterogeneity_tests_y_10_20, trend_tests_y_10_20, 
   
   main_results_y_10_20)


### Tables and figures ----
#### Table proteomic - als occurence - base and adjusted sd ----
results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_sd_ALS_table_y_10_20 <-                                                       # select both base and adjusted results
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis %in% c("main_y_10ans", "main_y_20ans")) |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted")) |>
  select(model, explanatory, analysis, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, model, contains("main_y_10ans"), contains("main_y_20ans")) |>
  rename("OR" = "OR_main_y_10ans", "95% CI" = "95% CI_main_y_10ans", "p-value" = "p_value_main_y_10ans", 
         "OR " = "OR_main_y_20ans", "95% CI " = "95% CI_main_y_20ans", "p-value " = "p_value_main_y_20ans") |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins within 10 years or 20 years.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "model" = "Models", 
    "OR" = "ALS diagnosis within 10 years", "95% CI" = "ALS diagnosis within 10 years", "p-value" = "ALS diagnosis within 10 years", 
    "OR " = "ALS diagnosis within 20 years", "95% CI " = "ALS diagnosis within 20 years", "p-value " = "ALS diagnosis within 20 years") |>
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
  merge_v(j = "model") |>
  bold(j = "model", part = "body") |>
  align(j = "model", align = "left", part = "all") |> 
  merge_at(j = "model", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |>
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")




#### Table proteomic - als occurence - base and adjusted quart ----
extra_rows <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis %in% c("main_y_10ans", "main_y_20ans")) |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory, model) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_main_y_10ans" = '-', "95% CI_main_y_10ans" = '-', "p_value_main_y_10ans" = '', "p_value_heterogeneity_main_y_10ans" = '', "p_value_trend_main_y_10ans" = '',
    "OR_main_y_20ans" = '-', "95% CI_main_y_20ans" = '-', "p_value_main_y_20ans" = '', "p_value_heterogeneity_main_y_20ans" = '', "p_value_trend_main_y_20ans" = '')

results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_quart_ALS_table_y_10_20 <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis %in% c("main_y_10ans", "main_y_20ans")) |>          
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, analysis, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, model, quartiles = term,  contains("main_y_10ans"), contains("main_y_20ans")) 


results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_quart_ALS_table_y_10_20 <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_quart_ALS_table_y_10_20 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_main_y_10ans = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_main_y_10ans[quartiles == 'Quartile 2'], ''), 
         p_value_trend_main_y_10ans = ifelse(quartiles == 'Quartile 1', p_value_trend_main_y_10ans[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_main_y_20ans = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_main_y_20ans[quartiles == 'Quartile 2'], ''), 
         p_value_trend_main_y_20ans = ifelse(quartiles == 'Quartile 1', p_value_trend_main_y_20ans[quartiles == 'Quartile 2'], ''), 
         model = fct_recode(model, "Base model" = "base", "Adjusted model" = "adjusted"), 
         model = fct_relevel(model, "Base model", "Adjusted model")) |>
  ungroup() |>
  rename("OR" = "OR_main_y_10ans", "95% CI" = "95% CI_main_y_10ans", "p-value" = "p_value_main_y_10ans", "Heterogeneity test" = "p_value_heterogeneity_main_y_10ans",  "Trend test" = "p_value_trend_main_y_10ans",
         "OR " = "OR_main_y_20ans", "95% CI " = "95% CI_main_y_20ans", "p-value " = "p_value_main_y_20ans", "Heterogeneity test " = "p_value_heterogeneity_main_y_20ans",  "Trend test " = "p_value_trend_main_y_20ans") |>
  arrange(protein_group, explanatory, model, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, adjusted on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, adjusted on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, adjusted on sex and birth year, smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, adjusted on sex and birth year, smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "model" = "Models", 
    "quartiles" = "Quartiles",
    "OR" = "ALS diagnosis within 10 years", "95% CI" = "ALS diagnosis within 10 years", "p-value" = "ALS diagnosis within 10 years",  "Heterogeneity test" = "ALS diagnosis within 10 years",   "Trend test" = "ALS diagnosis within 10 years", 
    "OR " = "ALS diagnosis within 20 years", "95% CI " = "ALS diagnosis within 20 years", "p-value " = "ALS diagnosis within 20 years", "Heterogeneity test " = "ALS diagnosis within 20 years",   "Trend test " = "ALS diagnosis within 20 years") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "model") |>
  bold(j = "model", part = "body") |>
  align(j = "model", align = "left", part = "all") |> 
  merge_at(j = "model", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

rm(extra_rows)


#### Figure proteomic - als occurence - base sd 10 y----
results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_sd_ALS_base_figure_y_10ans <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "main_y_10ans") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
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
    title = "Base logistic models - within 10 years",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

#### Figure proteomic - als occurence - base sd 20 y----
results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_sd_ALS_base_figure_y_20ans <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "main_y_20ans") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
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
    title = "Base logistic models - within 20 years",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))




#### Figure proteomic - als occurence - adjusted sd 10 y ----
results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_sd_ALS_adjusted_figure_y_10ans <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "main_y_10ans") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
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
    title = "Adjusted logistic models - within 10 years",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

#### Figure proteomic - als occurence - adjusted sd 20 y ----
results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$proteomic_sd_ALS_adjusted_figure_y_20ans <- 
  results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$proteome_wide$main_results_y_10_20 |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "main_y_20ans") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
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
    title = "Adjusted logistic models - within 20 years",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


## 2.2 Machine learning (XGboost) ----
### horizon 10 ans ----
params_binary <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",            # Notre métrique cible reste l'AUC
  eta              = 0.01,             # Apprentissage progressif pour éviter le surapprentissage
  max_depth        = 3,                # Permet de capter les interactions (ex: protéine x âge)
  subsample        = 0.8,              # 80% des lignes par arbre
  colsample_bytree = 0.6,              # 60% des variables testées par arbre (aide à l'ajustement)
  alpha            = 1,                # Régularisation L1 (Lasso)
  lambda           = 1)                # Régularisation L2 (Ridge)


X_matrix_10 <- bdd_danish |> 
  select(birth_year, sex, bmi, smoking_2cat_i, all_of(proteomic)) |> 
  mutate(
    sex_male = ifelse(sex == "Male", 1, 0),
    smoking_ever = ifelse(smoking_2cat_i == "Ever", 1, 0)) |> 
  select(-sex, -smoking_2cat_i) |> 
  as.matrix()

dtrain_10 <- xgb.DMatrix(data = X_matrix_10, label = bdd_danish$y_10ans)

# XGbosst avec validation croisée pour l'horizon 10 ans
set.seed(1996)
cv_10ans <- xgb.cv(
  params                = params_binary,
  data                  = dtrain_10,
  nrounds               = 1000,
  nfold                 = 5,
  early_stopping_rounds = 50,
  verbose              = FALSE)


# Récupération du nombre d'arbres optimal
best_nrounds_10 <- cv_10ans$best_iteration
if (is.null(best_nrounds_10) || length(best_nrounds_10) == 0 || best_nrounds_10 == 0) {
  best_nrounds_10 <- which.max(cv_10ans$evaluation_log$test_auc_mean)
}
cat("---> Nombre optimal d'arbres retenu (10 ans) :", best_nrounds_10, "\n")
cat("---> AUC maximum en cross-validation :", max(cv_10ans$evaluation_log$test_auc_mean), "\n")

# Entraînement final
final_xgb_10 <- xgb.train(
  params  = params_binary,
  data    = dtrain_10,
  nrounds = best_nrounds_10)

# Extraction et affichage de l'importance
importance_10 <- xgb.importance(feature_names = colnames(X_matrix_10), model = final_xgb_10)
cat("\nTop 10 des variables - Horizon 10 ans :\n")
print(head(importance_10, 10))

# Graphique
xgb.plot.importance(
  importance_matrix = importance_10[1:20, ], 
  main = "Top 20 de l'importance des variables (Horizon 10 ans)")

rm(params_binary, X_matrix_10, dtrain_10)

### horizon 20 ans ----
params_binary <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",            # Notre métrique cible reste l'AUC
  eta              = 0.01,             # Apprentissage progressif pour éviter le surapprentissage
  max_depth        = 3,                # Permet de capter les interactions (ex: protéine x âge)
  subsample        = 0.8,              # 80% des lignes par arbre
  colsample_bytree = 0.6,              # 60% des variables testées par arbre (aide à l'ajustement)
  alpha            = 1,                # Régularisation L1 (Lasso)
  lambda           = 1)                # Régularisation L2 (Ridge)

# Construction de la matrice X
X_matrix_20 <- bdd_danish |> 
  select(birth_year, sex, bmi, smoking_2cat_i, all_of(proteomic)) |> 
  mutate(
    sex_male = ifelse(sex == "Male", 1, 0),
    smoking_ever = ifelse(smoking_2cat_i == "Ever", 1, 0)) |> 
  select(-sex, -smoking_2cat_i) |> 
  as.matrix()

dtrain_20 <- xgb.DMatrix(data = X_matrix_20, label = bdd_danish$y_20ans)

# Validation croisée pour l'horizon 20 ans
set.seed(1996)
cv_20ans <- xgb.cv(
  params                = params_binary,
  data                  = dtrain_20,
  nrounds               = 1000,
  nfold                 = 5,
  early_stopping_rounds = 50,
  verbose               = FALSE)

# Récupération du nombre d'arbres optimal
best_nrounds_20 <- cv_20ans$best_iteration
if (is.null(best_nrounds_20) || length(best_nrounds_20) == 0 || best_nrounds_20 == 0) {
  best_nrounds_20 <- which.max(cv_20ans$evaluation_log$test_auc_mean)
}
cat("---> Nombre optimal d'arbres retenu (20 ans) :", best_nrounds_20, "\n")
cat("---> AUC maximum en cross-validation :", max(cv_20ans$evaluation_log$test_auc_mean), "\n")

# Entraînement final (20 ans)
final_xgb_20 <- xgb.train(
  params  = params_binary,
  data    = dtrain_20,
  nrounds = best_nrounds_20)

# Extraction et affichage de l'importance (20 ans)
importance_20 <- xgb.importance(feature_names = colnames(X_matrix_20), model = final_xgb_20)
cat("\nTop 10 des variables - Horizon 20 ans :\n")
print(head(importance_20, 10))

# Graphique
xgb.plot.importance(
  importance_matrix = importance_20[1:20, ], 
  main = "Top 20 de l'importance des variables (Horizon 20 ans)")

rm(params_binary, X_matrix_20, dtrain_20)

results_proteomic_ALS_occurrence_y_to_als$fixed_horizon_logistic_regression_analysis$machine_learning_xgboost <-  
  list(cv_10ans = cv_10ans, 
       final_xgb_10 = final_xgb_10, 
       importance_10 = importance_10, 
       best_nrounds_10 = best_nrounds_10,
       
       cv_20ans = cv_20ans, 
       final_xgb_20 = final_xgb_20, 
       importance_20 = importance_20, 
       best_nrounds_20 = best_nrounds_20)

rm(cv_10ans, final_xgb_10, importance_10, best_nrounds_10, 
  cv_20ans, final_xgb_20,  importance_20,  best_nrounds_20)


# 3. Stratified logistic models?----
# 4. Logistic models including time as a predictor?----
# 5. Analyse Cox tronqué (Landmark)?----

# 6. Cox + AUC dépendante du temps?----

# 7. Modèle AFT (survreg)? ----

# Aline Davias
# April 21, 2025 
# Analysis of survival after ALS diagnosis depending on fatty acids levels  

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_fattyacids_ALS_occurrence.R")

# Danish cohort ----
## Verif data ----
bdd_danish |> 
  filter(als == 1) |>
  select(time_diagnosis_death) |>
  tbl_summary()

densityplot(data = bdd_finnish, vars = "time_diagnosis_death")

bdd |> 
  filter(study %in% c("Danish")) |>
  filter(als == 1) |>
  select(sex, baseline_age, all_of(explanatory_raw), all_of(explanatory_quart)) |>
  tbl_summary( type = all_of(explanatory_raw) ~ "continuous",
               statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", 
                                all_categorical() ~ "{n} ({p}%)")) |>
  add_n()


## Base model ----
## Adjusted model ----

# Finnish cohort ----
## Verif data ----
bdd_finnish |> 
  filter(als == 1) |>
  select(follow_up_death, status_death) |>
  mutate(status_death = as.factor(as.character(status_death))) |>
  tbl_summary()

bdd_cases_finnish <- bdd_finnish |> 
  filter(als == 1)
densityplot(data = bdd_cases_finnish, vars = "follow_up_death")
descrip_num(data = bdd_cases_finnish, vars = "follow_up_death")

bdd |> 
  filter(study %in% c("FMC", "FMCF")) |>
  filter(als == 1) |>
  select(follow_up_death, status_death, sex, baseline_age, thawed, level_urbanization, all_of(covariates_finnish), education) |>
  tbl_summary(statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", 
                                all_categorical() ~ "{n} ({p}%)")) |>
  add_n()


bdd |> 
  filter(study == "FMCF") |>
  filter(als == 1) |>
  select(follow_up_death, status_death, sex, baseline_age, thawed, level_urbanization, all_of(explanatory_raw), all_of(explanatory_quart)) |>
  tbl_summary( type = all_of(explanatory_raw) ~ "continuous",
               statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", 
                                all_categorical() ~ "{n} ({p}%)")) |>
  add_n()

run_cox <- function(formula, data) {
  model <- coxph(formula, data = data)
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"])
}

## Base model sd ----
model1_cox_sd_finnish <- map_dfr(explanatory, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))
  
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_FMC <-                                                                # set the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox models 
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl) |>
    filter(explanatory == term)
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    model = "base", 
    explanatory = expl,
    term = expl,
    OR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Base model quart ----
model1_cox_quart_finnish <- map_dfr(explanatory_quart, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))
  
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_FMC <-                                                                # creation of the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox model
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl) |>
    filter(str_starts(term, explanatory)) |>
    mutate(term = case_when(
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"))
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        OR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

bdd |>
  filter(study == "FMCF") |>
  filter(als == 1) |>
  select(pufas_quart, baseline_age, sex, thawed) |>
  tbl_summary(by = "pufas_quart")

## Adjusted model sd ----
model2_cox_sd_finnish <- map_dfr(explanatory, function(expl) {
  
  bdd_cases_FMC <- bdd |> filter(als == 1) |> filter(study == "FMC")            # set the datasets
  bdd_cases_FMCF <- bdd |> filter(als == 1) |> filter(study == "FMCF")
  
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                              paste(covariates_finnish, collapse = " + "), 
                              "+ baseline_age + sex + thawed + level_urbanization"))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                  paste(covariates_finnish, collapse = " + "), 
                                  "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox models 
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl) |>
    filter(explanatory == term)
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    model = "adjusted", 
    explanatory = expl,
    term = expl,
    OR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Ajusted model quart ----
model2_cox_quart_finnish <- map_dfr(explanatory_quart, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of fatty acids quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of fatty acids quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                                  paste(covariates_finnish, collapse = " + "), 
                                  "+ baseline_age + sex + thawed + level_urbanization"))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                   paste(covariates_finnish, collapse = " + "), 
                                   "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox model
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl) |>
    filter(str_starts(term, explanatory)) |>
    mutate(term = case_when(
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"))
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        OR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

## Assemblage ----
main_results_fattyacids_ALS_finnish <- 
  bind_rows(
    model1_cox_sd_finnish, model2_cox_sd_finnish,
    model1_cox_quart_finnish, model2_cox_quart_finnish) |>
  mutate(
    OR = as.numeric(sprintf("%.1f", OR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    "95% CI" = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`, 
    `p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
    term = ifelse(str_detect(term, "_sd"), "continuous", term), 
    explanatory = gsub("_sd", "", explanatory), 
    explanatory = gsub("_quart", "", explanatory)) 

rm(model1_cox_sd_finnish, model2_cox_sd_finnish,
   model1_cox_quart_finnish, model2_cox_quart_finnish)

# Verif results sd ----
# run_cox <- function(formula, data) {
#   model <- coxph(formula, data = data)
#   model_summary <- summary(model)
#   coefs <- model_summary$coefficients
#   tibble(
#     term = rownames(coefs),
#     coef = coefs[, "coef"],
#     se = coefs[, "se(coef)"])
# }
# 
# bdd_cases_FMC <- bdd |>                                                       # set the datasets
#   filter(als == 1) |>                                                         # case selection
#   filter(study == "FMC") |>                                                   # cohort selection
#   mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
#                 scale,
#                 .names = "{.col}_sd"))  
# 
# bdd_cases_FMCF <- bdd |>                                                      # set the datasets
#   filter(als == 1) |>                                                         # case selection
#   filter(study == "FMCF") |>                                                  # cohort selection
#   mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
#                 scale,
#                 .names = "{.col}_sd"))
# 
# surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
#                      event = bdd_cases_FMC$status_death)
# surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
#                       event = bdd_cases_FMCF$status_death)
# 
# formula_FMC <- as.formula("surv_obj_FMC ~ pufas_sd + baseline_age + sex + thawed")  # formulas
# formula_FMCF <- as.formula("surv_obj_FMCF ~ pufas_sd + baseline_age + sex + thawed")
# 
# results <- list(                                                              # run of the simple cox models 
#   finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
#   finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
#   bind_rows(.id = "dataset") %>%
#   mutate(var = se^2, explanatory = "pufas_sd") |>
#   filter(explanatory == term)
# 
# rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
# 
# tibble(                                                                       # creation of the table of results
#   model = "base", 
#   explanatory = "pufas_sd",
#   term = "pufas_sd",
#   OR = exp(as.numeric(rma_fit$beta)),
#   lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
#   upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
#   `p-value` = as.numeric(rma_fit$pval))
# 
# ## Verif results quart ----
# run_cox <- function(formula, data) {
#   model <- coxph(formula, data = data)
#   model_summary <- summary(model)
#   coefs <- model_summary$coefficients
#   tibble(
#     term = rownames(coefs),
#     coef = coefs[, "coef"],
#     se = coefs[, "se(coef)"])
# }
# 
# bdd_cases_FMC <- bdd |>                                                       # set the datasets
#   filter(als == 1) |>                                                         # case selection
#   filter(study == "FMC") |>                                                   # cohort selection
#   mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of fatty acids quartiles (cohort specific)
#                                              labels = c("Q1", "Q2", "Q3", "Q4")),
#                 .names = "{.col}_quart"))
# 
# bdd_cases_FMCF <- bdd |>                                                      # set the datasets
#   filter(als == 1) |>                                                         # case selection
#   filter(study == "FMCF") |>                                                  # cohort selection
#   mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of fatty acids quartiles (cohort specific)
#                                              labels = c("Q1", "Q2", "Q3", "Q4")),
#                 .names = "{.col}_quart"))
# 
# surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
#                      event = bdd_cases_FMC$status_death)
# surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death,
#                       event = bdd_cases_FMCF$status_death)
# 
# formula_FMC <- as.formula("surv_obj_FMC ~ pufas_ω3_quart + baseline_age + sex + thawed")  # formulas
# formula_FMCF <- as.formula("surv_obj_FMCF ~ pufas_ω3_quart + baseline_age + sex + thawed")
# 
# results <- list(                                                              # run of the simple cox models
#   finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
#   finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
#   bind_rows(.id = "dataset") |>
#   mutate(var = se^2,
#          explanatory = "pufas_ω3_quart") |>
#   filter(str_starts(term, explanatory)) |>
#   mutate(term = case_when(
#            str_detect(term, "Q2") ~ "quartile 2",
#            str_detect(term, "Q3") ~ "quartile 3",
#            str_detect(term, "Q4") ~ "quartile 4"))
# 
# meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
#   group_by(explanatory, term) |>
#   group_modify(~ {
#     rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
#     tibble(                                                                   # results table creation
#       OR = exp(as.numeric(rma_fit$beta)),
#       lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
#       upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
#       `p-value` = as.numeric(rma_fit$pval))
#   }) |>
#   ungroup() |>
#   mutate(model = "base") |>
#   relocate(model, explanatory, term)
# 
# meta_results
# explanatory_quart


### table fatty acids (sd) - als survival ----
fattyacids_sd_als_table_finnish <- main_results_fattyacids_ALS_finnish |>
  select(model, explanatory, term, OR, "95% CI", "p-value") |>
  filter(term == "continuous") |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = fct_recode(explanatory, !!!explanatory_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age, sex, level of urbanization and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
  2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of PUFAs.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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

### table fatty acids (quart) - als survival ----
quartile1_rows <- main_results_fattyacids_ALS_finnish %>%
  distinct(model, explanatory) %>%
  mutate(
    term = "quartile 1",
    OR = "-",
    "95% CI" = "-",
    `p-value` = "")

fattyacids_quart_als_table_finnish <- main_results_fattyacids_ALS_finnish |>
  filter(!term == "continuous") |>
  select(model, explanatory, term, OR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

fattyacids_quart_als_table_finnish <- 
  bind_rows(quartile1_rows, fattyacids_quart_als_table_finnish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p-value")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = explanatory_labels), 
         explanatory = fct_recode(explanatory, !!!explanatory_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age, sex, level of urbanization and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
  2Estimated risk of ALS when pre-disease serum concentration of PUFAs compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Fatty acids", 
    term = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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


### figure fatty acids (sd) - als survival ----
fattyacids_sd_als_figure_finnish <- main_results_fattyacids_ALS_finnish |> 
  filter(term == "continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = fattyacids_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!fattyacids_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = OR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "PUFAs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

### figure fatty acids (quart) - als survival ----
fattyacids_quart_als_figure_finnish <- main_results_fattyacids_ALS_finnish |> 
  filter(!term == "continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = fattyacids_labels),
         explanatory = fct_recode(explanatory, !!!fattyacids_labels), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = OR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "PUFAs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

results_fattyacids_ALS_survival <- 
  list(finnish = list(main_results_fattyacids_ALS_finnish = main_results_fattyacids_ALS_finnish, 
                      # covar_als_table_finnish = covar_als_table_finnish, 
                      fattyacids_sd_als_table_finnish = fattyacids_sd_als_table_finnish, 
                      fattyacids_quart_als_table_finnish = fattyacids_quart_als_table_finnish, 
                      fattyacids_sd_als_figure_finnish = fattyacids_sd_als_figure_finnish, 
                      fattyacids_quart_als_figure_finnish = fattyacids_quart_als_figure_finnish))

rm(main_results_fattyacids_ALS_finnish, 
   fattyacids_sd_als_table_finnish, 
   fattyacids_quart_als_table_finnish, 
   fattyacids_sd_als_figure_finnish, 
   fattyacids_quart_als_figure_finnish)



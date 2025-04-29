# Aline Davias
# April 9, 2025 
# Analyses on fatty acids and ALS occurrence in both Danish and Finnish cohorts


# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")

# Occurrence analyses ----
## Danish cohort ----
### Base model ----
model1_sd_danish <- map_dfr(explanatory, function(expl) {                       # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  tibble(
    model = "base",
    explanatory = expl,
    term = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

model1_quart_danish <- map_dfr(explanatory_quart, function(expl) {              # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  tibble(
    model = "base",
    explanatory = expl,
    term = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

### Adjusted model ----
model2_sd_danish <- map_dfr(explanatory, function(expl) {                       # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", expl, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  tibble(
    model = "adjusted",
    explanatory = expl,
    term = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

model2_quart_danish <- map_dfr(explanatory_quart, function(expl) {              # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", expl, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) 
  tibble(
    model = "adjusted",
    explanatory = expl,
    term = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

### Assemblage ----
main_results_fattyacids_ALS_danish <- 
  bind_rows(model1_sd_danish, model2_sd_danish, model1_quart_danish, model2_quart_danish) |>
  filter(str_starts(term, explanatory))   |>  # to keep only the quartiles term without the covariates (str_starts() is checking whether the string in term starts with the string in expl)
  mutate(
    OR = as.numeric(sprintf("%.1f", OR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    "95% CI" = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw`= `p-value`, 
    `p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
    term = ifelse(str_detect(term, "_sd"), "continuous", term), 
    term = ifelse(str_detect(term, "_quartQ2"), "quartile 2", term), 
    term = ifelse(str_detect(term, "_quartQ3"), "quartile 3", term), 
    term = ifelse(str_detect(term, "_quartQ4"), "quartile 4", term), 
    explanatory = gsub("_sd", "", explanatory), 
    explanatory = gsub("_quart", "", explanatory)) 

rm(model1_sd_danish, model2_sd_danish, model1_quart_danish, model2_quart_danish)

## Finnish cohorts ----
run_clogit <- function(formula, data) {
  model <- clogit(formula, data = data)
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"])
}

### Base model ----
model1_sd_finnish <- map_dfr(explanatory, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))                # base formula: matched, not ajstuded
  
  bdd_finnish_FMC <- bdd |>                                                     # set the datasets
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_finnish_FMCF <- bdd |>                                                    # set the datasets
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_clogit(formula, bdd_finnish_FMC),
    finnish_2 = run_clogit(formula, bdd_finnish_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl)
  
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

model1_quart_finnish <- map_dfr(explanatory_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))                # base formula: matched, not ajstuded
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of quartiles cohort specific                      
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of quartiles cohort specific    
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_clogit(formula, bdd_finnish_FMC),
    finnish_2 = run_clogit(formula, bdd_finnish_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "quartile 2",
             str_detect(term, "Q3") ~ "quartile 3",
             str_detect(term, "Q4") ~ "quartile 4",
             TRUE ~ NA_character_)) 
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
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

### Adjusted model ----
model2_sd_finnish <- map_dfr(explanatory, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+",                               # adjusted formula: matched and ajstuded
                              paste(covariates_finnish, collapse = " + "), 
                              "+ strata(match)"))
  
  bdd_finnish_FMC <- bdd |>                                                     # set the datasets
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_finnish_FMCF <- bdd |>                                                    # set the datasets
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(fattyacids),                                           # create cohort specific scaled fatty acids variables 
                  scale,
                  .names = "{.col}_sd"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_clogit(formula, bdd_finnish_FMC),
    finnish_2 = run_clogit(formula, bdd_finnish_FMCF)) |>
    bind_rows(.id = "dataset") |>
    mutate(var = se^2, explanatory = expl) |>
    filter(str_starts(term, explanatory))    # to keep only the quartiles term without the covariates (str_starts() is checking whether the string in term starts with the string in expl)
  
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

model2_quart_finnish <- map_dfr(explanatory_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+",                               # adjusted formula: matched and adjusded
                              paste(covariates_finnish, collapse = " + "), 
                              "+ strata(match)"))
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of quartiles cohort specific                      
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                    # creation of quartiles cohort specific    
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_clogit(formula, bdd_finnish_FMC),
    finnish_2 = run_clogit(formula, bdd_finnish_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "quartile 2",
             str_detect(term, "Q3") ~ "quartile 3",
             str_detect(term, "Q4") ~ "quartile 4",
             TRUE ~ NA_character_)) |>
    filter(str_detect(term, "quartile"))                                        # filter to remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
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

### Assemblage ----
main_results_fattyacids_ALS_finnish <- 
  bind_rows(model1_sd_finnish, model2_sd_finnish, model1_quart_finnish, model2_quart_finnish) |>
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

rm(model1_sd_finnish, model2_sd_finnish, model1_quart_finnish, model2_quart_finnish)

## Sensitivity analysis - remove samples with Ca2+<2 ----
descrip_num(data = bdd_finnish, vars = "S_Ca")
bdd |>
  filter(study %in% c('FMC', 'FMCF', 'MFH')) |>
  filter(S_Ca>2) 

# bdd_finnish_sensi_ca <- bdd |>           # filtrer pout garder que les triplets de cas-control qui ont tous Ca2+>2
#   filter(study %in% c('FMC', 'FMCF')) |>
#   group_by(match) |>
#   filter(n() == 3, all(S_Ca > 2)) |>    
#   ungroup()                              # n = 111 au lieu de 257

bdd_finnish_sensi_ca <- bdd |>           # filtrer pout garder que les lignes qui ont Ca2+>2
  filter(study %in% c('FMC', 'FMCF')) |>
  filter(S_Ca>2)                         # n = 177 au lieu de 257


run_logit_sensi_Ca <- function(formula, data) {
  model <- glm(formula, data = data, family = binomial())
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "Estimate"],
    se = coefs[, "Std. Error"]
  )
}

#### Base model ----
model1_sd_sensi_ca_finnish <- map_dfr(explanatory, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ sex + baseline_age + municipality + thawed"))         # base formula: matched, not ajstuded
  
  bdd_finnish_1 <- bdd |> filter(study == "FMC") |> filter(S_Ca>2)              # creation of one dataset per finnish cohort
  bdd_finnish_2 <- bdd |> filter(study == "FMCF") |> filter(S_Ca>2) 
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_logit_sensi_Ca(formula, bdd_finnish_1),
    finnish_2 = run_logit_sensi_Ca(formula, bdd_finnish_2)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl)
  
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

model1_quart_sensi_ca_finnish <- map_dfr(explanatory_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ + sex + baseline_age + municipality + thawed"))     # base formula: matched, not ajstuded
  
  bdd_finnish_1 <- bdd |> filter(study == "FMC")                                # creation of one dataset per finnish cohort
  bdd_finnish_2 <- bdd |> filter(study == "FMCF")
  
  bdd_finnish_1 <- bdd_finnish_1 |>                                             # creation of quartiles cohort specific 
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_2 <- bdd_finnish_2 |>                                             
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_logit_sensi_Ca(formula, bdd_finnish_1),
    finnish_2 = run_logit_sensi_Ca(formula, bdd_finnish_2)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "quartile 2",
             str_detect(term, "Q3") ~ "quartile 3",
             str_detect(term, "Q4") ~ "quartile 4",
             TRUE ~ NA_character_)) 
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
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

#### Adjusted model ----
model2_sd_sensi_ca_finnish <- map_dfr(explanatory, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+",                               # adjusted formula: matched and ajstuded
                              paste(covariates_finnish, collapse = " + "), 
                              "+ sex + baseline_age + municipality + thawed"))
  
  bdd_finnish_1 <- bdd |> filter(study == "FMC")                                # creation of one dataset per finnish cohort
  bdd_finnish_2 <- bdd |> filter(study == "FMCF")
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_logit_sensi_Ca(formula, bdd_finnish_1),
    finnish_2 = run_logit_sensi_Ca(formula, bdd_finnish_2)) |>
    bind_rows(.id = "dataset") |>
    mutate(var = se^2, explanatory = expl) |>
    filter(str_starts(term, explanatory))    # to keep only the quartiles term without the covariates (str_starts() is checking whether the string in term starts with the string in expl)
  
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

model2_quart_sensi_ca_finnish <- map_dfr(explanatory_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+",                               # adjusted formula: matched and adjusded
                              paste(covariates_finnish, collapse = " + "), 
                              "+ sex + baseline_age + municipality + thawed"))
  
  bdd_finnish_1 <- bdd |> filter(study == "FMC")                          # creation of one dataset per finnish cohort
  bdd_finnish_2 <- bdd |> filter(study == "FMCF")
  
  bdd_finnish_1 <- bdd_finnish_1 |>                                             # creation of quartiles cohort specific 
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_2 <- bdd_finnish_2 |>                                             
    mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  
  results <- list(                                                              # run of the simple conditional logistic regression
    finnish_1 = run_logit_sensi_Ca(formula, bdd_finnish_1),
    finnish_2 = run_logit_sensi_Ca(formula, bdd_finnish_2)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "quartile 2",
             str_detect(term, "Q3") ~ "quartile 3",
             str_detect(term, "Q4") ~ "quartile 4",
             TRUE ~ NA_character_)) |>
    filter(str_detect(term, "quartile"))                                        # filter to remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
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

#### Assemblage ----
sensi_ca_results_fattyacids_ALS_finnish <- 
  bind_rows(model1_sd_sensi_ca_finnish, model2_sd_sensi_ca_finnish, 
            model1_quart_sensi_ca_finnish, model2_quart_sensi_ca_finnish) |>
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

rm(model1_sd_sensi_ca_finnish, model2_sd_sensi_ca_finnish, 
   model1_quart_sensi_ca_finnish, model2_quart_sensi_ca_finnish, 
   bdd_finnish_sensi_ca)




# Results presentation ----
## Danish cohort ----
### table covar description by als ----
descriptive_table_by_als_danish <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1)) |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  


### table covar - als  ----
covar_als_table_danish <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
      select(als, all_of(covariates_danish)) |>
      tbl_uvregression(
        y = als,
        method = glm,
        method.args = list(family = binomial), 
        exponentiate = TRUE, 
        estimate_fun = label_number(accuracy = 0.1, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun)|>
      bold_labels(), 
    tbl_2 = clogit(als ~ sex + baseline_age +
                     smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
                   data = bdd_danish) |>
      tbl_regression(exponentiate = TRUE, 
                     estimate_fun = label_number(accuracy = .1, decimal.mark = "."),
                     pvalue_fun = custom_pvalue_fun) |>
      bold_labels()), 
  tab_spanner = c("**Univariate**", "**Adjusted**")) |> 
  as_flex_table()|>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")

### table fatty acids (sd) - als ----
fattyacids_sd_als_table_danish <- main_results_fattyacids_ALS_danish |>
  filter(term == "continuous") |>
  select(model, explanatory, OR, "95% CI", "p-value") |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = fct_recode(explanatory, !!!explanatory_labels)) |>
  flextable() |>
  add_footer_lines(
  "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. 
  2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of PUFAs.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Fatty acids", 
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

### table fatty acids (quart) - als ----
quartile1_rows <- main_results_fattyacids_ALS_danish %>%
  distinct(model, explanatory) %>%
  mutate(
    term = "quartile 1",
    OR = "-",
    "95% CI" = "-",
    `p-value` = "")

fattyacids_quart_als_table_danish <- main_results_fattyacids_ALS_danish |>
  filter(!term == "continuous") |>
  select(model, explanatory, term, OR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

fattyacids_quart_als_table_danish <- 
  bind_rows(quartile1_rows, fattyacids_quart_als_table_danish) |>
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
    "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. 
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

### figure fattyacids description by als ----
fattyacids_by_als_figure_danish <- bdd_danish |>
  select(als, all_of(explanatory_raw)) |>
  pivot_longer(cols = -als, names_to = "PUFAs", values_to = "Values") |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1"), 
    PUFAs = factor(PUFAs, levels = fattyacids_labels),
    PUFAs = fct_rev(PUFAs),
    PUFAs = fct_recode(PUFAs, !!!fattyacids_labels)) |> 
  arrange(PUFAs) |>
  ggplot() +
  aes(x = Values, y = PUFAs, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_x_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(fill = "ALS ", x = "Values (%, log transformed)", y = "PUFAs") +
  theme_lucid() + 
  theme(legend.position = "bottom")

### figure fatty acids (sd) - als ----
fattyacids_sd_als_figure_danish <- main_results_fattyacids_ALS_danish |> 
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

### figure fatty acids (quart) - als ----
fattyacids_quart_als_figure_danish <- main_results_fattyacids_ALS_danish |> 
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

## Finnish cohorts ----
### table covar description by als ----
descriptive_table_by_als_finnish <- bdd_finnish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1)) |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

### table covar - als  ----

### table fatty acids (sd) - als ----
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
    "1All models are matched for age, sex, municipality, and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
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

### table fatty acids (quart) - als ----
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
    "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
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

### figure fattyacids description by als ----
fattyacids_by_als_figure_finnish <- bdd_finnish|>
  select(als, all_of(explanatory_raw)) |>
  pivot_longer(cols = -als, names_to = "PUFAs", values_to = "Values") |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1"), 
    PUFAs = factor(PUFAs, levels = fattyacids_labels),
    PUFAs = fct_rev(PUFAs),
    PUFAs = fct_recode(PUFAs, !!!fattyacids_labels)) |> 
  arrange(PUFAs) |>
  ggplot() +
  aes(x = Values, y = PUFAs, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_x_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(fill = "ALS ", x = "Values (%, log transformed)", y = "PUFAs") +
  theme_lucid() + 
  theme(legend.position = "bottom")


### figure fatty acids (sd) - als ----
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

### figure fatty acids (quart) - als ----
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

### table sensi ca fatty acids (sd) - als ----
fattyacids_sd_als_table_sensi_ca_finnish <- sensi_ca_results_fattyacids_ALS_finnish |>
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
    "1All models are matched for age, sex, municipality, and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
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

### table sensi ca fatty acids (quart) - als ----
quartile1_rows <- sensi_ca_results_fattyacids_ALS_finnish %>%
  distinct(model, explanatory) %>%
  mutate(
    term = "quartile 1",
    OR = "-",
    "95% CI" = "-",
    `p-value` = "")

fattyacids_quart_als_table_sensi_ca_finnish <- sensi_ca_results_fattyacids_ALS_finnish |>
  filter(term == c("quartile 2", "quartile 3", "quartile 4")) |>
  select(model, explanatory, term, OR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

fattyacids_quart_als_table_sensi_ca_finnish <- 
  bind_rows(quartile1_rows, fattyacids_quart_als_table_sensi_ca_finnish) |>
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
    "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
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

### figure sensi ca fatty acids (sd) - als ----
fattyacids_sd_als_figure_sensi_ca_finnish <- sensi_ca_results_fattyacids_ALS_finnish |> 
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

### figure sensi ca fatty acids (quart) - als ----
fattyacids_quart_als_figure_sensi_ca_finnish <- sensi_ca_results_fattyacids_ALS_finnish |> 
  filter(term == c("quartile 2", "quartile 3", "quartile 4")) |>
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



results_fattyacids_ALS_occurrence <- 
  list(danish = list(main_results_fattyacids_ALS_danish = main_results_fattyacids_ALS_danish, 
                     descriptive_table_by_als_danish = descriptive_table_by_als_danish, 
                     covar_als_table_danish = covar_als_table_danish, 
                     fattyacids_sd_als_table_danish = fattyacids_sd_als_table_danish, 
                     fattyacids_quart_als_table_danish = fattyacids_quart_als_table_danish, 
                     fattyacids_by_als_figure_danish = fattyacids_by_als_figure_danish, 
                     fattyacids_sd_als_figure_danish = fattyacids_sd_als_figure_danish, 
                     fattyacids_quart_als_figure_danish = fattyacids_quart_als_figure_danish), 
       finnish = list(main_results_fattyacids_ALS_finnish = main_results_fattyacids_ALS_finnish, 
                      descriptive_table_by_als_finnish = descriptive_table_by_als_finnish, 
                      # covar_als_table_finnish = covar_als_table_finnish, 
                      fattyacids_sd_als_table_finnish = fattyacids_sd_als_table_finnish, 
                      fattyacids_quart_als_table_finnish = fattyacids_quart_als_table_finnish, 
                      fattyacids_by_als_figure_finnish = fattyacids_by_als_figure_finnish, 
                      fattyacids_sd_als_figure_finnish = fattyacids_sd_als_figure_finnish, 
                      fattyacids_quart_als_figure_finnish = fattyacids_quart_als_figure_finnish), 
       sensi_ca = list(sensi_ca_results_fattyacids_ALS_finnish = sensi_ca_results_fattyacids_ALS_finnish, 
                       fattyacids_sd_als_table_sensi_ca_finnish = fattyacids_sd_als_table_sensi_ca_finnish, 
                       fattyacids_quart_als_table_sensi_ca_finnish = fattyacids_quart_als_table_sensi_ca_finnish, 
                       fattyacids_sd_als_figure_sensi_ca_finnish = fattyacids_sd_als_figure_sensi_ca_finnish, 
                       fattyacids_quart_als_figure_sensi_ca_finnish = fattyacids_quart_als_figure_sensi_ca_finnish))

rm(main_results_fattyacids_ALS_danish, 
   descriptive_table_by_als_danish, 
   covar_als_table_danish, 
   fattyacids_sd_als_table_danish, 
   fattyacids_quart_als_table_danish, 
   fattyacids_by_als_figure_danish, 
   fattyacids_sd_als_figure_danish, 
   fattyacids_quart_als_figure_danish, 
   main_results_fattyacids_ALS_finnish, 
   descriptive_table_by_als_finnish, 
   # covar_als_table_finnish, 
   fattyacids_sd_als_table_finnish, 
   fattyacids_quart_als_table_finnish, 
   fattyacids_by_als_figure_finnish, 
   fattyacids_sd_als_figure_finnish, 
   fattyacids_quart_als_figure_finnish, 
   
   sensi_ca_results_fattyacids_ALS_finnish, 
   fattyacids_sd_als_table_sensi_ca_finnish, 
   fattyacids_quart_als_table_sensi_ca_finnish, 
   fattyacids_sd_als_figure_sensi_ca_finnish, 
   fattyacids_quart_als_figure_sensi_ca_finnish)

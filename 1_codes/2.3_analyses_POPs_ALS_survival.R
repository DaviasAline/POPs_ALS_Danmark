# Aline Davias
# April 29, 2025 
# Analysis of survival after ALS diagnosis depending on POPs levels  

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")

# Danish cohort ----
## Covar model ----
bdd_cases_danish <- bdd_danish |>                                               # set the dataset
  filter(als == 1) |>                                                           # case selection
  select(follow_up_death, status_death, sex, baseline_age, 
         bmi, marital_status_2cat_i, smoking_2cat_i, education_i, cholesterol_i) 

surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                        event = bdd_cases_danish$status_death)

covar_danish <- tbl_merge(tbls = list(
  tbl_uvregression(
    data = bdd_cases_danish, 
    y = surv_obj_danish, 
    method = survival::coxph,  
    exponentiate = TRUE,
    include = c("sex", "baseline_age", 
                "bmi", "marital_status_2cat_i", "smoking_2cat_i", "education_i", "cholesterol_i")) |>
    bold_labels() |>
    bold_p(), 
  tbl_regression(
    coxph(surv_obj_danish ~ sex + baseline_age + bmi + marital_status_2cat_i + smoking_2cat_i + education_i + cholesterol_i, data = bdd_cases_danish),
    exponentiate = TRUE) |>
    bold_labels() |>
    bold_p()), 
  tab_spanner = c("**Crude**", "**Adjusted**"))
rm(bdd_cases_danish, surv_obj_danish)

## Main analysis 
## Base model sd ----
model1_cox_sd_danish <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- bdd_danish |>                                             # set the dataset
    filter(als == 1) |>                                                         # case selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <-                                                             # set the formula
    as.formula(paste("surv_obj_danish ~", expl, "+ baseline_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "base", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
})

## Base model quart ----
model1_cox_quart_danish <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_danish <- bdd_danish |>                                             # set the datasets
    filter(als == 1) |>                                                         # case selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ baseline_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "base", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

## Adjusted model sd ----
model2_cox_sd_danish <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- bdd_danish |>                                             # set the dataset
    filter(als == 1) |>                            
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <- as.formula(paste("surv_obj_danish ~", expl, "+",            # set the formulas              
                                     paste(covariates_danish, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "adjusted", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

## Ajusted model quart ----
model2_cox_quart_danish <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_danish <- bdd_danish |>                                             # set the dataset
    filter(als == 1) |>                                                         # case selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <- as.formula(paste("surv_obj_danish ~", expl, "+",            # set the formulas              
                                     paste(covariates_danish, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "adjusted", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})


# Finnish cohort ----
run_cox <- function(formula, data) {
  model <- coxph(formula, data = data)
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"])
}

## Covar model ----
### crude ----
covar_crude_finnish <- map_dfr(c("baseline_age", "sex", "thawed", "level_urbanization", 
                           covariates_finnish), function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC")                                                      # cohort selection
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF")                                                     # cohort selection

  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl))                      # set the formulas
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl))
  
  results <- list(                                                              # run of the simple cox models 
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") |>
    mutate(var = se^2, variable = expl)
  
  meta_results <- results |>                                                    # run metanalyse
    group_by(variable, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        HR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "base") |> 
    relocate(model, variable, term)
  return(meta_results)
  })


### adjusted ----
bdd_cases_FMC <- bdd |>                                                         # set the datasets
  filter(als == 1) |>                                                           # case selection
  filter(study == "FMC")                                                        # cohort selection
 
bdd_cases_FMCF <- bdd |>                                                        # set the datasets
  filter(als == 1) |>                                                           # case selection
  filter(study == "FMCF")                                                       # cohort selection

formula <- "Surv(follow_up_death, status_death) ~ baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol"
formula <- as.formula(formula)

covar_adjusted_finnish <- list(                                                 # run of the simple cox models
  finnish_FMC = run_cox(formula, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2)

covar_adjusted_finnish <- covar_adjusted_finnish |>                             # run metanalyse
  group_by(term) |>
  group_modify( ~ {
    rma_fit <- rma(yi = .x$coef,
                   vi = .x$var,
                   method = "DL")
    tibble(                                                                     # results table creation
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval))
  }) |>
  ungroup() |>
  mutate(model = "adjusted", 
         variable = fct_recode(
           term, 
             "level_urbanization" = "level_urbanization2",
             "level_urbanization" = "level_urbanization3",
             "level_urbanization" = "level_urbanization4",
             "marital_status_2cat" = "marital_status_2catMarried/cohabit",
             "sex" = "sexMale",
             "smoking_2cat" = "smoking_2catEver",
             "thawed" = "thawed1")) |>
  relocate(model, variable, term) 

rm(bdd_cases_FMC, bdd_cases_FMCF, formula)
covar_finnish <- bind_rows(covar_crude_finnish, covar_adjusted_finnish)
rm(covar_crude_finnish, covar_adjusted_finnish)

## Base model sd ----
model1_cox_sd_finnish <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
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
    filter(explanatory == term)                                                 # remove the covariates results                           
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    study = "Finnish", 
    model = "base", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Base model quart ----
model1_cox_quart_finnish <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
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
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        HR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(study = "Finnish", model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

## Adjusted model sd ----
model2_cox_sd_finnish <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # filter to get only the cases
    filter(study == "FMC") |>                                                   # filter to get only one cohort  
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # filter to get only the cases
    filter(study == "FMCF") |>                                                  # filter to get only one cohort 
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
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
    filter(explanatory == term)                                                 # remove the covariates results
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    study = "Finnish", 
    model = "adjusted", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Adjusted model quart ----
model2_cox_quart_finnish <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
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
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        study = "Finnish", 
        HR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})


# Metaanalysis ----
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
model1_cox_sd_metanalysis <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- bdd |>                                                    # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "Danish") |>                                                # cohort selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ baseline_age + sex"))  
  formula_FMC <-                                                                
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox models 
    danish = run_cox(formula_danish, bdd_cases_danish),
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl) |>
    filter(explanatory == term)                                                 # remove the covariates results
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    study = "Metanalysis", 
    model = "base", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Base model quart ----
model1_cox_quart_metanalysis <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_danish <- bdd |>                                                    # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "Danish") |>                                                # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                    # set the outcomes
                          event = bdd_cases_danish$status_death)
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ baseline_age + sex")) 
  formula_FMC <-                                                                # creation of the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox model
    danish = run_cox(formula_danish, bdd_cases_danish),
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        HR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(study = "Metanalysis", 
           model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

## Adjusted model sd ----
model2_cox_sd_metanalysis <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- bdd |>                                                    # set the datasets
    filter(als == 1) |>                                                         # filter to get only the cases
    filter(study == "Danish") |>                                                # filter to get only one cohort
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables
                  scale,
                  .names = "{.col}_sd"))
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # filter to get only the cases
    filter(study == "FMC") |>                                                   # filter to get only one cohort  
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # filter to get only the cases
    filter(study == "FMCF") |>                                                  # filter to get only one cohort 
    mutate(across(all_of(POPs_group),                                           # create cohort specific scaled POPs variables 
                  scale,
                  .names = "{.col}_sd"))  
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_danish <-
    as.formula(paste("surv_obj_danish ~", expl,                                 # set the formulas
                     "+ baseline_age + sex +
                     smoking_2cat + bmi + cholesterol + marital_status_2cat + education"))
  formula_FMC <- 
    as.formula(paste("surv_obj_FMC ~", expl, 
                     "+ baseline_age + sex + 
                     thawed + level_urbanization + 
                     smoking_2cat + bmi + cholesterol + marital_status_2cat"))
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, 
                     "+ baseline_age + sex + 
                     thawed + level_urbanization + 
                     smoking_2cat + bmi + cholesterol + marital_status_2cat + education"))
  
  results <- list(                                                              # run of the simple cox models 
    danish = run_cox(formula_danish, bdd_cases_danish),
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, explanatory = expl) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
  rma_fit <- rma(yi = coef, vi = var, data = results, method = "DL")            # run of the meta-analyse (metafor package as Ian did)
  
  tibble(                                                                       # creation of the table of results
    study = "Metanalysis", 
    model = "adjusted", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval)
    
  )
})

## Ajusted model quart ----
model2_cox_quart_metanalysis <- map_dfr(POPs_group_quart, function(expl) {
  
  bdd_cases_danish <- bdd |>                                                    # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "Danish") |>                                                # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMC <- bdd |>                                                       # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMC") |>                                                   # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_cases_FMCF <- bdd |>                                                      # set the datasets
    filter(als == 1) |>                                                         # case selection
    filter(study == "FMCF") |>                                                  # cohort selection
    mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                    # creation of POPs quartiles (cohort specific)                        
                                               labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,                    # set the outcomes
                       event = bdd_cases_FMC$status_death)
  surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death, 
                        event = bdd_cases_FMCF$status_death)
  
  formula_danish <- 
    as.formula(paste("surv_obj_danish ~", expl,                                 # set the formulas  
                     "+ baseline_age + sex + 
                     smoking_2cat + bmi + cholesterol + marital_status_2cat + education"))
  formula_FMC <- 
    as.formula(paste("surv_obj_FMC ~", expl, 
                     "+ baseline_age + sex + 
                     thawed + level_urbanization + 
                     smoking_2cat + bmi + cholesterol + marital_status_2cat"))  # education not available
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, 
                     "+ baseline_age + sex + 
                     thawed + level_urbanization + 
                     smoking_2cat + bmi + cholesterol + marital_status_2cat + education"))
  
  results <- list(                                                              # run of the simple cox model
    danish = run_cox(formula_danish, bdd_cases_danish),
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") %>%
    mutate(var = se^2, 
           explanatory = expl) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per explanatory variable)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        study = "Metanalysis", 
        HR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

# Assemblage ----
main_results_POPs_ALS_survival <-       
  bind_rows(
    model1_cox_sd_danish, model2_cox_sd_danish,
    model1_cox_quart_danish, model2_cox_quart_danish) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  
  bind_rows(model1_cox_sd_finnish, model2_cox_sd_finnish,
            model1_cox_quart_finnish, model2_cox_quart_finnish, 
            model1_cox_sd_metanalysis, model2_cox_sd_metanalysis,
            model1_cox_quart_metanalysis, model2_cox_quart_metanalysis) |>
  mutate(
    term = case_when(
      str_detect(term, "_sd") ~ "Continuous", 
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"), 
    explanatory = gsub("_sd", "", explanatory), 
    explanatory = gsub("_quart", "", explanatory), 
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`, 
    `p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI)

rm(model1_cox_sd_finnish, model2_cox_sd_finnish,
   model1_cox_quart_finnish, model2_cox_quart_finnish, 
   model1_cox_sd_danish, model2_cox_sd_danish,
   model1_cox_quart_danish, model2_cox_quart_danish, 
   model1_cox_sd_metanalysis, model2_cox_sd_metanalysis,
   model1_cox_quart_metanalysis, model2_cox_quart_metanalysis)

# Tables and figures ----
## Danish ----
### table covariates - als survival ----
covar_danish

### table POPs (sd) - als survival ----
POPs_sd_ALS_table_danish <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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

### table POPs (quart) - als survival ----
quartile1_rows <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "")

POPs_quart_ALS_table_danish <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status and education. 
  2Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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


### figure POPs (sd) - als survival ----
POPs_sd_ALS_figure_danish <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

### figure POPs (quart) - als survival ----
POPs_quart_ALS_figure_danish <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(!term == "Continuous") |>
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
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

### figure cumulative incidence POPs (quart) - als survival ----
# bdd_cases_danish <- bdd_danish |>                                               # set the datasets
#   filter(als == 1) |>                                                           # case selection
#   mutate(across(all_of(POPs_group), ~ factor(
#     ntile(.x, 4),                                                               # creation of POPs quartiles (cohort and cases specific)
#     labels = c("Q1", "Q2", "Q3", "Q4")), 
#     .names = "{.col}_quart")) 
# 
# create_surv_plot <- function(expl) {
#   formula <- as.formula(paste0("Surv(follow_up_death, status_death) ~ `", expl, "`"))
#   fit <- survfit(formula, data = bdd_cases_danish)
#   fit$call$formula <- formula
#   
#   plot <- ggsurvplot(
#     fit,
#     data = bdd_cases_danish,
#     fun = "event",
#     risk.table = TRUE,
#     pval = FALSE,
#     conf.int = FALSE,
#     palette = "Dark2",
#     xlab = "Follow-up (months)",
#     ylab = "Cumulative incidence",
#     legend.title = paste("Pre-disease", POPs_group_quart_labels[[expl]], "level"),
#     legend.labs = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")
#   )
#   
#   return(plot = plot)
# }
# 
# survival_plots_danish <- map(POPs_group_quart, create_surv_plot)
# names(survival_plots_danish) <- POPs_group_quart
# rm(create_surv_plot, bdd_cases_danish)

## Finnish ----
### table covariates - als survival ----
ref_rows <- covar_finnish |>
  distinct(model, variable) |>
  arrange(model, variable) |>
  mutate(
    term = c("continuous",  "continuous", "continuous", "1", "Other","Female","Never", "0",
             "continuous",  "continuous", "continuous", "1", "Other","Female","Never", "0"),
    HR = "-",
    "95% CI" = "-",
    `p-value` = "" ,
    lower_CI = "", 
    upper_CI = "")

covar_ALS_table_finnish <- covar_finnish |>
  mutate(
    term = fct_recode(term, 
                      "continuous" = "baseline_age",
                      "continuous" = "bmi",
                      "continuous" = "cholesterol",
                      "2" = "level_urbanization2",
                      "3" = "level_urbanization3",
                      "4" = "level_urbanization4",
                      "Married/cohabit" = "marital_status_2catMarried/cohabit",
                      "Male" = "sexMale",
                      "Ever" = "smoking_2catEver",
                      "1" = "thawed1"),
    HR = sprintf("%.1f", HR),
    lower_CI = sprintf("%.1f", lower_CI),
    upper_CI = sprintf("%.1f", upper_CI), 
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`, 
    `p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) 

covar_ALS_table_finnish <- 
  bind_rows(ref_rows, covar_ALS_table_finnish) |>
  select(model, variable, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI) |>
  mutate(
    term = fct_relevel(term, 
    "0", "1", "2", "3", "4", "continuous", "Never", "Ever", "Female", "Male", "Other", "Married/cohabit"), 
    model = fct_relevel(model, "base", "adjusted"))|>
  arrange(model, variable, term) |>
  filter(!(term == "continuous" & HR == "-"))

covar_ALS_table_finnish <- covar_ALS_table_finnish |>
  select(-`p-value_raw`, -`p-value_shape`, -lower_CI, -upper_CI) |>
  pivot_wider(names_from = model, values_from = c(HR, `95% CI`, `p-value`)) |>
  select(variable, term, ends_with("base"), ends_with("adjusted"))
rm(ref_rows)

### table POPs (sd) - als survival ----
POPs_sd_ALS_table_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age, sex, level of urbanization and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC study cohort). 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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

### table POPs (quart) - als survival ----
quartile1_rows <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "")

POPs_quart_ALS_table_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_finnish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_finnish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age, sex, level of urbanization and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC study cohort).
  2Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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


### figure POPs (sd) - als survival ----
POPs_sd_ALS_figure_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

### figure POPs (quart) - als survival ----
POPs_quart_ALS_figure_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
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
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

## Metanalysis ----
### table covariates - als survival ----

### table POPs (sd) - als survival ----
POPs_sd_ALS_table_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1Base models were all adjusted for age, sex. Finnish base models were futher adjusted for level of urbanization and serum freeze-thaw cycles. Adjusted models all further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC study cohort). 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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

### table POPs (quart) - als survival ----
quartile1_rows <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "")

POPs_quart_ALS_table_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_metanalysis <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_metanalysis) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, term, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
  "1Base models were all adjusted for age, sex. Finnish base models were futher adjusted for level of urbanization and serum freeze-thaw cycles. Adjusted models all further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC study cohort). 
  2Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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


### figure POPs (sd) - als survival ----
POPs_sd_ALS_figure_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  filter(term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

### figure POPs (quart) - als survival ----
POPs_quart_ALS_figure_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  filter(!term == "Continuous") |>
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
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

# Assemblage ----
results_POPs_ALS_survival <- 
  list(main_analysis = list(main_results_POPs_ALS_survival = main_results_POPs_ALS_survival), 
       danish = list(
         covar_danish = covar_danish, 
         POPs_sd_ALS_table_danish = POPs_sd_ALS_table_danish, 
         POPs_quart_ALS_table_danish = POPs_quart_ALS_table_danish, 
         POPs_sd_ALS_figure_danish = POPs_sd_ALS_figure_danish, 
         POPs_quart_ALS_figure_danish = POPs_quart_ALS_figure_danish), 
       finnish = list(
         covar_finnish = covar_finnish, 
         covar_ALS_table_finnish = covar_ALS_table_finnish,
         POPs_sd_ALS_table_finnish = POPs_sd_ALS_table_finnish, 
         POPs_quart_ALS_table_finnish = POPs_quart_ALS_table_finnish, 
         POPs_sd_ALS_figure_finnish = POPs_sd_ALS_figure_finnish, 
         POPs_quart_ALS_figure_finnish = POPs_quart_ALS_figure_finnish), 
       metanalysis = list(
         # covar_metanalysis = covar_metanalysis, 
         POPs_sd_ALS_table_metanalysis = POPs_sd_ALS_table_metanalysis, 
         POPs_quart_ALS_table_metanalysis = POPs_quart_ALS_table_metanalysis, 
         POPs_sd_ALS_figure_metanalysis = POPs_sd_ALS_figure_metanalysis, 
         POPs_quart_ALS_figure_metanalysis = POPs_quart_ALS_figure_metanalysis))

rm(main_results_POPs_ALS_survival, 
   covar_danish,
   POPs_sd_ALS_table_danish, 
   POPs_quart_ALS_table_danish, 
   POPs_sd_ALS_figure_danish, 
   POPs_quart_ALS_figure_danish, 
   
   covar_finnish,
   covar_ALS_table_finnish, 
   POPs_sd_ALS_table_finnish, 
   POPs_quart_ALS_table_finnish, 
   POPs_sd_ALS_figure_finnish, 
   POPs_quart_ALS_figure_finnish, 
   
   # covar_metanalysis,
   POPs_sd_ALS_table_metanalysis, 
   POPs_quart_ALS_table_metanalysis, 
   POPs_sd_ALS_figure_metanalysis, 
   POPs_quart_ALS_figure_metanalysis)

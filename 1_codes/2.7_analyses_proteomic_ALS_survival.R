# Aline Davias
# October 22, 2025 
# Analysis of als survival depending on proteomic profile

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# Creation of cases specific datasets ----
bdd_cases_danish <- 
  bdd_danish |>
  filter(als == 1) |>                                                           # remove controls 
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  select(sample, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic), 
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)                        
                labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)  
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)      
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
                  quart_meds[quartiles]       
                  },
                .names = "{.col}_quart_med"))

surv_obj <- Surv(time = bdd_cases_danish$follow_up_death,                       # set the outcomes
                 event = bdd_cases_danish$status_death)
covariates <-                                                                   # set the covariates 
  c("sex", "diagnosis_age", "smoking_2cat_i", "bmi")

var_label(bdd_cases_danish[proteomic]) <- NULL


# Effects of the covariates on ALS survival ----
covar <- tbl_merge(tbls = list(
  tbl_uvregression(
    data = bdd_cases_danish, 
    y = surv_obj, 
    method = survival::coxph,  
    exponentiate = TRUE,
    include = c("sex", "diagnosis_age", 
                "bmi", "marital_status_2cat_i", "smoking_i", "education_i", "cholesterol_i")) |>
    bold_labels() |>
    bold_p(), 
  tbl_regression(
    coxph(surv_obj ~ sex + diagnosis_age + bmi + marital_status_2cat_i + smoking_i + education_i + cholesterol_i, data = bdd_cases_danish),
    exponentiate = TRUE) |>
    bold_labels() |>
    bold_p()), 
  tab_spanner = c("**Crude**", "**Adjusted**"))

# Main analysis ----
## Cox model (sd) ----
### Base ----
model1_cox_sd <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj ~", expl, "+ diagnosis_age + sex"))              # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  
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
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj ~", expl, "+",  
                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
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
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj ~", expl, "+ diagnosis_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  
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
  
  formula_danish <- as.formula(paste("surv_obj ~", expl, "+",                   # set the formulas              
                                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
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
    
    data_sub <- bdd_cases_danish |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "main",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex"))  
    test_2 <- coxph(formula, data = data_sub)  
    
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
    
    data_sub <- bdd_cases_danish |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "main",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
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
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish) |> summary()
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
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish) |> summary()
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
# direct dans output 


# Sensi 1 + sensi 3 - Removing the oulier for NEFL + filtering cases with follow_up > 5 years ----
bdd_cases_danish_sensi_1_3 <- 
  bdd_danish |>
  filter (als == 1 & follow_up > 60) |>                                         # filtering cases with follow_up > 5 years
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  select(sample, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic), 
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)                        
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)  
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)      
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
                  quart_meds[quartiles]       
                },
                .names = "{.col}_quart_med"))

surv_obj_sensi_1_3 <- Surv(time = bdd_cases_danish_sensi_1_3$follow_up_death,   # set the outcomes
                           event = bdd_cases_danish_sensi_1_3$status_death)

## Cox model (sd) ----
### Base 
model1_cox_sd_sensi_1_3 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ diagnosis_age + sex"))              # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "sensi_1_3", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Adjusted 
model2_cox_sd_sensi_1_3 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+",  
                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

## Cox model (quart) ----
### Base ----
model1_cox_quart_sensi_1_3 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ diagnosis_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "sensi_1_3", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Adjusted ----
model2_cox_quart_sensi_1_3 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+",                   # set the formulas              
                                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Heterogeneity tests ----
#### base 

heterogeneity_base_sensi_1_3 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3",
      p_value_heterogeneity = p_lr)
  })



#### adjusted 
heterogeneity_adjusted_sensi_1_3 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                           # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3",
      p_value_heterogeneity = p_lr)
  })


heterogeneity_tests_sensi_1_3 <- 
  bind_rows(heterogeneity_base_sensi_1_3, 
            heterogeneity_adjusted_sensi_1_3) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base 
trend_base_sensi_1_3 <- data.frame(explanatory = character(),
                                   model = factor(), 
                                   p_value_trend = numeric(), 
                                   stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base_sensi_1_3 <- rbind(trend_base_sensi_1_3, 
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "sensi_1_3", 
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted 
trend_adjusted_sensi_1_3 <- data.frame(explanatory = character(),
                                       model = factor(), 
                                       p_value_trend = numeric(), 
                                       stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3 <- rbind(trend_adjusted_sensi_1_3, 
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "sensi_1_3", 
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests_sensi_1_3 <- 
  bind_rows(trend_base_sensi_1_3, trend_adjusted_sensi_1_3) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base_sensi_1_3, heterogeneity_adjusted_sensi_1_3,
   trend_base_sensi_1_3, trend_adjusted_sensi_1_3)


## Cox model (GAMs) ----
# direct dans output 


# Sensi 1 + sensi 3 + sensi 4 - Removing NEFL outlier + filtering cases with follow- up < 5 years + filtering follow-up <= 50%----
bdd_cases_danish_sensi_1_3_4 <- bdd_danish |>
  filter(als == 1 & follow_up > 60) |>                                          # remove controls and remove cases with follow-up<5years 
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 4 : we remove 50% of the highest values of follow-up (after removing the follow-up<5years)
  filter(follow_up <= seuil) |>
  select(sample, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic), 
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)                        
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)  
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)      
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
                  quart_meds[quartiles]       
                },
                .names = "{.col}_quart_med"))

surv_obj_sensi_1_3_4 <- Surv(time = bdd_cases_danish_sensi_1_3_4$follow_up_death,           # set the outcomes
                             event = bdd_cases_danish_sensi_1_3_4$status_death)

## Cox model (sd) ----
### Base ----
model1_cox_sd_sensi_1_3_4 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ diagnosis_age + sex"))              # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_4) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "sensi_1_3_4", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Adjusted ----
model2_cox_sd_sensi_1_3_4 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+",  
                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_4) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3_4", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

## Cox model (quart) ----
### Base ----
model1_cox_quart_sensi_1_3_4 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ diagnosis_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_4) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base",
    analysis = "sensi_1_3_4", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Adjusted ----
model2_cox_quart_sensi_1_3_4 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+",                   # set the formulas              
                                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_4) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3_4", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Heterogeneity tests ----
#### base 

heterogeneity_base_sensi_1_3_4 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3_4 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3_4",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3_4",
      p_value_heterogeneity = p_lr)
  })



#### adjusted 
heterogeneity_adjusted_sensi_1_3_4 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3_4 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3_4",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                           # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3_4",
      p_value_heterogeneity = p_lr)
  })

heterogeneity_tests_sensi_1_3_4 <- 
  bind_rows(heterogeneity_base_sensi_1_3_4, 
            heterogeneity_adjusted_sensi_1_3_4) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base 
trend_base_sensi_1_3_4 <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_4) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base_sensi_1_3_4 <- rbind(trend_base_sensi_1_3_4, 
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "sensi_1_3_4", 
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted 
trend_adjusted_sensi_1_3_4 <- data.frame(explanatory = character(),
                                         model = factor(), 
                                         p_value_trend = numeric(), 
                                         stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_4) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3_4 <- rbind(trend_adjusted_sensi_1_3_4, 
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "sensi_1_3_4", 
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests_sensi_1_3_4 <- 
  bind_rows(trend_base_sensi_1_3_4, trend_adjusted_sensi_1_3_4) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base_sensi_1_3_4, heterogeneity_adjusted_sensi_1_3_4,
   trend_base_sensi_1_3_4, trend_adjusted_sensi_1_3_4)


## Cox model (GAMs) ----
# direct dans output 


# Sensi 1 + sensi 3 + sensi 5 - Removing NEFL outlier + filtering cases with follow-up < 5 years + filtering follow-up > 50%----
bdd_cases_danish_sensi_1_3_5 <- bdd_danish |>
  filter(als == 1 & follow_up > 60) |>                                          # remove controls and remove cases with follow-up<5years 
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 4 : we remove 50% of the highest values of follow-up (after removing the follow-up<5years)
  filter(follow_up > seuil) |>
  select(sample, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic), 
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)                        
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)  
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)      
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
                  quart_meds[quartiles]       
                },
                .names = "{.col}_quart_med"))

surv_obj_sensi_1_3_5 <- Surv(time = bdd_cases_danish_sensi_1_3_5$follow_up_death,           # set the outcomes
                             event = bdd_cases_danish_sensi_1_3_5$status_death)

## Cox model (sd) ----
### Base ----
model1_cox_sd_sensi_1_3_5 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ diagnosis_age + sex"))              # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_5) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "sensi_1_3_5", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Adjusted ----
model2_cox_sd_sensi_1_3_5 <- map_dfr(proteomic_sd, function(expl) {
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+",  
                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_5) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3_5", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

## Cox model (quart) ----
### Base ----
model1_cox_quart_sensi_1_3_5 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ diagnosis_age + sex"))  
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_5) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base", 
    analysis = "sensi_1_3_5", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Adjusted ----
model2_cox_quart_sensi_1_3_5 <- map_dfr(proteomic_quart, function(expl) {
  
  formula_danish <- as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+",                   # set the formulas              
                                     paste(covariates, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_3_5) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted", 
    analysis = "sensi_1_3_5", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### Heterogeneity tests ----
#### base 
heterogeneity_base_sensi_1_3_5 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3_5 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3_5",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3_5",
      p_value_heterogeneity = p_lr)
  })



#### adjusted 
heterogeneity_adjusted_sensi_1_3_5 <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_3_5 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, sex, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3_5",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                           # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + sex + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3_5",
      p_value_heterogeneity = p_lr)
  })

heterogeneity_tests_sensi_1_3_5 <- 
  bind_rows(heterogeneity_base_sensi_1_3_5, 
            heterogeneity_adjusted_sensi_1_3_5) |>
  mutate(explanatory = gsub("_quart", "", explanatory))


### Trend tests ----
#### base ----
trend_base_sensi_1_3_5 <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_5) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base_sensi_1_3_5 <- rbind(trend_base_sensi_1_3_5, 
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "sensi_1_3_5", 
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted ----
trend_adjusted_sensi_1_3_5 <- data.frame(explanatory = character(),
                                         model = factor(), 
                                         p_value_trend = numeric(), 
                                         stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_5) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3_5 <- rbind(trend_adjusted_sensi_1_3_5, 
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "sensi_1_3_5", 
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests_sensi_1_3_5 <- 
  bind_rows(trend_base_sensi_1_3_5, trend_adjusted_sensi_1_3_5) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base_sensi_1_3_5, heterogeneity_adjusted_sensi_1_3_5,
   trend_base_sensi_1_3_5, trend_adjusted_sensi_1_3_5)


## Cox model (GAMs) ----
# direct dans output


# Sensi 1 + sensi 7 - Removing the oulier for NEFL + filtering to females ----
bdd_cases_danish_sensi_1_7_female <- 
  bdd_danish |>
  filter(als == 1) |>                                                           # remove controls
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  filter(sex == "Female") |>                                                    # sensi 7 female : we keep only female cases
  select(sample, als, als_date, follow_up_death, status_death, baseline_age, diagnosis_age, death_age, follow_up,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i,
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic),
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)
                  quart_meds[quartiles]
                },
                .names = "{.col}_quart_med"))

surv_obj_sensi_1_7_female <- Surv(time = bdd_cases_danish_sensi_1_7_female$follow_up_death,           # set the outcomes
                                  event = bdd_cases_danish_sensi_1_7_female$status_death)
covariates_sensi_7 <-                                                                   # set the covariates_sensi_7
  c("diagnosis_age", "smoking_2cat_i", "bmi")

## Cox model (sd) ----
### Base ----
model1_cox_sd_sensi_1_7_female <- map_dfr(proteomic_sd, function(expl) {

  formula_danish <-
    as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ diagnosis_age"))              # set the formulas
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_female) |> summary()  # run cox model

  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base",
    analysis = "sensi_1_7_female",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Adjusted ----
model2_cox_sd_sensi_1_7_female <- map_dfr(proteomic_sd, function(expl) {

  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+",
                     paste(covariates_sensi_7, collapse = " + ")))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_female) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted",
    analysis = "sensi_1_7_female",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

## Cox model (quart) ----
### Base ----
model1_cox_quart_sensi_1_7_female <- map_dfr(proteomic_quart, function(expl) {

  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ diagnosis_age"))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_female) |> summary()  # run cox model

  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base",
    analysis = "sensi_1_7_female",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Adjusted ----
model2_cox_quart_sensi_1_7_female <- map_dfr(proteomic_quart, function(expl) {

  formula_danish <- as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+",                   # set the formulas
                                     paste(covariates_sensi_7, collapse = " + ")))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_female) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted",
    analysis = "sensi_1_7_female",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Heterogeneity tests ----
#### base 

heterogeneity_base_sensi_1_7_female <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_7_female |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_7_female",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_7_female",
      p_value_heterogeneity = p_lr)
  })



#### adjusted 
heterogeneity_adjusted_sensi_1_7_female <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_7_female |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_7_female",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                           # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_7_female",
      p_value_heterogeneity = p_lr)
  })


heterogeneity_tests_sensi_1_7_female <-
  bind_rows(heterogeneity_base_sensi_1_7_female,
            heterogeneity_adjusted_sensi_1_7_female) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base 
trend_base_sensi_1_7_female <- data.frame(explanatory = character(),
                                          model = factor(),
                                          p_value_trend = numeric(),
                                          stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {

  formula <- as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ diagnosis_age"))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_female) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]

  trend_base_sensi_1_7_female <- rbind(trend_base_sensi_1_7_female,
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "sensi_1_7_female",
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted 
trend_adjusted_sensi_1_7_female <- data.frame(explanatory = character(),
                                              model = factor(),
                                              p_value_trend = numeric(),
                                              stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {

  formula <- as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ ", paste(covariates_sensi_7, collapse = "+")))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_female) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]

  trend_adjusted_sensi_1_7_female <- rbind(trend_adjusted_sensi_1_7_female,
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "sensi_1_7_female",
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests_sensi_1_7_female <-
  bind_rows(trend_base_sensi_1_7_female, trend_adjusted_sensi_1_7_female) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base_sensi_1_7_female, heterogeneity_adjusted_sensi_1_7_female,
   trend_base_sensi_1_7_female, trend_adjusted_sensi_1_7_female)


## Cox model (GAMs) ----
# direct output 

# Sensi 1 + sensi 7 - Removing the oulier for NEFL + filtering to male ----
bdd_cases_danish_sensi_1_7_male <- bdd_danish |>
  filter(als == 1) |>                                                           # remove controls
  mutate(                                                                       # remove NEFL outlier
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>   
  filter(sex == "Male") |>                                                    # sensi 7 male : we keep only male cases
  select(sample, als, als_date, follow_up_death, status_death, baseline_age, diagnosis_age, death_age, follow_up,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i,
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic),
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)
                  quart_meds[quartiles]
                },
                .names = "{.col}_quart_med"))

surv_obj_sensi_1_7_male <- Surv(time = bdd_cases_danish_sensi_1_7_male$follow_up_death,           # set the outcomes
                                event = bdd_cases_danish_sensi_1_7_male$status_death)
covariates_sensi_7 <-                                                                   # set the covariates_sensi_7
  c("diagnosis_age", "smoking_2cat_i", "bmi")

## Cox model (sd) ----
### Base ----
model1_cox_sd_sensi_1_7_male <- map_dfr(proteomic_sd, function(expl) {

  formula_danish <-
    as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ diagnosis_age"))              # set the formulas
  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_male) |> summary()  # run cox model

  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base",
    analysis = "sensi_1_7_male",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Adjusted ----
model2_cox_sd_sensi_1_7_male <- map_dfr(proteomic_sd, function(expl) {

  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+",
                     paste(covariates_sensi_7, collapse = " + ")))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_male) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted",
    analysis = "sensi_1_7_male",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

## Cox model (quart) ----
### Base ----
model1_cox_quart_sensi_1_7_male <- map_dfr(proteomic_quart, function(expl) {

  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ diagnosis_age"))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_male) |> summary()  # run cox model

  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "base",
    analysis = "sensi_1_7_male",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Adjusted ----
model2_cox_quart_sensi_1_7_male <- map_dfr(proteomic_quart, function(expl) {

  formula_danish <- as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+",                   # set the formulas
                                     paste(covariates_sensi_7, collapse = " + ")))

  model_summary <- coxph(formula_danish, data = bdd_cases_danish_sensi_1_7_male) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    model = "adjusted",
    analysis = "sensi_1_7_male",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    p_value = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates_sensi_7 results
})

### Heterogeneity tests ----
#### base 

heterogeneity_base_sensi_1_7_male <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_7_male |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_7_male",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                       # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_7_male",
      p_value_heterogeneity = p_lr)
  })



#### adjusted 
heterogeneity_adjusted_sensi_1_7_male <- map_dfr(
  proteomic_quart,
  function(var) {
    
    data_sub <- bdd_cases_danish_sensi_1_7_male |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, diagnosis_age, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_7_male",
        p_value_heterogeneity = NA_real_))
    }
    
    surv_obj <- Surv(time = data_sub$follow_up_death,                           # set the outcomes
                     event = data_sub$status_death)
    
    formula_raw <- as.formula("surv_obj ~ diagnosis_age + smoking_2cat_i + bmi")
    test_1 <- coxph(formula_raw, data = data_sub) 
    
    formula <- as.formula(paste("surv_obj ~", var, "+ diagnosis_age + smoking_2cat_i + bmi"))  
    test_2 <- coxph(formula, data = data_sub)  
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_7_male",
      p_value_heterogeneity = p_lr)
  })

heterogeneity_tests_sensi_1_7_male <-
  bind_rows(heterogeneity_base_sensi_1_7_male,
            heterogeneity_adjusted_sensi_1_7_male) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base 
trend_base_sensi_1_7_male <- data.frame(explanatory = character(),
                                        model = factor(),
                                        p_value_trend = numeric(),
                                        stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {

  formula <- as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ diagnosis_age"))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_male) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]

  trend_base_sensi_1_7_male <- rbind(trend_base_sensi_1_7_male,
                      data.frame(explanatory = expl,
                                 model = "base",
                                 analysis = "sensi_1_7_male",
                                 p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

#### adjusted 
trend_adjusted_sensi_1_7_male <- data.frame(explanatory = character(),
                                            model = factor(),
                                            p_value_trend = numeric(),
                                            stringsAsFactors = FALSE)

for (expl in proteomic_quart_med) {

  formula <- as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ ", paste(covariates_sensi_7, collapse = "+")))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_male) |> summary()
  p_value_trend <- model$coefficients[expl, "Pr(>|z|)"]

  trend_adjusted_sensi_1_7_male <- rbind(trend_adjusted_sensi_1_7_male,
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     analysis = "sensi_1_7_male",
                                     p_value_trend = p_value_trend))
}
rm(expl, model, formula, p_value_trend)

trend_tests_sensi_1_7_male <-
  bind_rows(trend_base_sensi_1_7_male, trend_adjusted_sensi_1_7_male) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_base_sensi_1_7_male, heterogeneity_adjusted_sensi_1_7_male,
   trend_base_sensi_1_7_male, trend_adjusted_sensi_1_7_male)


## Cox model (GAMs) ----
# direct output


# Assemblage main analyses ----
main_results <-       
  bind_rows(
    model1_cox_sd, model2_cox_sd,
    model1_cox_quart, model2_cox_quart, 
    model1_cox_sd_sensi_1_3, model2_cox_sd_sensi_1_3,
    model1_cox_quart_sensi_1_3, model2_cox_quart_sensi_1_3, 
    model1_cox_sd_sensi_1_3_4, model2_cox_sd_sensi_1_3_4,
    model1_cox_quart_sensi_1_3_4, model2_cox_quart_sensi_1_3_4, 
    model1_cox_sd_sensi_1_3_5, model2_cox_sd_sensi_1_3_5,
    model1_cox_quart_sensi_1_3_5, model2_cox_quart_sensi_1_3_5,
    model1_cox_sd_sensi_1_7_female, model2_cox_sd_sensi_1_7_female,
    model1_cox_quart_sensi_1_7_female, model2_cox_quart_sensi_1_7_female,
    model1_cox_sd_sensi_1_7_male, model2_cox_sd_sensi_1_7_male,
    model1_cox_quart_sensi_1_7_male, model2_cox_quart_sensi_1_7_male) |>
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
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    p_value_raw = p_value, 
    p_value_shape = ifelse(p_value_raw<0.05, "p_value<0.05", "p_value≥0.05"), 
    p_value = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
    p_value = ifelse(p_value == "1.00", ">0.99", p_value)) |>
  select(model, explanatory, analysis, term,  HR, HR_raw, `95% CI`, p_value, p_value_raw, p_value_shape, lower_CI, upper_CI)

heterogeneity_tests <- bind_rows(heterogeneity_tests, 
                                 heterogeneity_tests_sensi_1_3, 
                                 heterogeneity_tests_sensi_1_3_4, 
                                 heterogeneity_tests_sensi_1_3_5, 
                                 heterogeneity_tests_sensi_1_7_female,
                                 heterogeneity_tests_sensi_1_7_male)

trend_tests <- bind_rows(trend_tests, 
                                 trend_tests_sensi_1_3, 
                                 trend_tests_sensi_1_3_4, 
                                 trend_tests_sensi_1_3_5, 
                                 trend_tests_sensi_1_7_female,
                                 trend_tests_sensi_1_7_male)

main_results <- 
  main_results |>
  left_join(heterogeneity_tests, by = c("explanatory", "analysis", "model")) |>
  mutate(p_value_heterogeneity = ifelse(term == "Continuous", NA, p_value_heterogeneity),
         p_value_heterogeneity = ifelse(p_value_heterogeneity < 0.01, "<0.01", number(p_value_heterogeneity, accuracy = 0.01, decimal.mark = ".")),
         p_value_heterogeneity = ifelse(p_value_heterogeneity == "1.00", ">0.99", p_value_heterogeneity))

main_results <- 
  main_results |> 
  left_join(trend_tests, by = c("explanatory", "analysis", "model")) |>
  mutate(p_value_trend = ifelse(term == "Continuous", NA, p_value_trend),
         p_value_trend = ifelse(p_value_trend < 0.01, "<0.01", number(p_value_trend, accuracy = 0.01, decimal.mark = ".")),
         p_value_trend = ifelse(p_value_trend == "1.00", ">0.99", p_value_trend), 
         protein_group = case_when(str_detect(explanatory, 'proteomic_immun_res') ~ "Immune response", 
                                   str_detect(explanatory, 'proteomic_metabolism') ~ "Metabolism", 
                                   str_detect(explanatory, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
         explanatory = str_replace(explanatory, 'proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_', ""))

rm(model1_cox_sd, model2_cox_sd, 
   model1_cox_quart, model2_cox_quart, 
   model1_cox_sd_sensi_1_3, model2_cox_sd_sensi_1_3, 
   model1_cox_quart_sensi_1_3, model2_cox_quart_sensi_1_3, 
   model1_cox_sd_sensi_1_3_4, model2_cox_sd_sensi_1_3_4, 
   model1_cox_quart_sensi_1_3_4, model2_cox_quart_sensi_1_3_4, 
   model1_cox_sd_sensi_1_3_5, model2_cox_sd_sensi_1_3_5, 
   model1_cox_quart_sensi_1_3_5, model2_cox_quart_sensi_1_3_5, 
   model1_cox_sd_sensi_1_7_female, model2_cox_sd_sensi_1_7_female,
   model1_cox_quart_sensi_1_7_female, model2_cox_quart_sensi_1_7_female,
   model1_cox_sd_sensi_1_7_male, model2_cox_sd_sensi_1_7_male,
   model1_cox_quart_sensi_1_7_male, model2_cox_quart_sensi_1_7_male,
   heterogeneity_tests, trend_tests, 
   heterogeneity_tests_sensi_1_3, trend_tests_sensi_1_3, 
   heterogeneity_tests_sensi_1_3_4, trend_tests_sensi_1_3_4,
   heterogeneity_tests_sensi_1_3_5, trend_tests_sensi_1_3_5,
   heterogeneity_tests_sensi_1_7_female, trend_tests_sensi_1_7_female,
   heterogeneity_tests_sensi_1_7_male, trend_tests_sensi_1_7_male)


# Tables and Figures ----

## Table covariates - als survival ----
covar

## Table proteomic (sd) - als survival (main) ----
proteomic_sd_ALS_table <- 
  main_results |>
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (main) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table <- 
  proteomic_quart_ALS_table |>
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


## Figure proteomic - als survival - base sd (main) ----
proteomic_sd_ALS_base_figure <- 
  main_results |>
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

## Figure proteomic - als survival - adjusted sd (main) ----
proteomic_sd_ALS_adjusted_figure <- 
  main_results |>
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




## Table proteomic (sd) - als survival (sensi_1_3) ----
proteomic_sd_ALS_table_sensi_1_3 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "sensi_1_3") |>            # select continuous results
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (sensi_1_3) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_3 <- 
  proteomic_quart_ALS_table_sensi_1_3 |>
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


## Figure proteomic - als survival - base sd (sensi_1_3) ----
proteomic_sd_ALS_base_figure_sensi_1_3 <- 
  main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "sensi_1_3") |>
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

## Figure proteomic - als survival - adjusted sd (sensi_1_3) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3 <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "sensi_1_3") |>
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



## Table proteomic (sd) - als survival (sensi_1_3_4) ----
proteomic_sd_ALS_table_sensi_1_3_4 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "sensi_1_3_4") |>            # select continuous results
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (sensi_1_3_4) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3_4") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3_4 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3_4") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_3_4 <- 
  proteomic_quart_ALS_table_sensi_1_3_4 |>
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


## Figure proteomic - als survival - base sd (sensi_1_3_4) ----
proteomic_sd_ALS_base_figure_sensi_1_3_4 <- 
  main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "sensi_1_3_4") |>
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

## Figure proteomic - als survival - adjusted sd (sensi_1_3_4) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "sensi_1_3_4") |>
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




## Table proteomic (sd) - als survival (sensi_1_3_5) ----
proteomic_sd_ALS_table_sensi_1_3_5 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "sensi_1_3_5") |>            # select continuous results
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (sensi_1_3_5) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3_5") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3_5 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_3_5") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_3_5 <- 
  proteomic_quart_ALS_table_sensi_1_3_5 |>
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


## Figure proteomic - als survival - base sd (sensi_1_3_5) ----
proteomic_sd_ALS_base_figure_sensi_1_3_5 <- 
  main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "sensi_1_3_5") |>
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

## Figure proteomic - als survival - adjusted sd (sensi_1_3_5) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "sensi_1_3_5") |>
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


## Table proteomic (sd) - als survival (sensi_1_7_female) ----
proteomic_sd_ALS_table_sensi_1_7_female <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "sensi_1_7_female") |>            # select continuous results
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (sensi_1_7_female) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_7_female") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_7_female <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_7_female") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_7_female <- 
  proteomic_quart_ALS_table_sensi_1_7_female |>
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


## Figure proteomic - als survival - base sd (sensi_1_7_female) ----
proteomic_sd_ALS_base_figure_sensi_1_7_female <- 
  main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "sensi_1_7_female") |>
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

## Figure proteomic - als survival - adjusted sd (sensi_1_7_female) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_7_female <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "sensi_1_7_female") |>
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





## Table proteomic (sd) - als survival (sensi_1_7_male) ----
proteomic_sd_ALS_table_sensi_1_7_male <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous" & analysis == "sensi_1_7_male") |>            # select continuous results
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
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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


## Table proteomic (quart) - als survival (sensi_1_7_male) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_7_male") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_7_male <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous" & analysis == "sensi_1_7_male") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_7_male <- 
  proteomic_quart_ALS_table_sensi_1_7_male |>
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



## Figure proteomic - als survival - base sd (sensi_1_7_male) ----
proteomic_sd_ALS_base_figure_sensi_1_7_male <- 
  main_results |>
  filter(model == "base" & term == "Continuous" & analysis == "sensi_1_7_male") |>
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

## Figure proteomic - als survival - adjusted sd (sensi_1_7_male) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_7_male <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous" & analysis == "sensi_1_7_male") |>
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



## NfL only results ----
### Table NfL - als occurrence - base and adjusted sd ----
NfL_sd_ALS_table <- main_results |>
  filter(model %in% c("base", "adjusted") & 
           term == "Continuous" & 
           analysis == "main" &
           explanatory == "NEFL") |>            
  select(model, explanatory, protein_group, term, HR, "95% CI", "p_value") |>
  arrange(protein_group, explanatory) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1Base model is adjusted for age at diagnosis and sex. Adjusted model further accounts for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
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

### Table NfL - als occurrence - base and adjusted quart ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & 
           term != "Continuous" & 
           analysis == "main" & 
           explanatory == "NEFL") |>   
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "HR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

NfL_quart_ALS_table <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & 
           term != "Continuous" & 
           analysis == "main" & 
           explanatory == "NEFL") |>    
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

NfL_quart_ALS_table <- 
  NfL_quart_ALS_table |>
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


### Table NfL - als occurrence - base and adjusted sd (sensi follow up) ----
NfL_sd_ALS_table_sensi_follow_up <- main_results |>
  filter(analysis %in% c("main", "sensi_1_3", "sensi_1_3_4", "sensi_1_3_5"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(analysis = fct_recode(analysis, 
                               "Main analysis (n=165)" = "main",
                               "Filtered to follow-up > 5 years (n=148)" = "sensi_1_3",
                               "Filtered to follow-up between 5 and 14.6 years (n=74)" = "sensi_1_3_4",
                               "Filtered to follow-up > 14.6 years (n=74)" = "sensi_1_3_5"), 
         analysis = fct_relevel(analysis, 
                                "Main analysis (n=165)", 
                                "Filtered to follow-up > 5 years (n=148)", 
                                "Filtered to follow-up between 5 and 14.6 years (n=74)",
                                "Filtered to follow-up > 14.6 years (n=74)")) |> 
  select(analysis, model, explanatory, HR, "95% CI", p_value) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(HR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 'HR' = 'HR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'HR ' = 'HR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and birth year. Adjusted models further account for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "analysis" = "Analyses",
    "explanatory" = "Explanatory variable", 
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
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

### Figure NfL - als occurrence - base and adjusted sd (sensi follow up) ----
NfL_sd_ALS_figure_sensi_follow_up <- main_results |>
  filter(analysis %in% c("main", "sensi_1_3", "sensi_1_3_4", "sensi_1_3_5"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted models" = "adjusted",
                            "Base models" = "base"), 
         model = fct_relevel(model, "Base models", "Adjusted models"), 
         analysis = fct_recode(analysis, 
                               "Main analysis\n(n=165)" = "main",
                               "Filtered to\nfollow-up > 5 years\n (n=148)" = "sensi_1_3",
                               "Filtered to follow-up\nbetween 5 and 14.6 years\n(n=74)" = "sensi_1_3_4",
                               "Filtered to\nfollow-up > 14.6 years\n (n=74)" = "sensi_1_3_5"), 
         analysis = fct_relevel(analysis, 
                                "Main analysis\n(n=165)", 
                                "Filtered to\nfollow-up > 5 years\n (n=148)", 
                                "Filtered to follow-up\nbetween 5 and 14.6 years\n(n=74)",
                                "Filtered to\nfollow-up > 14.6 years\n (n=74)")) |> 
  ggplot(aes(x = explanatory, y = HR_raw, ymin = lower_CI, ymax = upper_CI, color = signif)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis), cols = dplyr::vars(model), switch = "y") +                         # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs( y = "Hazard ratios (HRs)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()

### Table NfL - als occurrence - base and adjusted sd (sensi sex) ----
NfL_sd_ALS_table_sensi_sex <- main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis %in% c("main", "sensi_1_7_female", "sensi_1_7_male") &
           explanatory == "NEFL") |>         
  select(analysis, explanatory, protein_group, term, HR, "95% CI", "p_value") |>
  arrange(protein_group, explanatory) |>
  pivot_wider(names_from = "analysis", values_from = c("HR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("main"), contains("sensi_1_7_female"), contains("sensi_1_7_male")) |>
  rename("HR" = "HR_main", "95% CI" = "95% CI_main", "p-value" = "p_value_main", 
         "HR " = "HR_sensi_1_7_female", "95% CI " = "95% CI_sensi_1_7_female", "p-value " = "p_value_sensi_1_7_female", 
         " HR " = "HR_sensi_1_7_male", " 95% CI " = "95% CI_sensi_1_7_male", " p-value " = "p_value_sensi_1_7_male") |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of proteins. 
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "HR" = "Main analysis\n(females and males)", "95% CI" = "Main analysis\n(females and males)", "p-value" = "Main analysis\n(females and males)", 
    "HR " = "Females", "95% CI " = "Females", "p-value " = "Females", 
    " HR " = "Males", " 95% CI " = "Males", " p-value " = "Males") |>
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

# Assemblage ----
results_proteomic_ALS_survival <- 
  list(
    main = list(covar = covar, 
                main_results= main_results, 
                
                proteomic_sd_ALS_table = proteomic_sd_ALS_table,
                proteomic_quart_ALS_table = proteomic_quart_ALS_table,
                
                proteomic_sd_ALS_base_figure = proteomic_sd_ALS_base_figure, 
                proteomic_sd_ALS_adjusted_figure = proteomic_sd_ALS_adjusted_figure), 
    
    sensi_1_3 = list(proteomic_sd_ALS_table_sensi_1_3 = proteomic_sd_ALS_table_sensi_1_3,
                     proteomic_quart_ALS_table_sensi_1_3 = proteomic_quart_ALS_table_sensi_1_3,
                     
                     proteomic_sd_ALS_base_figure_sensi_1_3 = proteomic_sd_ALS_base_figure_sensi_1_3, 
                     proteomic_sd_ALS_adjusted_figure_sensi_1_3 = proteomic_sd_ALS_adjusted_figure_sensi_1_3), 
    
    sensi_1_3_4 = list(proteomic_sd_ALS_table_sensi_1_3_4 = proteomic_sd_ALS_table_sensi_1_3_4,
                     proteomic_quart_ALS_table_sensi_1_3_4 = proteomic_quart_ALS_table_sensi_1_3_4,
                     
                     proteomic_sd_ALS_base_figure_sensi_1_3_4 = proteomic_sd_ALS_base_figure_sensi_1_3_4, 
                     proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_4), 
    
    sensi_1_3_5 = list(proteomic_sd_ALS_table_sensi_1_3_5 = proteomic_sd_ALS_table_sensi_1_3_5,
                       proteomic_quart_ALS_table_sensi_1_3_5 = proteomic_quart_ALS_table_sensi_1_3_5,
                       
                       proteomic_sd_ALS_base_figure_sensi_1_3_5 = proteomic_sd_ALS_base_figure_sensi_1_3_5, 
                       proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_5), 
    
    sensi_1_7_female = list(proteomic_sd_ALS_table_sensi_1_7_female = proteomic_sd_ALS_table_sensi_1_7_female,
                            proteomic_quart_ALS_table_sensi_1_7_female = proteomic_quart_ALS_table_sensi_1_7_female,
                            
                            proteomic_sd_ALS_base_figure_sensi_1_7_female = proteomic_sd_ALS_base_figure_sensi_1_7_female, 
                            proteomic_sd_ALS_adjusted_figure_sensi_1_7_female = proteomic_sd_ALS_adjusted_figure_sensi_1_7_female), 
    
    sensi_1_7_male = list(proteomic_sd_ALS_table_sensi_1_7_male = proteomic_sd_ALS_table_sensi_1_7_male,
                          proteomic_quart_ALS_table_sensi_1_7_male = proteomic_quart_ALS_table_sensi_1_7_male,
                          
                          proteomic_sd_ALS_base_figure_sensi_1_7_male = proteomic_sd_ALS_base_figure_sensi_1_7_male, 
                          proteomic_sd_ALS_adjusted_figure_sensi_1_7_male = proteomic_sd_ALS_adjusted_figure_sensi_1_7_male), 
    
    NfL_results = list(NfL_sd_ALS_table = NfL_sd_ALS_table, 
                       NfL_quart_ALS_table = NfL_quart_ALS_table, 
                       NfL_sd_ALS_table_sensi_follow_up = NfL_sd_ALS_table_sensi_follow_up, 
                       NfL_sd_ALS_figure_sensi_follow_up = NfL_sd_ALS_figure_sensi_follow_up, 
                       NfL_sd_ALS_table_sensi_sex = NfL_sd_ALS_table_sensi_sex))


saveRDS(results_proteomic_ALS_survival, file = "~/Documents/POP_ALS_2025_02_03/2_output/2.7_results_proteomic_ALS_survival.rds")


rm(bdd_cases_danish, 
   covariates, 
   covar, 
   main_results, 
   surv_obj, 
   proteomic_sd_ALS_table,
   proteomic_quart_ALS_table, 
   proteomic_sd_ALS_base_figure, 
   proteomic_sd_ALS_adjusted_figure, 
   NfL_sd_ALS_table,
   NfL_quart_ALS_table,
   NfL_sd_ALS_table_sensi_follow_up,
   NfL_sd_ALS_figure_sensi_follow_up,
   NfL_sd_ALS_table_sensi_sex)

rm(list = ls(pattern = "sensi_1_3"))
rm(list = ls(pattern = "sensi_1_3_4"))
rm(list = ls(pattern = "sensi_1_3_5"))
rm(list = ls(pattern = "sensi_1_7_female"))
rm(list = ls(pattern = "sensi_1_7_male"))

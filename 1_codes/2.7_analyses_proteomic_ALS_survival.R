# Aline Davias
# October 22, 2025 
# Analysis of als survival depending on proteomic profile

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.6_analyses_proteomic_ALS_occurrence.R", echo=TRUE)

# Creation of cases specific datasets ----
bdd_cases_danish <- bdd_danish |>
  filter(als == 1) |>   # remove controls 
  filter(match != 159) |>    # remove NEFL outlier
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
#### base ----
heterogeneity_base_quart <- data.frame(explanatory = character(),
                                       model = factor(),
                                       p_value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula("surv_obj ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish) 
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart <- rbind(heterogeneity_base_quart, 
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "main", 
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart <- data.frame(explanatory = character(),
                                           model = factor(),
                                           p_value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula(paste("surv_obj ~ ", paste(covariates, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish) 
  
  formula <- as.formula(paste("surv_obj ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "main", 
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_quart, 
            heterogeneity_adjusted_quart) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base ----
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

#### adjusted ----
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

rm(heterogeneity_base_quart, heterogeneity_adjusted_quart,
   trend_base, trend_adjusted)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)


### Base model ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish) {

  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,  
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted model ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var, 
    x_label = x_label)
}

cox_gam_results_adjusted <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels)

# Sensi 1 + sensi 3 - Removing the oulier for NEFL + filtering cases with follow_up > 5 years ----
bdd_cases_danish_sensi_1_3 <- 
  bdd_danish |>
  filter (als == 1 & follow_up > 60) |>                                         # filtering cases with follow_up > 5 years
  filter(match != 159) |>                                                       # removing NEFL outlier
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
heterogeneity_base_quart_sensi_1_3 <- data.frame(explanatory = character(),
                                                 model = factor(),
                                                 p_value_heterogeneity = numeric(), 
                                                 stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula("surv_obj_sensi_1_3 ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart_sensi_1_3 <- rbind(heterogeneity_base_quart_sensi_1_3, 
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "sensi_1_3", 
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted 
heterogeneity_adjusted_quart_sensi_1_3 <- data.frame(explanatory = character(),
                                                     model = factor(),
                                                     p_value_heterogeneity = numeric(), 
                                                     stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula(paste("surv_obj_sensi_1_3 ~ ", paste(covariates, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart_sensi_1_3 <- rbind(heterogeneity_adjusted_quart_sensi_1_3, 
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "sensi_1_3", 
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests_sensi_1_3 <- 
  bind_rows(heterogeneity_base_quart_sensi_1_3, 
            heterogeneity_adjusted_quart_sensi_1_3) |>
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

rm(heterogeneity_base_quart_sensi_1_3, heterogeneity_adjusted_quart_sensi_1_3,
   trend_base_sensi_1_3, trend_adjusted_sensi_1_3)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)


### Base  ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish_sensi_1_3) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,  
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base_sensi_1_3 <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted  ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish_sensi_1_3) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var, 
    x_label = x_label)
}

cox_gam_results_adjusted_sensi_1_3 <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels)



# Sensi 1 + sensi 3 + sensi 4 - Removing NEFL outlier + filtering cases with follow- up < 5 years + filtering follow-up <= 50%----
bdd_cases_danish_sensi_1_3_4 <- bdd_danish |>
  filter(als == 1 & follow_up > 60) |>                                          # remove controls and remove cases with follow-up<5years 
  filter(match != 159) |>                                                       # remove NEFL outlier
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
#### base ----
heterogeneity_base_quart_sensi_1_3_4 <- data.frame(explanatory = character(),
                                                   model = factor(),
                                                   p_value_heterogeneity = numeric(), 
                                                   stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula("surv_obj_sensi_1_3_4 ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3_4) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_4)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart_sensi_1_3_4 <- rbind(heterogeneity_base_quart_sensi_1_3_4, 
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "sensi_1_3_4", 
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart_sensi_1_3_4 <- data.frame(explanatory = character(),
                                                       model = factor(),
                                                       p_value_heterogeneity = numeric(), 
                                                       stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula(paste("surv_obj_sensi_1_3_4 ~ ", paste(covariates, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3_4) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_4 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_4)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart_sensi_1_3_4 <- rbind(heterogeneity_adjusted_quart_sensi_1_3_4, 
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "sensi_1_3_4", 
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests_sensi_1_3_4 <- 
  bind_rows(heterogeneity_base_quart_sensi_1_3_4, 
            heterogeneity_adjusted_quart_sensi_1_3_4) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base ----
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

#### adjusted ----
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

rm(heterogeneity_base_quart_sensi_1_3_4, heterogeneity_adjusted_quart_sensi_1_3_4,
   trend_base_sensi_1_3_4, trend_adjusted_sensi_1_3_4)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)


### Base model ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish_sensi_1_3_4) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3_4, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,  
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base_sensi_1_3_4 <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted model ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish_sensi_1_3_4) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3_4, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var, 
    x_label = x_label)
}

cox_gam_results_adjusted_sensi_1_3_4 <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels)



# Sensi 1 + sensi 3 + sensi 5 - Removing NEFL outlier + filtering cases with follow-up < 5 years + filtering follow-up > 50%----
bdd_cases_danish_sensi_1_3_5 <- bdd_danish |>
  filter(als == 1 & follow_up > 60) |>                                          # remove controls and remove cases with follow-up<5years 
  filter(match != 159) |>                                                       # remove NEFL outlier
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
#### base ----
heterogeneity_base_quart_sensi_1_3_5 <- data.frame(explanatory = character(),
                                                   model = factor(),
                                                   p_value_heterogeneity = numeric(), 
                                                   stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula("surv_obj_sensi_1_3_5 ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3_5) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_5)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart_sensi_1_3_5 <- rbind(heterogeneity_base_quart_sensi_1_3_5, 
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "sensi_1_3_5", 
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart_sensi_1_3_5 <- data.frame(explanatory = character(),
                                                       model = factor(),
                                                       p_value_heterogeneity = numeric(), 
                                                       stringsAsFactors = FALSE)

for (expl in proteomic_quart) {
  
  formula_raw <- as.formula(paste("surv_obj_sensi_1_3_5 ~ ", paste(covariates, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_3_5) 
  
  formula <- as.formula(paste("surv_obj_sensi_1_3_5 ~", expl, "+ ", paste(covariates, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_3_5)  
  
  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart_sensi_1_3_5 <- rbind(heterogeneity_adjusted_quart_sensi_1_3_5, 
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "sensi_1_3_5", 
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests_sensi_1_3_5 <- 
  bind_rows(heterogeneity_base_quart_sensi_1_3_5, 
            heterogeneity_adjusted_quart_sensi_1_3_5) |>
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

rm(heterogeneity_base_quart_sensi_1_3_5, heterogeneity_adjusted_quart_sensi_1_3_5,
   trend_base_sensi_1_3_5, trend_adjusted_sensi_1_3_5)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)


### Base model ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish_sensi_1_3_5) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3_5, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,  
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base_sensi_1_3_5 <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted model ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish_sensi_1_3_5) {
  
  outcome <- with(bdd_cases_danish_sensi_1_3_5, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var, 
    x_label = x_label)
}

cox_gam_results_adjusted_sensi_1_3_5 <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels)



# Sensi 1 + sensi 7 - Removing the oulier for NEFL + filtering to females ----
bdd_cases_danish_sensi_1_7_female <- 
  bdd_danish |>
  filter(als == 1) |>                                                           # remove controls
  filter(match != 159) |>                                                       # remove NEFL outlier
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
#### base ----
heterogeneity_base_quart_sensi_1_7_female <- data.frame(explanatory = character(),
                                                        model = factor(),
                                                        p_value_heterogeneity = numeric(),
                                                        stringsAsFactors = FALSE)

for (expl in proteomic_quart) {

  formula_raw <- as.formula("surv_obj_sensi_1_7_female ~ diagnosis_age")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_7_female)

  formula <- as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ diagnosis_age"))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_female)

  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]

  heterogeneity_base_quart_sensi_1_7_female <- rbind(heterogeneity_base_quart_sensi_1_7_female,
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "sensi_1_7_female",
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart_sensi_1_7_female <- data.frame(explanatory = character(),
                                                            model = factor(),
                                                            p_value_heterogeneity = numeric(),
                                                            stringsAsFactors = FALSE)

for (expl in proteomic_quart) {

  formula_raw <- as.formula(paste("surv_obj_sensi_1_7_female ~ ", paste(covariates_sensi_7, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_7_female)

  formula <- as.formula(paste("surv_obj_sensi_1_7_female ~", expl, "+ ", paste(covariates_sensi_7, collapse = "+")))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_female)

  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]

  heterogeneity_adjusted_quart_sensi_1_7_female <- rbind(heterogeneity_adjusted_quart_sensi_1_7_female,
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "sensi_1_7_female",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests_sensi_1_7_female <-
  bind_rows(heterogeneity_base_quart_sensi_1_7_female,
            heterogeneity_adjusted_quart_sensi_1_7_female) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base ----
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

#### adjusted ----
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

rm(heterogeneity_base_quart_sensi_1_7_female, heterogeneity_adjusted_quart_sensi_1_7_female,
   trend_base_sensi_1_7_female, trend_adjusted_sensi_1_7_female)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic,
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_",
                                     ""),
                         proteomic)


### Base model ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish_sensi_1_7_female) {

  outcome <- with(bdd_cases_danish_sensi_1_7_female, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())

  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1)
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2)
  pval <- case_when(pval < 0.01 ~ "< 0.01",
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]]

  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)

  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base_sensi_1_7_female <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted model ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish_sensi_1_7_female) {

  outcome <- with(bdd_cases_danish_sensi_1_7_female, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())

  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1)
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2)
  pval <- case_when(pval < 0.01 ~ "< 0.01",
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]]

  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)

  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_adjusted_sensi_1_7_female <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels)

# Sensi 1 + sensi 7 - Removing the oulier for NEFL + filtering to male ----
bdd_cases_danish_sensi_1_7_male <- bdd_danish |>
  filter(als == 1) |>                                                           # remove controls
  filter(match != 159) |>                                                       # remove NEFL outlier
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
#### base ----
heterogeneity_base_quart_sensi_1_7_male <- data.frame(explanatory = character(),
                                                      model = factor(),
                                                      p_value_heterogeneity = numeric(),
                                                      stringsAsFactors = FALSE)

for (expl in proteomic_quart) {

  formula_raw <- as.formula("surv_obj_sensi_1_7_male ~ diagnosis_age")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_7_male)

  formula <- as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ diagnosis_age"))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_male)

  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]

  heterogeneity_base_quart_sensi_1_7_male <- rbind(heterogeneity_base_quart_sensi_1_7_male,
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               analysis = "sensi_1_7_male",
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart_sensi_1_7_male <- data.frame(explanatory = character(),
                                                          model = factor(),
                                                          p_value_heterogeneity = numeric(),
                                                          stringsAsFactors = FALSE)

for (expl in proteomic_quart) {

  formula_raw <- as.formula(paste("surv_obj_sensi_1_7_male ~ ", paste(covariates_sensi_7, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish_sensi_1_7_male)

  formula <- as.formula(paste("surv_obj_sensi_1_7_male ~", expl, "+ ", paste(covariates_sensi_7, collapse = "+")))
  model <- coxph(formula, data = bdd_cases_danish_sensi_1_7_male)

  anova <- anova(model_raw, model, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]

  heterogeneity_adjusted_quart_sensi_1_7_male <- rbind(heterogeneity_adjusted_quart_sensi_1_7_male,
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   analysis = "sensi_1_7_male",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p_value_heterogeneity)

heterogeneity_tests_sensi_1_7_male <-
  bind_rows(heterogeneity_base_quart_sensi_1_7_male,
            heterogeneity_adjusted_quart_sensi_1_7_male) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### Trend tests ----
#### base ----
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

#### adjusted ----
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

rm(heterogeneity_base_quart_sensi_1_7_male, heterogeneity_adjusted_quart_sensi_1_7_male,
   trend_base_sensi_1_7_male, trend_adjusted_sensi_1_7_male)


## Cox model (GAMs) ----
vars_labels <- set_names(str_replace(proteomic,
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_",
                                     ""),
                         proteomic)


### Base model ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish_sensi_1_7_male) {

  outcome <- with(bdd_cases_danish_sensi_1_7_male, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())

  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1)
  pval_raw <- smry$s.table[1, "p-value"]
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2)
  pval <- case_when(pval < 0.01 ~ "< 0.01",
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]]

  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)

  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval_raw = pval_raw,
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base_sensi_1_7_male <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)



### Adjusted model ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish_sensi_1_7_male) {

  outcome <- with(bdd_cases_danish_sensi_1_7_male, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())

  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1)
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2)
  pval <- case_when(pval < 0.01 ~ "< 0.01",
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]]

  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)

  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_adjusted_sensi_1_7_male <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted, vars_labels, covariates_sensi_7)

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


extra_rows <- 
  main_results |>
  filter(model == "base" & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "HR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '')

proteomic_quart_ALS_table <- 
  main_results |>
  filter(model == "base" & term != "Continuous" & analysis == "main") |>            # select quartile results
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, HR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base")) 

proteomic_quart_ALS_table <- 
  proteomic_quart_ALS_table |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base") |>
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
    "Heterogeneity test" = "Base models",  "Trend test" = "Base models") |>
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



## Figure proteomic - als survival - base gam (main) ----
# pvals <- sapply(model1_gam, function(m) m$s.table[1, "p-value"])
# signif_vars <- names(pvals)[pvals < 0.05]
# signif_vars_labels <- set_names(str_replace(signif_vars, 
#                                             "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
#                                             ""), 
#                                 signif_vars)
# 
# plot_base_gam <- map(signif_vars, function(var) {
#   
#   outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
#   formula <- as.formula(paste("outcome ~ s(", var, ") + diagnosis_age + sex")) 
#   
#   model <- gam(formula,                                                         # run the cox-gam model
#                family = cox.ph(), 
#                method = "ML", 
#                data = bdd_cases_danish)            
#   
#   bdd_pred <- bdd_cases_danish |>                                               # création bdd avec protein + covariables ramenées à leur moyenne
#     mutate(
#       adj_diagnosis_age = mean(diagnosis_age, na.rm = TRUE),
#       adj_sex = names(which.max(table(sex)))) |>
#     select(all_of(var), starts_with("adj_")) |>
#     rename_with(~ gsub("adj_", "", .x)) 
#   
#   pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
#   
#   bdd_pred <- bdd_pred |>
#     mutate(
#       hazard_ratio = exp(pred$fit),
#       hazard_lower = exp(pred$fit - 1.96 * pred$se.fit),
#       hazard_upper = exp(pred$fit + 1.96 * pred$se.fit))
#   
#   model_summary <- summary(model)
#   edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
#   p_value <- model_summary$s.table[1, "p-value"]
#   p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
#   x_max <- max(bdd_cases_danish[[var]], na.rm = TRUE)
#   x_label <- signif_vars_labels[var] 
#   
#   p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
#     geom_line(color = "blue", size = 1) +
#     geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
#     labs(x = var, y = "Hazard Ratio (HR)") +
#     annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
#              hjust = 1, vjust = 1.2, size = 4, color = "black") +
#     theme_minimal() + 
#     scale_y_log10(limits = c(10, 15000)) +
#     theme(axis.text.x = element_text(color = 'white'),
#           axis.title.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.ticks.x = element_blank()) +
#     ggtitle("Base model")
#   # distribution of the non-scaled protein variables 
#   p2 <- bdd_cases_danish |>
#     ggplot() +
#     aes(x = "", y = .data[[var]]) +
#     geom_boxplot(fill = "blue") +
#     coord_flip() +
#     ylab(x_label) + 
#     xlab("") + 
#     theme_minimal()
#   
#   p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
#                            guides = 'collect') + 
#     theme_minimal()
#   p
# }) |>
#   set_names(signif_vars_labels)
# rm(pvals, signif_vars, signif_vars_labels)

vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish <- map(cox_gam_results_base, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base)

## Figure proteomic - als survival - adjusted gam (main) ----
# pvals <- sapply(model2_gam, function(m) m$s.table[1, "p-value"])
# signif_vars <- names(pvals)[pvals < 0.05]
# signif_vars_labels <- set_names(str_replace(signif_vars, 
#                                             "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
#                                             ""), 
#                                 signif_vars)
# 
# plot_adjusted_gam <- map(signif_vars, function(var) {
#   
#   outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
#   formula <- as.formula(paste("outcome ~ s(", var, ") + ", paste(covariates, collapse = "+"))) 
#   
#   model <- gam(formula,                                                         # run the cox-gam model
#                family = cox.ph(), 
#                method = "ML", 
#                data = bdd_cases_danish)     
#   
#   bdd_pred <- bdd_cases_danish |>                                               # création bdd avec protein + covariables ramenées à leur moyenne
#     mutate(
#       adj_diagnosis_age = mean(diagnosis_age, na.rm = TRUE),
#       adj_sex = names(which.max(table(sex))), 
#       adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i))), 
#       adj_bmi = mean(bmi, na.rm = TRUE)) |>
#     select(all_of(var), starts_with("adj_")) |>
#     rename_with(~ gsub("adj_", "", .x)) 
#   
#   pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
#   
#   bdd_pred <- bdd_pred |>
#     mutate(
#       hazard_ratio = exp(pred$fit),
#       hazard_lower = exp(pred$fit - 1.96 * pred$se.fit),
#       hazard_upper = exp(pred$fit + 1.96 * pred$se.fit))
#   
#   model_summary <- summary(model)
#   edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
#   p_value <- model_summary$s.table[1, "p-value"]
#   p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
#   x_max <- max(bdd_cases_danish[[var]], na.rm = TRUE)
#   x_label <- signif_vars_labels[var] 
#   
#   p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
#     geom_line(color = "blue", size = 1) +
#     geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
#     labs(x = var, y = "Hazard Ratio (HR)") +
#     annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
#              hjust = 1, vjust = 1.2, size = 4, color = "black") +
#     theme_minimal() + 
#     scale_y_log10(limits = c(100, 100000),
#                   labels = scales::label_number(accuracy = 1)) +
#     theme(axis.text.x = element_text(color = 'white'),
#           axis.title.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.ticks.x = element_blank()) +
#     ggtitle("Adjusted model")
#   # distribution of non-scaled protein variables 
#   p2 <- bdd_cases_danish |>
#     ggplot() +
#     aes(x = "", y = .data[[var]]) +
#     geom_boxplot(fill = "blue") +
#     coord_flip() +
#     ylab(x_label) + 
#     xlab("") + 
#     theme_minimal()
#   
#   p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
#                            guides = 'collect') + 
#     theme_minimal()
#   p
# }) |>
#   set_names(signif_vars_labels)
# rm(pvals, signif_vars, signif_vars_labels)


all_fits_adjusted <- map_dfr(cox_gam_results_adjusted, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish <- map(cox_gam_results_adjusted, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted)





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



## Figure proteomic - als survival - base gam (sensi_1_3) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base_sensi_1_3, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish_sensi_1_3 <- map(cox_gam_results_base_sensi_1_3, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base_sensi_1_3)

## Figure proteomic - als survival - adjusted gam (sensi_1_3) ----
all_fits_adjusted <- map_dfr(cox_gam_results_adjusted_sensi_1_3, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish_sensi_1_3 <- map(cox_gam_results_adjusted_sensi_1_3, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted_sensi_1_3)




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



## Figure proteomic - als survival - base gam (sensi_1_3_4) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base_sensi_1_3_4, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish_sensi_1_3_4 <- map(cox_gam_results_base_sensi_1_3_4, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3_4, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base_sensi_1_3_4)

## Figure proteomic - als survival - adjusted gam (sensi_1_3_4) ----
all_fits_adjusted <- map_dfr(cox_gam_results_adjusted_sensi_1_3_4, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish_sensi_1_3_4 <- map(cox_gam_results_adjusted_sensi_1_3_4, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3_4, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted_sensi_1_3_4)


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



## Figure proteomic - als survival - base gam (sensi_1_3_5) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base_sensi_1_3_5, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish_sensi_1_3_5 <- map(cox_gam_results_base_sensi_1_3_5, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3_5, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base_sensi_1_3_5)

## Figure proteomic - als survival - adjusted gam (sensi_1_3_5) ----
all_fits_adjusted <- map_dfr(cox_gam_results_adjusted_sensi_1_3_5, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish_sensi_1_3_5 <- map(cox_gam_results_adjusted_sensi_1_3_5, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_3_5, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted_sensi_1_3_5)



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



## Figure proteomic - als survival - base gam (sensi_1_7_female) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base_sensi_1_7_female, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish_sensi_1_7_female <- map(cox_gam_results_base_sensi_1_7_female, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_7_female, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base_sensi_1_7_female)

## Figure proteomic - als survival - adjusted gam (sensi_1_7_female) ----
all_fits_adjusted <- map_dfr(cox_gam_results_adjusted_sensi_1_7_female, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish_sensi_1_7_female <- map(cox_gam_results_adjusted_sensi_1_7_female, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_7_female, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted_sensi_1_7_female)



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



## Figure proteomic - als survival - base gam (sensi_1_7_male) ----
vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

all_fits_base <- map_dfr(cox_gam_results_base_sensi_1_7_male, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish_sensi_1_7_male <- map(cox_gam_results_base_sensi_1_7_male, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_7_male, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base_sensi_1_7_male)

## Figure proteomic - als survival - adjusted gam (sensi_1_7_male) ----
all_fits_adjusted <- map_dfr(cox_gam_results_adjusted_sensi_1_7_male, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish_sensi_1_7_male <- map(cox_gam_results_adjusted_sensi_1_7_male, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish_sensi_1_7_male, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)


rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted_sensi_1_7_male)





# Assemblage ----
results_proteomic_ALS_survival <- 
  list(
    main = list(covar = covar, 
                main_results= main_results, 
                
                proteomic_sd_ALS_table = proteomic_sd_ALS_table,
                proteomic_quart_ALS_table = proteomic_quart_ALS_table,
                
                proteomic_sd_ALS_base_figure = proteomic_sd_ALS_base_figure, 
                proteomic_sd_ALS_adjusted_figure = proteomic_sd_ALS_adjusted_figure, 
                
                plot_base_cox_gam_danish = plot_base_cox_gam_danish, 
                plot_adjusted_cox_gam_danish = plot_adjusted_cox_gam_danish), 
    
    sensi_1_3 = list(proteomic_sd_ALS_table_sensi_1_3 = proteomic_sd_ALS_table_sensi_1_3,
                     proteomic_quart_ALS_table_sensi_1_3 = proteomic_quart_ALS_table_sensi_1_3,
                     
                     proteomic_sd_ALS_base_figure_sensi_1_3 = proteomic_sd_ALS_base_figure_sensi_1_3, 
                     proteomic_sd_ALS_adjusted_figure_sensi_1_3 = proteomic_sd_ALS_adjusted_figure_sensi_1_3,
                     
                     plot_base_cox_gam_danish_sensi_1_3 = plot_base_cox_gam_danish_sensi_1_3, 
                     plot_adjusted_cox_gam_danish_sensi_1_3 = plot_adjusted_cox_gam_danish_sensi_1_3), 
    
    sensi_1_3_4 = list(proteomic_sd_ALS_table_sensi_1_3_4 = proteomic_sd_ALS_table_sensi_1_3_4,
                     proteomic_quart_ALS_table_sensi_1_3_4 = proteomic_quart_ALS_table_sensi_1_3_4,
                     
                     proteomic_sd_ALS_base_figure_sensi_1_3_4 = proteomic_sd_ALS_base_figure_sensi_1_3_4, 
                     proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_4,
                     
                     plot_base_cox_gam_danish_sensi_1_3_4 = plot_base_cox_gam_danish_sensi_1_3_4, 
                     plot_adjusted_cox_gam_danish_sensi_1_3_4 = plot_adjusted_cox_gam_danish_sensi_1_3_4), 
    
    sensi_1_3_5 = list(proteomic_sd_ALS_table_sensi_1_3_5 = proteomic_sd_ALS_table_sensi_1_3_5,
                       proteomic_quart_ALS_table_sensi_1_3_5 = proteomic_quart_ALS_table_sensi_1_3_5,
                       
                       proteomic_sd_ALS_base_figure_sensi_1_3_5 = proteomic_sd_ALS_base_figure_sensi_1_3_5, 
                       proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_5,
                       
                       plot_base_cox_gam_danish_sensi_1_3_5 = plot_base_cox_gam_danish_sensi_1_3_5, 
                       plot_adjusted_cox_gam_danish_sensi_1_3_5 = plot_adjusted_cox_gam_danish_sensi_1_3_5), 
    sensi_1_7_female = list(proteomic_sd_ALS_table_sensi_1_7_female = proteomic_sd_ALS_table_sensi_1_7_female,
                            proteomic_quart_ALS_table_sensi_1_7_female = proteomic_quart_ALS_table_sensi_1_7_female,
                            
                            proteomic_sd_ALS_base_figure_sensi_1_7_female = proteomic_sd_ALS_base_figure_sensi_1_7_female, 
                            proteomic_sd_ALS_adjusted_figure_sensi_1_7_female = proteomic_sd_ALS_adjusted_figure_sensi_1_7_female,
                            
                            plot_base_cox_gam_danish_sensi_1_7_female = plot_base_cox_gam_danish_sensi_1_7_female, 
                            plot_adjusted_cox_gam_danish_sensi_1_7_female = plot_adjusted_cox_gam_danish_sensi_1_7_female), 
    sensi_1_7_male = list(proteomic_sd_ALS_table_sensi_1_7_male = proteomic_sd_ALS_table_sensi_1_7_male,
                          proteomic_quart_ALS_table_sensi_1_7_male = proteomic_quart_ALS_table_sensi_1_7_male,
                          
                          proteomic_sd_ALS_base_figure_sensi_1_7_male = proteomic_sd_ALS_base_figure_sensi_1_7_male, 
                          proteomic_sd_ALS_adjusted_figure_sensi_1_7_male = proteomic_sd_ALS_adjusted_figure_sensi_1_7_male,
                          
                          plot_base_cox_gam_danish_sensi_1_7_male = plot_base_cox_gam_danish_sensi_1_7_male, 
                          plot_adjusted_cox_gam_danish_sensi_1_7_male = plot_adjusted_cox_gam_danish_sensi_1_7_male))

rm(bdd_cases_danish, 
   bdd_cases_danish_sensi_1_3, 
   bdd_cases_danish_sensi_1_3_4, 
   bdd_cases_danish_sensi_1_3_5, 
   bdd_cases_danish_sensi_1_7_female,
   bdd_cases_danish_sensi_1_7_male,
   surv_obj, 
   surv_obj_sensi_1_3, 
   surv_obj_sensi_1_3_4, 
   surv_obj_sensi_1_3_5, 
   surv_obj_sensi_1_7_female, 
   surv_obj_sensi_1_7_male, 
   covariates, 
   covar, 
   main_results, 
   proteomic_sd_ALS_table,
   proteomic_quart_ALS_table, 
   proteomic_sd_ALS_base_figure, 
   proteomic_sd_ALS_adjusted_figure, 
   plot_base_cox_gam_danish, 
   plot_adjusted_cox_gam_danish, 
   
   proteomic_sd_ALS_table_sensi_1_3,
   proteomic_quart_ALS_table_sensi_1_3,
   proteomic_sd_ALS_base_figure_sensi_1_3, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3,
   plot_base_cox_gam_danish_sensi_1_3, 
   plot_adjusted_cox_gam_danish_sensi_1_3, 
   
   proteomic_sd_ALS_table_sensi_1_3_4,
   proteomic_quart_ALS_table_sensi_1_3_4,
   proteomic_sd_ALS_base_figure_sensi_1_3_4, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_4,
   plot_base_cox_gam_danish_sensi_1_3_4, 
   plot_adjusted_cox_gam_danish_sensi_1_3_4,
   
   proteomic_sd_ALS_table_sensi_1_3_5,
   proteomic_quart_ALS_table_sensi_1_3_5,
   proteomic_sd_ALS_base_figure_sensi_1_3_5,
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_5,
   plot_base_cox_gam_danish_sensi_1_3_5,
   plot_adjusted_cox_gam_danish_sensi_1_3_5, 
   
   proteomic_sd_ALS_table_sensi_1_7_female,
   proteomic_quart_ALS_table_sensi_1_7_female,
   proteomic_sd_ALS_base_figure_sensi_1_7_female,
   proteomic_sd_ALS_adjusted_figure_sensi_1_7_female,
   plot_base_cox_gam_danish_sensi_1_7_female,
   plot_adjusted_cox_gam_danish_sensi_1_7_female, 
   
   proteomic_sd_ALS_table_sensi_1_7_male,
   proteomic_quart_ALS_table_sensi_1_7_male,
   proteomic_sd_ALS_base_figure_sensi_1_7_male,
   proteomic_sd_ALS_adjusted_figure_sensi_1_7_male,
   plot_base_cox_gam_danish_sensi_1_7_male,
   plot_adjusted_cox_gam_danish_sensi_1_7_male)



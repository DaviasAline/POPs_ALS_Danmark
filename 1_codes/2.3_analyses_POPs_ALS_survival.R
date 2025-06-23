# Aline Davias
# April 29, 2025 
# Analysis of survival after ALS diagnosis depending on POPs levels  

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")

# Creation of cases specific datasets ----
bdd_cases_danish <- bdd_danish |>
  filter (als == 1) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_group)) |>
  mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group),                                             # create cohort and cases specific scaled POPs variables 
                  ~as.numeric(scale(.x)),
                  .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(ΣPBDE, ΣPBDE_quart) 

bdd_cases_FMC <- bdd_finnish |>
  filter (als == 1) |>
  filter(study == "FMC") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol,   
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", 
         OCP_γ_HCH = "OCP_γ_HCH_raw", "Σchlordane", OCP_PeCB = "OCP_PeCB_raw") |>
  mutate(across(all_of(POPs_group_finnish), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_finnish),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) 

bdd_cases_FMCF <- bdd_finnish |>
  filter (als == 1) |>
  filter(study == "FMCF") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol,   
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", 
         OCP_γ_HCH = "OCP_γ_HCH_raw", "Σchlordane", OCP_PeCB = "OCP_PeCB_raw") |>
  mutate(across(all_of(POPs_group_finnish), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                     labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_finnish),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) 

bdd_cases_MFH <- bdd_finnish |>
  filter (als == 1) |>
  filter(study == "MFH") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol,   
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", 
         OCP_γ_HCH = "OCP_γ_HCH_raw", "Σchlordane", OCP_PeCB = "OCP_PeCB_raw") |>
  mutate(across(all_of(POPs_group_finnish), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                     labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_finnish),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) 

surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                        event = bdd_cases_danish$status_death)
surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,               
                        event = bdd_cases_FMC$status_death)
surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death,              
                        event = bdd_cases_FMCF$status_death)
surv_obj_MFH <- Surv(time = bdd_cases_MFH$follow_up_death,               
                        event = bdd_cases_MFH$status_death)

# Danish cohort ----
## Covar model ----
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

## Main analysis 
## Cox model (sd) ----
### Base ----
model1_cox_sd_danish <- map_dfr(POPs_group_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_danish ~", expl, "+ baseline_age + sex"))        # set the formulas                
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()          # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "base", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_starts(term, explanatory))                                         # remove the covariates results 
  
})

### Adjusted ----
model2_cox_sd_danish <- map_dfr(POPs_group_sd, function(expl) {

  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+",  
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

### Copollutant ----
POPs_group_sd_bis <- setdiff(POPs_group_sd, "PCB_4_sd")                         # remove the 4 most abundant PCB because they are already NDL-PCB
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_sd_bis)

formula_danish <-                                                               # set the formulas  
  as.formula(paste("surv_obj_danish ~",   
                   paste(c(POPs_group_sd_bis, covariates_danish), collapse = " + ")))

model_summary <- 
  coxph(formula_danish, data = bdd_cases_danish) |> summary() 
coefs <- model_summary$coefficients
model3_cox_sd_danish <- tibble(                                                 # creation of a table of results
  study = "Danish", 
  model = "copollutant", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))
rm(POPs_group_sd_bis, pollutant_labels_bis, formula_danish, model_summary, coefs)

## Cox model (quart) ----
### Base ----
model1_cox_quart_danish <- map_dfr(POPs_group_quart, function(expl) {
  
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
    filter(str_starts(term, explanatory))                                      # remove the covariates results
})

### Adjusted ----
model2_cox_quart_danish <- map_dfr(POPs_group_quart, function(expl) {
  
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


### Copollutant ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL <- 
  gam(outcome ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_danish) |>
  summary()
model3_quart_PCB_DL <- model3_quart_PCB_DL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_PCB_NDL <- 
  gam(outcome ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_PCB_NDL <- model3_quart_PCB_NDL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_HCB <- 
  gam(outcome ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_HCB <- model3_quart_HCB$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣDDT <- 
  gam(outcome ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_ΣDDT <- model3_quart_ΣDDT$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_β_HCH <- 
  gam(outcome ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_β_HCH <- model3_quart_β_HCH$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_Σchlordane <- 
  gam(outcome ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_Σchlordane <- model3_quart_Σchlordane$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣPBDE <- 
  gam(outcome ~ 
        ΣPBDE_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |>
  summary()
model3_quart_ΣPBDE <- model3_quart_ΣPBDE$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_cox_quart_danish <- bind_rows(
  model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane, model3_quart_ΣPBDE) |>
  filter(grepl("quart", variable)) |>
  mutate(
    study = "Danish", 
    model = "copollutant", 
    term = variable,
    explanatory = gsub("Q2", "", variable), 
    explanatory = gsub("Q3", "", explanatory), 
    explanatory = gsub("Q4", "", explanatory), 
    coef = Estimate, 
    se = `Std. Error`, 
    `p-value` =`Pr(>|z|)`) |> 
  select(study, model, term, explanatory, coef, se, `p-value`)

rm(model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane, model3_quart_ΣPBDE, 
   outcome)

## Cox-gam model (sd) ----
### Base  ----
plot_base_cox_gam_danish <- map(POPs_group_sd, function(var) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula <- as.formula(paste("outcome ~ s(", var, ") + baseline_age + sex")) 
  
  model <- gam(formula,                                                         # run the cox-gam model
               family = cox.ph(), 
               method = "ML", 
               data = bdd_cases_danish)            
  
  bdd_pred <- bdd_cases_danish |>                                               # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      hazard_ratio = exp(pred$fit),
      hazard_lower = exp(pred$fit - 1.96 * pred$se.fit),
      hazard_upper = exp(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_cases_danish[[var]], na.rm = TRUE)
  x_label <- POPs_group_sd_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Hazard Ratio (HR)") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() + 
    scale_y_log10(limits = c(0.1, 110)) +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Base model")
  
  var_orig <- gsub("_sd$", "", var)                                             # distribution of the non-scaled POP variables 
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var_orig]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |> 
  set_names(POPs_group)


### Adjusted ----
plot_adjusted_cox_gam_danish <- map(POPs_group_sd, function(var) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula <- as.formula(paste("outcome ~ s(", var, ") + baseline_age + sex +  smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i")) 
  
  model <- gam(formula,                                                         # run the cox-gam model
               family = cox.ph(), 
               method = "ML", 
               data = bdd_cases_danish)     
  
  bdd_pred <- bdd_cases_danish |>                                               # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i))), 
      adj_bmi = mean(bmi, na.rm = TRUE),  
      adj_cholesterol_i = mean(cholesterol_i, na.rm = TRUE), 
      adj_marital_status_2cat_i = names(which.max(table(marital_status_2cat_i))), 
      adj_education_i = names(which.max(table(education_i)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      hazard_ratio = exp(pred$fit),
      hazard_lower = exp(pred$fit - 1.96 * pred$se.fit),
      hazard_upper = exp(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_cases_danish[[var]], na.rm = TRUE)
  x_label <- POPs_group_sd_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Hazard Ratio (HR)") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() + 
    scale_y_log10(limits = c(0.1, 2000), 
                  labels = scales::label_number(accuracy = 1)) +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
  
  var_orig <- gsub("_sd$", "", var)                                             # distribution of non-scaled POP variables 
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var_orig]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |> 
  set_names(POPs_group)

### Copollutant ----
POPs_group_sd_bis <- setdiff(POPs_group_sd, "PCB_4_sd")                         # PCB4 not included because already present in NDL-PCBs
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
POPs_group_sd_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_sd_bis)

outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model <- gam(outcome ~ s(PCB_DL_sd) + s(PCB_NDL_sd) + s(OCP_HCB_sd) + s(ΣDDT_sd) + 
               s(OCP_β_HCH_sd) + s(Σchlordane_sd) + s(ΣPBDE_sd) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = cox.ph(), 
             method = "ML", 
             data = bdd_cases_danish)


plot_copollutant_cox_gam_danish <- map(POPs_group_sd_bis, function(var) {
  
  bdd_pred <- bdd_cases_danish |>
    mutate(across(all_of(covariates_danish),                                    # fixe toutes les covariables à leurs moyennes
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))))) |>
    mutate(across(setdiff(POPs_group_sd_bis, var), 
                  ~ mean(., na.rm = TRUE)))|>                                   # Fixe tous les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates_danish), all_of(setdiff(POPs_group_sd_bis, var)))  
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      hazard_ratio = exp(pred$fit),
      hazard_lower = exp(pred$fit - 1.96 * pred$se.fit),
      hazard_upper = exp(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_cases_danish[[var]], na.rm = TRUE)
  x_label <- POPs_group_sd_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Hazard Ratio (HR)") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_log10(limits = c(0.1, 2000), 
                  labels = scales::label_number(accuracy = 1)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Co-pollutant model")
  
  var_orig <- gsub("_sd$", "", var)          
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var_orig]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
}) |> set_names(POPs_group_bis)
rm(POPs_group_sd_bis, POPs_group_sd_labels_bis, model, outcome, POPs_group_bis)

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
### crude 
covar_crude_finnish <- map_dfr(c("baseline_age", "sex", "thawed", "level_urbanization", 
                           covariates_finnish), function(expl) {
  
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

### adjusted 
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

covar_finnish <- bind_rows(covar_crude_finnish, covar_adjusted_finnish)
rm(formula, covar_crude_finnish, covar_adjusted_finnish)

## Cox model (sd) ----
### Base ----
model1_cox_sd_finnish <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_FMC <-                                                                # set the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox models 
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") |>
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

### Adjusted ----
model2_cox_sd_finnish <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                                  paste(covariates_finnish, collapse = " + "), 
                                  "+ baseline_age + sex + thawed + level_urbanization"))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                   paste(covariates_finnish, collapse = " + "), 
                                   "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox models 
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") |>
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

### Copollutant ----
POPs_group_sd_finnish_bis <- setdiff(POPs_group_sd_finnish, c("OCP_β_HCH_sd",  "OCP_γ_HCH_sd", "PCB_4_sd"))

formula_FMC <-                                                                  # creation of the formulas
  as.formula(paste("surv_obj_FMC ~", paste(c(POPs_group_sd_finnish_bis, covariates_finnish, "baseline_age", "sex", "thawed", "level_urbanization"), collapse = "+")))  
formula_FMCF <- 
  as.formula(paste("surv_obj_FMCF ~", paste(c(POPs_group_sd_finnish_bis, covariates_finnish, "baseline_age", "sex", "thawed", "level_urbanization"), collapse = "+")))  

results <- list(                                                                # run of the simple cox models 
  finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2, 
         explanatory = term) |>
  filter(str_detect(term, "_sd"))                                               # remove the covariates results

model3_cox_sd_finnish <- results |>
  group_by(term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      study = "Finnish", 
      model = "copollutant", 
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval))
  }) |> 
  ungroup() |>
  mutate(explanatory = term)

rm(POPs_group_sd_finnish_bis, formula_FMC, formula_FMCF, results)

## Cox model (quart) ----
### Base ----
model1_cox_quart_finnish <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_FMC <-                                                                # creation of the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox model
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") |>
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

### Adjusted ----
model2_cox_quart_finnish <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                                  paste(covariates_finnish, collapse = " + "), 
                                  "+ baseline_age + sex + thawed + level_urbanization"))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                   paste(covariates_finnish, collapse = " + "), 
                                   "+ baseline_age + sex + thawed + level_urbanization"))
  
  results <- list(                                                              # run of the simple cox model
    finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
    finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
    bind_rows(.id = "dataset") |>
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

### Copollutant ---- 
# POPs_group_quart_finnish_bis <- setdiff(POPs_group_quart_finnish, c("OCP_β_HCH_quart",  "OCP_γ_HCH_quart", "PCB_4_quart"))
# 
# model3_quart_PCB_DL_FMC <-
#   gam(surv_obj_FMC ~
#         PCB_DL_quart +
#         s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',                                                            # maximum likelihood
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_PCB_DL_FMC <- model3_quart_PCB_DL_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_PCB_NDL_FMC <-
#   gam(surv_obj_FMC ~
#         PCB_NDL_quart +
#         s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_PCB_NDL_FMC <- model3_quart_PCB_NDL_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_HCB_FMC <-
#   gam(surv_obj_FMC ~
#         OCP_HCB_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_HCB_FMC <- model3_quart_HCB_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_ΣDDT_FMC <-
#   gam(surv_obj_FMC ~
#         ΣDDT_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_ΣDDT_FMC <- model3_quart_ΣDDT_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_ΣHCH_FMC <-
#   gam(surv_obj_FMC ~
#         ΣHCH_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_ΣHCH_FMC <- model3_quart_ΣHCH_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_Σchlordane_FMC <-
#   gam(surv_obj_FMC ~
#         Σchlordane_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_Σchlordane_FMC <- model3_quart_Σchlordane_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_PeCB_FMC <-
#   gam(surv_obj_FMC ~
#         OCP_PeCB_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMC) |>
#   summary()
# model3_quart_PeCB_FMC <- model3_quart_PeCB_FMC$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_PCB_DL_FMCF <-
#   gam(surv_obj_FMCF ~
#         PCB_DL_quart +
#         s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',                                                            # maximum likelihood
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_PCB_DL_FMCF <- model3_quart_PCB_DL_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_PCB_NDL_FMCF <-
#   gam(surv_obj_FMCF ~
#         PCB_NDL_quart +
#         s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_PCB_NDL_FMCF <- model3_quart_PCB_NDL_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_HCB_FMCF <-
#   gam(surv_obj_FMCF ~
#         OCP_HCB_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_HCB_FMCF <- model3_quart_HCB_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_ΣDDT_FMCF <-
#   gam(surv_obj_FMCF ~
#         ΣDDT_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_ΣDDT_FMCF <- model3_quart_ΣDDT_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_ΣHCH_FMCF <-
#   gam(surv_obj_FMCF ~
#         ΣHCH_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_ΣHCH_FMCF <- model3_quart_ΣHCH_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_Σchlordane_FMCF <-
#   gam(surv_obj_FMCF ~
#         Σchlordane_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_Σchlordane_FMCF <- model3_quart_Σchlordane_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_quart_PeCB_FMCF <-
#   gam(surv_obj_FMCF ~
#         OCP_PeCB_quart +
#         s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
#         baseline_age + sex + thawed + level_urbanization + marital_status_2cat + smoking_2cat + bmi + cholesterol,
#       family = cox.ph(),
#       method = 'ML',
#       data = bdd_cases_FMCF) |>
#   summary()
# model3_quart_PeCB_FMCF <- model3_quart_PeCB_FMCF$p.table |>
#   as.data.frame() |>
#   rownames_to_column("variable")
# 
# model3_cox_quart_FMC <- bind_rows(
#   model3_quart_PCB_DL_FMC, model3_quart_PCB_NDL_FMC, model3_quart_HCB_FMC, model3_quart_ΣDDT_FMC,
#   model3_quart_ΣHCH_FMC, model3_quart_Σchlordane_FMC, model3_quart_PeCB_FMC) |>
#   filter(grepl("quart", variable)) |>
#   mutate(
#     study = "FMC",
#     model = "copollutant",
#     term = variable,
#     explanatory = gsub("Q2", "", variable),
#     explanatory = gsub("Q3", "", explanatory),
#     explanatory = gsub("Q4", "", explanatory),
#     coef = Estimate,
#     se = `Std. Error`,
#     `p-value` =`Pr(>|z|)`) |>
#   select(study, model, term, explanatory, coef, se, `p-value`)
# 
# rm(model3_quart_PCB_DL_FMC, model3_quart_PCB_NDL_FMC, model3_quart_HCB_FMC, model3_quart_ΣDDT_FMC,
#    model3_quart_ΣHCH_FMC, model3_quart_Σchlordane_FMC, model3_quart_PeCB_FMC,
#    outcome)
# 
# results <- list(                                                                # run of the simple cox models
#   finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
#   finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
#   bind_rows(.id = "dataset") |>
#   mutate(var = se^2) |>
#   filter(str_detect(term, "_quart"))                                            # remove the covariates results
# 
# model3_cox_quart_finnish <- results |>
#   group_by(term) |>
#   group_modify(~ {
#     rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
#     tibble(                                                                     # results table creation
#       study = "Finnish",
#       model = "copollutant",
#       HR = exp(as.numeric(rma_fit$beta)),
#       lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
#       upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
#       `p-value` = as.numeric(rma_fit$pval))
#   }) |>
#   ungroup() |>
#   mutate(explanatory = gsub("Q2", "", term),
#          explanatory = gsub("Q3", "", explanatory),
#          explanatory = gsub("Q4", "", explanatory))
# 
# rm(POPs_group_quart_bis, formula_FMC, formula_FMCF, results)
# 

# Metaanalysis ----
## Cox model (sd) ----
### Base ----
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
    bind_rows(.id = "dataset") |>
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


### Adjusted ----
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
    bind_rows(.id = "dataset") |>
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

## Cox model (quart) ----
### Base ----
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
    bind_rows(.id = "dataset") |>
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

### Adjusted  ----
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
    bind_rows(.id = "dataset") |>
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
    model1_cox_sd_danish, model2_cox_sd_danish, model3_cox_sd_danish,
    model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  
  bind_rows(model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish,
            model1_cox_quart_finnish, model2_cox_quart_finnish, #model3_cox_quart_finnish,
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

rm(model1_cox_sd_danish, model2_cox_sd_danish, model3_cox_sd_danish, 
   model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish,
   model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish, 
   model1_cox_quart_finnish, model2_cox_quart_finnish, #model3_cox_quart_finnish,  
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
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI, serum total cholesterol and marital status. 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", 
    " HR " = "Copollutant Model", " 95% CI " = "Copollutant Model", " p-value " = "Copollutant Model") |>
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
  select(explanatory, term, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", 
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant") |>
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
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", 
    " HR " = "Copollutant Model", " 95% CI " = "Copollutant Model", " p-value " = "Copollutant Model") |>
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
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", 
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age, sex, level of urbanization and serum freeze-thaw cycles. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education (except for FMC study cohort). 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", 
    " HR " = "Copollutant Model", " 95% CI " = "Copollutant Model", " p-value " = "Copollutant Model") |>
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
  select(explanatory, term, contains("base"), contains("adjusted")
         #, contains("copollutant")
         ) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted"
         #, " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant"
         ) |>
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
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model"
    #,  " HR " = "Copollutant Model", " 95% CI " = "Copollutant Model", " p-value " = "Copollutant Model"
    ) |>
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
                            "Adjusted model" = "adjusted",
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
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
                            "Adjusted model" = "adjusted"
                            #, "Copollutant model" = "copollutant"
                            ),
         model = fct_relevel(model, 'Base model', 'Adjusted model'
                             #, 'Copollutant model'
                             ), 
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
         POPs_quart_ALS_figure_danish = POPs_quart_ALS_figure_danish, 
         plot_base_cox_gam_danish = plot_base_cox_gam_danish, 
         plot_adjusted_cox_gam_danish = plot_adjusted_cox_gam_danish, 
         plot_copollutant_cox_gam_danish = plot_copollutant_cox_gam_danish), 
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
   plot_base_cox_gam_danish, 
   plot_adjusted_cox_gam_danish, 
   plot_copollutant_cox_gam_danish,
   
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

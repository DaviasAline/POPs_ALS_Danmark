# Aline Davias
# April 29, 2025 
# Analysis of survival after ALS diagnosis depending on POPs levels  

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")

# Creation of cases specific datasets ----
POPs_group_bis <- setdiff(POPs_group, "ΣPBDE")
bdd_cases_tot <- bdd |>
  filter (als == 1) |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, diagnosis_age, baseline_age, marital_status_2cat, smoking_2cat, bmi, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_bis),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))
rm(POPs_group_bis)

bdd_cases_danish <- bdd_danish |>
  filter (als == 1) |>
  filter(study == "Danish") |>
  select(study, study_2cat, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_group), 
         all_of(POPs_included)) |>
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
  replace_with_median(ΣPBDE, ΣPBDE_quart) |>
  mutate(across(all_of(POPs_included), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_included),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_28, PCB_28_quart) |>
  replace_with_median(PCB_52, PCB_52_quart) |>
  replace_with_median(PCB_74, PCB_74_quart) |>
  replace_with_median(PCB_99, PCB_99_quart) |>
  replace_with_median(PCB_101, PCB_101_quart) |>
  replace_with_median(PCB_118, PCB_118_quart) |>
  replace_with_median(PCB_138, PCB_138_quart) |>
  replace_with_median(PCB_153, PCB_153_quart) |>
  replace_with_median(PCB_156, PCB_156_quart) |>
  replace_with_median(PCB_170, PCB_170_quart) |>
  replace_with_median(PCB_180, PCB_180_quart) |>
  replace_with_median(PCB_183, PCB_183_quart) |>
  replace_with_median(PCB_187, PCB_187_quart) |>
  replace_with_median(OCP_pp_DDT, OCP_pp_DDT_quart) |>
  replace_with_median(OCP_pp_DDE, OCP_pp_DDE_quart) |>
  replace_with_median(OCP_oxychlordane, OCP_oxychlordane_quart) |>
  replace_with_median(OCP_transnonachlor, OCP_transnonachlor_quart) |>
  replace_with_median(PBDE_47, PBDE_47_quart) |>
  replace_with_median(PBDE_99, PBDE_99_quart) |>
  replace_with_median(PBDE_153, PBDE_153_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))

bdd_cases_finnish <- bdd_finnish |>
  filter (als == 1) |>
  filter(study %in% c("FMC", "FMCF", "MFH")) |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
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
  replace_with_median(ΣHCH, ΣHCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

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
  replace_with_median(ΣHCH, ΣHCH_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

bdd_cases_FMCF <- bdd_finnish |>
  filter (als == 1) |>
  filter(study == "FMCF") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol, education_merged,  
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
  replace_with_median(ΣHCH, ΣHCH_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

bdd_cases_MFH <- bdd_finnish |>
  filter (als == 1) |>
  filter(study == "MFH") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol, education_merged,   
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
  replace_with_median(ΣHCH, ΣHCH_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

surv_obj_tot <- Surv(time = bdd_cases_tot$follow_up_death,                      # set the outcomes
                     event = bdd_cases_tot$status_death)
surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                        event = bdd_cases_danish$status_death)
surv_obj_finnish <- Surv(time = bdd_cases_finnish$follow_up_death,               
                     event = bdd_cases_finnish$status_death)
surv_obj_FMC <- Surv(time = bdd_cases_FMC$follow_up_death,               
                        event = bdd_cases_FMC$status_death)
surv_obj_FMCF <- Surv(time = bdd_cases_FMCF$follow_up_death,              
                        event = bdd_cases_FMCF$status_death)
surv_obj_MFH <- Surv(time = bdd_cases_MFH$follow_up_death,               
                        event = bdd_cases_MFH$status_death)

covariates_danish <- c("sex", "diagnosis_age", "smoking_2cat_i", "bmi", "marital_status_2cat_i")
covariates_finnish <- c("sex", "diagnosis_age", "smoking_2cat", "bmi", "marital_status_2cat")

# Danish cohort ----
## Covar model ----
covar_danish <- tbl_merge(tbls = list(
  tbl_uvregression(
    data = bdd_cases_danish, 
    y = surv_obj_danish, 
    method = survival::coxph,  
    exponentiate = TRUE,
    include = c("sex", "diagnosis_age", 
                "bmi", "marital_status_2cat_i", "smoking_i", "education_i", "cholesterol_i")) |>
    bold_labels() |>
    bold_p(), 
  tbl_regression(
    coxph(surv_obj_danish ~ sex + diagnosis_age + bmi + marital_status_2cat_i + smoking_i + education_i + cholesterol_i, data = bdd_cases_danish),
    exponentiate = TRUE) |>
    bold_labels() |>
    bold_p()), 
  tab_spanner = c("**Crude**", "**Adjusted**"))

## Main analysis 
## Cox model (sd) ----
### Base ----
model1_cox_sd_danish <- map_dfr(POPs_group_sd, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex"))       # set the formulas                
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
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
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

model <- coxph(formula_danish, data = bdd_cases_danish) 
cheking_model3_cox_sd_danish <- check_model(model, residual_type = "normal")

rm(POPs_group_sd_bis, pollutant_labels_bis, formula_danish, model_summary, model, coefs)


## Cox model (quart) ----
### Base ----
model1_cox_quart_danish <- map_dfr(POPs_group_quart, function(expl) {
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex"))  
  
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i,  
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
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

### Heterogeneity tests ----
#### base ----
heterogeneity_base_quart <- data.frame(explanatory = character(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (expl in POPs_group_quart) {
  
  formula_raw <- as.formula("surv_obj_danish ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_danish) 
  
  formula <- as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart <- rbind(heterogeneity_base_quart, 
                                    data.frame(explanatory = expl,
                                               model = "base",
                                               p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart <- data.frame(explanatory = character(),
                                           model = factor(),
                                           p.value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (expl in POPs_group_quart) {
  
  formula_raw <- as.formula(paste("surv_obj_danish ~ ", paste(covariates_danish, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_danish) 
  
  formula <- as.formula(paste("surv_obj_danish ~", expl, "+ ", paste(covariates_danish, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(explanatory = expl,
                                                   model = "adjusted",
                                                   p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

#### copollutant ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL_full <- 
  gam(outcome ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_danish) 

model3_quart_PCB_DL_raw <- 
  gam(outcome ~ 
        #PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_danish) 

anova <- anova(model3_quart_PCB_DL_raw, model3_quart_PCB_DL_full, test = "Chisq")
p.value_heterogeneity_PCB_DL <- tibble(explanatory = "PCB_DL_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_PCB_NDL_full <- 
  gam(outcome ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_PCB_NDL_raw <- 
  gam(outcome ~ 
        #PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_PCB_NDL_raw, model3_quart_PCB_NDL_full, test = "Chisq")
p.value_heterogeneity_PCB_NDL <- tibble(explanatory = "PCB_NDL_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_HCB_full <- 
  gam(outcome ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_HCB_raw <- 
  gam(outcome ~ 
        #OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_HCB_raw, model3_quart_HCB_full, test = "Chisq")
p.value_heterogeneity_HCB <- tibble(explanatory = "OCP_HCB_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])


model3_quart_ΣDDT_full <- 
  gam(outcome ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

model3_quart_ΣDDT_raw <- 
  gam(outcome ~ 
        #ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

anova <- anova(model3_quart_ΣDDT_raw, model3_quart_ΣDDT_full, test = "Chisq")
p.value_heterogeneity_ΣDDT <- tibble(explanatory = "ΣDDT_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_β_HCH_full <- 
  gam(outcome ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

model3_quart_β_HCH_raw <- 
  gam(outcome ~ 
        #OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

anova <- anova(model3_quart_β_HCH_raw, model3_quart_β_HCH_full, test = "Chisq")
p.value_heterogeneity_β_HCH <- tibble(explanatory = "OCP_β_HCH_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_Σchlordane_full <- 
  gam(outcome ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_Σchlordane_raw <- 
  gam(outcome ~ 
        #Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_Σchlordane_raw, model3_quart_Σchlordane_full, test = "Chisq")
p.value_heterogeneity_Σchlordane <- tibble(explanatory = "Σchlordane_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_ΣPBDE_full <- 
  gam(outcome ~ 
        ΣPBDE_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

model3_quart_ΣPBDE_raw <- 
  gam(outcome ~ 
        #ΣPBDE_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) 

anova <- anova(model3_quart_ΣPBDE_raw, model3_quart_ΣPBDE_full, test = "Chisq")
p.value_heterogeneity_ΣPBDE <- tibble(explanatory = "ΣPBDE_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

heterogeneity_copollutant_quart <- bind_rows(p.value_heterogeneity_PCB_DL, 
                                       p.value_heterogeneity_PCB_NDL, 
                                       p.value_heterogeneity_HCB, 
                                       p.value_heterogeneity_ΣDDT, 
                                       p.value_heterogeneity_β_HCH, 
                                       p.value_heterogeneity_Σchlordane, 
                                       p.value_heterogeneity_ΣPBDE) |>
  mutate(model = "copollutant")

rm(anova, outcome, 
   model3_quart_PCB_DL_full, model3_quart_PCB_DL_raw, 
   model3_quart_PCB_NDL_full, model3_quart_PCB_NDL_raw, 
   model3_quart_HCB_full, model3_quart_HCB_raw, 
   model3_quart_ΣDDT_full, model3_quart_ΣDDT_raw, 
   model3_quart_β_HCH_full, model3_quart_β_HCH_raw, 
   model3_quart_Σchlordane_full, model3_quart_Σchlordane_raw, 
   model3_quart_ΣPBDE_full, model3_quart_ΣPBDE_raw,
   p.value_heterogeneity_PCB_DL, 
   p.value_heterogeneity_PCB_NDL, 
   p.value_heterogeneity_HCB, 
   p.value_heterogeneity_ΣDDT, 
   p.value_heterogeneity_β_HCH, 
   p.value_heterogeneity_Σchlordane, 
   p.value_heterogeneity_ΣPBDE)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_quart, 
            heterogeneity_adjusted_quart, 
            heterogeneity_copollutant_quart) |>
  mutate(explanatory = gsub("_quart", "", explanatory), 
         study = "Danish")

### Trend tests ----
#### base ----
trend_base <- data.frame(explanatory = character(),
                         model = factor(), 
                         p.value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med) {
  
  formula <- as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_danish) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base <- rbind(trend_base, 
                      data.frame(explanatory = expl,
                                 model = "base",
                                 p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

#### adjusted ----
trend_adjusted <- data.frame(explanatory = character(),
                             model = factor(), 
                             p.value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med) {
  
  formula <- as.formula(paste("surv_obj_danish ~", expl, "+ ", paste(covariates_danish, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_danish) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(explanatory = expl,
                                     model = "adjusted",
                                     p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

#### copollutant ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL_trend <- 
  gam(outcome ~ 
        PCB_DL_quart_med + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_PCB_DL <- model3_quart_PCB_DL_trend$p.table["PCB_DL_quart_med", "Pr(>|z|)"]

model3_quart_PCB_NDL_trend <- 
  gam(outcome ~ 
        PCB_NDL_quart_med + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)|> 
  summary()
p.value_trend_PCB_NDL <- model3_quart_PCB_NDL_trend$p.table["PCB_NDL_quart_med", "Pr(>|z|)"]

model3_quart_HCB_trend <- 
  gam(outcome ~ 
        OCP_HCB_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)|> 
  summary()
p.value_trend_HCB <- model3_quart_HCB_trend$p.table["OCP_HCB_quart_med", "Pr(>|z|)"]

model3_quart_ΣDDT_trend <- 
  gam(outcome ~ 
        ΣDDT_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_ΣDDT <- model3_quart_ΣDDT_trend$p.table["ΣDDT_quart_med", "Pr(>|z|)"]

model3_quart_β_HCH_trend <- 
  gam(outcome ~ 
        OCP_β_HCH_quart_med +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_β_HCH <- model3_quart_β_HCH_trend$p.table["OCP_β_HCH_quart_med", "Pr(>|z|)"]

model3_quart_Σchlordane_trend <- 
  gam(outcome ~ 
        Σchlordane_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_Σchlordane <- model3_quart_Σchlordane_trend$p.table["Σchlordane_quart_med", "Pr(>|z|)"]

model3_quart_ΣPBDE_trend <- 
  gam(outcome ~ 
        ΣPBDE_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_ΣPBDE <- model3_quart_ΣPBDE_trend$p.table["ΣPBDE_quart_med", "Pr(>|z|)"]

trend_copollutant <- 
  data.frame(explanatory = c("PCB_DL_quart_med", "PCB_NDL_quart_med", "OCP_HCB_quart_med", 
                             "ΣDDT_quart_med", "OCP_β_HCH_quart_med", "Σchlordane_quart_med", 
                             "ΣPBDE_quart_med" ),
             model = "copollutant",
             p.value_trend = c(p.value_trend_PCB_DL, 
                               p.value_trend_PCB_NDL, 
                               p.value_trend_HCB, 
                               p.value_trend_ΣDDT, 
                               p.value_trend_β_HCH, 
                               p.value_trend_Σchlordane, 
                               p.value_trend_ΣPBDE))

rm(outcome, 
   model3_quart_PCB_DL_trend,
   model3_quart_PCB_NDL_trend, 
   model3_quart_HCB_trend, 
   model3_quart_ΣDDT_trend, 
   model3_quart_β_HCH_trend,
   model3_quart_Σchlordane_trend, 
   model3_quart_ΣPBDE_trend, 
   
   p.value_trend_PCB_DL, 
   p.value_trend_PCB_NDL, 
   p.value_trend_HCB, 
   p.value_trend_ΣDDT, 
   p.value_trend_β_HCH, 
   p.value_trend_Σchlordane, 
   p.value_trend_ΣPBDE)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted, trend_copollutant) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory), 
         study = "Danish")

rm(heterogeneity_base_quart, heterogeneity_adjusted_quart, heterogeneity_copollutant_quart,
   trend_base, trend_adjusted, trend_copollutant)


## Cox-gam model ----
### Base  ----
POPs_group_labels <- set_names(
  c("Most prevalent PCBs", "Dioxin-like PCBs","Non-dioxin-like PCBs", "HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)
plot_base_cox_gam_danish <- map(POPs_group, function(var) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula <- as.formula(paste("outcome ~ s(", var, ") + diagnosis_age + sex")) 
  
  model <- gam(formula,                                                         # run the cox-gam model
               family = cox.ph(), 
               method = "ML", 
               data = bdd_cases_danish)            
  
  bdd_pred <- bdd_cases_danish |>                                               # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_diagnosis_age = mean(diagnosis_age, na.rm = TRUE),
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
  x_label <- POPs_group_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Hazard Ratio (HR)") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() + 
    scale_y_log10(limits = c(10, 15000)) +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Base model")
                                              # distribution of the non-scaled POP variables 
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var]]) +
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
plot_adjusted_cox_gam_danish <- map(POPs_group, function(var) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula <- as.formula(paste("outcome ~ s(", var, ") + ", paste(covariates_danish, collapse = "+"))) 
  
  model <- gam(formula,                                                         # run the cox-gam model
               family = cox.ph(), 
               method = "ML", 
               data = bdd_cases_danish)     
  
  bdd_pred <- bdd_cases_danish |>                                               # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_diagnosis_age = mean(diagnosis_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i))), 
      adj_bmi = mean(bmi, na.rm = TRUE),  
      adj_marital_status_2cat_i = names(which.max(table(marital_status_2cat_i)))) |>
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
  x_label <- POPs_group_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = hazard_ratio)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Hazard Ratio (HR)") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() + 
    scale_y_log10(limits = c(100, 100000),
                  labels = scales::label_number(accuracy = 1)) +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
                                             # distribution of non-scaled POP variables 
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var]]) +
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
POPs_group_bis <- setdiff(POPs_group, "PCB_4")                         # PCB4 not included because already present in NDL-PCBs
POPs_group_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_bis)

outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model <- gam(outcome ~ s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + 
               s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) + 
               sex + diagnosis_age + 
               smoking_2cat_i + bmi + marital_status_2cat_i, 
             family = cox.ph(), 
             method = "ML", 
             data = bdd_cases_danish)


plot_copollutant_cox_gam_danish <- map(POPs_group_bis, function(var) {
  
  bdd_pred <- bdd_cases_danish |>
    mutate(across(all_of(covariates_danish),                                    # fixe toutes les covariables à leurs moyennes
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))))) |>
    mutate(across(setdiff(POPs_group_bis, var), 
                  ~ mean(., na.rm = TRUE)))|>                                   # Fixe tous les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates_danish), all_of(setdiff(POPs_group_bis, var)))  
  
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
  x_label <- POPs_group_labels_bis[var]
  
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
  
  p2 <- bdd_cases_danish |>
    ggplot() +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
}) |> set_names(POPs_group_bis)
rm(POPs_group_bis, POPs_group_labels_bis, model, outcome)

## Q-gcomp analysis ---- 
### All pollutants (grouped) ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")                                  # remove the 4 most abundant PCB because they are already NDL-PCB

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(covariates_danish, collapse = " + ")))
set.seed(1996)
qgcomp_boot_danish <-
  qgcomp.cox.boot(
    f = formula_danish,                                                         # formula
    data =   bdd_cases_danish, 
    q = 4,                                                                      # nb of quantiles
    expnms = POPs_group_bis,                                                    # exposures of interest 
    B = 1000,                                                                   # nb of boostrap
    MCsize = 5000,
    seed = 1996,
    parallel = TRUE,                                                            # shorter run time
    parplan = TRUE)                                                             # shorter run time
# print(qgcomp_boot_danish)
# plot(qgcomp_boot_danish)

# run the code without bootsrapping just to get weights even if the mixture is not significant 
formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(covariates_danish, POPs_group_bis), collapse = " + ")))
qgcomp_noboot_danish <-                                                         
  qgcomp.cox.noboot(
    f = formula_danish,                                                         # formula
    bdd_cases_danish[, c(POPs_group_bis, covariates_danish, 'follow_up_death', 'status_death')],
    q = 4,                                                                      # number of quantiles
    expnms = POPs_group_bis)                                                    # exposures of interest
# print(qgcomp_noboot_danish)
# qgcomp_noboot_danish$pos.weights
# qgcomp_noboot_danish$neg.weights
# plot(qgcomp_noboot_danish, suppressprint = TRUE)
rm(formula_danish, POPs_group_bis)

### Just the positives, adjusted for the others (grouped) ----
# run qgcomp just for the positive weigths and by adjusted for the negatives 
formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(covariates_danish,                                     # attention pour qgcomp.cox.boot() pas besoin de mettre les expo dans la formule sinon sont comptés 2 fois
                           # c("PCB_NDL", "OCP_β_HCH", "Σchlordane"), 
                           c("PCB_DL_quart",  "OCP_HCB_quart", "ΣDDT_quart", "ΣPBDE_quart")), 
                         collapse = " + ")))
set.seed(1996)
qgcomp_positive_boot_danish <-
  qgcomp.cox.boot(
    f = formula_danish,                                                         # formula
    data =   bdd_cases_danish, 
    q = 4,                                                                      # nb of quantiles
    expnms = c( "PCB_NDL", "OCP_β_HCH", "Σchlordane"),                          # exposures of interest 
    B = 1000,                                                                   # nb of boostrap
    MCsize = 5000,
    seed = 1996,
    parallel = TRUE,                                                            # shorter run time
    parplan = TRUE)                                                             # shorter run time
# print(qgcomp_positive_boot_danish)
# plot(qgcomp_positive_boot_danish)


formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
             paste(c(covariates_danish,                                         # attention pour qgcomp.cox.noboot() il faut mettre les expo dans la formule 
                     c("PCB_NDL", "OCP_β_HCH", "Σchlordane"), 
                     c("PCB_DL_quart",  "OCP_HCB_quart", "ΣDDT_quart", "ΣPBDE_quart")), 
                   collapse = " + ")))
qgcomp_positive_noboot_danish <-                                                         
  qgcomp.cox.noboot(
    f = formula_danish,                                                         # formula
    bdd_cases_danish[, c(
      c("PCB_NDL", "OCP_β_HCH", "Σchlordane"), 
      c("PCB_DL_quart",  "OCP_HCB_quart", "ΣDDT_quart", "ΣPBDE_quart"), 
      covariates_danish, 
      'follow_up_death', 'status_death')],
    q = 4,                                                                      # number of quantiles
    expnms = c("PCB_NDL", "OCP_β_HCH", "Σchlordane"))                           # exposures of interest
# print(qgcomp_positive_noboot_danish)
# qgcomp_positive_noboot_danish$pos.weights
# qgcomp_positive_noboot_danish$neg.weights
# plot(qgcomp_noboot_danish, suppressprint = TRUE)
rm(formula_danish, POPs_group_bis)


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
covar_crude_finnish <- map_dfr(
  c("diagnosis_age", "sex", "level_urbanization", "smoking_2cat", "bmi", "marital_status_2cat", "cholesterol"), function(expl) {
  
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
formula <- "Surv(follow_up_death, status_death) ~ diagnosis_age + sex + level_urbanization + smoking_2cat + bmi + marital_status_2cat + cholesterol"
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
             "marital_status_2cat" = "marital_status_2catOther",
             "sex" = "sexFemale",
             "smoking_2cat" = "smoking_2catNever")) |>
  relocate(model, variable, term) 

covar_finnish <- bind_rows(covar_crude_finnish, covar_adjusted_finnish)
rm(formula, covar_crude_finnish, covar_adjusted_finnish)

## Cox model (sd) - meta-analysis ----
### Base ----
model1_cox_sd_finnish <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_FMC <-                                                                # set the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ diagnosis_age + sex"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ diagnosis_age + sex"))
  
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
    study_design = "meta-analysis", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp)
  )
})

### Adjusted ----
model2_cox_sd_finnish <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                                  paste(covariates_finnish, collapse = " + ")))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                   paste(covariates_finnish, collapse = " + ")))
  
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
    study_design = "meta-analysis", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp)
    
  )
})

### Copollutant ----
POPs_group_sd_finnish_bis <- setdiff(POPs_group_sd_finnish, c("OCP_β_HCH_sd",  "OCP_γ_HCH_sd", "PCB_4_sd"))

formula_FMC <-                                                                  # creation of the formulas
  as.formula(paste("surv_obj_FMC ~", paste(c(POPs_group_sd_finnish_bis, covariates_finnish), collapse = "+")))  
formula_FMCF <- 
  as.formula(paste("surv_obj_FMCF ~", paste(c(POPs_group_sd_finnish_bis, covariates_finnish), collapse = "+")))  

model3_cox_sd_finnish <- list(                                                                # run of the simple cox models 
  finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2, 
         explanatory = term) |>
  filter(str_detect(term, "_sd"))                                               # remove the covariates results

model3_cox_sd_finnish <- model3_cox_sd_finnish |>
  group_by(term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      study = "Finnish", 
      model = "copollutant", 
      study_design = "meta-analysis", 
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |>
  mutate(explanatory = term)

rm(POPs_group_sd_finnish_bis, formula_FMC, formula_FMCF)

## Cox model (quart) - meta-analysis ----
### Base ----
model1_cox_quart_finnish <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_FMC <-                                                                # creation of the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ diagnosis_age + sex"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ diagnosis_age + sex"))
  
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
        `p-value` = as.numeric(rma_fit$pval),
        p.value_heterogeneity = as.numeric(rma_fit$QEp))
    }) |> 
    ungroup() |> 
    mutate(study = "Finnish", model = "base", study_design = "meta-analysis") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

### Adjusted ----
model2_cox_quart_finnish <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_FMC <- as.formula(paste("surv_obj_FMC ~", expl, "+",                  # set the formulas              
                                  paste(covariates_finnish, collapse = " + ")))
  formula_FMCF <- as.formula(paste("surv_obj_FMCF ~", expl, "+",                            
                                   paste(covariates_finnish, collapse = " + ")))
  
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
        `p-value` = as.numeric(rma_fit$pval),
        p.value_heterogeneity = as.numeric(rma_fit$QEp))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted", study_design = "meta-analysis") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

### Copollutant ---- 
POPs_group_quart_finnish_bis <- setdiff(POPs_group_quart_finnish, c("OCP_β_HCH_quart",  "OCP_γ_HCH_quart", "PCB_4_quart"))

formula_FMC <- as.formula(paste("surv_obj_FMC ~ ",                  # set the formulas              
                                paste(c(POPs_group_quart_finnish_bis, covariates_finnish), 
                                      collapse = " + ")))
formula_FMCF <- as.formula(paste("surv_obj_FMCF ~ ",                             
                                 paste(c(POPs_group_quart_finnish_bis, covariates_finnish), 
                                       collapse = " + ")))

model3_cox_quart_finnish <- list(                                               # run of the simple cox model
  finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2, 
         explanatory = term) |>
  filter(str_detect(term, "_quart"))                                            # remove the covariates results

model3_cox_quart_finnish <- model3_cox_quart_finnish |>                         # run metanalyse (one per quartile per explanatory variable)
  group_by(explanatory, term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      study = "Finnish", 
      model = "copollutant", 
      study_design = "meta-analysis",
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
      p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |> 
  mutate(
    explanatory = gsub("Q2", "", explanatory),
    explanatory = gsub("Q3", "", explanatory),
    explanatory = gsub("Q4", "", explanatory)) |>
  relocate(model, explanatory, term)

rm(POPs_group_quart_finnish_bis, formula_FMC, formula_FMCF)

## Cox model (sd) - raw - not adjusted on cohort ----
### Base ----
model1_cox_sd_finnish_raw <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <- 
    as.formula(paste("surv_obj_finnish ~", expl, "+ diagnosis_age + sex"))      # set the formulas                
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "base", 
    study_design = "raw unadjusted on cohort", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Adjusted ----
model2_cox_sd_finnish_raw <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <-                                                            # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "+",  
                     paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "raw unadjusted on cohort", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})

### Copollutant ----
POPs_group_sd_finnish_bis <- setdiff(POPs_group_sd_finnish, c("OCP_β_HCH_sd",  "OCP_γ_HCH_sd", "PCB_4_sd"))

formula_finnish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(POPs_group_sd_finnish_bis, covariates_finnish), collapse = " + ")))

model_summary <- 
  coxph(formula_finnish, data = bdd_cases_finnish) |> summary() 
coefs <- model_summary$coefficients
model3_cox_sd_finnish_raw <- tibble(                                                 # creation of a table of results
  study = "Finnish", 
  model = "copollutant", 
  study_design = "raw unadjusted on cohort", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

model <- coxph(formula_finnish, data = bdd_cases_finnish) 
cheking_model3_cox_sd_finnish <- check_model(model, residual_type = "normal")

rm(POPs_group_sd_finnish_bis, formula_finnish, model_summary, model, coefs)


## Cox model (quart) - raw - not adjusted on cohort ----
### Base ----
model1_cox_quart_finnish_raw <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_finnish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "+ diagnosis_age + sex"))  
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()  # run cox model
  
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "base", 
    study_design = "raw unadjusted on cohort",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                      # remove the covariates results
})

### Adjusted ----
model2_cox_quart_finnish_raw <- map_dfr(POPs_group_quart_finnish, function(expl) {
  
  formula_finnish <- as.formula(paste("surv_obj_finnish ~", expl, "+",            # set the formulas              
                                      paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "raw unadjusted on cohort",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})


### Copollutant ----
outcome <- with(bdd_cases_finnish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL <- 
  gam(outcome ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_finnish) |>
  summary()
model3_quart_PCB_DL <- model3_quart_PCB_DL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_PCB_NDL <- 
  gam(outcome ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_PCB_NDL <- model3_quart_PCB_NDL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_HCB <- 
  gam(outcome ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_HCB <- model3_quart_HCB$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣDDT <- 
  gam(outcome ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_ΣDDT <- model3_quart_ΣDDT$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣHCH <- 
  gam(outcome ~ 
        ΣHCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_ΣHCH <- model3_quart_ΣHCH$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_Σchlordane <- 
  gam(outcome ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_Σchlordane <- model3_quart_Σchlordane$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_OCP_PeCB <- 
  gam(outcome ~ 
        OCP_PeCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |>
  summary()
model3_quart_OCP_PeCB <- model3_quart_OCP_PeCB$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_cox_quart_finnish_raw <- bind_rows(
  model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_ΣHCH, model3_quart_Σchlordane, model3_quart_OCP_PeCB) |>
  filter(grepl("quart", variable)) |>
  mutate(
    study = "Finnish", 
    model = "copollutant", 
    study_design = "raw unadjusted on cohort", 
    term = variable,
    explanatory = gsub("Q2", "", variable), 
    explanatory = gsub("Q3", "", explanatory), 
    explanatory = gsub("Q4", "", explanatory), 
    coef = Estimate, 
    se = `Std. Error`, 
    `p-value` =`Pr(>|z|)`) |> 
  select(study, model, study_design, term, explanatory, coef, se, `p-value`)

rm(model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_ΣHCH, model3_quart_Σchlordane, model3_quart_OCP_PeCB, 
   outcome)

### Heterogeneity tests ----
#### base ----
heterogeneity_base_quart_finnish <- data.frame(explanatory = character(),
                                               model = factor(),
                                               p.value_heterogeneity = numeric(), 
                                               stringsAsFactors = FALSE)

for (expl in POPs_group_quart_finnish) {
  
  formula_raw <- as.formula("surv_obj_finnish ~ diagnosis_age + sex")
  model_raw <- coxph(formula_raw, data = bdd_cases_finnish) 
  
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_finnish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart_finnish <- rbind(heterogeneity_base_quart_finnish, 
                                            data.frame(explanatory = expl,
                                                       model = "base",
                                                       p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

#### adjusted ----
heterogeneity_adjusted_quart_finnish <- data.frame(explanatory = character(),
                                                   model = factor(),
                                                   p.value_heterogeneity = numeric(), 
                                                   stringsAsFactors = FALSE)

for (expl in POPs_group_quart_finnish) {
  
  formula_raw <- as.formula(paste("surv_obj_finnish ~ ", paste(covariates_finnish, collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_finnish) 
  
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ ", paste(covariates_finnish, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_finnish)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart_finnish <- rbind(heterogeneity_adjusted_quart_finnish, 
                                                data.frame(explanatory = expl,
                                                           model = "adjusted",
                                                           p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

#### copollutant ----
outcome <- with(bdd_cases_finnish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL_full <- 
  gam(outcome ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_finnish) 

model3_quart_PCB_DL_raw <- 
  gam(outcome ~ 
        #PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_finnish) 

anova <- anova(model3_quart_PCB_DL_raw, model3_quart_PCB_DL_full, test = "Chisq")
p.value_heterogeneity_PCB_DL <- tibble(explanatory = "PCB_DL_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_PCB_NDL_full <- 
  gam(outcome ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

model3_quart_PCB_NDL_raw <- 
  gam(outcome ~ 
        #PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

anova <- anova(model3_quart_PCB_NDL_raw, model3_quart_PCB_NDL_full, test = "Chisq")
p.value_heterogeneity_PCB_NDL <- tibble(explanatory = "PCB_NDL_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_HCB_full <- 
  gam(outcome ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

model3_quart_HCB_raw <- 
  gam(outcome ~ 
        #OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

anova <- anova(model3_quart_HCB_raw, model3_quart_HCB_full, test = "Chisq")
p.value_heterogeneity_HCB <- tibble(explanatory = "OCP_HCB_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])


model3_quart_ΣDDT_full <- 
  gam(outcome ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

model3_quart_ΣDDT_raw <- 
  gam(outcome ~ 
        #ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

anova <- anova(model3_quart_ΣDDT_raw, model3_quart_ΣDDT_full, test = "Chisq")
p.value_heterogeneity_ΣDDT <- tibble(explanatory = "ΣDDT_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_ΣHCH_full <- 
  gam(outcome ~ 
        ΣHCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

model3_quart_ΣHCH_raw <- 
  gam(outcome ~ 
        #ΣHCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

anova <- anova(model3_quart_ΣHCH_raw, model3_quart_ΣHCH_full, test = "Chisq")
p.value_heterogeneity_ΣHCH <- tibble(explanatory = "ΣHCH_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_Σchlordane_full <- 
  gam(outcome ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

model3_quart_Σchlordane_raw <- 
  gam(outcome ~ 
        #Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)

anova <- anova(model3_quart_Σchlordane_raw, model3_quart_Σchlordane_full, test = "Chisq")
p.value_heterogeneity_Σchlordane <- tibble(explanatory = "Σchlordane_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_OCP_PeCB_full <- 
  gam(outcome ~ 
        OCP_PeCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

model3_quart_OCP_PeCB_raw <- 
  gam(outcome ~ 
        #OCP_PeCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) 

anova <- anova(model3_quart_OCP_PeCB_raw, model3_quart_OCP_PeCB_full, test = "Chisq")
p.value_heterogeneity_OCP_PeCB <- tibble(explanatory = "OCP_PeCB_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

heterogeneity_copollutant_quart_finnish <- 
  bind_rows(p.value_heterogeneity_PCB_DL, 
            p.value_heterogeneity_PCB_NDL, 
            p.value_heterogeneity_HCB, 
            p.value_heterogeneity_ΣDDT, 
            p.value_heterogeneity_ΣHCH, 
            p.value_heterogeneity_Σchlordane, 
            p.value_heterogeneity_OCP_PeCB) |>
  mutate(model = "copollutant")

rm(anova, outcome, 
   model3_quart_PCB_DL_full, model3_quart_PCB_DL_raw, 
   model3_quart_PCB_NDL_full, model3_quart_PCB_NDL_raw, 
   model3_quart_HCB_full, model3_quart_HCB_raw, 
   model3_quart_ΣDDT_full, model3_quart_ΣDDT_raw, 
   model3_quart_ΣHCH_full, model3_quart_ΣHCH_raw, 
   model3_quart_Σchlordane_full, model3_quart_Σchlordane_raw, 
   model3_quart_OCP_PeCB_full, model3_quart_OCP_PeCB_raw,
   p.value_heterogeneity_PCB_DL, 
   p.value_heterogeneity_PCB_NDL, 
   p.value_heterogeneity_HCB, 
   p.value_heterogeneity_ΣDDT, 
   p.value_heterogeneity_ΣHCH, 
   p.value_heterogeneity_Σchlordane, 
   p.value_heterogeneity_OCP_PeCB)

heterogeneity_tests_finnish <- 
  bind_rows(heterogeneity_base_quart_finnish, 
            heterogeneity_adjusted_quart_finnish, 
            heterogeneity_copollutant_quart_finnish) |>
  mutate(explanatory = gsub("_quart", "", explanatory), 
         study = "Finnish")

### Trend tests ----
#### base ----
trend_base_finnish <- data.frame(explanatory = character(),
                                 model = factor(), 
                                 p.value_trend = numeric(), 
                                 stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med_finnish) {
  
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ diagnosis_age + sex"))  
  model <- coxph(formula, data = bdd_cases_finnish) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_base_finnish <- rbind(trend_base_finnish, 
                              data.frame(explanatory = expl,
                                         model = "base",
                                         p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

#### adjusted ----
trend_adjusted_finnish <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p.value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med_finnish) {
  
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ ", paste(covariates_finnish, collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_finnish) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_adjusted_finnish <- rbind(trend_adjusted_finnish, 
                                  data.frame(explanatory = expl,
                                             model = "adjusted",
                                             p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

#### copollutant ----
outcome <- with(bdd_cases_finnish, cbind(follow_up_death, status_death))

model3_quart_PCB_DL_trend <- 
  gam(outcome ~ 
        PCB_DL_quart_med + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',                                                            # maximum likelihood
      data = bdd_cases_finnish) |> 
  summary()
p.value_trend_PCB_DL <- model3_quart_PCB_DL_trend$p.table["PCB_DL_quart_med", "Pr(>|z|)"]

model3_quart_PCB_NDL_trend <- 
  gam(outcome ~ 
        PCB_NDL_quart_med + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)|> 
  summary()
p.value_trend_PCB_NDL <- model3_quart_PCB_NDL_trend$p.table["PCB_NDL_quart_med", "Pr(>|z|)"]

model3_quart_HCB_trend <- 
  gam(outcome ~ 
        OCP_HCB_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish)|> 
  summary()
p.value_trend_HCB <- model3_quart_HCB_trend$p.table["OCP_HCB_quart_med", "Pr(>|z|)"]

model3_quart_ΣDDT_trend <- 
  gam(outcome ~ 
        ΣDDT_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣHCH) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |> 
  summary()
p.value_trend_ΣDDT <- model3_quart_ΣDDT_trend$p.table["ΣDDT_quart_med", "Pr(>|z|)"]

model3_quart_ΣHCH_trend <- 
  gam(outcome ~ 
        ΣHCH_quart_med +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat,  
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |> 
  summary()
p.value_trend_ΣHCH <- model3_quart_ΣHCH_trend$p.table["ΣHCH_quart_med", "Pr(>|z|)"]

model3_quart_Σchlordane_trend <- 
  gam(outcome ~ 
        Σchlordane_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(OCP_PeCB) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |> 
  summary()
p.value_trend_Σchlordane <- model3_quart_Σchlordane_trend$p.table["Σchlordane_quart_med", "Pr(>|z|)"]

model3_quart_OCP_PeCB_trend <- 
  gam(outcome ~ 
        OCP_PeCB_quart_med + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(ΣHCH) + s(Σchlordane) +
        sex + diagnosis_age + smoking_2cat + bmi + marital_status_2cat, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_finnish) |> 
  summary()
p.value_trend_OCP_PeCB <- model3_quart_OCP_PeCB_trend$p.table["OCP_PeCB_quart_med", "Pr(>|z|)"]

trend_copollutant_finnish <- 
  data.frame(explanatory = c("PCB_DL_quart_med", "PCB_NDL_quart_med", "OCP_HCB_quart_med", 
                             "ΣDDT_quart_med", "ΣHCH_quart_med", "Σchlordane_quart_med", 
                             "OCP_PeCB_quart_med" ),
             model = "copollutant",
             p.value_trend = c(p.value_trend_PCB_DL, 
                               p.value_trend_PCB_NDL, 
                               p.value_trend_HCB, 
                               p.value_trend_ΣDDT, 
                               p.value_trend_ΣHCH, 
                               p.value_trend_Σchlordane, 
                               p.value_trend_OCP_PeCB))

rm(outcome, 
   model3_quart_PCB_DL_trend,
   model3_quart_PCB_NDL_trend, 
   model3_quart_HCB_trend, 
   model3_quart_ΣDDT_trend, 
   model3_quart_ΣHCH_trend,
   model3_quart_Σchlordane_trend, 
   model3_quart_OCP_PeCB_trend, 
   
   p.value_trend_PCB_DL, 
   p.value_trend_PCB_NDL, 
   p.value_trend_HCB, 
   p.value_trend_ΣDDT, 
   p.value_trend_ΣHCH, 
   p.value_trend_Σchlordane, 
   p.value_trend_OCP_PeCB)

trend_tests_finnish <- 
  bind_rows(trend_base_finnish, trend_adjusted_finnish, trend_copollutant_finnish) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory), 
         study = "Finnish")

rm(heterogeneity_base_quart_finnish, heterogeneity_adjusted_quart_finnish, heterogeneity_copollutant_quart_finnish,
   trend_base_finnish, trend_adjusted_finnish, trend_copollutant_finnish)



## Q-gcomp analysis ---- 
POPs_group_finnish_bis <-                                                 # remove the 4 most abundant PCB because they are already NDL-PCB
  setdiff(POPs_group_finnish, 
          c("OCP_β_HCH",  "OCP_γ_HCH", "PCB_4"))                    

formula_finnish <-                                                              # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(covariates_finnish, "study"), collapse = " + ")))
set.seed(1996)
qgcomp_boot_finnish <-
  qgcomp.cox.boot(
    f = formula_finnish,                                                        # formula
    data =   bdd_cases_finnish, 
    q = 4,                                                                      # nb of quantiles
    expnms = POPs_group_finnish_bis,                                                 # exposures of interest 
    B = 1000,                                                                   # nb of boostrap
    MCsize = 5000,
    seed = 1996,
    parallel = TRUE,                                                            # shorter run time
    parplan = TRUE)                                                             # shorter run time
# print(qgcomp_boot_finnish)
# qgcomp_boot_finnish$pos.weights                                                  # NULL because the model is not significant 
# qgcomp_boot_finnish$neg.weights                                                  # NULL because the model is not significant 
# plot(qgcomp_boot_finnish)

# run the code without bootsrapping just to get weights even if the mixture is not significant 
formula_finnish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(covariates_finnish, POPs_group_finnish_bis, "study"), collapse = " + ")))
qgcomp_noboot_finnish <-                                                         
  qgcomp.cox.noboot(
    f = formula_finnish,                                                         # formula
    bdd_cases_finnish[, c(POPs_group_finnish_bis, covariates_finnish, 'follow_up_death', 'status_death', "study")],
    q = 4,                                                                      # number of quantiles
    expnms = POPs_group_finnish_bis)                                                 # exposures of interest
# print(qgcomp_noboot_finnish)
# qgcomp_noboot_finnish$pos.weights
# qgcomp_noboot_finnish$neg.weights
# plot(qgcomp_noboot_finnish, suppressprint = TRUE)

rm(formula_finnish, POPs_group_finnish_bis)

# Metaanalysis ----
## Cox model (sd) ----
POPs_group_sd_bis <- setdiff(POPs_group_sd, "ΣPBDE_sd")

### Base ----
model1_cox_sd_metanalysis <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex"))  
  formula_FMC <-                                                                
    as.formula(paste("surv_obj_FMC ~", expl, "+ diagnosis_age + sex"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ diagnosis_age + sex"))
  
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
    `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp)
    
  )
})

### Adjusted ----
model2_cox_sd_metanalysis <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+",
                     paste(covariates_danish, collapse = "+")))  
  formula_FMC <-                                                                
    as.formula(paste("surv_obj_FMC ~", expl,  "+",
                     paste(covariates_finnish, collapse = "+")))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl,  "+",
                     paste(covariates_finnish, collapse = "+")))
  
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
    model = "adjusted", 
    explanatory = expl,
    term = expl,
    HR = exp(as.numeric(rma_fit$beta)),
    lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
    upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
    `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp)
    
  )
})

### Copollutant ----
POPs_group_sd_bis <- setdiff(POPs_group_sd, c("PCB_4_sd", "ΣPBDE_sd"))

formula_danish <-                                                                  # creation of the formulas
  as.formula(paste("surv_obj_danish ~", paste(c(POPs_group_sd_bis, 
                                                covariates_danish), 
                                              collapse = "+")))  
formula_FMC <-                                                                  # creation of the formulas
  as.formula(paste("surv_obj_FMC ~", paste(c(POPs_group_sd_bis, 
                                             covariates_finnish), 
                                           collapse = "+")))  
formula_FMCF <- 
  as.formula(paste("surv_obj_FMCF ~", paste(c(POPs_group_sd_bis, covariates_finnish), 
                                            collapse = "+")))  

model3_cox_sd_metanalysis <- list(                                                                # run of the simple cox models 
  danish = run_cox(formula_danish, bdd_cases_danish),
  finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2, 
         explanatory = term) |>
  filter(str_detect(term, "_sd"))                                               # remove the covariates results

model3_cox_sd_metanalysis <- model3_cox_sd_metanalysis |>
  group_by(term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      study = "Metanalysis", 
      model = "copollutant", 
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
      p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |>
  mutate(explanatory = term)

rm(POPs_group_sd_bis, formula_danish, formula_FMC, formula_FMCF)


## Cox model (quart) ----
POPs_group_quart_bis <- setdiff(POPs_group_quart, "ΣPBDE_quart")

### Base ----
model1_cox_quart_metanalysis <- map_dfr(POPs_group_quart_bis, function(expl) {
  
  formula_danish <-                                                             # creation of the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+ diagnosis_age + sex")) 
  formula_FMC <-                                                                # creation of the formulas
    as.formula(paste("surv_obj_FMC ~", expl, "+ diagnosis_age + sex"))  
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+ diagnosis_age + sex"))
  
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
        `p-value` = as.numeric(rma_fit$pval),
        p.value_heterogeneity = as.numeric(rma_fit$QEp))
    }) |> 
    ungroup() |> 
    mutate(study = "Metanalysis", 
           model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

### Adjusted  ----
model2_cox_quart_metanalysis <- map_dfr(POPs_group_quart_bis, function(expl) {
  
  formula_danish <- 
    as.formula(paste("surv_obj_danish ~", expl, "+",                            # set the formulas 
                      paste(covariates_danish, collapse = "+")))
  formula_FMC <- 
    as.formula(paste("surv_obj_FMC ~", expl, "+", 
                     paste(covariates_finnish, collapse = "+")))
  formula_FMCF <- 
    as.formula(paste("surv_obj_FMCF ~", expl, "+", 
                     paste(covariates_finnish, collapse = "+")))
  
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
        `p-value` = as.numeric(rma_fit$pval),
        p.value_heterogeneity = as.numeric(rma_fit$QEp))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

### Copollutant ---- 
POPs_group_quart_bis <- setdiff(POPs_group_quart, c("PCB_4_quart", "ΣPBDE_quart"))

formula_danish <- as.formula(paste("surv_obj_danish ~ ",                        # set the formulas              
                                paste(c(POPs_group_quart_bis, covariates_danish), collapse = " + ")))
formula_FMC <- as.formula(paste("surv_obj_FMC ~ ",                               
                                paste(c(POPs_group_quart_bis, covariates_finnish), collapse = " + ")))
formula_FMCF <- as.formula(paste("surv_obj_FMCF ~ ",                             
                                 paste(c(POPs_group_quart_bis, covariates_finnish), collapse = " + ")))

model3_cox_quart_metanalysis <- list(                                               # run of the simple cox model
  danish = run_cox(formula_danish, bdd_cases_danish),
  finnish_FMC = run_cox(formula_FMC, bdd_cases_FMC),
  finnish_FMCF = run_cox(formula_FMCF, bdd_cases_FMCF)) |>
  bind_rows(.id = "dataset") |>
  mutate(var = se^2, 
         explanatory = term) |>
  filter(str_detect(term, "_quart"))                                            # remove the covariates results

model3_cox_quart_metanalysis <- model3_cox_quart_metanalysis |>                         # run metanalyse (one per quartile per explanatory variable)
  group_by(explanatory, term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      study = "Metanalysis", 
      model = "copollutant",
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
      p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |> 
  mutate(
    explanatory = gsub("Q2", "", explanatory),
    explanatory = gsub("Q3", "", explanatory),
    explanatory = gsub("Q4", "", explanatory)) |>
  relocate(model, explanatory, term)

rm(POPs_group_quart_bis, formula_danish, formula_FMC, formula_FMCF)

# Assemblage main analyses ----
main_results_POPs_ALS_survival <-       
  bind_rows(
    model1_cox_sd_danish, model2_cox_sd_danish, model3_cox_sd_danish,
    model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish, 
    model1_cox_sd_finnish_raw, model2_cox_sd_finnish_raw, model3_cox_sd_finnish_raw,  # raw analysis of the finnish cohorts (not adjusting for cohort)
    model1_cox_quart_finnish_raw, model2_cox_quart_finnish_raw, model3_cox_quart_finnish_raw) |> # raw analysis of the finnish cohorts (not adjusting for cohort)
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, study_design, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  
  bind_rows(model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish,  # meta-anlaysis of the finnish cohorts
            model1_cox_quart_finnish, model2_cox_quart_finnish, model3_cox_quart_finnish,  # meta-anlaysis of the finnish cohorts
            model1_cox_sd_metanalysis, model2_cox_sd_metanalysis, model3_cox_sd_metanalysis, # meta-anlaysis of all cohorts 
            model1_cox_quart_metanalysis, model2_cox_quart_metanalysis, model3_cox_quart_metanalysis) |> # meta-anlaysis of all cohorts 
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
  select(study, model, study_design, explanatory, term,  HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI, "p.value_heterogeneity")

heterogeneity_tests <- 
  bind_rows(heterogeneity_tests, heterogeneity_tests_finnish) |>
  mutate(study_design = ifelse(study == "Danish", NA, "raw unadjusted on cohort"))

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, heterogeneity_tests, 
            by = c("explanatory", "model", "study", "study_design")) |>
  mutate(p.value_heterogeneity = ifelse(is.na(p.value_heterogeneity.x), p.value_heterogeneity.y, p.value_heterogeneity.x), 
         p.value_heterogeneity = ifelse(term == "Continuous" & study == "Danish", NA, p.value_heterogeneity), 
         p.value_heterogeneity = ifelse(term == "Continuous" & study == "Finnish" & study_design == "raw unadjusted on cohort", NA, p.value_heterogeneity), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity)) |>
  select(-p.value_heterogeneity.x, -p.value_heterogeneity.y)

trend_tests <- bind_rows(trend_tests, trend_tests_finnish)|>
  mutate(study_design = ifelse(study == "Danish", NA, "raw unadjusted on cohort"))

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, trend_tests, 
            by = c("explanatory", "model", "study", "study_design")) |>
  mutate(p.value_trend = ifelse(term == "Continuous", NA, p.value_trend), 
         p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
         p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(model1_cox_sd_danish, 
   # model2_cox_sd_danish, 
   model3_cox_sd_danish, 
   model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish,
   model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish, 
   model1_cox_sd_finnish_raw, #model2_cox_sd_finnish_raw, 
   model3_cox_sd_finnish_raw, 
   model1_cox_quart_finnish_raw, model2_cox_quart_finnish_raw, model3_cox_quart_finnish_raw, 
   model1_cox_quart_finnish, model2_cox_quart_finnish, model3_cox_quart_finnish,  
   model1_cox_sd_metanalysis, model2_cox_sd_metanalysis, model3_cox_sd_metanalysis, 
   model1_cox_quart_metanalysis, model2_cox_quart_metanalysis, model3_cox_quart_metanalysis, 
   bdd_cases_FMC, bdd_cases_FMCF, bdd_cases_MFH, 
   heterogeneity_tests, trend_tests, 
   heterogeneity_tests_finnish, trend_tests_finnish)


# Sensitivity analysis 1 - effects of baseline age and duration baseline - diagnosis among all cohorts ----
## data prep ----
POPs_group_bis <- setdiff(POPs_group, "ΣPBDE")
POPs_group_quart_bis <- setdiff(POPs_group_quart, "ΣPBDE_quart")
POPs_group_sd_bis <- setdiff(POPs_group_sd, "ΣPBDE_sd")

bdd_cases_tot_sauf_MFH <- bdd |>
  filter (als == 1) |>
  filter(!study == 'MFH') |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, diagnosis_age, baseline_age, marital_status_2cat, smoking_2cat, bmi, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_bis),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

bdd_cases_red_purple <- bdd |>
  filter (als == 1) |>
  # filter(baseline_age>45) |>
  filter(follow_up<300) |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, diagnosis_age, marital_status_2cat, smoking_2cat, bmi, baseline_age, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_bis),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

bdd_cases_red_green <- bdd |>
  filter (als == 1) |>
  filter(baseline_age>45) |>
  # filter(follow_up<300) |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, diagnosis_age, marital_status_2cat, smoking_2cat, bmi, baseline_age, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_bis),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

bdd_cases_red_orange <- bdd |>
  filter (als == 1) |>
  filter(baseline_age>40) |>
  filter(follow_up<350) |>
  select(study, study_2cat, als, follow_up_death, status_death, sex, diagnosis_age, marital_status_2cat, smoking_2cat, bmi, baseline_age, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group_bis),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

surv_obj_tot_sauf_MFH <- Surv(time = bdd_cases_tot_sauf_MFH$follow_up_death,    # set the outcomes
                     event = bdd_cases_tot_sauf_MFH$status_death)

surv_obj_red_purple <- Surv(time = bdd_cases_red_purple$follow_up_death,        # set the outcomes
                            event = bdd_cases_red_purple$status_death)

surv_obj_red_green <- Surv(time = bdd_cases_red_green$follow_up_death,          # set the outcomes
                           event = bdd_cases_red_green$status_death)

surv_obj_red_orange <- Surv(time = bdd_cases_red_orange$follow_up_death,        # set the outcomes
                     event = bdd_cases_red_orange$status_death)

sample_size_check <- 
  tbl_merge(
    list(bdd_cases_tot |> select(study) |> tbl_summary(),
         bdd_cases_tot_sauf_MFH |> select(study) |> tbl_summary(),
         bdd_cases_red_purple |> select(study) |> tbl_summary(),
         bdd_cases_red_green |> select(study) |> tbl_summary(),
         bdd_cases_red_orange |> select(study) |> tbl_summary()), 
    tab_spanner = c('**All subjects**', 
                    '**All subjects except MFH**', 
                    '**follow up < 300 months**', 
                    '**baseline age > 45 years**', 
                    '**follow up < 350 months and baseline age > 40 years**'))

sample_size_check_2cat <- 
  tbl_merge(
    list(bdd_cases_tot |> select(study_2cat) |> tbl_summary(),
         bdd_cases_tot_sauf_MFH |> select(study_2cat) |> tbl_summary(),
         bdd_cases_red_purple |> select(study_2cat) |> tbl_summary(),
         bdd_cases_red_green |> select(study_2cat) |> tbl_summary(),
         bdd_cases_red_orange |> select(study_2cat) |> tbl_summary()), 
    tab_spanner = c('**All subjects**', 
                    '**All subjects except MFH**', 
                    '**follow up < 300 months**', 
                    '**baseline age > 45 years**', 
                    '**follow up < 350 months and baseline age > 40 years**'))

## justification ----
### 4 cohorts ----
plot_follow_up <- ggplot(bdd_cases_tot, aes(x = "", y = follow_up)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Follow-up (months)", title = "Distribution of the duration \nbetween baseline and diagnosis")

plot_baseline_age <- ggplot(bdd_cases_tot, aes(x = "", y = baseline_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Age at baseline (years)", title = "Distribution of the subjects \nage at baseline")

plot_diagnosis_age <- ggplot(bdd_cases_tot, aes(x = "", y = diagnosis_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Age at diagnosis (years)", title = "Distribution of the subjects \nage at diagnosis")

plot_baseline_follow_up <- 
  ggplot(bdd_cases_tot) +
  aes(x = baseline_age, y = follow_up, colour = study) +
  geom_point(alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_rect(aes(xmin = 45, xmax = 66, ymin = 0, ymax = 575),
            fill = NA, color = "chartreuse3", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 45.5, y = 575, label = "n=204",
           hjust = 0, vjust = -0.5, size = 5, color = "chartreuse3") +
  
  geom_rect(aes(xmin = 15, xmax = 66, ymin = 0, ymax = 300),
            fill = NA, color = "darkorchid", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 16, y = 300, label = "n=208",
           hjust = 0, vjust = -0.5, size = 5, color = "darkorchid") +
  
  geom_rect(aes(xmin = 40, xmax = 66, ymin = 0, ymax = 350),
            fill = NA, color = "darksalmon", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 53, y = 350, label = "n=207",
           hjust = 0, vjust = -0.5, size = 5, color = "darksalmon") +
  
  labs(x = "Age at baseline (years)", 
       y = "Duration between baseline and \ndiagnosis (months)")

plot_justif <- (plot_baseline_age + plot_follow_up + plot_diagnosis_age)/plot_baseline_follow_up + plot_layout()

### 2 cohorts ----
plot_follow_up_2cat <- ggplot(bdd_cases_tot, aes(x = "", y = follow_up)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study_2cat), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Follow-up (months)", 
       title = "Follow-up duration \nbetween baseline and diagnosis",
       colour = "Cohorts")

plot_baseline_age_2cat <- ggplot(bdd_cases_tot, aes(x = "", y = baseline_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study_2cat), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, 
       y = "Age at baseline (years)", 
       title = "Age at baseline",
       colour = "Cohorts")

plot_diagnosis_age_2cat <- ggplot(bdd_cases_tot, aes(x = "", y = diagnosis_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study_2cat), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, 
       y = "Age at diagnosis (years)", 
       title = "Age at diagnosis",
       colour = "Cohorts")

plot_baseline_follow_up_2cat <- 
  bdd_cases_tot |>
  ggplot() +
  aes(x = baseline_age, y = follow_up, colour = study_2cat) +
  geom_point(alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_rect(aes(xmin = 45, xmax = 66, ymin = 0, ymax = 575),
            fill = NA, color = "chartreuse3", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 45.5, y = 575, label = "Stratified to baseline age > 45 years",
           hjust = 0, vjust = -0.5, size = 4, color = "chartreuse3") +
  
  geom_rect(aes(xmin = 15, xmax = 66, ymin = 0, ymax = 300),
            fill = NA, color = "darkorchid", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 16, y = 300, label = "Stratified to follow-up < 300 months",
           hjust = 0, vjust = -0.5, size = 4, color = "darkorchid") +
  
  geom_rect(aes(xmin = 40, xmax = 66, ymin = 0, ymax = 350),
            fill = NA, color = "darksalmon", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 41, y = 350, label = "Stratified to baseline age > 40 years and follow-up < 350 months",
           hjust = 0, vjust = -0.5, size = 4, color = "darksalmon") +
  
  labs(x = "Age at baseline (years)", 
       y = "Follow-up (months)", 
       colour = "Cohorts")

plot_justif_2cat <- (plot_baseline_age_2cat + plot_follow_up_2cat + plot_diagnosis_age_2cat)/plot_baseline_follow_up_2cat + plot_layout()

cor_test_sensi1 <- 
  bdd_cases_tot |> 
  select(baseline_age, diagnosis_age, follow_up)
cor_test_sensi1 <- 
  cor(cor_test_sensi1, 
      use = "pairwise.complete.obs", 
      method = "pearson")

rm(plot_baseline_age, plot_follow_up, plot_diagnosis_age, plot_baseline_follow_up, 
   plot_baseline_age_2cat, plot_follow_up_2cat,  plot_diagnosis_age_2cat, plot_baseline_follow_up_2cat)
         

## analysis ----
### with adjustment for study (4 cat) ----
model2_cox_sd_sensi1 <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+ study +",              # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot) |> summary()            # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot", 
    model = "adjusted", 
    cohort_adj = "adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi1_sauf_MFH <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot_sauf_MFH ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot_sauf_MFH) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot_sauf_MFH", 
    model = "adjusted", 
    cohort_adj = "adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi1_red_purple <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_purple ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_purple) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_purple", 
    model = "adjusted", 
    cohort_adj = "adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi1_red_green <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_green ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_green) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_green", 
    model = "adjusted", 
    cohort_adj = "adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi1_red_orange <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_orange ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_orange) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_orange", 
    model = "adjusted", 
    cohort_adj = "adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### with adjustment for study (2 cat) ----
model2_cox_sd_sensi0 <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+ study_2cat +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot", 
    model = "adjusted", 
    cohort_adj = "adjusted (2 cat)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi0_sauf_MFH <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot_sauf_MFH ~", expl, "+ study_2cat +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot_sauf_MFH) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot_sauf_MFH", 
    model = "adjusted", 
    cohort_adj = "adjusted (2 cat)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi0_red_purple <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_purple ~", expl, "+ study_2cat +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_purple) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_purple", 
    model = "adjusted", 
    cohort_adj = "adjusted (2 cat)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi0_red_green <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_green ~", expl, "+ study_2cat +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_green) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_green", 
    model = "adjusted", 
    cohort_adj = "adjusted (2 cat)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi0_red_orange <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_orange ~", expl, "+ study_2cat +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_orange) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_orange", 
    model = "adjusted", 
    cohort_adj = "adjusted (2 cat)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})


### without adjustment for study ----
model2_cox_sd_sensi2 <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+",           # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot", 
    model = "adjusted", 
    cohort_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_sauf_MFH <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot_sauf_MFH ~", expl, "+ ",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot_sauf_MFH) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot_sauf_MFH", 
    model = "adjusted", 
    cohort_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_red_purple <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_purple ~", expl,  "+",          # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_purple) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_purple", 
    model = "adjusted", 
    cohort_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_red_green <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_green ~", expl,  "+",          # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_green) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_green", 
    model = "adjusted", 
    cohort_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_red_orange <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red_orange ~", expl,  "+",           # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red_orange) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red_orange", 
    model = "adjusted", 
    cohort_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

### with interaction on cohort ----
model2_cox_sd_sensi1_interac <- map_dfr(POPs_group_sd_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "* study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot_interaction", 
    model = "adjusted", 
    cohort_adj = "interaction",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})


### assemblage ----
results_sensi1 <-       
  bind_rows(
    model2_cox_sd_danish, # Danish main analysis
    model2_cox_sd_finnish_raw, 
    
    model2_cox_sd_sensi0, 
    model2_cox_sd_sensi0_sauf_MFH, 
    model2_cox_sd_sensi0_red_purple, 
    model2_cox_sd_sensi0_red_green, 
    model2_cox_sd_sensi0_red_orange, 
    
    model2_cox_sd_sensi1, 
    model2_cox_sd_sensi1_sauf_MFH, 
    model2_cox_sd_sensi1_red_purple, 
    model2_cox_sd_sensi1_red_green, 
    model2_cox_sd_sensi1_red_orange, 
    
    model2_cox_sd_sensi2, 
    model2_cox_sd_sensi2_sauf_MFH, 
    model2_cox_sd_sensi2_red_purple, 
    model2_cox_sd_sensi2_red_green, 
    model2_cox_sd_sensi2_red_orange, 
    model2_cox_sd_sensi1_interac) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, cohort_adj, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  mutate(
    term = case_when(
      str_detect(term, "FMC") ~ "FMC",
      str_detect(term, "FMCF") ~ "FMCF",
      str_detect(term, "MFH") ~ "MFH", 
      .default = "Danish"), 
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
  select(study, model, cohort_adj, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI)

POPs_sd_ALS_figure_sensi1_not_adjusted <- 
  results_sensi1 |>
  filter(!explanatory == 'ΣPBDE') |>
  filter(study %in% c("Danish", "Finnish", "sensi1_tot", "sensi1_red_orange")) |>
  mutate(
    study =  fct_relevel(study, 
                         "Danish", "Finnish", "sensi1_tot", "sensi1_red_orange"
                         # "sensi1_tot_sauf_MFH", "sensi1_red_purple", "sensi1_red_green", 
                         ), 
    study = fct_recode(study, 
                       "Danish (n=166)\n not adjusted on cohort" = "Danish",
                       "Finnish (n=97)\n not adjusted on cohort" = "Finnish",
                       "All cohorts (n=263)\n not adjusted on cohort" = "sensi1_tot", 
                       "Restricted to follow\nup < 350 months and \nbaseline age > 40 years \n(n=208)\n not adjusted on cohort" = "sensi1_red_orange"
                       # "Restricted to \nfollow up < 300 months \n(n=204)" = "sensi1_red_purple",
                       # "Restricted to \nbaseline age > 45 years \n(n=207)" = "sensi1_red_green",
                       # "All cohorts expect MFH (n=252)" = "sensi1_tot_sauf_MFH"
                       ),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
    cohort_adj = ifelse(is.na(cohort_adj), "not adjusted", cohort_adj), 
    cohort_adj = factor(cohort_adj, levels = c("adjusted", "not adjusted"))) |>
  filter(cohort_adj == "not adjusted") |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid(~study)

POPs_sd_ALS_figure_sensi1_adjusted <- 
  results_sensi1 |>
  filter(!explanatory == 'ΣPBDE') |>
  filter(study %in% c("sensi1_tot", "sensi1_red_orange")) |>
  mutate(
    study =  fct_relevel(study, 
                          "sensi1_tot", "sensi1_red_orange"
                         # "sensi1_tot_sauf_MFH", "sensi1_red_purple", "sensi1_red_green", 
    ), 
    study = fct_recode(study, 
                       "All cohorts (n=263)\n adjusted on cohort (4 cat)" = "sensi1_tot", 
                       "Restricted to \nfollow up < 350 months and \nbaseline age > 40 years \n(n=208)\n adjusted on cohort (4 cat)" = "sensi1_red_orange"
                       # "Restricted to \nfollow up < 300 months \n(n=204)" = "sensi1_red_purple",
                       # "Restricted to \nbaseline age > 45 years \n(n=207)" = "sensi1_red_green",
                       # "All cohorts expect MFH (n=252)" = "sensi1_tot_sauf_MFH"
    ),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
    cohort_adj = ifelse(is.na(cohort_adj), "adjusted", cohort_adj), 
    cohort_adj = factor(cohort_adj, levels = c("adjusted", "not adjusted"))) |>
  filter(cohort_adj == "adjusted") |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid(~study)

POPs_sd_ALS_figure_sensi1_adjusted_2cat <- 
  results_sensi1 |>
  filter(!explanatory == 'ΣPBDE') |>
  filter(study %in% c( "sensi1_tot", "sensi1_red_orange")) |>
  mutate(
    study =  fct_relevel(study, 
                        "sensi1_tot", "sensi1_red_orange"
                         # "sensi1_tot_sauf_MFH", "sensi1_red_purple", "sensi1_red_green", 
    ), 
    study = fct_recode(study, 
                       "All cohorts (n=263)\n adjusted on cohort (2 cat)" = "sensi1_tot", 
                       "Restricted to \nfollow up < 350 months and \nbaseline age > 40 years \n(n=208)\n adjusted on cohort (2 cat)" = "sensi1_red_orange"
                       # "Restricted to \nfollow up < 300 months \n(n=204)" = "sensi1_red_purple",
                       # "Restricted to \nbaseline age > 45 years \n(n=207)" = "sensi1_red_green",
                       # "All cohorts expect MFH (n=252)" = "sensi1_tot_sauf_MFH"
    ),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
    cohort_adj = ifelse(is.na(cohort_adj), "adjusted (2 cat)", cohort_adj), 
    cohort_adj = factor(cohort_adj, levels = c("adjusted (2 cat)", "not adjusted"))) |>
  filter(cohort_adj == "adjusted (2 cat)") |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid(~study)


rm(POPs_group_bis, POPs_group_quart_bis, POPs_group_sd_bis,  
   bdd_cases_red_purple, bdd_cases_red_green, bdd_cases_red_orange, bdd_cases_tot_sauf_MFH,
   surv_obj_tot, surv_obj_red_purple, surv_obj_red_green, surv_obj_red_orange, surv_obj_tot_sauf_MFH, 
   model2_cox_sd_sensi0, model2_cox_sd_sensi0_red_purple, model2_cox_sd_sensi0_red_green, model2_cox_sd_sensi0_red_orange, 
   model2_cox_sd_sensi1, model2_cox_sd_sensi1_red_purple, model2_cox_sd_sensi1_red_green, model2_cox_sd_sensi1_red_orange, 
   model2_cox_sd_sensi2, model2_cox_sd_sensi2_red_purple, model2_cox_sd_sensi2_red_green, model2_cox_sd_sensi2_red_orange, 
   model2_cox_sd_sensi0_sauf_MFH, model2_cox_sd_sensi1_interac, model2_cox_sd_sensi1_sauf_MFH, model2_cox_sd_sensi2_sauf_MFH)

POPs_group_bis <- setdiff(POPs_group_sd, c("PCB_4_sd", 'ΣPBDE_sd'))
pollutant_labels_bis <- set_names(
  POPs_group_bis,
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane"))

POPs_heatmap_cases_group <- 
  bdd_cases_danish |> 
  select(all_of(POPs_group_bis)) |>
  rename(!!!pollutant_labels_bis) 
POPs_heatmap_cases_group <- 
  cor(POPs_heatmap_cases_group, 
      use = "pairwise.complete.obs", 
      method = "pearson")
rm(POPs_group_bis, pollutant_labels_bis)


# Sensitivity analysis 2 - effects of baseline age and duration baseline - diagnosis just among the Finnish cohorts ----
# Testing the effects of age at baseline, diagnosis and duration between the two, just among the finnish cohorts
bdd_cases_finnish <- bdd_cases_finnish |>
  mutate(box_orange = ifelse(baseline_age > 40 &
                        follow_up < 350, "in", "out"), 
         box_purple = ifelse(follow_up<300, "in", "out"), 
         box_green = ifelse(baseline_age>45, "in", "out"), 
         box_orange = fct_relevel(box_orange, "out", "in"),
         box_purple = fct_relevel(box_purple, "out", "in"), 
         box_green = fct_relevel(box_green, "out", "in"))

bdd_cases_danish <- bdd_cases_danish |>                                         # all the danish cases are in the boxs
  mutate(box_orange = ifelse(baseline_age > 40 &
                               follow_up < 350, "in", "out"), 
         box_purple = ifelse(follow_up<300, "in", "out"), 
         box_green = ifelse(baseline_age>45, "in", "out"), 
         box_orange = fct_relevel(box_orange, "out", "in"),
         box_purple = fct_relevel(box_purple, "out", "in"), 
         box_green = fct_relevel(box_green, "out", "in"))

tbl_distrib_among_green_box <- tbl_merge(
  tbls = list(bdd_cases_finnish |>
                filter(box_green == "in") |>
                select(baseline_age) |>
                tbl_summary(), 
              bdd_cases_danish |>
                filter(box_green == "in") |>
                select(baseline_age) |>
                tbl_summary()), 
  tab_spanner = c("**Finnish cases in the green box (baseline age > 45 years)**", 
                  "**Danish cases among the green box (baseline age > 45 years)**"))

tbl_distrib_among_purple_box <- tbl_merge(
  tbls = list(bdd_cases_finnish |>
                filter(box_purple == "in") |>
                select(follow_up) |>
                tbl_summary(), 
              bdd_cases_danish |>
                filter(box_purple == "in") |>
                select(follow_up) |>
                tbl_summary()),
  tab_spanner = c("**Finnish cases in the purple box (follow-up duration < 300 months)**", 
                  "**Danish cases in the purple box (follow-up duration < 300 months)**"))

tbl_distrib_among_orange_box <- tbl_merge(
  tbls = list(
    tbl_stack(
      tbls = list(bdd_cases_finnish |>
                    filter(box_orange == "in") |>
                    select(baseline_age) |>
                    tbl_summary(), 
                  bdd_cases_finnish |>
                    filter(box_orange == "in") |>
                    select(follow_up) |>
                    tbl_summary())), 
    bdd_cases_danish |>
      filter(box_orange == "in") |>
      select(baseline_age, follow_up) |>
      tbl_summary()), 
  tab_spanner = c("**Finnish cases in the orange box (baseline age > 40 and follow-up duration < 350 months)**", 
                  "**Danish cases in the orange box (baseline age > 40 and follow-up duration < 350 months)**"))


model2_cox_sd_sensi2_finnish <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ ",                 # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_finnish) |> summary()        # run cox model
  coefs <- model_summary$coefficients
  tibble(
    # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "not adjusted",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

## orange box ----
model2_cox_sd_sensi2_finnish_adj_orange <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ box_orange +",                    # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_finnish) |> summary()        # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "adj_orange",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_in_orange <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_in <- bdd_cases_finnish |> filter(box_orange == "in")
  
  surv_obj <- Surv(time = bdd_cases_finnish_in$follow_up_death,                 # set the outcomes
                   event = bdd_cases_finnish_in$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_in) |> summary()     # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "in_orange",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_out_orange <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_out <- bdd_cases_finnish |> filter(box_orange == "out")
  
  surv_obj <- Surv(time = bdd_cases_finnish_out$follow_up_death,                # set the outcomes
                   event = bdd_cases_finnish_out$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_out) |> summary()    # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "out_orange",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"] ) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

## purple box ----
model2_cox_sd_sensi2_finnish_adj_purple <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ box_purple +",                    # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_finnish) |> summary()        # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "adj_purple",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_in_purple <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_in <- bdd_cases_finnish |> filter(box_purple == "in")
  
  surv_obj <- Surv(time = bdd_cases_finnish_in$follow_up_death,                 # set the outcomes
                   event = bdd_cases_finnish_in$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_in) |> summary()     # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "in_purple",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_out_purple <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_out <- bdd_cases_finnish |> filter(box_purple == "out")
  
  surv_obj <- Surv(time = bdd_cases_finnish_out$follow_up_death,                # set the outcomes
                   event = bdd_cases_finnish_out$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_out) |> summary()    # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "out_purple",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"] ) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

## green box ----
model2_cox_sd_sensi2_finnish_adj_green <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj_finnish ~", expl, "+ box_green +",                    # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_finnish) |> summary()        # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "adj_green",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_in_green <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_in <- bdd_cases_finnish |> filter(box_green == "in")
  
  surv_obj <- Surv(time = bdd_cases_finnish_in$follow_up_death,                 # set the outcomes
                   event = bdd_cases_finnish_in$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_in) |> summary()     # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "in_green",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_sd_sensi2_finnish_out_green <- map_dfr(POPs_group_sd_finnish, function(expl) {
  formula <- as.formula(paste("surv_obj ~", expl, "+",                          # set the formulas
                              paste(covariates_finnish, collapse = " + ")))
  
  bdd_cases_finnish_out <- bdd_cases_finnish |> filter(box_green == "out")
  
  surv_obj <- Surv(time = bdd_cases_finnish_out$follow_up_death,                # set the outcomes
                   event = bdd_cases_finnish_out$status_death)
  
  model_summary <- coxph(formula, data = bdd_cases_finnish_out) |> summary()    # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi2_finnish",
    model = "adjusted",
    box_adj = "out_green",
    term = rownames(coefs),
    explanatory = expl,
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"],
    `p-value` = coefs[, "Pr(>|z|)"] ) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

results_sensi2 <-
  bind_rows(
    model2_cox_sd_sensi2_finnish,
    model2_cox_sd_sensi2_finnish_adj_orange,
    model2_cox_sd_sensi2_finnish_in_orange,
    model2_cox_sd_sensi2_finnish_out_orange,
    model2_cox_sd_sensi2_finnish_adj_purple,
    model2_cox_sd_sensi2_finnish_in_purple,
    model2_cox_sd_sensi2_finnish_out_purple,
    model2_cox_sd_sensi2_finnish_adj_green,
    model2_cox_sd_sensi2_finnish_in_green,
    model2_cox_sd_sensi2_finnish_out_green) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    explanatory = gsub("_sd", "", explanatory),
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`,
    `p-value_shape` = ifelse(`p-value_raw` < 0.05, "p-value<0.05", "p-value≥0.05"),
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(
    study, model, box_adj, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`,
    lower_CI,  upper_CI)

rm(model2_cox_sd_sensi2_finnish,
  model2_cox_sd_sensi2_finnish_adj_orange,
  model2_cox_sd_sensi2_finnish_in_orange,
  model2_cox_sd_sensi2_finnish_out_orange,
  model2_cox_sd_sensi2_finnish_adj_purple,
  model2_cox_sd_sensi2_finnish_in_purple,
  model2_cox_sd_sensi2_finnish_out_purple,
  model2_cox_sd_sensi2_finnish_adj_green,
  model2_cox_sd_sensi2_finnish_in_green,
  model2_cox_sd_sensi2_finnish_out_green)

POPs_sd_ALS_figure_sensi2_finnish_orange <-
  results_sensi2 |>
  filter(box_adj %in% c("not adjusted", "adj_orange", "in_orange", "out_orange")) |>
  mutate(
    box_adj =  fct_relevel(box_adj, "not adjusted", "adj_orange", "in_orange", "out_orange"),
    box_adj = fct_recode(
      box_adj,
      "Not adjusted on baseline age>40\n& follow up<350 (Y/N, n=97)\nMain analysis" = "not adjusted",
      "Adjusted on baseline age>40 &\nfollow up<350 (Y/N, n=97)" = "adj_orange",
      "Stratified to baseline age>40 &\nfollow up<350 (n=41)" = "in_orange",
      "Stratified to baseline age<40 &\nfollow up>350 (n=56)" = "out_orange"),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid( ~ box_adj)
  
POPs_sd_ALS_figure_sensi2_finnish_purple <-
  results_sensi2 |>
  filter(box_adj %in% c("not adjusted", "adj_purple", "in_purple", "out_purple")) |>
  mutate(
    box_adj =  fct_relevel(box_adj, "not adjusted", "adj_purple", "in_purple", "out_purple"),
    box_adj = fct_recode(
      box_adj,
      "Not adjusted on follow up<300\n(Y/N, n=97)\nMain analysis" = "not adjusted",
      "Adjusted on follow up<300\n(Y/N, n=97)" = "adj_purple",
      "Stratified to follow up<300 (n=42)" = "in_purple",
      "Stratified to follow up>300 (n=55)" = "out_purple"),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid( ~ box_adj)

POPs_sd_ALS_figure_sensi2_finnish_green <-
  results_sensi2 |>
  filter(box_adj %in% c("not adjusted", "adj_green", "in_green", "out_green")) |>
  mutate(
    box_adj =  fct_relevel(box_adj, "not adjusted", "adj_green", "in_green", "out_green"),
    box_adj = fct_recode(
      box_adj,
      "Not adjusted on baseline\nage>45 (Y/N, n=97)\nMain analysis" = "not adjusted",
      "Adjusted on baseline\nage>45 (Y/N, n=97)" = "adj_green",
      "Stratified to baseline age>45\n(n=38)" = "in_green",
      "Stratified to baseline age>45\n(n=59)" = "out_green"),
    explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid( ~ box_adj)

# Sensitivity analysis 3 - age at baseline as interaction term in the finnish cohort ----
sensi3_interac_baseline_age <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <-                                                             # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "*baseline_age +",  
                     paste(covariates_finnish, collapse = " + ")))
  
  
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "interac baseline age continuous", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory) | term == "baseline_age")              # remove the covariates results 
})

sensi3_interac_follow_up <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  bdd_cases_finnish <- bdd_cases_finnish |>
    mutate(follow_up_years = follow_up / 12)
  
  formula_finnish <-                                                             # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "*follow_up_years +",  
                     paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "interac follow up continuous", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory) | term == "follow_up")                 # remove the covariates results 
})

sensi3_interac_box_orange <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <-                                                             # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "*box_orange +",  
                     paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "interac orange", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory) | term == "box_orangeout")             # remove the covariates results 
})

sensi3_interac_box_purple <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <-                                                             # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "*box_purple +",  
                     paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "interac purple", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory) | term == "box_purpleout")             # remove the covariates results 
})

sensi3_interac_box_green <- map_dfr(POPs_group_sd_finnish, function(expl) {
  
  formula_finnish <-                                                             # set the formulas
    as.formula(paste("surv_obj_finnish ~", expl, "*box_green +",  
                     paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula_finnish, data = bdd_cases_finnish) |> summary()# run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Finnish", 
    model = "adjusted", 
    study_design = "interac green", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory) | term == "box_greenout")              # remove the covariates results 
})

results_sensi3 <-
  bind_rows(
    sensi3_interac_baseline_age,
    sensi3_interac_follow_up,
    sensi3_interac_box_orange,
    sensi3_interac_box_purple,
    sensi3_interac_box_green) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    explanatory = gsub("_sd", "", explanatory),
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`,
    `p-value_shape` = ifelse(`p-value_raw` < 0.05, "p-value<0.05", "p-value≥0.05"),
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(
    study, model, study_design, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`,
    lower_CI,  upper_CI)

rm(sensi3_interac_baseline_age,
  sensi3_interac_follow_up,
  sensi3_interac_box_orange,
  sensi3_interac_box_purple,
  sensi3_interac_box_green)

POPs_sd_ALS_table_sensi3 <- 
  results_sensi3 |>
  filter(study == "Finnish") |>
  select(study_design, explanatory, term, HR, "95% CI", "p-value_raw") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
         term = term |>
           str_replace_all(c(
             "baseline_age" = "Baseline age",
             "follow_up" = "Follow-up",
             "box_orangein" = "In of the orange box",
             "box_purplein" = "In of the purple box",
             "box_greenin" = "In of the green box", 
             "_sd" = "")) |>
           str_replace_all(setNames(names(POPs_group_labels_finnish), POPs_group_labels_finnish)), 
         study_design = fct_recode(study_design, 
                                   "Interaction with baseline age (continuous)" = "interac baseline age continuous",
                                   "Interaction with follow up (continuous)" = "interac follow up continuous",
                                   "Interaction with baseline age (categorical)" = "interac green",
                                   "Interaction with both baseline age and follow-up duration (categorical)" = "interac orange",
                                   "Interaction with follow-up (categorical)" = "interac purple"), 
         `p-value` = if_else(sprintf("%.2f", `p-value_raw`) == "0.00", "<0.01", sprintf("%.2f", `p-value_raw`))) 

POPs_sd_ALS_table_sensi3_baseline_age <- POPs_sd_ALS_table_sensi3 |>
  filter(study_design %in% c("Interaction with baseline age (continuous)", 
                             "Interaction with baseline age (categorical)")) |>
  flextable(col_keys = c("study_design", "explanatory", "term",  "HR", "95% CI", "p-value")) |>
  add_footer_lines(
    "1All models are adjusted for age at diagnosis, sex, smoking, BMI and marital status. 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "study_design" = "Study design",
    "explanatory" = "Exposures", 
    "term" = "Term") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  merge_v(j = "study_design") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  bold(j = "study_design", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "term", part = "header") |>
  merge_at(j = "study_design", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")|> 
  color(i = ~ `p-value_raw` < 0.05, color = "red", part = "body")

POPs_sd_ALS_table_sensi3_follow_up <- POPs_sd_ALS_table_sensi3 |>
  filter(study_design %in% c("Interaction with follow up (continuous)", 
                             "Interaction with follow-up (categorical)")) |>
  flextable(col_keys = c("study_design", "explanatory", "term",  "HR", "95% CI", "p-value")) |>
  add_footer_lines(
    "1All models are adjusted for age at diagnosis, sex, smoking, BMI and marital status. 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "study_design" = "Study design",
    "explanatory" = "Exposures", 
    "term" = "Term") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  merge_v(j = "study_design") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  bold(j = "study_design", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "term", part = "header") |>
  merge_at(j = "study_design", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")|> 
  color(i = ~ `p-value_raw` < 0.05, color = "red", part = "body")


POPs_sd_ALS_table_sensi3_baseline_age_follow_up <- POPs_sd_ALS_table_sensi3 |>
  filter(study_design == "Interaction with both baseline age and follow-up duration (categorical)") |>
  flextable(col_keys = c("study_design", "explanatory", "term",  "HR", "95% CI", "p-value")) |>
  add_footer_lines(
    "1All models are adjusted for age at diagnosis, sex, smoking, BMI and marital status. 
  2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease plasma concentration of POPs.
  3CI: Confidence interval.") |>
  add_header(
    "study_design" = "Study design",
    "explanatory" = "Exposures", 
    "term" = "Term") |>
  merge_h(part = "header") |>
  merge_v(j = "explanatory") |>
  merge_v(j = "study_design") |>
  theme_vanilla() |>
  bold(j = "explanatory", part = "body") |>
  bold(j = "study_design", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_at(j = "term", part = "header") |>
  merge_at(j = "study_design", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")|> 
  color(i = ~ `p-value_raw` < 0.05, color = "red", part = "body")
rm(POPs_sd_ALS_table_sensi3)

# Sensitivity analysis 4 - effects of survival duration among the Danish population ----
## justification ----
bdd_cases_tot |>
  select(study_2cat, follow_up_death) |>
  tbl_summary(by = study_2cat)

plot_follow_up_death <- ggplot(bdd_cases_tot, aes(x = "", y = follow_up_death)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Follow-up from diagnosis to death or end of study (months)", title = "Distribution of the duration\nbetween diagnosis and death\nor end fo the study")

plot_diagnosis_age <- ggplot(bdd_cases_tot, aes(x = "", y = diagnosis_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Age at diagnosis (years)", title = "Distribution of the subjects \nage at diagnosis")

plot_diagnosis_follow_up_death <- 
  ggplot(bdd_cases_tot) +
  aes(x = diagnosis_age, y = follow_up_death, colour = study) +
  geom_point(alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_rect(aes(xmin = 43, xmax = 93, ymin = 0, ymax = 24),
            fill = NA, color = "darkorchid", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 45, y = 24, label = "n=85 (51%)",
           hjust = 0, vjust = -0.5, size = 5, color = "darkorchid") +
  
  labs(x = "Age at diagnosis (years)", 
       y = "Duration between diagnosis and \ndeath or end of study (months)")

plot_justif_sensi4 <- (plot_follow_up_death + plot_diagnosis_age)/plot_diagnosis_follow_up_death

rm(plot_follow_up_death, plot_diagnosis_age, plot_diagnosis_follow_up_death)

bdd_cases_danish <- 
  bdd_cases_danish |>
  mutate(follow_up_death_cat = ifelse(follow_up_death>23, "survival ≥ 24 months", "survival <  24 months"))

bdd_cases_danish |> select(follow_up_death_cat) |> tbl_summary()

## analysis ----
model2_cox_sd_danish_sensi4_in <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- 
    bdd_cases_danish |> 
    mutate(follow_up_death_cat = ifelse(follow_up_death>23, "survival ≥ 24 months", "survival <  24 months")) |>
    filter(follow_up_death_cat == "survival <  24 months")
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,              # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+",  
                     paste(covariates_danish, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "adjusted", 
    study_design = "stratified to survival duration\n<  24 months (n=85, 51%)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})


model2_cox_sd_danish_sensi4_out <- map_dfr(POPs_group_sd, function(expl) {
  
  bdd_cases_danish <- 
    bdd_cases_danish |> 
    mutate(follow_up_death_cat = ifelse(follow_up_death>23, "survival ≥ 24 months", "survival <  24 months")) |>
    filter(follow_up_death_cat == "survival ≥ 24 months")
  
  surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                          event = bdd_cases_danish$status_death)
  
  formula_danish <-                                                             # set the formulas
    as.formula(paste("surv_obj_danish ~", expl, "+",  
                     paste(covariates_danish, collapse = " + ")))
  
  model_summary <- coxph(formula_danish, data = bdd_cases_danish) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "Danish", 
    model = "adjusted", 
    study_design = "stratified to survival duration\n≥ 24 months (n=81, 49%)",
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results 
})


results_sensi4 <-
  bind_rows(
    model2_cox_sd_danish_sensi4_in,
    model2_cox_sd_danish_sensi4_out) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    explanatory = gsub("_sd", "", explanatory),
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`,
    `p-value_shape` = ifelse(`p-value_raw` < 0.05, "p-value<0.05", "p-value≥0.05"),
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(
    study, model, study_design, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`,
    lower_CI,  upper_CI)

rm(model2_cox_sd_danish_sensi4_in,
  model2_cox_sd_danish_sensi4_out)

## figure ----
POPs_sd_ALS_figure_sensi4_danish <-
  results_sensi4 |>
  mutate(
    explanatory = factor(explanatory, levels = POPs_group_labels),
    explanatory = fct_rev(explanatory),
    explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() +
  facet_grid( ~ study_design)


# Sensitivity analysis 5 - mixture model in the Danish EPIC cohort ----
## Data prep ----
bdd_cases_danish_bis <- bdd_danish |>
  filter (als == 1) |>
  filter(follow_up_death>0) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
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
  replace_with_median(ΣPBDE, ΣPBDE_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))

POPs_group_sd_bis <- setdiff(POPs_group_sd, "PCB_4_sd")    
POPs_group_quart_bis <- setdiff(POPs_group_quart, "PCB_4_quart")    

X_matrix_sd <- model.matrix(~ ., data = bdd_cases_danish_bis[, c(covariates_danish, POPs_group_sd_bis)])[, -1] 
X_matrix_quart <- model.matrix(~ ., data = bdd_cases_danish_bis[, c(covariates_danish, POPs_group_quart_bis)])[, -1] 

penalty_factor_sd <-                                                            # création d'un penalty factor pour forcer le modele lasso à inclure les covariables dans le modele 
  ifelse(colnames(X_matrix_sd) %in% 
           colnames(model.matrix(~ ., data = bdd_cases_danish_bis[, covariates_danish])[, -1]), 
         0, 1)

penalty_factor_quart <-                                                         # création d'un penalty factor pour forcer le modele lasso à inclure les covariables dans le modele 
  ifelse(colnames(X_matrix_quart) %in% 
           colnames(model.matrix(~ ., data = bdd_cases_danish_bis[, covariates_danish])[, -1]), 
         0, 1)


## LASSO POPs selection ----
set.seed(1996)
sensi5_lasso_sd_danish <- cv.glmnet(                                             # lasso + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_sd)                                           # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi5_lasso_sd_danish)
coef(sensi5_lasso_sd_danish, s = "lambda.min")                                   # lasso keeps only chlordane 
coef(sensi5_lasso_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi5_lasso_quart_danish <- cv.glmnet(                                          # lasso + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_quart)                                        # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi5_lasso_quart_danish)
coef(sensi5_lasso_quart_danish, s = "lambda.min")                                # lasso doesn't keep any pollutant when they are quartiles
coef(sensi5_lasso_quart_danish, s = "lambda.1se")  

## Ridge POPs selection ----
set.seed(1996)
sensi5_ridge_sd_danish <- cv.glmnet(                                             # ridge + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                                             

plot(sensi5_ridge_sd_danish)
coef(sensi5_ridge_sd_danish, s = "lambda.min")  
coef(sensi5_ridge_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi5_ridge_quart_danish <- cv.glmnet(                                          # ridge + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                                             

plot(sensi5_ridge_quart_danish)
coef(sensi5_ridge_quart_danish, s = "lambda.min")  
coef(sensi5_ridge_quart_danish, s = "lambda.1se")  

## Elastic net selection -----
### choose the best alpha ----
alpha_grid <- seq(0.3, 1, by = 0.1)  # de 0 (Ridge) à 1 (Lasso)                 # starting at least at 0.3 because we want to select variables, not just run a ridge 

cv_results <- list()
cvm_min <- numeric(length(alpha_grid))

for (i in seq_along(alpha_grid)) {   # Boucle pour trouver le meilleur alpha
  set.seed(1996)
  cv_model <- cv.glmnet(
    x = X_matrix_sd,
    y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),
    family = "cox",
    alpha = alpha_grid[i],
    type.measure = "deviance",
    nfolds = 10,
    standardize = FALSE,
    penalty.factor = penalty_factor_sd
  )
  cv_results[[i]] <- cv_model
  cvm_min[i] <- min(cv_model$cvm)  # stocke la deviance minimale
}

best_alpha <- alpha_grid[which.min(cvm_min)]
best_alpha

rm(alpha_grid, cv_results, cvm_min, i, cv_model, best_alpha)

### analysis ----
set.seed(1996)
sensi5_elastic_net_sd_danish <- cv.glmnet(                                       # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.4,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                          

plot(sensi5_elastic_net_sd_danish)
coef(sensi5_elastic_net_sd_danish, s = "lambda.min")                             # Elastic net keeps HCB and Σchlordane
coef(sensi5_elastic_net_sd_danish, s = "lambda.1se")

set.seed(1996)
sensi5_elastic_net_quart_danish <- cv.glmnet(                                    # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.3,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                          

plot(sensi5_elastic_net_quart_danish)
coef(sensi5_elastic_net_quart_danish, s = "lambda.min")                          # elastic net doesn't keep any pollutant when they are quartiles
coef(sensi5_elastic_net_quart_danish, s = "lambda.1se")
rm(bdd_cases_danish_bis)

### new copollutant model with selected POPs ----
#### sd ----
POPs_sd_selected <-                                                             # selected POPs by LASSO
  c("OCP_HCB_sd", "Σchlordane_sd")                 
POPs_sd_selected_labels <- 
  set_names(c("OCP_HCB", "Σchlordane"), c("HCB", "Σchlordane"))

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(POPs_sd_selected, covariates_danish), collapse = " + ")))

model_summary <- 
  coxph(formula_danish, data = bdd_cases_danish) |> summary() 
coefs <- model_summary$coefficients
sensi5_model3_cox_sd_elastic_net_danish <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "copollutant_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

#### quartiles ----
POPs_quart_selected <-                                                          # selected POPs by LASSO
  c("OCP_HCB_quart", "Σchlordane_quart")                 
POPs_quart_selected_labels <- 
  set_names(c("OCP_HCB", "Σchlordane"), c("HCB", "Σchlordane"))

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(POPs_quart_selected, covariates_danish), collapse = " + ")))

model_summary <- 
  coxph(formula_danish, data = bdd_cases_danish) |> summary() 
coefs <- model_summary$coefficients
sensi5_model3_cox_quart_elastic_net_danish <- tibble(                                  # creation of a table of results
  study = "Danish", 
  model = "copollutant_quart", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_quart"))


##### heterogeneity tests ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_HCB_full <- 
  gam(outcome ~ 
        OCP_HCB_quart + s(Σchlordane) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_HCB_raw <- 
  gam(outcome ~ 
        #OCP_HCB_quart + 
        s(Σchlordane) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_HCB_raw, model3_quart_HCB_full, test = "Chisq")
p.value_heterogeneity_HCB <- tibble(explanatory = "OCP_HCB_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_Σchlordane_full <- 
  gam(outcome ~ 
        Σchlordane_quart + 
        s(OCP_HCB) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_Σchlordane_raw <- 
  gam(outcome ~ 
        #Σchlordane_quart + 
        s(OCP_HCB) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_Σchlordane_raw, model3_quart_Σchlordane_full, test = "Chisq")
p.value_heterogeneity_Σchlordane <- tibble(explanatory = "Σchlordane_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

heterogeneity_copollutant_quart <- 
  bind_rows(p.value_heterogeneity_HCB, 
            p.value_heterogeneity_Σchlordane) |>
  mutate(model = "copollutant_quart", 
         explanatory = gsub("_quart", "", explanatory), 
         study = "Danish")

rm(anova, outcome, 
   model3_quart_HCB_full, model3_quart_HCB_raw, 
   model3_quart_Σchlordane_full, model3_quart_Σchlordane_raw, 
   p.value_heterogeneity_HCB, 
   p.value_heterogeneity_Σchlordane)

##### trend tests ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_HCB_trend <- 
  gam(outcome ~ 
        OCP_HCB_quart_med + 
         s(Σchlordane) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)|> 
  summary()
p.value_trend_HCB <- model3_quart_HCB_trend$p.table["OCP_HCB_quart_med", "Pr(>|z|)"]

model3_quart_Σchlordane_trend <- 
  gam(outcome ~ 
        Σchlordane_quart_med + 
        s(OCP_HCB) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_Σchlordane <- model3_quart_Σchlordane_trend$p.table["Σchlordane_quart_med", "Pr(>|z|)"]

trend_copollutant <- 
  data.frame(explanatory = c("OCP_HCB_quart_med", 
                             "Σchlordane_quart_med"),
             model = "copollutant_quart",
             p.value_trend = c(p.value_trend_HCB, 
                               p.value_trend_Σchlordane)) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory), 
         study = "Danish")


### ERS model avec les POP selectionnés par elastic net 0.4 ----
coef(sensi5_elastic_net_sd_danish, s = "lambda.min")                            # Elastic net keeps HCB and Σchlordane
cox_model <-                                                                    # regular co-pollutant cox model to get unpenalized coefficients 
  coxph(Surv(follow_up_death, status_death) ~ 
          OCP_HCB_sd + Σchlordane_sd + 
          sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
        data = bdd_cases_danish)  |> summary()

weight_HCB <- cox_model$coefficients["OCP_HCB_sd", "coef"] / (cox_model$coefficients["OCP_HCB_sd", "coef"]  + cox_model$coefficients["Σchlordane_sd", "coef"])
weight_Σchlordane <- cox_model$coefficients["Σchlordane_sd", "coef"] / (cox_model$coefficients["OCP_HCB_sd", "coef"]  + cox_model$coefficients["Σchlordane_sd", "coef"])

bdd_cases_danish <-
  bdd_cases_danish |>
  mutate(
    ERS_score_from_elastic_net = 
      OCP_HCB_sd * weight_HCB +                                                 # somme des 2 POPs qui ont ete multiplié par leur poids respectif
      Σchlordane_sd * weight_Σchlordane,   
    
    ERS_score_from_elastic_net_sd = 
      scale(ERS_score_from_elastic_net), 
    
    ERS_score_from_elastic_net_quart = 
      factor(ntile(ERS_score_from_elastic_net, 4),
             labels = c("Q1", "Q2", "Q3", "Q4")))


sensi5_cox_model_ERS_from_elastic_net_sd <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_sd + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi5_cox_model_ERS_from_elastic_net_sd)
coefs <- model_summary$coefficients
sensi5_cox_model_ERS_from_elastic_net_sd <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "ERS_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

sensi5_cox_model_ERS_from_elastic_net_quart <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_quart + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi5_cox_model_ERS_from_elastic_net_quart)
coefs <- model_summary$coefficients
sensi5_cox_model_ERS_from_elastic_net_quart <- tibble(                                  # creation of a table of results
  study = "Danish", 
  model = "ERS_quart", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_quart"))

rm(cox_model, coefs, model_summary, weight_HCB, weight_Σchlordane)


### assemblage ----
results_sensi5 <-
  bind_rows(sensi5_model3_cox_sd_elastic_net_danish, 
            sensi5_model3_cox_quart_elastic_net_danish, 
            sensi5_cox_model_ERS_from_elastic_net_sd, 
            sensi5_cox_model_ERS_from_elastic_net_quart) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    explanatory = gsub("_quartQ2", "", explanatory),
    explanatory = gsub("_quartQ3", "", explanatory),
    explanatory = gsub("_quartQ4", "", explanatory),
    explanatory = gsub("_sd", "", explanatory),
    term = case_when(
      str_detect(term, "_sd") ~ "Continuous", 
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"), 
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`,
    `p-value_shape` = ifelse(`p-value_raw` < 0.05, "p-value<0.05", "p-value≥0.05"),
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(
    study, model, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`,
    lower_CI,  upper_CI)

results_sensi5 <- left_join(results_sensi5, heterogeneity_copollutant_quart, 
                            by = c("study", "model", "explanatory"))

results_sensi5 <- 
  left_join(results_sensi5, trend_copollutant, 
            by = c("study", "model", "explanatory")) 

results_sensi5 <- results_sensi5 |>
  mutate(
    p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity), 
    p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
    p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(outcome, 
   model3_quart_HCB_trend, 
   model3_quart_Σchlordane_trend, 
   
   p.value_trend_HCB, 
   p.value_trend_Σchlordane, 
   
   heterogeneity_copollutant_quart, 
   trend_copollutant, 
   
   sensi5_cox_model_ERS_from_elastic_net_sd, 
   sensi5_cox_model_ERS_from_elastic_net_quart)

### table sd copollutant ----
POPs_sd_ALS_table_sensi5_danish <- results_sensi5 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_sd") |>
  select(explanatory, HR, "95% CI", "p-value") |>
  mutate(
    explanatory = factor(explanatory, levels = POPs_group_labels), 
    explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1Σchlordane corresponds to trans-nonanchlor and oxychlordane.
     2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
     3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Copollutant model", "95% CI" = "Copollutant model", "p-value" = "Copollutant model") |>
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

### table quart copollutant ----
quartile1_rows <- 
  results_sensi5 |>
  filter(model == "copollutant_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '') |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") 

POPs_quart_ALS_table_sensi5_danish <- 
  results_sensi5 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi5_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi5_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  group_by(explanatory) |>
  mutate(
    `Heterogeneity test` = ifelse(term == "quartile 1", first(`Heterogeneity test`[term == "quartile 2"]), ""),
    `Trend test` = ifelse(term == "quartile 1", first(`Trend test`[term == "quartile 2"]), "")) |>
  ungroup() |>
  mutate(explanatory = factor(explanatory, levels = POPs_quart_selected_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_quart_selected_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.
    5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", term = "Quartiles",
    "HR" = "Copollutant model", "95% CI" = "Copollutant model", "p-value" = "Copollutant model",  
    "Heterogeneity test" = "Copollutant model",  "Trend test" = "Copollutant model") |>
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


### figure sd copollutant ----
POPs_sd_ALS_figure_sensi5_danish <-
  results_sensi5 |>
  filter(model == "copollutant_sd") |>
  mutate(
    # explanatory = factor(explanatory, levels = POPs_sd_selected_labels),
    explanatory = fct_recode(explanatory, !!!POPs_sd_selected_labels), 
    explanatory = fct_rev(explanatory)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() 

### figure quart copollutant ----
POPs_quart_ALS_figure_sensi5_danish <-
  results_sensi5 |>
  filter(model == "copollutant_quart") |>
  mutate(
    # explanatory = factor(explanatory, levels = POPs_sd_selected_labels),
    explanatory = fct_recode(explanatory, !!!POPs_quart_selected_labels), 
    term = str_replace(term, "quartile", "Quartile"), 
    term = fct_rev(term)
    ) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = term,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  facet_grid(rows = dplyr::vars(explanatory), switch = "y") +  
  coord_flip() 

### table sd ERS ----
POPs_sd_ALS_table_sensi5_ERS_danish <- results_sensi5 |>
  filter(study == "Danish") |>
  filter(model == "ERS_sd") |>
  select(explanatory, HR, "95% CI", "p-value") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net")) |> 
  flextable() |>
  add_footer_lines(
    "1Environmental risk score is the weighted sum of relevant pollutants (HCB and Σchlordane) selected with elastic net regularization.
     2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
     3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase of the environmental risk score.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Mixture model", "95% CI" = "Mixture model", "p-value" = "Mixture model") |>
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

### table quart ERS ----
quartile1_rows <- 
  results_sensi5 |>
  filter(model == "ERS_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") 

POPs_quart_ALS_table_sensi5_ERS_danish <- 
  results_sensi5 |>
  filter(study == "Danish") |>
  filter(model == "ERS_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi5_ERS_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi5_ERS_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99"), 
         explanatory = fct_recode(explanatory, 
                                  "Environmental risk score" = "ERS_score_from_elastic_net")) |>
  arrange(explanatory, term) |>
  flextable() |>
  add_footer_lines(
    "1Environmental risk score is the weighted sum of relevant pollutants (HCB and Σchlordane) selected with elastic net regularization.
    2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    3Estimated risk of ALS death when the ERS is at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", term = "Quartiles",
    "HR" = "Mixture model", "95% CI" = "Mixture model", "p-value" = "Mixture model") |>
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

### figure sd ERS ----
POPs_sd_ALS_figure_sensi5_ERS_danish <-
  results_sensi5 |>
  filter(model == "ERS_sd") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net")) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() 


### figure quart ERS ----
POPs_quart_ALS_figure_sensi5_ERS_danish <-
  results_sensi5 |>
  filter(model == "ERS_quart") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net"),
    term = str_replace(term, "quartile", "Quartile"), 
    term = fct_rev(term)) |>
  ggplot(aes(
    x = term,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  facet_grid(rows = dplyr::vars(explanatory), switch = "y") +  
  coord_flip() 



rm(POPs_group_sd_bis, POPs_group_quart_bis, 
   X_matrix_sd, X_matrix_quart, 
   penalty_factor_sd, penalty_factor_quart, 
   POPs_sd_selected, POPs_sd_selected_labels, 
   POPs_quart_selected, POPs_quart_selected_labels,
   formula_danish, model_summary, coefs, 
   sensi5_lasso_quart_danish, 
   sensi5_ridge_quart_danish, 
   sensi5_elastic_net_quart_danish, 
   sensi5_model3_cox_sd_elastic_net_danish, 
   sensi5_model3_cox_quart_elastic_net_danish)


# Sensitivity analysis 6 - elastic net on all the POPs (ungrouped) ----
## Data prep ----
bdd_cases_danish_bis <- bdd_danish |>
  filter (als == 1) |>
  filter(follow_up_death>0) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_included)) |>
  mutate(across(all_of(POPs_included), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_included),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_28, PCB_28_quart) |>
  replace_with_median(PCB_52, PCB_52_quart) |>
  replace_with_median(PCB_74, PCB_74_quart) |>
  replace_with_median(PCB_99, PCB_99_quart) |>
  replace_with_median(PCB_101, PCB_101_quart) |>
  replace_with_median(PCB_118, PCB_118_quart) |>
  replace_with_median(PCB_138, PCB_138_quart) |>
  replace_with_median(PCB_153, PCB_153_quart) |>
  replace_with_median(PCB_156, PCB_156_quart) |>
  replace_with_median(PCB_170, PCB_170_quart) |>
  replace_with_median(PCB_180, PCB_180_quart) |>
  replace_with_median(PCB_183, PCB_183_quart) |>
  replace_with_median(PCB_187, PCB_187_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_pp_DDT, OCP_pp_DDT_quart) |>
  replace_with_median(OCP_pp_DDE, OCP_pp_DDE_quart) |>
  replace_with_median(OCP_oxychlordane, OCP_oxychlordane_quart) |>
  replace_with_median(OCP_transnonachlor, OCP_transnonachlor_quart) |>
  replace_with_median(PBDE_47, PBDE_47_quart) |>
  replace_with_median(PBDE_99, PBDE_99_quart) |>
  replace_with_median(PBDE_153, PBDE_153_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))

X_matrix_sd <- model.matrix(~ ., data = bdd_cases_danish_bis[, c(covariates_danish, POPs_included_sd)])[, -1] 
X_matrix_quart <- model.matrix(~ ., data = bdd_cases_danish_bis[, c(covariates_danish, POPs_included_quart)])[, -1] 

penalty_factor_sd <-                                                            # création d'un penalty factor pour forcer le modele lasso à inclure les covariables dans le modele 
  ifelse(colnames(X_matrix_sd) %in% 
           colnames(model.matrix(~ ., data = bdd_cases_danish_bis[, covariates_danish])[, -1]), 
         0, 1)

penalty_factor_quart <-                                                         # création d'un penalty factor pour forcer le modele lasso à inclure les covariables dans le modele 
  ifelse(colnames(X_matrix_quart) %in% 
           colnames(model.matrix(~ ., data = bdd_cases_danish_bis[, covariates_danish])[, -1]), 
         0, 1)


## LASSO POPs selection ----
set.seed(1996)
sensi6_lasso_sd_danish <- cv.glmnet(                                             # lasso + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_sd)                                           # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi6_lasso_sd_danish)
coef(sensi6_lasso_sd_danish, s = "lambda.min")                                   # lasso keeps only chlordane 
coef(sensi6_lasso_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi6_lasso_quart_danish <- cv.glmnet(                                          # lasso + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_quart)                                        # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi6_lasso_quart_danish)
coef(sensi6_lasso_quart_danish, s = "lambda.min")                                # lasso doesn't keep any pollutant when they are quartiles
coef(sensi6_lasso_quart_danish, s = "lambda.1se")  


## Elastic net selection -----
### choose the best alpha ----
alpha_grid <- seq(0.3, 1, by = 0.1)  # de 0 (Ridge) à 1 (Lasso)                 # starting at least at 0.3 because we want to select variables, not just run a ridge 

cv_results <- list()
cvm_min <- numeric(length(alpha_grid))

for (i in seq_along(alpha_grid)) {   # Boucle pour trouver le meilleur alpha
  set.seed(1996)
  cv_model <- cv.glmnet(
    x = X_matrix_sd,
    y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),
    family = "cox",
    alpha = alpha_grid[i],
    type.measure = "deviance",
    nfolds = 10,
    standardize = FALSE,
    penalty.factor = penalty_factor_sd
  )
  cv_results[[i]] <- cv_model
  cvm_min[i] <- min(cv_model$cvm)  # stocke la deviance minimale
}

best_alpha <- alpha_grid[which.min(cvm_min)]
best_alpha

rm(alpha_grid, cv_results, cvm_min, i, cv_model, best_alpha)

### analysis ----
set.seed(1996)
sensi6_elastic_net_sd_danish <- cv.glmnet(                                       # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.4,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                          

plot(sensi6_elastic_net_sd_danish)
coef(sensi6_elastic_net_sd_danish, s = "lambda.min")                             # Elastic net keeps HCB, β-HCH and Σchlordane
coef(sensi6_elastic_net_sd_danish, s = "lambda.1se")

set.seed(1996)
sensi6_elastic_net_quart_danish <- cv.glmnet(                                    # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.3,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                          

plot(sensi6_elastic_net_quart_danish)
coef(sensi6_elastic_net_quart_danish, s = "lambda.min")                          # elastic net doesn't keep any pollutant when they are quartiles
coef(sensi6_elastic_net_quart_danish, s = "lambda.1se")
rm(bdd_cases_danish_bis)

### new copollutant model with selected POPs ----
#### sd ----
POPs_sd_selected <-                                                             # selected POPs by LASSO
  c("PCB_28_sd", "PCB_52_sd", "OCP_oxychlordane_sd", "OCP_transnonachlor_sd", "PBDE_47_sd")                 
POPs_sd_selected_labels <- 
  set_names(c("PCB_28", "PCB_52", "OCP_oxychlordane", "OCP_transnonachlor", "PBDE_47"), 
            c("PCB-28", "PCB-52", "Oxychlordane", "Transnonachlor", "PBDE-47"))

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(POPs_sd_selected, covariates_danish), collapse = " + ")))

model_summary <- 
  coxph(formula_danish, data = bdd_cases_danish) |> summary() 
coefs <- model_summary$coefficients
sensi6_model3_cox_sd_elastic_net_danish <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "copollutant_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

#### quartiles ----
POPs_quart_selected <-                                                             # selected POPs by LASSO
  c("PCB_28_quart", "PCB_52_quart", "OCP_oxychlordane_quart", "OCP_transnonachlor_quart", "PBDE_47_quart")                 
POPs_quart_selected_labels <- 
  set_names(c("PCB_28", "PCB_52", "OCP_oxychlordane", "OCP_transnonachlor", "PBDE_47"), 
            c("PCB-28", "PCB-52", "Oxychlordane", "Transnonachlor", "PBDE-47"))

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(POPs_quart_selected, covariates_danish), collapse = " + ")))

model_summary <- 
  coxph(formula_danish, data = bdd_cases_danish) |> summary() 
coefs <- model_summary$coefficients
sensi6_model3_cox_quart_elastic_net_danish <- tibble(                                  # creation of a table of results
  study = "Danish", 
  model = "copollutant_quart", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_quart"))

##### heterogeneity tests ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_PCB_28_full <- 
  gam(outcome ~ 
        PCB_28_quart + 
        s(PCB_52) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_PCB_28_raw <- 
  gam(outcome ~ 
        #PCB_28_quart + 
        s(PCB_52) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_PCB_28_raw, model3_quart_PCB_28_full, test = "Chisq")
p.value_heterogeneity_PCB_28 <- tibble(explanatory = "PCB_28_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_PCB_52_full <- 
  gam(outcome ~ 
        PCB_52_quart + 
        s(PCB_28) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_PCB_52_raw <- 
  gam(outcome ~ 
        #PCB_52_quart + 
        s(PCB_28) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_PCB_52_raw, model3_quart_PCB_52_full, test = "Chisq")
p.value_heterogeneity_PCB_52 <- tibble(explanatory = "PCB_52_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_oxychlordane_full <- 
  gam(outcome ~ 
        OCP_oxychlordane_quart + 
        s(PCB_28) + s(PCB_52) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_oxychlordane_raw <- 
  gam(outcome ~ 
        #OCP_oxychlordane_quart + 
        s(PCB_28) + s(PCB_52) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_oxychlordane_raw, model3_quart_oxychlordane_full, test = "Chisq")
p.value_heterogeneity_oxychlordane <- tibble(explanatory = "OCP_oxychlordane_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_transnonachlor_full <- 
  gam(outcome ~ 
        OCP_transnonachlor_quart +
        s(PCB_28) + s(PCB_52) + s(OCP_oxychlordane) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_transnonachlor_raw <- 
  gam(outcome ~ 
        #OCP_transnonachlor_quart +
        s(PCB_28) + s(PCB_52) + s(OCP_oxychlordane) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_transnonachlor_raw, model3_quart_transnonachlor_full, test = "Chisq")
p.value_heterogeneity_transnonachlor <- tibble(explanatory = "OCP_transnonachlor_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

model3_quart_PBDE_47_full <- 
  gam(outcome ~ 
        PBDE_47_quart +
        s(PCB_28) + s(PCB_52) + s(OCP_oxychlordane) + s(OCP_transnonachlor) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

model3_quart_PBDE_47_raw <- 
  gam(outcome ~ 
        #PBDE_47_quart +
        s(PCB_28) + s(PCB_52) + s(OCP_oxychlordane) + s(OCP_transnonachlor) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)

anova <- anova(model3_quart_PBDE_47_raw, model3_quart_PBDE_47_full, test = "Chisq")
p.value_heterogeneity_PBDE_47 <- tibble(explanatory = "PBDE_47_quart", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

heterogeneity_copollutant_quart <- 
  bind_rows(p.value_heterogeneity_PCB_28, 
            p.value_heterogeneity_PCB_52,
            p.value_heterogeneity_oxychlordane,
            p.value_heterogeneity_transnonachlor,
            p.value_heterogeneity_PBDE_47) |>
  mutate(model = "copollutant_quart", 
         explanatory = gsub("_quart", "", explanatory), 
         study = "Danish")

rm(anova, outcome, 
   
   model3_quart_PCB_28_full, model3_quart_PCB_28_raw, 
   model3_quart_PCB_52_full, model3_quart_PCB_52_raw, 
   model3_quart_oxychlordane_full, model3_quart_oxychlordane_raw, 
   model3_quart_transnonachlor_full, model3_quart_transnonachlor_raw, 
   model3_quart_PBDE_47_full, model3_quart_PBDE_47_raw, 
   
   p.value_heterogeneity_PCB_28, 
   p.value_heterogeneity_PCB_52,
   p.value_heterogeneity_oxychlordane,
   p.value_heterogeneity_transnonachlor,
   p.value_heterogeneity_PBDE_47)

##### trend tests ----
outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))

model3_quart_PCB_28_trend <- 
  gam(outcome ~ 
        PCB_28_quart_med + 
        s(PCB_52) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)|> 
  summary()
p.value_trend_PCB_28 <- model3_quart_PCB_28_trend$p.table["PCB_28_quart_med", "Pr(>|z|)"]

model3_quart_PCB_52_trend <- 
  gam(outcome ~ 
        PCB_52_quart_med + 
        s(PCB_28) + s(OCP_oxychlordane) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish)|> 
  summary()
p.value_trend_PCB_52 <- model3_quart_PCB_52_trend$p.table["PCB_52_quart_med", "Pr(>|z|)"]

model3_quart_oxychlordane_trend <- 
  gam(outcome ~ 
        OCP_oxychlordane_quart_med + 
        s(PCB_28) + s(PCB_52) + s(OCP_transnonachlor) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_oxychlordane <- model3_quart_oxychlordane_trend$p.table["OCP_oxychlordane_quart_med", "Pr(>|z|)"]

model3_quart_transnonachlor_trend <- 
  gam(outcome ~ 
        OCP_transnonachlor_quart_med + 
        s(PCB_28) + s(PCB_52) + s(OCP_transnonachlor) + s(OCP_oxychlordane) + 
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_transnonachlor <- model3_quart_transnonachlor_trend$p.table["OCP_transnonachlor_quart_med", "Pr(>|z|)"]

model3_quart_PBDE_47_trend <- 
  gam(outcome ~ 
        PBDE_47_quart_med + 
        s(PCB_28) + s(PCB_52) + s(OCP_oxychlordane) + s(PBDE_47) +
        sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
      family = cox.ph(), 
      method = 'ML',
      data = bdd_cases_danish) |> 
  summary()
p.value_trend_PBDE_47 <- model3_quart_PBDE_47_trend$p.table["PBDE_47_quart_med", "Pr(>|z|)"]

trend_copollutant <- 
  data.frame(explanatory = c("PCB_28_quart_med", 
                             "PCB_52_quart_med", 
                             "OCP_oxychlordane_quart_med", 
                             "OCP_transnonachlor_quart_med", 
                             "PBDE_47_quart_med"),
             model = "copollutant_quart",
             p.value_trend = c(p.value_trend_PCB_28,
                               p.value_trend_PCB_52, 
                               p.value_trend_oxychlordane, 
                               p.value_trend_transnonachlor,
                               p.value_trend_PBDE_47)) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory), 
         study = "Danish")

### ERS model avec les POP selectionnés par elastic net 0.4 ----
coef(sensi6_elastic_net_sd_danish, s = "lambda.min")                            # Elastic net keeps HCB and Σchlordane
cox_model <-                                                                    # regular co-pollutant cox model to get unpenalized coefficients 
  coxph(Surv(follow_up_death, status_death) ~ 
          PCB_28_sd + PCB_52_sd + OCP_oxychlordane_sd + OCP_transnonachlor_sd + PBDE_47_sd + 
          sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, 
        data = bdd_cases_danish)
summary(cox_model)
pollutants <- c("PCB_28_sd", "PCB_52_sd", "OCP_oxychlordane_sd", "OCP_transnonachlor_sd", "PBDE_47_sd")                                        # selected POPs
scaled_betas <- c(PCB_28_sd = 0.17657 / (0.17657 + 0.23046 + 0.17271 + 0.07953 + 0.14389),                      # weigths of the POPs depending on the scaled coefficients of the regular cox copollutant model
                  PCB_52_sd = 0.23046 / (0.17657 + 0.23046 + 0.17271 + 0.07953 + 0.14389), 
                  OCP_oxychlordane_sd = 0.17271 / (0.17657 + 0.23046 + 0.17271 + 0.07953 + 0.14389), 
                  OCP_transnonachlor_sd = 0.07953 / (0.17657 + 0.23046 + 0.17271 + 0.07953 + 0.14389), 
                  PBDE_47_sd = 0.14389 / (0.17657 + 0.23046 + 0.17271 + 0.07953 + 0.14389))
Z <- scale(bdd_cases_danish[, pollutants])   
bdd_cases_danish$ERS_score_from_elastic_net <- as.numeric(Z %*% scaled_betas)
bdd_cases_danish <-
  bdd_cases_danish |>
  mutate(ERS_score_from_elastic_net_sd = scale(ERS_score_from_elastic_net), 
         ERS_score_from_elastic_net_quart = factor(ntile(ERS_score_from_elastic_net, 4),
                                                   labels = c("Q1", "Q2", "Q3", "Q4")))


sensi6_cox_model_ERS_from_elastic_net_sd <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_sd + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <-  summary(sensi6_cox_model_ERS_from_elastic_net_sd)
coefs <- model_summary$coefficients
sensi6_cox_model_ERS_from_elastic_net_sd <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "ERS_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

sensi6_cox_model_ERS_from_elastic_net_quart <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_quart + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi6_cox_model_ERS_from_elastic_net_quart)
coefs <- model_summary$coefficients
sensi6_cox_model_ERS_from_elastic_net_quart <- tibble(                                  # creation of a table of results
  study = "Danish", 
  model = "ERS_quart", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_quart"))

rm(cox_model, pollutants, scaled_betas, Z, model_summary, coefs)


### assemblage ----
results_sensi6 <-
  bind_rows(sensi6_model3_cox_sd_elastic_net_danish, 
            sensi6_model3_cox_quart_elastic_net_danish,
            sensi6_cox_model_ERS_from_elastic_net_sd, 
            sensi6_cox_model_ERS_from_elastic_net_quart) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    explanatory = gsub("_quartQ2", "", explanatory),
    explanatory = gsub("_quartQ3", "", explanatory),
    explanatory = gsub("_quartQ4", "", explanatory),
    explanatory = gsub("_sd", "", explanatory),
    term = case_when(
      str_detect(term, "_sd") ~ "Continuous", 
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"), 
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`,
    `p-value_shape` = ifelse(`p-value_raw` < 0.05, "p-value<0.05", "p-value≥0.05"),
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(
    study, model, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`,
    lower_CI,  upper_CI)


results_sensi6 <- left_join(results_sensi6, heterogeneity_copollutant_quart, 
                            by = c("study", "model", "explanatory"))

results_sensi6 <- 
  left_join(results_sensi6, trend_copollutant, 
            by = c("study", "model", "explanatory")) 

results_sensi6 <- results_sensi6 |>
  mutate(
    p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity), 
    p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
    p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(outcome, 
   model3_quart_PCB_28_trend, 
   model3_quart_PCB_52_trend, 
   model3_quart_oxychlordane_trend, 
   model3_quart_transnonachlor_trend, 
   model3_quart_PBDE_47_trend, 
   
   p.value_trend_PCB_28,
   p.value_trend_PCB_52, 
   p.value_trend_oxychlordane, 
   p.value_trend_transnonachlor,
   p.value_trend_PBDE_47,
   
   heterogeneity_copollutant_quart, 
   trend_copollutant, 
   sensi6_model3_cox_sd_elastic_net_danish, 
   sensi6_model3_cox_quart_elastic_net_danish,
   sensi6_cox_model_ERS_from_elastic_net_sd, 
   sensi6_cox_model_ERS_from_elastic_net_quart)


### table sd copollutant ----
POPs_sd_ALS_table_sensi6_danish <- 
  results_sensi6 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_sd") |>
  select(explanatory, HR, "95% CI", "p-value") |>
  mutate(
    # explanatory = factor(explanatory, levels = POPs_sd_selected_labels), 
    explanatory = fct_recode(explanatory, !!!POPs_sd_selected_labels)) |> 
  flextable() |>
  add_footer_lines(
    "1The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
     2Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Copollutant model", "95% CI" = "Copollutant model", "p-value" = "Copollutant model") |>
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

### table quart copollutant ----
quartile1_rows <- 
  results_sensi6 |>
  filter(model == "copollutant_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '') |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") 

POPs_quart_ALS_table_sensi6_danish <- 
  results_sensi6 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi6_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi6_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  group_by(explanatory) |>
  mutate(
    `Heterogeneity test` = ifelse(term == "quartile 1", first(`Heterogeneity test`[term == "quartile 2"]), ""),
    `Trend test` = ifelse(term == "quartile 1", first(`Trend test`[term == "quartile 2"]), "")) |>
  ungroup() |>
  mutate(explanatory = factor(explanatory, levels = POPs_sd_selected_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_sd_selected_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    2Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    5Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", term = "Quartiles",
    "HR" = "Copollutant model", "95% CI" = "Copollutant model", "p-value" = "Copollutant model",  
    "Heterogeneity test" = "Copollutant model",  "Trend test" = "Copollutant model") |>
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


### figure sd copollutant ----
POPs_sd_ALS_figure_sensi6_danish <-
  results_sensi6 |>
  filter(model == "copollutant_sd") |>
  mutate(
    explanatory = factor(explanatory, levels = POPs_sd_selected_labels),
    explanatory = fct_recode(explanatory, !!!POPs_sd_selected_labels), 
    explanatory = fct_rev(explanatory)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() 

### figure quart copollutant ----
POPs_quart_ALS_figure_sensi6_danish <-
  results_sensi6 |>
  filter(model == "copollutant_quart") |>
  mutate(
    explanatory = factor(explanatory, levels = POPs_quart_selected_labels),
    explanatory = fct_recode(explanatory, !!!POPs_quart_selected_labels), 
    explanatory = str_replace(explanatory, "Transnonachlor", "Trans-\nnonachlor"), 
    explanatory = fct_relevel(explanatory, "PCB-28", "PCB-52", "Oxychlordane", "Trans-\nnonachlor", "PBDE-47"),
    term = str_replace(term, "quartile", "Quartile"), 
    term = fct_rev(term)) |>
  arrange(explanatory) |>
  ggplot(aes(
    x = term,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  facet_grid(rows = dplyr::vars(explanatory), switch = "y") +  
  coord_flip() 

rm(POPs_sd_selected, POPs_sd_selected_labels,
   POPs_quart_selected, POPs_quart_selected_labels,
   X_matrix_sd, X_matrix_quart, 
   penalty_factor_sd, penalty_factor_quart, 
   formula_danish, model_summary, coefs, 
   sensi6_lasso_quart_danish, 
   #sensi6_model_ridge_quart_danish, 
   sensi6_elastic_net_quart_danish)


### table sd ERS ----
POPs_sd_ALS_table_sensi6_ERS_danish <- results_sensi6 |>
  filter(study == "Danish") |>
  filter(model == "ERS_sd") |>
  select(explanatory, HR, "95% CI", "p-value") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net")) |> 
  flextable() |>
  add_footer_lines(
    "1Environmental risk score is the weighted sum of relevant pollutants (PCB-28, PCB-52, oxychlordane, transnonachlor, PBDE-47) selected with elastic net regularization.
     2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
     3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase of the environmental risk score.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Mixture model", "95% CI" = "Mixture model", "p-value" = "Mixture model") |>
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

### table quart ERS ----
quartile1_rows <- 
  results_sensi6 |>
  filter(model == "ERS_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") 

POPs_quart_ALS_table_sensi6_ERS_danish <- 
  results_sensi6 |>
  filter(study == "Danish") |>
  filter(model == "ERS_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi6_ERS_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi6_ERS_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99"), 
         explanatory = fct_recode(explanatory, 
                                  "Environmental risk score" = "ERS_score_from_elastic_net")) |>
  arrange(explanatory, term) |>
  flextable() |>
  add_footer_lines(
    "1Environmental risk score is the weighted sum of relevant pollutants (PCB-28, PCB-52, oxychlordane, transnonachlor, PBDE-47) selected with elastic net regularization.
    2The model was adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    3Estimated risk of ALS death when the ERS is at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", term = "Quartiles",
    "HR" = "Mixture model", "95% CI" = "Mixture model", "p-value" = "Mixture model") |>
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

### figure sd ERS ----
POPs_sd_ALS_figure_sensi6_ERS_danish <-
  results_sensi6 |>
  filter(model == "ERS_sd") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net")) |>
  ggplot(aes(
    x = explanatory,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  coord_flip() 


### figure quart ERS ----
POPs_quart_ALS_figure_sensi6_ERS_danish <-
  results_sensi6 |>
  filter(model == "ERS_quart") |>
  mutate(
    explanatory = fct_recode(explanatory, 
                             "Environmental risk score" = "ERS_score_from_elastic_net"),
    term = str_replace(term, "quartile", "Quartile"), 
    term = fct_rev(term)) |>
  ggplot(aes(
    x = term,
    y = HR,
    ymin = lower_CI,
    ymax = upper_CI,
    color = `p-value_shape`)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +                        # , scales = "free_x"
  scale_color_manual(values = c(
    "p-value<0.05" = "red",
    "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.y = element_text(hjust = 0.5)) +
  facet_grid(rows = dplyr::vars(explanatory), switch = "y") +  
  coord_flip() 


# Sensitivity analysis 7 - follow-up duration ----
# investigation du probleme de follow-up duration très courte par rapport à la réalité (3-5 ans)
sensi7_figure <- bdd_cases_danish |>
  mutate(
    als_year = year(als_date), 
    status_death = as.factor(as.character(status_death)),
    status_death_rec = fct_recode(status_death, 
                                  "Alive" = "0",
                                  "Deceased" = "1")) |>
  ggplot() +
  aes(
    y = als_year,
    x = follow_up_death,
    colour = status_death_rec) +
  geom_point() +
  scale_color_hue(direction = 1) +
  labs(
    y = "Year of ALS diagnosis",
    x = " Follow-up between ALS diagnosis and death or end of study",
    color = "Death statut at the end of the study") +
  theme_minimal()

sensi7_table <- tbl_merge(
  tbls = list(
    bdd_cases_danish |> 
      select(follow_up_death, status_death) |>
      tbl_summary(by = status_death) |>
      add_overall(), 
    bdd_cases_finnish |> 
      select(follow_up_death, status_death) |>
      tbl_summary(by = status_death) |>
      add_overall()), 
  tab_spanner = c('**Danish EPIC**', "**Finnish Health Surveys**"))

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
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.
    2All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking, BMI and marital status. 
    3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    4CI: Confidence interval.") |>
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
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '')

POPs_quart_ALS_table_danish <- main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_danish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend")) |>
  select(explanatory, term, contains("base"), contains("adjusted"), contains("copollutant")) |>
  group_by(explanatory) |>
  mutate(p.value_heterogeneity_base = ifelse(term == 'quartile 1', p.value_heterogeneity_base[term == 'quartile 2'], ''), 
         p.value_trend_base = ifelse(term == 'quartile 1', p.value_trend_base[term == 'quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(term == 'quartile 1', p.value_heterogeneity_adjusted[term == 'quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(term == 'quartile 1', p.value_trend_adjusted[term == 'quartile 2'], ''),
         p.value_heterogeneity_copollutant = ifelse(term == 'quartile 1', p.value_heterogeneity_copollutant[term == 'quartile 2'], ''), 
         p.value_trend_copollutant = ifelse(term == 'quartile 1', p.value_trend_copollutant[term == 'quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Heterogeneity test" = "p.value_heterogeneity_base", "Trend test" = "p.value_trend_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted",  "Heterogeneity test " = "p.value_heterogeneity_adjusted", "Trend test " = "p.value_trend_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant",  " Heterogeneity test " = "p.value_heterogeneity_copollutant", " Trend test " = "p.value_trend_copollutant") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, Σchlordane corresponds to trans-nonanchlor and oxychlordane and finally ΣPBDE corresponds to PBDEs 47, 99, 153.
    2All models are adjusted for sex and age at diagnosis. Adjusted models further account for smoking, BMI and marital status.
    3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.
    5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex and age at diagnosis.
    6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex and age at diagnosis.
    7Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", 
    term = "Quartiles",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",  "Trend test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model",  "Heterogeneity test " = "Adjusted model",  "Trend test " = "Adjusted model", 
    " HR " = "Copollutant model", " 95% CI " = "Copollutant model", " p-value " = "Copollutant model",  " Heterogeneity test " = "Copollutant model",  " Trend test " = "Copollutant model") |>
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
                            "Adjusted model" = "adjusted", 
                            "Co-pollutant model" = 'copollutant'),
         model = fct_relevel(model, 'Base model', 'Adjusted model', "Co-pollutant model"), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y") +                         # , scales = "free_x"
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
                            "Adjusted model" = "adjusted", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()


### table POPs - ALS survival (qgcomp analysis all POPs grouped) ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")                                  # remove the 4 most abundant PCB because they are already NDL-PCB
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_bis) 

p <- summary(qgcomp_boot_danish)
POPs_ALS_qgcomp_table_danish <-                                                 # overall results
  tibble(
    study = "Danish", 
    model = "copollutant", 
    HR = exp(qgcomp_boot_danish$psi),
    lower_CI = exp(qgcomp_boot_danish$ci[1]), 
    upper_CI = exp(qgcomp_boot_danish$ci[2]), 
    p_value = p$coefficients[1, "Pr(>|z|)"]) |>
  mutate(
    HR = sprintf("%.1f", HR),
    lower_CI = sprintf("%.2f", lower_CI),
    upper_CI = sprintf("%.2f", upper_CI),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''), 
    `p-value` = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, HR, `95% CI`, `p-value`)

### figure POPs - ALS survival (qgcomp analysis all POPs grouped) ----
POPs_ALS_qgcomp_figure_danish <- 
  tibble(
    pollutant = c(names(qgcomp_noboot_danish$pos.weights), names(qgcomp_noboot_danish$neg.weights)),
    weight = c(qgcomp_noboot_danish$pos.weights, - qgcomp_noboot_danish$neg.weights)) |>
  mutate(
    pollutant_label = pollutant_labels_bis[pollutant] %||% pollutant,
    pollutant_label = factor(pollutant_label, levels = rev(pollutant_labels_bis))) |>
  ggplot(
    aes(x = weight, y = pollutant_label, fill = weight > 0)) +
  geom_col(show.legend = FALSE) +
  labs(y = "Exposures", x = "       Negative weights                             Positive weights") +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_lucid()


### table POPs - ALS survival (qgcomp analysis positive POPs grouped) ----
p <- summary(qgcomp_positive_boot_danish)
POPs_ALS_qgcomp_positive_table_danish <-                                                 # overall results
  tibble(
    study = "Danish", 
    model = "copollutant", 
    HR = exp(qgcomp_positive_boot_danish$psi),
    lower_CI = exp(qgcomp_positive_boot_danish$ci[1]), 
    upper_CI = exp(qgcomp_positive_boot_danish$ci[2]), 
    p_value = p$coefficients[1, "Pr(>|z|)"]) |>
  mutate(
    HR = sprintf("%.1f", HR),
    lower_CI = sprintf("%.2f", lower_CI),
    upper_CI = sprintf("%.2f", upper_CI),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''), 
    `p-value` = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, HR, `95% CI`, `p-value`)

### figure POPs - ALS survival (qgcomp analysis positive POPs grouped) ----
POPs_ALS_qgcomp_positive_figure_danish <- 
  tibble(
    pollutant = c(names(qgcomp_positive_noboot_danish$pos.weights), names(qgcomp_positive_noboot_danish$neg.weights)),
    weight = c(qgcomp_positive_noboot_danish$pos.weights, - qgcomp_positive_noboot_danish$neg.weights)) |>
  mutate(
    pollutant_label = pollutant_labels_bis[pollutant] %||% pollutant,
    pollutant_label = factor(pollutant_label, levels = rev(pollutant_labels_bis))) |>
  ggplot(
    aes(x = weight, y = pollutant_label, fill = weight > 0)) +
  geom_col(show.legend = FALSE) +
  labs(y = "Exposures", x = "Positive weights") +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_lucid()

rm(p, POPs_group_bis, pollutant_labels_bis)

## Finnish ----
### table covariates - als survival ----
ref_rows <- covar_finnish |>
  distinct(model, variable) |>
  arrange(model, variable) |>
  mutate(
    term = c("continuous",  "continuous", "continuous", "1", "Married/cohabit","Male","Ever",
             "continuous",  "continuous", "continuous", "1", "Married/cohabit","Male","Ever"),
    HR = "-",
    "95% CI" = "-",
    `p-value` = "" ,
    lower_CI = "", 
    upper_CI = "")

covar_ALS_table_finnish <- covar_finnish |>
  mutate(
    term = fct_recode(term, 
                      "continuous" = "diagnosis_age",
                      "continuous" = "bmi",
                      "continuous" = "cholesterol",
                      "2" = "level_urbanization2",
                      "3" = "level_urbanization3",
                      "4" = "level_urbanization4",
                      "Other" = "marital_status_2catOther",
                      "Female" = "sexFemale",
                      "Never" = "smoking_2catNever"),
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
  mutate(
    variable = fct_recode(variable, 
                          "Sex" = "sex", 
                          "Age at diagnosis (years)" = "diagnosis_age", 
                          "Boby mass index (kg/m²)" = "bmi", 
                          "Marital status" = "marital_status_2cat", 
                          "Smoking status" = "smoking_2cat", 
                          "Serum cholesterol (mmol/L)" = "cholesterol", 
                          "Level of urbanization" = "level_urbanization")) |>
  select(model, variable, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI) |>
  mutate(
    term = fct_relevel(term, 
    "1", "2", "3", "4", "continuous", "Ever", "Never", "Male", "Female", "Married/cohabit",  "Other"), 
    model = fct_relevel(model, "base", "adjusted"))|>
  arrange(model, variable, term) |>
  filter(!(term == "continuous" & HR == "-"))

covar_ALS_table_finnish <- covar_ALS_table_finnish |>
  select(-`p-value_raw`, -`p-value_shape`, -lower_CI, -upper_CI) |>
  pivot_wider(names_from = model, values_from = c(HR, `95% CI`, `p-value`)) |>
  select(variable, term, ends_with("base"), ends_with("adjusted")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted") |>
  flextable() |>
  add_header(
    "variable" = "Characteristics", 
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
  merge_h(part = "header") |>
  merge_v(j = "variable") |>
  theme_vanilla() |>
  bold(j = "variable", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "variable", align = "left", part = "all") |> 
  merge_at(j = "variable", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
  
rm(ref_rows)

### table POPs (sd) - als survival (meta-analaysis) ----
POPs_sd_ALS_table_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "meta-analysis") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity")) |>
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Hetero-geneity test" = "p.value_heterogeneity_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "Hetero-geneity test " = "p.value_heterogeneity_adjusted", 
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant", " Hetero-geneity test " = "p.value_heterogeneity_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |> 
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking, BMI and marital status. 
    3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "Hetero-geneity test" = "Base model", 
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Hetero-geneity test " = "Adjusted model", 
    " HR " = "Copollutant model", " 95% CI " = "Copollutant model", " p-value " = "Copollutant model", " Hetero-geneity test " = "Copollutant model") |>
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

### table POPs (quart) - als survival (meta-analaysis) ----
quartile1_rows <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "meta-analysis") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = "")

POPs_quart_ALS_table_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "meta-analysis") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_finnish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_finnish) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity")) |>
  select(explanatory, term, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Hetero-geneity test" = "p.value_heterogeneity_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "Hetero-geneity test " = "p.value_heterogeneity_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant", " Hetero-geneity test " = "p.value_heterogeneity_copollutant") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels_finnish), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2All models are adjusted for age at diagnosis and sex. Adjusted models further account for smoking, BMI and marital status.
    3Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "Hetero-geneity test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model",  "Hetero-geneity test " = "Adjusted model",
    " HR " = "Copollutant model", " 95% CI " = "Copollutant model", " p-value " = "Copollutant model", " Hetero-geneity test " = "Copollutant model") |>
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

### table POPs (sd) - als survival (raw analysis not adjusted on cohort) ----
POPs_sd_ALS_table_finnish_raw <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!study_design == "meta-analysis") |>
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", 
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |> 
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2All models are adjusted for age and sex. Adjusted models further account for smoking, BMI and marital status. 
    3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    4CI: Confidence interval.") |>
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


### table POPs (quart) - als survival (raw analysis not adjusted on cohort) ----
quartile1_rows <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '')

POPs_quart_ALS_table_finnish_raw <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
  filter(study_design == "raw unadjusted on cohort") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_finnish_raw <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_finnish_raw) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend")) |>
  select(explanatory, term, contains("base"), contains("adjusted"), contains("copollutant")) |>
  group_by(explanatory) |>
  mutate(p.value_heterogeneity_base = ifelse(term == 'quartile 1', p.value_heterogeneity_base[term == 'quartile 2'], ''), 
         p.value_trend_base = ifelse(term == 'quartile 1', p.value_trend_base[term == 'quartile 2'], ''),
         p.value_heterogeneity_adjusted = ifelse(term == 'quartile 1', p.value_heterogeneity_adjusted[term == 'quartile 2'], ''), 
         p.value_trend_adjusted = ifelse(term == 'quartile 1', p.value_trend_adjusted[term == 'quartile 2'], ''),
         p.value_heterogeneity_copollutant = ifelse(term == 'quartile 1', p.value_heterogeneity_copollutant[term == 'quartile 2'], ''), 
         p.value_trend_copollutant = ifelse(term == 'quartile 1', p.value_trend_copollutant[term == 'quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Heterogeneity test" = "p.value_heterogeneity_base", "Trend test" = "p.value_trend_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted",  "Heterogeneity test " = "p.value_heterogeneity_adjusted", "Trend test " = "p.value_trend_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant",  " Heterogeneity test " = "p.value_heterogeneity_copollutant", " Trend test " = "p.value_trend_copollutant") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels_finnish), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2All models are adjusted for sex and age at diagnosis. Adjusted models further account for smoking, BMI and marital status.
    3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
    4CI: Confidence interval.
    5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex and age at diagnosis.
    6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex and age at diagnosis.
    7Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
    8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", 
    term = "Quartiles",
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",  "Trend test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model",  "Heterogeneity test " = "Adjusted model",  "Trend test " = "Adjusted model", 
    " HR " = "Copollutant model", " 95% CI " = "Copollutant model", " p-value " = "Copollutant model",  " Heterogeneity test " = "Copollutant model",  " Trend test " = "Copollutant model") |>
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

### figure POPs (sd) - als survival (meta-analaysis) ----
POPs_sd_ALS_figure_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "meta-analysis") |>
  filter(term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted",
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
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

### figure POPs (quart) - als survival (meta-analaysis) ----
POPs_quart_ALS_figure_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "meta-analysis") |>
  filter(!term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
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

### figure POPs (sd) - als survival (raw analysis not adjusted on cohort) ----
POPs_sd_ALS_figure_finnish_raw <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(study_design == "raw unadjusted on cohort") |>
  filter(term == "Continuous") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted",
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |>
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


### figure POPs (quart) - als survival (raw analysis not adjusted on cohort) ----
POPs_quart_ALS_figure_finnish_raw <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
  filter(!term == "Continuous") |>
  filter(study_design == "raw unadjusted on cohort") |>
  mutate(model = fct_recode(model, 
                            "Base model" = "base",
                            "Adjusted model" = "adjusted", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_finnish),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()


### table POPs (sd) - ALS survival (qgcomp analysis) ----
POPs_group_finnish_bis <-                                                 # remove the 4 most abundant PCB because they are already NDL-PCB
  setdiff(POPs_group_finnish, 
          c("OCP_β_HCH",  "OCP_γ_HCH", "PCB_4"))                    
pollutant_labels_finnish_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "ΣHCH", "Σchlordane", "PeCB"), 
  POPs_group_finnish_bis)

p <- summary(qgcomp_boot_finnish)
POPs_ALS_qgcomp_table_finnish <-                                                 # overall results
  tibble(
    study = "Finnish", 
    model = "copollutant", 
    HR = exp(qgcomp_boot_finnish$psi),
    lower_CI = exp(qgcomp_boot_finnish$ci[1]), 
    upper_CI = exp(qgcomp_boot_finnish$ci[2]), 
    p_value = p$coefficients[1, "Pr(>|z|)"]) |>
  mutate(
    HR = sprintf("%.1f", HR),
    lower_CI = sprintf("%.2f", lower_CI),
    upper_CI = sprintf("%.2f", upper_CI),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''), 
    `p-value` = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, HR, `95% CI`, `p-value`)

### figure POPs (sd) - ALS survival (qgcomp analysis) ----
POPs_ALS_qgcomp_figure_finnish <- 
  tibble(
    pollutant = c(names(qgcomp_noboot_finnish$pos.weights), names(qgcomp_noboot_finnish$neg.weights)),
    weight = c(qgcomp_noboot_finnish$pos.weights, - qgcomp_noboot_finnish$neg.weights)) |>
  mutate(
    pollutant_label = pollutant_labels_finnish_bis[pollutant] %||% pollutant,
    pollutant_label = factor(pollutant_label, levels = rev(pollutant_labels_finnish_bis))) |>
  ggplot(
    aes(x = weight, y = pollutant_label, fill = weight > 0)) +
  geom_col(show.legend = FALSE) +
  labs(y = "Exposures", x = "Negative weights                             Positive weights") +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_lucid()

rm(qgcomp_boot_finnish, qgcomp_noboot_finnish, p, POPs_group_finnish_bis, pollutant_labels_finnish_bis)


## Metanalysis ----
POPs_group_labels_metanalysis <- c(
  "Most prevalent PCBs" = "PCB_4",
  "Dioxin-like PCBs" = "PCB_DL",
  "Non dioxin-like PCBs" = "PCB_NDL",
  "HCB" = "OCP_HCB",
  "ΣDDT" = "ΣDDT",
  "β-HCH" = "OCP_β_HCH",
  "Σchlordane" = "Σchlordane")

### table POPs (sd) - als survival ----
POPs_sd_ALS_table_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity")) |>
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Hetero-geneity test" = "p.value_heterogeneity_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "Hetero-geneity test " = "p.value_heterogeneity_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant", " Hetero-geneity test " = "p.value_heterogeneity_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_metanalysis)) |> 
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH, ΣPBDE corresponds to PBDEs 47, 99, 153 and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
    2Base models were all adjusted for age at diagnosis and sex. Adjusted models all further account for smoking, BMI and marital status. 
    3Estimated risk of death after ALS diagnosis associated with a one standard deviation increase in pre-disease serum concentration of POPs.
    4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Exposures", 
    "HR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", "Hetero-geneity test" = "Base model",
    "HR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Hetero-geneity test " = 'Adjusted model',
    " HR " = "Copollutant model", " 95% CI " = "Copollutant model", " p-value " = "Copollutant model", " Hetero-geneity test " = 'Copollutant model') |>
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
    `p-value` = "", 
    "p.value_heterogeneity" = "")

POPs_quart_ALS_table_metanalysis <- main_results_POPs_ALS_survival |>
  filter(study == "Metanalysis") |>
  filter(!term == "Continuous") |>
  select(model, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_metanalysis <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_metanalysis) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity")) |>
  select(explanatory, term, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", "Hetero-geneity test" = "p.value_heterogeneity_base",
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", "Hetero-geneity test " = "p.value_heterogeneity_adjusted",
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant", " Hetero-geneity test " = "p.value_heterogeneity_copollutant") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels_metanalysis), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_metanalysis)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
  "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE, ΣHCH corresponds to β-HCH and γ-HCH, ΣPBDE corresponds to PBDEs 47, 99, 153 and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2Base models were all adjusted for age at diagnosis and sex. Adjusted models all further account for smoking, BMI and marital status. 
  3Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  4CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", "Hetero-geneity test" = "Base Model", 
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model", "Hetero-geneity test " = "Adjusted Model", 
    " HR " = "Copollutant Model", " 95% CI " = "Copollutant Model", " p-value " = "Copollutant Model", " Hetero-geneity test " = "Copollutant Model") |>
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
                            "Adjusted model" = "adjusted", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_metanalysis),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_metanalysis)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y") +  
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
                            "Adjusted model" = "adjusted", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         explanatory = factor(explanatory, levels = POPs_group_labels_metanalysis),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels_metanalysis), 
         term = fct_rev(term)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = term, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()

rm(POPs_group_labels_metanalysis)

# Assemblage ----
results_POPs_ALS_survival <- 
  list(main_analysis = list(main_results_POPs_ALS_survival = main_results_POPs_ALS_survival), 
       danish = list(
         covar_danish = covar_danish, 
         POPs_sd_ALS_table_danish = POPs_sd_ALS_table_danish, 
         POPs_quart_ALS_table_danish = POPs_quart_ALS_table_danish, 
         POPs_sd_ALS_figure_danish = POPs_sd_ALS_figure_danish, 
         cheking_model3_cox_sd_danish = cheking_model3_cox_sd_danish, 
         POPs_quart_ALS_figure_danish = POPs_quart_ALS_figure_danish, 
         plot_base_cox_gam_danish = plot_base_cox_gam_danish, 
         plot_adjusted_cox_gam_danish = plot_adjusted_cox_gam_danish, 
         plot_copollutant_cox_gam_danish = plot_copollutant_cox_gam_danish, 
         qgcomp = list(qgcomp_boot_danish = qgcomp_boot_danish, 
                       qgcomp_noboot_danish = qgcomp_noboot_danish, 
                       qgcomp_positive_boot_danish = qgcomp_positive_boot_danish, 
                       qgcomp_positive_noboot_danish = qgcomp_positive_noboot_danish, 
                       POPs_ALS_qgcomp_table_danish = POPs_ALS_qgcomp_table_danish,
                       POPs_ALS_qgcomp_figure_danish = POPs_ALS_qgcomp_figure_danish, 
                       POPs_ALS_qgcomp_positive_table_danish = POPs_ALS_qgcomp_positive_table_danish,
                       POPs_ALS_qgcomp_positive_figure_danish = POPs_ALS_qgcomp_positive_figure_danish)), 
       finnish = list(
         covar_finnish = covar_finnish, 
         covar_ALS_table_finnish = covar_ALS_table_finnish,
         POPs_sd_ALS_table_finnish = POPs_sd_ALS_table_finnish, 
         POPs_quart_ALS_table_finnish = POPs_quart_ALS_table_finnish, 
         POPs_sd_ALS_table_finnish_raw = POPs_sd_ALS_table_finnish_raw, 
         POPs_quart_ALS_table_finnish_raw = POPs_quart_ALS_table_finnish_raw, 
         cheking_model3_cox_sd_finnish = cheking_model3_cox_sd_finnish, 
         POPs_sd_ALS_figure_finnish = POPs_sd_ALS_figure_finnish, 
         POPs_quart_ALS_figure_finnish = POPs_quart_ALS_figure_finnish, 
         POPs_sd_ALS_figure_finnish_raw = POPs_sd_ALS_figure_finnish_raw, 
         POPs_quart_ALS_figure_finnish_raw = POPs_quart_ALS_figure_finnish_raw, 
         POPs_ALS_qgcomp_table_finnish = POPs_ALS_qgcomp_table_finnish, 
         POPs_ALS_qgcomp_figure_finnish = POPs_ALS_qgcomp_figure_finnish), 
       metanalysis = list(
         POPs_sd_ALS_table_metanalysis = POPs_sd_ALS_table_metanalysis, 
         POPs_quart_ALS_table_metanalysis = POPs_quart_ALS_table_metanalysis, 
         POPs_sd_ALS_figure_metanalysis = POPs_sd_ALS_figure_metanalysis, 
         POPs_quart_ALS_figure_metanalysis = POPs_quart_ALS_figure_metanalysis), 
       sensi1 = list(
         plot_justif = plot_justif, 
         plot_justif_2cat = plot_justif_2cat, 
         cor_test_sensi1 = cor_test_sensi1, 
         results_sensi1 = results_sensi1, 
         sample_size_check = sample_size_check, 
         sample_size_check_2cat = sample_size_check_2cat, 
         POPs_heatmap_cases_group = POPs_heatmap_cases_group, 
         POPs_sd_ALS_figure_sensi1_not_adjusted = POPs_sd_ALS_figure_sensi1_not_adjusted, 
         POPs_sd_ALS_figure_sensi1_adjusted = POPs_sd_ALS_figure_sensi1_adjusted, 
         POPs_sd_ALS_figure_sensi1_adjusted_2cat = POPs_sd_ALS_figure_sensi1_adjusted_2cat), 
       sensi2 = list(
         results_sensi2 = results_sensi2, 
         POPs_sd_ALS_figure_sensi2_finnish_orange = POPs_sd_ALS_figure_sensi2_finnish_orange, 
         POPs_sd_ALS_figure_sensi2_finnish_purple = POPs_sd_ALS_figure_sensi2_finnish_purple, 
         POPs_sd_ALS_figure_sensi2_finnish_green = POPs_sd_ALS_figure_sensi2_finnish_green, 
         tbl_distrib_among_green_box = tbl_distrib_among_green_box, 
         tbl_distrib_among_purple_box = tbl_distrib_among_purple_box, 
         tbl_distrib_among_orange_box = tbl_distrib_among_orange_box), 
       sensi3 = list(
         results_sensi3 = results_sensi3, 
         POPs_sd_ALS_table_sensi3_baseline_age = POPs_sd_ALS_table_sensi3_baseline_age, 
         POPs_sd_ALS_table_sensi3_follow_up = POPs_sd_ALS_table_sensi3_follow_up, 
         POPs_sd_ALS_table_sensi3_baseline_age_follow_up = POPs_sd_ALS_table_sensi3_baseline_age_follow_up), 
       sensi4 = list(
         plot_justif_sensi4 = plot_justif_sensi4, 
         results_sensi4 = results_sensi4, 
         POPs_sd_ALS_figure_sensi4_danish = POPs_sd_ALS_figure_sensi4_danish), 
       sensi5 = list(
         sensi5_lasso_sd_danish = sensi5_lasso_sd_danish, 
         sensi5_ridge_sd_danish = sensi5_ridge_sd_danish, 
         sensi5_elastic_net_sd_danish = sensi5_elastic_net_sd_danish, 
         results_sensi5 = results_sensi5, 
         POPs_sd_ALS_figure_sensi5_danish = POPs_sd_ALS_figure_sensi5_danish, 
         POPs_quart_ALS_figure_sensi5_danish = POPs_quart_ALS_figure_sensi5_danish, 
         POPs_sd_ALS_table_sensi5_danish = POPs_sd_ALS_table_sensi5_danish, 
         POPs_quart_ALS_table_sensi5_danish = POPs_quart_ALS_table_sensi5_danish, 
         POPs_sd_ALS_figure_sensi5_ERS_danish = POPs_sd_ALS_figure_sensi5_ERS_danish, 
         POPs_quart_ALS_figure_sensi5_ERS_danish = POPs_quart_ALS_figure_sensi5_ERS_danish, 
         POPs_sd_ALS_table_sensi5_ERS_danish = POPs_sd_ALS_table_sensi5_ERS_danish, 
         POPs_quart_ALS_table_sensi5_ERS_danish = POPs_quart_ALS_table_sensi5_ERS_danish), 
       sensi6 = list(
         sensi6_lasso_sd_danish = sensi6_lasso_sd_danish, 
         # sensi6_ridge_sd_danish = sensi6_ridge_sd_danish, 
         sensi6_elastic_net_sd_danish = sensi6_elastic_net_sd_danish, 
         results_sensi6 = results_sensi6, 
         POPs_sd_ALS_figure_sensi6_danish = POPs_sd_ALS_figure_sensi6_danish, 
         POPs_quart_ALS_figure_sensi6_danish = POPs_quart_ALS_figure_sensi6_danish, 
         POPs_sd_ALS_table_sensi6_danish = POPs_sd_ALS_table_sensi6_danish, 
         POPs_quart_ALS_table_sensi6_danish = POPs_quart_ALS_table_sensi6_danish, 
         POPs_sd_ALS_figure_sensi6_ERS_danish = POPs_sd_ALS_figure_sensi6_ERS_danish, 
         POPs_quart_ALS_figure_sensi6_ERS_danish = POPs_quart_ALS_figure_sensi6_ERS_danish, 
         POPs_sd_ALS_table_sensi6_ERS_danish = POPs_sd_ALS_table_sensi6_ERS_danish, 
         POPs_quart_ALS_table_sensi6_ERS_danish = POPs_quart_ALS_table_sensi6_ERS_danish), 
       sensi7 = list(
         sensi7_figure = sensi7_figure, 
         sensi7_table = sensi7_table))

rm(bdd_cases_tot, 
   bdd_cases_danish, 
   bdd_cases_finnish, 
   main_results_POPs_ALS_survival, 
   covar_danish,
   POPs_sd_ALS_table_danish, 
   cheking_model3_cox_sd_danish, 
   POPs_quart_ALS_table_danish, 
   POPs_sd_ALS_figure_danish, 
   POPs_quart_ALS_figure_danish, 
   plot_base_cox_gam_danish, 
   plot_adjusted_cox_gam_danish, 
   plot_copollutant_cox_gam_danish,
   POPs_ALS_qgcomp_table_danish, 
   POPs_ALS_qgcomp_figure_danish,
   POPs_ALS_qgcomp_positive_table_danish, 
   POPs_ALS_qgcomp_positive_figure_danish,
   qgcomp_boot_danish, 
   qgcomp_noboot_danish, 
   qgcomp_positive_boot_danish, 
   qgcomp_positive_noboot_danish, 
   
   covar_finnish,
   covar_ALS_table_finnish, 
   POPs_sd_ALS_table_finnish, 
   POPs_quart_ALS_table_finnish, 
   POPs_sd_ALS_table_finnish_raw,
   POPs_quart_ALS_table_finnish_raw, 
   POPs_sd_ALS_figure_finnish, 
   POPs_quart_ALS_figure_finnish,
   POPs_sd_ALS_figure_finnish_raw,
   POPs_quart_ALS_figure_finnish_raw, 
   cheking_model3_cox_sd_finnish, 
   POPs_ALS_qgcomp_table_finnish,
   POPs_ALS_qgcomp_figure_finnish, 
   model2_cox_sd_finnish_raw, 
   
   POPs_sd_ALS_table_metanalysis, 
   POPs_quart_ALS_table_metanalysis, 
   POPs_sd_ALS_figure_metanalysis, 
   POPs_quart_ALS_figure_metanalysis, 
   
   covariates_danish, covariates_finnish, 
   surv_obj_danish, surv_obj_FMC, surv_obj_FMCF, surv_obj_MFH, surv_obj_finnish, 
   
   plot_justif, plot_justif_2cat, results_sensi1, sample_size_check, sample_size_check_2cat, 
   POPs_heatmap_cases_group, cor_test_sensi1, 
   model2_cox_sd_danish, 
   POPs_sd_ALS_figure_sensi1_not_adjusted, 
   POPs_sd_ALS_figure_sensi1_adjusted, 
   POPs_sd_ALS_figure_sensi1_adjusted_2cat,
   
   results_sensi2, 
   POPs_sd_ALS_figure_sensi2_finnish_orange, 
   POPs_sd_ALS_figure_sensi2_finnish_purple, 
   POPs_sd_ALS_figure_sensi2_finnish_green, 
   tbl_distrib_among_green_box, 
   tbl_distrib_among_purple_box, 
   tbl_distrib_among_orange_box, 
   
   results_sensi3, 
   POPs_sd_ALS_table_sensi3_baseline_age,
   POPs_sd_ALS_table_sensi3_follow_up, 
   POPs_sd_ALS_table_sensi3_baseline_age_follow_up,
   
   results_sensi4, 
   plot_justif_sensi4, 
   POPs_sd_ALS_figure_sensi4_danish, 
   
   results_sensi5, 
   sensi5_lasso_sd_danish, 
   sensi5_ridge_sd_danish, 
   sensi5_elastic_net_sd_danish, 
   POPs_sd_ALS_figure_sensi5_danish, 
   POPs_quart_ALS_figure_sensi5_danish, 
   POPs_sd_ALS_table_sensi5_danish, 
   POPs_quart_ALS_table_sensi5_danish, 
   POPs_sd_ALS_figure_sensi5_ERS_danish, 
   POPs_quart_ALS_figure_sensi5_ERS_danish, 
   POPs_sd_ALS_table_sensi5_ERS_danish, 
   POPs_quart_ALS_table_sensi5_ERS_danish,
   
   results_sensi6, 
   sensi6_lasso_sd_danish, 
   #sensi6_ridge_sd_danish, 
   sensi6_elastic_net_sd_danish, 
   POPs_sd_ALS_figure_sensi6_danish, 
   POPs_quart_ALS_figure_sensi6_danish, 
   POPs_sd_ALS_table_sensi6_danish, 
   POPs_quart_ALS_table_sensi6_danish, 
   
   POPs_sd_ALS_figure_sensi6_danish, 
   POPs_quart_ALS_figure_sensi6_danish, 
   POPs_sd_ALS_table_sensi6_danish, 
   POPs_quart_ALS_table_sensi6_danish, 
   POPs_sd_ALS_figure_sensi6_ERS_danish, 
   POPs_quart_ALS_figure_sensi6_ERS_danish, 
   POPs_sd_ALS_table_sensi6_ERS_danish, 
   POPs_quart_ALS_table_sensi6_ERS_danish,
   
   sensi7_figure, 
   sensi7_table)


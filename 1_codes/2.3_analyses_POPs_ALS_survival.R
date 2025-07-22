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
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(OCP_γ_HCH, OCP_γ_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                        event = bdd_cases_danish$status_death)
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
rm(POPs_group_sd_bis, pollutant_labels_bis, formula_danish, model_summary, coefs)

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


## Cox-gam model (sd) ----
### Base  ----
plot_base_cox_gam_danish <- map(POPs_group_sd, function(var) {
  
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
  x_label <- POPs_group_sd_labels[var] 
  
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
  x_label <- POPs_group_sd_labels[var] 
  
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
               sex + diagnosis_age + 
               smoking_2cat_i + bmi + marital_status_2cat_i, 
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

## Q-gcomp analysis ---- 
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
# qgcomp_boot_danish$pos.weights                                                  # NULL because the model is not significant 
# qgcomp_boot_danish$neg.weights                                                  # NULL because the model is not significant 
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
    expnms = POPs_group_bis)                                                 # exposures of interest
# print(qgcomp_noboot_danish)
# qgcomp_noboot_danish$pos.weights
# qgcomp_noboot_danish$neg.weights
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

## Cox model (sd) ----
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
      HR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
    p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |>
  mutate(explanatory = term)

rm(POPs_group_sd_finnish_bis, formula_FMC, formula_FMCF)

## Cox model (quart) ----
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
    mutate(study = "Finnish", model = "base") |> 
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
    mutate(model = "adjusted") |> 
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

## Investigation d'autres modeles pour les cohortes finnish ----



## Q-gcomp analysis ---- 
POPs_group_finnish_bis <-                                                 # remove the 4 most abundant PCB because they are already NDL-PCB
  setdiff(POPs_group_finnish, 
          c("OCP_β_HCH",  "OCP_γ_HCH", "PCB_4"))                    

bdd_cases_finnish <- bdd_finnish |>
  filter (als == 1) |>
  filter(study %in% c("FMC", "FMCF")) |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, 
         thawed, level_urbanization, marital_status_2cat, smoking_2cat, bmi, cholesterol, study, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", 
         OCP_γ_HCH = "OCP_γ_HCH_raw", "Σchlordane", OCP_PeCB = "OCP_PeCB_raw") |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"))

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

rm(formula_finnish, bdd_cases_finnish, POPs_group_finnish_bis)

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
    model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  
  bind_rows(model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish,
            model1_cox_quart_finnish, model2_cox_quart_finnish, model3_cox_quart_finnish,
            model1_cox_sd_metanalysis, model2_cox_sd_metanalysis, model3_cox_sd_metanalysis, 
            model1_cox_quart_metanalysis, model2_cox_quart_metanalysis, model3_cox_quart_metanalysis) |>
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
  select(study, model, explanatory, term, HR, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI, "p.value_heterogeneity")

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, heterogeneity_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_heterogeneity = ifelse(is.na(p.value_heterogeneity.x), p.value_heterogeneity.y, p.value_heterogeneity.x), 
         p.value_heterogeneity = ifelse(term == "Continuous" & study == "Danish", NA, p.value_heterogeneity), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity)) |>
  select(-p.value_heterogeneity.x, -p.value_heterogeneity.y)

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, trend_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_trend = ifelse(term == "Continuous", NA, p.value_trend), 
         p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
         p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(model1_cox_sd_danish, model2_cox_sd_danish, model3_cox_sd_danish, 
   model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish,
   model1_cox_sd_finnish, model2_cox_sd_finnish, model3_cox_sd_finnish, 
   model1_cox_quart_finnish, model2_cox_quart_finnish, model3_cox_quart_finnish,  
   model1_cox_sd_metanalysis, model2_cox_sd_metanalysis, model3_cox_sd_metanalysis, 
   model1_cox_quart_metanalysis, model2_cox_quart_metanalysis, model3_cox_quart_metanalysis, 
   bdd_cases_danish, bdd_cases_FMC, bdd_cases_FMCF, bdd_cases_MFH, 
   heterogeneity_tests, trend_tests)


# Sensitivity analysis 1 ----
## data prep ----
POPs_group_bis <- setdiff(POPs_group, "ΣPBDE")
POPs_group_quart_bis <- setdiff(POPs_group_quart, "ΣPBDE_quart")

bdd_cases_tot <- bdd |>
  filter (als == 1) |>
  select(study, als, follow_up_death, status_death, sex, diagnosis_age, baseline_age, marital_status_2cat, smoking_2cat, bmi, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
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

bdd_cases_red <- bdd |>
  filter (als == 1) |>
  filter(baseline_age>49) |>
  filter(follow_up<300) |>
  select(study, als, follow_up_death, status_death, sex, diagnosis_age, marital_status_2cat, smoking_2cat, bmi, baseline_age, follow_up, 
         "PCB_4", "PCB_DL", "PCB_NDL",  "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "Σchlordane") |>
  mutate(across(all_of(POPs_group_bis), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                                 labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
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

surv_obj_tot <- Surv(time = bdd_cases_tot$follow_up_death,                # set the outcomes
                     event = bdd_cases_tot$status_death)

surv_obj_red <- Surv(time = bdd_cases_red$follow_up_death,                # set the outcomes
                     event = bdd_cases_red$status_death)

## justification ----
plot_follow_up <- ggplot(bdd_cases_tot, aes(x = "", y = follow_up)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Follow-up (months)", title = "Distribution of the duration between baseline and diagnosis")

plot_baseline_age <- ggplot(bdd_cases_tot, aes(x = "", y = baseline_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Age at baseline (years)", title = "Distribution of the subjects age at baseline")

plot_diagnosis_age <- ggplot(bdd_cases_tot, aes(x = "", y = diagnosis_age)) +
  geom_violin(fill = "gray90", color = "gray50", width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = "black") +
  geom_jitter(aes(color = study), width = 0.2, alpha = 0.6, size = 1.8) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(x = NULL, y = "Age at diagnosis (years)", title = "Distribution of the subjects age at diagnosis")

plot_baseline_follow_up <- 
  ggplot(bdd_cases_tot) +
  aes(x = baseline_age, y = follow_up, colour = study) +
  geom_point(alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() + 
  geom_rect(aes(xmin = 49, xmax = 66, ymin = 0, ymax = 300),
            fill = NA, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Age at baseline (years)", 
       y = "Duration between baseline and diagnosis (months)")

plot_justif <- (plot_baseline_age + plot_follow_up + plot_diagnosis_age)/plot_baseline_follow_up + plot_layout()

## analysis ----
model2_cox_quart_sensi1 <- map_dfr(POPs_group_quart_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_tot) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_tot", 
    model = "adjusted", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})

model2_cox_quart_sensi1_red <- map_dfr(POPs_group_quart_bis, function(expl) {
  
  formula <- as.formula(paste("surv_obj_red ~", expl, "+ study +",            # set the formulas              
                              paste(covariates_finnish, collapse = " + ")))
  
  model_summary <- coxph(formula, data = bdd_cases_red) |> summary()  # run cox model
  coefs <- model_summary$coefficients
  tibble(                                                                       # creation of a table of results
    study = "sensi1_red", 
    model = "adjusted", 
    term = rownames(coefs),
    explanatory = expl, 
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"], 
    `p-value` = coefs[, "Pr(>|z|)"]) |>
    filter(str_starts(term, explanatory))                                       # remove the covariates results
})


## Heterogeneity tests ----
heterogeneity_sensi1_tot <- data.frame(explanatory = character(),
                                       study = factor(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (expl in POPs_group_quart_bis) {
  
  formula_raw <- as.formula(paste("surv_obj_tot ~ ", paste(c(covariates_finnish, "study"), collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_tot) 
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+ ", paste(c(covariates_finnish, "study"), collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_tot)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_sensi1_tot <- rbind(heterogeneity_sensi1_tot, 
                                    data.frame(explanatory = expl,
                                               study = "sensi1_tot",
                                               model = "adjusted",
                                               p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

heterogeneity_sensi1_red <- data.frame(explanatory = character(),
                                       study = factor(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (expl in POPs_group_quart_bis) {
  
  formula_raw <- as.formula(paste("surv_obj_red ~ ", paste(c(covariates_finnish, "study"), collapse = "+")))
  model_raw <- coxph(formula_raw, data = bdd_cases_red) 
  
  formula <- as.formula(paste("surv_obj_red ~", expl, "+ ", paste(c(covariates_finnish, "study"), collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_red)  
  
  anova <- anova(model_raw, model, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_sensi1_red <- rbind(heterogeneity_sensi1_red, 
                                    data.frame(explanatory = expl,
                                               study = "sensi1_red",
                                               model = "adjusted",
                                               p.value_heterogeneity = p.value_heterogeneity))
}
rm(expl, formula_raw, model_raw, formula, model, anova, p.value_heterogeneity)

heterogeneity_tests <- 
  bind_rows(heterogeneity_sensi1_tot, 
            heterogeneity_sensi1_red) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

## Trend tests ----
POPs_group_quart_med_bis <- setdiff(POPs_group_quart_med, "ΣPBDE_quart_med")

trend_sensi1_tot <- data.frame(explanatory = character(),
                               model = factor(), 
                               study = factor(),
                               p.value_trend = numeric(), 
                               stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med_bis) {
  
  formula <- as.formula(paste("surv_obj_tot ~", expl, "+ ", paste(c(covariates_finnish, "study"), collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_tot) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_sensi1_tot <- rbind(trend_sensi1_tot, 
                            data.frame(explanatory = expl,
                                       study = "sensi1_tot",
                                       model = "adjusted",
                                       p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

trend_sensi1_red <- data.frame(explanatory = character(),
                               model = factor(), 
                               study = factor(),
                               p.value_trend = numeric(), 
                               stringsAsFactors = FALSE)

for (expl in POPs_group_quart_med_bis) {
  
  formula <- as.formula(paste("surv_obj_red ~", expl, "+ ", paste(c(covariates_finnish, "study"), collapse = "+")))  
  model <- coxph(formula, data = bdd_cases_red) |> summary()
  p.value_trend <- model$coefficients[expl, "Pr(>|z|)"]
  
  trend_sensi1_red <- rbind(trend_sensi1_red, 
                            data.frame(explanatory = expl,
                                       study = "sensi1_red", 
                                       model = "adjusted",
                                       p.value_trend = p.value_trend))
}
rm(expl, model, formula, p.value_trend)

trend_tests <- 
  bind_rows(trend_sensi1_tot, trend_sensi1_red) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

rm(heterogeneity_sensi1_tot, heterogeneity_sensi1_red, 
   trend_sensi1_tot, trend_sensi1_red)


results_sensi1 <-       
  bind_rows(
    model2_cox_quart_sensi1, model2_cox_quart_sensi1_red) |>
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  select(study, model, explanatory, term, HR, lower_CI, upper_CI, "p-value") |>
  mutate(
    term = case_when(
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


results_sensi1 <- 
  left_join(results_sensi1, heterogeneity_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity)) 

results_sensi1 <- 
  left_join(results_sensi1, trend_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
         p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(POPs_group_bis, POPs_group_quart_bis, bdd_cases_tot, bdd_cases_red, surv_obj_tot, surv_obj_red, 
   model2_cox_quart_sensi1, model2_cox_quart_sensi1_red, 
   heterogeneity_tests, trend_tests)

## presentation of the results ----
quartile1_rows <- results_sensi1 |>
  distinct(study, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "", 
    "p.value_heterogeneity" = '', 
    "p.value_trend" = '')

sensi1_table <- results_sensi1 |>
  select(study, explanatory, term, HR, "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend") |>
  mutate(across(everything(), as.character))

sensi1_table <- 
  bind_rows(quartile1_rows, sensi1_table) |>
  mutate(`p-value` = str_replace(`p-value`, "1.00", ">0.99")) |>
  arrange(explanatory, term) |>
  pivot_wider(names_from = "study", values_from = c("HR", "95% CI", "p-value", "p.value_heterogeneity", "p.value_trend")) |>
  select(explanatory, term, contains("sensi1_tot"), contains("sensi1_red")) |>
  group_by(explanatory) |>
  mutate(p.value_heterogeneity_sensi1_tot = ifelse(term == 'quartile 1', p.value_heterogeneity_sensi1_tot[term == 'quartile 2'], ''), 
         p.value_trend_sensi1_tot = ifelse(term == 'quartile 1', p.value_trend_sensi1_tot[term == 'quartile 2'], ''),
         p.value_heterogeneity_sensi1_red = ifelse(term == 'quartile 1', p.value_heterogeneity_sensi1_red[term == 'quartile 2'], ''), 
         p.value_trend_sensi1_red = ifelse(term == 'quartile 1', p.value_trend_sensi1_red[term == 'quartile 2'], '')) |>
  ungroup() |>
  rename("HR" = "HR_sensi1_tot", "95% CI" = "95% CI_sensi1_tot", "p-value" = "p-value_sensi1_tot", "Heterogeneity test" = "p.value_heterogeneity_sensi1_tot", "Trend test" = "p.value_trend_sensi1_tot",
         "HR " = "HR_sensi1_red", "95% CI " = "95% CI_sensi1_red", "p-value " = "p-value_sensi1_red",  "Heterogeneity test " = "p.value_heterogeneity_sensi1_red", "Trend test " = "p.value_trend_sensi1_red") |>
  mutate(explanatory = factor(explanatory, levels = POPs_group_labels), 
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |>
  flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are adjusted for sex, age at diagnosis, smoking, BMI and marital status.
  3Estimated risk of ALS death when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.
  5Heterogeneity tests in outcome value across POP quartiles, adjusted for sex and age at diagnosis.
  6Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex and age at diagnosis.
  7Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at diagnosis, smoking, BMI and marital status.
  8Trend tests using continuous variables whose values corresponded to the quartile specific median POP levels, adjusted for sex, age at diagnosis, smoking, BMI and marital status.") |>
  add_header(
    "explanatory" = "Exposures", 
    term = "Quartiles",
    "HR" = "Total sample size (n=263)", "95% CI" = "Total sample size (n=263)", "p-value" = "Total sample size (n=263)",  "Heterogeneity test" = "Total sample size (n=263)",  "Trend test" = "Total sample size (n=263)",
    "HR " = "Reduced sample size (n=184)", "95% CI " = "Reduced sample size (n=184)", "p-value " = "Reduced sample size (n=184)",  "Heterogeneity test " = "Reduced sample size (n=184)",  "Trend test " = "Reduced sample size (n=184)") |>
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
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI and marital status. 
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
                            "Adjusted model" = "adjusted"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
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


### table POPs - ALS survival (qgcomp analysis) ----
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

### figure POPs - ALS survival (qgcomp analysis) ----
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

rm(qgcomp_boot_danish, qgcomp_noboot_danish, p, POPs_group_bis, pollutant_labels_bis)


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
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_finnish)) |> 
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI and marital status. 
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
    `p-value` = "", 
    "p.value_heterogeneity" = "")

POPs_quart_ALS_table_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
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
    "1All models are adjusted for age and sex. Adjusted models further account for smoking, BMI and marital status.
  2Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "POPs", 
    term = "Quartiles",
    "HR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", "Hetero-geneity test" = "Base Model",
    "HR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model",  "Hetero-geneity test " = "Adjusted Model",
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
POPs_sd_ALS_figure_finnish <- main_results_POPs_ALS_survival |>
  filter(study == "Finnish") |>
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
  facet_grid(cols = dplyr::vars(model), switch = "y") +  
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
    study = "finnish", 
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
  select(model, explanatory, term, HR, "95% CI", "p-value") |>
  filter(term == "Continuous") |>
  pivot_wider(names_from = "model", values_from = c("HR", "95% CI", "p-value")) |>
  select(explanatory, contains("base"), contains("adjusted"), contains("copollutant")) |>
  rename("HR" = "HR_base", "95% CI" = "95% CI_base", "p-value" = "p-value_base", 
         "HR " = "HR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p-value_adjusted", 
         " HR " = "HR_copollutant", " 95% CI " = "95% CI_copollutant", " p-value " = "p-value_copollutant") |>
  mutate(explanatory = fct_recode(explanatory, !!!POPs_group_labels_metanalysis)) |> 
  flextable() |>
  add_footer_lines(
    "1Base models were all adjusted for age and sex. Adjusted models all further account for smoking, BMI and marital status. 
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
  "1Base models were all adjusted for age and sex. Adjusted models all further account for smoking, BMI and marital status. 
  2Estimated risk of death after ALS diagnosis when pre-disease serum concentration of POPs compared to quartile 1.
  3CI: Confidence interval.") |>
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
         POPs_quart_ALS_figure_danish = POPs_quart_ALS_figure_danish, 
         plot_base_cox_gam_danish = plot_base_cox_gam_danish, 
         plot_adjusted_cox_gam_danish = plot_adjusted_cox_gam_danish, 
         plot_copollutant_cox_gam_danish = plot_copollutant_cox_gam_danish, 
         POPs_ALS_qgcomp_table_danish = POPs_ALS_qgcomp_table_danish,
         POPs_ALS_qgcomp_figure_danish = POPs_ALS_qgcomp_figure_danish), 
       finnish = list(
         covar_finnish = covar_finnish, 
         covar_ALS_table_finnish = covar_ALS_table_finnish,
         POPs_sd_ALS_table_finnish = POPs_sd_ALS_table_finnish, 
         POPs_quart_ALS_table_finnish = POPs_quart_ALS_table_finnish, 
         POPs_sd_ALS_figure_finnish = POPs_sd_ALS_figure_finnish, 
         POPs_quart_ALS_figure_finnish = POPs_quart_ALS_figure_finnish, 
         POPs_ALS_qgcomp_table_finnish = POPs_ALS_qgcomp_table_finnish, 
         POPs_ALS_qgcomp_figure_finnish = POPs_ALS_qgcomp_figure_finnish), 
       metanalysis = list(
         POPs_sd_ALS_table_metanalysis = POPs_sd_ALS_table_metanalysis, 
         POPs_quart_ALS_table_metanalysis = POPs_quart_ALS_table_metanalysis, 
         POPs_sd_ALS_figure_metanalysis = POPs_sd_ALS_figure_metanalysis, 
         POPs_quart_ALS_figure_metanalysis = POPs_quart_ALS_figure_metanalysis), 
       sensi1 = list(
         plot_justif = plot_justif, 
         results_sensi1 = results_sensi1, 
         sensi1_table = sensi1_table))

rm(main_results_POPs_ALS_survival, 
   covar_danish,
   POPs_sd_ALS_table_danish, 
   POPs_quart_ALS_table_danish, 
   POPs_sd_ALS_figure_danish, 
   POPs_quart_ALS_figure_danish, 
   plot_base_cox_gam_danish, 
   plot_adjusted_cox_gam_danish, 
   plot_copollutant_cox_gam_danish,
   POPs_ALS_qgcomp_table_danish, 
   POPs_ALS_qgcomp_figure_danish,
   
   covar_finnish,
   covar_ALS_table_finnish, 
   POPs_sd_ALS_table_finnish, 
   POPs_quart_ALS_table_finnish, 
   POPs_sd_ALS_figure_finnish, 
   POPs_quart_ALS_figure_finnish,
   POPs_ALS_qgcomp_table_finnish,
   POPs_ALS_qgcomp_figure_finnish, 
   
   POPs_sd_ALS_table_metanalysis, 
   POPs_quart_ALS_table_metanalysis, 
   POPs_sd_ALS_figure_metanalysis, 
   POPs_quart_ALS_figure_metanalysis, 
   
   covariates_danish, covariates_finnish, 
   surv_obj_danish, surv_obj_FMC, surv_obj_FMCF, surv_obj_MFH, 
   
   plot_justif, results_sensi1, sensi1_table)

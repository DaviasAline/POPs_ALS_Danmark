# Aline Davias
# April 29, 2025 
# Analysis of survival after ALS diagnosis depending on POPs levels  

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")

# Creation of cases specific datasets ----
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

surv_obj_danish <- Surv(time = bdd_cases_danish$follow_up_death,                # set the outcomes
                        event = bdd_cases_danish$status_death)

covariates_danish <- c("sex", "diagnosis_age", "smoking_2cat_i", "bmi", "marital_status_2cat_i")

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


## Cox-gam model ----
### Base  ----
POPs_group_labels_cox_gam <- set_names(
  POPs_group, 
  c("Most prevalent PCBs", "Dioxin-like PCBs","Non-dioxin-like PCBs", "HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"))
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
  x_label <- POPs_group_labels_cox_gam[var] 
  
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
  set_names(POPs_group_labels_cox_gam)


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
  x_label <- POPs_group_labels_cox_gam[var] 
  
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
  set_names(POPs_group_labels_cox_gam)
rm(POPs_group_labels_cox_gam)

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
rm(formula_danish)



# Assemblage main analyses ----
main_results_POPs_ALS_survival <-       
  bind_rows(
    model1_cox_sd_danish, model2_cox_sd_danish, model3_cox_sd_danish,
    model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish) |> 
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * se),
    upper_CI = exp(coef + 1.96 * se)) |>
  mutate(
    term = case_when(
      str_detect(term, "_sd") ~ "Continuous", 
      str_detect(term, "Q2") ~ "quartile 2",
      str_detect(term, "Q3") ~ "quartile 3",
      str_detect(term, "Q4") ~ "quartile 4"), 
    explanatory = gsub("_sd", "", explanatory), 
    explanatory = gsub("_quart", "", explanatory), 
    HR_raw = HR, 
    HR = as.numeric(sprintf("%.1f", HR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw` = `p-value`, 
    `p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, explanatory, term,  HR, HR_raw, `95% CI`, `p-value`, `p-value_raw`, `p-value_shape`, lower_CI, upper_CI)

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, heterogeneity_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_heterogeneity = ifelse(term == "Continuous" & study == "Danish", NA, p.value_heterogeneity), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
         p.value_heterogeneity = ifelse(p.value_heterogeneity == "1.00", ">0.99", p.value_heterogeneity))

main_results_POPs_ALS_survival <- 
  left_join(main_results_POPs_ALS_survival, trend_tests, 
            by = c("explanatory", "model", "study")) |>
  mutate(p.value_trend = ifelse(term == "Continuous", NA, p.value_trend), 
         p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")), 
         p.value_trend = ifelse(p.value_trend == "1.00", ">0.99", p.value_trend)) 

rm(model1_cox_sd_danish, 
   model2_cox_sd_danish, 
   model3_cox_sd_danish, 
   model1_cox_quart_danish, model2_cox_quart_danish, model3_cox_quart_danish,
   heterogeneity_tests, trend_tests)


# Sensitivity analysis 1 - mixture model ----
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
sensi1_lasso_sd_danish <- cv.glmnet(                                            # lasso + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_sd)                                           # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi1_lasso_sd_danish)
coef(sensi1_lasso_sd_danish, s = "lambda.min")                                   # lasso keeps only chlordane 
coef(sensi1_lasso_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi1_lasso_quart_danish <- cv.glmnet(                                          # lasso + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_quart)                                        # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi1_lasso_quart_danish)
coef(sensi1_lasso_quart_danish, s = "lambda.min")                                # lasso doesn't keep any pollutant when they are quartiles
coef(sensi1_lasso_quart_danish, s = "lambda.1se")  

## Ridge POPs selection ----
set.seed(1996)
sensi1_ridge_sd_danish <- cv.glmnet(                                             # ridge + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                                             

plot(sensi1_ridge_sd_danish)
coef(sensi1_ridge_sd_danish, s = "lambda.min")  
coef(sensi1_ridge_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi1_ridge_quart_danish <- cv.glmnet(                                          # ridge + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                                             

plot(sensi1_ridge_quart_danish)
coef(sensi1_ridge_quart_danish, s = "lambda.min")  
coef(sensi1_ridge_quart_danish, s = "lambda.1se")  

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
sensi1_elastic_net_sd_danish <- cv.glmnet(                                       # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.4,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                          

plot(sensi1_elastic_net_sd_danish)
coef(sensi1_elastic_net_sd_danish, s = "lambda.min")                             # Elastic net keeps HCB and Σchlordane
coef(sensi1_elastic_net_sd_danish, s = "lambda.1se")

set.seed(1996)
sensi1_elastic_net_quart_danish <- cv.glmnet(                                    # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.3,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                          

plot(sensi1_elastic_net_quart_danish)
coef(sensi1_elastic_net_quart_danish, s = "lambda.min")                          # elastic net doesn't keep any pollutant when they are quartiles
coef(sensi1_elastic_net_quart_danish, s = "lambda.1se")
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
sensi1_model3_cox_sd_elastic_net_danish <- tibble(                                     # creation of a table of results
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
sensi1_model3_cox_quart_elastic_net_danish <- tibble(                                  # creation of a table of results
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
coef(sensi1_elastic_net_sd_danish, s = "lambda.min")                            # Elastic net keeps HCB and Σchlordane
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


sensi1_cox_model_ERS_from_elastic_net_sd <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_sd + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi1_cox_model_ERS_from_elastic_net_sd)
coefs <- model_summary$coefficients
sensi1_cox_model_ERS_from_elastic_net_sd <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "ERS_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

sensi1_cox_model_ERS_from_elastic_net_quart <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_quart + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi1_cox_model_ERS_from_elastic_net_quart)
coefs <- model_summary$coefficients
sensi1_cox_model_ERS_from_elastic_net_quart <- tibble(                                  # creation of a table of results
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
results_sensi1 <-
  bind_rows(sensi1_model3_cox_sd_elastic_net_danish, 
            sensi1_model3_cox_quart_elastic_net_danish, 
            sensi1_cox_model_ERS_from_elastic_net_sd, 
            sensi1_cox_model_ERS_from_elastic_net_quart) |>
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

results_sensi1 <- left_join(results_sensi1, heterogeneity_copollutant_quart, 
                            by = c("study", "model", "explanatory"))

results_sensi1 <- 
  left_join(results_sensi1, trend_copollutant, 
            by = c("study", "model", "explanatory")) 

results_sensi1 <- results_sensi1 |>
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
   
   sensi1_cox_model_ERS_from_elastic_net_sd, 
   sensi1_cox_model_ERS_from_elastic_net_quart)

### table sd copollutant ----
POPs_sd_ALS_table_sensi1_danish <- results_sensi1 |>
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
  results_sensi1 |>
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

POPs_quart_ALS_table_sensi1_danish <- 
  results_sensi1 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi1_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi1_danish) |>
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
POPs_sd_ALS_figure_sensi1_danish <-
  results_sensi1 |>
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
POPs_quart_ALS_figure_sensi1_danish <-
  results_sensi1 |>
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
POPs_sd_ALS_table_sensi1_ERS_danish <- results_sensi1 |>
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
  results_sensi1 |>
  filter(model == "ERS_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") 

POPs_quart_ALS_table_sensi1_ERS_danish <- 
  results_sensi1 |>
  filter(study == "Danish") |>
  filter(model == "ERS_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi1_ERS_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi1_ERS_danish) |>
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
POPs_sd_ALS_figure_sensi1_ERS_danish <-
  results_sensi1 |>
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
POPs_quart_ALS_figure_sensi1_ERS_danish <-
  results_sensi1 |>
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
   sensi1_lasso_quart_danish, 
   sensi1_ridge_quart_danish, 
   sensi1_elastic_net_quart_danish, 
   sensi1_model3_cox_sd_elastic_net_danish, 
   sensi1_model3_cox_quart_elastic_net_danish)


# Sensitivity analysis 2 - elastic net on all the POPs (ungrouped) ----
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
sensi2_lasso_sd_danish <- cv.glmnet(                                             # lasso + cross validation to choose the lamda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_sd)                                           # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi2_lasso_sd_danish)
coef(sensi2_lasso_sd_danish, s = "lambda.min")                                   # lasso keeps only chlordane 
coef(sensi2_lasso_sd_danish, s = "lambda.1se")  

set.seed(1996)
sensi2_lasso_quart_danish <- cv.glmnet(                                          # lasso + cross validation to choose the lamda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 1,                                                                    # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10.
  penalty.factor = penalty_factor_quart)                                        # pre-set penalty factor because we want to force the model to select at least the covariates

plot(sensi2_lasso_quart_danish)
coef(sensi2_lasso_quart_danish, s = "lambda.min")                                # lasso doesn't keep any pollutant when they are quartiles
coef(sensi2_lasso_quart_danish, s = "lambda.1se")  


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
sensi2_elastic_net_sd_danish <- cv.glmnet(                                       # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_sd,                                                              # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.4,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_sd)                          

plot(sensi2_elastic_net_sd_danish)
coef(sensi2_elastic_net_sd_danish, s = "lambda.min")                             # Elastic net keeps HCB, β-HCH and Σchlordane
coef(sensi2_elastic_net_sd_danish, s = "lambda.1se")

set.seed(1996)
sensi2_elastic_net_quart_danish <- cv.glmnet(                                    # elastic net + cross validation to choose the lambda parameter
  x = X_matrix_quart,                                                           # matrice of explanatory variables 
  y = with(bdd_cases_danish_bis, Surv(follow_up_death, status_death)),          # survival outcome                   
  family = "cox",                                                               # cox regression 
  standardize = FALSE,                                                          # standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. 
  alpha = 0.3,                                                                  # alpha is for the elastic net mixing parameter 𝛼, with range 𝛼∈[0,1]. 𝛼=1 is lasso regression (default) and 𝛼=0 is ridge regression.
  type.measure = "deviance",                                                    # in cox models, we can choose between C and deviance. type.measure = "deviance" cv.glmnet choisira λ qui minimise la partial likelihood deviance (≈ maximise la log-vraisemblance).
  nfolds = 10,                                                                  # number of folds for the cross validation process. Default is 10. 
  penalty.factor = penalty_factor_quart)                          

plot(sensi2_elastic_net_quart_danish)
coef(sensi2_elastic_net_quart_danish, s = "lambda.min")                          # elastic net doesn't keep any pollutant when they are quartiles
coef(sensi2_elastic_net_quart_danish, s = "lambda.1se")
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
sensi2_model3_cox_sd_elastic_net_danish <- tibble(                                     # creation of a table of results
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
sensi2_model3_cox_quart_elastic_net_danish <- tibble(                                  # creation of a table of results
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
coef(sensi2_elastic_net_sd_danish, s = "lambda.min")                            # Elastic net keeps HCB and Σchlordane
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


sensi2_cox_model_ERS_from_elastic_net_sd <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_sd + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <-  summary(sensi2_cox_model_ERS_from_elastic_net_sd)
coefs <- model_summary$coefficients
sensi2_cox_model_ERS_from_elastic_net_sd <- tibble(                                     # creation of a table of results
  study = "Danish", 
  model = "ERS_sd", 
  term = rownames(coefs),
  explanatory = rownames(coefs),
  coef = coefs[, "coef"],
  se = coefs[, "se(coef)"], 
  `p-value` = coefs[, "Pr(>|z|)"]) |>
  filter(str_detect(term, "_sd"))

sensi2_cox_model_ERS_from_elastic_net_quart <- coxph(Surv(follow_up_death, status_death) ~ ERS_score_from_elastic_net_quart + sex + diagnosis_age + smoking_2cat_i + bmi + marital_status_2cat_i, data = bdd_cases_danish)
model_summary <- summary(sensi2_cox_model_ERS_from_elastic_net_quart)
coefs <- model_summary$coefficients
sensi2_cox_model_ERS_from_elastic_net_quart <- tibble(                                  # creation of a table of results
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
results_sensi2 <-
  bind_rows(sensi2_model3_cox_sd_elastic_net_danish, 
            sensi2_model3_cox_quart_elastic_net_danish,
            sensi2_cox_model_ERS_from_elastic_net_sd, 
            sensi2_cox_model_ERS_from_elastic_net_quart) |>
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


results_sensi2 <- left_join(results_sensi2, heterogeneity_copollutant_quart, 
                            by = c("study", "model", "explanatory"))

results_sensi2 <- 
  left_join(results_sensi2, trend_copollutant, 
            by = c("study", "model", "explanatory")) 

results_sensi2 <- results_sensi2 |>
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
   sensi2_model3_cox_sd_elastic_net_danish, 
   sensi2_model3_cox_quart_elastic_net_danish,
   sensi2_cox_model_ERS_from_elastic_net_sd, 
   sensi2_cox_model_ERS_from_elastic_net_quart)


### table sd copollutant ----
POPs_sd_ALS_table_sensi2_danish <- 
  results_sensi2 |>
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
  results_sensi2 |>
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

POPs_quart_ALS_table_sensi2_danish <- 
  results_sensi2 |>
  filter(study == "Danish") |>
  filter(model == "copollutant_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value", 
         "Heterogeneity test" = "p.value_heterogeneity", 
         "Trend test" = "p.value_trend") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi2_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi2_danish) |>
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
POPs_sd_ALS_figure_sensi2_danish <-
  results_sensi2 |>
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
POPs_quart_ALS_figure_sensi2_danish <-
  results_sensi2 |>
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
   sensi2_lasso_quart_danish, 
   #sensi2_model_ridge_quart_danish, 
   sensi2_elastic_net_quart_danish)


### table sd ERS ----
POPs_sd_ALS_table_sensi2_ERS_danish <- results_sensi2 |>
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
  results_sensi2 |>
  filter(model == "ERS_quart") |>
  distinct(model, explanatory) |>
  mutate(
    term = "quartile 1",
    HR = "-",
    "95% CI" = "-",
    `p-value` = "") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") 

POPs_quart_ALS_table_sensi2_ERS_danish <- 
  results_sensi2 |>
  filter(study == "Danish") |>
  filter(model == "ERS_quart") |>
  select("explanatory", "term", "HR", "95% CI", "p-value") |>
  mutate(across(everything(), as.character))

POPs_quart_ALS_table_sensi2_ERS_danish <- 
  bind_rows(quartile1_rows, POPs_quart_ALS_table_sensi2_ERS_danish) |>
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
POPs_sd_ALS_figure_sensi2_ERS_danish <-
  results_sensi2 |>
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
POPs_quart_ALS_figure_sensi2_ERS_danish <-
  results_sensi2 |>
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


# Sensitivity analysis 3 - follow-up duration ----
# investigation du probleme de follow-up duration très courte par rapport à la réalité (3-5 ans)
sensi3_figure <- bdd_cases_danish |>
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

sensi3_table <- bdd_cases_danish |> 
      select(follow_up_death, status_death) |>
      tbl_summary(by = status_death) |>
      add_overall()

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
  ggplot(aes(x = explanatory, y = HR_raw, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
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
  ggplot(aes(x = term, y = HR_raw, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
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
         qgcomp = list(qgcomp_boot_danish = qgcomp_boot_danish, 
                       qgcomp_noboot_danish = qgcomp_noboot_danish, 
                       qgcomp_positive_boot_danish = qgcomp_positive_boot_danish, 
                       qgcomp_positive_noboot_danish = qgcomp_positive_noboot_danish, 
                       POPs_ALS_qgcomp_table_danish = POPs_ALS_qgcomp_table_danish,
                       POPs_ALS_qgcomp_figure_danish = POPs_ALS_qgcomp_figure_danish, 
                       POPs_ALS_qgcomp_positive_table_danish = POPs_ALS_qgcomp_positive_table_danish,
                       POPs_ALS_qgcomp_positive_figure_danish = POPs_ALS_qgcomp_positive_figure_danish)), 
       sensi1 = list(
         sensi1_lasso_sd_danish = sensi1_lasso_sd_danish, 
         sensi1_ridge_sd_danish = sensi1_ridge_sd_danish, 
         sensi1_elastic_net_sd_danish = sensi1_elastic_net_sd_danish, 
         results_sensi1 = results_sensi1, 
         POPs_sd_ALS_figure_sensi1_danish = POPs_sd_ALS_figure_sensi1_danish, 
         POPs_quart_ALS_figure_sensi1_danish = POPs_quart_ALS_figure_sensi1_danish, 
         POPs_sd_ALS_table_sensi1_danish = POPs_sd_ALS_table_sensi1_danish, 
         POPs_quart_ALS_table_sensi1_danish = POPs_quart_ALS_table_sensi1_danish, 
         POPs_sd_ALS_figure_sensi1_ERS_danish = POPs_sd_ALS_figure_sensi1_ERS_danish, 
         POPs_quart_ALS_figure_sensi1_ERS_danish = POPs_quart_ALS_figure_sensi1_ERS_danish, 
         POPs_sd_ALS_table_sensi1_ERS_danish = POPs_sd_ALS_table_sensi1_ERS_danish, 
         POPs_quart_ALS_table_sensi1_ERS_danish = POPs_quart_ALS_table_sensi1_ERS_danish), 
       sensi2 = list(
         sensi2_lasso_sd_danish = sensi2_lasso_sd_danish, 
         # sensi2_ridge_sd_danish = sensi2_ridge_sd_danish, 
         sensi2_elastic_net_sd_danish = sensi2_elastic_net_sd_danish, 
         results_sensi2 = results_sensi2, 
         POPs_sd_ALS_figure_sensi2_danish = POPs_sd_ALS_figure_sensi2_danish, 
         POPs_quart_ALS_figure_sensi2_danish = POPs_quart_ALS_figure_sensi2_danish, 
         POPs_sd_ALS_table_sensi2_danish = POPs_sd_ALS_table_sensi2_danish, 
         POPs_quart_ALS_table_sensi2_danish = POPs_quart_ALS_table_sensi2_danish, 
         POPs_sd_ALS_figure_sensi2_ERS_danish = POPs_sd_ALS_figure_sensi2_ERS_danish, 
         POPs_quart_ALS_figure_sensi2_ERS_danish = POPs_quart_ALS_figure_sensi2_ERS_danish, 
         POPs_sd_ALS_table_sensi2_ERS_danish = POPs_sd_ALS_table_sensi2_ERS_danish, 
         POPs_quart_ALS_table_sensi2_ERS_danish = POPs_quart_ALS_table_sensi2_ERS_danish), 
       sensi3 = list(
         sensi3_figure = sensi3_figure, 
         sensi3_table = sensi3_table))

rm(bdd_cases_danish, 
   main_results_POPs_ALS_survival, 
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
   POPs_ALS_qgcomp_positive_table_danish, 
   POPs_ALS_qgcomp_positive_figure_danish,
   qgcomp_boot_danish, 
   qgcomp_noboot_danish, 
   qgcomp_positive_boot_danish, 
   qgcomp_positive_noboot_danish, 
   
   covariates_danish, 
   surv_obj_danish, 
   
   results_sensi1, 
   sensi1_lasso_sd_danish, 
   sensi1_ridge_sd_danish, 
   sensi1_elastic_net_sd_danish, 
   POPs_sd_ALS_figure_sensi1_danish, 
   POPs_quart_ALS_figure_sensi1_danish, 
   POPs_sd_ALS_table_sensi1_danish, 
   POPs_quart_ALS_table_sensi1_danish, 
   POPs_sd_ALS_figure_sensi1_ERS_danish, 
   POPs_quart_ALS_figure_sensi1_ERS_danish, 
   POPs_sd_ALS_table_sensi1_ERS_danish, 
   POPs_quart_ALS_table_sensi1_ERS_danish,
   
   results_sensi2, 
   sensi2_lasso_sd_danish, 
   #sensi2_ridge_sd_danish, 
   sensi2_elastic_net_sd_danish, 
   POPs_sd_ALS_figure_sensi2_danish, 
   POPs_quart_ALS_figure_sensi2_danish, 
   POPs_sd_ALS_table_sensi2_danish, 
   POPs_quart_ALS_table_sensi2_danish, 
   POPs_sd_ALS_figure_sensi2_ERS_danish, 
   POPs_quart_ALS_figure_sensi2_ERS_danish, 
   POPs_sd_ALS_table_sensi2_ERS_danish, 
   POPs_quart_ALS_table_sensi2_ERS_danish,
   
   sensi3_figure, 
   sensi3_table)


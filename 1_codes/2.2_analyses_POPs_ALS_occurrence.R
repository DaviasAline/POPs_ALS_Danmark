# Aline Davias
# 03/02/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.1_analyses_descriptive.R")

covariates <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'fS_Kol', 'marital_status_2cat_i', 'education_i')
covariates_danish <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'fS_Kol', 'marital_status_2cat_i', 'education_i')
covariates_finnish <- c("marital_status_2cat", 'smoking_2cat', 'bmi', 'fS_Kol')     # education removed because missing in one finnish cohort 

# effects of the covariates on ALS ----
covar <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
  select(als, all_of(covariates)) |>
  tbl_uvregression(
    y = als,
    method = glm,
    method.args = list(family = binomial), 
    exponentiate = TRUE, 
    estimate_fun = label_number(accuracy = 0.1, decimal.mark = "."),
    pvalue_fun = custom_pvalue_fun)|>
  bold_labels(), 
  tbl_2 = clogit(als ~ sex + baseline_age +
           smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tbl_regression(exponentiate = TRUE, 
                 estimate_fun = label_number(accuracy = .1, decimal.mark = "."),
                 pvalue_fun = custom_pvalue_fun) |>
  bold_labels()), 
  tab_spanner = c("**Univariate**", "**Adjusted**"))

# main analysis ----
### model 1 ----
#### quartiles ----
model1_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quart <- rbind(model1_quart, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model1_quart <- model1_quart |> 
  mutate(
    df = case_when(
      grepl("_quartQ2", df) ~ "Quartile 2",
      grepl("_quartQ3", df) ~ "Quartile 3",
      grepl("_quartQ4", df) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base_quart") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


#### gams ----
model1_gam <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### model 2 ----
# matched on sex and age, adjusted on for smoking_2cat_i, BMI, serum total fS_Kol, marital status and education

#### quartiles ----
model2_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quart <- rbind(model2_quart, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model2_quart <- model2_quart |> 
  mutate(
    df = case_when(
      grepl("_quartQ2", df) ~ "Quartile 2",
      grepl("_quartQ3", df) ~ "Quartile 3",
      grepl("_quartQ4", df) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted_quart") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### gams ----
model2_gam <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### models 3 ----
#### quartiles ----
##### pollutant of interest as quartile while the others are s() transformed in a gam model 
model3_quart_PCB_DL <- 
  gam(als ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_PCB_DL <- model3_quart_PCB_DL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_PCB_NDL <- 
  gam(als ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_PCB_NDL <- model3_quart_PCB_NDL$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_HCB <- 
  gam(als ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_HCB <- model3_quart_HCB$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣDDT <- 
  gam(als ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_ΣDDT <- model3_quart_ΣDDT$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_β_HCH <- 
  gam(als ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_β_HCH <- model3_quart_β_HCH$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_Σchlordane <- 
  gam(als ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_Σchlordane <- model3_quart_Σchlordane$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart_ΣPBDE <- 
  gam(als ~ 
        ΣPBDE_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
      family = binomial, 
      method = 'REML',
      data = bdd_danish) |>
  summary()
model3_quart_ΣPBDE <- model3_quart_ΣPBDE$p.table |>
  as.data.frame() |>
  rownames_to_column("variable") 

model3_quart <- bind_rows(
  model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane, model3_quart_ΣPBDE) |>
  filter(grepl("quart", variable)) |>
  mutate(
    df = case_when(
      grepl("_quartQ2", variable) ~ "Quartile 2",
      grepl("_quartQ3", variable) ~ "Quartile 3",
      grepl("_quartQ4", variable) ~ "Quartile 4",
      TRUE ~ NA_character_),
    variable = gsub("_quartQ2", "", variable), 
    variable = gsub("_quartQ3", "", variable), 
    variable = gsub("_quartQ4", "", variable), 
    model = "copollutant_quart", 
    OR = exp(Estimate), 
    lower_CI = exp(Estimate - 1.96 * `Std. Error`), 
    upper_CI = exp(Estimate + 1.96 * `Std. Error`), 
    p.value = `Pr(>|z|)`) |> 
  select(model, variable, df, OR, lower_CI, upper_CI, p.value)

rm(model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_HCB, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane, model3_quart_ΣPBDE)

#### qgcomp ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
set.seed(1996)
qgcomp_boot_danish <- 
  qgcomp.glm.boot(
    f = as.formula(paste("als ~", paste(c(POPs_group_bis, covariates_danish), collapse = " + "))),
    expnms = POPs_group_bis,                                                    # pollutants of interest
    data = bdd_danish, 
    family = binomial(), 
    q = 4, 
    B = 1000,                                                                   # nb of boostrap
    seed = 1996, 
    rr = FALSE)                                                                 # rr=FALSE to allow estimation of ORs when using qgcomp.glm.boot
#qgcomp_boot_danish$pos.weights

qgcomp_noboot_danish <-                                                         
  qgcomp.glm.noboot(
    f = as.formula(paste("als ~", paste(c(POPs_group_bis, covariates_danish), collapse = " + "))),                                                   # formula
    data = bdd_danish, 
    q = 4,                                                                      # number of quantiles
    expnms = POPs_group_bis)   
rm(POPs_group_bis)

### heterogeneity tests ----
#### model 1 quartile ----
heterogeneity_base_quart <- data.frame(variable = character(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart <- rbind(heterogeneity_base_quart, 
                                    data.frame(variable = var,
                                               model = "base_quart",
                                               p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

#### model 2 quartile ----
heterogeneity_adjusted_quart <- data.frame(variable = character(),
                                           model = factor(), 
                                           p.value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(variable = var,
                                                   model = "adjusted_quart",
                                                   p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)


#### model 3 quartile ----
test_1 <- gam(als ~ 
                s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                PCB_DL_quart + 
                s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_PCB_DL <- tibble(variable = "PCB_DL", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                PCB_NDL_quart + 
                s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_PCB_NDL <- tibble(variable = "PCB_NDL", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(PCB_DL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                OCP_HCB_quart + 
                s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_OCP_HCB <- tibble(variable = "OCP_HCB", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(PCB_DL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                ΣDDT_quart + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_ΣDDT <- tibble(variable = "ΣDDT", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                OCP_β_HCH_quart + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_OCP_β_HCH <- tibble(variable = "OCP_β_HCH", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                Σchlordane_quart + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_Σchlordane <- tibble(variable = "Σchlordane", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

test_1 <- gam(als ~ 
                s(PCB_DL) + s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)
test_2 <- gam(als ~ 
                ΣPBDE_quart + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish)

anova <- anova(test_1, test_2, test = "Chisq")
p.value_heterogeneity_ΣPBDE <- tibble(variable = "ΣPBDE", p.value_heterogeneity = anova$`Pr(>Chi)`[2])

heterogeneity_copollutant <- bind_rows(p.value_heterogeneity_PCB_DL, 
                                               p.value_heterogeneity_PCB_NDL, 
                                               p.value_heterogeneity_OCP_HCB, 
                                               p.value_heterogeneity_ΣDDT, 
                                               p.value_heterogeneity_OCP_β_HCH, 
                                               p.value_heterogeneity_Σchlordane, 
                                               p.value_heterogeneity_ΣPBDE) |>
  mutate(model = "copollutant_quart")

rm(test_1, test_2, anova, 
   p.value_heterogeneity_PCB_DL, 
   p.value_heterogeneity_PCB_NDL, 
   p.value_heterogeneity_OCP_HCB, 
   p.value_heterogeneity_ΣDDT, 
   p.value_heterogeneity_OCP_β_HCH, 
   p.value_heterogeneity_Σchlordane, 
   p.value_heterogeneity_ΣPBDE)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_quart, 
            heterogeneity_adjusted_quart, 
            heterogeneity_copollutant) |>
  mutate(variable = gsub("_quart", "", variable))

### trend tests ----
#### model 1 quartile ----
trend_base <- data.frame(variable = character(),
                         model = factor(), 
                         p.value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish) |>
    summary() 
  p.value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base <- rbind(trend_base, 
                      data.frame(variable = var,
                                 model = "base_quart",
                                 p.value_trend = p.value_trend))
}
rm(var, test, formula, p.value_trend)

#### model 2 quartile ----
trend_adjusted <- data.frame(variable = character(),
                             model = factor(), 
                             p.value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  test <- clogit(formula, data = bdd_danish) |> summary()
  p.value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(variable = var,
                                     model = "adjusted_quart",
                                     p.value_trend = p.value_trend))
}
rm(var, test, formula, p.value_trend)

#### model 3 quartile ----
test <- gam(als ~ 
                PCB_DL_quart_med + 
                s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_PCB_DL <- test$p.table["PCB_DL_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                PCB_NDL_quart_med + 
                s(PCB_DL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_PCB_NDL <- test$p.table["PCB_NDL_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                OCP_HCB_quart_med + 
                s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_OCP_HCB <- test$p.table["OCP_HCB_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                ΣDDT_quart_med + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_ΣDDT <- test$p.table["ΣDDT_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                OCP_β_HCH_quart_med + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_OCP_β_HCH <- test$p.table["OCP_β_HCH_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                Σchlordane_quart_med + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(ΣPBDE) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_Σchlordane <- test$p.table["Σchlordane_quart_med", "Pr(>|z|)"]

test <- gam(als ~ 
                ΣPBDE_quart_med + 
                s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
                sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
              family = binomial, 
              method = 'ML',
              data = bdd_danish) |>
  summary()
p.value_trend_ΣPBDE <- test$p.table["ΣPBDE_quart_med", "Pr(>|z|)"]

trend_copollutant <- 
  data.frame(variable = c("PCB_DL_quart_med", "PCB_NDL_quart_med", "OCP_HCB_quart_med", 
                          "ΣDDT_quart_med", "OCP_β_HCH_quart_med", "Σchlordane_quart_med", 
                          "ΣPBDE_quart_med" ),
           model = "copollutant_quart",
           p.value_trend = c(p.value_trend_PCB_DL, 
                             p.value_trend_PCB_NDL, 
                             p.value_trend_OCP_HCB, 
                             p.value_trend_ΣDDT, 
                             p.value_trend_OCP_β_HCH, 
                             p.value_trend_Σchlordane, 
                             p.value_trend_ΣPBDE))

rm(test, 
   p.value_trend_PCB_DL, 
   p.value_trend_PCB_NDL, 
   p.value_trend_OCP_HCB, 
   p.value_trend_ΣDDT, 
   p.value_trend_OCP_β_HCH, 
   p.value_trend_Σchlordane, 
   p.value_trend_ΣPBDE)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted, trend_copollutant) |>
  mutate(variable = gsub("_quart_med", "", variable))

### metanalysis (quartiles) ----
run_clogit <- function(formula, data) {
  model <- clogit(formula, data = data)
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"])
}
POPs_group_metanalysis <- setdiff(POPs_group, "ΣPBDE")                          # we don't include ΣPBDE in the metanalysis because low levels 
POPs_group_metanalysis_quart <- paste0(POPs_group_metanalysis, "_quart")

#### Base model ----
metanalysis_base_quart <- map_dfr(POPs_group_metanalysis_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))                # base formula: matched, not ajstuded
  
  bdd_danish <- bdd |>                                                     
    filter(study == "Danish") |>                                                # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  bdd_finnish_MFH <- bdd |> 
    filter(study == "MFH") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    danish = run_clogit(formula, bdd_danish),
    finnish_FMC = run_clogit(formula, bdd_finnish_FMC),
    finnish_FMCF = run_clogit(formula, bdd_finnish_FMCF), 
    finnish_MFH = run_clogit(formula, bdd_finnish_MFH)) |>
    bind_rows(.id = "dataset") |>
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "Quartile 2",
             str_detect(term, "Q3") ~ "Quartile 3",
             str_detect(term, "Q4") ~ "Quartile 4",
             TRUE ~ NA_character_)) 
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        OR = exp(as.numeric(rma_fit$beta)),
        lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
        upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
        `p-value` = as.numeric(rma_fit$pval), 
        p.value_heterogeneity = as.numeric(rma_fit$QEp))
    }) |> 
    ungroup() |> 
    mutate(model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

#### Adjusted model ----
metanalysis_adjusted_quart <- map_dfr(POPs_group_metanalysis_quart, function(expl) {
  formula_educ <- 
    as.formula(paste("als ~", expl, 
                     "+ strata(match) + marital_status_2cat + smoking_2cat + bmi + fS_Kol + education"))
  formula_no_educ <- 
    as.formula(paste("als ~", expl,                                             # education mot available in one finnish cohort
                     "+ strata(match) + marital_status_2cat + smoking_2cat + bmi + fS_Kol")) 
  
  bdd_danish <- bdd |>                                                     
    filter(study == "Danish") |>                                                # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  bdd_finnish_MFH <- bdd |> 
    filter(study == "MFH") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metanalysis), ~ factor(ntile(.x, 4),        # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    danish = run_clogit(formula_educ, bdd_danish),
    finnish_FMC = run_clogit(formula_no_educ, bdd_finnish_FMC),                 # no education data for this cohort 
    finnish_FMCF = run_clogit(formula_educ, bdd_finnish_FMCF)
    # , 
    # finnish_MFH = run_clogit(formula_educ, bdd_finnish_MFH)
    ) |>
    bind_rows(.id = "dataset") |>
    mutate(var = se^2, 
           explanatory = expl,
           term = case_when(
             str_detect(term, "Q2") ~ "Quartile 2",
             str_detect(term, "Q3") ~ "Quartile 3",
             str_detect(term, "Q4") ~ "Quartile 4",
             TRUE ~ NA_character_)) |>
    filter(str_detect(term, "Quartile"))                                        # filter to remove the covariates results
  
  meta_results <- results |>                                                    # run metanalyse (one per quartile per expl var)
    group_by(explanatory, term) |> 
    group_modify(~ {
      rma_fit <- rma(yi = .x$coef, vi = .x$var, method = "DL")
      tibble(                                                                   # results table creation 
        OR = exp(as.numeric(rma_fit$beta)),
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

#### Copollutant model ----
POPs_group_matanalysis <- setdiff(POPs_group, "ΣPBDE")  

bdd_metanalysis_danish <- bdd |>                                                     
  filter(study == "Danish") |>                                                  # creation of one dataset per finnish cohort
  mutate(across(all_of(POPs_group_matanalysis), ~ factor(ntile(.x, 4),          # creation of quartiles cohort specific                      
                                                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) 

bdd_metanalysis_FMC <- bdd |>                                                     
  filter(study == "FMC") |>                                                     # creation of one dataset per finnish cohort
  mutate(across(all_of(POPs_group_matanalysis), ~ factor(ntile(.x, 4),          # creation of quartiles cohort specific                      
                                                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) 
bdd_metanalysis_FMCF <- bdd |> 
  filter(study == "FMCF") |>                                                    # creation of one dataset per finnish cohort
  mutate(across(all_of(POPs_group_matanalysis), ~ factor(ntile(.x, 4),          # creation of quartiles cohort specific    
                                                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart"))
bdd_metanalysis_MFH <- bdd |> 
  filter(study == "MFH") |>                                                     # creation of one dataset per finnish cohort
  mutate(across(all_of(POPs_group_matanalysis), ~ factor(ntile(.x, 4),          # creation of quartiles cohort specific    
                                                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart"))
rm(POPs_group_matanalysis)

##### danish ----
model3_quart_PCB_DL_danish <- 
  gam(als ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_PCB_DL_danish <- model3_quart_PCB_DL_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_PCB_NDL_danish <- 
  gam(als ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_PCB_NDL_danish <- model3_quart_PCB_NDL_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_HCB_danish <- 
  gam(als ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_HCB_danish <- model3_quart_HCB_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_ΣDDT_danish <- 
  gam(als ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_ΣDDT_danish <- model3_quart_ΣDDT_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_β_HCH_danish <- 
  gam(als ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_β_HCH_danish <- model3_quart_β_HCH_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_Σchlordane_danish <- 
  gam(als ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) +
        sex + baseline_age + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_danish) |>
  summary()
model3_quart_Σchlordane_danish <- model3_quart_Σchlordane_danish$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

##### FMC ----
model3_quart_PCB_DL_FMC <- 
  gam(als ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_PCB_DL_FMC <- model3_quart_PCB_DL_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_PCB_NDL_FMC <- 
  gam(als ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_PCB_NDL_FMC <- model3_quart_PCB_NDL_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_HCB_FMC <- 
  gam(als ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_HCB_FMC <- model3_quart_HCB_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_ΣDDT_FMC <- 
  gam(als ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_ΣDDT_FMC <- model3_quart_ΣDDT_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_β_HCH_FMC <- 
  gam(als ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_β_HCH_FMC <- model3_quart_β_HCH_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_Σchlordane_FMC <- 
  gam(als ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat, # no education data forFMC cohort
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMC) |>
  summary()
model3_quart_Σchlordane_FMC <- model3_quart_Σchlordane_FMC$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

##### FMCF ----
model3_quart_PCB_DL_FMCF <- 
  gam(als ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_PCB_DL_FMCF <- model3_quart_PCB_DL_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_PCB_NDL_FMCF <- 
  gam(als ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_PCB_NDL_FMCF <- model3_quart_PCB_NDL_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_HCB_FMCF <- 
  gam(als ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_HCB_FMCF <- model3_quart_HCB_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_ΣDDT_FMCF <- 
  gam(als ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_ΣDDT_FMCF <- model3_quart_ΣDDT_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_β_HCH_FMCF <- 
  gam(als ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_β_HCH_FMCF <- model3_quart_β_HCH_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_Σchlordane_FMCF <- 
  gam(als ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_FMCF) |>
  summary()
model3_quart_Σchlordane_FMCF <- model3_quart_Σchlordane_FMCF$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

##### MFH -----
model3_quart_PCB_DL_MFH <- 
  gam(als ~ 
        PCB_DL_quart + 
        s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_PCB_DL_MFH <- model3_quart_PCB_DL_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_PCB_NDL_MFH <- 
  gam(als ~ 
        PCB_NDL_quart + 
        s(PCB_DL)  + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_PCB_NDL_MFH <- model3_quart_PCB_NDL_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_HCB_MFH <- 
  gam(als ~ 
        OCP_HCB_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(ΣDDT) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_HCB_MFH <- model3_quart_HCB_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_ΣDDT_MFH <- 
  gam(als ~ 
        ΣDDT_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(OCP_β_HCH) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education, 
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_ΣDDT_MFH <- model3_quart_ΣDDT_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_β_HCH_MFH <- 
  gam(als ~ 
        OCP_β_HCH_quart +
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(Σchlordane) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_β_HCH_MFH <- model3_quart_β_HCH_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

model3_quart_Σchlordane_MFH <- 
  gam(als ~ 
        Σchlordane_quart + 
        s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + s(OCP_β_HCH) +
        sex + baseline_age + thawed + level_urbanization + smoking_2cat + bmi + fS_Kol + marital_status_2cat + education,
      family = binomial, 
      method = 'REML',
      data = bdd_metanalysis_MFH) |>
  summary()
model3_quart_Σchlordane_MFH <- model3_quart_Σchlordane_MFH$p.table |>
  as.data.frame() |>
  rownames_to_column("explanatory") 

metanalysis_copollutant_quart <- 
  bind_rows(
    model3_quart_PCB_DL_danish, model3_quart_PCB_NDL_danish, model3_quart_HCB_danish, model3_quart_ΣDDT_danish, model3_quart_β_HCH_danish, model3_quart_Σchlordane_danish, 
    model3_quart_PCB_DL_FMC, model3_quart_PCB_NDL_FMC, model3_quart_HCB_FMC, model3_quart_ΣDDT_FMC, model3_quart_β_HCH_FMC, model3_quart_Σchlordane_FMC,  
    model3_quart_PCB_DL_FMCF, model3_quart_PCB_NDL_FMCF, model3_quart_HCB_FMCF, model3_quart_ΣDDT_FMCF, model3_quart_β_HCH_FMCF, model3_quart_Σchlordane_FMCF
    # , 
    # model3_quart_PCB_DL_MFH, model3_quart_PCB_NDL_MFH, model3_quart_HCB_MFH, model3_quart_ΣDDT_MFH, model3_quart_β_HCH_MFH, model3_quart_Σchlordane_MFH
  )

metanalysis_copollutant_quart <- metanalysis_copollutant_quart |>
  filter(grepl("quart", explanatory)) |>
  mutate(
    term = case_when(
      grepl("_quartQ2", explanatory) ~ "Quartile 2",
      grepl("_quartQ3", explanatory) ~ "Quartile 3",
      grepl("_quartQ4", explanatory) ~ "Quartile 4",
      TRUE ~ NA_character_),
    explanatory = gsub("_quartQ2", "", explanatory), 
    explanatory = gsub("_quartQ3", "", explanatory), 
    explanatory = gsub("_quartQ4", "", explanatory), 
    var = `Std. Error`^2) 

metanalysis_copollutant_quart <- metanalysis_copollutant_quart |>               # run metanalyse (one per quartile per expl var)
  group_by(explanatory, term) |> 
  group_modify(~ {
    rma_fit <- rma(yi = .x$Estimate, vi = .x$var, method = "DL")
    tibble(                                                                     # results table creation 
      OR = exp(as.numeric(rma_fit$beta)),
      lower_CI = exp(as.numeric(rma_fit$beta) - 1.96 * as.numeric(rma_fit$se)),
      upper_CI = exp(as.numeric(rma_fit$beta) + 1.96 * as.numeric(rma_fit$se)),
      `p-value` = as.numeric(rma_fit$pval),
      p.value_heterogeneity = as.numeric(rma_fit$QEp))
  }) |> 
  ungroup() |> 
  mutate(model = "copollutant") |> 
  relocate(model, explanatory, term)

rm(model3_quart_PCB_DL_danish, model3_quart_PCB_NDL_danish, model3_quart_HCB_danish, model3_quart_ΣDDT_danish, model3_quart_β_HCH_danish, model3_quart_Σchlordane_danish, 
   model3_quart_PCB_DL_FMC, model3_quart_PCB_NDL_FMC, model3_quart_HCB_FMC, model3_quart_ΣDDT_FMC, model3_quart_β_HCH_FMC, model3_quart_Σchlordane_FMC,  
   model3_quart_PCB_DL_FMCF, model3_quart_PCB_NDL_FMCF, model3_quart_HCB_FMCF, model3_quart_ΣDDT_FMCF, model3_quart_β_HCH_FMCF, model3_quart_Σchlordane_FMCF, 
   model3_quart_PCB_DL_MFH, model3_quart_PCB_NDL_MFH, model3_quart_HCB_MFH, model3_quart_ΣDDT_MFH, model3_quart_β_HCH_MFH, model3_quart_Σchlordane_MFH, 
   bdd_metanalysis_danish, bdd_metanalysis_FMC, bdd_metanalysis_FMCF, bdd_metanalysis_MFH)

metanalysis_quart <- bind_rows(metanalysis_base_quart, metanalysis_adjusted_quart, metanalysis_copollutant_quart) |> 
  mutate(explanatory = gsub("_quart", "", explanatory), 
         OR = sprintf("%.1f", OR),
         lower_CI = sprintf("%.1f", lower_CI),
         upper_CI = sprintf("%.1f", upper_CI),
         `p-value_raw` = `p-value`, 
         `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
         `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = ''),
         `p-value_heterogeneity` = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
         `p-value_heterogeneity` = ifelse(`p-value_heterogeneity` == "1.00", ">0.99", `p-value_heterogeneity`)) |>
  select(model,
         explanatory, 
         term,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p-value"), 
         lower_CI, upper_CI) 

rm(metanalysis_base_quart, metanalysis_adjusted_quart, metanalysis_copollutant_quart, run_clogit, POPs_group_metanalysis, POPs_group_metanalysis_quart)

### merging the main results ----
main_results <- bind_rows(model1_quart, 
                          model2_quart,
                          model3_quart) |> 
  mutate(variable = gsub("_quart", "", variable), 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  arrange(variable) |>
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value"), 
         lower_CI, upper_CI) 

main_results <- left_join(main_results, heterogeneity_tests, by = c("variable", "model"))
main_results <- left_join(main_results, trend_tests, by = c("variable", "model"))
main_results <- main_results |>
  mutate(
    p.value_heterogeneity = ifelse(p.value_heterogeneity < 0.01, "<0.01", number(p.value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p.value_trend = ifelse(p.value_trend < 0.01, "<0.01", number(p.value_trend, accuracy = 0.01, decimal.mark = ".")))


results_quart <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI, - p.value_raw) |>
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart"), 
         contains("copollutant_quart")) 
colnames(results_quart) <- gsub('_quart', '', colnames(results_quart))

rm(model1_quart, 
   model2_quart, 
   model3_quart, 
   heterogeneity_base_quart, heterogeneity_adjusted_quart, heterogeneity_copollutant, 
   trend_base, trend_adjusted, trend_copollutant, 
   heterogeneity_tests, trend_tests)

# sensitivity analyses pollutants not summed ----
### quartiles ----
model1_quart_not_summed <- data.frame(variable = character(),
                                      df = integer(),
                                      OR = numeric(),
                                      lower_CI = numeric(),
                                      upper_CI = numeric(),
                                      p_value = numeric(),
                                      stringsAsFactors = FALSE)

for (var in POPs_included_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quart_not_summed <- rbind(model1_quart_not_summed, data.frame(variable = var,
                                                                       df = df_value, 
                                                                       OR = OR,
                                                                       lower_CI = lower_CI,
                                                                       upper_CI = upper_CI,
                                                                       "p-value" = p_value))
}

model1_quart_not_summed <- model1_quart_not_summed |> 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quart_not_summed") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

model2_quart_not_summed <- data.frame(variable = character(),
                                      df = integer(),
                                      OR = numeric(),
                                      lower_CI = numeric(),
                                      upper_CI = numeric(),
                                      p_value = numeric(),
                                      stringsAsFactors = FALSE)

for (var in POPs_included_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quart_not_summed <- rbind(model2_quart_not_summed, data.frame(variable = var,
                                                                       df = df_value, 
                                                                       OR = OR,
                                                                       lower_CI = lower_CI,
                                                                       upper_CI = upper_CI,
                                                                       "p-value" = p_value))
}

model2_quart_not_summed <- model2_quart_not_summed |> 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quart_not_summed") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

sensitivity_results_not_summed_quart <- 
  bind_rows(model1_quart_not_summed, model2_quart_not_summed) |>
  mutate(variable = gsub("_quart", "", variable), 
         variable = gsub('_', '-', variable),
         variable = gsub('OCP-', '', variable),
         variable = fct_recode(variable, "p,p'-DDE" = 'pp-DDE',  "p,p'-DDT" ="pp-DDT"),
         OR = format(OR, nsmall = 1, digits = 1),
         lower_CI = format(lower_CI, nsmall = 1, digits = 1),
         upper_CI =  format(upper_CI, nsmall = 1, digits = 1),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", format(p.value, nsmall = 1, digits = 1)), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  arrange(variable) |>
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value"), 
         lower_CI, upper_CI) 

### gam ----
model1_gam_not_summed <- list()

for (var in POPs_included) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gam_not_summed[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

model2_gam_not_summed <- list()

for (var in POPs_included) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + 
                              smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gam_not_summed[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

# figures ----
## quartiles ----
plot_quart <- main_results |> 
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) |>
  mutate(p.value_shape = ifelse(p.value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted_quart",
                            "Base model" = "base_quart",
                            "Copollutant model" = "copollutant_quart"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         df = fct_relevel(df, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         variable = fct_recode(variable, 
                               "Most\nprevalent\nPCBs" = "PCB_4",
                               "Dioxin-like\nPCBs" = "PCB_DL",
                               "Non-dioxin-\nlike PCBs" = "PCB_NDL",
                               "β-HCH" = "OCP_β_HCH", 
                               "HCB" = "OCP_HCB"), 
        variable = fct_relevel(variable, 
                               "Most\nprevalent\nPCBs", "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE")) |>
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(variable), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

## quartiles not summed ----
plot_quart_sensi_not_summed <- sensitivity_results_not_summed_quart |> 
  mutate(df = fct_recode(df, "Quartile 2" = "2", "Quartile 3" = "3", "Quartile 4" = "4"), 
         df = fct_relevel(df, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         p.value_shape = ifelse(p.value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted_quart_not_summed",
                            "Base model" = "base_quart_not_summed"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         variable = fct_relevel(variable, 
                                "PCB-118", "PCB-156", "PCB-28", "PCB-52", "PCB-74", "PCB-99",
                                "PCB-101", "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183",
                                "PCB-187", "OCP-HCB", "p,p'-DDE", "p,p'-DDT", "OCP-β-HCH", "OCP-transnonachlor",
                                "OCP-oxychlordane", "PBDE-47", "PBDE-99", "PBDE-153"), 
         OR = as.numeric(as.character(OR)), 
         lower_CI = as.numeric(as.character(lower_CI)), 
         upper_CI = as.numeric(as.character(upper_CI))) |>
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(variable), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

## model 1 ----
### gam ----
pollutant_labels <- set_names(
  c("Most prevalent PCBs", "Dioxin-like PCBs","Non-dioxin-like PCBs", "HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

plot_base_gam <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Base model")
  
  p2 <- ggplot(bdd_pred) +
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
rm(pollutant_labels)

### gam not summed ----
POPs_included_labels <- gsub("_", "-", POPs_included)
POPs_included_labels <- gsub("BDE", "PBDE", POPs_included_labels)
POPs_included_labels <- gsub("pp", "p,p'", POPs_included_labels)
pollutant_labels <- set_names(c(POPs_included_labels, POPs_included))

plot_base_gam_not_summed <- map(POPs_included, function(var) {
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula,
               family = binomial,
               method = "REML",
               data = bdd_danish)
  
  bdd_pred <- bdd_danish |>                                                    # création bdd avec covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with( ~ gsub("adj_", "", .x))
  
  pred <- predict(model,
                  newdata = bdd_pred,
                  type = "link",
                  se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      prob = plogis(pred$fit),                                       # plogit does exp(pred$fit) / (1 + exp(pred$fit))
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"],
                nsmall = 1,
                digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =
                                                            2, digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "blue",
                alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate(
      "text",
      x = x_min,
      y = Inf,
      label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
      hjust = 1,
      vjust = 1.2,
      size = 4,
      color = "black"
    ) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(
      axis.text.x = element_text(color = 'white'),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ggtitle("Base model")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) +
    xlab("") +
    theme_minimal()
  
  p <- p1 / p2 + plot_layout(
    ncol = 1,
    nrow = 2,
    heights = c(10, 1),
    guides = 'collect'
  ) +
    theme_minimal()
  p
}) |>
  set_names(POPs_included_labels)
rm(pollutant_labels)

## model 2 ----
### gam ----
pollutant_labels <- set_names(
  c("Most prevalent PCBs", "Dioxin-like PCBs","Non-dioxin-like PCBs", "HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

plot_adjusted_gam <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish |>                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
  
  p2 <- ggplot(bdd_pred) +
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
rm(pollutant_labels)

### gam not summed ----
POPs_included_labels <- gsub("_", "-", POPs_included)
POPs_included_labels <- gsub("BDE", "PBDE", POPs_included_labels)
POPs_included_labels <- gsub("pp", "p,p'", POPs_included_labels)
pollutant_labels <- set_names(c(POPs_included_labels, POPs_included))

plot_adjusted_gam_not_summed <- map(POPs_included, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + 
                              smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish |>                                                    # création bdd avec covariables ramenées à leur moyenne
    mutate(across(all_of(c("sex", "baseline_age", 'smoking_2cat_i', 'bmi', 'fS_Kol', 'marital_status_2cat_i', 'education_i')), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
  
  p2 <- ggplot(bdd_pred) +
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
  set_names(POPs_included_labels)
rm(pollutant_labels)


## model 3 ----
### gam ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_bis)

model <- gam(als ~ s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + 
               s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = "REML", 
             data = bdd_danish)


plot_copollutant_gam <- map(POPs_group_bis, function(var) {
  
  bdd_pred <- bdd_danish|>
    mutate(across(all_of(covariates),                                           # fixe toutes les covariables à leurs moyennes
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.)))))|>
    mutate(across(setdiff(POPs_group_bis, var), 
                  ~ mean(., na.rm = TRUE)))|>                                   # Fixe tous les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates), all_of(setdiff(POPs_group_bis, var)))  
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred|>
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Co-pollutant model")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
})|> set_names(POPs_group_bis)
rm(POPs_group_bis, pollutant_labels_bis, model)


### table POPs - ALS occurence (qgcomp analysis) ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")                                  # remove the 4 most abundant PCB because they are already NDL-PCB
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs",
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"),
  POPs_group_bis)
qgcomp_boot_danish
p <- summary(qgcomp_boot_danish)
POPs_ALS_qgcomp_table_danish <-                                                 # overall results
  tibble(
    study = "Danish",
    model = "copollutant",
    OR = exp(p$coefficients["psi1", "Estimate"]),
    lower_CI = exp(p$coefficients["psi1", "Lower CI"]),
    upper_CI = exp(p$coefficients["psi1", "Upper CI"]),
    p_value = p$coefficients["psi1", "Pr(>|z|)"] ) |>
  mutate(
    OR = sprintf("%.1f", OR),
    lower_CI = sprintf("%.2f", lower_CI),
    upper_CI = sprintf("%.2f", upper_CI),
    `95% CI` = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value` = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")),
    `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`)) |>
  select(study, model, OR, `95% CI`, `p-value`)

### figure POPs - ALS occurence (qgcomp analysis) ----
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

## metanalysis ----
plot_metanalysis_quart <- metanalysis_quart |> 
  mutate(`p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted",
                            "Base model" = "base", 
                            "Copollutant model" = "copollutant"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         term = fct_relevel(term, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         explanatory = fct_recode(explanatory, 
                                  "Most\nprevalent\nPCBs" = "PCB_4",
                                  "Dioxin-like\nPCBs" = "PCB_DL",
                                  "Non-dioxin-\nlike PCBs" = "PCB_NDL",
                                  "β-HCH" = "OCP_β_HCH", 
                                  "HCB" = "OCP_HCB"), 
         explanatory = fct_relevel(explanatory, 
                                   "Most\nprevalent\nPCBs", "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane"), 
         OR = as.numeric(OR), 
         lower_CI = as.numeric(lower_CI), 
         upper_CI = as.numeric(upper_CI)) |>
  ggplot(aes(x = term, y = OR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()


# Assemblage -----
results_POPs_ALS_occurrence <- 
  list(main = list(main_results = main_results, 
                   covar = covar, 
                   results_quart = results_quart, 
                   model1_gam = model1_gam, 
                   model2_gam = model2_gam, 
                   plot_quart = plot_quart,                                     # co-pollutant model with only the POP of interest as quartile and the other as s() in a GAM model
                   plot_base_gam = plot_base_gam, 
                   plot_adjusted_gam = plot_adjusted_gam, 
                   plot_copollutant_gam = plot_copollutant_gam, 
                   POPs_ALS_qgcomp_table_danish = POPs_ALS_qgcomp_table_danish, 
                   POPs_ALS_qgcomp_figure_danish = POPs_ALS_qgcomp_figure_danish), 
       sensitivity_not_summed = list(sensitivity_results_not_summed_quart = sensitivity_results_not_summed_quart, 
                                     model1_gam_not_summed = model1_gam_not_summed, 
                                     model2_gam_not_summed = model2_gam_not_summed,
                                     model1_quart_not_summed = model1_quart_not_summed, 
                                     model2_quart_not_summed = model2_quart_not_summed, 
                                     plot_quart_sensi_not_summed = plot_quart_sensi_not_summed, 
                                     plot_base_gam_not_summed = plot_base_gam_not_summed, 
                                     plot_adjusted_gam_not_summed = plot_adjusted_gam_not_summed), 
       metanalysis = list(metanalysis_quart = metanalysis_quart, 
                          plot_metanalysis_quart = plot_metanalysis_quart))

rm(main_results, covar, results_quart, model1_gam, model2_gam, 
   sensitivity_results_not_summed_quart, model1_gam_not_summed, model2_gam_not_summed, model1_quart_not_summed, model2_quart_not_summed, 
   plot_quart, plot_quart_sensi_not_summed, 
   plot_base_gam, plot_adjusted_gam, 
   plot_base_gam_not_summed, plot_adjusted_gam_not_summed, 
   plot_copollutant_gam, 
   metanalysis_quart, plot_metanalysis_quart, 
   covariates_danish, covariates_finnish, 
   POPs_ALS_qgcomp_table_danish, 
   POPs_ALS_qgcomp_figure_danish)


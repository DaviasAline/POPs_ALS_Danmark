# Aline Davias
# 03/02/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.1_analyses_descriptive.R")

covariates <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'cholesterol_i', 'marital_status_2cat_i', 'education_i')


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
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tbl_regression(exponentiate = TRUE, 
                 estimate_fun = label_number(accuracy = .1, decimal.mark = "."),
                 pvalue_fun = custom_pvalue_fun) |>
  bold_labels()), 
  tab_spanner = c("**Univariate**", "**Adjusted**"))

# main analysis ----
### model 1 ----
#### spline transformation ----
model1_spline <- data.frame(variable = character(),
                            df = integer(),
                            OR = numeric(),
                            lower_CI = numeric(),
                            upper_CI = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_spline <- rbind(model1_spline, data.frame(variable = var,
                                                   df = df_value, 
                                                   OR = OR,
                                                   lower_CI = lower_CI,
                                                   upper_CI = upper_CI,
                                                   "p-value" = p_value))
}

model1_spline <- model1_spline %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quartile transformation ----
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

model1_quart <- model1_quart %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quart") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model1_quadratic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quadratic <- rbind(model1_quadratic, data.frame(variable = var,
                                                         df = df_value, 
                                                         OR = OR,
                                                         lower_CI = lower_CI,
                                                         upper_CI = upper_CI,
                                                         "p-value" = p_value))
}

model1_quadratic <- model1_quadratic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


#### cubic transformation ----
model1_cubic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_cubic <- rbind(model1_cubic, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model1_cubic <- model1_cubic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### gamm ----
model1_gamm <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gamm[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### model 2 ----
# matched on sex and age, adjusted on for smoking_2cat_i, BMI, serum total cholesterol_i, marital status and education

#### spline transformation ----
model2_spline <- data.frame(variable = character(),
                            df = integer(),
                            OR = numeric(),
                            lower_CI = numeric(),
                            upper_CI = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_spline <- rbind(model2_spline, data.frame(variable = var,
                                                   df = df_value, 
                                                   OR = OR,
                                                   lower_CI = lower_CI,
                                                   upper_CI = upper_CI,
                                                   "p-value" = p_value))
}

model2_spline <- model2_spline %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quartile transformation ----
model2_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl(paste0("^", var), term))
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

model2_quart <- model2_quart %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quart") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model2_quadratic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quadratic <- rbind(model2_quadratic, 
                            data.frame(variable = var,
                                       df = df_value, 
                                       OR = OR,
                                       lower_CI = lower_CI,
                                       upper_CI = upper_CI,
                                       "p-value" = p_value))
}

model2_quadratic <- model2_quadratic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### cubic transformation ----
model2_cubic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_cubic <- rbind(model2_cubic, 
                        data.frame(variable = var,
                                   df = df_value, 
                                   OR = OR,
                                   lower_CI = lower_CI,
                                   upper_CI = upper_CI,
                                   "p-value" = p_value))
}

model2_cubic <- model2_cubic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### gamm ----
model2_gamm <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### models 3 ----
#### spline transformation ----
model3_spline_HCB <- 
  clogit(als ~ 
           ns(OCP_HCB, df = 4) + 
           strata(match) + 
           PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_HCB", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_PCB_DL <- 
  clogit(als ~ 
           ns(PCB_DL, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_PCB_NDL <- 
  clogit(als ~ 
           ns(PCB_NDL, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_DL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_ΣPBDE <- 
  clogit(als ~ 
           ns(ΣPBDE, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_ΣDDT <- 
  clogit(als ~ 
           ns(ΣDDT, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + OCP_β_HCH + Σchlordane +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_β_HCH <- 
  clogit(als ~ 
           ns(OCP_β_HCH, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_β_HCH", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_Σchlordane <- 
  clogit(als ~ 
           ns(Σchlordane, df = 4) + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline <- bind_rows(
  model3_spline_HCB, model3_spline_PCB_DL, model3_spline_PCB_NDL, model3_spline_ΣPBDE, model3_spline_ΣDDT, model3_spline_β_HCH, model3_spline_Σchlordane) |>
  mutate(
    variable = gsub('ns\\(', '', term), 
    variable = gsub(', df = 4)1', '', variable),
    variable = gsub(', df = 4)2', '', variable),
    variable = gsub(', df = 4)3', '', variable),
    variable = gsub(', df = 4)4', '', variable),
    model = "copollutant_spline", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)

rm(model3_spline_HCB, model3_spline_PCB_DL, model3_spline_PCB_NDL, model3_spline_ΣPBDE, model3_spline_ΣDDT, model3_spline_β_HCH, model3_spline_Σchlordane)

#### quartile tranformation ----
model3_quart_PCB_DL <- 
  clogit(als ~ 
           PCB_DL_quart + 
           strata(match) + 
           OCP_HCB + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_PCB_NDL <- 
  clogit(als ~ 
           PCB_NDL_quart + 
           strata(match) + 
           OCP_HCB + PCB_DL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_HCB <- 
  clogit(als ~ 
           OCP_HCB_quart + 
           strata(match) + 
           PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_HCB_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_ΣDDT <- 
  clogit(als ~ 
           ΣDDT_quart + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_β_HCH <- 
  clogit(als ~ 
           OCP_β_HCH_quart + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_β_HCH_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_Σchlordane <- 
  clogit(als ~ 
           Σchlordane_quart + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_ΣPBDE <- 
  clogit(als ~ 
           ΣPBDE_quart + 
           strata(match) + 
           OCP_HCB + PCB_DL + PCB_NDL + ΣDDT + OCP_β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart <- bind_rows(
  model3_quart_HCB, model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_ΣPBDE, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane) |>
  mutate(
    variable = gsub("_quartQ2", "", term), 
    variable = gsub("_quartQ3", "", variable), 
    variable = gsub("_quartQ4", "", variable), 
    model = "copollutant_quart", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)

model3_quart_bis <- 
  clogit(als ~ 
           PCB_DL_quart + 
           strata(match) + 
           OCP_HCB_quart + PCB_NDL_quart + ΣPBDE_quart + ΣDDT_quart + OCP_β_HCH_quart + Σchlordane_quart + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("quart", term)) |>
  mutate(df = str_extract(term, "Q[2-4]"), 
         variable = str_remove(term, "Q[2-4]"), 
         variable = str_remove(variable, "_quart"), 
         model = "copollutant_quart_bis", 
         OR = exp(estimate), 
         lower_CI = exp(estimate - 1.96 * std.error), 
         upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)

rm(model3_quart_HCB, model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_ΣPBDE, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane)

#### quadratic tranformation ----
model3_quadra_HCB <- 
  glm(als ~ 
        poly(OCP_HCB, 2) + 
        baseline_age + sex + 
        PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_HCB", term)) |>
  mutate(df = c(1, 2))

model3_quadra_PCB_DL <- 
  glm(als ~ 
        poly(PCB_DL, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term)) |>
  mutate(df = c(1, 2))

model3_quadra_PCB_NDL <- 
  glm(als ~ 
        poly(PCB_NDL, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term)) |>
  mutate(df = c(1, 2))

model3_quadra_ΣPBDE <- 
  glm(als ~ 
        poly(ΣPBDE, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2))

model3_quadra_ΣDDT <- 
  glm(als ~ 
        poly(ΣDDT, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term)) |>
  mutate(df = c(1, 2))

model3_quadra_β_HCH <- 
  glm(als ~ 
        poly(OCP_β_HCH, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_β_HCH", term)) |>
  mutate(df = c(1, 2))

model3_quadra_Σchlordane <- 
  glm(als ~ 
        poly(Σchlordane, 2) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term)) |>
  mutate(df = c(1, 2))

model3_quadratic <- bind_rows(
  model3_quadra_HCB, model3_quadra_PCB_DL, model3_quadra_PCB_NDL, model3_quadra_ΣPBDE, model3_quadra_ΣDDT, model3_quadra_β_HCH, model3_quadra_Σchlordane) |>
  mutate(
    variable = gsub("poly\\(", "", term),
    variable = gsub(", 2)1", "", variable),
    variable = gsub(", 2)2", "", variable),
    model = "copollutant_quadra", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)

rm(model3_quadra_HCB, model3_quadra_PCB_DL, model3_quadra_PCB_NDL, model3_quadra_ΣPBDE, model3_quadra_ΣDDT, model3_quadra_β_HCH, model3_quadra_Σchlordane)

#### cubic tranformation ----
model3_cubic_HCB <- 
  glm(als ~ 
        poly(OCP_HCB, 3) + 
        baseline_age + sex + 
        PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_HCB", term)) |>
  mutate(df = c(1, 2, 3))

model3_cubic_PCB_DL <- 
  glm(als ~ 
        poly(PCB_DL, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_PCB_NDL <- 
  glm(als ~ 
        poly(PCB_NDL, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + ΣPBDE + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_ΣPBDE <- 
  glm(als ~ 
        poly(ΣPBDE, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣDDT + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2, 3))

model3_cubic_ΣDDT <- 
  glm(als ~ 
        poly(ΣDDT, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + OCP_β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_β_HCH <- 
  glm(als ~ 
        poly(OCP_β_HCH, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("OCP_β_HCH", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_Σchlordane <- 
  glm(als ~ 
        poly(Σchlordane, 3) + 
        baseline_age + sex + 
        OCP_HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + OCP_β_HCH + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic <- bind_rows(
  model3_cubic_HCB, model3_cubic_PCB_DL, model3_cubic_PCB_NDL, model3_cubic_ΣPBDE, model3_cubic_ΣDDT, model3_cubic_β_HCH, model3_cubic_Σchlordane) |>
  mutate(
    variable = gsub("poly\\(", "", term),
    variable = gsub(", 3)1", "", variable),
    variable = gsub(", 3)2", "", variable),
    variable = gsub(", 3)3", "", variable),
    model = "copollutant_cubic", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)

rm(model3_cubic_HCB, model3_cubic_PCB_DL, model3_cubic_PCB_NDL, model3_cubic_ΣPBDE, model3_cubic_ΣDDT, model3_cubic_β_HCH, model3_cubic_Σchlordane)

### heterogeneity tests ----
#### model 1 spline ----
heterogeneity_base_spline <- data.frame(variable = character(),
                                        model = factor(), 
                                        p.value_heterogeneity = numeric(), 
                                        stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_spline <- rbind(heterogeneity_base_spline, 
                                     data.frame(variable = var,
                                                model = "base_spline",
                                                p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

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

#### model 2 spline ----
heterogeneity_adjusted_spline <- data.frame(variable = character(),
                                            model = factor(), 
                                            p.value_heterogeneity = numeric(), 
                                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_spline <- rbind(heterogeneity_adjusted_spline, 
                                         data.frame(variable = var,
                                                    model = "adjusted_spline",
                                                    p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

#### model 2 quartile ----
heterogeneity_adjusted_quart <- data.frame(variable = character(),
                                           model = factor(), 
                                           p.value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(variable = var,
                                                   model = "adjusted_quart",
                                                   p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_spline, heterogeneity_base_quart, heterogeneity_adjusted_spline, heterogeneity_adjusted_quart) %>%
  mutate(variable = gsub("_quart", "", variable))

### trend tests ----
#### model 1 quartile ----
trend_base <- data.frame(variable = character(),
                         model = factor(), 
                         p.value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_trend <- anova$`Pr(>|Chi|)`[2]
  
  trend_base <- rbind(trend_base, 
                      data.frame(variable = var,
                                 model = "base_quart",
                                 p.value_trend = p.value_trend))
}
rm(var, test_1, test_2, formula, anova, p.value_trend)

#### model 2 quartile ----
trend_adjusted <- data.frame(variable = character(),
                             model = factor(), 
                             p.value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_trend <- anova$`Pr(>|Chi|)`[2]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(variable = var,
                                     model = "adjusted_quart",
                                     p.value_trend = p.value_trend))
}
rm(var, test_1, test_2, formula, anova, p.value_trend)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted) %>%
  mutate(variable = gsub("_quart_med", "", variable))

### merging the main results ----
main_results <- bind_rows(model1_spline, 
                          model1_quart, 
                          model1_quadratic,
                          model1_cubic, 
                          model2_spline, 
                          model2_quart,
                          model2_quadratic,
                          model2_cubic,
                          model3_spline, 
                          model3_quart, 
                          model3_quart_bis,
                          model3_quadratic, 
                          model3_cubic) %>% 
  mutate(variable = gsub("_quart", "", variable), 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) %>%
  arrange(variable) %>%
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

results_spline <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI, -p.value_raw) |>
  filter(model %in% c('base_spline', 'adjusted_spline', 'copollutant_spline')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_spline"), 
         contains("adjusted_spline"), 
         contains("copollutant_spline")) 
colnames(results_spline) <- gsub('_spline', '', colnames(results_spline))

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

results_quadratic <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI, -p.value_raw) |>
  filter(model %in% c('base_quadra', 'adjusted_quadra', 'copollutant_quadra')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         degree = df,
         contains("base_quadra"), 
         contains("adjusted_quadra"), 
         contains("copollutant_quadra")) 
colnames(results_quadratic) <- gsub('_quadra', '', colnames(results_quadratic))

results_cubic <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI, -p.value_raw) |>
  filter(model %in% c('base_cubic', 'adjusted_cubic', 'copollutant_cubic')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         degree = df,
         contains("base_cubic"), 
         contains("adjusted_cubic"), 
         contains("copollutant_cubic")) 
colnames(results_cubic) <- gsub('_cubic', '', colnames(results_cubic))

rm(model1_quart, model1_spline, model1_quadratic, model1_cubic,
   model2_quart, model2_spline, model2_quadratic, model2_cubic,
   model3_quart, model3_spline, model3_quadratic, model3_cubic,
   model3_quart_bis, 
   heterogeneity_base_quart, heterogeneity_base_spline,
   heterogeneity_adjusted_spline, heterogeneity_adjusted_quart, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests)

# sensitivity analyses without outliers ----
### model 1  ----
# adjusted for sex and age

#### spline transformation ----
model1_spline_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_spline_outlier <- rbind(model1_spline_outlier, data.frame(variable = var,
                                                                   df = df_value, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   "p-value" = p_value))
}

model1_spline_outlier <- model1_spline_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### quadratic transformation ----
model1_quadratic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quadratic_outlier <- rbind(model1_quadratic_outlier, data.frame(variable = var,
                                                                         df = df_value, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         "p-value" = p_value))
}

model1_quadratic_outlier <- model1_quadratic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model1_cubic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_cubic_outlier <- rbind(model1_cubic_outlier, data.frame(variable = var,
                                                                 df = df_value, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 "p-value" = p_value))
}

model1_cubic_outlier <- model1_cubic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)


#### gamm ----
model1_gamm_outliers <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)
### model 2 ----
# adjusted for sex, age, smoking_2cat_i, BMI, serum total cholesterol_i, marital status and education

#### spline transformation ----
model2_spline_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_spline_outlier <- rbind(model2_spline_outlier, data.frame(variable = var,
                                                                   df = df_value, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   "p-value" = p_value))
}

model2_spline_outlier <- model2_spline_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model2_quadratic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quadratic_outlier <- rbind(model2_quadratic_outlier, 
                                    data.frame(variable = var,
                                               df = df_value, 
                                               OR = OR,
                                               lower_CI = lower_CI,
                                               upper_CI = upper_CI,
                                               "p-value" = p_value))
}

model2_quadratic_outlier <- model2_quadratic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model2_cubic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_cubic_outlier <- rbind(model2_cubic_outlier, 
                                data.frame(variable = var,
                                           df = df_value, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           "p-value" = p_value))
}

model2_cubic_outlier <- model2_cubic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### gamm ----
model2_gamm_outliers <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + 
                              smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

#### merging  ----
sensitivity_results_outlier <- bind_rows(model1_spline_outlier,
                                  model1_quadratic_outlier,
                                 model1_cubic_outlier, 
                                 model2_spline_outlier,
                                model2_quadratic_outlier,
                                 model2_cubic_outlier) %>% 
  arrange(variable) %>%
  mutate(
    OR = as.numeric(sprintf("%.1f", OR)), 
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
    p.value = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")),
    "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) %>%
  select(-starts_with("lower"), -starts_with("upper")) %>%
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value")) 

results_spline_outliers <- 
  sensitivity_results_outlier |>
  filter(model %in% c('base_spline', 'adjusted_spline', 'copollutant_spline')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_spline"), 
         contains("adjusted_spline"), 
         contains("copollutant_spline")) 
colnames(results_spline_outliers) <- gsub('_spline', '', colnames(results_spline_outliers))

results_quadratic_outliers <- 
  sensitivity_results_outlier |>
  filter(model %in% c('base_quadra', 'adjusted_quadra', 'copollutant_quadra')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_quadra"), 
         contains("adjusted_quadra"), 
         contains("copollutant_quadra")) 
colnames(results_quadratic_outliers) <- gsub('_quadra', '', colnames(results_quadratic_outliers))

results_cubic_outliers <- 
  sensitivity_results_outlier |>
  filter(model %in% c('base_cubic', 'adjusted_cubic', 'copollutant_cubic')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_cubic"), 
         contains("adjusted_cubic"), 
         contains("copollutant_cubic")) 
colnames(results_cubic_outliers) <- gsub('_cubic', '', colnames(results_cubic_outliers))

rm(model1_spline_outlier,
   model1_quadratic_outlier,
   model1_cubic_outlier, 
   model2_spline_outlier,
   model2_quadratic_outlier,
   model2_cubic_outlier)


## sensitivity analyses pollutants not summed ----
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

model1_quart_not_summed <- model1_quart_not_summed %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quart_not_summed") %>%
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
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl(paste0("^", var), term))
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

model2_quart_not_summed <- model2_quart_not_summed %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quart_not_summed") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


sensitivity_results_not_summed_quart <- 
  bind_rows(model1_quart_not_summed, model2_quart_not_summed) |>
  mutate(variable = gsub("_quart", "", variable), 
         variable = gsub("BDE", "PBDE", variable),
         variable = gsub('_', '-', variable),
         variable = fct_recode(variable, "p,p'-DDE" = 'pp-DDE',  "p,p'-DDT" ="pp-DDT"),
         OR = format(OR, nsmall = 1, digits = 1),
         lower_CI = format(lower_CI, nsmall = 1, digits = 1),
         upper_CI =  format(upper_CI, nsmall = 1, digits = 1),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", format(p.value, nsmall = 1, digits = 1)), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) %>%
  arrange(variable) %>%
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value"), 
         lower_CI, upper_CI) 

### gamm ----
model1_gamm_not_summed <- list()

for (var in POPs_included) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

model2_gamm_not_summed <- list()

for (var in POPs_included) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + 
                              smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

# figures ----
## quartiles ----
plot_quart <- main_results %>% 
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) %>%
  mutate(df = fct_recode(df, "Quartile 2" = "2", "Quartile 3" = "3", "Quartile 4" = "4"), 
         p.value_shape = ifelse(p.value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
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
                               "β-HCH" = "OCP_β_HCH"), 
        variable = fct_relevel(variable, 
                               "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "Most\nprevalent\nPCBs", "OCP_HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE")) %>%
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


plot_quart_bis <- main_results %>% 
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart_bis')) %>%
  mutate(df = fct_recode(df, "Quartile 2" = "2", "Quartile 3" = "3", "Quartile 4" = "4", "Quartile 2" = "Q2", "Quartile 3" = "Q3", "Quartile 4" = "Q4"), 
         p.value_shape = ifelse(p.value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted_quart",
                            "Base model" = "base_quart",
                            "Copollutant model" = "copollutant_quart_bis"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         df = fct_relevel(df, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         variable = fct_recode(variable, 
                               "Most\nprevalent\nPCBs" = "PCB_4",
                               "Dioxin-like\nPCBs" = "PCB_DL",
                               "Non-dioxin-\nlike PCBs" = "PCB_NDL",
                               "β-HCH" = "OCP_β_HCH"), 
         variable = fct_relevel(variable, 
                                "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "Most\nprevalent\nPCBs", "OCP_HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE")) %>%
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(variable), cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

## quartiles not summed ----
plot_quart_sensi_not_summed <- sensitivity_results_not_summed_quart %>% 
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
                                "PCB-187", "OCP_HCB", "p,p'-DDE", "p,p'-DDT", "β-HCH", "Transnonachlor",
                                "Oxychlordane", "PBDE-47", "PBDE-99", "PBDE-153"), 
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
### spline ----
plot_base_spline <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  pred <- predict(model, newdata = new_data, type = "lp", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) + 
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = "spline") +
    theme_minimal()
  
  plot_base_spline[[var]] <- plot
}
rm(var, formula, model, new_data, pred, plot)

### spline outliers ----
plot_base_spline_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(c('sex', 'baseline_age'), collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in c('sex', 'baseline_age')) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_base_spline_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)


### quadratic ----
plot_base_quadratic <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in c("sex", "baseline_age")) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic') +
    theme_minimal()
  
  plot_base_quadratic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### quadratic outliers ----
plot_base_quadratic_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + sex + baseline_age"))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in c("sex", "baseline_age")) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_base_quadratic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic ----
plot_base_cubic <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in c("sex", "baseline_age")) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic') +
    theme_minimal()
  
  plot_base_cubic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### cubic outliers ----
plot_base_cubic_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + sex + baseline_age"))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in c("sex", "baseline_age")) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_base_cubic_outlier[[var]] <- plot
}
rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","OCP_HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

plot_base_gamm <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(c("sex", "baseline_age")), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels)

### gamm outliers ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","OCP_HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group_outlier)

plot_base_gamm_outlier <- map(POPs_group_outlier, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(c("sex", "baseline_age")), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels)

### gamm not summed ----
POPs_included_labels <- gsub("_", "-", POPs_included)
POPs_included_labels <- gsub("BDE", "PBDE", POPs_included_labels)
POPs_included_labels <- gsub("pp", "p,p'", POPs_included_labels)
pollutant_labels <- set_names(c(POPs_included_labels, POPs_included))

plot_base_gamm_not_summed <- map(POPs_included, function(var) {
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula,
               family = binomial,
               method = "REML",
               data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec covariables ramenées à leur moyenne
    mutate(across(
      all_of(c("sex", "baseline_age")),
      ~ if (is.numeric(.))
        mean(., na.rm = TRUE)
      else
        names(which.max(table(.))),
      .names = "adj_{.col}"
    )) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with( ~ gsub("adj_", "", .x))
  
  pred <- predict(model,
                  newdata = bdd_pred,
                  type = "link",
                  se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(
      prob = plogis(pred$fit),
      # plogit does exp(pred$fit) / (1 + exp(pred$fit))
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit)
    )
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"],
                nsmall = 1,
                digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =
                                                            2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "blue",
                alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate(
      "text",
      x = x_max,
      y = Inf,
      label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
      hjust = 1,
      vjust = 1.2,
      size = 4,
      color = "black"
    ) +
    theme_minimal() +
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
}) %>%
  set_names(POPs_included_labels)
rm(pollutant_labels)

## model 2 ----
### spline ----
plot_adjusted_spline <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + 
                               strata(match) + 
                               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  
  model <- clogit(formula, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "lp", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = "spline") +
    theme_minimal()
  
  plot_adjusted_spline[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### spline outliers ----
plot_adjusted_spline_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_adjusted_spline_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)


### quadratic  ----
plot_adjusted_quadratic <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic') +
    theme_minimal()
  
  plot_adjusted_quadratic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### quadratic outliers ----
plot_adjusted_quadratic_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_adjusted_quadratic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic ----
plot_adjusted_cubic <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic') +
    theme_minimal()
  
  plot_adjusted_cubic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### cubic outliers -----
plot_adjusted_cubic_outlier <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_adjusted_cubic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)


### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","OCP_HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

plot_adjusted_gamm <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels)

### gamm outliers ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","OCP_HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group_outlier)

plot_adjusted_gamm_outlier <- map(POPs_group_outlier, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels)

### gamm not summed ----
POPs_included_labels <- gsub("_", "-", POPs_included)
POPs_included_labels <- gsub("BDE", "PBDE", POPs_included_labels)
POPs_included_labels <- gsub("pp", "p,p'", POPs_included_labels)
pollutant_labels <- set_names(c(POPs_included_labels, POPs_included))

plot_adjusted_gamm_not_summed <- map(POPs_included, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + 
                              smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec covariables ramenées à leur moyenne
    mutate(across(all_of(c("sex", "baseline_age", 'smoking_2cat_i', 'bmi', 'cholesterol_i', 'marital_status_2cat_i', 'education_i')), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
}) %>% 
  set_names(POPs_included_labels)
rm(pollutant_labels)


## model 3 ----
### gamm ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "OCP_HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_bis)

model <- gam(als ~ s(PCB_DL) + s(PCB_NDL) + s(OCP_HCB) + s(ΣDDT) + 
               s(OCP_β_HCH) + s(Σchlordane) + s(ΣPBDE) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = "REML", 
             data = bdd_danish)


plot_copollutant_gamm <- map(POPs_group_bis, function(var) {
  
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
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
rm(POPs_group_bis, pollutant_labels_bis)


### gamm outliers ----
POPs_group_outlier_bis <- setdiff(POPs_group_outlier, "PCB_4_outlier")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "OCP_HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_outlier_bis)

model <- gam(als ~ s(PCB_DL_outlier) + s(PCB_NDL_outlier) + s(OCP_HCB_outlier) + s(ΣDDT_outlier) + 
               s(OCP_β_HCH_outlier) + s(Σchlordane_outlier) + s(ΣPBDE_outlier) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = "REML", 
             data = bdd_danish)

plot_copollutant_gamm_outlier <- map(POPs_group_outlier_bis, function(var) {
  
  bdd_pred <- bdd_danish|>
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.)))))|>
    mutate(across(setdiff(POPs_group_outlier_bis, var), 
                  ~ mean(., na.rm = TRUE)))|>                                 # Fixe les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates), all_of(setdiff(POPs_group_outlier_bis, var)))  
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred|>
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  
  # edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  # p_value <- model_summary$s.table[1, "p-value"]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
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
})|> set_names(POPs_group_outlier_bis)
rm(POPs_group_outlier_bis, pollutant_labels_bis)



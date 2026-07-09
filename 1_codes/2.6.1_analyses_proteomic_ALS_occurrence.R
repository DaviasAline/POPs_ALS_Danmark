# Aline Davias
# October 20, 2025 
# Analysis of als risk depending on proteomic profile

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# Data preparation ----
## Sensi_1 ----
proteomic_sensi_1 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1")
proteomic_sd_sensi_1 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_sd")
proteomic_quart_sensi_1 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_quart")
proteomic_quart_med_sensi_1 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_quart_med")

bdd_danish_sensi_1 <- 
  bdd_danish  |> 
  mutate(
    proteomic_neuro_explo_NEFL_sensi_1 = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic_sensi_1),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(
    all_of(proteomic_sensi_1),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med")) |>
  mutate(across(all_of(proteomic_sensi_1), ~ factor(ntile(.x, 4),                           
                                                    labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) 

## Sensi_2 ----
bdd_danish_sensi_2 <- 
  bdd_danish |>
  group_by(match) |>                                                             
  filter(any(als == 1 & follow_up < 60)) |>                                     # sensi 2 : we remove cases (and their controls) with follow-up > 60 months 
  ungroup()|>
  mutate(across(all_of(proteomic),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic), ~ factor(ntile(.x, 4),                           
                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"
  ))

## Sensi_3 ----
bdd_danish_sensi_3 <- 
  bdd_danish |>
  group_by(match) |>
  filter(any(als == 1 & follow_up > 60)) |>                                     # had to remove 17 cases with follow-up<60 months (and their controls) -> 51 people
  ungroup() |>
  mutate(across(all_of(proteomic),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic), ~ factor(ntile(.x, 4),                           
                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"
  ))

## Sensi_1_3 ----

proteomic_sensi_1_3 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_3")
proteomic_sd_sensi_1_3 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_3_sd")
proteomic_quart_sensi_1_3 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_3_quart")
proteomic_quart_med_sensi_1_3 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_3_quart_med")


bdd_danish_sensi_1_3 <- 
  bdd_danish |>
  group_by(match) |>                                                            # sensi 3 : we remove cases (and their controls) with follow-up <60 months to see if the association with NfL really happens a long time before ALS diagnosis 
  filter(any(als == 1 & follow_up > 60)) |>
  ungroup()|>
  mutate(                                                                       # sensi 1 : we remove NfL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sensi_1_3 =
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic_sensi_1_3),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic_sensi_1_3), ~ factor(ntile(.x, 4),                           
                                                      labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic_sensi_1_3),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"
  ))

## Sensi 1_3_4 ----
proteomic_sensi_1_3_4 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_3_4")
proteomic_sd_sensi_1_3_4 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_3_4_sd")
proteomic_quart_sensi_1_3_4 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_3_4_quart")
proteomic_quart_med_sensi_1_3_4 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_3_4_quart_med")

bdd_danish_sensi_1_3_4 <- 
  bdd_danish |>
  group_by(match) |>                                                              
  filter(any(als == 1 & follow_up > 60)) |>                                     # sensi 3 : we remove cases (and their controls) with follow-up < 60 months
  ungroup() |>
  mutate(                                                                       # sensi 1 : we remove NfL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sensi_1_3_4 = 
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 4 : we remove 50% of the highest values of follow-up (after removing the follow-up<5years)
  group_by(match) |>
  filter(any(als == 1 & follow_up <= seuil)) |>
  ungroup() |>
  select(-seuil) |>
  mutate(across(all_of(proteomic_sensi_1_3_4),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic_sensi_1_3_4), ~ factor(ntile(.x, 4),                           
                                                        labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic_sensi_1_3_4),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"
  ))

## Sensi_1_3_5 ----
proteomic_sensi_1_3_5 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_3_5")
proteomic_sd_sensi_1_3_5 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_3_5_sd")
proteomic_quart_sensi_1_3_5 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_3_5_quart")
proteomic_quart_med_sensi_1_3_5 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_3_5_quart_med")


bdd_danish_sensi_1_3_5 <- 
  bdd_danish |>
  group_by(match) |>                                                            # sensi 3 : we remove cases (and their controls) with follow-up < 60 months  
  filter(any(als == 1 & follow_up > 60)) |>
  ungroup()|>
  mutate(                                                                       # sensi 1 : we remove NfL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sensi_1_3_5 = 
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 5 : we remove 50% of the lowest values of follow-up (after removing the follow-up<5years)
  group_by(match) |>
  filter(any(als == 1 & follow_up > seuil)) |>
  ungroup() |>
  select(-seuil) |>
  mutate(across(all_of(proteomic_sensi_1_3_5),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic_sensi_1_3_5), ~ factor(ntile(.x, 4),                           
                                                        labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic_sensi_1_3_5),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"
  ))

## Sensi_6 ----
bdd_danish_sensi_6 <- bdd_danish |>                                             # just for visualisation of follow up distribution 
  group_by(match) |>
  mutate(follow_up_years = follow_up/12, 
         follow_up_ter = follow_up_years[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE))

## Sensi_1_6 ----

bdd_danish_sensi_1_6 <- 
  bdd_danish  |>
  mutate(
    proteomic_neuro_explo_NEFL_sd_sensi_1_6 = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL_sd)) |>
  group_by(match) |>
  mutate(follow_up_ter = follow_up[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE, 
      labels = c("Tertile 1", "Tertile 2", "Tertile 3")))
bdd_danish_sensi_1_6_T1 <- 
  bdd_danish_sensi_1_6 |>
  filter(follow_up_ter == "Tertile 1")
bdd_danish_sensi_1_6_T2 <- 
  bdd_danish_sensi_1_6 |>
  filter(follow_up_ter == "Tertile 2")
bdd_danish_sensi_1_6_T3 <- 
  bdd_danish_sensi_1_6 |>
  filter(follow_up_ter == "Tertile 3")

proteomic_sd_sensi_1_6 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sd_sensi_1_6")

## Sensi_1_7_female ----
proteomic_sensi_1_7_female <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_7_female")
proteomic_sd_sensi_1_7_female <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_7_female_sd")
proteomic_quart_sensi_1_7_female <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_7_female_quart")
proteomic_quart_med_sensi_1_7_female <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_7_female_quart_med")


bdd_danish_sensi_1_7_female <- 
  bdd_danish |>
  group_by(match) |>                                                            # sensi 7 : we keep cases (and their controls) that are females 
  filter(any(als == 1 & sex == "Female")) |>
  ungroup()|>
  mutate(                                                                       # sensi 1 : we remove NfL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sensi_1_7_female =
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic_sensi_1_7_female),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic_sensi_1_7_female), ~ factor(ntile(.x, 4),                           
                                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic_sensi_1_7_female),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med"))

## Sensi_1_7_male ----
proteomic_sensi_1_7_male <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_7_male")
proteomic_sd_sensi_1_7_male <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sensi_1_7_male_sd")
proteomic_quart_sensi_1_7_male <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_sensi_1_7_male_quart")
proteomic_quart_med_sensi_1_7_male <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_sensi_1_7_male_quart_med")


bdd_danish_sensi_1_7_male <- 
  bdd_danish |>
  group_by(match) |>                                                            # sensi 7 : we keep cases (and their controls) that are male 
  filter(any(als == 1 & sex == "Male")) |>
  ungroup()|>
  mutate(                                                                       # sensi 1 : we remove NfL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sensi_1_7_male =
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic_sensi_1_7_male),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(proteomic_sensi_1_7_male), ~ factor(ntile(.x, 4),                           
                                                           labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(
    all_of(proteomic_sensi_1_7_male),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med")) 




# Effects of the covariates on ALS ----
covar <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
      select(als, 'sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'fS_Kol', 'marital_status_2cat_i', 'education_i') |>
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


# Main analysis ----
## model 1 sd ----
model1_sd <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd <- rbind(model1_sd, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model1_sd <- model1_sd |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----
model1_quart <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart <- rbind(model1_quart, data.frame(explanatory = var,
                                                 term = term, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 p_value = p_value))
}

model1_quart <- model1_quart |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "main") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base <- data.frame(explanatory = character(),
                                       model = factor(),
                                       p_value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base <- rbind(heterogeneity_base, 
                                    data.frame(explanatory = var,
                                               model = "base", 
                                               analysis = "main",
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base <- data.frame(explanatory = character(),
                         model = factor(), 
                         p_value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base <- rbind(trend_base, 
                      data.frame(explanatory = var,
                                 model = "base", 
                                 analysis = "main",
                                 p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)

## model 1 gams ----
model1_gam <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

## model 2 sd ----
model2_sd <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_sd <- rbind(model2_sd, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model2_sd <- model2_sd |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

## model 2 quartiles ----
model2_quart <- data.frame(explanatory = character(),
                           term = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart <- rbind(model2_quart, data.frame(explanatory = var,
                                                 term = term, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 p_value = p_value))
}

model2_quart <- model2_quart |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "main") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted <- data.frame(explanatory = character(),
                                           model = factor(), 
                                           p_value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted <- rbind(heterogeneity_adjusted, 
                                        data.frame(explanatory = var,
                                                   model = "adjusted", 
                                                   analysis = "main",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted <- data.frame(explanatory = character(),
                             model = factor(), 
                             p_value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(explanatory = var,
                                     model = "adjusted", 
                                     analysis = "main",
                                     p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)


## model 2 gams ----
model2_gam <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


# Sensi 1 - Removing NfL outlier ----
### model 1 sd ----
model1_sd_sensi_1 <- data.frame(explanatory = character(),
                        term = integer(),
                        OR = numeric(),
                        lower_CI = numeric(),
                        upper_CI = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1 <- rbind(model1_sd_sensi_1, data.frame(explanatory = var,
                                           term = term, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           p_value = p_value))
}

model1_sd_sensi_1 <- model1_sd_sensi_1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1") |>
  select(explanatory, model, everything())


rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 1 quartiles ----
model1_quart_sensi_1 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1 <- rbind(model1_quart_sensi_1, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_quart_sensi_1 <- model1_quart_sensi_1 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1 <- map_dfr(
  proteomic_quart_sensi_1,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1",
      p_value_heterogeneity = p_lr)
  })



#### trend test 
trend_base_sensi_1 <- data.frame(explanatory = character(),
                                model = factor(), 
                                p_value_trend = numeric(), 
                                stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1 <- rbind(trend_base_sensi_1, 
                             data.frame(explanatory = var,
                                        model = "base", 
                                        analysis = "sensi_1",
                                        p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)




### model 1 gams ----

model1_gam_sensi_1 <- list()

for (var in proteomic_sensi_1) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1)
  model_summary <- summary(model)
  model1_gam_sensi_1[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### model 2 sd ----

model2_sd_sensi_1 <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1)

  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1 <- rbind(model2_sd_sensi_1, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model2_sd_sensi_1 <- model2_sd_sensi_1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1") |>
  select(explanatory, model, everything())


rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1 <- rbind(model2_quart_sensi_1, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model2_quart_sensi_1 <- model2_quart_sensi_1 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests
heterogeneity_adjusted_sensi_1 <- map_dfr(
  proteomic_quart_sensi_1,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1",
      p_value_heterogeneity = p_lr)
  })


### trend tests

trend_adjusted_sensi_1 <- data.frame(explanatory = character(),
                                    model = factor(), 
                                    p_value_trend = numeric(), 
                                    stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1 <- rbind(trend_adjusted_sensi_1, 
                                 data.frame(explanatory = var,
                                            model = "adjusted", 
                                            analysis = "sensi_1",
                                            p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)



### model 2 gams ----

model2_gam_sensi_1 <- list()

for (var in proteomic_sensi_1) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1)
  model_summary <- summary(model)
  model2_gam_sensi_1[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



# Sensi 2 - Filtering to cases and their controls with follow-up < 5 years ----
### model 1 sd ----
model1_sd_sensi_2 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_2 <- rbind(model1_sd_sensi_2, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model1_sd_sensi_2 <- model1_sd_sensi_2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----
model1_quart_sensi_2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_2 <- rbind(model1_quart_sensi_2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_quart_sensi_2 <- model1_quart_sensi_2 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_2Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_2Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_2Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_2 <- data.frame(explanatory = character(),
                                           model = factor(),
                                           p_value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish_sensi_2)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_2)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_sensi_2 <- rbind(heterogeneity_base_sensi_2, 
                                        data.frame(explanatory = var,
                                                   model = "base", 
                                                   analysis = "sensi_2",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base_sensi_2 <- data.frame(explanatory = character(),
                                   model = factor(), 
                                   p_value_trend = numeric(), 
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_2) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_2 <- rbind(trend_base_sensi_2, 
                                data.frame(explanatory = var,
                                           model = "base", 
                                           analysis = "sensi_2",
                                           p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 


model1_gam_sensi_2 <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_2)
  model_summary <- summary(model)
  model1_gam_sensi_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



### model 2 sd ----


model2_sd_sensi_2 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_2 <- rbind(model2_sd_sensi_2, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model2_sd_sensi_2 <- model2_sd_sensi_2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_2)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_2 <- rbind(model2_quart_sensi_2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_quart_sensi_2 <- model2_quart_sensi_2 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_2Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_2Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_2Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_sensi_2 <- data.frame(explanatory = character(),
                                               model = factor(), 
                                               p_value_heterogeneity = numeric(), 
                                               stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi, data = bdd_danish_sensi_2)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_2)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_sensi_2 <- rbind(heterogeneity_adjusted_sensi_2, 
                                            data.frame(explanatory = var,
                                                       model = "adjusted", 
                                                       analysis = "sensi_2",
                                                       p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted_sensi_2 <- data.frame(explanatory = character(),
                                       model = factor(), 
                                       p_value_trend = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_2) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_2 <- rbind(trend_adjusted_sensi_2, 
                                    data.frame(explanatory = var,
                                               model = "adjusted", 
                                               analysis = "sensi_2",
                                               p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 2 gams ---- 

model2_gam_sensi_2 <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_2)
  model_summary <- summary(model)
  model2_gam_sensi_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


# Sensi 3 - Filtering to cases and their controls with follow_up > 5 years ----
### model 1 sd ----

model1_sd_sensi_3 <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_3 <- rbind(model1_sd_sensi_3, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model1_sd_sensi_3 <- model1_sd_sensi_3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_3 <- rbind(model1_quart_sensi_3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_quart_sensi_3 <- model1_quart_sensi_3 |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_3 <- data.frame(explanatory = character(),
                                         model = factor(),
                                         p_value_heterogeneity = numeric(), 
                                         stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish_sensi_3)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_3)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_sensi_3 <- rbind(heterogeneity_base_sensi_3, 
                                      data.frame(explanatory = var,
                                                 model = "base", 
                                                 analysis = "sensi_3",
                                                 p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base_sensi_3 <- data.frame(explanatory = character(),
                                 model = factor(), 
                                 p_value_trend = numeric(), 
                                 stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_3) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_3 <- rbind(trend_base_sensi_3, 
                              data.frame(explanatory = var,
                                         model = "base", 
                                         analysis = "sensi_3",
                                         p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)



### model 1 gams ----
model1_gam_sensi_3 <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_3)
  model_summary <- summary(model)
  model1_gam_sensi_3[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




### model 2 sd ----

model2_sd_sensi_3 <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_3 <- rbind(model2_sd_sensi_3, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model2_sd_sensi_3 <- model2_sd_sensi_3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 2 quartiles ----
model2_quart_sensi_3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_3)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_3 <- rbind(model2_quart_sensi_3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_quart_sensi_3 <- model2_quart_sensi_3 |> 
  mutate(
    term = case_when(
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_sensi_3 <- data.frame(explanatory = character(),
                                             model = factor(), 
                                             p_value_heterogeneity = numeric(), 
                                             stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi, data = bdd_danish_sensi_3)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_3)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_sensi_3 <- rbind(heterogeneity_adjusted_sensi_3, 
                                          data.frame(explanatory = var,
                                                     model = "adjusted", 
                                                     analysis = "sensi_3",
                                                     p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted_sensi_3 <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_3) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_3 <- rbind(trend_adjusted_sensi_3, 
                                  data.frame(explanatory = var,
                                             model = "adjusted", 
                                             analysis = "sensi_3",
                                             p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 2 gams ----
model2_gam_sensi_3 <- list()

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_3)
  model_summary <- summary(model)
  model2_gam_sensi_3[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



# Sensi 1 + sensi 3 - Removing the oulier for NfL + filtering cases and their controls with follow_up > 5 years ----
### model 1 sd ----

model1_sd_sensi_1_3 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3 <- rbind(model1_sd_sensi_1_3, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model1_sd_sensi_1_3 <- model1_sd_sensi_1_3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_3 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_3 <- rbind(model1_quart_sensi_1_3, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_quart_sensi_1_3 <- model1_quart_sensi_1_3 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_3 <- map_dfr(
  proteomic_quart_sensi_1_3,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3",
      p_value_heterogeneity = p_lr)
  })


#### trend test 
trend_base_sensi_1_3 <- data.frame(explanatory = character(),
                                   model = factor(), 
                                   p_value_trend = numeric(), 
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_3) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_3 <- rbind(trend_base_sensi_1_3, 
                                data.frame(explanatory = var,
                                           model = "base", 
                                           analysis = "sensi_1_3",
                                           p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 

model1_gam_sensi_1_3 <- list()

for (var in proteomic_sensi_1_3) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3)
  model_summary <- summary(model)
  model1_gam_sensi_1_3[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




### model 2 sd ----

model2_sd_sensi_1_3 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3 <- rbind(model2_sd_sensi_1_3, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model2_sd_sensi_1_3 <- model2_sd_sensi_1_3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_3 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_3)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_3 <- rbind(model2_quart_sensi_1_3, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_quart_sensi_1_3 <- model2_quart_sensi_1_3 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_sensi_1_3 <- map_dfr(
  proteomic_quart_sensi_1_3,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3",
      p_value_heterogeneity = p_lr)
  })



### trend tests

trend_adjusted_sensi_1_3 <- data.frame(explanatory = character(),
                                       model = factor(), 
                                       p_value_trend = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_3) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3 <- rbind(trend_adjusted_sensi_1_3, 
                                    data.frame(explanatory = var,
                                               model = "adjusted", 
                                               analysis = "sensi_1_3",
                                               p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)



### model 2 gams ---- 

model2_gam_sensi_1_3 <- list()

for (var in proteomic_sensi_1_3) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3)
  model_summary <- summary(model)
  model2_gam_sensi_1_3[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




# Sensi 1 + sensi 3 + sensi 4 - Removing NfL outlier + filtering cases and their controls with follow- up < 5 years + filtering follow-up <= 50%----
### model 1 sd ----

model1_sd_sensi_1_3_4 <- data.frame(explanatory = character(),
                                    term = integer(),
                                    OR = numeric(),
                                    lower_CI = numeric(),
                                    upper_CI = numeric(),
                                    p_value = numeric(),
                                    stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_4)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3_4 <- rbind(model1_sd_sensi_1_3_4, data.frame(explanatory = var,
                                                                   term = term, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   p_value = p_value))
}

model1_sd_sensi_1_3_4 <- model1_sd_sensi_1_3_4 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_4") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_3_4 <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_4)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_3_4 <- rbind(model1_quart_sensi_1_3_4, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model1_quart_sensi_1_3_4 <- model1_quart_sensi_1_3_4 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3_4Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3_4Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3_4Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_4") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_3_4 <- map_dfr(
  proteomic_quart_sensi_1_3_4,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3_4 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3_4")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3_4",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3_4",
      p_value_heterogeneity = p_lr)
  })


#### trend test 
trend_base_sensi_1_3_4 <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_3_4) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_3_4 <- rbind(trend_base_sensi_1_3_4, 
                                  data.frame(explanatory = var,
                                             model = "base", 
                                             analysis = "sensi_1_3_4",
                                             p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 

model1_gam_sensi_1_3_4 <- list()

for (var in proteomic_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3_4)
  model_summary <- summary(model)
  model1_gam_sensi_1_3_4[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




### model 2 sd ----

model2_sd_sensi_1_3_4 <- data.frame(explanatory = character(),
                                    term = integer(),
                                    OR = numeric(),
                                    lower_CI = numeric(),
                                    upper_CI = numeric(),
                                    p_value = numeric(),
                                    stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_4)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3_4 <- rbind(model2_sd_sensi_1_3_4, data.frame(explanatory = var,
                                                                   term = term, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   p_value = p_value))
}

model2_sd_sensi_1_3_4 <- model2_sd_sensi_1_3_4 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_4") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_3_4 <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_3_4)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_3_4 <- rbind(model2_quart_sensi_1_3_4, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model2_quart_sensi_1_3_4 <- model2_quart_sensi_1_3_4 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3_4Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3_4Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3_4Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_4") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests
heterogeneity_adjusted_sensi_1_3_4 <- map_dfr(
  proteomic_quart_sensi_1_3_4,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3_4 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3_4")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3_4",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3_4",
      p_value_heterogeneity = p_lr)
  })



### trend tests

trend_adjusted_sensi_1_3_4 <- data.frame(explanatory = character(),
                                         model = factor(), 
                                         p_value_trend = numeric(), 
                                         stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_3_4) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3_4 <- rbind(trend_adjusted_sensi_1_3_4, 
                                      data.frame(explanatory = var,
                                                 model = "adjusted", 
                                                 analysis = "sensi_1_3_4",
                                                 p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 2 gams ---- 

model2_gam_sensi_1_3_4 <- list()

for (var in proteomic_sensi_1_3_4) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3_4)
  model_summary <- summary(model)
  model2_gam_sensi_1_3_4[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



# Sensi 1 + sensi 3 + sensi 5 - Removing NfL outlier + filtering cases and their controls with follow-up < 5 years + filtering follow-up > 50%----
### model 1 sd ----

model1_sd_sensi_1_3_5 <- data.frame(explanatory = character(),
                                    term = integer(),
                                    OR = numeric(),
                                    lower_CI = numeric(),
                                    upper_CI = numeric(),
                                    p_value = numeric(),
                                    stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_5)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3_5 <- rbind(model1_sd_sensi_1_3_5, data.frame(explanatory = var,
                                                                   term = term, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   p_value = p_value))
}

model1_sd_sensi_1_3_5 <- model1_sd_sensi_1_3_5 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_5") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_3_5 <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_5)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_3_5 <- rbind(model1_quart_sensi_1_3_5, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model1_quart_sensi_1_3_5 <- model1_quart_sensi_1_3_5 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3_5Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3_5Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3_5Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_5") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_3_5 <- map_dfr(
  proteomic_quart_sensi_1_3_5,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3_5 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3_5")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_3_5",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_3_5",
      p_value_heterogeneity = p_lr)
  })

#### trend test 
trend_base_sensi_1_3_5 <- data.frame(explanatory = character(),
                                     model = factor(), 
                                     p_value_trend = numeric(), 
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_3_5) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_3_5 <- rbind(trend_base_sensi_1_3_5, 
                                  data.frame(explanatory = var,
                                             model = "base", 
                                             analysis = "sensi_1_3_5",
                                             p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 


model1_gam_sensi_1_3_5 <- list()

for (var in proteomic_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3_5)
  model_summary <- summary(model)
  model1_gam_sensi_1_3_5[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



### model 2 sd ----


model2_sd_sensi_1_3_5 <- data.frame(explanatory = character(),
                                    term = integer(),
                                    OR = numeric(),
                                    lower_CI = numeric(),
                                    upper_CI = numeric(),
                                    p_value = numeric(),
                                    stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_5)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3_5 <- rbind(model2_sd_sensi_1_3_5, data.frame(explanatory = var,
                                                                   term = term, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   p_value = p_value))
}

model2_sd_sensi_1_3_5 <- model2_sd_sensi_1_3_5 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_5") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_3_5 <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_3_5)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_3_5 <- rbind(model2_quart_sensi_1_3_5, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model2_quart_sensi_1_3_5 <- model2_quart_sensi_1_3_5 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_3_5Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_3_5Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_3_5Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_5") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_sensi_1_3_5 <- map_dfr(
  proteomic_quart_sensi_1_3_5,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_3_5 |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_3_5")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_3_5",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_3_5",
      p_value_heterogeneity = p_lr)
  })



### trend tests

trend_adjusted_sensi_1_3_5 <- data.frame(explanatory = character(),
                                         model = factor(), 
                                         p_value_trend = numeric(), 
                                         stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_3_5) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_3_5 <- rbind(trend_adjusted_sensi_1_3_5, 
                                      data.frame(explanatory = var,
                                                 model = "adjusted", 
                                                 analysis = "sensi_1_3_5",
                                                 p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)



### model 2 gams ---- 

model2_gam_sensi_1_3_5 <- list()

for (var in proteomic_sensi_1_3_5) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_3_5)
  model_summary <- summary(model)
  model2_gam_sensi_1_3_5[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




# Sensi 6 - stratifiyng the analysis in tertiles of follow-up duration  ----
### Visualization ----
sensi_6_table_follow_up <- bdd_danish_sensi_6 |>                                # for visualisation of follow up distribution 
  select(follow_up_years, follow_up_ter) |>
  tbl_summary(by = follow_up_ter)

sensi_6_densityplot_follow_up <-                                                # for visualisation of follow up distribution 
  bdd_danish_sensi_6 |>
  filter(als ==1) |>
  ggplot() +
  aes(x = follow_up_years) +
  geom_density(fill = "lightgray") +
  geom_vline(
    xintercept = quantile(bdd_danish_sensi_6$follow_up_years, probs = c(1/3, 2/3), na.rm = TRUE),
    color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = min(bdd_danish_sensi_6$follow_up_years, na.rm = TRUE),
           y = Inf, label = "Tertile 1", color = "red", hjust = 0, vjust = 2) +
  annotate("text",
           x = median(bdd_danish_sensi_6$follow_up_years, na.rm = TRUE),  # centre entre les deux lignes
           y = Inf, label = "Tertile 2", color = "red", vjust = 2) +
  annotate("text",
           x = max(bdd_danish_sensi_6$follow_up_years, na.rm = TRUE),
           y = Inf, label = "Tertile 3", color = "red", hjust = 1, vjust = 2) +
  labs(x = "Follow-up duration from baseline to diagnosis (years)")

bdd_danish_sensi_6 <- bdd_danish |>
  group_by(match) |>
  mutate(follow_up_ter = follow_up[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE, 
      labels = c("Tertile 1", "Tertile 2", "Tertile 3")))
bdd_danish_sensi_6_T1 <-                                                        # stratifiyng the data on tertiles of follow up duration 
  bdd_danish_sensi_6 |>
  filter(follow_up_ter == "Tertile 1")
bdd_danish_sensi_6_T2 <- 
  bdd_danish_sensi_6 |>
  filter(follow_up_ter == "Tertile 2")
bdd_danish_sensi_6_T3 <- 
  bdd_danish_sensi_6 |>
  filter(follow_up_ter == "Tertile 3")

### model 1 T1 ----
model1_sd_sensi_6_T1 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_6_T1 <- rbind(model1_sd_sensi_6_T1, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_6_T1 <- model1_sd_sensi_6_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_6_T1") |>
  select(explanatory, model, everything())

### model 1 T2 ----
model1_sd_sensi_6_T2 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_6_T2 <- rbind(model1_sd_sensi_6_T2, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_6_T2 <- model1_sd_sensi_6_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_6_T2") |>
  select(explanatory, model, everything())


### model 1 T3 ----
model1_sd_sensi_6_T3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_6_T3 <- rbind(model1_sd_sensi_6_T3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_6_T3 <- model1_sd_sensi_6_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_6_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



### model 2 T1 ----
model2_sd_sensi_6_T1 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T1)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_6_T1 <- rbind(model2_sd_sensi_6_T1, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_6_T1 <- model2_sd_sensi_6_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_6_T1") |>
  select(explanatory, model, everything())

### model 2 T2 ----
model2_sd_sensi_6_T2 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_6_T2 <- rbind(model2_sd_sensi_6_T2, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_6_T2 <- model2_sd_sensi_6_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_6_T2") |>
  select(explanatory, model, everything())


### model 2 T3 ----
model2_sd_sensi_6_T3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_6_T3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_6_T3 <- rbind(model2_sd_sensi_6_T3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_6_T3 <- model2_sd_sensi_6_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_6_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var, 
   bdd_danish_sensi_6, bdd_danish_sensi_6_T1, bdd_danish_sensi_6_T2, bdd_danish_sensi_6_T3)




# Sensi 1 + sensi 6 - Removing the oulier for NfL + stratifiyng the analysis in tertiles of follow-up duration  ----
### model 1 T1 ----
model1_sd_sensi_1_6_T1 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_6_T1 <- rbind(model1_sd_sensi_1_6_T1, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_sd_sensi_1_6_T1 <- model1_sd_sensi_1_6_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_6_T1") |>
  select(explanatory, model, everything())


### model 1 T2 ----
model1_sd_sensi_1_6_T2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_6_T2 <- rbind(model1_sd_sensi_1_6_T2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_sd_sensi_1_6_T2 <- model1_sd_sensi_1_6_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_6_T2") |>
  select(explanatory, model, everything())



### model 1 T3 ----
model1_sd_sensi_1_6_T3 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)


for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_6_T3 <- rbind(model1_sd_sensi_1_6_T3, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_sd_sensi_1_6_T3 <- model1_sd_sensi_1_6_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_6_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



### model 2 T1 ----
model2_sd_sensi_1_6_T1 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T1)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_6_T1 <- rbind(model2_sd_sensi_1_6_T1, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_6_T1 <- model2_sd_sensi_1_6_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_6_T1") |>
  select(explanatory, model, everything())


### model 2 T2 ----
model2_sd_sensi_1_6_T2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_6_T2 <- rbind(model2_sd_sensi_1_6_T2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_6_T2 <- model2_sd_sensi_1_6_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_6_T2") |>
  select(explanatory, model, everything())



### model 2 T3 ----
model2_sd_sensi_1_6_T3 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)


for (var in proteomic_sd_sensi_1_6) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_6_T3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_6_T3 <- rbind(model2_sd_sensi_1_6_T3, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_6_T3 <- model2_sd_sensi_1_6_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_6_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var, 
   bdd_danish_sensi_1_6, bdd_danish_sensi_1_6_T1, bdd_danish_sensi_1_6_T2, bdd_danish_sensi_1_6_T3, 
   proteomic_sd_sensi_1_6)


# Sensi 1 + sensi 7 - Removing the oulier for NfL + filtering to females ----
### model 1 sd ----

model1_sd_sensi_1_7_female <- data.frame(explanatory = character(),
                                         term = integer(),
                                         OR = numeric(),
                                         lower_CI = numeric(),
                                         upper_CI = numeric(),
                                         p_value = numeric(),
                                         stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_female)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_7_female <- rbind(model1_sd_sensi_1_7_female, data.frame(explanatory = var,
                                                                             term = term, 
                                                                             OR = OR,
                                                                             lower_CI = lower_CI,
                                                                             upper_CI = upper_CI,
                                                                             p_value = p_value))
}

model1_sd_sensi_1_7_female <- model1_sd_sensi_1_7_female |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_7_female") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_7_female <- data.frame(explanatory = character(),
                                            term = integer(),
                                            OR = numeric(),
                                            lower_CI = numeric(),
                                            upper_CI = numeric(),
                                            p_value = numeric(),
                                            stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_female)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_7_female <- rbind(model1_quart_sensi_1_7_female, data.frame(explanatory = var,
                                                                                   term = term, 
                                                                                   OR = OR,
                                                                                   lower_CI = lower_CI,
                                                                                   upper_CI = upper_CI,
                                                                                   p_value = p_value))
}

model1_quart_sensi_1_7_female <- model1_quart_sensi_1_7_female |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_7_femaleQ2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_7_femaleQ3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_7_femaleQ4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_7_female") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_7_female <- map_dfr(
  proteomic_quart_sensi_1_7_female,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_7_female |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_7_female")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_7_female",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_7_female",
      p_value_heterogeneity = p_lr)
  })



#### trend test 
trend_base_sensi_1_7_female <- data.frame(explanatory = character(),
                                          model = factor(), 
                                          p_value_trend = numeric(), 
                                          stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_7_female) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_7_female <- rbind(trend_base_sensi_1_7_female, 
                                       data.frame(explanatory = var,
                                                  model = "base", 
                                                  analysis = "sensi_1_7_female",
                                                  p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 
model1_gam_sensi_1_7_female <- list()

for (var in proteomic_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_7_female)
  model_summary <- summary(model)
  model1_gam_sensi_1_7_female[[var]] <- model_summary
}

rm(var, formula, model, model_summary)



### model 2 sd ----
model2_sd_sensi_1_7_female <- data.frame(explanatory = character(),
                                         term = integer(),
                                         OR = numeric(),
                                         lower_CI = numeric(),
                                         upper_CI = numeric(),
                                         p_value = numeric(),
                                         stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_female)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_7_female <- rbind(model2_sd_sensi_1_7_female, data.frame(explanatory = var,
                                                                             term = term, 
                                                                             OR = OR,
                                                                             lower_CI = lower_CI,
                                                                             upper_CI = upper_CI,
                                                                             p_value = p_value))
}

model2_sd_sensi_1_7_female <- model2_sd_sensi_1_7_female |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_7_female") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_7_female <- data.frame(explanatory = character(),
                                            term = integer(),
                                            OR = numeric(),
                                            lower_CI = numeric(),
                                            upper_CI = numeric(),
                                            p_value = numeric(),
                                            stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_7_female)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_7_female <- rbind(model2_quart_sensi_1_7_female, data.frame(explanatory = var,
                                                                                   term = term, 
                                                                                   OR = OR,
                                                                                   lower_CI = lower_CI,
                                                                                   upper_CI = upper_CI,
                                                                                   p_value = p_value))
}

model2_quart_sensi_1_7_female <- model2_quart_sensi_1_7_female |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_7_femaleQ2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_7_femaleQ3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_7_femaleQ4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_7_female") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests
heterogeneity_adjusted_sensi_1_7_female <- map_dfr(
  proteomic_quart_sensi_1_7_female,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_7_female |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_7_female")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_7_female",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_7_female",
      p_value_heterogeneity = p_lr)
  })



### trend tests

trend_adjusted_sensi_1_7_female <- data.frame(explanatory = character(),
                                              model = factor(), 
                                              p_value_trend = numeric(), 
                                              stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_7_female) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_7_female <- rbind(trend_adjusted_sensi_1_7_female, 
                                           data.frame(explanatory = var,
                                                      model = "adjusted", 
                                                      analysis = "sensi_1_7_female",
                                                      p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)



### model 2 gams ---- 

model2_gam_sensi_1_7_female <- list()

for (var in proteomic_sensi_1_7_female) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_7_female)
  model_summary <- summary(model)
  model2_gam_sensi_1_7_female[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




# Sensi 1 + sensi 7 - Removing the oulier for NfL + filtering to male ----
### model 1 sd ----

model1_sd_sensi_1_7_male <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_male)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_7_male <- rbind(model1_sd_sensi_1_7_male, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model1_sd_sensi_1_7_male <- model1_sd_sensi_1_7_male |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_7_male") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_7_male <- data.frame(explanatory = character(),
                                          term = integer(),
                                          OR = numeric(),
                                          lower_CI = numeric(),
                                          upper_CI = numeric(),
                                          p_value = numeric(),
                                          stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_male)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_7_male <- rbind(model1_quart_sensi_1_7_male, data.frame(explanatory = var,
                                                                               term = term, 
                                                                               OR = OR,
                                                                               lower_CI = lower_CI,
                                                                               upper_CI = upper_CI,
                                                                               p_value = p_value))
}

model1_quart_sensi_1_7_male <- model1_quart_sensi_1_7_male |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_7_maleQ2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_7_maleQ3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_7_maleQ4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_7_male") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_7_male <- map_dfr(
  proteomic_quart_sensi_1_7_male,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_7_male |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_7_male")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "base",
        analysis = "sensi_1_7_male",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match),
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "base",
      analysis = "sensi_1_7_male",
      p_value_heterogeneity = p_lr)
  })



#### trend test 
trend_base_sensi_1_7_male <- data.frame(explanatory = character(),
                                        model = factor(), 
                                        p_value_trend = numeric(), 
                                        stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_7_male) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_7_male <- rbind(trend_base_sensi_1_7_male, 
                                     data.frame(explanatory = var,
                                                model = "base", 
                                                analysis = "sensi_1_7_male",
                                                p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 

model1_gam_sensi_1_7_male <- list()

for (var in proteomic_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + baseline_age"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_7_male)
  model_summary <- summary(model)
  model1_gam_sensi_1_7_male[[var]] <- model_summary
}

rm(var, formula, model, model_summary)




### model 2 sd ----

model2_sd_sensi_1_7_male <- data.frame(explanatory = character(),
                                       term = integer(),
                                       OR = numeric(),
                                       lower_CI = numeric(),
                                       upper_CI = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_7_male)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_7_male <- rbind(model2_sd_sensi_1_7_male, data.frame(explanatory = var,
                                                                         term = term, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         p_value = p_value))
}

model2_sd_sensi_1_7_male <- model2_sd_sensi_1_7_male |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_7_male") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_7_male <- data.frame(explanatory = character(),
                                          term = integer(),
                                          OR = numeric(),
                                          lower_CI = numeric(),
                                          upper_CI = numeric(),
                                          p_value = numeric(),
                                          stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_7_male)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_7_male <- rbind(model2_quart_sensi_1_7_male, data.frame(explanatory = var,
                                                                               term = term, 
                                                                               OR = OR,
                                                                               lower_CI = lower_CI,
                                                                               upper_CI = upper_CI,
                                                                               p_value = p_value))
}

model2_quart_sensi_1_7_male <- model2_quart_sensi_1_7_male |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_7_maleQ2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_7_maleQ3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_7_maleQ4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_7_male") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



### heterogeneity tests

heterogeneity_adjusted_sensi_1_7_male <- map_dfr(
  proteomic_quart_sensi_1_7_male,
  function(var) {
    
    data_sub <- bdd_danish_sensi_1_7_male |>
      filter(
        if (var == "proteomic_neuro_explo_NEFL_sensi_1_7_male")
          match != 159
        else
          TRUE) |>
      
      filter(
        complete.cases(
          pick(als, match, smoking_2cat_i, bmi, all_of(var))))
    
    if (nrow(data_sub) == 0) {
      return(tibble(
        explanatory = var,
        model = "adjusted",
        analysis = "sensi_1_7_male",
        p_value_heterogeneity = NA_real_))
    }
    
    test_1 <- clogit(
      als ~ strata(match) + smoking_2cat_i + bmi,
      data = data_sub)
    
    test_2 <- clogit(
      reformulate(
        termlabels = c(var, "strata(match)", "smoking_2cat_i", "bmi"),
        response = "als"),
      data = data_sub)
    
    p_lr <- anova(test_1, test_2, test = "LR")$`Pr(>|Chi|)`[2]
    
    tibble(
      explanatory = var,
      model = "adjusted",
      analysis = "sensi_1_7_male",
      p_value_heterogeneity = p_lr)
  })




### trend tests

trend_adjusted_sensi_1_7_male <- data.frame(explanatory = character(),
                                            model = factor(), 
                                            p_value_trend = numeric(), 
                                            stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_7_male) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_7_male <- rbind(trend_adjusted_sensi_1_7_male, 
                                         data.frame(explanatory = var,
                                                    model = "adjusted", 
                                                    analysis = "sensi_1_7_male",
                                                    p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 2 gams ---- 

model2_gam_sensi_1_7_male <- list()

for (var in proteomic_sensi_1_7_male) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + baseline_age + smoking_2cat_i + bmi"))
  model <- mgcv::gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_7_male)
  model_summary <- summary(model)
  model2_gam_sensi_1_7_male[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


# Merging the main results ----
main_results <- bind_rows(model1_sd,                                            # main analyses
                          model2_sd, 
                          model1_quart, 
                          model2_quart, 
                         
                          model1_sd_sensi_1,                                    # sensi 1 
                          model2_sd_sensi_1, 
                          model1_quart_sensi_1, 
                          model2_quart_sensi_1, 
                          
                          model1_sd_sensi_2,                                    # sensi 2
                          model2_sd_sensi_2,
                          model1_quart_sensi_2,
                          model2_quart_sensi_2,
                          
                          model1_sd_sensi_3,                                    # sensi 3
                          model2_sd_sensi_3,
                          model1_quart_sensi_3,
                          model2_quart_sensi_3,
                          
                          model1_sd_sensi_1_3,                                  # sensi 1_3
                          model2_sd_sensi_1_3,
                          model1_quart_sensi_1_3,
                          model2_quart_sensi_1_3,
                          
                          model1_sd_sensi_1_3_4,                                # sensi 1_3_4
                          model2_sd_sensi_1_3_4,
                          model1_quart_sensi_1_3_4,
                          model2_quart_sensi_1_3_4,
                          
                          model1_sd_sensi_1_3_5,                                # sensi 1_3_5
                          model2_sd_sensi_1_3_5,
                          model1_quart_sensi_1_3_5,
                          model2_quart_sensi_1_3_5,

                          model1_sd_sensi_6_T1,                                 # sensi 6
                          model1_sd_sensi_6_T2,
                          model1_sd_sensi_6_T3,
                          model2_sd_sensi_6_T1,
                          model2_sd_sensi_6_T2,
                          model2_sd_sensi_6_T3,

                          model1_sd_sensi_1_6_T1,                               # sensi 1_6
                          model1_sd_sensi_1_6_T2,
                          model1_sd_sensi_1_6_T3,
                          model2_sd_sensi_1_6_T1,
                          model2_sd_sensi_1_6_T2,
                          model2_sd_sensi_1_6_T3, 
                          
                          model1_sd_sensi_1_7_female,                           # sensi 1_7 female 
                          model2_sd_sensi_1_7_female,
                          model1_quart_sensi_1_7_female,
                          model2_quart_sensi_1_7_female, 
                          
                          model1_sd_sensi_1_7_male,                             # sensi 1_7 male 
                          model2_sd_sensi_1_7_male,
                          model1_quart_sensi_1_7_male,
                          model2_quart_sensi_1_7_male) |>
  
  mutate(explanatory = if_else(str_detect(explanatory, "NEFL"), "proteomic_neuro_explo_NEFL", explanatory), 
         
         explanatory = gsub("_quart_med", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory), 
         explanatory = gsub("_sd", "", explanatory), 
         
         OR_raw = OR, 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p_value_raw = p_value, 
         p_value = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
         "95% CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  group_by(model) |>                               
  mutate(
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      NA_real_), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_1" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_2" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_1_3" & model == "base",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "main" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_1" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_2" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr), 
    p_value_fdr = if_else(
      term == "Continuous" & analysis == "sensi_1_3" & model == "adjusted",                            
      p.adjust(p_value_raw, method = "fdr"),
      p_value_fdr)) |>
  ungroup() |>
  arrange(explanatory) |>
  select(analysis, 
         model,
         explanatory, 
         term,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p_value"), 
         lower_CI, upper_CI) 



heterogeneity_tests <-                                                          # main
  bind_rows(heterogeneity_base, heterogeneity_adjusted)  |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests <-      
  bind_rows(trend_base, trend_adjusted) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1 <-                                                  # sensi 1 
  bind_rows(heterogeneity_base_sensi_1, heterogeneity_adjusted_sensi_1) |>
  mutate(explanatory = gsub("_quart_sensi_1", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1 <- 
  bind_rows(trend_base_sensi_1, trend_adjusted_sensi_1) |>
  mutate(explanatory = gsub("_quart_med_sensi_1", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_2 <-                                                  # sensi 2 
  bind_rows(heterogeneity_base_sensi_2, heterogeneity_adjusted_sensi_2) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_2 <- 
  bind_rows(trend_base_sensi_2, trend_adjusted_sensi_2) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_3 <-                                                  # sensi 3
  bind_rows(heterogeneity_base_sensi_3, heterogeneity_adjusted_sensi_3) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_3 <- 
  bind_rows(trend_base_sensi_3, trend_adjusted_sensi_3) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_3 <-                                                # sensi 1 3 
  bind_rows(heterogeneity_base_sensi_1_3, heterogeneity_adjusted_sensi_1_3) |>
  mutate(explanatory = gsub("_quart_sensi_1_3", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1_3 <- 
  bind_rows(trend_base_sensi_1_3, trend_adjusted_sensi_1_3) |>
  mutate(explanatory = gsub("_quart_med_sensi_1_3", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_3_4 <-                                              # sensi 1 3 4 
  bind_rows(heterogeneity_base_sensi_1_3_4, heterogeneity_adjusted_sensi_1_3_4) |>
  mutate(explanatory = gsub("_quart_sensi_1_3_4", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1_3_4 <- 
  bind_rows(trend_base_sensi_1_3_4, trend_adjusted_sensi_1_3_4) |>
  mutate(explanatory = gsub("_quart_med_sensi_1_3_4", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_3_5 <-                                              # sensi 1 3 5 
  bind_rows(heterogeneity_base_sensi_1_3_5, heterogeneity_adjusted_sensi_1_3_5) |>
  mutate(explanatory = gsub("_quart_sensi_1_3_5", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1_3_5 <- 
  bind_rows(trend_base_sensi_1_3_5, trend_adjusted_sensi_1_3_5) |>
  mutate(explanatory = gsub("_quart_med_sensi_1_3_5", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_7_female <-                                         # sensi 1 7 female
  bind_rows(heterogeneity_base_sensi_1_7_female, heterogeneity_adjusted_sensi_1_7_female) |>
  mutate(explanatory = gsub("_quart_sensi_1_7_female", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory)) 

trend_tests_sensi_1_7_female <- 
  bind_rows(trend_base_sensi_1_7_female, trend_adjusted_sensi_1_7_female) |>
  mutate(explanatory = gsub("_quart_med_sensi_1_7_female", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_7_male <-                                           # sensi 1 7 male
  bind_rows(heterogeneity_base_sensi_1_7_male, heterogeneity_adjusted_sensi_1_7_male) |>
  mutate(explanatory = gsub("_quart_sensi_1_7_male", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1_7_male <- 
  bind_rows(trend_base_sensi_1_7_male, trend_adjusted_sensi_1_7_male) |>
  mutate(explanatory = gsub("_quart_med_sensi_1_7_male", "", explanatory), 
         explanatory = gsub("_quart_med", "", explanatory))


heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_2)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_3)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_3)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_3_4)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_3_5)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_7_female)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_7_male)

trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_2)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_3)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_3)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_3_4)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_3_5)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_7_female)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_7_male)


main_results <- left_join(main_results, heterogeneity_tests, by = c("explanatory", "model", "analysis"))
main_results <- left_join(main_results, trend_tests, by = c("explanatory", "model", "analysis"))

main_results <- main_results |>
  mutate(
    p_value_heterogeneity = ifelse(term == "Continuous", NA, p_value_heterogeneity), 
    p_value_trend = ifelse(term == "Continuous", NA, p_value_trend), 
    p_value_heterogeneity = ifelse(p_value_heterogeneity < 0.01, "<0.01", number(p_value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p_value_trend = ifelse(p_value_trend < 0.01, "<0.01", number(p_value_trend, accuracy = 0.01, decimal.mark = ".")), 
    protein_group = case_when(str_detect(explanatory, 'proteomic_immun_res') ~ "Immune response", 
                              str_detect(explanatory, 'proteomic_metabolism') ~ "Metabolism", 
                              str_detect(explanatory, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
    explanatory = str_replace(explanatory, 'proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_', ""))


rm(model1_sd, model1_quart, model2_sd, model2_quart, 
   
   model1_sd_sensi_1, model1_quart_sensi_1, model2_sd_sensi_1, model2_quart_sensi_1, 
   model1_sd_sensi_2, model1_quart_sensi_2, model2_sd_sensi_2, model2_quart_sensi_2, 
   model1_sd_sensi_3, model1_quart_sensi_3, model2_sd_sensi_3, model2_quart_sensi_3, 
   model1_sd_sensi_1_3, model1_quart_sensi_1_3, model2_sd_sensi_1_3, model2_quart_sensi_1_3, 
   model1_sd_sensi_1_3_4, model1_quart_sensi_1_3_4, model2_sd_sensi_1_3_4, model2_quart_sensi_1_3_4, 
   model1_sd_sensi_1_3_5, model1_quart_sensi_1_3_5, model2_sd_sensi_1_3_5, model2_quart_sensi_1_3_5, 
   model1_sd_sensi_1_7_female, model1_quart_sensi_1_7_female, model2_sd_sensi_1_7_female, model2_quart_sensi_1_7_female, 
   model1_sd_sensi_1_7_male, model1_quart_sensi_1_7_male, model2_sd_sensi_1_7_male, model2_quart_sensi_1_7_male, 

   model1_sd_sensi_6_T1, 
   model1_sd_sensi_6_T2, 
   model1_sd_sensi_6_T3, 
   model1_sd_sensi_1_6_T1, 
   model1_sd_sensi_1_6_T2, 
   model1_sd_sensi_1_6_T3, 
   
   model2_sd_sensi_6_T1, 
   model2_sd_sensi_6_T2, 
   model2_sd_sensi_6_T3, 
   model2_sd_sensi_1_6_T1, 
   model2_sd_sensi_1_6_T2, 
   model2_sd_sensi_1_6_T3, 
   
   heterogeneity_base, heterogeneity_adjusted, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests, 
   
   heterogeneity_base_sensi_1, heterogeneity_adjusted_sensi_1, 
   trend_base_sensi_1, trend_adjusted_sensi_1, 
   heterogeneity_tests_sensi_1, trend_tests_sensi_1, 
   
   heterogeneity_base_sensi_2, heterogeneity_adjusted_sensi_2, 
   trend_base_sensi_2, trend_adjusted_sensi_2, 
   heterogeneity_tests_sensi_2, trend_tests_sensi_2, 
   
   heterogeneity_base_sensi_3, heterogeneity_adjusted_sensi_3, 
   trend_base_sensi_3, trend_adjusted_sensi_3, 
   heterogeneity_tests_sensi_3, trend_tests_sensi_3, 
   
   heterogeneity_base_sensi_1_3, heterogeneity_adjusted_sensi_1_3, 
   trend_base_sensi_1_3, trend_adjusted_sensi_1_3, 
   heterogeneity_tests_sensi_1_3, trend_tests_sensi_1_3, 
   
   heterogeneity_base_sensi_1_3_4, heterogeneity_adjusted_sensi_1_3_4, 
   trend_base_sensi_1_3_4, trend_adjusted_sensi_1_3_4, 
   heterogeneity_tests_sensi_1_3_4, trend_tests_sensi_1_3_4, 
   
   heterogeneity_base_sensi_1_3_5, heterogeneity_adjusted_sensi_1_3_5, 
   trend_base_sensi_1_3_5, trend_adjusted_sensi_1_3_5, 
   heterogeneity_tests_sensi_1_3_5, trend_tests_sensi_1_3_5, 
   
   heterogeneity_base_sensi_1_7_female, heterogeneity_adjusted_sensi_1_7_female, 
   trend_base_sensi_1_7_female, trend_adjusted_sensi_1_7_female, 
   heterogeneity_tests_sensi_1_7_female, trend_tests_sensi_1_7_female, 
   
   heterogeneity_base_sensi_1_7_male, heterogeneity_adjusted_sensi_1_7_male, 
   trend_base_sensi_1_7_male, trend_adjusted_sensi_1_7_male, 
   heterogeneity_tests_sensi_1_7_male, trend_tests_sensi_1_7_male)


# Additionnal analysis 1 - Testing the normality of the distribution of all the ORs -----
OR_distribution <- 
  main_results |>
  filter(model == "adjusted") |>
  filter(analysis == "sensi_1") |>
  filter(term == "Continuous") |>
  mutate(OR_log = log(OR_raw), 
         weight = -log(p_value_raw))

w_median <- wtd.quantile(OR_distribution$OR_log, weights = OR_distribution$weight, probs = 0.5, na.rm = TRUE)


boxplot_OR <- OR_distribution|>
  ggplot() +
  aes(x = "", y = OR_log) +
  geom_boxplot(shape = "circle", fill = "#112446") +
  theme_bw() +
  theme(axis.title = element_blank()) 

densityplot_OR <- OR_distribution |>
  ggplot() +
  aes(x = OR_log) +
  geom_density(fill = "lightgray") +
  theme_bw() +
  #stat_overlay_normal_density(color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +              # normal distribution                         
  geom_vline(xintercept = median(OR_distribution$OR_log, na.rm = TRUE),         # mediane
             color = "darkgreen", 
             linetype = "dashed") +
  geom_vline(xintercept = w_median, na.rm = TRUE,                               # weighted mediane
             color = "darkblue", 
             linetype = "dashed") +
  theme(axis.title = element_blank())

exp(w_median)
median(OR_distribution$OR_raw)
shapiro.test(OR_distribution$OR_log)
rm(w_median)



# Additionnal analysis 2 - Test for multiple testing -----
# Fasano–Franceschini and Jackknife Leave-One-Out 

# Function to compute FF statistic on numeric columns only
ff_stat <- function(mat1, mat2) {
  res <- fasano.franceschini.test(mat1[, c("OR_log2", "p_value_log")],
                                  mat2[, c("OR_log2", "p_value_log")],
                                  seed = 0)
  return(res$statistic)
}


## All proteins ----
### data prep ----
data_matrix1_all <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw < 1) |>
  mutate(
    OR_log2 = abs(log2(OR_raw)),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)

data_matrix2_all <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw >= 1) |>
  mutate(
    OR_log2 = log2(OR_raw),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)



### Fasano Franceschini test ----
additional_analysis_2_all_results <- 
  fasano.franceschini.test(data_matrix1_all[, c("OR_log2", "p_value_log")],
                           data_matrix2_all[, c("OR_log2", "p_value_log")],
                           seed = 0)

### Jackknife approach ----
T0 <- additional_analysis_2_all_results$statistic                               # d statistic from the inital ff test with all observations 

influence_1 <- numeric(nrow(data_matrix1_all))                                  # Initialize vectors
influence_2 <- numeric(nrow(data_matrix2_all))

for (i in seq_len(nrow(data_matrix1_all))) {                                    # Jackknife for data_matrix1 OR_raw < 1
  mat1_minus <- data_matrix1_all[-i, ]
  Ti <- ff_stat(mat1_minus, data_matrix2_all)
  influence_1[i] <- T0 - Ti
}

for (j in seq_len(nrow(data_matrix2_all))) {                                    # Jackknife for data_matrix2 OR_raw >= 1
  mat2_minus <- data_matrix2_all[-j, ]
  Tj <- ff_stat(data_matrix1_all, mat2_minus)
  influence_2[j] <- T0 - Tj
}

jackknife_all_results <- tibble(                                                # Assemble results 
  group = c(rep("OR<1", nrow(data_matrix1_all)),
            rep("OR>=1", nrow(data_matrix2_all))),
  
  explanatory = c(data_matrix1_all$explanatory,
                  data_matrix2_all$explanatory),
  
  influence = c(influence_1, influence_2), 
  
  d_stat_percent = influence/T0*100) |>
  arrange(desc(influence))


proteomic_influence <- 
  jackknife_all_results |>                                                      # selection of the prot with highest influence to the ff test
  filter(d_stat_percent>2) |>                                                   # filter to prot that if removed, modify the D stat of at least 2%
  pull(explanatory) 
proteomic_influence <- sort(proteomic_influence)

proteomic_selected <-                                                           # selection des prot avec p<0.05
  main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         p_value_raw<0.05) |> 
  pull(explanatory) 
proteomic_selected <- sort(proteomic_selected)

valeurs_communes <- intersect(proteomic_influence, proteomic_selected)          # check the commune proteins 

### plots ----
plot1_all <- 
  data_matrix1_all |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR<1", 
       y = "-log10(p-value)", 
       x = "abslog2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

plot2_all <- 
  data_matrix2_all |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR>=1", 
       y = "-log10(p-value)", 
       x = "log2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

additional_analysis_2_all_figure <- plot1_all + plot2_all


rm(data_matrix1_all, data_matrix2_all, 
   plot1_all, plot2_all, 
   T0, influence_1, influence_2, i, j, mat1_minus, mat2_minus, Ti, Tj, 
   proteomic_influence, proteomic_selected, valeurs_communes)

## Immune response ----
### data prep ----
data_matrix1_immune <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Immune response") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw < 1) |>
  mutate(
    OR_log2 = abs(log2(OR_raw)),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)

data_matrix2_immune <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Immune response") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw >= 1) |>
  mutate(
    OR_log2 = log2(OR_raw),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)


### Fasano Franceschini test ----
additional_analysis_2_immune_results <- 
  fasano.franceschini.test(data_matrix1_immune[, c("OR_log2", "p_value_log")],
                           data_matrix2_immune[, c("OR_log2", "p_value_log")],
                           seed = 0)

### Jackknife approach ----
T0 <- additional_analysis_2_immune_results$statistic                               # d statistic from the inital ff test with all observations 

influence_1 <- numeric(nrow(data_matrix1_immune))                                      # Initialize vectors
influence_2 <- numeric(nrow(data_matrix2_immune))

for (i in seq_len(nrow(data_matrix1_immune))) {                                        # Jackknife for data_matrix1 OR_raw < 1
  mat1_minus <- data_matrix1_immune[-i, ]
  Ti <- ff_stat(mat1_minus, data_matrix2_immune)
  influence_1[i] <- T0 - Ti
}

for (j in seq_len(nrow(data_matrix2_immune))) {                                        # Jackknife for data_matrix2 OR_raw >= 1
  mat2_minus <- data_matrix2_immune[-j, ]
  Tj <- ff_stat(data_matrix1_immune, mat2_minus)
  influence_2[j] <- T0 - Tj
}

jackknife_immune_results <- tibble(                                                # Assemble results 
  group = c(rep("OR<1", nrow(data_matrix1_immune)),
            rep("OR>=1", nrow(data_matrix2_immune))),
  
  explanatory = c(data_matrix1_immune$explanatory,
                  data_matrix2_immune$explanatory),
  
  influence = c(influence_1, influence_2), 
  
  d_stat_percent = influence/T0*100) |>
  
  arrange(desc(influence))


proteomic_influence <- 
  jackknife_immune_results |>                              # selection of the prot with highest influence to the ff test
  filter(d_stat_percent>2) |>   
  pull(explanatory) 
proteomic_influence <- sort(proteomic_influence)

proteomic_selected <-                                                           # selection des prot avec p<0.05
  main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         protein_group == "Immune response",
         p_value_raw<0.05) |> 
  pull(explanatory) 
proteomic_selected <- sort(proteomic_selected)

valeurs_communes <- intersect(proteomic_influence, proteomic_selected)          # check the commune proteins 


### plots ----

plot1_immune <- 
  data_matrix1_immune |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR<1", 
       y = "-log10(p-value)", 
       x = "abslog2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

plot2_immune <- 
  data_matrix2_immune |> 
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR>=1", 
       y = "-log10(p-value)", 
       x = "log2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

additional_analysis_2_immune_figure <- plot1_immune + plot2_immune


rm(data_matrix1_immune, data_matrix2_immune, 
   plot1_immune, plot2_immune, 
   T0, influence_1, influence_2, i, j, mat1_minus, mat2_minus, Ti, Tj, 
   proteomic_influence, proteomic_selected, valeurs_communes)



## Metabolism ----
### data prep ----
data_matrix1_metabolism <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Metabolism") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw < 1) |>
  mutate(
    OR_log2 = abs(log2(OR_raw)),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)

data_matrix2_metabolism <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Metabolism") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw >= 1) |>
  mutate(
    OR_log2 = log2(OR_raw),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)


### Fasano Franceschini test ----
additional_analysis_2_metabolism_results <- 
  fasano.franceschini.test(data_matrix1_metabolism[, c("OR_log2", "p_value_log")],
                           data_matrix2_metabolism[, c("OR_log2", "p_value_log")],
                           seed = 0)

### Jackknife approach ----
T0 <- additional_analysis_2_metabolism_results$statistic                        # d statistic from the inital ff test with all observations 

influence_1 <- numeric(nrow(data_matrix1_metabolism))                           # Initialize vectors
influence_2 <- numeric(nrow(data_matrix2_metabolism))

for (i in seq_len(nrow(data_matrix1_metabolism))) {                             # Jackknife for data_matrix1 OR_raw < 1
  mat1_minus <- data_matrix1_metabolism[-i, ]
  Ti <- ff_stat(mat1_minus, data_matrix2_metabolism)
  influence_1[i] <- T0 - Ti
}

for (j in seq_len(nrow(data_matrix2_metabolism))) {                             # Jackknife for data_matrix2 OR_raw >= 1
  mat2_minus <- data_matrix2_metabolism[-j, ]
  Tj <- ff_stat(data_matrix1_metabolism, mat2_minus)
  influence_2[j] <- T0 - Tj
}

jackknife_metabolism_results <- tibble(                                         # Assemble results 
  group = c(rep("OR<1", nrow(data_matrix1_metabolism)),
            rep("OR>=1", nrow(data_matrix2_metabolism))),
  
  explanatory = c(data_matrix1_metabolism$explanatory,
                  data_matrix2_metabolism$explanatory),
  
  influence = c(influence_1, influence_2), 
  
  d_stat_percent = influence/T0*100) |>
  
  arrange(desc(influence))


proteomic_influence <- 
  jackknife_metabolism_results |>                          # selection of the prot with highest influence to the ff test
  filter(d_stat_percent>2) |>   
  pull(explanatory) 
proteomic_influence <- sort(proteomic_influence)

proteomic_selected <-                                                           # selection des prot avec p<0.05
  main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         protein_group == "Metabolism",
         p_value_raw<0.05) |> 
  pull(explanatory) 
proteomic_selected <- sort(proteomic_selected)

valeurs_communes <- intersect(proteomic_influence, proteomic_selected)          # check the commune proteins 


### plots ----

plot1_metabolism <- 
  data_matrix1_metabolism |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR<1", 
       y = "-log10(p-value)", 
       x = "abslog2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

plot2_metabolism <- 
  data_matrix2_metabolism |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR>=1", 
       y = "-log10(p-value)", 
       x = "log2(OR)") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

additional_analysis_2_metabolism_figure <- plot1_metabolism + plot2_metabolism

rm(data_matrix1_metabolism, data_matrix2_metabolism, 
   plot1_metabolism, plot2_metabolism, 
   T0, influence_1, influence_2, i, j, mat1_minus, mat2_minus, Ti, Tj, 
   proteomic_influence, proteomic_selected, valeurs_communes)

## Neuro-exploratory ----
### data prep ----
data_matrix1_neuro <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Neuro-exploratory") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw < 1) |>
  mutate(
    OR_log2 = abs(log2(OR_raw)),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)

data_matrix2_neuro <- 
  main_results |> 
  filter(analysis == "main",
         model == "adjusted",
         term == "Continuous", 
         protein_group == "Neuro-exploratory") |>
  select(explanatory, OR_raw, p_value_raw) |>
  filter(OR_raw >= 1) |>
  mutate(
    OR_log2 = log2(OR_raw),
    p_value_log = -log10(p_value_raw)) |>
  select(explanatory, OR_log2, p_value_log)

### Fasano Franceschini test ----
additional_analysis_2_neuro_results <- 
  fasano.franceschini.test(data_matrix1_neuro[, c("OR_log2", "p_value_log")],
                           data_matrix2_neuro[, c("OR_log2", "p_value_log")],
                           seed = 0)

### Jackknife approach ----
T0 <- additional_analysis_2_neuro_results$statistic                             # d statistic from the inital ff test with all observations 

influence_1 <- numeric(nrow(data_matrix1_neuro))                                # Initialize vectors
influence_2 <- numeric(nrow(data_matrix2_neuro))

for (i in seq_len(nrow(data_matrix1_neuro))) {                                  # Jackknife for data_matrix1 OR_raw < 1
  mat1_minus <- data_matrix1_neuro[-i, ]
  Ti <- ff_stat(mat1_minus, data_matrix2_neuro)
  influence_1[i] <- T0 - Ti
}

for (j in seq_len(nrow(data_matrix2_neuro))) {                                  # Jackknife for data_matrix2 OR_raw >= 1
  mat2_minus <- data_matrix2_neuro[-j, ]
  Tj <- ff_stat(data_matrix1_neuro, mat2_minus)
  influence_2[j] <- T0 - Tj
}

jackknife_neuro_results <- tibble(                                              # Assemble results 
  group = c(rep("OR<1", nrow(data_matrix1_neuro)),
            rep("OR>=1", nrow(data_matrix2_neuro))),
  
  explanatory = c(data_matrix1_neuro$explanatory,
                  data_matrix2_neuro$explanatory),
  
  influence = c(influence_1, influence_2), 
  
  d_stat_percent = influence/T0*100) |>
  
  arrange(desc(influence))


proteomic_influence <- 
  jackknife_neuro_results |>                               # selection of the prot with highest influence to the ff test
  filter(d_stat_percent>2) |>   
  pull(explanatory) 
proteomic_influence <- sort(proteomic_influence)

proteomic_selected <-                                                           # selection des prot avec p<0.05
  main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         protein_group == "Neuro-exploratory",
         p_value_raw<0.05) |> 
  pull(explanatory) 
proteomic_selected <- sort(proteomic_selected)

valeurs_communes <- intersect(proteomic_influence, proteomic_selected)          # check the commune proteins 


### plots ----
plot1_neuro <- 
  data_matrix1_neuro |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR<1", 
       y = "-log10(p-value)", 
       x = "abslog2(OR)", 
       color = "Influence") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

plot2_neuro <- 
  data_matrix2_neuro |>
  mutate(influence_flag = ifelse(explanatory %in% proteomic_influence,
                                 "influence > 2%",
                                 "influence < 2%")) |>
  ggplot() +
  aes(x = OR_log2, y = p_value_log, color = influence_flag) +
  geom_point() +
  scale_color_manual(values = c("influence > 2%" = "red",
                                "influence < 2%" = "black")) +
  theme_minimal() +
  labs(title = "Associations with OR>=1", 
       y = "-log10(p-value)", 
       x = "log2(OR)", 
       color = "Influence") + 
  xlim(0, 0.6) +
  ylim(0, 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 

additional_analysis_2_neuro_figure <- plot1_neuro + plot2_neuro


rm(data_matrix1_neuro, data_matrix2_neuro, 
   plot1_neuro, plot2_neuro, 
   T0, influence_1, influence_2, i, j, mat1_minus, mat2_minus, Ti, Tj, 
   proteomic_influence, proteomic_selected, valeurs_communes)

# Additionnal analysis 3 - SL learner ----
# cf specific code 

# Additionnal analysis 4 - tidymodels ----
# cf specific code 


# Analyses special NfL ----
## NfL x sex interaction ----
model_interac <- 
  clogit(als ~  sex*proteomic_neuro_explo_NEFL_sd + 
           baseline_age + smoking_2cat_i + bmi, data = bdd_danish)
summary(model_interac)


model_interac_sensi_2 <- 
  clogit(als ~  sex*proteomic_neuro_explo_NEFL_sd +
           baseline_age + smoking_2cat_i + bmi, data = bdd_danish_sensi_2)
summary(model_interac_sensi_2)


model_interac_sensi_1_3_4 <- 
  clogit(als ~  sex*proteomic_neuro_explo_NEFL_sd + 
           baseline_age + smoking_2cat_i + bmi, data = bdd_danish_sensi_1_3_4)
summary(model_interac_sensi_1_3_4)

model_interac_sensi_1_3_5 <- 
  clogit(als ~  sex*proteomic_neuro_explo_NEFL_sd + 
           baseline_age + smoking_2cat_i + bmi, data = bdd_danish_sensi_1_3_5)
summary(model_interac_sensi_1_3_5)

rm(model_interac, 
   model_interac_sensi_2, 
   model_interac_sensi_1_3_4, 
   model_interac_sensi_1_3_5)


## NfL over time (LOESS) ----
# The locally estimated scatterplot smoothing (LOESS) curve illustrates the changes in the ratio of NfL levels in each case–control pair/triplet over time. Each point corresponds to 1 case–control pair/triplet.

# Calcul du ratio : proteomic_neuro_explo_NEFL_case / moyenne(proteomic_neuro_explo_NEFL_controls)
ratios <- 
  bdd_danish |>
  group_by(match) |>
  summarise(
    ratio_proteomic_neuro_explo_NEFL = proteomic_neuro_explo_NEFL[als == 1] / mean(proteomic_neuro_explo_NEFL[als == 0]),
    follow_up = unique(follow_up[als == 1])) |>
  ungroup() |> 
  mutate(follow_up_neg = -follow_up, 
         follow_up_neg = follow_up_neg/12) 

# Calcul du ratio : proteomic_neuro_explo_NEFL_case / moyenne(proteomic_neuro_explo_NEFL_controls)
ratios_sensi_1 <- 
  bdd_danish |>
  mutate(
    proteomic_neuro_explo_NEFL = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  group_by(match) |>
  summarise(
    ratio_proteomic_neuro_explo_NEFL = proteomic_neuro_explo_NEFL[als == 1] / mean(proteomic_neuro_explo_NEFL[als == 0]),
    follow_up = unique(follow_up[als == 1])) |>
  ungroup() |> 
  mutate(follow_up_neg = -follow_up, 
         follow_up_neg = follow_up_neg/12)


## NfL AUC, sensi, speci ----
### All data, unadjusted ----
roc_NfL_all <- roc(
  response = bdd_danish_sensi_1$als,
  predictor = bdd_danish_sensi_1$proteomic_neuro_explo_NEFL,
  levels = c(0, 1), 
  direction = "<")

auc(roc_NfL_all)
ci.auc(roc_NfL_all)

youden_NfL_all <- coords(             # youden = J=Se+Sp−1 donc rapporte 3 memes valeurs youden. donc il vaut mieux choisir celle avec la meilleur sensibilité (ne pas rater des ALS) 
  roc_NfL_all,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")) |>
  as.data.frame()

youden_best_NfL_all <- 
  youden_NfL_all |>
  slice_max(sensitivity, n = 1)

### All data, adjusted ----
model_adjusted_all <- clogit(
  als ~ proteomic_neuro_explo_NEFL + strata(match) + smoking_2cat_i + bmi,
  data = bdd_danish_sensi_1)

bdd_danish_sensi_1$pred_adjusted_all <- predict(model_adjusted_all, type = "lp")

roc_NfL_all_adjusted <- 
  roc(
    response = bdd_danish_sensi_1$als, 
    predictor = bdd_danish_sensi_1$pred_adjusted_all, 
    direction = "<")

auc(roc_NfL_all_adjusted)
ci.auc(roc_NfL_all_adjusted)

youden_NfL_all_adjusted <- coords(
  roc_NfL_all_adjusted,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

youden_NfL_all_adjusted



### Follow-up < 5 years, unadjusted ----
roc_NfL_sensi_2 <- roc(
  response = bdd_danish_sensi_2$als,
  predictor = bdd_danish_sensi_2$proteomic_neuro_explo_NEFL,
  levels = c(0, 1), 
  direction = "<")

auc(roc_NfL_sensi_2)
ci.auc(roc_NfL_sensi_2)

youden_NfL_sensi_2 <- coords(             # youden = J=Se+Sp−1 donc rapporte 3 memes valeurs youden. donc il vaut mieux choisir celle avec la meilleur sensibilité (ne pas rater des ALS) 
  roc_NfL_sensi_2,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")) |>
  as.data.frame()

youden_best_NfL_sensi_2 <- 
  youden_NfL_sensi_2 |>
  slice_max(sensitivity, n = 1)

### Follow-up < 5 years, adjusted ----
model_adjusted_sensi_2 <- clogit(
  als ~ proteomic_neuro_explo_NEFL + strata(match) + smoking_2cat_i + bmi,
  data = bdd_danish_sensi_2)

bdd_danish_sensi_2$pred_adjusted_sensi_2 <- predict(model_adjusted_sensi_2, type = "lp")

roc_NfL_sensi_2_adjusted <- 
  roc(
    response = bdd_danish_sensi_2$als, 
    predictor = bdd_danish_sensi_2$pred_adjusted_sensi_2, 
    direction = "<")

auc(roc_NfL_sensi_2_adjusted)
ci.auc(roc_NfL_sensi_2_adjusted)

youden_NfL_sensi_2_adjusted <- coords(
  roc_NfL_sensi_2_adjusted,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

youden_NfL_sensi_2_adjusted


### Follow-up > 5 years, unadjusted ----
roc_NfL_sensi_1_3 <- roc(
  response = bdd_danish_sensi_1_3$als,
  predictor = bdd_danish_sensi_1_3$proteomic_neuro_explo_NEFL,
  levels = c(0, 1), 
  direction = "<")

auc(roc_NfL_sensi_1_3)
ci.auc(roc_NfL_sensi_1_3)

youden_NfL_sensi_1_3 <- coords(             # youden = J=Se+Sp−1 donc rapporte 3 memes valeurs youden. donc il vaut mieux choisir celle avec la meilleur sensibilité (ne pas rater des ALS) 
  roc_NfL_sensi_1_3,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")) |>
  as.data.frame()

youden_best_NfL_sensi_1_3 <- 
  youden_NfL_sensi_1_3 |>
  slice_max(sensitivity, n = 1)

### Follow-up > 5 years, adjusted ----
model_adjusted_sensi_1_3 <- clogit(
  als ~ proteomic_neuro_explo_NEFL + strata(match) + smoking_2cat_i + bmi,
  data = bdd_danish_sensi_1_3)

bdd_danish_sensi_1_3$pred_adjusted_sensi_1_3 <- predict(model_adjusted_sensi_1_3, type = "lp")

roc_NfL_sensi_1_3_adjusted <- 
  roc(
    response = bdd_danish_sensi_1_3$als, 
    predictor = bdd_danish_sensi_1_3$pred_adjusted_sensi_1_3, 
    direction = "<")

auc(roc_NfL_sensi_1_3_adjusted)
ci.auc(roc_NfL_sensi_1_3_adjusted)

youden_NfL_sensi_1_3_adjusted <- coords(
  roc_NfL_sensi_1_3_adjusted,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

youden_NfL_sensi_1_3_adjusted


### Follow-up > 5 years and < 14.6 years, unadjusted ----
roc_NfL_sensi_1_3_4 <- roc(
  response = bdd_danish_sensi_1_3_4$als,
  predictor = bdd_danish_sensi_1_3_4$proteomic_neuro_explo_NEFL,
  levels = c(0, 1), 
  direction = "<")

auc(roc_NfL_sensi_1_3_4)
ci.auc(roc_NfL_sensi_1_3_4)

youden_NfL_sensi_1_3_4 <- coords(             # youden = J=Se+Sp−1 donc rapporte 3 memes valeurs youden. donc il vaut mieux choisir celle avec la meilleur sensibilité (ne pas rater des ALS) 
  roc_NfL_sensi_1_3_4,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")) |>
  as.data.frame()

youden_best_NfL_sensi_1_3_4 <- 
  youden_NfL_sensi_1_3_4 |>
  slice_max(sensitivity, n = 1)


### Follow-up > 5 years and < 14.6 years, adjusted ----
model_adjusted_sensi_1_3_4 <- clogit(
  als ~ proteomic_neuro_explo_NEFL + strata(match) + smoking_2cat_i + bmi,
  data = bdd_danish_sensi_1_3_4)

bdd_danish_sensi_1_3_4$pred_adjusted_sensi_1_3_4 <- predict(model_adjusted_sensi_1_3_4, type = "lp")

roc_NfL_sensi_1_3_4_adjusted <- 
  roc(
    response = bdd_danish_sensi_1_3_4$als, 
    predictor = bdd_danish_sensi_1_3_4$pred_adjusted_sensi_1_3_4, 
    direction = "<")

auc(roc_NfL_sensi_1_3_4_adjusted)
ci.auc(roc_NfL_sensi_1_3_4_adjusted)

youden_NfL_sensi_1_3_4_adjusted <- coords(
  roc_NfL_sensi_1_3_4_adjusted,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

youden_NfL_sensi_1_3_4_adjusted


### Follow-up > 14.6 years, unadjusted ----
roc_NfL_sensi_1_3_5 <- roc(
  response = bdd_danish_sensi_1_3_5$als,
  predictor = bdd_danish_sensi_1_3_5$proteomic_neuro_explo_NEFL,
  levels = c(0, 1), 
  direction = "<")

auc(roc_NfL_sensi_1_3_5)
ci.auc(roc_NfL_sensi_1_3_5)

youden_NfL_sensi_1_3_5 <- coords(             # youden = J=Se+Sp−1 donc rapporte 3 memes valeurs youden. donc il vaut mieux choisir celle avec la meilleur sensibilité (ne pas rater des ALS) 
  roc_NfL_sensi_1_3_5,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")) |>
  as.data.frame()

youden_best_NfL_sensi_1_3_5 <- 
  youden_NfL_sensi_1_3_5 |>
  slice_max(sensitivity, n = 1)



### Follow-up > 14.6 years, adjusted ----
model_adjusted_sensi_1_3_5 <- clogit(
  als ~ proteomic_neuro_explo_NEFL + strata(match) + smoking_2cat_i + bmi,
  data = bdd_danish_sensi_1_3_5)

bdd_danish_sensi_1_3_5$pred_adjusted_sensi_1_3_5 <- predict(model_adjusted_sensi_1_3_5, type = "lp")

roc_NfL_sensi_1_3_5_adjusted <- 
  roc(
    response = bdd_danish_sensi_1_3_5$als, 
    predictor = bdd_danish_sensi_1_3_5$pred_adjusted_sensi_1_3_5, 
    direction = "<")

auc(roc_NfL_sensi_1_3_5_adjusted)
ci.auc(roc_NfL_sensi_1_3_5_adjusted)

youden_NfL_sensi_1_3_5_adjusted <- coords(
  roc_NfL_sensi_1_3_5_adjusted,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

youden_NfL_sensi_1_3_5_adjusted



# Figures and Tables ----


## Main ----
### Table covariates - als survival ----
covar

### Table proteomic - als occurence - base and adjusted sd ----
proteomic_sd_ALS_table <-                                                       # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "main") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

# proteomic_sd_ALS_table <-                                                     # select only base results
#   main_results |>
#   filter(model == "base" &                                                    # select only continuous results
#            term == "Continuous" & 
#            p_value_raw < 0.05 & 
#            analysis == "main") |>        
#   select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
#   pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
#   select(protein_group, explanatory, contains("base")) |>
#   rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base") |>
#   flextable() |>
#   add_footer_lines(
#     "1All models are matched on birth year and sex. 
#     2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
#     3CI: Confidence interval.") |>
#   add_header(
#     "explanatory" = "Pre-disease plasma proteins", 
#     "protein_group" = "Protein group", 
#     "OR" = "Base models", "95% CI" = "Base models", "p-value" = "Base models") |>
#   theme_vanilla() |>
#   merge_h(part = "header") |>
#   align(align = "center", part = "all") |>
#   merge_v(j = "explanatory") |>
#   bold(j = "explanatory", part = "body") |>
#   align(j = "explanatory", align = "left", part = "all") |> 
#   merge_at(j = "explanatory", part = "header") |>
#   merge_v(j = "protein_group") |>
#   bold(j = "protein_group", part = "body") |>
#   align(j = "protein_group", align = "left", part = "all") |> 
#   merge_at(j = "protein_group", part = "header") |>
#   flextable::font(fontname = "Calibri", part = "all") |> 
#   fontsize(size = 10, part = "all") |>
#   padding(padding.top = 0, padding.bottom = 0, part = "all")



### Table proteomic - als occurence - base and adjusted quart ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "main") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "main") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
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
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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

# extra_rows <- 
#   main_results |>
#   filter(model == "base" &                                                      # select only the base results                                 
#            term != "Continuous" &                                               # select only quartile results
#            analysis == "main") |>            
#   group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
#   filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
#   distinct(protein_group, explanatory) |> 
#   mutate(
#     quartiles = "Quartile 1",
#     "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '')
# 
# proteomic_quart_ALS_table <- 
#   main_results |>
#   filter(model == "base" &                                                      # select quartile results
#            term != "Continuous" & 
#            analysis == "main") |>            
#   group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
#   filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
#   ungroup() |>
#   select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
#   pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
#   select(protein_group, explanatory, quartiles = term, contains("base")) 
# 
# proteomic_quart_ALS_table <- 
#   proteomic_quart_ALS_table |>
#   mutate_if(is.numeric, as.character) |>
#   bind_rows(extra_rows) |>
#   group_by(explanatory) |>
#   mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
#          p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], '')) |>
#   ungroup() |>
#   rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base") |>
#   arrange(protein_group, explanatory, quartiles) |>
#   flextable() |>
#   add_footer_lines(
#     "1All models are matched on birth year and sex. 
#     2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
#     3CI: Confidence interval.
#     4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
#     5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year.") |>
#   add_header(
#     "explanatory" = "Pre-disease plasma proteins", 
#     "protein_group" = "Protein group", 
#     "quartiles" = "Quartiles",
#     "OR" = "Base models", "95% CI" = "Base models", "p-value" = "Base models", 
#     "Heterogeneity test" = "Base models",  "Trend test" = "Base models") |>
#   theme_vanilla() |>
#   merge_h(part = "header") |>
#   align(align = "center", part = "all") |>
#   merge_v(j = "protein_group") |>
#   bold(j = "protein_group", part = "body") |>
#   align(j = "protein_group", align = "left", part = "all") |> 
#   merge_at(j = "protein_group", part = "header") |>
#   merge_v(j = "explanatory") |>
#   bold(j = "explanatory", part = "body") |>
#   align(j = "explanatory", align = "left", part = "all") |> 
#   merge_at(j = "explanatory", part = "header") |>
#   merge_at(j = "quartiles", part = "header") |>
#   flextable::font(fontname = "Calibri", part = "all") |> 
#   fontsize(size = 10, part = "all") |>
#   padding(padding.top = 0, padding.bottom = 0, part = "all")
# 
# rm(extra_rows)



### Figure proteomic - als occurence - base sd ----
proteomic_sd_ALS_base_figure <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "main") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd ----
proteomic_sd_ALS_adjusted_figure <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "main") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))




## Sensi_1 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_1) ----
proteomic_sd_ALS_table_sensi_1 <-                                                       # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_1") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")





### Table proteomic - als occurence - base and adjusted quart (sensi_1) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1 <- 
  proteomic_quart_ALS_table_sensi_1 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p_value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p_value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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

### Figure proteomic - als occurence - base sd (sensi_1) ----
proteomic_sd_ALS_base_figure_sensi_1 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


### Figure proteomic - als occurence - adjusted sd (sensi_1) ----
proteomic_sd_ALS_adjusted_figure_sensi_1 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



rm(all_vars_labels, 
   bdd_danish_sensi_1, 
   proteomic_sensi_1, 
   proteomic_sd_sensi_1, 
   proteomic_quart_sensi_1, 
   proteomic_quart_med_sensi_1)


## Sensi_2 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_2) ----
proteomic_sd_ALS_table_sensi_2 <-                                               # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_2") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - base and adjusted quart (sensi_2) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_2") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_2 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_2") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_2 <- 
  proteomic_quart_ALS_table_sensi_2 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p_value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p_value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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


### Figure proteomic - als occurence - base sd (sensi_2) ----
proteomic_sd_ALS_base_figure_sensi_2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd (sensi_2) ----
proteomic_sd_ALS_adjusted_figure_sensi_2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


rm(bdd_danish_sensi_2, all_vars_labels)


## Sensi_3 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_3) ----
proteomic_sd_ALS_table_sensi_3 <-                                               # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_3") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - base and adjusted quart (sensi_3) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_3") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_3 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_3") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_3 <- 
  proteomic_quart_ALS_table_sensi_3 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p_value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p_value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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


### Figure proteomic - als occurence - base sd (sensi_3) ----
proteomic_sd_ALS_base_figure_sensi_3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd (sensi_3) ----
proteomic_sd_ALS_adjusted_figure_sensi_3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

rm(bdd_danish_sensi_3)






## Sensi_1_3 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_1_3) ----
proteomic_sd_ALS_table_sensi_1_3 <-                                               # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_1_3") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - base and adjusted quart (sensi_1_3) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
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
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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


### Figure proteomic - als occurence - base sd (sensi_1_3) ----
proteomic_sd_ALS_base_figure_sensi_1_3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd (sensi_1_3) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


rm(bdd_danish_sensi_1_3, 
   proteomic_sensi_1_3,
   proteomic_sd_sensi_1_3,
   proteomic_quart_sensi_1_3, 
   proteomic_quart_med_sensi_1_3)



## Sensi_1_3_4 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_1_3_4) ----
proteomic_sd_ALS_table_sensi_1_3_4 <-                                               # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_1_3_4") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - base and adjusted quart (sensi_1_3_4) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3_4") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3_4 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3_4") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
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
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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


### Figure proteomic - als occurence - base sd (sensi_1_3_4) ----
proteomic_sd_ALS_base_figure_sensi_1_3_4 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_4") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd (sensi_1_3_4) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_4") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

rm(proteomic_sensi_1_3_4,
   proteomic_sd_sensi_1_3_4,
   proteomic_quart_sensi_1_3_4, 
   proteomic_quart_med_sensi_1_3_4, 
   bdd_danish_sensi_1_3_4)


## Sensi_1_3_5 ----
### Table proteomic - als occurence - base and adjusted sd (sensi_1_3_5) ----
proteomic_sd_ALS_table_sensi_1_3_5 <-                                               # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_1_3_5") |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - base and adjusted quart (sensi_1_3_5) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3_5") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_3_5 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1_3_5") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
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
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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


### Figure proteomic - als occurence - base sd (sensi_1_3_5) ----
proteomic_sd_ALS_base_figure_sensi_1_3_5 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_5") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



### Figure proteomic - als occurence - adjusted sd (sensi_1_3_5) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_5") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


rm(proteomic_sensi_1_3_5,
   proteomic_sd_sensi_1_3_5,
   proteomic_quart_sensi_1_3_5, 
   proteomic_quart_med_sensi_1_3_5, 
   bdd_danish_sensi_1_3_5)





## Sensi_6 ----
### Figure proteomic - als occurence - base sd (sensi_6) ----
proteomic_sd_ALS_base_figure_sensi_6_T1 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_6_T1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T1)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_base_figure_sensi_6_T2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_6_T2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T2)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_base_figure_sensi_6_T3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_6_T3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T3)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


### Figure proteomic - als occurence - adjusted sd (sensi_6) ----
proteomic_sd_ALS_adjusted_figure_sensi_6_T1 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_6_T1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T1)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_adjusted_figure_sensi_6_T2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_6_T2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T2)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_adjusted_figure_sensi_6_T3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_6_T3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T3)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))


## Sensi_1_6 ----
### Figure proteomic - als occurence - base sd (sensi_1_6) ----
proteomic_sd_ALS_base_figure_sensi_1_6_T1 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T1)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_base_figure_sensi_1_6_T2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T2)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_base_figure_sensi_1_6_T3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models (filtered to follow-up T3)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))





### Figure proteomic - als occurence - adjusted sd (sensi_1_6) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_6_T1 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T1") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T1)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_adjusted_figure_sensi_1_6_T2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T2") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T2)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))

proteomic_sd_ALS_adjusted_figure_sensi_1_6_T3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_6_T3") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models (filtered to follow-up T3)",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5))



## Sensi_1_7 ----
### Table proteomic - als occurence - adjusted sd (sensi_1_7) ----
proteomic_sd_ALS_table_sensi_1_7 <-                                             
  main_results |>   
  filter(model == "adjusted" &                                                  # select only adjusted results
           term == "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1", "sensi_1_7_female", "sensi_1_7_male")) |>            
  group_by(explanatory) |>                                                      # select explanatory vars significant                
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                     
  ungroup() |>
  select(analysis, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, 
         "OR_sensi_1", "95% CI_sensi_1", "p_value_sensi_1", 
         "OR_sensi_1_7_female", "95% CI_sensi_1_7_female", "p_value_sensi_1_7_female", 
         "OR_sensi_1_7_male", "95% CI_sensi_1_7_male", "p_value_sensi_1_7_male") |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", 
         "OR " = "OR_sensi_1_7_female", "95% CI " = "95% CI_sensi_1_7_female", "p-value " = "p_value_sensi_1_7_female", 
         " OR " = "OR_sensi_1_7_male", " 95% CI " = "95% CI_sensi_1_7_male", " p-value " = "p_value_sensi_1_7_male") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and adjusted for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Females and males (main analyses)", "95% CI" = "Females and males (main analyses)", "p-value" = "Females and males (main analyses)", 
    "OR " = "Females", "95% CI " = "Females", "p-value " = "Females", 
    " OR " = "Males", " 95% CI " = "Males", " p-value " = "Males") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### Table proteomic - als occurence - adjusted quart (sensi_1_7) ----
extra_rows <- 
  main_results |>
  filter(model == "adjusted" &                                                  # select only adjusted results
           term != "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1", "sensi_1_7_female", "sensi_1_7_male")) |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_sensi_1" = '-', "95% CI_sensi_1" = '-', "p_value_sensi_1" = '', "p_value_heterogeneity_sensi_1" = '', "p_value_trend_sensi_1" = '',
    "OR_sensi_1_7_female" = '-', "95% CI_sensi_1_7_female" = '-', "p_value_sensi_1_7_female" = '', "p_value_heterogeneity_sensi_1_7_female" = '', "p_value_trend_sensi_1_7_female" = '',
    "OR_sensi_1_7_male" = '-', "95% CI_sensi_1_7_male" = '-', "p_value_sensi_1_7_male" = '', "p_value_heterogeneity_sensi_1_7_male" = '', "p_value_trend_sensi_1_7_male" = '')

proteomic_quart_ALS_table_sensi_1_7 <- 
  main_results |>
  filter(model == "adjusted" &                                                  # select only adjusted results
           term != "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1", "sensi_1_7_female", "sensi_1_7_male")) |>              
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(analysis, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, matches("_sensi_1$"), matches("_sensi_1_7_female$"), matches("_sensi_1_7_male$")) 

proteomic_quart_ALS_table_sensi_1_7 <- 
  proteomic_quart_ALS_table_sensi_1_7 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_sensi_1 = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_sensi_1[quartiles == 'Quartile 2'], ''), 
         p_value_trend_sensi_1 = ifelse(quartiles == 'Quartile 1', p_value_trend_sensi_1[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_sensi_1_7_female = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_sensi_1_7_female[quartiles == 'Quartile 2'], ''), 
         p_value_trend_sensi_1_7_female = ifelse(quartiles == 'Quartile 1', p_value_trend_sensi_1_7_female[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_sensi_1_7_male = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_sensi_1_7_male[quartiles == 'Quartile 2'], ''), 
         p_value_trend_sensi_1_7_male = ifelse(quartiles == 'Quartile 1', p_value_trend_sensi_1_7_male[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", "Heterogeneity test" = "p_value_heterogeneity_sensi_1",  "Trend test" = "p_value_trend_sensi_1",
         "OR " = "OR_sensi_1_7_female", "95% CI " = "95% CI_sensi_1_7_female", "p-value " = "p_value_sensi_1_7_female", "Heterogeneity test " = "p_value_heterogeneity_sensi_1_7_female",  "Trend test " = "p_value_trend_sensi_1_7_female",
         " OR " = "OR_sensi_1_7_male", " 95% CI " = "95% CI_sensi_1_7_male", " p-value " = "p_value_sensi_1_7_male", " Heterogeneity test " = "p_value_heterogeneity_sensi_1_7_male",  " Trend test " = "p_value_trend_sensi_1_7_male") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and adjusted for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across POP quartiles, matched on birth year and adjusted for smoking and body mass index.
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on birth year and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Females and males (main analyses)", "95% CI" = "Females and males (main analyses)", "p-value" = "Females and males (main analyses)",  "Heterogeneity test" = "Females and males (main analyses)",   "Trend test" = "Females and males (main analyses)", 
    "OR " = "Females", "95% CI " = "Females", "p-value " = "Females",  "Heterogeneity test " = "Females",   "Trend test " = "Females", 
    " OR " = "Males", " 95% CI " = "Males", " p-value " = "Males", " Heterogeneity test " = "Males",   " Trend test " = "Males") |>
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


### Figure proteomic - als occurence - base sd (sensi_1_7) ----
proteomic_sd_ALS_base_figure_sensi_1_7 <- 
  main_results |>
  filter(model == "base" &                                                      # select only adjusted results
           term != "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1_7_female", "sensi_1_7_male")) |>     
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  facet_grid(cols = vars(analysis))



### Figure proteomic - als occurence - adjusted sd (sensi_1_7) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_7 <- 
  main_results |>
  filter(model == "adjusted" &                                                  # select only adjusted results
           term != "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1_7_female", "sensi_1_7_male")) |>     
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p_value_raw),
    significance = case_when(
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      p_value_raw < 0.05 & OR > 1 ~ "OR>1 & p-value<0.05",
      p_value_raw < 0.05 & OR < 1 ~ "OR<1 & p-value<0.05",
      TRUE ~ "p-value>0.05")) |>
  ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
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
    values = c("OR<1 & p-value<0.05" = "blue", 
               "p-value>0.05" = "grey70", 
               "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  facet_grid(cols = vars(analysis))


rm(proteomic_sensi_1_7_female,
   proteomic_sd_sensi_1_7_female,
   proteomic_quart_sensi_1_7_female, 
   proteomic_quart_med_sensi_1_7_female, 
   
   proteomic_sensi_1_7_male,
   proteomic_sd_sensi_1_7_male,
   proteomic_quart_sensi_1_7_male, 
   proteomic_quart_med_sensi_1_7_male, 
   
   bdd_danish_sensi_1_7_female,
   bdd_danish_sensi_1_7_male)

## NfL only results ----
### Main analysis ----
#### Table NfL - als occurrence - base and adjusted sd (sensi_1) ----
NfL_sd_ALS_table_sensi_1 <- 
  main_results |>
  filter(analysis == "sensi_1", 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(analysis = fct_recode(analysis, "Main analysis (n=495)" = "sensi_1")) |> 
  select(analysis, model, explanatory, OR, "95% CI", p_value) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 'OR' = 'OR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'OR ' = 'OR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1Base model is matched for sex and birth year. Adjusted model further accounts for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "analysis" = "Analyses",
    "explanatory" = "Explanatory variable", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
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

#### Table NfL - als occurrence - base and adjusted quart (sensi_1) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1" &
           explanatory == "NEFL") |>           
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

NfL_quart_ALS_table_sensi_1 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1" & 
           explanatory == "NEFL") |>             
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

NfL_quart_ALS_table_sensi_1 <- 
  NfL_quart_ALS_table_sensi_1 |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  group_by(explanatory) |>
  mutate(p_value_heterogeneity_base = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_base[quartiles == 'Quartile 2'], ''), 
         p_value_trend_base = ifelse(quartiles == 'Quartile 1', p_value_trend_base[quartiles == 'Quartile 2'], ''),
         p_value_heterogeneity_adjusted = ifelse(quartiles == 'Quartile 1', p_value_heterogeneity_adjusted[quartiles == 'Quartile 2'], ''), 
         p_value_trend_adjusted = ifelse(quartiles == 'Quartile 1', p_value_trend_adjusted[quartiles == 'Quartile 2'], '')) |>
  ungroup() |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", "Heterogeneity test" = "p_value_heterogeneity_base",  "Trend test" = "p_value_trend_base",
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted", "Heterogeneity test " = "p_value_heterogeneity_adjusted",  "Trend test " = "p_value_trend_adjusted") |>
  arrange(protein_group, explanatory, quartiles) |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and birth year. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and birth year, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and birth year, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model",  "Heterogeneity test" = "Base model",   "Trend test" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model", "Heterogeneity test " = "Adjusted model",   "Trend test " = "Adjusted model") |>
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

#### Figure NfL - als occurrence - base and adjusted sd (sensi_1) ----
NfL_sd_ALS_figure_sensi_1 <- main_results |>
  filter(analysis == "sensi_1", 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted models" = "adjusted",
                            "Base models" = "base"), 
         model = fct_relevel(model, "Base models", "Adjusted models"), 
         analysis = fct_recode(analysis, "Main analysis\n(n=495)" = "sensi_1")) |> 
  ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI, color = signif)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis), cols = dplyr::vars(model), switch = "y") +                         # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs( y = "Odd Ratios (ORs)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()

#### Table NfL - als occurrence - base and adjusted sd (sensi_follow-up) ----
NfL_sd_ALS_table_sensi_follow_up_base_adj <- main_results |>
  filter(analysis %in% c("sensi_1", "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(    analysis = fct_recode(analysis, 
                                   "All cases and\ncontrols (n=495)" = "sensi_1",
                                   "Years to ALS < 5 years\n (n=51)" = "sensi_2",
                                   #"Filtered to\nfollow-up > 5 years\n (n=444)" = "sensi_1_3",
                                   "Years to ALS\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
                                   "Years to ALS > 14.6 years\n (n=219)" = "sensi_1_3_5"), 
             analysis = fct_relevel(analysis, 
                                    "All cases and\ncontrols (n=495)", 
                                    "Years to ALS < 5 years\n (n=51)", 
                                    #"Filtered to\nfollow-up > 5 years\n (n=444)", 
                                    "Years to ALS\nbetween 5 and 14.6 years\n(n=225)",
                                    "Years to ALS > 14.6 years\n (n=219)")) |> 
  select(analysis, model, explanatory, OR, "95% CI", p_value) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 'OR' = 'OR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'OR ' = 'OR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and birth year. Adjusted models further account for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "analysis" = "Analyses",
    "explanatory" = "Explanatory variable", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
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


#### Figure NfL - als occurrence - base and adjusted sd (sensi_follow-up) ----
NfL_sd_ALS_figure_sensi_follow_up_base_adj <- main_results |>
  filter(analysis %in% c("sensi_1", "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted models" = "adjusted",
                            "Base models" = "base"), 
         model = fct_relevel(model, "Base models", "Adjusted models"), 
         analysis = fct_recode(analysis, 
                               "All cases and\ncontrols (n=495)" = "sensi_1",
                               "Years to ALS < 5 years\n (n=51)" = "sensi_2",
                               #"Filtered to\nfollow-up > 5 years\n (n=444)" = "sensi_1_3",
                               "Years to ALS\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
                               "Years to ALS > 14.6 years\n (n=219)" = "sensi_1_3_5"), 
         analysis = fct_relevel(analysis, 
                                "All cases and\ncontrols (n=495)", 
                                "Years to ALS < 5 years\n (n=51)", 
                                #"Filtered to\nfollow-up > 5 years\n (n=444)", 
                                "Years to ALS\nbetween 5 and 14.6 years\n(n=225)",
                                "Years to ALS > 14.6 years\n (n=219)")) |> 
  ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI, color = signif)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis), cols = dplyr::vars(model), switch = "y") +                         # , scales = "free_x"
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs( y = "Odd Ratios (ORs)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()


#### Table NfL - als occurrence - adjusted sd (sensi_sex) ----
NfL_sd_ALS_table_sensi_sex_adj <- 
  main_results |>   
  filter(model == "adjusted" &                                                  # select only adjusted results
           term == "Continuous" &                                               # select only continuous results
           analysis %in% c("sensi_1", "sensi_1_7_female", "sensi_1_7_male") &
           explanatory == "NEFL") |>            
  select(analysis, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "analysis", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, 
         "OR_sensi_1", "95% CI_sensi_1", "p_value_sensi_1", 
         "OR_sensi_1_7_female", "95% CI_sensi_1_7_female", "p_value_sensi_1_7_female", 
         "OR_sensi_1_7_male", "95% CI_sensi_1_7_male", "p_value_sensi_1_7_male") |>
  rename("OR" = "OR_sensi_1", "95% CI" = "95% CI_sensi_1", "p-value" = "p_value_sensi_1", 
         "OR " = "OR_sensi_1_7_female", "95% CI " = "95% CI_sensi_1_7_female", "p-value " = "p_value_sensi_1_7_female", 
         " OR " = "OR_sensi_1_7_male", " 95% CI " = "95% CI_sensi_1_7_male", " p-value " = "p_value_sensi_1_7_male") |>
  flextable() |>
  add_footer_lines(
    "1All models are matched on birth year and adjusted for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease plasma proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Females and males (main analyses)", "95% CI" = "Females and males (main analyses)", "p-value" = "Females and males (main analyses)", 
    "OR " = "Females", "95% CI " = "Females", "p-value " = "Females", 
    " OR " = "Males", " 95% CI " = "Males", " p-value " = "Males") |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  align(align = "center", part = "all") |>
  merge_v(j = "explanatory") |>
  bold(j = "explanatory", part = "body") |>
  align(j = "explanatory", align = "left", part = "all") |> 
  merge_at(j = "explanatory", part = "header") |>
  merge_v(j = "protein_group") |>
  bold(j = "protein_group", part = "body") |>
  align(j = "protein_group", align = "left", part = "all") |> 
  merge_at(j = "protein_group", part = "header") |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")

#### Table NfL - als occurrence - base and adjusted sd (sensi_follow-up and sex, for article) ----
NfL_sd_ALS_table_sensi_follow_up_sex_base_adj <- 
  main_results |>
  filter(analysis %in% c("sensi_1", "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", "sensi_1_7_female", "sensi_1_7_male"), 
         term == "Continuous", 
         explanatory == "NEFL") |>
  mutate(analysis = fct_recode(analysis, 
                               "All cases and\ncontrols (n=495)" = "sensi_1",
                               "Years to ALS < 5 years\n (n=51)" = "sensi_2",
                               #"Filtered to\nfollow-up > 5 years\n (n=444)" = "sensi_1_3",
                               "Years to ALS\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
                               "Years to ALS > 14.6 years\n (n=219)" = "sensi_1_3_5", 
                               "Females (n=192)" = "sensi_1_7_female", 
                               "Males (n=303)" = "sensi_1_7_male"), 
         analysis = fct_relevel(analysis, 
                                "All cases and\ncontrols (n=495)", 
                                "Years to ALS < 5 years\n (n=51)", 
                                #"Filtered to\nfollow-up > 5 years\n (n=444)", 
                                "Years to ALS\nbetween 5 and 14.6 years\n(n=225)",
                                "Years to ALS > 14.6 years\n (n=219)", 
                                "Females (n=192)", 
                                "Males (n=303)")) |> 
  select(analysis, model, explanatory, OR, "95% CI", p_value) |> 
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95% CI`, p_value)) |> 
  select(analysis, explanatory, 'OR' = 'OR_base', '95% CI' = '95% CI_base', 'p-value' = 'p_value_base', 
         'OR ' = 'OR_adjusted', '95% CI ' = '95% CI_adjusted', 'p-value ' = 'p_value_adjusted') |>
  flextable() |>
  add_footer_lines(
    "1All models are matched for sex and birth year. Adjusted models further account for smoking, body mass index and marital status. 
  2Estimated risk of ALS for a one standard deviation increase of pre-disease NEFL (NPX). 
  3CI: Confidence interval.") |>
  add_header(
    "analysis" = "Analyses",
    "explanatory" = "Explanatory variable", 
    "OR" = "Base model", "95% CI" = "Base model", "p-value" = "Base model", 
    "OR " = "Adjusted model", "95% CI " = "Adjusted model", "p-value " = "Adjusted model") |>
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



#### Figure NfL - als occurrence - adjusted sd (sensi_follow-up and sex, pour article) ----
NfL_sd_ALS_figure_sensi_follow_up_sex_adj <- 
  main_results |>
  filter(analysis %in% c("sensi_1", 
                         "sensi_2", 
                         #"sensi_1_3", 
                         "sensi_1_3_4", "sensi_1_3_5", 
                         "sensi_1_7_female", "sensi_1_7_male"), 
         term == "Continuous", 
         explanatory == "NEFL", 
         model == "adjusted") |> 
  mutate(signif = ifelse(p_value_raw<0.05, "p-value<0.05", "p-value≥0.05"), 
         analysis = fct_recode(analysis, 
                               "All cases and\ncontrols (n=495)" = "sensi_1",
                               "Years to ALS < 5 years\n (n=51)" = "sensi_2",
                               #"Filtered to\nfollow-up > 5 years\n (n=444)" = "sensi_1_3",
                               "Years to ALS\nbetween 5 and 14.6 years\n(n=225)" = "sensi_1_3_4",
                               "Years to ALS > 14.6 years\n (n=219)" = "sensi_1_3_5", 
                               "Females (n=192)" = "sensi_1_7_female", 
                               "Males (n=303)" = "sensi_1_7_male"), 
         analysis = fct_relevel(analysis, 
                                "All cases and\ncontrols (n=495)", 
                                "Years to ALS < 5 years\n (n=51)", 
                                #"Filtered to\nfollow-up > 5 years\n (n=444)", 
                                "Years to ALS\nbetween 5 and 14.6 years\n(n=225)",
                                "Years to ALS > 14.6 years\n (n=219)", 
                                "Females (n=192)",
                                "Males (n=303)")) |> 
  ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(analysis),  switch = "y") +   
  labs( y = "Odd Ratios") +
  scale_x_discrete(labels = NULL) +
  theme_lucid() +
  theme(axis.text.x = element_text(size = 14,  color = "black", face = "bold"), 
        axis.title.x = element_text(size = 14,  color = "black", face = "bold"), 
        legend.position = "bottom", 
        strip.text.y.left = element_text(hjust = 0.5, vjust = 0.5, angle = 0,  size = 12, color = "black", face = "bold"), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()



### NfL over time (LOESS) ----
figure_NfL_over_time <- 
  ggplot(ratios, aes(x = follow_up_neg, y = ratio_proteomic_neuro_explo_NEFL)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.75) +
  labs(
    x = "Time to ALS diagnosis (years)",
    y = "Ratio of NfL in case and matched controls") +
  theme_lucid() +
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"), 
        axis.text = element_text(size = 16, color = "black", face = "bold"))
rm(ratios)

figure_NfL_over_time_sensi_1 <- 
  ggplot(ratios_sensi_1, aes(x = follow_up_neg, y = ratio_proteomic_neuro_explo_NEFL)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "loess", se = TRUE, color = "red", span = 0.75) +
  labs(
    x = "Time to ALS diagnosis (years)",
    y = "Ratio of NfL in case and matched controls") +
  theme_lucid() +
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"), 
        axis.text = element_text(size = 16, color = "black", face = "bold"))
rm(ratios_sensi_1)


### NfL AUC, sensi, speci ----
label_main_analysis <- 
  paste0("All cases and controls (n=495)\n", 
         "AUC = ", round(auc(roc_NfL_all), 2), 
         "\nOptimal NfL cut-off: ", round(youden_best_NfL_all["threshold"], 2), 
         "\n(sensitivity: ", round(youden_best_NfL_all["sensitivity"], 2), 
         " and specificity: ", round(youden_best_NfL_all["specificity"], 2), ")\n")

label_sensi_2 <- 
  paste0("Years to ALS < 5 years (n=51)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_2), 2), 
         "\nOptimal NfL cut-off: ", round(youden_best_NfL_sensi_2["threshold"], 2), 
         "\n(sensitivity: ", round(youden_best_NfL_sensi_2["sensitivity"], 2), 
         " and specificity: ", round(youden_best_NfL_sensi_2["specificity"], 2), ")\n")

label_sensi_1_3_4 <- 
  paste0("Years to ALS between 5 and 14.6 years (n=225)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_1_3_4), 2), 
         "\nOptimal NfL cut-off: ", round(youden_best_NfL_sensi_1_3_4["threshold"], 2), 
         "\n(sensitivity: ", round(youden_best_NfL_sensi_1_3_4["sensitivity"], 2), 
         " and specificity: ", round(youden_best_NfL_sensi_1_3_4["specificity"], 2), ")\n")

label_sensi_1_3_5 <- 
  paste0("Years to ALS > 14.6 years (n=219)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_1_3_5), 2), 
         "\nOptimal NfL cut-off: ", round(youden_best_NfL_sensi_1_3_5["threshold"], 2), 
         "\n(sensitivity: ", round(youden_best_NfL_sensi_1_3_5["sensitivity"], 2), 
         " and specificity: ", round(youden_best_NfL_sensi_1_3_5["specificity"], 2), ")\n")

roc_patterns <- c(
  label_main_analysis = "solid",
  label_sensi_2 = "dotdash",
  label_sensi_1_3_4 = "dashed",
  label_sensi_1_3_5 = "dotted")

names(roc_patterns) <- c(
  label_main_analysis,
  label_sensi_2,
  label_sensi_1_3_4,
  label_sensi_1_3_5)

#### unadjusted_pattern ----
AUC_figure_unadjusted_pattern <- ggroc(
  list(
    label_main_analysis = roc_NfL_all,
    label_sensi_2 = roc_NfL_sensi_2,
    label_sensi_1_3_4 = roc_NfL_sensi_1_3_4,
    label_sensi_1_3_5 = roc_NfL_sensi_1_3_5) |>
    setNames(c(
      label_main_analysis,
      label_sensi_2,
      label_sensi_1_3_4,
      label_sensi_1_3_5)),
  legacy.axes = TRUE,
  aes = c("linetype"), 
  linewidth = 0.5) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50") +
  scale_linetype_manual(values = roc_patterns) +
  labs(
    x = "1 − Specificity",
    y = "Sensitivity",
    title = "ROC curves for NfL pre-disease biomarker (unmatched, unadjusted)",
    linetype = "") +
  theme_lucid() +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold"),
        axis.text = element_text(size = 18, color = "black", face = "bold"), 
        legend.text = element_text(size = 16, color = "black"), 
        legend.position = "bottom") 

#### unadjusted_color ----
AUC_figure_unadjusted_color <- ggroc(
  list(
    label_main_analysis = roc_NfL_all,
    label_sensi_2 = roc_NfL_sensi_2,
    label_sensi_1_3_4 = roc_NfL_sensi_1_3_4,
    label_sensi_1_3_5 = roc_NfL_sensi_1_3_5) |>
    setNames(c(
      label_main_analysis,
      label_sensi_2,
      label_sensi_1_3_4,
      label_sensi_1_3_5)),
  legacy.axes = TRUE,
  aes = c("color"), 
  linewidth = 2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50") +
  #scale_linetype_manual(values = roc_patterns) +
  labs(
    x = "1 − Specificity",
    y = "Sensitivity",
    title = "ROC curves for NfL pre-disease biomarker (unmatched, unadjusted)",
    color = "") +
  theme_lucid() +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold"),
        axis.text = element_text(size = 18, color = "black", face = "bold"), 
        legend.text = element_text(size = 16, color = "black"), 
        title = element_text(size = 16, color = "black"), 
        legend.position = "bottom") 

#### adjusted_pattern ----
label_main_analysis <- 
  paste0("All cases and controls (n=495)\n", 
         "AUC = ", round(auc(roc_NfL_all_adjusted), 2), "\n")

label_sensi_2 <- 
  paste0("Years to ALS < 5 years (n=51)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_2_adjusted), 2), "\n")

label_sensi_1_3_4 <- 
  paste0("Years to ALS between 5 and 14.6 years (n=225)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_1_3_4_adjusted), 2), "\n")

label_sensi_1_3_5 <- 
  paste0("Years to ALS > 14.6 years (n=219)\n", 
         "AUC = ", round(auc(roc_NfL_sensi_1_3_5_adjusted), 2), "\n")

roc_patterns <- c(
  label_main_analysis = "solid",
  label_sensi_2 = "dotdash",
  label_sensi_1_3_4 = "dashed",
  label_sensi_1_3_5 = "dotted")

names(roc_patterns) <- c(
  label_main_analysis,
  label_sensi_2,
  label_sensi_1_3_4,
  label_sensi_1_3_5)


AUC_figure_adjusted_pattern <- ggroc(
  list(
    label_main_analysis = roc_NfL_all_adjusted,
    label_sensi_2 = roc_NfL_sensi_2_adjusted,
    label_sensi_1_3_4 = roc_NfL_sensi_1_3_4_adjusted,
    label_sensi_1_3_5 = roc_NfL_sensi_1_3_5_adjusted) |>
    setNames(c(
      label_main_analysis,
      label_sensi_2,
      label_sensi_1_3_4,
      label_sensi_1_3_5)),
  legacy.axes = TRUE,
  aes = c("linetype")) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50") +
  scale_linetype_manual(values = roc_patterns) +
  labs(
    x = "1 − Specificity",
    y = "Sensitivity",
    title = "ROC curves for NfL pre-disease biomarker (matched and adjusted)",
    linetype = "") +
  theme_lucid() +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold"),
        axis.text = element_text(size = 18, color = "black", face = "bold"), 
        legend.text = element_text(size = 16, color = "black"), 
        title = element_text(size = 16, color = "black"), 
        legend.position = "bottom") 

rm(roc_NfL_all, youden_NfL_all, youden_best_NfL_all, 
   model_adjusted_all, roc_NfL_all_adjusted, youden_NfL_all_adjusted, 
   roc_NfL_sensi_2, youden_NfL_sensi_2, youden_best_NfL_sensi_2, 
   model_adjusted_sensi_2, roc_NfL_sensi_2_adjusted, youden_NfL_sensi_2_adjusted, 
   roc_NfL_sensi_1_3, youden_NfL_sensi_1_3, youden_best_NfL_sensi_1_3, 
   model_adjusted_sensi_1_3, roc_NfL_sensi_1_3_adjusted, youden_NfL_sensi_1_3_adjusted, 
   roc_NfL_sensi_1_3_4, youden_NfL_sensi_1_3_4, youden_best_NfL_sensi_1_3_4, 
   model_adjusted_sensi_1_3_4, roc_NfL_sensi_1_3_4_adjusted, youden_NfL_sensi_1_3_4_adjusted, 
   roc_NfL_sensi_1_3_5, youden_NfL_sensi_1_3_5, youden_best_NfL_sensi_1_3_5, 
   model_adjusted_sensi_1_3_5, roc_NfL_sensi_1_3_5_adjusted, youden_NfL_sensi_1_3_5_adjusted, 
   roc_patterns, 
   label_main_analysis, label_sensi_2, label_sensi_1_3_4, label_sensi_1_3_5)


### Removing cases and controls with other neurological diseases ----
bdd_danish_sensi_1 <- 
  bdd_danish  |> 
  filter(!match == 159) |>
  mutate(
    proteomic_neuro_explo_NEFL_sd = scale(proteomic_neuro_explo_NEFL)) 

bdd_danish_sensi_2 <- 
  bdd_danish |>
  group_by(match) |>                                                             
  filter(any(als == 1 & follow_up < 60)) |>                                     # sensi 2 : we remove cases (and their controls) with follow-up > 60 months 
  ungroup()|>
  mutate(proteomic_neuro_explo_NEFL_sd = scale(proteomic_neuro_explo_NEFL)) 

bdd_danish_sensi_1_3_4 <- 
  bdd_danish |>
  filter(!match == 159) |>
  group_by(match) |>                                                              
  filter(any(als == 1 & follow_up > 60)) |>                                     # sensi 3 : we remove cases (and their controls) with follow-up < 60 months
  ungroup() |>
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 4 : we remove 50% of the highest values of follow-up (after removing the follow-up<5years)
  group_by(match) |>
  filter(any(als == 1 & follow_up <= seuil)) |>
  ungroup() |>
  select(-seuil) |>
  mutate(proteomic_neuro_explo_NEFL_sd = scale(proteomic_neuro_explo_NEFL))

bdd_danish_sensi_1_3_5 <- 
  bdd_danish |>
  filter(!match == 159) |>
  group_by(match) |>                                                            # sensi 3 : we remove cases (and their controls) with follow-up < 60 months  
  filter(any(als == 1 & follow_up > 60)) |>
  ungroup()|>
  mutate(seuil = quantile(follow_up, 0.5, na.rm = TRUE)) |>                     # sensi 5 : we remove 50% of the lowest values of follow-up (after removing the follow-up<5years)
  group_by(match) |>
  filter(any(als == 1 & follow_up > seuil)) |>
  ungroup() |>
  select(-seuil) |>
  mutate(proteomic_neuro_explo_NEFL_sd = scale(proteomic_neuro_explo_NEFL))

model_list <- list(
  # --- Sensi 1 ---
  "Conditional and adjusted_sensi1" = 
    clogit(als ~ proteomic_neuro_explo_NEFL_sd + strata(match) + bmi + smoking_2cat_i, data = bdd_danish_sensi_1),
  "Adjusted_sensi1" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = bdd_danish_sensi_1),
  "Adjusted removed individuals with other neuro conditions_sensi1neuro" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = filter(bdd_danish_sensi_1, c_diag_name == "None")),
  
  # --- Sensi 2 ---
  "Conditional and adjusted_sensi2" = 
    clogit(als ~ proteomic_neuro_explo_NEFL_sd + strata(match) + bmi + smoking_2cat_i, data = bdd_danish_sensi_2),
  "Adjusted_sensi2" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = bdd_danish_sensi_2),
  "Adjusted removed individuals with other neuro conditions_sensi2neuro" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = filter(bdd_danish_sensi_2, c_diag_name == "None")),
  
  # --- Sensi 1_3_4 ---
  "Conditional and adjusted_sensi134" = 
    clogit(als ~ proteomic_neuro_explo_NEFL_sd + strata(match) + bmi + smoking_2cat_i, data = bdd_danish_sensi_1_3_4),
  "Adjusted_sensi134" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = bdd_danish_sensi_1_3_4),
  "Adjusted removed individuals with other neuro conditions_sensi134neuro" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = filter(bdd_danish_sensi_1_3_4, c_diag_name == "None")),
  
  # --- Sensi 1_3_5 ---
  "Conditional and adjusted_sensi135" = 
    clogit(als ~ proteomic_neuro_explo_NEFL_sd + strata(match) + bmi + smoking_2cat_i, data = bdd_danish_sensi_1_3_5),
  "Adjusted_sensi135" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = bdd_danish_sensi_1_3_5),
  "Adjusted removed individuals with other neuro conditions_sensi135neuro" = 
    glm(als ~ proteomic_neuro_explo_NEFL_sd + birth_year + sex + bmi + smoking_2cat_i, family = binomial, data = filter(bdd_danish_sensi_1_3_5, c_diag_name == "None")))

format_model <- function(mod, id) {
  tidy(mod) |> 
    filter(str_detect(term, "_sd")) |>
    mutate(
      analysis    = str_split_i(id, "_", 2), 
      model       = str_split_i(id, "_", 1), 
      explanatory = "proteomic_neuro_explo_NEFL",
      term        = "Continuous",
      N           = nobs(mod),
      OR_raw      = exp(estimate),
      lower_raw   = exp(estimate - 1.96 * std.error),
      upper_raw   = exp(estimate + 1.96 * std.error),
      OR          = as.numeric(sprintf("%.1f", OR_raw)),
      lower_CI    = as.numeric(sprintf("%.1f", lower_raw)),
      upper_CI    = as.numeric(sprintf("%.1f", upper_raw)),
      `95% CI`    = paste0(sprintf("%.1f", lower_raw), ", ", sprintf("%.1f", upper_raw)),
      p_value_raw = p.value,
      p_value     = if_else(p.value < 0.01, "<0.01", scales::number(p.value, accuracy = 0.01)))
}

all_models_formatted <- imap_dfr(model_list, format_model) |> 
  select(analysis, model, N, explanatory, term, starts_with("OR"), starts_with("95%"), starts_with("p_value"))

table_without_other_diseases <- all_models_formatted |>
  mutate(
    sensi_type = case_when(
      str_detect(analysis, "sensi1$|sensi1neuro") ~ "All timing period",
      str_detect(analysis, "sensi2") ~ "Years to ALS < 5 years",
      str_detect(analysis, "sensi134") ~ "Years to ALS between 5 and 14.6 years",
      str_detect(analysis, "sensi135") ~ "Years to ALS > 14.6 years"),
    sensi_type = fct_relevel(sensi_type, 
                             "All timing period", "Years to ALS < 5 years", "Years to ALS between 5 and 14.6 years",
                             "Years to ALS > 14.6 years"), 
    model_type = fct_recode(model, 
                            "Conditional" = "Conditional and adjusted",
                            "Adjusted_All" = "Adjusted",
                            "Adjusted_NoNeuro" = "Adjusted removed individuals with other neuro conditions")) |>
  complete(sensi_type, model_type) |>
  select(sensi_type, model_type, N, OR, `95% CI`, p_value) |>
  pivot_wider(
    id_cols = sensi_type,
    names_from = model_type,
    values_from = c(N, OR, `95% CI`, p_value),
    names_glue = "{model_type}_{.value}") |>
  mutate(Conditional_N = as.character(Conditional_N), 
         Conditional_N = fct_recode(Conditional_N, "51" = "17", "222" = "74", "495" = "165")) |>
  select(
    Analysis = sensi_type,
    Conditional_N, Conditional_OR, `Conditional_95% CI`, Conditional_p_value,
    Adjusted_All_N, Adjusted_All_OR, `Adjusted_All_95% CI`, Adjusted_All_p_value,
    Adjusted_NoNeuro_N, Adjusted_NoNeuro_OR, `Adjusted_NoNeuro_95% CI`, Adjusted_NoNeuro_p_value)

header_df <- data.frame(
  col_keys = colnames(table_without_other_diseases),
  top_header = c(
    "Analysis",
    rep("Conditional and adjusted logistic regressions (main analyses)", 4),
    rep("Adjusted logistic regressions", 4),
    rep("Adjusted logistic regressions (after removing individuals with neurological condition diagnoses other than ALS", 4)),
  bottom_header = c(
    "Analysis",
    "N", "OR", "95% CI", "p-value",
    "N", "OR", "95% CI", "p-value",
    "N", "OR", "95% CI", "p-value"),
  stringsAsFactors = FALSE)

table_without_other_diseases <- table_without_other_diseases |>
  flextable() |>
  set_header_df(mapping = header_df, key = "col_keys") |>
  merge_h(part = "header") |>
  merge_v(part = "header", j = "Analysis") |> 
  theme_vanilla() |>
  add_footer_lines(
    "Model 1: Conditional logistic regression matched on birth year and sex, adjusted on BMI and smoking status.
     Model 2: Logistic regression adjusted on birth year, sex, BMI and smoking status.
     Model 3: Logistic regression adjusted on birth year, sex, BMI and smoking status (among individuals with no record of neurological condition diagnosis other than ALS).
     CI: Confidence interval. OR estimated for a one standard deviation increase of pre-disease plasma NfL level.") |>
  bold(j = "Analysis", part = "body") |>
  bold(part = "header") |>
  align(align = "center", part = "all") |>
  align(j = "Analysis", align = "left", part = "all") |> 
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 9.5, part = "all") |>
  padding(padding.top = 5, padding.bottom = 5, part = "all") |>
  autofit()

rm(bdd_danish_sensi_1, bdd_danish_sensi_2, bdd_danish_sensi_1_3_4, bdd_danish_sensi_1_3_5, model_list, all_models_formatted, header_df)

# Assemblage ----
results_proteomic_ALS_occurrence <- 
  list(
    main = list(
      covar = covar, 
      main_results= main_results, 
    
      proteomic_sd_ALS_table = proteomic_sd_ALS_table,
      proteomic_quart_ALS_table = proteomic_quart_ALS_table,
    
      proteomic_sd_ALS_base_figure = proteomic_sd_ALS_base_figure, 
      proteomic_sd_ALS_adjusted_figure = proteomic_sd_ALS_adjusted_figure, 
    
      model1_gam = model1_gam, 
      model2_gam = model2_gam),
    
    sensi_1 = list(
      proteomic_sd_ALS_table_sensi_1 = proteomic_sd_ALS_table_sensi_1, 
      proteomic_quart_ALS_table_sensi_1 = proteomic_quart_ALS_table_sensi_1, 
      
      proteomic_sd_ALS_base_figure_sensi_1 = proteomic_sd_ALS_base_figure_sensi_1, 
      model1_gam_sensi_1 = model1_gam_sensi_1, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_1 = proteomic_sd_ALS_adjusted_figure_sensi_1, 
      model2_gam_sensi_1 = model2_gam_sensi_1), 
    
    sensi_2 = list(
      proteomic_sd_ALS_table_sensi_2 = proteomic_sd_ALS_table_sensi_2, 
      proteomic_quart_ALS_table_sensi_2 = proteomic_quart_ALS_table_sensi_2, 
      
      proteomic_sd_ALS_base_figure_sensi_2 = proteomic_sd_ALS_base_figure_sensi_2, 
      model1_gam_sensi_2 = model1_gam_sensi_2, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_2 = proteomic_sd_ALS_adjusted_figure_sensi_2, 
      model2_gam_sensi_2 = model2_gam_sensi_2), 
    
    sensi_3 = list(
      proteomic_sd_ALS_table_sensi_3 = proteomic_sd_ALS_table_sensi_3, 
      proteomic_quart_ALS_table_sensi_3 = proteomic_quart_ALS_table_sensi_3, 
      
      proteomic_sd_ALS_base_figure_sensi_3 = proteomic_sd_ALS_base_figure_sensi_3, 
      model1_gam_sensi_3 = model1_gam_sensi_3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_3 = proteomic_sd_ALS_adjusted_figure_sensi_3, 
      model2_gam_sensi_3 = model2_gam_sensi_3), 
    
    sensi_1_3 = list(
      proteomic_sd_ALS_table_sensi_1_3 = proteomic_sd_ALS_table_sensi_1_3, 
      proteomic_quart_ALS_table_sensi_1_3 = proteomic_quart_ALS_table_sensi_1_3, 
      
      proteomic_sd_ALS_base_figure_sensi_1_3 = proteomic_sd_ALS_base_figure_sensi_1_3, 
      model1_gam_sensi_1_3 = model1_gam_sensi_1_3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_3 = proteomic_sd_ALS_adjusted_figure_sensi_1_3, 
      model2_gam_sensi_1_3 = model2_gam_sensi_1_3), 
    
    sensi_1_3_4 = list(
      proteomic_sd_ALS_table_sensi_1_3_4 = proteomic_sd_ALS_table_sensi_1_3_4, 
      proteomic_quart_ALS_table_sensi_1_3_4 = proteomic_quart_ALS_table_sensi_1_3_4, 
      
      proteomic_sd_ALS_base_figure_sensi_1_3_4 = proteomic_sd_ALS_base_figure_sensi_1_3_4, 
      model1_gam_sensi_1_3_4 = model1_gam_sensi_1_3_4,
      
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_3_4 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_4, 
      model2_gam_sensi_1_3_4 = model2_gam_sensi_1_3_4), 
    
    sensi_1_3_5 = list(
      proteomic_sd_ALS_table_sensi_1_3_5 = proteomic_sd_ALS_table_sensi_1_3_5, 
      proteomic_quart_ALS_table_sensi_1_3_5 = proteomic_quart_ALS_table_sensi_1_3_5, 
      
      proteomic_sd_ALS_base_figure_sensi_1_3_5 = proteomic_sd_ALS_base_figure_sensi_1_3_5, 
      model1_gam_sensi_1_3_5 = model1_gam_sensi_1_3_5, 
      
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_3_5 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_5, 
      model2_gam_sensi_1_3_5 = model2_gam_sensi_1_3_5), 
  
    
    sensi_6 = list(
      sensi_6_table_follow_up = sensi_6_table_follow_up, 
      sensi_6_densityplot_follow_up = sensi_6_densityplot_follow_up,
      
      proteomic_sd_ALS_base_figure_sensi_6_T1 = proteomic_sd_ALS_base_figure_sensi_6_T1, 
      proteomic_sd_ALS_base_figure_sensi_6_T2 = proteomic_sd_ALS_base_figure_sensi_6_T2, 
      proteomic_sd_ALS_base_figure_sensi_6_T3 = proteomic_sd_ALS_base_figure_sensi_6_T3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_6_T1 = proteomic_sd_ALS_adjusted_figure_sensi_6_T1, 
      proteomic_sd_ALS_adjusted_figure_sensi_6_T2 = proteomic_sd_ALS_adjusted_figure_sensi_6_T2, 
      proteomic_sd_ALS_adjusted_figure_sensi_6_T3 = proteomic_sd_ALS_adjusted_figure_sensi_6_T3),
    
    
    sensi_1_6 = list(
      proteomic_sd_ALS_base_figure_sensi_1_6_T1 = proteomic_sd_ALS_base_figure_sensi_1_6_T1, 
      proteomic_sd_ALS_base_figure_sensi_1_6_T2 = proteomic_sd_ALS_base_figure_sensi_1_6_T2, 
      proteomic_sd_ALS_base_figure_sensi_1_6_T3 = proteomic_sd_ALS_base_figure_sensi_1_6_T3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_6_T1 = proteomic_sd_ALS_adjusted_figure_sensi_1_6_T1, 
      proteomic_sd_ALS_adjusted_figure_sensi_1_6_T2 = proteomic_sd_ALS_adjusted_figure_sensi_1_6_T2, 
      proteomic_sd_ALS_adjusted_figure_sensi_1_6_T3 = proteomic_sd_ALS_adjusted_figure_sensi_1_6_T3), 
    
    sensi_1_7 = list(
      proteomic_sd_ALS_table_sensi_1_7 = proteomic_sd_ALS_table_sensi_1_7, 
      proteomic_quart_ALS_table_sensi_1_7 = proteomic_quart_ALS_table_sensi_1_7, 
      proteomic_sd_ALS_base_figure_sensi_1_7 = proteomic_sd_ALS_base_figure_sensi_1_7, 
      proteomic_sd_ALS_adjusted_figure_sensi_1_7 = proteomic_sd_ALS_adjusted_figure_sensi_1_7), 
    
    additional_analysis_1 = list(
      OR_distribution = OR_distribution, 
      boxplot_OR = boxplot_OR, 
      densityplot_OR = densityplot_OR), 
    
    additional_analysis_2 = list(
      additional_analysis_2_all_figure = additional_analysis_2_all_figure, 
      additional_analysis_2_all_results = additional_analysis_2_all_results, 
      jackknife_all_results = jackknife_all_results, 
      
      additional_analysis_2_immune_figure = additional_analysis_2_immune_figure, 
      additional_analysis_2_immune_results = additional_analysis_2_immune_results, 
      jackknife_immune_results = jackknife_immune_results, 
      
      additional_analysis_2_metabolism_figure = additional_analysis_2_metabolism_figure, 
      additional_analysis_2_metabolism_results = additional_analysis_2_metabolism_results, 
      jackknife_metabolism_results = jackknife_metabolism_results, 
      
      additional_analysis_2_neuro_figure = additional_analysis_2_neuro_figure, 
      additional_analysis_2_neuro_results = additional_analysis_2_neuro_results, 
      jackknife_neuro_results = jackknife_neuro_results), 
    
    #additional_analysis_3 = SL (check specific code)
    
    NfL_results = list(NfL_sd_ALS_table_sensi_1 = NfL_sd_ALS_table_sensi_1, 
                       NfL_quart_ALS_table_sensi_1 = NfL_quart_ALS_table_sensi_1, 
                       NfL_sd_ALS_figure_sensi_1 = NfL_sd_ALS_figure_sensi_1, 
                       NfL_sd_ALS_table_sensi_follow_up_base_adj = NfL_sd_ALS_table_sensi_follow_up_base_adj, 
                       NfL_sd_ALS_table_sensi_follow_up_sex_base_adj = NfL_sd_ALS_table_sensi_follow_up_sex_base_adj, 
                       NfL_sd_ALS_figure_sensi_follow_up_base_adj = NfL_sd_ALS_figure_sensi_follow_up_base_adj, 
                       NfL_sd_ALS_figure_sensi_follow_up_sex_adj = NfL_sd_ALS_figure_sensi_follow_up_sex_adj, 
                       NfL_sd_ALS_table_sensi_sex_adj = NfL_sd_ALS_table_sensi_sex_adj, 
                       NfL_over_time = list(
                         figure_NfL_over_time = figure_NfL_over_time ,
                         figure_NfL_over_time_sensi_1 = figure_NfL_over_time_sensi_1), 
                       AUC_NfL = list(
                         AUC_figure_adjusted_pattern = AUC_figure_adjusted_pattern, 
                         AUC_figure_unadjusted_pattern = AUC_figure_unadjusted_pattern, 
                         AUC_figure_unadjusted_color = AUC_figure_unadjusted_color), 
                       other_disease = list(
                         table_without_other_diseases = table_without_other_diseases)))


saveRDS(results_proteomic_ALS_occurrence, file = "~/Documents/POP_ALS_2025_02_03/2_output/2.6.1_results_proteomic_ALS_occurrence.rds")

rm(covar,                                   # main results 
   main_results, 
   proteomic_sd_ALS_table,
   proteomic_quart_ALS_table, 
   proteomic_sd_ALS_base_figure, 
   proteomic_sd_ALS_adjusted_figure, 
   model1_gam, 
   model2_gam,
   
   proteomic_sd_ALS_table_sensi_1,         # sensi_1
   proteomic_quart_ALS_table_sensi_1, 
   proteomic_sd_ALS_base_figure_sensi_1, 
   model1_gam_sensi_1, 
   proteomic_sd_ALS_adjusted_figure_sensi_1, 
   model2_gam_sensi_1, 
   
   proteomic_quart_ALS_table_sensi_2,       # sensi_2
   proteomic_sd_ALS_table_sensi_2, 
   proteomic_sd_ALS_base_figure_sensi_2, 
   model1_gam_sensi_2, 
   proteomic_sd_ALS_adjusted_figure_sensi_2, 
   model2_gam_sensi_2, 
   
   proteomic_sd_ALS_table_sensi_3,         # sensi_3
   proteomic_quart_ALS_table_sensi_3, 
   proteomic_sd_ALS_base_figure_sensi_3,
   model1_gam_sensi_3, 
   proteomic_sd_ALS_adjusted_figure_sensi_3,
   model2_gam_sensi_3, 
   
   proteomic_sd_ALS_table_sensi_1_3,          # sensi_1_3
   proteomic_quart_ALS_table_sensi_1_3, 
   proteomic_sd_ALS_base_figure_sensi_1_3,
   model1_gam_sensi_1_3, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3,
   model2_gam_sensi_1_3, 
   
   proteomic_sd_ALS_table_sensi_1_3_4,       # sensi_1_3_4
   proteomic_quart_ALS_table_sensi_1_3_4, 
   proteomic_sd_ALS_base_figure_sensi_1_3_4,
   model1_gam_sensi_1_3_4, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_4,
   model2_gam_sensi_1_3_4, 
   
   proteomic_sd_ALS_table_sensi_1_3_5,        # sensi_1_3_5
   proteomic_quart_ALS_table_sensi_1_3_5, 
   proteomic_sd_ALS_base_figure_sensi_1_3_5,
   model1_gam_sensi_1_3_5, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_5,
   model2_gam_sensi_1_3_5, 
   
   sensi_6_table_follow_up,                       # sensi_6  
   sensi_6_densityplot_follow_up, 
   proteomic_sd_ALS_base_figure_sensi_6_T1, 
   proteomic_sd_ALS_base_figure_sensi_6_T2, 
   proteomic_sd_ALS_base_figure_sensi_6_T3, 
   proteomic_sd_ALS_adjusted_figure_sensi_6_T1, 
   proteomic_sd_ALS_adjusted_figure_sensi_6_T2, 
   proteomic_sd_ALS_adjusted_figure_sensi_6_T3, 
   
   proteomic_sd_ALS_base_figure_sensi_1_6_T1,      # sensi_1_6
   proteomic_sd_ALS_base_figure_sensi_1_6_T2, 
   proteomic_sd_ALS_base_figure_sensi_1_6_T3, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_6_T1, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_6_T2, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_6_T3, 
   
   proteomic_sd_ALS_table_sensi_1_7,               # sensi_1_7
   proteomic_quart_ALS_table_sensi_1_7, 
   proteomic_sd_ALS_base_figure_sensi_1_7, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_7, 
   model1_gam_sensi_1_7_female, 
   model1_gam_sensi_1_7_male, 
   model2_gam_sensi_1_7_female, 
   model2_gam_sensi_1_7_male,
   
   OR_distribution,                        # additionnal analysis 1
   boxplot_OR, 
   densityplot_OR, 
   
   additional_analysis_2_all_figure,      # additionnal analysis 2
   additional_analysis_2_all_results, 
   jackknife_all_results, 
   
   additional_analysis_2_immune_figure, 
   additional_analysis_2_immune_results, 
   jackknife_immune_results, 
   
   additional_analysis_2_metabolism_figure, 
   additional_analysis_2_metabolism_results, 
   jackknife_metabolism_results, 
   
   additional_analysis_2_neuro_figure, 
   additional_analysis_2_neuro_results, 
   jackknife_neuro_results, 
   
   # additionnal analysis 3 : check specfic SL code 
   
   NfL_sd_ALS_table_sensi_1,                      # NfL main results
   NfL_quart_ALS_table_sensi_1, 
   NfL_sd_ALS_figure_sensi_1, 
   NfL_sd_ALS_table_sensi_follow_up_base_adj, 
   NfL_sd_ALS_table_sensi_follow_up_sex_base_adj, 
   NfL_sd_ALS_figure_sensi_follow_up_base_adj, 
   NfL_sd_ALS_figure_sensi_follow_up_sex_adj, 
   NfL_sd_ALS_table_sensi_sex_adj,
   
   figure_NfL_over_time,                      # NfL results over time
   figure_NfL_over_time_sensi_1, 
   
   AUC_figure_unadjusted_pattern,          # NfL AUC
   AUC_figure_unadjusted_color, 
   AUC_figure_adjusted_pattern, 
   
   table_without_other_diseases          # NfL other diseases
   
)
  

# Aline Davias
# October 20, 2025 
# Analysis of als risk ALS depending on proteomic profile


# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.5_analyses_fattyacids_ALS_survival.R", echo=TRUE)
covariates <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'fS_Kol', 'marital_status_2cat_i', 'education_i')


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
#### sd ----
model1_sd <- data.frame(variable = character(),
                           df = integer(),
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
  df_value <- model_summary$term
  
  model1_sd <- rbind(model1_sd, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model1_sd <- model1_sd |> 
  mutate(
    df = case_when(
      grepl("_sd", df) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base_sd") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


#### quartiles ----
model1_quart <- data.frame(variable = character(),
                           df = integer(),
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

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### model 2 ----
# matched on sex and age, adjusted on for smoking_2cat_i, BMI, serum total fS_Kol, marital status and education

#### sd ----
model2_sd <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_sd <- rbind(model2_sd, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model2_sd <- model2_sd |> 
  mutate(
    df = case_when(
      grepl("_sd", df) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted_sd") |>
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


#### quartiles ----
model2_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
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

for (var in proteomic) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


### heterogeneity tests ----
#### model 1 quartile ----
heterogeneity_base_quart <- data.frame(variable = character(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
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

for (var in proteomic_quart) {
  
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


heterogeneity_tests <- 
  bind_rows(heterogeneity_base_quart, 
            heterogeneity_adjusted_quart) |>
  mutate(variable = gsub("_quart", "", variable))

### trend tests ----
#### model 1 quartile ----
trend_base <- data.frame(variable = character(),
                         model = factor(), 
                         p.value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
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

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  test <- clogit(formula, data = bdd_danish) |> summary()
  p.value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(variable = var,
                                     model = "adjusted_quart",
                                     p.value_trend = p.value_trend))
}
rm(var, test, formula, p.value_trend)


trend_tests <- 
  bind_rows(trend_base, trend_adjusted) |>
  mutate(variable = gsub("_quart_med", "", variable))



### merging the main results ----
main_results <- bind_rows(model1_quart, 
                          model2_quart, 
                          model1_sd, 
                          model2_sd) |> 
  mutate(variable = gsub("_quart", "", variable), 
         variable = gsub("_sd", "", variable), 
         OR_raw = OR, 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", number(p.value, accuracy = 0.01, decimal.mark = ".")), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  group_by(model) %>%                                  # correction séparée par modèle
  mutate(
    p.value_fdr = if_else(
      df == "Continuous",                             # seulement pour sd == "Continuous"
      p.adjust(p.value_raw, method = "fdr"),
      NA_real_                                        # NA sinon
    )
  ) |>
  ungroup() |>
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


rm(model1_quart, model2_quart, 
   model1_sd, model2_sd,
   heterogeneity_base_quart, heterogeneity_adjusted_quart, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests)



## Figures and Tables ----
main_results |>
  filter(model == "base_sd") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p.value_raw),
    significance = case_when(
      p.value_raw < 0.05 & OR > 1 ~ "Up",
      p.value_raw < 0.05 & OR < 1 ~ "Down",
      TRUE ~ "NS"
    )
  ) |>
  ggplot(aes(x = log2OR, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot des résultats logistiques",
    x = "log2(OR)",
    y = "-log10(p-value)",
    color = "Effet"
  )


main_results |>
  filter(model == "base_sd") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p.value_raw),
    significance = case_when(
      p.value_fdr < 0.05 & OR_raw > 1 ~ "Up_FDR",
      p.value_fdr < 0.05 & OR_raw < 1 ~ "Down_FDR",
      p.value_raw < 0.05 & OR_raw > 1 ~ "Up_raw",
      p.value_raw < 0.05 & OR_raw < 1 ~ "Down_raw",
      TRUE ~ "NS"
    )
  ) |>
  ggplot(aes(x = log2OR, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(
    values = c(
      "Down_FDR" = "blue",
      "Up_FDR" = "red",
      "Down_raw" = "lightblue",
      "Up_raw" = "lightcoral",
      "NS" = "grey70"
    ),
    labels = c(
      "Down_FDR" = "OR < 1 (FDR < 0.05)",
      "Up_FDR" = "OR > 1 (FDR < 0.05)",
      "Down_raw" = "OR < 1 (p < 0.05)",
      "Up_raw" = "OR > 1 (p < 0.05)",
      "NS" = "Non significatif"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot des résultats logistiques",
    x = "log2(OR)",
    y = "-log10(p-value)",
    color = "Effet"
  )


main_results |>
  filter(model == "adjusted_sd") |>
  mutate(
    log2OR = log2(OR_raw),
    neg_log10_p = -log10(p.value_raw),
    significance = case_when(
      p.value_fdr < 0.05 & OR_raw > 1 ~ "Up_FDR",
      p.value_fdr < 0.05 & OR_raw < 1 ~ "Down_FDR",
      p.value_raw < 0.05 & OR_raw > 1 ~ "Up_raw",
      p.value_raw < 0.05 & OR_raw < 1 ~ "Down_raw",
      TRUE ~ "NS"
    )
  ) |>
  ggplot(aes(x = log2OR, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(
    values = c(
      "Down_FDR" = "blue",
      "Up_FDR" = "red",
      "Down_raw" = "lightblue",
      "Up_raw" = "lightcoral",
      "NS" = "grey70"
    ),
    labels = c(
      "Down_FDR" = "OR < 1 (FDR < 0.05)",
      "Up_FDR" = "OR > 1 (FDR < 0.05)",
      "Down_raw" = "OR < 1 (p < 0.05)",
      "Up_raw" = "OR > 1 (p < 0.05)",
      "NS" = "Non significatif"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot des résultats logistiques",
    x = "log2(OR)",
    y = "-log10(p-value)",
    color = "Effet"
  )



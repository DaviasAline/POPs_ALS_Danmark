# Aline Davias
# October 20, 2025 
# Analysis of als risk depending on proteomic profile


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
    model = "base") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### quartiles ----
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
    model = "base") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


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
model2_sd <- data.frame(explanatory = character(),
                           term = integer(),
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
    model = "adjusted") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

#### quartiles ----
model2_quart <- data.frame(explanatory = character(),
                           term = integer(),
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
    model = "adjusted") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

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
heterogeneity_base_quart <- data.frame(explanatory = character(),
                                       model = factor(),
                                       p_value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart <- rbind(heterogeneity_base_quart, 
                                    data.frame(explanatory = var,
                                               model = "base",
                                               p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)

#### model 2 quartile ----
heterogeneity_adjusted_quart <- data.frame(explanatory = character(),
                                           model = factor(), 
                                           p_value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in proteomic_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(explanatory = var,
                                                   model = "adjusted",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_quart, 
            heterogeneity_adjusted_quart)  |>
  mutate(explanatory = gsub("_quart", "", explanatory))

### trend tests ----
#### model 1 quartile ----
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
                                 p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)

#### model 2 quartile ----
trend_adjusted <- data.frame(explanatory = character(),
                             model = factor(), 
                             p_value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in proteomic_quart_med) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + fS_Kol + marital_status_2cat_i + education_i"))
  test <- clogit(formula, data = bdd_danish) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(explanatory = var,
                                     model = "adjusted",
                                     p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))



### merging the main results ----
main_results <- bind_rows(model1_quart, 
                          model2_quart, 
                          model1_sd, 
                          model2_sd) |> 
  mutate(explanatory = gsub("_quart", "", explanatory), 
         explanatory = gsub("_sd", "", explanatory), 
         OR_raw = OR, 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         p_value_raw = p_value, 
         p_value = ifelse(p_value < 0.01, "<0.01", number(p_value, accuracy = 0.01, decimal.mark = ".")), 
         "95% CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  group_by(model) %>%                                  # correction séparée par modèle
  mutate(
    p_value_fdr = if_else(
      term == "Continuous",                             # seulement pour sd == "Continuous"
      p.adjust(p_value_raw, method = "fdr"),
      NA_real_                                        # NA sinon
    )
  ) |>
  ungroup() |>
  arrange(explanatory) |>
  select(explanatory, 
         model,
         term,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p_value"), 
         lower_CI, upper_CI) 

main_results <- left_join(main_results, heterogeneity_tests, by = c("explanatory", "model"))
main_results <- left_join(main_results, trend_tests, by = c("explanatory", "model"))
main_results <- main_results |>
  mutate(
    p_value_heterogeneity = ifelse(p_value_heterogeneity < 0.01, "<0.01", number(p_value_heterogeneity, accuracy = 0.01, decimal.mark = ".")), 
    p_value_trend = ifelse(p_value_trend < 0.01, "<0.01", number(p_value_trend, accuracy = 0.01, decimal.mark = ".")), 
    protein_group = case_when(str_detect(explanatory, 'proteomic_immun_res') ~ "Immune response", 
                              str_detect(explanatory, 'proteomic_metabolism') ~ "Metabolism", 
                              str_detect(explanatory, 'proteomic_neuro_explo') ~ "Neuro-exploratory"), 
    explanatory = str_replace(explanatory, 'proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_', ""))


rm(model1_quart, model2_quart, 
   model1_sd, model2_sd,
   heterogeneity_base_quart, heterogeneity_adjusted_quart, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests)



# Figures and Tables ----
### Table covariates - als survival ----
covar

## Table proteomic - als occurence - base and adjusted sd ----
proteomic_sd_ALS_table <- main_results |>
  filter(model %in% c("base", "adjusted") & term == "Continuous") |>
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, explanatory, protein_group, OR, "95% CI", "p_value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value")) |>
  select(protein_group, explanatory, contains("base"), contains("adjusted")) |>
  rename("OR" = "OR_base", "95% CI" = "95% CI_base", "p-value" = "p_value_base", 
         "OR " = "OR_adjusted", "95% CI " = "95% CI_adjusted", "p-value " = "p_value_adjusted") |>
  flextable() |>
  add_footer_lines(
    "1All models are adjusted for age at baseline and sex. Adjusted models further account for smoking, BMI, cholesterol, marital status and education. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
    "protein_group" = "Protein group", 
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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


## Table proteomic - als occurence - base and adjusted quart ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous") |>              # select quartile results
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table <- 
  main_results |>
  filter(model %in% c("base", "adjusted") & term != "Continuous") |>              # select quartile results
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
    "1All models are adjusted for age at baseline and sex. Adjusted models further account for smoking, BMI, cholesterol, marital status and education. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, adjusted for sex and age at baseline. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, adjusted for sex and age at baseline. 
    6Heterogeneity tests in outcome value across POP quartiles, adjusted for sex, age at baseline, smoking, BMI, cholesterol, marital status and education.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, adjusted for sex, age at baseline, smoking, BMI, cholesterol, marital status and education.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
    "protein_group" = "Protein group", 
    "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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



## Figure proteomic - als occurence - base sd ----
proteomic_sd_ALS_base_figure <- 
  main_results |>
  filter(model == "base" & term == "Continuous") |>
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
    color = "")

## Figure proteomic - als occurence - adjusted sd ----
proteomic_sd_ALS_adjusted_figure <- 
  main_results |>
  filter(model == "adjusted" & term == "Continuous") |>
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
    values = c("OR<1 & p-value<0.05" = "blue", "p-value>0.05" = "grey70", "OR>1 & p-value<0.05" = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Adjusted logistic models",
    x = "OR",
    y = "-log10(p-value)", 
    color = "")

## Figure proteomic - als occurrence - base gam ----
pvals <- sapply(model1_gam, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]

# pollutant_labels <- set_names(
#   c("Most prevalent PCBs", "Dioxin-like PCBs","Non-dioxin-like PCBs", "HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
#   POPs_group)
plot_base_gam <- map(signif_vars, function(var) {
  
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
  # x_label <- pollutant_labels[var] 
  
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
    # ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |>
  set_names(signif_vars)
# rm(pollutant_labels)

rm(pvals, signif_vars)

## Figure proteomic - als occurrence adjusted gam ----
pvals <- sapply(model2_gam, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]

plot_adjusted_gam <- map(signif_vars, function(var) {
  
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
  # x_label <- pollutant_labels[var] 
  
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
    # ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |> 
  set_names(signif_vars)
# rm(pollutant_labels)

rm(pvals, signif_vars)


# Assemblage ----
results_proteomic_ALS_occurrence <- 
  list(
    covar = covar, 
    main_results= main_results, 
    
    proteomic_sd_ALS_table = proteomic_sd_ALS_table,
    proteomic_quart_ALS_table = proteomic_quart_ALS_table,
    
    proteomic_sd_ALS_base_figure = proteomic_sd_ALS_base_figure, 
    proteomic_sd_ALS_adjusted_figure = proteomic_sd_ALS_adjusted_figure, 
    
    model1_gam = model1_gam, 
    model2_gam = model2_gam,
    plot_base_gam = plot_base_gam, 
    plot_adjusted_gam = plot_adjusted_gam)

rm(covariates, 
   covar, 
   main_results, 
   proteomic_sd_ALS_table,
   proteomic_quart_ALS_table, 
   proteomic_sd_ALS_base_figure, 
   proteomic_sd_ALS_adjusted_figure, 
   model1_gam, 
   model2_gam,
   plot_base_gam, 
   plot_adjusted_gam)
  

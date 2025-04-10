# Aline Davias
# April 9, 2025 
# Analyses on fatty acids and PUFAs


## data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2_analyses_POPs_ALS.R")

## main analysis ----
model1 <- map_dfr(explanatory, function(var) {                                   # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  
  tibble(
    model = "model 1",
    variable = var,
    df = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

model2 <- map_dfr(explanatory, function(var) {                                   # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  
  tibble(
    model = "model 2",
    variable = var,
    df = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

main_results_fattyacids_ALS <- bind_rows(model1, model2) |>
  filter(variable == df) |>
  mutate(
    OR = as.numeric(sprintf("%.1f", OR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    "95%CI" = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw`= `p-value`, 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = "."))) 

rm(model1, model2)

## results presentation ----
### table 1 ----
table_1 <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1)) |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  


### table 2 ----
table_2 <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
      select(als, all_of(covariates_danish)) |>
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
  tab_spanner = c("**Univariate**", "**Adjusted**")) |> 
  as_flex_table()|>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")

### table 3 ----
table_3 <- main_results_fattyacids_ALS |>
  select(model, variable, OR, "95%CI", "p-value") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95%CI", "p-value")) |>
  select(variable, contains("model 1"), contains("model 2")) |>
  rename("OR" = "OR_model 1", "95% CI" = "95%CI_model 1", "p-value" = "p-value_model 1", 
         "OR " = "OR_model 2", "95% CI " = "95%CI_model 2", "p-value " = "p-value_model 2") |>
  mutate(variable = fct_recode(variable, 
    "Adrenic acid (ω6)" = "adrenic_acid_ω6",
    "Arachidonic acid (ω6)" = "arachidonic_acid_ω6",
    "Cervonic acid (DHA, ω3)" = "cervonic_acid_ω3",
    "Clupanodonic acid (DPA, ω3)" = "clupanodonic_acid_ω3",
    "Dihomo-γ-linolenic acid (ω6)" = "dihomo_γ_linolenic_acid_ω6",
    "Linoleic acid (LA, ω6)" = "linoleic_acid_ω6",
    "Total PUFAs" = "pufas",
    "Total ω3 PUFAs" = "pufas_ω3",
    "Total ω6 PUFAs" = "pufas_ω6",
    "Total ω9 PUFAs" = "pufas_ω9",
    "Rumenic acid (ω6)" = "rumenic_acid_ω6",
    "Timnodonic acid (EPA, ω3)" = "timnodonic_acid_ω3",
    "α-linolenic acid (ALA, ω3)" = "α_linolenic_acid_ω3"))

table_3 <- table_3 |> flextable() |>
  add_footer_lines(
    "1All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. 
  2Estimated risk of ALS when pre-disease serum concentration of PUFAs is increasing by one unit (pg/ml).
  3CI: Confidence interval.") |>
  add_header(
    "variable" = "Exposures", 
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
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



### figure 1 ----
figure_1 <- bdd_danish |>
  select(als, all_of(explanatory)) |>
  pivot_longer(cols = -als, names_to = "PUFAs", values_to = "Values") |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1"), 
    PUFAs =   fct_relevel(PUFAs, 
                          "pufas", "pufas_ω9", "pufas_ω6", "pufas_ω3", "linoleic_acid_ω6",
                          "dihomo_γ_linolenic_acid_ω6", "arachidonic_acid_ω6", "adrenic_acid_ω6",
                          "rumenic_acid_ω6", "α_linolenic_acid_ω3", "timnodonic_acid_ω3",
                          "clupanodonic_acid_ω3", "cervonic_acid_ω3"), 
    PUFAs = fct_rev(PUFAs)) |> 
  arrange(PUFAs) |>
  ggplot() +
  aes(x = Values, y = PUFAs, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_x_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(fill = "ALS ", x = "Values (pg/ml, log transformed)", y = "PUFAs") +
  theme_lucid() + 
  theme(legend.position = "bottom")


### figure 2 ----
figure_2 <- main_results_fattyacids_ALS %>% 
  filter(!variable == "adrenic_acid_ω6") |>
  mutate(p.value_shape = ifelse('p-value_raw'<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Base model" = "model 1",
                            "Adjusted model" = "model 2"),
         model = fct_relevel(model, 'Base model', 'Adjusted model')) %>%
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "PUFAs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()


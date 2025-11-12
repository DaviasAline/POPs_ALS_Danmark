# Aline Davias
# October 20, 2025 
# Analysis of als risk depending on proteomic profile


# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.5_analyses_fattyacids_ALS_survival.R", echo=TRUE)


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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gam[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


# Sensitivity 1 - Removing the oulier for NEFL ----
bdd_danish_sensi_1 <- 
  bdd_danish  |> 
  mutate(
    proteomic_neuro_explo_NEFL_sensi_1 = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd_sensi_1")) |>
  mutate(across(all_of(proteomic), ~ factor(ntile(.x, 4),                           
                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart_sensi_1")) |>
  mutate(across(
    all_of(proteomic),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med_sensi_1"
  ))

proteomic_sensi_1 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1")
proteomic_sd_sensi_1 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sd_sensi_1")
proteomic_quart_sensi_1 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_quart_sensi_1")
proteomic_quart_med_sensi_1 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_quart_med_sensi_1")


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
heterogeneity_base_sensi_1 <- data.frame(explanatory = character(),
                                        model = factor(),
                                        p_value_heterogeneity = numeric(), 
                                        stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish_sensi_1)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_1)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_sensi_1 <- rbind(heterogeneity_base_sensi_1, 
                                     data.frame(explanatory = var,
                                                model = "base", 
                                                analysis = "sensi_1",
                                                p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1)
  model_summary <- summary(model)
  model1_gam_sensi_1[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model1_gam_sensi_1, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_base_gam_sensi_1 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_1)
  
  bdd_pred <- bdd_danish_sensi_1 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
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
  x_min <- min(bdd_danish_sensi_1[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Base model excluding NEFL oulier")
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, signif_vars_labels)



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

heterogeneity_adjusted_sensi_1 <- data.frame(explanatory = character(),
                                            model = factor(), 
                                            p_value_heterogeneity = numeric(), 
                                            stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi, data = bdd_danish_sensi_1)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_1)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_sensi_1 <- rbind(heterogeneity_adjusted_sensi_1, 
                                         data.frame(explanatory = var,
                                                    model = "adjusted", 
                                                    analysis = "sensi_1",
                                                    p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1)
  model_summary <- summary(model)
  model2_gam_sensi_1[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model2_gam_sensi_1, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_adjusted_gam_sensi_1 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + smoking_2cat_i + bmi"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_1)
  
  bdd_pred <- bdd_danish_sensi_1 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
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
  x_min <- min(bdd_danish_sensi_1[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Adjusted model excluding NEFL oulier")
  
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
  set_names(signif_vars)


rm(pvals, signif_vars, 
   bdd_danish_sensi_1, proteomic_sensi_1, proteomic_sd_sensi_1, proteomic_quart_sensi_1, proteomic_quart_med_sensi_1, 
   signif_vars_labels)


# Sensi 2 - Removing cases and their controls with follow_up < 5 years ----
# bdd_danish |>                                              # have to remove 17 cases with follow-up<60 months (and their controls) -> 51 people
#   filter(follow_up<60) |>
#   select(sample, match, follow_up, everything()) |>
#   arrange(follow_up) 

bdd_danish_sensi_2 <- 
  bdd_danish |>
  group_by(match) |>
  filter(any(als == 1 & follow_up > 60)) |>
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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_2)
  model_summary <- summary(model)
  model1_gam_sensi_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model1_gam_sensi_2, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_base_gam_sensi_2 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_2)
  
  bdd_pred <- bdd_danish_sensi_2 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
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
  x_min <- min(bdd_danish_sensi_2[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Base model stratified to follow-up>5years")
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, 
   signif_vars_labels)




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
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_2)
  model_summary <- summary(model)
  model2_gam_sensi_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model2_gam_sensi_2, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_adjusted_gam_sensi_2 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + smoking_2cat_i + bmi"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_2)
  
  bdd_pred <- bdd_danish_sensi_2 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
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
  x_min <- min(bdd_danish_sensi_2[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Adjusted model stratified to follow-up>5years")
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, 
   bdd_danish_sensi_2, 
   signif_vars_labels)



# Sensi 1 + sensi 2 - Removing the oulier for NEFL + removing cases and their controls with follow_up < 5 years ----

# bdd_danish |>                                              # have to remove 17 cases with follow-up<60 months (and their controls) -> 51 people
#   filter(follow_up<60) |>
#   select(sample, match, follow_up, everything()) |>
#   arrange(follow_up) 

bdd_danish_sensi_1_2 <- 
  bdd_danish |>
  group_by(match) |>                                                            # sensi 2 : we remove cases (and their controls) with follow-up <60 months to see if the association with NEFL really happens a long time before ALS diagnosis 
  filter(any(als == 1 & follow_up > 60)) |>
  ungroup()|>
  mutate(                                                                       # sensi 1 : we remove NEFL in match 159 because it's an outlier
    proteomic_neuro_explo_NEFL_sd_sensi_1_2 = 
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL_sd), 
    proteomic_neuro_explo_NEFL_sensi_1_2 = 
      ifelse(match == 159, NA, proteomic_neuro_explo_NEFL)) |>
  mutate(across(all_of(proteomic),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd_sensi_1_2")) |>
  mutate(across(all_of(proteomic), ~ factor(ntile(.x, 4),                           
                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart_sensi_1_2")) |>
  mutate(across(
    all_of(proteomic),
    ~ {
      cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)               
      quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
      quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
      quart_meds[quartiles]                                                     
    },
    .names = "{.col}_quart_med_sensi_1_2"
  ))


proteomic_sensi_1_2 <-
  proteomic |> 
  str_replace("proteomic_neuro_explo_NEFL", "proteomic_neuro_explo_NEFL_sensi_1_2")
proteomic_sd_sensi_1_2 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sd_sensi_1_2")
proteomic_quart_sensi_1_2 <-
  proteomic_quart |> 
  str_replace("proteomic_neuro_explo_NEFL_quart", "proteomic_neuro_explo_NEFL_quart_sensi_1_2")
proteomic_quart_med_sensi_1_2 <-
  proteomic_quart_med |> 
  str_replace("proteomic_neuro_explo_NEFL_quart_med", "proteomic_neuro_explo_NEFL_quart_med_sensi_1_2")



### model 1 sd ----


model1_sd_sensi_1_2 <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_2 <- rbind(model1_sd_sensi_1_2, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model1_sd_sensi_1_2 <- model1_sd_sensi_1_2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


## model 1 quartiles ----

model1_quart_sensi_1_2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_quart_sensi_1_2 <- rbind(model1_quart_sensi_1_2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model1_quart_sensi_1_2 <- model1_quart_sensi_1_2 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_2Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_2Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_2Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)


#### heterogeneity test 
heterogeneity_base_sensi_1_2 <- data.frame(explanatory = character(),
                                           model = factor(),
                                           p_value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_2) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish_sensi_1_2)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_1_2)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_sensi_1_2 <- rbind(heterogeneity_base_sensi_1_2, 
                                        data.frame(explanatory = var,
                                                   model = "base", 
                                                   analysis = "sensi_1_2",
                                                   p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)


#### trend test 
trend_base_sensi_1_2 <- data.frame(explanatory = character(),
                                   model = factor(), 
                                   p_value_trend = numeric(), 
                                   stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test <- 
    clogit(formula, data = bdd_danish_sensi_1_2) |>
    summary() 
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_base_sensi_1_2 <- rbind(trend_base_sensi_1_2, 
                                data.frame(explanatory = var,
                                           model = "base", 
                                           analysis = "sensi_1_2",
                                           p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 1 gams ---- 


model1_gam_sensi_1_2 <- list()

for (var in proteomic_sensi_1_2) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_2)
  model_summary <- summary(model)
  model1_gam_sensi_1_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model1_gam_sensi_1_2, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_base_gam_sensi_1_2 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_1_2)
  
  bdd_pred <- bdd_danish_sensi_1_2 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
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
  x_min <- min(bdd_danish_sensi_1_2[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Base model stratified to follow-up>5years")
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, 
   signif_vars_labels)




### model 2 sd ----


model2_sd_sensi_1_2 <- data.frame(explanatory = character(),
                                  term = integer(),
                                  OR = numeric(),
                                  lower_CI = numeric(),
                                  upper_CI = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_2 <- rbind(model2_sd_sensi_1_2, data.frame(explanatory = var,
                                                               term = term, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               p_value = p_value))
}

model2_sd_sensi_1_2 <- model2_sd_sensi_1_2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



## model 2 quartiles ----
model2_quart_sensi_1_2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  model <- clogit(formula, data = bdd_danish_sensi_1_2)
  model_summary <- tidy(model) |> filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  model2_quart_sensi_1_2 <- rbind(model2_quart_sensi_1_2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_quart_sensi_1_2 <- model2_quart_sensi_1_2 |> 
  mutate(
    term = case_when(
      grepl("_quart_sensi_1_2Q2", term) ~ "Quartile 2",
      grepl("_quart_sensi_1_2Q3", term) ~ "Quartile 3",
      grepl("_quart_sensi_1_2Q4", term) ~ "Quartile 4",
      grepl("_quartQ2", term) ~ "Quartile 2",
      grepl("_quartQ3", term) ~ "Quartile 3",
      grepl("_quartQ4", term) ~ "Quartile 4",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_2") |>
  select(explanatory, model, everything())
rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)

### heterogeneity tests

heterogeneity_adjusted_sensi_1_2 <- data.frame(explanatory = character(),
                                               model = factor(), 
                                               p_value_heterogeneity = numeric(), 
                                               stringsAsFactors = FALSE)

for (var in proteomic_quart_sensi_1_2) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi, data = bdd_danish_sensi_1_2)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test_2 <- clogit(formula, data = bdd_danish_sensi_1_2)
  
  anova <- anova(test_1, test_2, test = "LR")
  p_value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_sensi_1_2 <- rbind(heterogeneity_adjusted_sensi_1_2, 
                                            data.frame(explanatory = var,
                                                       model = "adjusted", 
                                                       analysis = "sensi_1_2",
                                                       p_value_heterogeneity = p_value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p_value_heterogeneity)



### trend tests

trend_adjusted_sensi_1_2 <- data.frame(explanatory = character(),
                                       model = factor(), 
                                       p_value_trend = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in proteomic_quart_med_sensi_1_2) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  test <- clogit(formula, data = bdd_danish_sensi_1_2) |> summary()
  p_value_trend <- test$coefficients[var, "Pr(>|z|)"]
  
  trend_adjusted_sensi_1_2 <- rbind(trend_adjusted_sensi_1_2, 
                                    data.frame(explanatory = var,
                                               model = "adjusted", 
                                               analysis = "sensi_1_2",
                                               p_value_trend = p_value_trend))
}
rm(var, test, formula, p_value_trend)





### model 2 gams ---- 

model2_gam_sensi_1_2 <- list()

for (var in proteomic_sensi_1_2) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish_sensi_1_2)
  model_summary <- summary(model)
  model2_gam_sensi_1_2[[var]] <- model_summary
}

rm(var, formula, model, model_summary)


pvals <- sapply(model2_gam_sensi_1_2, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_adjusted_gam_sensi_1_2 <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + smoking_2cat_i + bmi"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_sensi_1_2)
  
  bdd_pred <- bdd_danish_sensi_1_2 |>                                                     # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
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
  x_min <- min(bdd_danish_sensi_1_2[[var]], na.rm = TRUE)
  x_label <- signif_vars_labels[var] 
  
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
    ggtitle("Adjusted model stratified to follow-up>5years")
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, 
   bdd_danish_sensi_1_2, 
   proteomic_sensi_1_2,
   proteomic_sd_sensi_1_2,
   proteomic_quart_sensi_1_2, 
   proteomic_quart_med_sensi_1_2,
   signif_vars_labels)





# Sensi 3 - stratifiyng the analysis in tertiles of follow-up duration  ----
### Visualization ----
bdd_danish_sensi_3 <- bdd_danish |>                                             # just for visualisation of follow up distribution 
  group_by(match) |>
  mutate(follow_up_years = follow_up/12, 
         follow_up_ter = follow_up_years[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE))

sensi_3_table_follow_up <- bdd_danish_sensi_3 |>                                # for visualisation of follow up distribution 
  select(follow_up_years, follow_up_ter) |>
  tbl_summary(by = follow_up_ter)

sensi_3_densityplot_follow_up <-                                                # for visualisation of follow up distribution 
  bdd_danish_sensi_3 |>
  filter(als ==1) |>
  ggplot() +
  aes(x = follow_up_years) +
  geom_density(fill = "lightgray") +
  geom_vline(
    xintercept = quantile(bdd_danish_sensi_3$follow_up_years, probs = c(1/3, 2/3), na.rm = TRUE),
    color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = min(bdd_danish_sensi_3$follow_up_years, na.rm = TRUE),
           y = Inf, label = "Tertile 1", color = "red", hjust = 0, vjust = 2) +
  annotate("text",
           x = median(bdd_danish_sensi_3$follow_up_years, na.rm = TRUE),  # centre entre les deux lignes
           y = Inf, label = "Tertile 2", color = "red", vjust = 2) +
  annotate("text",
           x = max(bdd_danish_sensi_3$follow_up_years, na.rm = TRUE),
           y = Inf, label = "Tertile 3", color = "red", hjust = 1, vjust = 2) +
  labs(x = "Follow-up duration from baseline to diagnosis (years)")

bdd_danish_sensi_3 <- bdd_danish |>
  group_by(match) |>
  mutate(follow_up_ter = follow_up[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE, 
      labels = c("Tertile 1", "Tertile 2", "Tertile 3")))
bdd_danish_sensi_3_T1 <-                                                        # stratifiyng the data on tertiles of follow up duration 
  bdd_danish_sensi_3 |>
  filter(follow_up_ter == "Tertile 1")
bdd_danish_sensi_3_T2 <- 
  bdd_danish_sensi_3 |>
  filter(follow_up_ter == "Tertile 2")
bdd_danish_sensi_3_T3 <- 
  bdd_danish_sensi_3 |>
  filter(follow_up_ter == "Tertile 3")

### model 1 T1 ----
model1_sd_sensi_3_T1 <- data.frame(explanatory = character(),
                                term = integer(),
                                OR = numeric(),
                                lower_CI = numeric(),
                                upper_CI = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_3_T1 <- rbind(model1_sd_sensi_3_T1, data.frame(explanatory = var,
                                                           term = term, 
                                                           OR = OR,
                                                           lower_CI = lower_CI,
                                                           upper_CI = upper_CI,
                                                           p_value = p_value))
}

model1_sd_sensi_3_T1 <- model1_sd_sensi_3_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_3_T1") |>
  select(explanatory, model, everything())

### model 1 T2 ----
model1_sd_sensi_3_T2 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_3_T2 <- rbind(model1_sd_sensi_3_T2, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_3_T2 <- model1_sd_sensi_3_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_3_T2") |>
  select(explanatory, model, everything())


### model 1 T3 ----
model1_sd_sensi_3_T3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_3_T3 <- rbind(model1_sd_sensi_3_T3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_3_T3 <- model1_sd_sensi_3_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_3_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



### model 2 T1 ----
model2_sd_sensi_3_T1 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T1)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_3_T1 <- rbind(model2_sd_sensi_3_T1, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_3_T1 <- model2_sd_sensi_3_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_3_T1") |>
  select(explanatory, model, everything())

### model 2 T2 ----
model2_sd_sensi_3_T2 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_3_T2 <- rbind(model2_sd_sensi_3_T2, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_3_T2 <- model2_sd_sensi_3_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_3_T2") |>
  select(explanatory, model, everything())


### model 2 T3 ----
model2_sd_sensi_3_T3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_3_T3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_3_T3 <- rbind(model2_sd_sensi_3_T3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model2_sd_sensi_3_T3 <- model2_sd_sensi_3_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_3_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var, 
   bdd_danish_sensi_3, bdd_danish_sensi_3_T1, bdd_danish_sensi_3_T2, bdd_danish_sensi_3_T3)



# Sensi 1 + sensi 3 - Removing the oulier for NEFL + stratifiyng the analysis in tertiles of follow-up duration  ----

bdd_danish_sensi_1_3 <- 
  bdd_danish  |>
  mutate(
    proteomic_neuro_explo_NEFL_sd_sensi_1_3 = ifelse(match == 159, NA, proteomic_neuro_explo_NEFL_sd)) |>
  group_by(match) |>
  mutate(follow_up_ter = follow_up[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE, 
      labels = c("Tertile 1", "Tertile 2", "Tertile 3")))
bdd_danish_sensi_1_3_T1 <- 
  bdd_danish_sensi_1_3 |>
  filter(follow_up_ter == "Tertile 1")
bdd_danish_sensi_1_3_T2 <- 
  bdd_danish_sensi_1_3 |>
  filter(follow_up_ter == "Tertile 2")
bdd_danish_sensi_1_3_T3 <- 
  bdd_danish_sensi_1_3 |>
  filter(follow_up_ter == "Tertile 3")

proteomic_sd_sensi_1_3 <-
  proteomic_sd |> 
  str_replace("proteomic_neuro_explo_NEFL_sd", "proteomic_neuro_explo_NEFL_sd_sensi_1_3")

### model 1 T1 ----
model1_sd_sensi_1_3_T1 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T1)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3_T1 <- rbind(model1_sd_sensi_1_3_T1, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_1_3_T1 <- model1_sd_sensi_1_3_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_T1") |>
  select(explanatory, model, everything())


### model 1 T2 ----
model1_sd_sensi_1_3_T2 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T2)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3_T2 <- rbind(model1_sd_sensi_1_3_T2, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_1_3_T2 <- model1_sd_sensi_1_3_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_T2") |>
  select(explanatory, model, everything())



### model 1 T3 ----
model1_sd_sensi_1_3_T3 <- data.frame(explanatory = character(),
                                   term = integer(),
                                   OR = numeric(),
                                   lower_CI = numeric(),
                                   upper_CI = numeric(),
                                   p_value = numeric(),
                                   stringsAsFactors = FALSE)


for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T3)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model1_sd_sensi_1_3_T3 <- rbind(model1_sd_sensi_1_3_T3, data.frame(explanatory = var,
                                                                 term = term, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 p_value = p_value))
}

model1_sd_sensi_1_3_T3 <- model1_sd_sensi_1_3_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "base", 
    analysis = "sensi_1_3_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var)



### model 2 T1 ----
model2_sd_sensi_1_3_T1 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T1)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3_T1 <- rbind(model2_sd_sensi_1_3_T1, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_3_T1 <- model2_sd_sensi_1_3_T1 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_T1") |>
  select(explanatory, model, everything())


### model 2 T2 ----
model2_sd_sensi_1_3_T2 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T2)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3_T2 <- rbind(model2_sd_sensi_1_3_T2, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_3_T2 <- model2_sd_sensi_1_3_T2 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_T2") |>
  select(explanatory, model, everything())



### model 2 T3 ----
model2_sd_sensi_1_3_T3 <- data.frame(explanatory = character(),
                                     term = integer(),
                                     OR = numeric(),
                                     lower_CI = numeric(),
                                     upper_CI = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)


for (var in proteomic_sd_sensi_1_3) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi"))
  
  model <- clogit(formula, data = bdd_danish_sensi_1_3_T3)
  
  model_summary <- tidy(model)  |> filter(grepl(paste0("^", var), term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  term <- model_summary$term
  
  model2_sd_sensi_1_3_T3 <- rbind(model2_sd_sensi_1_3_T3, data.frame(explanatory = var,
                                                                     term = term, 
                                                                     OR = OR,
                                                                     lower_CI = lower_CI,
                                                                     upper_CI = upper_CI,
                                                                     p_value = p_value))
}

model2_sd_sensi_1_3_T3 <- model2_sd_sensi_1_3_T3 |> 
  mutate(
    term = case_when(
      grepl("_sd", term) ~ "Continuous",
      TRUE ~ NA_character_),
    model = "adjusted", 
    analysis = "sensi_1_3_T3") |>
  select(explanatory, model, everything())

rm(model, lower_CI, upper_CI, term, formula, p_value, OR, model_summary, var, 
   bdd_danish_sensi_1_3, bdd_danish_sensi_1_3_T1, bdd_danish_sensi_1_3_T2, bdd_danish_sensi_1_3_T3, 
   proteomic_sd_sensi_1_3)




# Merging the results ----
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

                          model1_sd_sensi_1_2,                                  # sensi 1_2
                          model2_sd_sensi_1_2,
                          model1_quart_sensi_1_2,
                          model2_quart_sensi_1_2,

                          model1_sd_sensi_3_T1,                                 # sensi 3
                          model1_sd_sensi_3_T2,
                          model1_sd_sensi_3_T3,
                          model2_sd_sensi_3_T1,
                          model2_sd_sensi_3_T2,
                          model2_sd_sensi_3_T3,

                          model1_sd_sensi_1_3_T1,                               # sensi 1_3
                          model1_sd_sensi_1_3_T2,
                          model1_sd_sensi_1_3_T3,
                          model2_sd_sensi_1_3_T1,
                          model2_sd_sensi_1_3_T2,
                          model2_sd_sensi_1_3_T3
                          ) |> 
  mutate(explanatory = gsub("_quart_med_sensi_1_2", "", explanatory), 
         explanatory = gsub("_quart_sensi_1_2", "", explanatory), 
         explanatory = gsub("_quart_med_sensi_1", "", explanatory), 
         explanatory = gsub("_quart_sensi_1", "", explanatory), 
         explanatory = gsub("_quart", "", explanatory), 
         
         explanatory = gsub("_sd_sensi_1_2", "", explanatory),
         explanatory = gsub("_sd_sensi_1", "", explanatory),
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
      term == "Continuous" & analysis == "sensi_1_2" & model == "base",                            
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
      term == "Continuous" & analysis == "sensi_1_2" & model == "adjusted",                            
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



heterogeneity_tests <- 
  bind_rows(heterogeneity_base, heterogeneity_adjusted)  |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests <- 
  bind_rows(trend_base, trend_adjusted) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1 <- 
  bind_rows(heterogeneity_base_sensi_1, heterogeneity_adjusted_sensi_1) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1 <- 
  bind_rows(trend_base_sensi_1, trend_adjusted_sensi_1) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_2 <- 
  bind_rows(heterogeneity_base_sensi_2, heterogeneity_adjusted_sensi_2) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_2 <- 
  bind_rows(trend_base_sensi_2, trend_adjusted_sensi_2) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests_sensi_1_2 <- 
  bind_rows(heterogeneity_base_sensi_1_2, heterogeneity_adjusted_sensi_1_2) |>
  mutate(explanatory = gsub("_quart", "", explanatory))

trend_tests_sensi_1_2 <- 
  bind_rows(trend_base_sensi_1_2, trend_adjusted_sensi_1_2) |>
  mutate(explanatory = gsub("_quart_med", "", explanatory))

heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_2)
heterogeneity_tests <- bind_rows(heterogeneity_tests, heterogeneity_tests_sensi_1_2)

trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_2)
trend_tests <- bind_rows(trend_tests, trend_tests_sensi_1_2)

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
   model1_sd_sensi_1_2, model1_quart_sensi_1_2, model2_sd_sensi_1_2, model2_quart_sensi_1_2, 

   model1_sd_sensi_3_T1, 
   model1_sd_sensi_3_T2, 
   model1_sd_sensi_3_T3, 
   model1_sd_sensi_1_3_T1, 
   model1_sd_sensi_1_3_T2, 
   model1_sd_sensi_1_3_T3, 
   
   model2_sd_sensi_3_T1, 
   model2_sd_sensi_3_T2, 
   model2_sd_sensi_3_T3, 
   model2_sd_sensi_1_3_T1, 
   model2_sd_sensi_1_3_T2, 
   model2_sd_sensi_1_3_T3, 
   
   heterogeneity_base, heterogeneity_adjusted, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests, 
   
   heterogeneity_base_sensi_1, heterogeneity_adjusted_sensi_1, 
   trend_base_sensi_1, trend_adjusted_sensi_1, 
   heterogeneity_tests_sensi_1, trend_tests_sensi_1, 
   
   heterogeneity_base_sensi_2, heterogeneity_adjusted_sensi_2, 
   trend_base_sensi_2, trend_adjusted_sensi_2, 
   heterogeneity_tests_sensi_2, trend_tests_sensi_2, 
   
   heterogeneity_base_sensi_1_2, heterogeneity_adjusted_sensi_1_2, 
   trend_base_sensi_1_2, trend_adjusted_sensi_1_2, 
   heterogeneity_tests_sensi_1_2, trend_tests_sensi_1_2)


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


# Additional analysis 2 - NEFL over time ----
## LOESS ----
# The locally estimated scatterplot smoothing (LOESS) curve illustrates the changes in the ratio of NEFL levels in each case–control pair/triplet over time. Each point corresponds to 1 case–control pair/triplet.

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

# Visualisation avec courbe LOESS
figure_NEFL_over_time <- 
  ggplot(ratios, aes(x = follow_up_neg, y = ratio_proteomic_neuro_explo_NEFL)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.75) +
  labs(
    x = "Time to ALS diagnosis (years)",
    y = "Ratio of NEFL in case and matched controls") +
  theme_minimal(base_size = 14)
rm(ratio)

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

# Visualisation avec courbe LOESS
figure_NEFL_over_time_sensi_1 <- 
  ggplot(ratios_sensi_1, aes(x = follow_up_neg, y = ratio_proteomic_neuro_explo_NEFL)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.75) +
  labs(
    x = "Time to ALS diagnosis (years)",
    y = "Ratio of NEFL in case and matched controls") +
  theme_lucid(base_size = 14)
rm(ratio_sensi_1)

## GAM ----

# Figures and Tables ----
### Table covariates - als survival ----
covar

## Table proteomic - als occurence - base and adjusted sd ----
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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
#     "1All models are matched on age at baseline and sex. 
#     2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
#     3CI: Confidence interval.") |>
#   add_header(
#     "explanatory" = "Pre-disease serum proteins", 
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



## Table proteomic - als occurence - base and adjusted quart ----
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and age at baseline. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and age at baseline, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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
#     "1All models are matched on age at baseline and sex. 
#     2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
#     3CI: Confidence interval.
#     4Heterogeneity tests in outcome value across protein quartiles, matched on sex and age at baseline. 
#     5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline.") |>
#   add_header(
#     "explanatory" = "Pre-disease serum proteins", 
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



## Figure proteomic - als occurence - base sd ----
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



## Figure proteomic - als occurence - adjusted sd ----
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



## Figure proteomic - als occurrence - base gam ----
pvals <- sapply(model1_gam, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

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
  x_label <- signif_vars_labels[var] 
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, signif_vars_labels)

## Figure proteomic - als occurrence - adjusted gam ----
pvals <- sapply(model2_gam, function(m) m$s.table[1, "p-value"])
signif_vars <- names(pvals)[pvals < 0.05]
signif_vars_labels <- set_names(str_replace(signif_vars, 
                                            "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                            ""), 
                                signif_vars)

plot_adjusted_gam <- map(signif_vars, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age + smoking_2cat_i + bmi"))   # ou {paste(covariates, collapse = ' + ')}
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish |>                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
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
  x_label <- signif_vars_labels[var] 
  
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
  set_names(signif_vars)

rm(pvals, signif_vars, signif_vars_labels)


## Table proteomic - als occurence - base and adjusted sd (sensi_1) ----
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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





## Table proteomic - als occurence - base and adjusted quart (sensi_1) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and age at baseline. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and age at baseline, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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

## Figure proteomic - als occurence - base sd (sensi_1) ----
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

## Table proteomic - als occurence - base and adjusted sd (sensi_2) ----
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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


## Table proteomic - als occurence - base and adjusted quart (sensi_2) ----
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and age at baseline. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and age at baseline, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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


## Figure proteomic - als occurence - base sd (sensi_2) ----
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



## Table proteomic - als occurence - base and adjusted sd (sensi_1_2) ----
proteomic_sd_ALS_table_sensi_1_2 <-                                                       # select both base and adjusted results
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only continuous results
           term == "Continuous" & 
           analysis == "sensi_1_2") |>            
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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




## Table proteomic - als occurence - base and adjusted quart (sensi_1_2) ----
extra_rows <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select only quartile results
           term != "Continuous" & 
           analysis == "sensi_1_2") |>           
  group_by(explanatory) |>                                                      # select explanatory vars with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  distinct(protein_group, explanatory) |> 
  mutate(
    quartiles = "Quartile 1",
    "OR_base" = '-', "95% CI_base" = '-', "p_value_base" = '', "p_value_heterogeneity_base" = '', "p_value_trend_base" = '',
    "OR_adjusted" = '-', "95% CI_adjusted" = '-', "p_value_adjusted" = '', "p_value_heterogeneity_adjusted" = '', "p_value_trend_adjusted" = '')

proteomic_quart_ALS_table_sensi_1_2 <- 
  main_results |>
  filter(model %in% c("base", "adjusted") &                                     # select quartile results
           term != "Continuous" & 
           analysis == "sensi_1_2") |>             
  group_by(explanatory) |>                                                      # select explanatory var s with at least one quartile significant 
  filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>  
  ungroup() |>
  select(model, protein_group, explanatory, term, OR, "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend") |>
  pivot_wider(names_from = "model", values_from = c("OR", "95% CI", "p_value", "p_value_heterogeneity", "p_value_trend")) |>
  select(protein_group, explanatory, quartiles = term, contains("base"), contains("adjusted")) 

proteomic_quart_ALS_table_sensi_1_2 <- 
  proteomic_quart_ALS_table_sensi_1_2 |>
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
    "1All models are matched on age at baseline and sex. Adjusted models further account for smoking and body mass index. 
    2Estimated risk of ALS associated with a one standard deviation increase in pre-disease serum concentration of proteins.
    3CI: Confidence interval.
    4Heterogeneity tests in outcome value across protein quartiles, matched on sex and age at baseline. 
    5Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline. 
    6Heterogeneity tests in outcome value across POP quartiles, matched on sex and age at baseline, and adjusted for smoking and body mass index.
    7Trend tests using continuous variables whose values corresponded to the quartile specific median protein levels, matched on sex and age at baseline, and adjusted for smoking and body mass index.") |>
  add_header(
    "explanatory" = "Pre-disease serum proteins", 
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



## Figure proteomic - als occurence - base sd (sensi_1_2) ----
proteomic_sd_ALS_base_figure_sensi_1_2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_2") |>
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



## Figure proteomic - als occurence - base sd (sensi_3) ----
proteomic_sd_ALS_base_figure_sensi_3_T1 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_3_T1") |>
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

proteomic_sd_ALS_base_figure_sensi_3_T2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_3_T2") |>
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

proteomic_sd_ALS_base_figure_sensi_3_T3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_3_T3") |>
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



## Figure proteomic - als occurence - base sd (sensi_1_3) ----
proteomic_sd_ALS_base_figure_sensi_1_3_T1 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T1") |>
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

proteomic_sd_ALS_base_figure_sensi_1_3_T2 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T2") |>
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

proteomic_sd_ALS_base_figure_sensi_1_3_T3 <- 
  main_results |>
  filter(model == "base" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T3") |>
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



## Figure proteomic - als occurence - adjusted sd (sensi_1) ----
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



## Figure proteomic - als occurence - adjusted sd (sensi_2) ----
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



## Figure proteomic - als occurence - adjusted sd (sensi_1_2) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_2") |>
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



## Figure proteomic - als occurence - adjusted sd (sensi_3) ----
proteomic_sd_ALS_adjusted_figure_sensi_3_T1 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_3_T1") |>
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

proteomic_sd_ALS_adjusted_figure_sensi_3_T2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_3_T2") |>
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

proteomic_sd_ALS_adjusted_figure_sensi_3_T3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_3_T3") |>
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



## Figure proteomic - als occurence - adjusted sd (sensi_1_3) ----
proteomic_sd_ALS_adjusted_figure_sensi_1_3_T1 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T1") |>
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

proteomic_sd_ALS_adjusted_figure_sensi_1_3_T2 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T2") |>
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

proteomic_sd_ALS_adjusted_figure_sensi_1_3_T3 <- 
  main_results |>
  filter(model == "adjusted" & 
           term == "Continuous" & 
           analysis == "sensi_1_3_T3") |>
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
      model2_gam = model2_gam,
      plot_base_gam = plot_base_gam, 
      plot_adjusted_gam = plot_adjusted_gam),
    
    sensi_1 = list(
      proteomic_sd_ALS_table_sensi_1 = proteomic_sd_ALS_table_sensi_1, 
      proteomic_quart_ALS_table_sensi_1 = proteomic_quart_ALS_table_sensi_1, 
      
      proteomic_sd_ALS_base_figure_sensi_1 = proteomic_sd_ALS_base_figure_sensi_1, 
      model1_gam_sensi_1 = model1_gam_sensi_1, 
      plot_base_gam_sensi_1 = plot_base_gam_sensi_1, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_1 = proteomic_sd_ALS_adjusted_figure_sensi_1, 
      model2_gam_sensi_1 = model2_gam_sensi_1, 
      plot_adjusted_gam_sensi_1 = plot_adjusted_gam_sensi_1), 
    
    sensi_2 = list(
      proteomic_sd_ALS_table_sensi_2 = proteomic_sd_ALS_table_sensi_2, 
      proteomic_quart_ALS_table_sensi_2 = proteomic_quart_ALS_table_sensi_2, 
      
      proteomic_sd_ALS_base_figure_sensi_2 = proteomic_sd_ALS_base_figure_sensi_2, 
      model1_gam_sensi_2 = model1_gam_sensi_2, 
      plot_base_gam_sensi_2 = plot_base_gam_sensi_2, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_2 = proteomic_sd_ALS_adjusted_figure_sensi_2, 
      model2_gam_sensi_2 = model2_gam_sensi_2, 
      plot_adjusted_gam_sensi_2 = plot_adjusted_gam_sensi_2), 
    
    sensi_1_2 = list(
      proteomic_sd_ALS_table_sensi_1_2 = proteomic_sd_ALS_table_sensi_1_2, 
      proteomic_quart_ALS_table_sensi_1_2 = proteomic_quart_ALS_table_sensi_1_2, 
      
      proteomic_sd_ALS_base_figure_sensi_1_2 = proteomic_sd_ALS_base_figure_sensi_1_2, 
      model1_gam_sensi_1_2 = model1_gam_sensi_1_2, 
      plot_base_gam_sensi_1_2 = plot_base_gam_sensi_1_2, 
      
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_2 = proteomic_sd_ALS_adjusted_figure_sensi_1_2, 
      model2_gam_sensi_1_2 = model2_gam_sensi_1_2, 
      plot_adjusted_gam_sensi_1_2 = plot_adjusted_gam_sensi_1_2), 
    
    sensi_3 = list(
      sensi_3_table_follow_up = sensi_3_table_follow_up, 
      sensi_3_densityplot_follow_up = sensi_3_densityplot_follow_up,
      
      proteomic_sd_ALS_base_figure_sensi_3_T1 = proteomic_sd_ALS_base_figure_sensi_3_T1, 
      proteomic_sd_ALS_base_figure_sensi_3_T2 = proteomic_sd_ALS_base_figure_sensi_3_T2, 
      proteomic_sd_ALS_base_figure_sensi_3_T3 = proteomic_sd_ALS_base_figure_sensi_3_T3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_3_T1 = proteomic_sd_ALS_adjusted_figure_sensi_3_T1, 
      proteomic_sd_ALS_adjusted_figure_sensi_3_T2 = proteomic_sd_ALS_adjusted_figure_sensi_3_T2, 
      proteomic_sd_ALS_adjusted_figure_sensi_3_T3 = proteomic_sd_ALS_adjusted_figure_sensi_3_T3),
    
    
    sensi_1_3 = list(
      proteomic_sd_ALS_base_figure_sensi_1_3_T1 = proteomic_sd_ALS_base_figure_sensi_1_3_T1, 
      proteomic_sd_ALS_base_figure_sensi_1_3_T2 = proteomic_sd_ALS_base_figure_sensi_1_3_T2, 
      proteomic_sd_ALS_base_figure_sensi_1_3_T3 = proteomic_sd_ALS_base_figure_sensi_1_3_T3, 
      
      proteomic_sd_ALS_adjusted_figure_sensi_1_3_T1 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_T1, 
      proteomic_sd_ALS_adjusted_figure_sensi_1_3_T2 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_T2, 
      proteomic_sd_ALS_adjusted_figure_sensi_1_3_T3 = proteomic_sd_ALS_adjusted_figure_sensi_1_3_T3), 
    
    additional_analysis_1 = list(
      OR_distribution = OR_distribution, 
      boxplot_OR = boxplot_OR, 
      densityplot_OR = densityplot_OR), 
    additional_analysis_2 = list(
      figure_NEFL_over_time = figure_NEFL_over_time ,
      figure_NEFL_over_time_sensi_1 = figure_NEFL_over_time_sensi_1))

rm(covar, 
   main_results, 
   proteomic_sd_ALS_table,
   proteomic_quart_ALS_table, 
   proteomic_sd_ALS_base_figure, 
   proteomic_sd_ALS_adjusted_figure, 
   model1_gam, 
   model2_gam,
   plot_base_gam, 
   plot_adjusted_gam, 
   
   proteomic_sd_ALS_table_sensi_1, 
   proteomic_quart_ALS_table_sensi_1, 
   proteomic_sd_ALS_base_figure_sensi_1, 
   model1_gam_sensi_1, 
   plot_base_gam_sensi_1, 
   proteomic_sd_ALS_adjusted_figure_sensi_1, 
   model2_gam_sensi_1, 
   plot_adjusted_gam_sensi_1, 
   
   proteomic_quart_ALS_table_sensi_2, 
   proteomic_sd_ALS_table_sensi_2, 
   proteomic_sd_ALS_base_figure_sensi_2, 
   model1_gam_sensi_2, 
   plot_base_gam_sensi_2, 
   proteomic_sd_ALS_adjusted_figure_sensi_2, 
   model2_gam_sensi_2, 
   plot_adjusted_gam_sensi_2, 
   
   proteomic_sd_ALS_table_sensi_1_2, 
   proteomic_quart_ALS_table_sensi_1_2, 
   proteomic_sd_ALS_base_figure_sensi_1_2,
   model1_gam_sensi_1_2, 
   plot_base_gam_sensi_1_2, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_2,
   model2_gam_sensi_1_2, 
   plot_adjusted_gam_sensi_1_2, 
   
   sensi_3_table_follow_up, 
   sensi_3_densityplot_follow_up, 
   proteomic_sd_ALS_base_figure_sensi_3_T1, 
   proteomic_sd_ALS_base_figure_sensi_3_T2, 
   proteomic_sd_ALS_base_figure_sensi_3_T3, 
   proteomic_sd_ALS_adjusted_figure_sensi_3_T1, 
   proteomic_sd_ALS_adjusted_figure_sensi_3_T2, 
   proteomic_sd_ALS_adjusted_figure_sensi_3_T3, 
   
   proteomic_sd_ALS_base_figure_sensi_1_3_T1, 
   proteomic_sd_ALS_base_figure_sensi_1_3_T2, 
   proteomic_sd_ALS_base_figure_sensi_1_3_T3, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_T1, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_T2, 
   proteomic_sd_ALS_adjusted_figure_sensi_1_3_T3, 
   
   OR_distribution, 
   boxplot_OR, 
   densityplot_OR, 
   figure_NEFL_over_time, 
   figure_NEFL_over_time_sensi_1)
  

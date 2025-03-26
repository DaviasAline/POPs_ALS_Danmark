
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


sensitivity_results_not_summed_quart <- bind_rows(model1_quart_not_summed, model2_quart_not_summed) 

sensitivity_results_not_summed_quart <- sensitivity_results_not_summed_quart |>
  mutate(variable = gsub("_quart", "", variable), 
         variable = gsub("BDE", "PBDE", variable),
         variable = gsub('_', '-', variable),
         variable = fct_recode(variable, "p,p'-DDE" = 'pp-DDE',  "p,p'-DDT" ="pp-DDT"),
         OR = format(OR, nsmall = 1, digits = 1),
         lower_CI = format(lower_CI, nsmall = 1, digits = 1),
         upper_CI =  format(upper_CI, nsmall = 1, digits = 1),
         p.value_raw = p.value, 
         p.value = ifelse(p.value < 0.01, "<0.01", format(p.value, nsmall = 1, digits = 1), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) %>%
  arrange(variable) %>%
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value"), 
         lower_CI, upper_CI) 


extra_rows <- sensitivity_results_not_summed_quart %>%
  distinct(variable) %>% 
  mutate(
    quartiles = "1",
    "OR_base_not_summed" = '-', "95%CI_base_not_summed" = '-', "p.value_base_not_summed" = '',
    "OR_adjusted_not_summed" = '-', "95%CI_adjusted_not_summed" = '-', "p.value_adjusted_not_summed" = '')

sensitivity_results_not_summed_quart <- 
  sensitivity_results_not_summed_quart |>
  select(-lower_CI, -upper_CI, - p.value_raw) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart"), 
         contains("copollutant_quart")) 
colnames(sensitivity_results_not_summed_quart) <- gsub('_quart', '', colnames(sensitivity_results_not_summed_quart))

sensitivity_results_not_summed_quart <- sensitivity_results_not_summed_quart |>
  mutate_if(is.numeric, as.character) |>
  bind_rows(extra_rows) |>
  mutate(
    # variable = fct_relevel(variable, "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE" ),
    variable = gsub("_", "-", variable)) |>
  arrange(variable, quartiles) |>
  rename('OR' = 'OR_base_not_summed', '95% CI' = '95%CI_base_not_summed', 'p-value' = 'p.value_base_not_summed', 
         'OR ' = 'OR_adjusted_not_summed', '95% CI ' = '95%CI_adjusted_not_summed', 'p-value ' = 'p.value_adjusted_not_summed')
sensitivity_results_not_summed_quart <- sensitivity_results_not_summed_quart |> flextable() |>
  add_footer_lines(
    "1POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
  2All models are matched for sex and age. Adjusted models further account for smoking, BMI, serum total cholesterol, marital status, and education. Co-pollutant models further include all chemicals, where the chemical of interest is categorized into quartiles, while the others are treated as continuous variables.
  3Estimated risk of ALS when exposures to POP are at quartiles 2, 3, and 4, compared to quartile 1.
  4CI: Confidence interval.") |>
  add_header(
    "variable" = "Exposures", "quartiles" = "Quartiles",
    "OR" = "Base Model", "95% CI" = "Base Model", "p-value" = "Base Model", 
    "OR " = "Adjusted Model", "95% CI " = "Adjusted Model", "p-value " = "Adjusted Model") |>
  merge_h(part = "header") |>
  merge_v(j = "variable") |>
  theme_vanilla() |>
  bold(j = "variable", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = "variable", align = "left", part = "all") |> 
  merge_at(j = "variable", part = "header") |>
  merge_at(j = "quartiles", part = "header") |>
  font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")
rm(extra_rows)



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

for (var in POPSs_included) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + 
                              smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

POPs_included_labels <- gsub("_", "-", POPs_included)
POPs_included_labels <- gsub("BDE", "PBDE", POPs_included_labels)
POPs_included_labels <- gsub("pp", "p,p'", POPs_included_labels)
pollutant_labels <- set_names(c(POPs_labels,  POPs_included))

plot_base_gamm_not_summed <- map(POPs_included, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec covariables ramenées à leur moyenne
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
  set_names(POPs_included_labels)

wrap_plots(plot_base_gamm_not_summed)


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

wrap_plots(plot_base_gamm_not_summed)


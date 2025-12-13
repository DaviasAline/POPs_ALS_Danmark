# Aline Davias
# 10/11/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.7_analyses_proteomic_ALS_occurrence.R")


# Table 1 - Subject characteristics description ---- 
# Description of the subject characteristics of the Danish Diet, Cancer and Health study cohort (sample size: 498).
table_1 <- 
  bdd_danish |>
  filter(match != 159) |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"),
    follow_up = follow_up/12) |>
  select(
    als, birth_year, baseline_age, diagnosis_age, follow_up,  
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, fS_Kol, 
    proteomic_neuro_explo_NEFL) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 0, 
                            bmi ~ 1, 
                            fS_Kol ~ 1, 
                            proteomic_neuro_explo_NEFL ~ 1), 
              label = list(proteomic_neuro_explo_NEFL ~ "Neurofilament light polypeptide (NEFL)")) |>
  bold_labels() |>
  # add_p(include = -c('diagnosis_age', 'follow_up')) |>
  add_overall() |>
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")




# Figure 1 - Descriptive figure of NEFL distribution (density and boxplots) ----
figure_1 <- results_descriptive$danish$figure_NEFL


# Figure 2 - Base and adjusted GAMs (als risk) ----
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


## main base ----

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
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "steelblue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Base model")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab("Neurofilament light polypeptide (NPX)") + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |>
  set_names(signif_vars)

rm(pvals, signif_vars, signif_vars_labels)


## main adjusted ----

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
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "steelblue", alpha = 0.2) +
    labs(x = var, y = "") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    ggtitle("Adjusted model")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab("Neurofilament light polypeptide (NPX)") + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |>
  set_names(signif_vars)


figure_2 <- 
  wrap_plots(plot_base_gam_sensi_1$proteomic_neuro_explo_NEFL_sensi_1, 
             plot_adjusted_gam_sensi_1$proteomic_neuro_explo_NEFL_sensi_1)

rm(pvals, signif_vars, 
   bdd_danish_sensi_1, proteomic_sensi_1, 
   signif_vars_labels, 
   model1_gam_sensi_1, model2_gam_sensi_1, plot_base_gam_sensi_1, plot_adjusted_gam_sensi_1)


# Figure 3 - Sensitivity analysis - Main adjusted model vs. filtered to > 5 years of follow-up adjusted model (ALS risk) ----



## main adjusted ----
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
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "steelblue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model (main analysis, n = 495)")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab("Neurofilament light polypeptide (NPX)") + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) |>
  set_names(signif_vars)


rm(pvals, signif_vars, 
   bdd_danish_sensi_1, proteomic_sensi_1, 
   signif_vars_labels, model2_gam_sensi_1)



## sensi adjusted ----
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
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "steelblue", alpha = 0.2) +
    labs(x = "Neurofilament light polypeptide (NPX)", y = "Predicted probability of ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    ggtitle("Adjusted model stratified to follow-up duration > 5 years\n(sensitivity analysis, n = 447)")
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab("Neurofilament light polypeptide (NPX)") + 
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
   signif_vars_labels, 
   model2_gam_sensi_1_2)

figure_3 <- wrap_plots(plot_adjusted_gam_sensi_1$proteomic_neuro_explo_NEFL_sensi_1, 
                       plot_adjusted_gam_sensi_1_2$proteomic_neuro_explo_NEFL_sensi_1_2)
rm(plot_adjusted_gam_sensi_1, plot_adjusted_gam_sensi_1_2)



# Figure 4 - Additional analysis (NEFL level over follow-up time) ----
figure_4 <- results_proteomic_ALS_occurrence$additional_analysis_2$figure_NEFL_over_time_sensi_1

# Figure 5 - Base and adjusted GAMs (ALS survival) ----
bdd_cases_danish <- bdd_danish |>
  filter (als == 1) |>
  filter(match != 159) |>
  select(sample, als, als_date, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age, follow_up, 
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(proteomic)) |>
  mutate(across(all_of(proteomic), 
                ~ factor(ntile(.x, 4),                                          # creation of proteomic quartiles variables (cohort and cases specific)                        
                         labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(proteomic),                                              # creation of proteomic standardized variables (cohort and cases specific)  
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  mutate(across(all_of(proteomic),
                ~ {
                  cuts <- quantile(.x, probs = seq(0, 1, 0.25), na.rm = TRUE)      
                  quartiles <- cut(.x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
                  quart_meds <- tapply(.x, quartiles, median, na.rm = TRUE)                 
                  quart_meds[quartiles]       
                },
                .names = "{.col}_quart_med"))

surv_obj <- Surv(time = bdd_cases_danish$follow_up_death,                       # set the outcomes
                 event = bdd_cases_danish$status_death)
covariates <-                                                                   # set the covariates 
  c("sex", "diagnosis_age", "smoking_2cat_i", "bmi")


vars_labels <- set_names(str_replace(proteomic, 
                                     "proteomic_immun_res_|proteomic_neuro_explo_|proteomic_metabolism_", 
                                     ""), 
                         proteomic)

## main base ----
fit_cox_gam_base <- function(var, data = bdd_cases_danish) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var,
    x_label = x_label)
}

cox_gam_results_base <- map(proteomic, fit_cox_gam_base)
rm(fit_cox_gam_base)


## main adjusted ----
fit_cox_gam_adjusted <- function(var, data = bdd_cases_danish) {
  
  outcome <- with(bdd_cases_danish, cbind(follow_up_death, status_death))
  formula_str <- paste0("outcome ~ s(", var, ") + sex + diagnosis_age + smoking_2cat_i + bmi")
  model <- gam(as.formula(formula_str), data = data, family = cox.ph())
  
  smry <- summary(model)
  edf <- format(smry$s.table[1, "edf"], nsmall = 1, digits = 1) 
  pval <- format(smry$s.table[1, "p-value"], nsmall = 2, digits = 2) 
  pval <- case_when(pval < 0.01 ~ "< 0.01", 
                    pval >0.99 ~ "> 0.99",
                    .default = format(pval, nsmall =2, digits = 2))
  x_label <- vars_labels[[var]] 
  
  plot_data <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE, pages = 0)[[1]]
  smooth_df <- data.frame(
    x = plot_data$x,
    fit = plot_data$fit,
    se = plot_data$se)
  
  list(
    model = model,
    plot_data = smooth_df,
    edf = edf,
    pval = pval,
    var = var, 
    x_label = x_label)
}

cox_gam_results_adjusted <- map(proteomic, fit_cox_gam_adjusted)
rm(fit_cox_gam_adjusted)

## plot----

all_fits_base <- map_dfr(cox_gam_results_base, "plot_data", .id = "var")
y_range_base <- range(all_fits_base$fit - 2 * all_fits_base$se,
                      all_fits_base$fit + 2 * all_fits_base$se, na.rm = TRUE)

plot_base_cox_gam_danish <- map(cox_gam_results_base, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Base model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +
    scale_x_continuous(limits = c(1, 5))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_cases_danish, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(all_fits_base, y_range_base, cox_gam_results_base)


all_fits_adjusted <- map_dfr(cox_gam_results_adjusted, "plot_data", .id = "var")
y_range_adjusted <- range(all_fits_adjusted$fit - 2 * all_fits_adjusted$se,
                          all_fits_adjusted$fit + 2 * all_fits_adjusted$se, na.rm = TRUE)


plot_adjusted_cox_gam_danish <- map(cox_gam_results_adjusted, function(res) {
  
  p1 <- ggplot(res$plot_data, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                alpha = 0.2, fill = "steelblue") +
    labs(
      title = "Adjusted model", 
      x = NULL,
      y = "LogHR (smooth estimate)") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("EDF: ", res$edf, "\np-value: ", res$pval),
      size = 4.2,
      color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(-12, 12)) +  
    scale_x_continuous(limits = c(1, 5))  +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p2 <- ggplot(bdd_cases_danish, aes_string(x = res$var, y = 1)) +
    geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Neurofilament light polypeptide (NPX)") +
    scale_x_continuous(limits = c(1, 5)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  p1 / p2 + plot_layout(heights = c(10, 1))
}) |> 
  set_names(vars_labels)

rm(vars_labels, all_fits_adjusted, y_range_adjusted, cox_gam_results_adjusted)


figure_5 <- 
  wrap_plots(plot_base_cox_gam_danish$NEFL, 
             plot_adjusted_cox_gam_danish$NEFL, 
             ncol = 2)

rm(plot_base_cox_gam_danish, plot_adjusted_cox_gam_danish)


# Export ----
## tables ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/table_1.docx")


## figures ----
ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_1.tiff",
  figure_1,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_2.tiff",
  figure_2,
  height = 5,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_3.tiff",
  figure_3,
  height = 5,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_4.tiff",
  figure_4,
  height = 4,
  width = 8,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/3.Article_NEFL_ALS/figure_5.tiff",
  figure_5,
  height = 5,
  width = 10,
  units = "in")

# figure 5 to add

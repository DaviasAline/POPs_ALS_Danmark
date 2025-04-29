# Aline Davias 
# 04.25.2025
# code for team meeting presentation 


## data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.2_analyses_POPs_ALS_occurrence.R")

## descriptif ----
covar_comp <- bdd |>
  select("study", "sex", "marital_status_2cat", "smoking_2cat", "education", 
         "bmi", "cholesterol") |>
  mutate(education = fct_recode(education,
           "Low" = "<7 years of primary school",
           "Medium" = "7-10 years of primary school",
           "High" = ">10 years of primary school",
           "Low" = "<7 years",
           "Medium" = "7-12 years",
           "High" = ">12 years"), 
         study = fct_recode(study, "EPIC" = "Danish"), 
         study = fct_relevel(study,  "FMC", "FMCF", "MFH", "EPIC")) |>
  tbl_summary(by = "study", 
              missing = "no") |>
  bold_labels()

covar_comp_cases <- bdd |>
  filter(als == 1) |>
  select(study, "baseline_age", "diagnosis_age", death_age, 
         follow_up, follow_up_death) |>
  mutate(         study = fct_recode(study, "EPIC" = "Danish"), 
                  study = fct_relevel(study,  "FMC", "FMCF", "MFH", "EPIC")) |>
  tbl_summary(by = "study", missing = "no") |>
  bold_labels()

POPs_boxplot_comp <- bdd |>
  select(study, all_of(POPs_group)) |>
  select(-"ΣPBDE") |>
  pivot_longer(cols = -study, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = study) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml)", fill = "Study") +
  coord_flip() +
  theme_lucid()

## model 1 ----
### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
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
    scale_y_continuous(limits = c(0, 1)) +  
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
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
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
  x_label <- paste(pollutant_labels[var], "without outliers")
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = paste(var, "without outliers"), y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    scale_y_continuous(limits = c(0, 1)) +  
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
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


## model 2 ----
### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
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
    labs(x = var, y = "") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
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
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels)

### gamm outliers ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
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
  x_label <- paste(pollutant_labels[var], "without outliers")
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = paste(var, "without outliers"), y = "") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
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


## model 3 ----
### gamm ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
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
    labs(x = var, y = "") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
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


### gamm outliers ----
POPs_group_outlier_bis <- setdiff(POPs_group_outlier, "PCB_4_outlier")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
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
  x_label <- paste(pollutant_labels_bis[var], "without outliers")
  
  
  # edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  # p_value <- model_summary$s.table[1, "p-value"]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = paste(var, "without outliers"), y = "") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
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
rm(POPs_group_outlier_bis, pollutant_labels_bis, model)




gaam_results_PCB_DL <- wrap_plots(plot_base_gamm$PCB_DL, 
           plot_adjusted_gamm$PCB_DL, 
           plot_copollutant_gamm$PCB_DL, 
           plot_base_gamm_outlier$PCB_DL,
           plot_adjusted_gamm_outlier$PCB_DL, 
           plot_copollutant_gamm_outlier$PCB_DL,
           nrow = 2)

gaam_results_PCB_NDL <- wrap_plots(plot_base_gamm$PCB_NDL, 
           plot_adjusted_gamm$PCB_NDL, 
           plot_copollutant_gamm$PCB_NDL, 
           plot_base_gamm_outlier$PCB_NDL,
           plot_adjusted_gamm_outlier$PCB_NDL, 
           plot_copollutant_gamm_outlier$PCB_NDL,
           nrow = 2)

gaam_results_PCB_4 <- wrap_plots(plot_base_gamm$PCB_4, 
           plot_adjusted_gamm$PCB_4, 
           plot_base_gamm_outlier$PCB_4,
           plot_adjusted_gamm_outlier$PCB_4,
           nrow = 2)

gaam_results_OCP_HCB <- wrap_plots(plot_base_gamm$OCP_HCB, 
           plot_adjusted_gamm$OCP_HCB, 
           plot_copollutant_gamm$OCP_HCB, 
           plot_base_gamm_outlier$OCP_HCB,
           plot_adjusted_gamm_outlier$OCP_HCB, 
           plot_copollutant_gamm_outlier$OCP_HCB,
           nrow = 2)

gaam_results_ΣDDT <- wrap_plots(plot_base_gamm$ΣDDT, 
           plot_adjusted_gamm$ΣDDT, 
           plot_copollutant_gamm$ΣDDT, 
           plot_base_gamm_outlier$ΣDDT,
           plot_adjusted_gamm_outlier$ΣDDT, 
           plot_copollutant_gamm_outlier$ΣDDT,
           nrow = 2)

gaam_results_OCP_β_HCH <- wrap_plots(plot_base_gamm$OCP_β_HCH, 
           plot_adjusted_gamm$OCP_β_HCH, 
           plot_copollutant_gamm$OCP_β_HCH, 
           plot_base_gamm_outlier$OCP_β_HCH,
           plot_adjusted_gamm_outlier$OCP_β_HCH, 
           plot_copollutant_gamm_outlier$OCP_β_HCH,
           nrow = 2)

gaam_results_Σchlordane <- wrap_plots(plot_base_gamm$Σchlordane, 
           plot_adjusted_gamm$Σchlordane, 
           plot_copollutant_gamm$Σchlordane, 
           plot_base_gamm_outlier$Σchlordane,
           plot_adjusted_gamm_outlier$Σchlordane, 
           plot_copollutant_gamm_outlier$Σchlordane,
           nrow = 2)

gaam_results_ΣPBDE <- wrap_plots(plot_base_gamm$ΣPBDE, 
           plot_adjusted_gamm$ΣPBDE, 
           plot_copollutant_gamm$ΣPBDE, 
           plot_base_gamm_outlier$ΣPBDE,
           plot_adjusted_gamm_outlier$ΣPBDE, 
           plot_copollutant_gamm_outlier$ΣPBDE,
           nrow = 2)


## metaanalysis (quart) ----
run_clogit <- function(formula, data) {
  model <- clogit(formula, data = data)
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  tibble(
    term = rownames(coefs),
    coef = coefs[, "coef"],
    se = coefs[, "se(coef)"])
}
POPs_group_metaanalysis <- setdiff(POPs_group, "ΣPBDE")                         # we don't include ΣPBDE in the metaanalysis because low levels 
POPs_group_metaanalysis_quart <- paste0(POPs_group_metaanalysis, "_quart")

### Base model ----
metaanalysis_base_quart <- map_dfr(POPs_group_metaanalysis_quart, function(expl) {
  formula <- as.formula(paste("als ~", expl, "+ strata(match)"))                # base formula: matched, not ajstuded
  
  bdd_danish <- bdd |>                                                     
    filter(study == "Danish") |>                                                # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  bdd_finnish_MFH <- bdd |> 
    filter(study == "MFH") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific    
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
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "base") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})

### Adjusted model ----
metaanalysis_adjusted_quart <- map_dfr(POPs_group_metaanalysis_quart, function(expl) {
  formula_educ <- 
    as.formula(paste("als ~", expl, 
                     "+ strata(match) + marital_status_2cat + smoking_2cat + bmi + cholesterol + education"))
  formula_no_educ <- 
    as.formula(paste("als ~", expl,                                             # education mot available in one finnish cohort
                     "+ strata(match) + marital_status_2cat + smoking_2cat + bmi + cholesterol")) 
  
  bdd_danish <- bdd |>                                                     
    filter(study == "Danish") |>                                                # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  
  bdd_finnish_FMC <- bdd |>                                                     
    filter(study == "FMC") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific                      
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) 
  bdd_finnish_FMCF <- bdd |> 
    filter(study == "FMCF") |>                                                  # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  bdd_finnish_MFH <- bdd |> 
    filter(study == "MFH") |>                                                   # creation of one dataset per finnish cohort
    mutate(across(all_of(POPs_group_metaanalysis), ~ factor(ntile(.x, 4),       # creation of quartiles cohort specific    
                                                            labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart"))
  
  results <- list(                                                              # run of the simple conditional logistic regression
    danish = run_clogit(formula_educ, bdd_danish),
    finnish_FMC = run_clogit(formula_no_educ, bdd_finnish_FMC),                 # no education data for this cohort 
    finnish_FMCF = run_clogit(formula_educ, bdd_finnish_FMCF), 
    finnish_MFH = run_clogit(formula_educ, bdd_finnish_MFH)) |>
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
        `p-value` = as.numeric(rma_fit$pval))
    }) |> 
    ungroup() |> 
    mutate(model = "adjusted") |> 
    relocate(model, explanatory, term)
  return(meta_results)
})


metaanalysis_quart <- bind_rows(metaanalysis_base_quart, metaanalysis_adjusted_quart) |> 
  mutate(explanatory = gsub("_quart", "", explanatory), 
         OR = as.numeric(sprintf("%.1f", OR)),
         lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
         upper_CI = as.numeric(sprintf("%.1f", upper_CI)),
         `p-value_raw` = `p-value`, 
         `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = ".")), 
         `p-value` = ifelse(`p-value` == "1.00", ">0.99", `p-value`), 
         "95%CI" = paste(lower_CI, ", ", upper_CI, sep = '')) |>
  select(model,
         explanatory, 
         term,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p-value"), 
         lower_CI, upper_CI) 

rm(metaanalysis_base_quart, metaanalysis_adjusted_quart, run_clogit, POPs_group_metaanalysis, POPs_group_metaanalysis_quart)

plot_metaanalysis_quart <- metaanalysis_quart %>% 
  mutate(`p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted",
                            "Base model" = "base"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         term = fct_relevel(term, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         explanatory = fct_recode(explanatory, 
                                  "Most\nprevalent\nPCBs" = "PCB_4",
                                  "Dioxin-like\nPCBs" = "PCB_DL",
                                  "Non-dioxin-\nlike PCBs" = "PCB_NDL",
                                  "β-HCH" = "OCP_β_HCH", 
                                  "HCB" = "OCP_HCB"), 
         explanatory = fct_relevel(explanatory, 
                                   "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "Most\nprevalent\nPCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane")) %>%
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

## plot metanalysis + danish ----
results_danish <- 
  results_POPs_ALS_occurrence$main$main_results |> 
  filter(model %in% c("base_quart", "adjusted_quart")) |> 
  select(explanatory = variable, model, term = df, OR, "95%CI", "p-value" = p.value, 
         "p-value_raw" = "p.value_raw", "lower_CI", "upper_CI") |>
  mutate(model = str_replace(model, "_quart", ""), 
         term = fct_recode(term, 
                           "Quartile 2" = "2",
                           "Quartile 3" = "3",
                           "Quartile 4" = "4"), 
         analysis = "Danish")
metaanalysis_quart <- metaanalysis_quart |>
  mutate(analysis = "Metanalysis")

total <- bind_rows(results_danish, metaanalysis_quart)


plot_quart_totale <- total %>% 
  filter(!explanatory == "ΣPBDE") |>
  mutate(`p-value_shape` = ifelse(`p-value_raw`<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted",
                            "Base model" = "base"),
         model = fct_relevel(model, 'Base model', 'Adjusted model'), 
         term = fct_relevel(term, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         explanatory = fct_recode(explanatory, 
                                  "Most\nprevalent\nPCBs" = "PCB_4",
                                  "Dioxin-like\nPCBs" = "PCB_DL",
                                  "Non-dioxin-\nlike PCBs" = "PCB_NDL",
                                  "β-HCH" = "OCP_β_HCH", 
                                  "HCB" = "OCP_HCB"), 
         explanatory = fct_relevel(explanatory, 
                                   "Dioxin-like\nPCBs", "Non-dioxin-\nlike PCBs", "Most\nprevalent\nPCBs", "HCB", "ΣDDT", "β-HCH", "Σchlordane")) |>
  ggplot(aes(x = term, y = OR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`, shape = analysis)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.5) +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(explanatory), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value", shape = "Analysis") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()



## export ----
ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_PCB_DL.tiff",
  gaam_results_PCB_DL,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_PCB_NDL.tiff",
  gaam_results_PCB_NDL,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_PCB_4.tiff",
  gaam_results_PCB_4,
  height = 5,
  width = 7,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_OCP_HCB.tiff",
  gaam_results_OCP_HCB,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_ΣDDT.tiff",
  gaam_results_ΣDDT,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_OCP_β_HCH.tiff",
  gaam_results_OCP_β_HCH,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_Σchlordane.tiff",
  gaam_results_Σchlordane,
  height = 8,
  width = 10,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/gaam_results_ΣPBDE.tiff",
  gaam_results_ΣPBDE,
  height = 8,
  width = 10,
  units = "in")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/plot_metaanalysis_quart.tiff",
  plot_metaanalysis_quart,
  height = 8,
  width = 6,
  units = "in")


ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/plot_quart_totale.tiff",
  plot_quart_totale,
  height = 8,
  width = 8,
  units = "in")



ggsave(
  "~/Documents/POP_ALS_2025_02_03/team meeting presentation 04.28.2025/POPs_boxplot_comp.tiff",
  POPs_boxplot_comp,
  height = 5,
  width = 8,
  units = "in")

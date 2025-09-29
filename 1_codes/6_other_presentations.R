# Harvard EH department retreat 2025 ----
# Figures for poster presentation Harvard EH department retreat 2025
# Aline Davias
# September, 29, 2025

## data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")

## figure 1 (als risk, mixture GAMs) ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "Hexachlorobenzene", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
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
  # p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", formatC(p_value, format = "f", digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = "", y = "Predicted probability of ALS", title = x_label) +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 5, color = "black") +
    theme_lucid() +
    scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +  
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 14))
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab("") + 
    xlab("") + 
    theme_lucid() +
    theme(axis.text.x = element_text(size = 14), 
          plot.background  = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
})|> set_names(POPs_group_bis)


plot_copollutant_gam_bis <- map(POPs_group_bis, function(var) {
  
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
  # p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", formatC(p_value, format = "f", digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = "", y = "", title = x_label) +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 5, color = "black") +
    theme_lucid() +
    scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +  
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 14))
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab("") + 
    xlab("") + 
    theme_lucid() +
    theme(axis.text.x = element_text(size = 14), 
          plot.background  = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
})|> set_names(POPs_group_bis)

Figure_1 <- wrap_plots(
  plot_copollutant_gam_bis$OCP_HCB, 
  plot_copollutant_gam$PCB_NDL, 
  plot_copollutant_gam_bis$Σchlordane,
  ncol = 1) & 
  theme(
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA))

rm(plot_copollutant_gam, POPs_group_bis, pollutant_labels_bis, model)

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Harvard_EH_retreat_2025/Figure_1.tiff",
  Figure_1,
  height = 7.5,
  width = 6,
  units = "in")

## figure 2 (als survival, single pollutant cox models) ----
Figure_2 <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(term == "Continuous") |>
  filter(model == "adjusted") |>
  mutate(model = fct_recode(model, 
                            "Adjusted model" = "adjusted"),
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = 'right',  
        legend.direction = "vertical",
        strip.text.y = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Harvard_EH_retreat_2025/Figure_2.tiff",
  Figure_2,
  dpi = 300, 
  height = 5,
  width = 7,
  units = "in")

## figure 3 (als survival, mixture cox model) ----
Figure_3 <- 
  results_POPs_ALS_survival$sensi5$POPs_sd_ALS_figure_sensi5_danish +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = 'right',  
        legend.direction = "vertical",
        strip.text.y = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA)) 
    
ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Harvard_EH_retreat_2025/Figure_3.tiff",
  Figure_3,
  dpi = 300, 
  height = 2,
  width = 6.5,
  units = "in")


# ISEE 2025 ----
# Figures for ISEE oral presentation 2025 Atlanta
# Aline Davias
# August, 14, 2025

## data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")

## result slide 1 (als risk) ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "Hexachlorobenzene", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
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
  # p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", formatC(p_value, format = "f", digits = 2))
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = "", y = "Predicted probability\nof ALS", title = x_label) +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 5, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +  
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 12))
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab("") + 
    xlab("") + 
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12))
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
})|> set_names(POPs_group_bis)

results_slide_1 <- wrap_plots(
  plot_copollutant_gam$OCP_HCB, 
  plot_copollutant_gam$PCB_NDL, 
  plot_copollutant_gam$Σchlordane,
  ncol = 1)

rm(plot_copollutant_gam, POPs_group_bis, pollutant_labels_bis, model)

## results slide 2 (als survival) ----
### a -----
results_slide_2_a <- results_POPs_ALS_survival$main_analysis$main_results_POPs_ALS_survival |>
  filter(study == "Danish") |>
  filter(term == "Continuous") |>
  filter(model == "adjusted") |>
  mutate(model = fct_recode(model, 
                            "Adjusted model" = "adjusted"),
         explanatory = factor(explanatory, levels = POPs_group_labels),
         explanatory = fct_rev(explanatory),
         explanatory = fct_recode(explanatory, !!!POPs_group_labels)) |>
  arrange(explanatory) |> 
  ggplot(aes(x = explanatory, y = HR, ymin = lower_CI, ymax = upper_CI, color = `p-value_shape`)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Hazard Ratio (HR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = 'right',  
        legend.direction = "vertical",
        strip.text.y = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) +
  coord_flip()

### b -----
bdd_cases_danish <- bdd_danish |>
  filter (als == 1) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_group)) |>
  mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group),                                             # create cohort and cases specific scaled POPs variables 
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd"))  |>
  replace_with_median(PCB_4, PCB_4_quart) |>
  replace_with_median(PCB_DL, PCB_DL_quart) |>
  replace_with_median(PCB_NDL, PCB_NDL_quart) |>
  replace_with_median(OCP_HCB, OCP_HCB_quart) |>
  replace_with_median(ΣDDT, ΣDDT_quart) |>
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) |>
  replace_with_median(Σchlordane, Σchlordane_quart) |>
  replace_with_median(ΣPBDE, ΣPBDE_quart) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  POPs_group_bis,
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"))

results_slide_2_b <- 
  bdd_cases_danish |> 
  select(all_of(POPs_group_bis)) |>
  rename(!!!pollutant_labels_bis) 
results_slide_2_b <- 
  cor(results_slide_2_b, 
      use = "pairwise.complete.obs", 
      method = "pearson")

rm(bdd_cases_danish, POPs_group_bis, pollutant_labels_bis)

## export ----
ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/ISEE_2025/results_slide_1.tiff",
  results_slide_1,
  dpi = 300, 
  height = 8,
  width = 5,
  units = "in")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/ISEE_2025/results_slide_2_a.tiff",
  results_slide_2_a,
  dpi = 300, 
  height = 5,
  width = 7,
  units = "in")

tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/ISEE_2025/results_slide_2_b.tiff", units = "mm", width = 130, height = 120, res = 300)
corrplot(results_slide_2_b, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         number.digits = 1,
         tl.cex = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()
 
ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/ISEE_2025/results_slide_2_c.tiff",
  results_slide_2_c,
  dpi = 300, 
  height = 3,
  width = 5,
  units = "in")



# TARGET ALS  ----
# Figures for the Target ALS grant application
# Aline Davias
# August, 6, 2025


## data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")

## figure POPs and ALS occurence ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "Hexachlorobenzene", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
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
  # p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall =2, digits = 2))
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(round(p_value, 2), nsmall = 2))
  
  x_min <- min(bdd_danish[[var]], na.rm = TRUE)
  x_label <- paste(pollutant_labels_bis[var], "(pg/ml)")
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability\nof ALS") +
    annotate("text", x = x_min, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 0, vjust = 1.2, size = 3, color = "black", family = "Arial") +
    theme_minimal(base_family = "Arial") +
    scale_y_continuous(limits = c(0, 1)) +  
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal(base_family = "Arial") +
    theme(axis.title.x = element_text(size = 10))
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
})|> set_names(POPs_group_bis)

rm(POPs_group_bis, pollutant_labels_bis, model)


figure_1 <- wrap_plots(
  plot_copollutant_gam$PCB_NDL, 
  plot_copollutant_gam$OCP_HCB, 
  nrow = 2) + 
  plot_annotation(caption = 
                    "Figure 1. Estimated risk of ALS occurrence\nattributed to pre-disease POP plasma\nconcentrations in the Danish Diet, Cancer\nand Health study cohort (generalized\nadditive models, sample size: 166 ALS case\nand 332 controls).", 
                  theme = theme(plot.caption = element_text(hjust = 0, size = 10)))


ggsave(
  "~/Documents/Appels à projets/Target ALS/figure_1.tiff",
  figure_1,
  dpi = 300, 
  height = 5,
  width = 3,
  units = "in")



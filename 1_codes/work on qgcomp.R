
# work on qgcomp models in the danish cohort ----
POPs_group_sd_bis <- setdiff(POPs_group_sd, "PCB_4_sd")                         # remove the 4 most abundant PCB because they are already NDL-PCB
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_sd_bis)

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(covariates_danish, collapse = " + ")))



## unadjusted model ----
set.seed(1996)
fit_qgcomp_boot_base <-
  qgcomp.cox.boot(
    f = Surv(follow_up_death, status_death) ~ .,                                # formula
    bdd_cases_danish[, c(POPs_group_sd_bis, 'follow_up_death', 'status_death')],# no covariates (unadjusted model)
    q = 4,                                                                      # nb of quantiles
    expnms = POPs_group_sd_bis,                                                 # exposures of interest 
    B = 1000,                                                                   # nb of boostrap
    MCsize = 5000,
    seed = 1996,
    parallel = TRUE,                                                            # shorter run time
    parplan = TRUE)                                                             # shorter run time
  
print(fit_qgcomp_boot_base)
fit_qgcomp_boot_base$pos.weights                                                # NULL because model not significant
fit_qgcomp_boot_base$neg.weights                                                # NULL because model not significant
plot(fit_qgcomp_boot_base, suppressprint = TRUE)

p <- summary(fit_qgcomp_boot_base)
results <-                                                                      # overall results
  tibble(
    study = "Danish", 
    model = "copollutant", 
    HR = round(exp(fit_qgcomp_boot_base$psi),1),
    lower_CI = round(exp(fit_qgcomp_boot_base$ci[1]), 1),
    upper_CI = round(exp(fit_qgcomp_boot_base$ci[2]), 1),
    p_value = round(p$coefficients[1, "Pr(>|z|)"], 2))
rm(p, fit_qgcomp_boot_base)

fit_qgcomp_noboot_base <-                                                       # run the code without bootsrapping just to get weights even if the mixture is not significant 
  qgcomp.cox.noboot(
    f = Surv(follow_up_death, status_death) ~ .,                                # formula
    bdd_cases_danish[, c(POPs_group_sd_bis, 'follow_up_death', 'status_death')],
    q = 4,                                                                      # number of quantiles
    expnms = POPs_group_sd_bis)                                                 # exposures of interest
print(fit_qgcomp_noboot_base)
fit_qgcomp_noboot_base$pos.weights
fit_qgcomp_noboot_base$neg.weights
plot(fit_qgcomp_noboot_base, suppressprint = TRUE)


weights_df <- tibble(
  pollutant = c(names(fit_qgcomp_noboot_base$pos.weights), names(fit_qgcomp_noboot_base$neg.weights)),
  weight = c(fit_qgcomp_noboot_base$pos.weights, -fit_qgcomp_noboot_base$neg.weights)) %>%
  mutate(
    pollutant_label = pollutant_labels_bis[pollutant] %||% pollutant,
    pollutant_label = factor(pollutant_label, levels = rev(pollutant_labels_bis))
  )

gg_weights <- ggplot(weights_df, aes(x = weight, y = pollutant_label, fill = weight > 0)) +
  geom_col(show.legend = FALSE) +
  labs(y = "Exposures", x = "                                                     Negative weights                                                               Positive weights") +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_lucid()



## adjusted model ----
formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(covariates_danish, collapse = " + ")))
set.seed(1996)
fit_qgcomp_boot_adj <-
  qgcomp.cox.boot(
    f = formula_danish,                                                         # formula
    data =   bdd_cases_danish, 
    q = 4,                                                                      # nb of quantiles
    expnms = POPs_group_sd_bis,                                                 # exposures of interest 
    B = 1000,                                                                   # nb of boostrap
    MCsize = 5000,
    seed = 1996,
    parallel = TRUE,                                                            # shorter run time
    parplan = TRUE)                                                             # shorter run time
print(fit_qgcomp_boot_adj)
fit_qgcomp_boot_adj$pos.weights
fit_qgcomp_boot_adj$neg.weights
plot(fit_qgcomp_boot_adj)


p <- summary(fit_qgcomp_boot_adj)
results_adj <-                                                                      # overall results
  tibble(
    study = "Danish", 
    model = "copollutant", 
    HR = round(exp(fit_qgcomp_boot_adj$psi),1),
    lower_CI = round(exp(fit_qgcomp_boot_adj$ci[1]), 1),
    upper_CI = round(exp(fit_qgcomp_boot_adj$ci[2]), 1),
    p_value = round(p$coefficients[1, "Pr(>|z|)"], 2))
rm(p, fit_qgcomp_boot_adj)

formula_danish <-                                                               # set the formulas  
  as.formula(paste("Surv(follow_up_death, status_death) ~",   
                   paste(c(covariates_danish, POPs_group_sd_bis), collapse = " + ")))
fit_qgcomp_noboot_adj <-                                                       # run the code without bootsrapping just to get weights even if the mixture is not significant 
  qgcomp.cox.noboot(
    f = formula_danish,                                                         # formula
    bdd_cases_danish[, c(POPs_group_sd_bis, covariates_danish, 'follow_up_death', 'status_death')],
    q = 4,                                                                      # number of quantiles
    expnms = POPs_group_sd_bis)                                                 # exposures of interest
print(fit_qgcomp_noboot_adj)
fit_qgcomp_noboot_adj$pos.weights
fit_qgcomp_noboot_adj$neg.weights
plot(fit_qgcomp_noboot_adj, suppressprint = TRUE)


weights_df <- tibble(
  pollutant = c(names(fit_qgcomp_noboot_adj$pos.weights), names(fit_qgcomp_noboot_adj$neg.weights)),
  weight = c(fit_qgcomp_noboot_adj$pos.weights, - fit_qgcomp_noboot_adj$neg.weights)) %>%
  mutate(
    pollutant_label = pollutant_labels_bis[pollutant] %||% pollutant,
    pollutant_label = factor(pollutant_label, levels = rev(pollutant_labels_bis)))

gg_weights <- ggplot(weights_df, aes(x = weight, y = pollutant_label, fill = weight > 0)) +
  geom_col(show.legend = FALSE) +
  labs(y = "Exposures", x = "                                                     Negative weights                                                               Positive weights") +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_lucid()







# test (old)----
model <- gam(als ~ s(ΣDDT) +
               sex + baseline_age +
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = 'REML', 
             data = bdd_danish)
test <- ggpredict(model, term = 'ΣDDT') |> plot() +
  geom_boxplot(aes(data = bdd_danish, y = ΣDDT)) + plot_layout(ncol = 1)
gam.check(model)

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", 
                "cholesterol_i", "marital_status_2cat_i", "education_i")

bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
  mutate(across(all_of(covariates), 
                ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                .names = "adj_{.col}")) %>%
  select(ΣDDT, starts_with("adj_")) %>%
  rename_with(~ gsub("adj_", "", .x)) 

pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)

bdd_pred <- bdd_pred %>%
  mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
         prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
         prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))


p1 <- ggplot(bdd_pred, aes(x = ΣDDT, y = prob)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
  labs(x = var, y = "Predicted probability of ALS") +
  theme_minimal()
p2 <- ggplot(bdd_pred) +
  aes(x = "", y = ΣDDT) +
  geom_boxplot(fill = "#112446") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(color = 'white'),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(2, 1),
                   guides = 'collect')





model <- gam(als ~ s(HCB_outlier, k = 4) +
               sex + baseline_age +
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = 'REML', 
             data = bdd_danish)

summary(model)
ggpredict(model, term = 'HCB_outlier') |> plot()
gam.check(model)
g


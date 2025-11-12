## sensi1 + sensi4 - Removing the oulier for NEFL + showing the NEFL results depending on follow-up duration  ----


s(follow_up, proteomic )

mgcv

or quadratic : 
  
  poly(time*proteomic )


formula <- as.formula(paste("als ~ s(follow_up, by = proteomic_neuro_explo_NEFL_sd)"))
model <- gam(formula, 
             family = binomial, 
             method = 'REML', 
             data = bdd_danish)
model_summary <- summary(model)


library(ggplot2)

# Choose representative NEFL values (e.g., 25th, 50th, 75th percentile)
nefl_vals <- quantile(bdd_danish$proteomic_neuro_explo_NEFL_sd, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Create prediction grid
newdat <- expand.grid(
  follow_up = seq(min(bdd_danish$follow_up), max(bdd_danish$follow_up), length = 100),
  proteomic_neuro_explo_NEFL_sd = nefl_vals
)

# Predict probabilities
newdat$pred <- predict(model, newdata = newdat, type = "response")

# Plot
ggplot(newdat, aes(x = follow_up, y = pred, color = factor(proteomic_neuro_explo_NEFL_sd))) +
  geom_line(size = 1) +
  labs(x = "Follow-up time",
       y = "Predicted probability of ALS",
       color = "NEFL (sd)") +
  theme_minimal()




# test GAM 12/11/2025 ----- 
library(mgcv)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glue)


bdd_danish_additional_1 <- 
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
    .names = "{.col}_quart_med_sensi_1" )) |>
  group_by(match) %>%
  mutate(
    follow_up_bis = if_else(is.na(follow_up), follow_up[als == 1], follow_up), 
    follow_up_bis = -follow_up_bis/12) |>
  ungroup()


# Modèle GAM
formula <- as.formula(glue("als ~ s(follow_up_bis, by = proteomic_neuro_explo_NEFL_sensi_1) + sex + baseline_age"))
model <- gam(formula, family = binomial, method = "REML", data = bdd_danish_additional_1)

# Plot all smooth terms
plot(model, pages = 1, 
     rug = TRUE,         # add data points
     se = TRUE)          # add CI

draw(model)



# Extract smooth for follow_up_bis by your by-variable
sm <- smooth_estimates(model, smooth = "s(follow_up_bis):proteomic_neuro_explo_NEFL_sensi_1")
head(sm)
ggplot(sm, aes(x = follow_up_bis, y = .estimate, color = proteomic_neuro_explo_NEFL_sensi_1)) +
  geom_line() +
  geom_ribbon(aes(ymin = .estimate - 2*.se, ymax = .estimate + 2*.se), alpha = 0.2) +
  theme_minimal() +
  labs(y = "Effect on ALS log-odds", x = "Follow-up BIS")



# Résumé du modèle
model_summary <- summary(model)
edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1)
p_value <- model_summary$s.table[1, "p-value"]
p_value_text <- ifelse(p_value < 0.01, "< 0.01", format(p_value, nsmall = 2, digits = 2))

# Création d'une BDD pour prédiction avec covariables ramenées à leur moyenne
bdd_pred <- bdd_danish_additional_1 |>
  mutate(
    adj_baseline_age = mean(baseline_age, na.rm = TRUE),
    adj_sex = names(which.max(table(sex)))) |>
  select(proteomic_neuro_explo_NEFL_sensi_1, starts_with("adj_")) |>
  rename_with(~ gsub("adj_", "", .x))

# Prédiction
pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
bdd_pred <- bdd_pred |>
  mutate(
    prob = plogis(pred$fit),
    prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
    prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))

# Labels
x_min <- min(bdd_danish_additional_1$proteomic_neuro_explo_NEFL_sensi_1, na.rm = TRUE)
x_label <- "NEFL"  # si tu veux enlever le préfixe

# Plot 1 : prédiction
p1 <- ggplot(bdd_pred, aes(x = .data$proteomic_neuro_explo_NEFL_sensi_1, y = prob)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
  labs(x = var, y = "Predicted probability of ALS") +
  annotate("text", x = x_min, y = Inf, 
           label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
           hjust = 0, vjust = 1.2, size = 4, color = "black") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(color = 'white'),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Base model excluding NEFL outlier")

# Plot 2 : boxplot
p2 <- ggplot(bdd_pred, aes(x = "", y = .data$proteomic_neuro_explo_NEFL_sensi_1)) +
  geom_boxplot(fill = "blue") +
  coord_flip() +
  ylab(x_label) +
  xlab("") +
  theme_minimal()

# Combinaison
plot_base_gam_sensi_1 <- p1 / p2 + 
  plot_layout(ncol = 1, nrow = 2, heights = c(10, 1), guides = 'collect') +
  theme_minimal()


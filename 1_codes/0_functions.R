library(adelie)
library(broom)
library(car)
library(corrplot)
library(distill)
library(dplyr)
library(esquisse)
library(expss)
library(earth)        # MARS
library(flextable)
library(fasano.franceschini.test)
library(ggpubr)
library(grid)
library(GGally)
library(gtsummary)
library(grDevices)
library(ggplot2)
library(ggstance)
library(gratia)
library(glue)
library(gridExtra)  
library(ggeffects)
library(gWQS)
library(glmnet)       # lasso/ridge/elastic net
library(ggrepel)
library(haven)
library(Hmisc)
library(knitr)
library(kableExtra)
library(knitr)
library(kernlab)      # SVM
library(labelled)
library(lazyeval)
library(mice)
library(mgcv)         # GAM
library(metafor)
library(officer)
library(patchwork)
library(plotly)
library(parameters)
library(purrr)
library(performance)
library(pROC)
library(probably)     # calibration des probabilités
library(questionr)
library(qgcomp)
library(ranger)       # random forest
library(rmarkdown)
library(rmdformats)
library(RColorBrewer)
library(rstatix)
library(readxl)
library(reshape2)
library(randomForest)
library(survival)
library(survminer)
library(splines)
library(see)
library(scales)
library(shapviz)      # SHAP values
library(SuperLearner)
library(sva)
library(tidyverse)
library(tidyr)
library(tidymodels)
library(tune)
library(vip)
library(writexl)
library(xgboost)    # gradient boosting



## Tableau descriptif variables numériques
descrip_num <- function(data, vars) {                                           
  table_1 <- data %>%
    select(all_of({{vars}})) %>%
    tbl_summary(
      missing = "no", 
      type = list(where(is.numeric)~ "continuous"), 
      statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}", 
      digits = all_continuous() ~ 2
    ) %>%
    bold_labels() %>% 
    as_gt()  %>% 
    as.data.frame()%>% 
    select(variable, stat_0) %>%
    separate(
      col = stat_0, 
      into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
      sep = "/", 
      remove = TRUE) %>%
    mutate(N = as.integer(N))
  
  table_2 <- data.frame(
    variable = {{vars}},
    "Zero count" = sapply({{vars}}, function(var) sum(data[[var]] == 0, na.rm = TRUE))
  )
  
  left_join(table_1, table_2, by = "variable")
}

## Tableau comparaison d'effectifs
comp_effectifs <- function(data, vars_col1, vars_col2, name_col1, name_col2){ 
  table_col1 <- data %>% select(all_of({{vars_col1}})) 
  table_col2 <- data %>% select(all_of({{vars_col2}}))
  colnames(table_col2) <- colnames(table_col1)
  
  table_col1 <- table_col1 %>% tbl_summary()
  table_col2 <- table_col2 %>% tbl_summary()
  comp <- tbl_merge(
    tbls = list(table_col1, table_col2), 
    tab_spanner = c({{name_col1}}, {{name_col2}}))
  return(comp)
}

## Vérif distrib variables continues
verif_distrib <- function(data, var) {                                               # data = df large dans lequel récuperer les variables                                                                                         # vars = vecteur de noms de variables (toutes numériques)
  
  boxplot <- data %>%
    ggplot() +
    aes(x = "", y = {{var}}) +
    geom_boxplot(shape = "circle", fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) 
  densityplot <- data %>%
    ggplot() +
    aes(x = {{var}}) +
    geom_density(fill = "lightgray") +
    theme_bw() +
    stat_overlay_normal_density(color = "red", linetype = "dashed") +
    theme(axis.title = element_blank())
  qqnorm <- data %>%
    ggplot(aes(sample = {{var}})) +
    stat_qq() +
    stat_qq_line()+
    theme_bw() + 
    theme(axis.title = element_blank())
  results <- boxplot + densityplot + qqnorm
  
  return(results)
}

## Scatterplots
scatterplot <- function(data, outcome, vars) {
  data_long <-
    data %>%
    select(sample, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric)
  data_long <- data_long %>%
    pivot_longer(cols = -sample,
                 names_to = "Variable",
                 values_to = "Value")
  bdd_outcome <-
    data %>% 
    select(sample, {{outcome}})
  data_long <- 
    left_join(data_long,
              bdd_outcome,
              by = "sample") %>%
    rename(Outcome = {{outcome}})
  
  scatterplot <- data_long %>%
    ggplot() +
    aes(x = Outcome, y = Value) +
    geom_point(
      shape = "circle",
      size = 1.55,
      colour = "#112446") +
    labs(y ="") +
    theme_bw() +
    facet_wrap(~Variable, scales = "free", ncol = 4L)
  return(scatterplot)
  
}

## Boxplots
boxplot <- function(data, vars) {                                               # data = df large dans lequel récuperer les variables                                                                                         # vars = vecteur de noms de variables (toutes numériques)
  data_long <-                                                     
    data %>% 
    select(sample, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -sample, names_to = "Variable", values_to = "Value")
  
  boxplot <- data_long %>%
    ggplot() +
    aes(x = "", y = Value) +
    geom_boxplot(shape = "circle", fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap(~Variable, scales = "free", ncol = 2L, nrow = 5L)
  
  return(boxplot)
}

## Histogrammes
histogram <- function(data, vars, order_var) {                                  # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(sample, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -sample, names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = fct_relevel(Variable, order_var))
  
  histogram <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_histogram(bins = 30L, fill = "#112446") +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = 2L)
  
  return(histogram)
} 

## Barplots 
barplot <- function(data, vars) {                                               # data = df large dans lequel récuperer les variables 
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes catégorielles)
    data %>% 
    select(sample, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], factor) 
  data_long <- data_long %>% 
    pivot_longer(cols = -sample, names_to = "Variable", values_to = "Value")
  
  barplot <- data_long %>%
    filter(!is.na(Value)) %>%
    ggplot() +
    aes(x = Value) +
    geom_bar(fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = 3L) 
  
  return(barplot)
}

## Densityplots 
densityplot <- function(data, vars) {                                           # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(sample, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -sample, names_to = "Variable", values_to = "Value")
  
  densityplot <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_density(fill = "lightgray") +
    theme_bw() +
    facet_wrap( ~ Variable, scales = "free", ncol = 3L, nrow = 3L) +
    stat_overlay_normal_density(color = "red", linetype = "dashed") +
    theme(axis.title = element_blank())
  
  return(densityplot)
} 


## Heatmap de correlation 
heatmap_cor <- function(cormat, decimal) {                  # cormat = df avec seulement les variables à mettre dans l'heatmap + des noms raccourcis
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){              # Utiliser la corrélation entre les variables comme mesure de distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cormat <- round(cor(cormat, 
                      use = "pairwise.complete.obs", 
                      method = "pearson"), decimal)
  cormat <- reorder_cormat(cormat)                 # réordonner les coef de cor
  upper_tri <- get_upper_tri(cormat)               # obtenir que le triangle sup
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)   # passer en df long rapidement 
  
  heatmap <-                                       # heatmap
    ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value),
              color = "black",
              size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 7,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    ))
  
  return(heatmap)
}

heatmap_cor_pairwise <- function(data, vars_1, vars_2, decimal){
  
  
  bdd_cormat <- data %>% select(all_of({{vars_1}}), all_of({{vars_2}}))
  
  cormat <- round(cor(bdd_cormat, 
                      use = "pairwise.complete.obs", 
                      method = "spearman"), decimal)
  
  cormat <- cormat %>% 
    as.data.frame() %>% 
    select(all_of({{vars_1}})) %>% 
    t() %>%
    as.data.frame() %>%
    select(all_of({{vars_2}})) %>%
    as.matrix()
  
  cormat_long <- 
    reshape2::melt(cormat, na.rm = TRUE)  # passer en df long rapidement 
  
  heatmap <-                                       # faire la heatmap
    ggplot(cormat_long, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Spearman\nCorrelation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), 
      axis.text.y = element_text(size = 12)) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value),
              color = "black",
              size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1.9, 0.5),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 5,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    ))
  return(heatmap)
}

# Fonction pour adapter les décimales des p-values
custom_pvalue_fun <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) {
      return(NA) # Retourner NA si p est NA
    } else if (p>0.99) {
      return(">0.99")
    } else if (p < 0.001) {
      # Pour p < 0.001, utiliser la notation scientifique pour afficher toutes les décimales
      return(format(p, scientific = TRUE))
    } else if (p >= 0.001 & p < 0.01) {
      # Pour 0.001 <= p < 0.01, afficher avec 3 décimales
      return(sprintf("%.3f", p))
    } else {
      # Pour p >= 0.01, afficher avec 2 décimales
      return(sprintf("%.2f", p))
    }
  })
}


replace_with_median <- function(data, var, quartile_var) {                       
  new_var_name <- paste0(rlang::as_name(ensym(var)), "_quart_med")  
  data %>%
    group_by({{quartile_var}}) %>%
    mutate(!!new_var_name := median({{var}}, na.rm = TRUE)) %>%
    ungroup()
}


make_gam_plot_base <- function(var, data) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age"))
  
  model <- mgcv::gam(formula, family = binomial, method = "REML", data = data)
  
  bdd_pred <- data |>
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x))
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      prob = plogis(pred$fit),
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01",
                         format(p_value, nsmall = 2, digits = 2))
  
  x_min <- min(data[[var]], na.rm = TRUE)
  x_label <- all_vars_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "steelblue", alpha = 0.2) +
    labs(y = "Predicted probability of ALS") +
    annotate(
      "text", x = x_min, y = Inf,
      label = paste("EDF:", edf, "\np-value:", p_value_text),
      hjust = 0, vjust = 1.2, size = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()) +
    ggtitle("Base model")
  
  p2 <- ggplot(bdd_pred, aes(x = "", y = .data[[var]])) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab(x_label) +
    xlab("") +
    theme_minimal()
  
  p <- p1 / p2 +
    plot_layout(heights = c(10, 1), guides = "collect")
  
  list(
    plot = p,
    p_value = p_value)
}


make_gam_plot_adjusted <- function(var, data) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + sex + baseline_age  + smoking_2cat_i + bmi"))
  
  model <- mgcv::gam(formula, family = binomial, method = "REML", data = data)
  
  bdd_pred <- data |>
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_sex = names(which.max(table(sex))), 
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x))
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      prob = plogis(pred$fit),
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01",
                         format(p_value, nsmall = 2, digits = 2))
  
  x_min <- min(data[[var]], na.rm = TRUE)
  x_label <- all_vars_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "steelblue", alpha = 0.2) +
    labs(y = "Predicted probability of ALS") +
    annotate(
      "text", x = x_min, y = Inf,
      label = paste("EDF:", edf, "\np-value:", p_value_text),
      hjust = 0, vjust = 1.2, size = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
  
  p2 <- ggplot(bdd_pred, aes(x = "", y = .data[[var]])) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab(x_label) +
    xlab("") +
    theme_minimal()
  
  p <- p1 / p2 +
    plot_layout(heights = c(10, 1), guides = "collect")
  
  list(
    plot = p,
    p_value = p_value)
}

make_gam_plot_base_sex <- function(var, data) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + baseline_age"))
  
  model <- mgcv::gam(formula, family = binomial, method = "REML", data = data)
  
  bdd_pred <- data |>
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE)) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x))
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      prob = plogis(pred$fit),
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01",
                         format(p_value, nsmall = 2, digits = 2))
  
  x_min <- min(data[[var]], na.rm = TRUE)
  x_label <- all_vars_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "steelblue", alpha = 0.2) +
    labs(y = "Predicted probability of ALS") +
    annotate(
      "text", x = x_min, y = Inf,
      label = paste("EDF:", edf, "\np-value:", p_value_text),
      hjust = 0, vjust = 1.2, size = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()) +
    ggtitle("Base model")
  
  p2 <- ggplot(bdd_pred, aes(x = "", y = .data[[var]])) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab(x_label) +
    xlab("") +
    theme_minimal()
  
  p <- p1 / p2 +
    plot_layout(heights = c(10, 1), guides = "collect")
  
  list(
    plot = p,
    p_value = p_value)
}

make_gam_plot_adjusted_sex <- function(var, data) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + baseline_age  + smoking_2cat_i + bmi"))
  
  model <- mgcv::gam(formula, family = binomial, method = "REML", data = data)
  
  bdd_pred <- data |>
    mutate(
      adj_baseline_age = mean(baseline_age, na.rm = TRUE),
      adj_bmi = mean(bmi, na.rm = TRUE),
      adj_smoking_2cat_i = names(which.max(table(smoking_2cat_i)))) |>
    select(all_of(var), starts_with("adj_")) |>
    rename_with(~ gsub("adj_", "", .x))
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred |>
    mutate(
      prob = plogis(pred$fit),
      prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
      prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1)
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01",
                         format(p_value, nsmall = 2, digits = 2))
  
  x_min <- min(data[[var]], na.rm = TRUE)
  x_label <- all_vars_labels[var]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper),
                fill = "steelblue", alpha = 0.2) +
    labs(y = "Predicted probability of ALS") +
    annotate(
      "text", x = x_min, y = Inf,
      label = paste("EDF:", edf, "\np-value:", p_value_text),
      hjust = 0, vjust = 1.2, size = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()) +
    ggtitle("Adjusted model")
  
  p2 <- ggplot(bdd_pred, aes(x = "", y = .data[[var]])) +
    geom_boxplot(fill = "steelblue") +
    coord_flip() +
    ylab(x_label) +
    xlab("") +
    theme_minimal()
  
  p <- p1 / p2 +
    plot_layout(heights = c(10, 1), guides = "collect")
  
  list(
    plot = p,
    p_value = p_value)
}


fit_and_plot_cox_gam <- function(
    data,
    vars = proteomic,
    time = "follow_up",
    status = "als", 
    covariates = c("sex", "baseline_age"),  # default
    vars_labels = NULL,
    y_limits = c(-12, 12),
    title = "Base model"                     # default
) {
  
  # 1. Gestion des labels
  if (is.null(vars_labels)) {
    vars_labels <- set_names(vars, vars)
  }
  
  # Sous-fonction pour traiter une protéine
  fit_one <- function(var) {
    
    outcome <- cbind(data[[time]], data[[status]])
    formula_str <- paste0(
      "outcome ~ s(", var, ") + ",
      paste(covariates, collapse = " + "))
    
    model <- mgcv::gam(as.formula(formula_str),
                       data = data,
                       family = cox.ph())
    
    smry <- summary(model)
    edf_val <- smry$s.table[1, "edf"]
    pval_raw <- smry$s.table[1, "p-value"]
    
    edf_label <- format(edf_val, nsmall = 1, digits = 1)
    pval_label <- case_when(
      pval_raw < 0.01 ~ "< 0.01",
      pval_raw > 0.99 ~ "> 0.99",
      TRUE ~ format(pval_raw, nsmall = 2, digits = 2))
    
    pdf(NULL)
    plot_data_raw <- plot(model, select = 1, seWithMean = TRUE, rug = FALSE)
    dev.off()
    
    if (length(plot_data_raw) == 0) return(NULL)
    
    plot_data <- plot_data_raw[[1]]
    smooth_df <- data.frame(
      x = plot_data$x,
      fit = plot_data$fit,
      se = plot_data$se)
    
    p1 <- ggplot(smooth_df, aes(x = x, y = fit)) +
      geom_line(linewidth = 1.2, color = "steelblue") +
      geom_ribbon(aes(ymin = fit - 2 * se, 
                      ymax = fit + 2 * se),
                  alpha = 0.2, fill = "steelblue") +
      labs(
        title = title,
        y = "LogHR (smooth estimate)",
        x = NULL) +
      annotate(
        "text",
        x = -Inf, y = Inf,
        hjust = -0.1, vjust = 1.5,
        label = paste0("EDF: ", edf_label, "\np-value: ", pval_label),
        size = 4.2) +
      scale_y_continuous(limits = y_limits) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
    
    p2 <- ggplot(data, aes(x = .data[[var]], y = 1)) +
      geom_boxplot(width = 0.3, fill = "steelblue", color = "black") +
      labs(x = vars_labels[[var]]) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
    
    list(
      model = model,
      pval_raw = pval_raw,
      plot = p1 / p2 + plot_layout(heights = c(10, 1)))
  }
  
  results <- map(vars, fit_one)
  names(results) <- vars_labels[vars]
  
  return(results)
}

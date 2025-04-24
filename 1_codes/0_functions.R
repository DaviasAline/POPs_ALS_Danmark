library(car)
library(corrplot)
library(distill)
library(dplyr)
library(esquisse)
library(expss)
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
library(haven)
library(Hmisc)
library(knitr)
library(kableExtra)
library(knitr)
library(labelled)
library(lazyeval)
library(mice)
library(mgcv)
library(officer)
library(patchwork)
library(plotly)
library(questionr)
library(rmarkdown)
library(parameters)
library(rmdformats)
library(RColorBrewer)
library(rstatix)
library(readxl)
library(reshape2)
library(purrr)
library(survival)
library(survminer)
library(splines)
library(writexl)
library(see)
library(scales)
library(tidyverse)
library(tidyr)
library(flextable)
library(metafor)

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
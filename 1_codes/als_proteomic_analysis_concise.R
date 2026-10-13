# Aline Davias
# Analysis of als risk depending on proteomic profile (version condensée + sensibilité outliers)

# Data and package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
results_proteomic_ALS_occurrence_y_to_als <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/2.6.4_results_proteomic_ALS_occurrence_y_to_als.rds") 

# Data preparation ----
# Merge des vecteurs protéomiques et EVs 
proteomic <- c(proteomic, EVs_log10)                                            

# Paramètres outliers 
# outlier si |z| > outlier_k, z étant calculé sur les sets complets (pour un prédicteur donné)
outlier_k <- 4  
models <- list(base     = NULL,
               adjusted = c("smoking_2cat_i", "bmi"))

# Garde les sets appariés dont le CAS vérifie la condition (les contrôles suivent)
keep_cases <- function(data, cond) {
  data |>
    group_by(match) |>
    filter(any(als == 1 & {{ cond }})) |>
    ungroup()
}
q50 <- function(data) quantile(data$follow_up, 0.5, na.rm = TRUE)

d3 <- keep_cases(bdd_danish, follow_up > 60)                                    # sensi 3 : on retire les cas (et leurs contrôles) avec follow-up < 60 mois

datasets <- list(
  main           = bdd_danish,                                                  # variables non standardisées
  sensi_2        = keep_cases(bdd_danish, follow_up < 60),                      # sensi 2 : follow-up < 60 mois
  sensi_3        = d3,
  sensi_3_4      = keep_cases(d3, follow_up <= q50(d3)),                        # sensi 4 : 50% des follow-up les plus courts (après sensi 3)
  sensi_3_5      = keep_cases(d3, follow_up >  q50(d3)),                        # sensi 5 : 50% des follow-up les plus longs (après sensi 3)
  sensi_7_female = keep_cases(bdd_danish, sex == "Female"),                     # sensi 7
  sensi_7_male   = keep_cases(bdd_danish, sex == "Male"))

# Sensibilité EV : normalisation des marqueurs EV (valeurs brutes EVs_raw) 
# Marqueur X / somme des dénominateurs -> log10 (ici). La standardisation est faite dans fit_one,
# après retrait des sets avec valeur manquante / non finie pour ce marqueur.
# Un marqueur qui fait partie du dénominateur est exclu de sa propre normalisation.
ev_norms <- list(tetra3      = c("CD9", "CD63", "CD81"),                        # X / (CD9+CD63+CD81)
                 CD9         = "CD9",                                           # X / CD9
                 CD63        = "CD63",                                          # X / CD63
                 CD81        = "CD81",                                          # X / CD81
                 tetra3_CD29 = c("CD9", "CD63", "CD81", "CD29"))                # X / (CD9+CD63+CD81+CD29)

analysis_vars <- map(datasets, ~proteomic)                                      # variables testées dans chaque analyse

for (nm in names(ev_norms)) {
  den_cols <- paste0("ev_", ev_norms[[nm]], "_raw")
  markers  <- setdiff(EVs_raw, den_cols)
  new_vars <- paste0(str_remove(markers, "_raw$"), "_norm_", nm)
  den      <- rowSums(bdd_danish[den_cols])                                     # NA si un des dénominateurs est manquant
  
  if (any(bdd_danish[c(markers, den_cols)] <= 0, na.rm = TRUE))
    warning("Valeurs brutes <= 0 pour ", nm, " : ratio non fini, les sets concernés seront retirés pour ces marqueurs")
  
  a <- paste0("sensi_ev_", nm)
  datasets[[a]] <- bind_cols(bdd_danish,
                             map(markers, ~suppressWarnings(log10(bdd_danish[[.x]] / den))) |>
                               set_names(new_vars) |> as_tibble())
  analysis_vars[[a]] <- new_vars
}
rm(nm, den_cols, markers, new_vars, den, a)

n_analysis <- map_int(datasets, nrow)                                           # N par analyse (avant exclusions propres à chaque prédicteur)
cut_years  <- round(unname(q50(d3)) / 12, 1)                                    # coupure du follow-up (années) pour sensi_3_4 / sensi_3_5


# Effects of the covariates on ALS risk ----
covar <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
      select(als, sex, baseline_age, smoking_2cat_i, bmi, fS_Kol, marital_status_2cat_i, education_i) |>
      tbl_uvregression(y = als, method = glm, method.args = list(family = binomial),
                       exponentiate = TRUE,
                       estimate_fun = label_number(accuracy = 0.1, decimal.mark = "."),
                       pvalue_fun = custom_pvalue_fun) |>
      bold_labels(),
    tbl_2 = clogit(als ~ sex + baseline_age + smoking_2cat_i + bmi + fS_Kol +
                     marital_status_2cat_i + education_i, data = bdd_danish) |>
      tbl_regression(exponentiate = TRUE,
                     estimate_fun = label_number(accuracy = .1, decimal.mark = "."),
                     pvalue_fun = custom_pvalue_fun) |>
      bold_labels()),
  tab_spanner = c("**Univariate**", "**Adjusted**"))


# Preparation des fonctions ----
# Fonction d'ajustement : 1 protéine, 1 dataset, 1 modèle 
# Pour CE prédicteur :
#  1) on retire TOUJOURS les sets appariés (cas + contrôles) contenant au moins une valeur
#     manquante ou non finie (cas OU contrôle) ;
#  2) scale() est calculé sur les sets restants (OR = par SD de l'échantillon analysé) ;
#  3) si outlier = TRUE : on retire en plus les sets contenant au moins un outlier (|z| > k),
#     sans re-standardiser (l'échelle reste celle de l'analyse sans outliers).
prep_data <- function(data, var, outlier = FALSE, k = outlier_k) {
  na_sets <- unique(data$match[!is.finite(data[[var]])])
  data <- filter(data, !match %in% na_sets)
  data$x_sd <- as.numeric(scale(data[[var]]))                                   # prédicteur standardisé
  
  out_sets <- if (outlier) unique(data$match[abs(data$x_sd) > k])
  list(data = filter(data, !match %in% out_sets), na_sets = na_sets, out_sets = out_sets)
}

fit_one <- function(var, analysis, model, outlier, k = outlier_k) {
  prep <- prep_data(datasets[[analysis]], var, outlier, k)
  
  f <- reformulate(c("x_sd", models[[model]], "strata(match)"), response = "als")
  
  tryCatch(
    clogit(f, data = prep$data) |>
      tidy(exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "x_sd") |>
      transmute(explanatory = var, term = "Continuous",
                OR = estimate, lower_CI = conf.low, upper_CI = conf.high, p_value = p.value,
                model = model, analysis = paste0(analysis, if (outlier) "_outlier"),
                n_sets_na = length(prep$na_sets), n_sets_outlier = length(prep$out_sets)),
    error = function(e) { message("Echec : ", var, " / ", analysis, " / ", model, " : ", conditionMessage(e)); NULL })
}


# Analyses ----
#Toutes les analyses (principale + sensibilités, avec et sans outliers) 
grid <- tibble(analysis = names(analysis_vars), var = analysis_vars) |>         # chaque analyse a ses propres variables
  unnest(var) |>
  cross_join(expand_grid(model = names(models), outlier = c(FALSE, TRUE)))
# -> pour limiter les outliers à l'analyse principale : filter(grid, !outlier | analysis == "main")

results <- pmap_dfr(grid, fit_one)

# Mise en forme des résultats 
main_results <- results |>
  rename(OR_raw = OR, p_value_raw = p_value) |>
  group_by(analysis, model, is_ev = str_detect(explanatory, "^ev_")) |>
  mutate(p_value_fdr = p.adjust(p_value_raw, method = "fdr")) |>               # FDR par analyse x modèle x (protéines vs EVs)
  ungroup() |>
  select(-is_ev) |>
  mutate(OR = number(OR_raw, accuracy = 0.01, decimal.mark = "."),
         "95% CI" = paste0(number(lower_CI, accuracy = 0.01, decimal.mark = "."), ", ",
                           number(upper_CI, accuracy = 0.01, decimal.mark = ".")),
         p_value = ifelse(p_value_raw < 0.01, "<0.01", number(p_value_raw, accuracy = 0.01, decimal.mark = ".")),
         protein_group = case_when(str_detect(explanatory, "proteomic_immun_res") ~ "Immune response",
                                   str_detect(explanatory, "proteomic_metabolism") ~ "Metabolism",
                                   str_detect(explanatory, "proteomic_neuro_explo") ~ "Neuro-exploratory",
                                   str_detect(explanatory, "ev_") ~ "Extracellular Vesicles"),
         explanatory = explanatory |>
           str_remove("_log10$") |>                                             # EV principaux : même nom que les EV normalisés
           str_remove("_norm_.*$") |>                                           # EV normalisés : le dénominateur est dans `analysis`
           str_remove("^(proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_|ev_)")) |>
  arrange(explanatory) |>
  select(analysis, model, explanatory, protein_group, term,
         starts_with("OR"), starts_with("95%"), starts_with("p_value"),
         lower_CI, upper_CI, starts_with("n_sets"))

rm(grid, results)


# Figures and Tables ----

# Paramètres  
std_analyses <- c("main", "sensi_2", "sensi_3", "sensi_3_4", "sensi_3_5")
ev_analyses  <- paste0("sensi_ev_", names(ev_norms))                            # EV normalisés (script d'analyses)
all_analyses <- c(std_analyses, ev_analyses, "sensi_7_female", "sensi_7_male")
with_outlier <- function(x) c(x, paste0(x, "_outlier"))
suffixes     <- c(standard = "", outlier = "_outlier")
models_lab   <- c(base = "Base model", adjusted = "Adjusted model")

footer_std <- c("1All models are matched on birth year and sex. Adjusted models further account for smoking and body mass index.",
                "2Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.",
                "3CI: Confidence interval.")
footer_sex <- c("1All models are matched on birth year and adjusted for smoking and body mass index.", footer_std[-1])
footer_for <- function(a) {                                                     # footer adapté à l'analyse (outliers, normalisation EV)
  extra <- c(
    if (str_starts(a, "sensi_ev_"))
      sprintf("EV marker intensities were divided by %s before log10 transformation and standardization.",
              paste(ev_norms[[a |> str_remove("^sensi_ev_") |> str_remove("_outlier$")]], collapse = " + ")),
    if (str_ends(a, "_outlier"))
      sprintf("For each marker, matched sets (cases and their controls) with at least one participant beyond %s SD from the mean were excluded.", outlier_k))
  c(footer_std, paste0(length(footer_std) + seq_along(extra), extra))
}


# Fonction : tableau flextable 
# data   : main_results déjà filtré (term, model et/ou analysis)
# by     : colonne qui définit les blocs de colonnes ("model" ou "analysis")
# groups : vecteur nommé c(valeur_de_by = "Libellé du bloc")
make_ft <- function(data, by, groups, footer = footer_std) {
  stats <- c("OR", "95% CI", "p_value")
  keys  <- c("protein_group", "explanatory", as.vector(outer(stats, names(groups), paste, sep = "_")))
  
  wide <- data |>
    group_by(explanatory) |>
    filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                            # protéines significatives (dans au moins un bloc)
    ungroup()
  
  if (nrow(wide) == 0) return(NULL)                                             # aucune protéine significative : pas de tableau
  
  wide <- wide |>
    select(all_of(by), protein_group, explanatory, all_of(stats)) |>
    pivot_wider(names_from = all_of(by), values_from = all_of(stats)) |>
    arrange(protein_group, explanatory) |>
    select(all_of(keys))
  
  header <- data.frame(
    col_keys = keys,
    top      = c("Protein group", "Pre-disease plasma proteins", rep(unname(groups), each = 3)),
    bottom   = c("Protein group", "Pre-disease plasma proteins", rep(c("OR", "95% CI", "p-value"), length(groups))))
  
  flextable(wide) |>
    set_header_df(mapping = header, key = "col_keys") |>
    add_footer_lines(footer) |>
    theme_vanilla() |>
    merge_h(part = "header") |>
    merge_v(part = "header") |>
    merge_v(j = "protein_group") |>
    align(align = "center", part = "all") |>
    align(j = c("protein_group", "explanatory"), align = "left", part = "all") |>
    bold(j = c("protein_group", "explanatory"), part = "body") |>
    flextable::font(fontname = "Calibri", part = "all") |>
    fontsize(size = 10, part = "all") |>
    padding(padding.top = 0, padding.bottom = 0, part = "all")
}


# Fonction : volcano plot 
volcano <- function(data, title) {
  p <- data |>
    mutate(neg_log10_p = -log10(p_value_raw),
           significance = case_when(p_value_raw < 0.05 & OR_raw > 1 ~ "OR>1 & p-value<0.05",
                                    p_value_raw < 0.05 & OR_raw < 1 ~ "OR<1 & p-value<0.05",
                                    TRUE ~ "p-value>0.05")) |>
    ggplot(aes(x = OR_raw, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_text_repel(data = ~filter(.x, p_value_raw < 0.05), aes(label = explanatory),
                    size = 3.5, max.overlaps = 20, box.padding = 0.4, point.padding = 0.2,
                    segment.color = "grey20", color = "black") +
    scale_color_manual(values = c("OR<1 & p-value<0.05" = "blue",
                                  "p-value>0.05" = "grey70",
                                  "OR>1 & p-value<0.05" = "red")) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(limits = c(0, or_max), breaks = seq(0, or_max, by = 0.5)) +   # bornes calculées sur toutes les analyses :
    scale_y_continuous(limits = c(0, p_max),  breaks = seq(0, p_max, by = 1)) +      # aucun point tronqué, graphiques comparables
    theme_minimal(base_size = 14) +
    labs(title = title, x = "OR", y = "-log10(p-value)", color = "")
  
  if (n_distinct(data$analysis) > 1) p <- p + facet_grid(cols = vars(analysis))  # ex : femmes / hommes
  p
}


# Données de base 
res <- main_results |> filter(term == "Continuous")

# Bornes des axes des volcano plots, communes à toutes les figures (OR max et -log10(p) max observés)
or_max <- ceiling(max(res$OR_raw, na.rm = TRUE) * 2) / 2
p_max  <- max(-log10(res$p_value_raw[res$p_value_raw > 0]), na.rm = TRUE) |> (\(x) ceiling(x * 10) / 10)()


## Tableaux : base vs ajusté, pour chaque analyse (avec et sans outliers) ----
# -> tables_std$main, tables_std$sensi_2, tables_std$sensi_2_outlier, ...
tables_std <- map(set_names(with_outlier(c(std_analyses, ev_analyses))), \(a)
                  make_ft(filter(res, analysis == a), by = "model", groups = models_lab, footer = footer_for(a)))


## Tableaux : comparaison AVEC vs SANS outliers ----
# -> tables_outlier$main$base, tables_outlier$main$adjusted, tables_outlier$sensi_3$adjusted, ...
tables_outlier <- map(set_names(all_analyses), \(a)
                      map(set_names(names(models_lab)), \(m)
                          make_ft(filter(res, model == m, analysis %in% c(a, paste0(a, "_outlier"))),
                                  by = "analysis",
                                  groups = set_names(c("All participants", "Outliers excluded"), c(a, paste0(a, "_outlier"))),
                                  footer = footer_for(paste0(a, "_outlier")))))


## Tableaux : sensi 7 (femmes / hommes, modèle ajusté) ----
# -> tables_sex$standard, tables_sex$outlier
sex_groups <- function(sfx) set_names(c("Females and males (main analyses)", "Females", "Males"),
                                      paste0(c("main", "sensi_7_female", "sensi_7_male"), sfx))
tables_sex <- map(suffixes, \(sfx)
                  make_ft(filter(res, model == "adjusted", analysis %in% names(sex_groups(sfx))),
                          by = "analysis", groups = sex_groups(sfx), footer = footer_sex))


## Figures : volcano plots ----
# -> figs$main$base, figs$main$adjusted, figs$sensi_2_outlier$base, ...
figs <- map(set_names(with_outlier(c(std_analyses, ev_analyses))), \(a)
            map(set_names(names(models_lab)), \(m)
                volcano(filter(res, analysis == a, model == m), paste(str_to_title(m), "logistic models"))))

# -> figs_sex$standard$base, figs_sex$outlier$adjusted, ... (facettes femmes / hommes)
figs_sex <- map(suffixes, \(sfx)
                map(set_names(names(models_lab)), \(m)
                    volcano(filter(res, model == m, analysis %in% paste0(c("sensi_7_female", "sensi_7_male"), sfx)),
                            paste(str_to_title(m), "logistic models"))))


## Figures : forest plot (femmes vs hommes, modèle ajusté) ----
forest_sex <- function(sfx = "") {
  res |>
    filter(model == "adjusted", analysis %in% paste0(c("sensi_7_female", "sensi_7_male"), sfx)) |>
    group_by(explanatory) |>
    filter(any(p_value_raw < 0.05)) |>
    ungroup() |>
    mutate(analysis = if_else(str_detect(analysis, "female"), "Females", "Males")) |>
    group_by(explanatory) |>
    mutate(or_male = first(OR_raw[analysis == "Males"])) |>
    ungroup() |>
    mutate(explanatory = fct_reorder(explanatory, or_male)) |>
    ggplot(aes(x = explanatory, y = OR_raw, ymin = lower_CI, ymax = upper_CI, color = analysis)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.4) +
    scale_color_manual(values = c("Females" = "gray70", "Males" = "black")) +
    labs(x = "Proteins", y = "Odds Ratio", color = "") +
    theme_lucid() +
    coord_flip()
}
forest_sex_main    <- forest_sex()
forest_sex_outlier <- forest_sex("_outlier")


### Table covariates - als risk ----
covar

rm(res, datasets, d3)


## Figures : différence cas - moyenne des témoins appariés en fonction du temps avant diagnostic ----
# Variables en échelle log (log2 pour les protéines, log10 pour les EVs) : la différence = log fold change
# Nécessite (script d'analyses) : main_results, proteomic (incluant EVs_log10), bdd_danish, prep_data()

# Nom "propre" d'une variable : même nettoyage que dans main_results
clean_name <- function(x) {
  x |>
    str_remove("_log10$") |>
    str_remove("^(proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_|ev_)")
}

# Prédicteurs avec p < 0.05 "à un moment donné" (modèle ajusté) 
# = dans au moins une analyse, avec ou sans outliers (hors EV normalisés) + fenêtres y_to_als
y_to_als_analyses <- c("main_y_5_10ans", "main_y_10_15ans", "main_y_15_20ans")

# (les deux objets n'ont pas les mêmes types de colonnes, ex. lower_CI numérique vs texte :
#  on extrait donc les noms séparément au lieu de les empiler avec bind_rows)
sig_from <- function(results) {
  results |>
    filter(model == "adjusted", term == "Continuous", as.numeric(p_value_raw) < 0.05) |>
    pull(explanatory)
}

sig_names <- unique(c(
  sig_from(filter(main_results, !str_starts(analysis, "sensi_ev_"))),
  sig_from(filter(results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results,
                  analysis %in% y_to_als_analyses))))

prot_signi <- proteomic[clean_name(proteomic) %in% sig_names]


# Graphique pour une variable, avec ou sans retrait des sets contenant un outlier
# Mêmes sets que dans les analyses (prep_data) : sets avec valeur manquante retirés ;
# si outlier = TRUE, sets contenant un outlier (|z| > outlier_k) retirés pour ce prédicteur.
ratio_plot <- function(col, outlier = FALSE) {
  name <- clean_name(col)
  prep <- prep_data(bdd_danish, col, outlier)
  
  prep$data |>
    group_by(match) |>
    summarise(diff_val  = .data[[col]][als == 1] - mean(.data[[col]][als == 0]),   # différence sur la variable non standardisée
              follow_up = unique(follow_up[als == 1]),
              .groups = "drop") |>
    ggplot(aes(x = -follow_up / 12, y = diff_val)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.75) +
    labs(x = "Time to ALS diagnosis (years)",
         y = paste0("Difference in ", name, " between case and matched controls"),
         title = name,
         subtitle = if (outlier) sprintf("Outliers excluded (%d matched sets removed)", length(prep$out_sets))
         else "All matched sets with complete data") +
    theme_lucid() +
    theme(text = element_text(size = 16, color = "black", face = "bold"),
          axis.text = element_text(size = 16, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5))
}

# -> plots_list$CD9$all_sets, plots_list$CD9$no_outliers, ...
plots_list <- set_names(prot_signi, clean_name(prot_signi)) |>
  map(\(col) list(all_sets    = ratio_plot(col, outlier = FALSE),
                  no_outliers = ratio_plot(col, outlier = TRUE)))

## Figures : ALS risk score simple ----
# Score = somme pondérée des z-scores des prédicteurs significatifs de l'analyse principale (poids = log OR par SD).
# Sorties : table OR (clogit), plot LOESS, courbes ROC par strate de temps avant diagnostic.
# Nécessite (script d'analyses) : main_results, proteomic, bdd_danish, models, keep_cases, q50, cut_years
# Nécessite aussi le package pROC.
# ATTENTION : prédicteurs et poids viennent des mêmes données que l'évaluation -> performance apparente (optimiste).

# Paramètres 
score_model <- "adjusted"   # modèle dont viennent les poids : "base" ou "adjusted"
p_thr       <- 0.05         # seuil de sélection (p_value_raw, analyse principale)
winsor      <- 3            # |z| plafonné à 3 (NULL = aucun)

clean_name <- function(x) {                                                     # même nettoyage que dans main_results
  x |>
    str_remove("_log10$") |>
    str_remove("^(proteomic_immun_res_|proteomic_metabolism_|proteomic_neuro_explo_|ev_)")
}
group_of <- \(x) case_when(str_detect(x, "proteomic_immun_res") ~ "Immune response",
                           str_detect(x, "proteomic_metabolism") ~ "Metabolism",
                           str_detect(x, "proteomic_neuro_explo") ~ "Neuro-exploratory",
                           str_detect(x, "ev_") ~ "Extracellular Vesicles")
std_z <- function(x, mu = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE), w = winsor) {   # z-score winsorisé
  x[!is.finite(x)] <- NA
  z <- (x - mu) / s
  if (!is.null(w)) z <- pmax(pmin(z, w), -w)
  z
}
na_by_set <- function(x, match) ave(x, match, FUN = \(v) if (anyNA(v)) rep(NA_real_, length(v)) else v)


# 1. Construction du score 
var_map <- tibble(var = proteomic, explanatory = clean_name(proteomic), protein_group = group_of(proteomic))
stopifnot(!anyDuplicated(var_map[c("explanatory", "protein_group")]))

weights <- main_results |>
  filter(analysis == "main", model == score_model, term == "Continuous", p_value_raw < p_thr) |>
  transmute(explanatory, protein_group, beta = log(OR_raw)) |>
  left_join(var_map, by = c("explanatory", "protein_group"))
stopifnot(!anyNA(weights$var))

Z <- sapply(weights$var, \(v) std_z(bdd_danish[[v]]))                           # 1 colonne par prédicteur retenu
bdd_score <- bdd_danish |>
  mutate(score = as.numeric(Z %*% weights$beta),                                # somme pondérée
         score = as.numeric(scale(na_by_set(score, match))))                    # 1 SD ; NA pour tout le set si un membre a une valeur manquante

# Strates de temps avant diagnostic (mêmes que la table S1 / la figure 1)
d3 <- keep_cases(bdd_score, follow_up > 60)
strata_labels <- c("All cases and controls", "Years to ALS diagnosis < 5 years",
                   sprintf("Years to ALS diagnosis between 5 and %s years", cut_years),
                   sprintf("Years to ALS diagnosis > %s years", cut_years))
strata <- list(bdd_score,
               keep_cases(bdd_score, follow_up < 60),
               keep_cases(d3, follow_up <= q50(d3)),
               keep_cases(d3, follow_up >  q50(d3))) |>
  set_names(strata_labels) |>
  map(\(x) filter(x, !is.na(score)))                                            # sets complets uniquement


### 1. Table : OR de la relation ALS ~ score (régression logistique conditionnelle) ----
fit_or <- function(data, covars) {
  clogit(reformulate(c("score", covars, "strata(match)"), response = "als"), data = data) |>
    tidy(exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "score")
}

score_or <- imap_dfr(strata, \(data, stratum) {
  map_dfr(c("base", "adjusted"), \(m) mutate(fit_or(data, models[[m]]), model = m)) |>
    mutate(stratum = stratum, n = nrow(data), n_sets = n_distinct(data$match))
})

score_tab <- score_or |>
  mutate(OR_ci = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p     = if_else(p.value < 0.01, "<0.01", sprintf("%.2f", p.value)),
         N     = sprintf("%d (%d)", n, n_sets)) |>
  select(stratum, N, model, OR_ci, p) |>
  pivot_wider(names_from = model, values_from = c(OR_ci, p))

keys <- c("stratum", "N", "OR_ci_base", "p_base", "OR_ci_adjusted", "p_adjusted")
header <- data.frame(
  col_keys = keys,
  top      = c("Time to ALS diagnosis", "N participants (sets)", "Base model", "Base model", "Adjusted model", "Adjusted model"),
  bottom   = c("Time to ALS diagnosis", "N participants (sets)", "OR (95% CI)", "p-value", "OR (95% CI)", "p-value"))

table_ALS_score_simple <- flextable(score_tab, col_keys = keys) |>
  set_header_df(mapping = header, key = "col_keys") |>
  add_footer_lines(c(
    "1Conditional logistic regressions matched on birth year and sex. Adjusted models further account for smoking and body mass index.",
    sprintf("2Odds ratio for a one standard deviation increase in the ALS risk score (weighted sum of the standardized levels of the %d proteins/markers with p<0.05 in the main analysis; weights = log OR).", nrow(weights)),
    "3Predictors and weights were derived from the same data: apparent performance, likely optimistic.")) |>
  theme_vanilla() |>
  merge_h(part = "header") |>
  merge_v(part = "header") |>
  align(align = "center", part = "all") |>
  align(j = "stratum", align = "left", part = "all") |>
  bold(j = "stratum", part = "body") |>
  flextable::font(fontname = "Calibri", part = "all") |>
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


### 2. Plot LOESS : différence de score cas - témoins appariés vs temps avant diagnostic ----
plot_theme <- theme(axis.title = element_text(size = 16, color = "black", face = "bold"),
                    axis.text  = element_text(size = 16, color = "black", face = "bold"))

plot_loess <- bdd_score |>
  filter(!is.na(score)) |>
  group_by(match) |>
  summarise(diff_val  = score[als == 1] - mean(score[als == 0]),                # score standardisé : on trace une différence, pas un ratio
            follow_up = unique(follow_up[als == 1]),
            .groups = "drop") |>
  ggplot(aes(x = -follow_up / 12, y = diff_val)) +
  geom_point(alpha = 0.6, color = "grey30") +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.75, color = "black", fill = "grey70") +
  labs(x = "Time to ALS diagnosis (years)",
       y = "Difference in ALS risk score\nbetween case and matched controls") +
  theme_lucid() +
  plot_theme


### 3. Courbes ROC (une par strate), avec seuil optimal (Youden) ----
roc_list <- map2(strata, strata_labels, \(data, lab) {
  r    <- pROC::roc(data$als, data$score, levels = c(0, 1), direction = "<", quiet = TRUE)   # cas = score plus élevé
  best <- pROC::coords(r, "best", best.method = "youden",
                       ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)[1, ]
  list(curve = pROC::coords(r, "all", ret = c("specificity", "sensitivity"), transpose = FALSE),
       label = sprintf("%s (n=%d)\nAUC = %.2f\nOptimal score cut-off: %.2f\n(sensitivity: %.2f and specificity: %.2f)",
                       lab, nrow(data), as.numeric(pROC::auc(r)), best$threshold, best$sensitivity, best$specificity))
})

roc_df <- imap_dfr(roc_list, \(x, nm) mutate(x$curve, fpr = 1 - specificity, label = x$label)) |>
  mutate(label = factor(label, levels = map_chr(roc_list, "label")))

plot_roc <- ggplot(roc_df, aes(x = fpr, y = sensitivity, linetype = label)) +
  geom_path(color = "black") +
  geom_abline(linetype = "dashed", color = "grey50") +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed", "dotted")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity") +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_lucid() +
  plot_theme +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 12, hjust = 0), legend.key.width = unit(1.5, "cm"))

### 4. Heatmap de corrélation entre les prédicteurs du score ----
# Nécessite (script als_risk_score_simple.R) : Z (z-scores winsorisés des prédicteurs) et weights

cor_method <- "pearson"                                                         # "spearman" pour une version robuste aux outliers

cm <- cor(Z, use = "pairwise.complete.obs", method = cor_method)
dimnames(cm) <- list(weights$explanatory, weights$explanatory)                  # noms courts (IL32, NEFL, ...)

ord <- hclust(as.dist(1 - cm), method = "average")$order                        # regroupe les protéines les plus corrélées
cm  <- cm[ord, ord]

heat_df <- as.data.frame(as.table(cm)) |>
  setNames(c("x", "y", "r")) |>
  mutate(x = factor(x, levels = rownames(cm)),
         y = factor(y, levels = rev(rownames(cm))),                             # diagonale en haut à gauche
         txt_color = if_else(abs(r) > 0.6, "white", "black"))

plot_heatmap <- ggplot(heat_df, aes(x = x, y = y, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r), color = txt_color), size = 2.3) +
  scale_color_identity() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        panel.grid = element_blank())




# Publication ----
## Table 1 ----
table_1 <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"),
    follow_up = follow_up/12) |>
  select(
    als, birth_year, baseline_age, diagnosis_age, follow_up,  
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            follow_up ~ 0, 
                            bmi ~ 1)) |>
  bold_labels() |>
  add_p(include = -c('diagnosis_age', 'follow_up')) |>
  add_overall() |>
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**N missing**") |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")


get_res <- function(analyses, sfx = "") {
  main_results |>
    filter(term == "Continuous", model == "adjusted", analysis %in% paste0(analyses, sfx)) |>
    mutate(analysis = str_remove(analysis, "_outlier$"))
}

fig1_analyses <- c("main", "sensi_2", "sensi_3_4", "sensi_3_5")


## Figure 1 - Main analysis (volcano plot global et stratifié sur le temps avant diagnostic) ----
strip_lab <- function(title, n) {
  sprintf("<b>%s</b><br><span style='font-size:9pt; color:grey40;'>Conditional logistic regressions (N=%d)<br>Matched on birth year and sex and adjusted on BMI and smoking status</span>",
          title, n)
}
fig1_labels <- c(
  main      = strip_lab("All cases and controls", n_analysis[["main"]]),
  sensi_2   = strip_lab("Years to ALS diagnostic < 5 years", n_analysis[["sensi_2"]]),
  sensi_3_4 = strip_lab(sprintf("Years to ALS diagnostic between 5 and %s years", cut_years), n_analysis[["sensi_3_4"]]),
  sensi_3_5 = strip_lab(sprintf("Years to ALS diagnostic > %s years", cut_years), n_analysis[["sensi_3_5"]]))

volcano_facets <- function(data, labels) {                                      # volcano plot à facettes (labels = vecteur nommé analyse -> titre HTML)
  data |>
    mutate(Significance = case_when(p_value_raw < 0.05 & OR_raw > 1 ~ "OR>1 & p-value<0.05",
                                    p_value_raw < 0.05 & OR_raw < 1 ~ "OR<1 & p-value<0.05",
                                    TRUE ~ "p-value>0.05"),
           analysis = factor(analysis, levels = names(labels), labels = labels)) |>
    ggplot(aes(x = OR_raw, y = -log10(p_value_raw), color = Significance)) +
    geom_point(alpha = 0.8, size = 2, show.legend = FALSE) +
    geom_text_repel(data = ~filter(.x, p_value_raw < 0.05), aes(label = explanatory),
                    size = 3.5, max.overlaps = 20, box.padding = 0.4, point.padding = 0.2,
                    segment.color = "grey20", color = "black") +
    scale_color_manual(values = c("OR<1 & p-value<0.05" = "blue",
                                  "p-value>0.05" = "grey70",
                                  "OR>1 & p-value<0.05" = "red")) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(x = "OR", y = "-log10(p-value)", color = "") +
    facet_wrap(vars(analysis), ncol = 2L, nrow = 3L) +
    theme(strip.text = element_markdown(hjust = 0, size = 12, lineheight = 1.1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1.5, "lines"))
}
make_figure_1 <- \(sfx = "") volcano_facets(get_res(fig1_analyses, sfx), fig1_labels)
figure_1 <- map(suffixes, make_figure_1)  


## Figure 2 - Forest plot stratifié sur le sexe (défini dans le script de figures) ----
figure_2 <- list(standard = forest_sex_main, outlier = forest_sex_outlier)

## Figure 3 - Machine learning tidymodels ----

## Figure 4 option 1 - ALS risk score simple à partir des resultats des regression logistique ----
figure_4_opt1 = plot_loess / plot_roc

## Figure 4 option 2 - ALS risk score complexe à partir des resultats machine learning tidymodels ----


## Table S1 - Main results (global + stratifié sur le temps avant diagnostic) ----
footer_paragraph <- function(footer) {                                          # notes numérotées (numéro en exposant), une par ligne
  do.call(as_paragraph,
          do.call(c, imap(footer, \(x, i) list(as_sup(as.character(i)),
                                               paste0(x, if (i < length(footer)) "\n" else "")))))
}
when_rows <- function(ft, i, f) if (length(i) > 0) f(ft, i) else ft            # n'applique f que si des lignes sont concernées

make_table_S <- function(data, groups, top_notes, or_note, ci_note, footer) {
  stats <- c("OR", "95% CI", "p_value")
  keys  <- c("protein_group", "explanatory", as.vector(outer(stats, names(groups), paste, sep = "_")))
  
  wide <- data |>
    group_by(explanatory) |>
    filter(any(p_value_raw < 0.05, na.rm = TRUE)) |>                            # prédicteurs significatifs dans au moins une analyse
    ungroup()
  
  if (nrow(wide) == 0) return(NULL)                                             # aucun prédicteur significatif : pas de tableau
  
  wide <- wide |>
    select(analysis, protein_group, explanatory, all_of(stats), p_value_raw, OR_raw) |>
    pivot_wider(names_from = analysis, values_from = c(all_of(stats), p_value_raw, OR_raw)) |>
    arrange(protein_group, explanatory)
  
  header <- data.frame(
    col_keys = keys,
    top      = c("Protein group", "Pre-disease plasma proteins", rep(unname(groups), each = 3)),
    bottom   = c("Protein group", "Pre-disease plasma proteins", rep(c("OR", "95% CI", "p-value"), length(groups))))
  
  ft <- flextable(wide, col_keys = keys) |>
    set_header_df(mapping = header, key = "col_keys")
  
  # Mise en forme conditionnelle : p < 0.05 en gras ; OR>1 orange, OR<1 vert
  for (g in names(groups)) {
    cols <- paste(stats, g, sep = "_")
    p <- wide[[paste0("p_value_raw_", g)]]
    o <- wide[[paste0("OR_raw_", g)]]
    ft <- ft |>
      when_rows(which(p < 0.05),           \(x, i) bold(x, i = i, j = paste0("p_value_", g))) |>
      when_rows(which(p < 0.05 & o > 1),   \(x, i) bg(x, i = i, j = cols, bg = "#FFE5D9")) |>
      when_rows(which(p < 0.05 & o < 1),   \(x, i) bg(x, i = i, j = cols, bg = "#D8F3DC"))
  }
  
  # Notes en exposant dans l'en-tête
  for (note in unique(top_notes)) {
    gs <- names(top_notes)[top_notes == note]
    ft <- append_chunks(ft, i = 1, j = as.vector(outer(stats, gs, paste, sep = "_")),
                        part = "header", as_sup(note))
  }
  ft <- ft |>
    append_chunks(i = 2, j = paste0("OR_",     names(groups)), part = "header", as_sup(or_note)) |>
    append_chunks(i = 2, j = paste0("95% CI_", names(groups)), part = "header", as_sup(ci_note))
  
  ft |>
    add_footer_row(top = FALSE, values = footer_paragraph(footer), colwidths = length(keys)) |>
    theme_vanilla() |>
    merge_h(part = "header") |>
    merge_v(part = "header") |>
    merge_v(j = "protein_group") |>
    align(align = "center", part = "all") |>
    align(j = c("protein_group", "explanatory"), align = "left", part = "all") |>
    bold(j = c("protein_group", "explanatory"), part = "body") |>
    flextable::font(fontname = "Calibri", part = "all") |>
    fontsize(size = 10, part = "all") |>
    padding(padding.top = 0, padding.bottom = 0, part = "all")
}

S1_groups <- c(
  main      = sprintf("Main analyses\nAll cases and controls (N=%d)", n_analysis[["main"]]),
  sensi_2   = sprintf("Filtered to years to ALS diagnosis < 5 years (N=%d)", n_analysis[["sensi_2"]]),
  sensi_3_4 = sprintf("Filtered to years to ALS diagnosis between 5 and %s years (N=%d)", cut_years, n_analysis[["sensi_3_4"]]),
  sensi_3_5 = sprintf("Filtered to years to ALS diagnosis > %s years (N=%d)", cut_years, n_analysis[["sensi_3_5"]]))

table_S1 <- map(suffixes, \(sfx)
                make_table_S(get_res(fig1_analyses, sfx),
                             groups    = S1_groups,
                             top_notes = c(main = "1", sensi_2 = "1", sensi_3_4 = "1", sensi_3_5 = "1"),
                             or_note   = "2", ci_note = "3",
                             footer    = c("All models are matched on birth year and sex, and adjusted for smoking and body mass index.",
                                           "Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.",
                                           "CI: Confidence interval.")))

## Table S2 - Sensibilité : fenêtres fixes de temps avant diagnostic (la stratification par sexe est dans S3) ----
# Les fenêtres 5-10 / 10-15 / 15-20 ans viennent de l'objet results_proteomic_ALS_occurrence_y_to_als (autre pipeline).
# Ses colonnes n'ont pas les mêmes types que main_results (ex. lower_CI en texte) : on n'en garde que les
# colonnes utiles à la table, converties au bon type, avant d'empiler les deux sources.
S2_groups <- c(
  main            = sprintf("Main analyses\nAll cases and controls (N=%d)", n_analysis[["main"]]),
  main_y_5_10ans  = "Sensitivity analyses\nCases diagnosed within 5 to 10 years (27 cases - 468 controls)",    # N repris de l'ancien code
  main_y_10_15ans = "Sensitivity analyses\nCases diagnosed within 10 to 15 years (54 cases - 441 controls)",
  main_y_15_20ans = "Sensitivity analyses\nCases diagnosed within 15 to 20 years (57 cases - 438 controls)")

S2_cols <- c("analysis", "protein_group", "explanatory", "OR", "95% CI", "p_value", "p_value_raw", "OR_raw")

S2_data <- bind_rows(
  get_res("main") |>
    select(all_of(S2_cols)),
  results_proteomic_ALS_occurrence_y_to_als$logistic_fixed_horizon$proteome_wide$main_results |>
    filter(term == "Continuous", model == "adjusted", analysis %in% names(S2_groups)) |>
    mutate(OR_raw      = as.numeric(OR),
           p_value_raw = as.numeric(p_value_raw),
           across(c(OR, `95% CI`, p_value), as.character),
           explanatory = str_remove(explanatory, "_log10$")) |>                 # même nom que les EV de main_results
    select(all_of(S2_cols)))

table_S2 <- make_table_S(
  S2_data, groups = S2_groups,
  top_notes = c(main = "1", main_y_5_10ans = "2", main_y_10_15ans = "2", main_y_15_20ans = "2"),
  or_note = "3", ci_note = "4",
  footer = c("Models are matched on birth year and sex, and adjusted for smoking and body mass index.",
             "Models are adjusted on birth year, sex, smoking and body mass index.",
             "Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.",
             "CI: Confidence interval."))



## Table S3 - Sensibilité stratifiée sur le sexe ----
S3_groups <- c(
  main           = sprintf("Main analyses\nAll cases and controls (N=%d)", n_analysis[["main"]]),
  sensi_7_male   = sprintf("Sensitivity analyses\nFiltered to males (N=%d)", n_analysis[["sensi_7_male"]]),
  sensi_7_female = sprintf("Sensitivity analyses\nFiltered to females (N=%d)", n_analysis[["sensi_7_female"]]))

table_S3 <- map(suffixes, \(sfx)
                make_table_S(get_res(names(S3_groups), sfx),
                             groups    = S3_groups,
                             top_notes = c(main = "1", sensi_7_male = "2", sensi_7_female = "2"),
                             or_note   = "3", ci_note = "4",
                             footer    = c("Models are matched on birth year and sex, and adjusted for smoking and body mass index.",
                                           "Models are adjusted on birth year, smoking and body mass index (sex-stratified models).",
                                           "Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.",
                                           "CI: Confidence interval.")))


## Table S4 - Analyses de la table S1 avant / après retrait des outliers ----
# pairs : tibble(ref, alt, label) ; ref / alt = codes d'analyses de main_results, label = groupe de lignes
# Un (prédicteur, paire) est retenu si :
#  - le sens de l'OR change (OR de part et d'autre de 1), seulement si p < 0.05 dans au moins l'une des 2 analyses
#    (sinon : simple bruit autour de 1) ;
#  - et/ou la significativité change (p < 0.05 -> p >= 0.05, ou l'inverse).
# Tous les prédicteurs sont examinés (pas seulement ceux de la table S1).
res_adj <- main_results |> filter(term == "Continuous", model == "adjusted")

make_table_diff <- function(pairs, group_title, ref_title, alt_title, footer, dir_requires_sig = TRUE) {
  stats <- c("OR", "95% CI", "p_value")
  keys  <- c("group", "protein_group", "explanatory", paste0(stats, "_ref"), paste0(stats, "_alt"), "n_extra", "change")
  
  diff <- pmap_dfr(pairs, \(ref, alt, label) {
    res_adj |>
      filter(analysis %in% c(ref, alt)) |>
      mutate(side = if_else(analysis == ref, "ref", "alt")) |>
      select(side, protein_group, explanatory, all_of(stats), p_value_raw, OR_raw, n_sets_na, n_sets_outlier) |>
      pivot_wider(names_from = side,
                  values_from = c(all_of(stats), p_value_raw, OR_raw, n_sets_na, n_sets_outlier)) |>
      mutate(group = label)
  }) |>
    filter(!is.na(OR_raw_ref), !is.na(OR_raw_alt)) |>                           # prédicteurs présents dans les 2 analyses
    mutate(dir  = (sign(log(OR_raw_ref)) != sign(log(OR_raw_alt)) &
                     (!dir_requires_sig | p_value_raw_ref < 0.05 | p_value_raw_alt < 0.05)) %in% TRUE,
           loss = (p_value_raw_ref <  0.05 & p_value_raw_alt >= 0.05) %in% TRUE,
           gain = (p_value_raw_ref >= 0.05 & p_value_raw_alt <  0.05) %in% TRUE,
           change = paste0(if_else(dir, "OR direction reversed", ""),
                           if_else(dir & (loss | gain), "; ", ""),
                           case_when(loss ~ "significant \u2192 not significant",
                                     gain ~ "not significant \u2192 significant",
                                     TRUE ~ "")),
           n_extra = (n_sets_na_alt + n_sets_outlier_alt) - (n_sets_na_ref + n_sets_outlier_ref),  # sets exclus en plus dans `alt`
           group = factor(group, levels = pairs$label)) |>
    filter(dir | loss | gain) |>
    arrange(group, protein_group, explanatory) |>
    mutate(group = as.character(group))
  
  if (nrow(diff) == 0) return(NULL)                                             # aucune différence majeure
  
  header <- data.frame(
    col_keys = keys,
    top      = c(group_title, "Protein group", "Pre-disease plasma proteins",
                 rep(ref_title, 3), rep(alt_title, 3), "Additional sets excluded", "Change"),
    bottom   = c(group_title, "Protein group", "Pre-disease plasma proteins",
                 rep(c("OR", "95% CI", "p-value"), 2), "Additional sets excluded", "Change"))
  
  ft <- flextable(diff, col_keys = keys) |>
    set_header_df(mapping = header, key = "col_keys")
  
  for (side in c("ref", "alt")) {                                               # p < 0.05 en gras ; OR>1 orange, OR<1 vert
    cols <- paste(stats, side, sep = "_")
    p <- diff[[paste0("p_value_raw_", side)]]
    o <- diff[[paste0("OR_raw_", side)]]
    ft <- ft |>
      when_rows(which(p < 0.05),         \(x, i) bold(x, i = i, j = paste0("p_value_", side))) |>
      when_rows(which(p < 0.05 & o > 1), \(x, i) bg(x, i = i, j = cols, bg = "#FFE5D9")) |>
      when_rows(which(p < 0.05 & o < 1), \(x, i) bg(x, i = i, j = cols, bg = "#D8F3DC"))
  }
  
  ft |>
    append_chunks(i = 1, j = c(paste0(stats, "_ref"), paste0(stats, "_alt")), part = "header", as_sup("1")) |>
    append_chunks(i = 2, j = c("OR_ref", "OR_alt"),         part = "header", as_sup("2")) |>
    append_chunks(i = 2, j = c("95% CI_ref", "95% CI_alt"), part = "header", as_sup("3")) |>
    append_chunks(i = 1, j = "change",                      part = "header", as_sup("4")) |>
    add_footer_row(top = FALSE, values = footer_paragraph(footer), colwidths = length(keys)) |>
    theme_vanilla() |>
    merge_h(part = "header") |>
    merge_v(part = "header") |>
    merge_v(j = "group") |>
    align(align = "center", part = "all") |>
    align(j = c("group", "protein_group", "explanatory", "change"), align = "left", part = "all") |>
    bold(j = c("group", "explanatory"), part = "body") |>
    flextable::font(fontname = "Calibri", part = "all") |>
    fontsize(size = 10, part = "all") |>
    padding(padding.top = 0, padding.bottom = 0, part = "all")
}

footer_common <- c("Estimated risk of ALS associated with a one standard deviation increase in pre-disease plasma concentration of proteins.",
                   "CI: Confidence interval.")
footer_change <- "Only major changes are shown: OR direction reversed (OR on the other side of 1, with p<0.05 in at least one of the two analyses) and/or change in statistical significance (p<0.05). Additional sets excluded: matched sets (cases and their controls) removed for the marker in the right-hand analysis compared with the left-hand one."


S4_labels <- c(main      = "All cases and controls",
               sensi_2   = "Years to diagnosis < 5",
               sensi_3_4 = sprintf("Years to diagnosis between 5 and %s", cut_years),
               sensi_3_5 = sprintf("Years to diagnosis > %s", cut_years))
S4_pairs <- tibble(ref = names(S4_labels), alt = paste0(names(S4_labels), "_outlier"), label = unname(S4_labels))

table_S4 <- make_table_diff(
  S4_pairs, group_title = "Analysis", ref_title = "All data", alt_title = "Outliers excluded",
  footer = c(sprintf("All models are matched on birth year and sex, and adjusted for smoking and body mass index. Outliers excluded: for each marker, matched sets (cases and their controls) with at least one participant beyond %s SD from the mean were removed.", outlier_k),
             footer_common, footer_change))


## Table S5 - Analyses principales EV vs EV normalisés (tétraspanines) ----
ev_label <- map_chr(ev_norms, \(d) if (length(d) == 1) d else paste0("(", paste(d, collapse = " + "), ")"))
S5_pairs <- function(sfx) tibble(ref = paste0("main", sfx), alt = paste0(ev_analyses, sfx),
                                 label = paste("Marker X /", unname(ev_label)))
S5_footer <- function(sfx) c(
  paste("All models are matched on birth year and sex, and adjusted for smoking and body mass index.",
        "Main analyses use log10-transformed EV intensities; in the normalized analyses, raw EV intensities were divided by the indicated marker(s) before log10 transformation.",
        if (sfx != "") sprintf("In both analyses, for each marker, matched sets with at least one participant beyond %s SD from the mean were removed.", outlier_k)),
  footer_common, footer_change)

table_S5 <- map(suffixes, \(sfx)
                make_table_diff(S5_pairs(sfx), group_title = "Normalization", ref_title = "Main analysis",
                                alt_title = "Normalized", footer = S5_footer(sfx)))
all_sets_plots <- map(plots_list, "all_sets")   # garde seulement all_sets


## Figure S1 Analyse principale + fenêtres fixes de temps avant diagnostic (mêmes données que la table S2) ----
strip_lab_win <- function(title, n_cases, n_controls) {                         # N des fenêtres : repris de l'ancien code
  sprintf("<b>%s</b><br><span style='font-size:9pt; color:grey40;'>Logistic regressions (%d cases - %d controls)<br>Adjusted on birth year, sex, BMI and smoking status</span>",
          title, n_cases, n_controls)
}
figS1_labels <- c(
  main            = strip_lab("All cases and controls", n_analysis[["main"]]),
  main_y_5_10ans  = strip_lab_win("Cases diagnosed within 5 to 10 years", 27, 468),
  main_y_10_15ans = strip_lab_win("Cases diagnosed within 10 to 15 years", 54, 441),
  main_y_15_20ans = strip_lab_win("Cases diagnosed within 15 to 20 years", 57, 438))

figure_S1 <- volcano_facets(S2_data, figS1_labels)

## Figure S2 - Machine learning tidymodels ----

## Figure S3 - biological network? ----

# Regroupement de tous les résultats ----
results_proteomic_ALS_occurrence <- list(
  covar = covar,                                                                # effets des covariables
  main = list(
    main_results = main_results,                                                # tous les modèles : analyses x modèles x (avec/sans outliers)
    tables = list(
      by_analysis        = tables_std,                                          # $main, $sensi_2, $sensi_2_outlier, $sensi_ev_CD9, ...
      outlier_comparison = tables_outlier,                                      # $<analyse>$base / $adjusted : avec vs sans outliers
      sex                = tables_sex,                                          # $standard / $outlier
      table_ALS_score_simple = table_ALS_score_simple),                         # ALS risk score from logistic regression results                  
    figures = list(
      volcano      = figs,                                                       # $<analyse>$base / $adjusted
      volcano_sex  = figs_sex,                                                   # $standard / $outlier -> $base / $adjusted
      forest_sex   = figure_2,                                                   # $standard / $outlier
      ratio        = plots_list,                                                 # $<prédicteur>$all_sets / $no_outliers
      plot_roc     = plot_roc, 
      plot_loess   = plot_loess, 
      plot_heatmap = plot_heatmap)),                                               
  publication = list(
      table_1       = table_1, 
      figure_1      = figure_1,                                                 # $standard / $outlier
      figure_2      = figure_2,
      figure_4_opt1 = figure_4_opt1, 
      table_S1      = table_S1,
      table_S2      = table_S2, 
      table_S3      = table_S3,
      table_S4      = table_S4,                                                 # S1 avant / après retrait des outliers
      table_S5      = table_S5,                                                 # EV principaux vs EV normalisés ($standard / $outlier)
      figure_S1     = figure_S1))                                               

rm(tables_std, tables_outlier, tables_sex, figs, figs_sex, forest_sex_main, forest_sex_outlier,
   plots_list, plot_roc, plot_loess, plot_heatmap,
   table_1, figure_1, figure_2, figure_4_opt1, 
   table_S1, table_S2, table_S3, table_S4, table_S5, figure_S1)                       



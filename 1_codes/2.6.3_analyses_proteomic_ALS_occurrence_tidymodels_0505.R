# Aline Davias
# Juin 2026
# Prédiction ALS par protéomique – tidymodels (v2)
# 
# Structure :
#   Test 1 : protéines + covariables (birth_year, sex, bmi, smoking_2cat_i)
#   Test 2 : idem + follow_up_no_na_y
#   Chaque test :
#     Phase A : tuning hyperparamètres modèle par modèle (tune_grid / tune_bayes)
#     Phase B : comparaison finale workflow_set, hyperparamètres fixés
#     Phase C : analyse détaillée du meilleur modèle (SHAP, calibration, etc.)
#
# Points clés :
#   - Early stopping xgboost via xgb.cv (même approche que ton script xgboost)
#   - Matching 1:2 géré par group_vfold_cv sur `match`
#   - Filtrage univarié (top-N) via step custom basé sur pROC (CRAN stable)
#   - tune_bayes pour glmnet et RF (plus efficace que grid_regular sur p >> n)
#   - grid_latin_hypercube pour les espaces continus (xgb, svm, mars)
#   - Vérification anti-leakage des folds



# Packages and data ----
library(tidymodels)
library(tidyverse)
library(ranger)       # random forest
library(xgboost)      # gradient boosting
library(glmnet)       # lasso/ridge/elastic net
library(kernlab)      # SVM
library(earth)        # MARS
library(mgcv)         # GAM
library(pROC)         # AUC univariate (pour le step de filtrage custom)
library(shapviz)      # SHAP values
library(probably)     # calibration des probabilités
library(vip)          # variable importance

tidymodels_prefer()
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")


# Fonctions ----
step_select_auc_top <- function(
    recipe, ..., outcome, top_n = 20, role = NA,
    trained = FALSE, selected_vars = NULL, skip = FALSE,
    id = recipes::rand_id("select_auc_top")) {
  
  recipes::add_step(
    recipe,
    structure(
      list(
        terms       = rlang::enquos(...),
        outcome     = outcome,
        top_n       = top_n,
        role        = role,
        trained     = trained,
        selected_vars = selected_vars,
        skip        = skip,
        id          = id),
      class = c("step_select_auc_top", "step")))
}

prep.step_select_auc_top <- function(x, training, info = NULL, ...) {
  
  # Variables candidates (numériques)
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  
  # Outcome binaire (0/1 ou factor)
  y_raw <- training[[x$outcome]]
  y_num <- as.numeric(y_raw) - ifelse(is.factor(y_raw), 1, 0)
  # → si factor avec levels 1/2 → devient 0/1
  # → si numeric 0/1 → reste 0/1
  
  # AUC univariée pour chaque prédicteur candidat
  auc_scores <- vapply(col_names, function(v) {
    tryCatch(
      as.numeric(pROC::auc(y_num, training[[v]], quiet = TRUE)),
      error = function(e) 0.5)
  }, numeric(1))
  
  # Symétriser l'AUC (AUC < 0.5 → 1 - AUC, direction inversée mais info réelle)
  auc_scores <- pmax(auc_scores, 1 - auc_scores)
  
  # Garder les top_n
  n_keep        <- min(x$top_n, length(auc_scores))
  top_idx       <- order(auc_scores, decreasing = TRUE)[seq_len(n_keep)]
  selected_vars <- col_names[top_idx]
  
  x$selected_vars <- selected_vars
  x$trained       <- TRUE
  x
}

bake.step_select_auc_top <- function(object, new_data, ...) {
  # Garder l'outcome + les variables non-candidates + les variables sélectionnées
  all_vars <- names(new_data)
  # Variables candidates qui ne sont PAS sélectionnées → à supprimer
  # (on ne peut pas récupérer les candidats directement en bake, 
  #  donc on supprime tout ce qui n'est ni sélectionné ni non-candidat)
  # Approche simple : on ne supprime rien ici — le filtrage est fait en prep
  # En pratique, on utilise recipes::check_new_data + select
  new_data |> dplyr::select(dplyr::all_of(
    union(object$selected_vars,
          setdiff(all_vars, 
                  # on retire les vars numériques qui ne sont pas sélectionnées
                  # SAUF l'outcome et les variables de rôle spécial
                  setdiff(all_vars, object$selected_vars)))))
}

# Version simplifiée du bake (plus robuste) :
bake.step_select_auc_top <- function(object, new_data, ...) {
  # Colonnes à conserver : tout sauf les prédicteurs non sélectionnés
  # On identifie les prédicteurs non sélectionnés en regardant les numeric
  # qui ne sont pas dans selected_vars et ne sont pas l'outcome
  non_selected_preds <- setdiff(
    names(new_data)[sapply(new_data, is.numeric)],
    object$selected_vars)
  
  # On ne retire que les colonnes candidates non sélectionnées
  # MAIS on doit savoir quelles colonnes étaient candidates → problème de bake
  # Solution : stocker aussi les rejected_vars en prep
  cols_to_keep <- setdiff(names(new_data), object$rejected_vars)
  new_data |> dplyr::select(dplyr::all_of(cols_to_keep))
}

# Version finale propre avec rejected_vars :
prep.step_select_auc_top <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  y_raw <- training[[x$outcome]]
  y_num <- if (is.factor(y_raw)) as.numeric(y_raw) - 1 else as.numeric(y_raw)
  
  auc_scores <- vapply(col_names, function(v) {
    tryCatch(
      as.numeric(pROC::auc(y_num, training[[v]], quiet = TRUE)),
      error = function(e) 0.5)
  }, numeric(1))
  auc_scores <- pmax(auc_scores, 1 - auc_scores)
  
  n_keep        <- min(x$top_n, length(auc_scores))
  top_idx       <- order(auc_scores, decreasing = TRUE)[seq_len(n_keep)]
  selected_vars <- col_names[top_idx]
  
  x$selected_vars <- selected_vars
  x$rejected_vars  <- setdiff(col_names, selected_vars)
  x$trained        <- TRUE
  x
}

# Méthode print pour lisibilité
print.step_select_auc_top <- function(x, width = max(20, options()$width - 35), ...) {
  cat("AUC-based top-", x$top_n, " feature selection\n", sep = "")
  invisible(x)
}





# Folds groupés sur `match` + vérification anti-leakage
make_folds_checked <- function(data, v = 10, seed = 1996) {
  data_with_match <- data  # match est déjà dans data
  
  set.seed(seed)
  folds <- group_vfold_cv(data_with_match, group = match, v = v)
  
  # Vérification leakage
  has_leakage <- FALSE
  for (i in seq_len(nrow(folds))) {
    tr_ids   <- training(folds$splits[[i]])$match
    te_ids   <- testing(folds$splits[[i]])$match
    leakage  <- intersect(tr_ids, te_ids)
    if (length(leakage) > 0) {
      cat("⚠️  Fold", i, ": LEAKAGE détecté ! Triplets:", leakage, "\n")
      has_leakage <- TRUE
    }
  }
  if (!has_leakage) cat("✅ Vérification leakage OK — aucun triplet splitté\n")
  
  cat("Nb folds:", nrow(folds), "\n")
  cat("Taille fold 1 — train:", nrow(training(folds$splits[[1]])),
      "| test:", nrow(testing(folds$splits[[1]])), "\n")
  folds
}

# Résumé comparatif d'un workflow_set
summarise_wf_results <- function(wf_results) {
  collect_metrics(wf_results) |>
    filter(.metric == "roc_auc") |>
    arrange(desc(mean)) |>
    select(wflow_id, mean, std_err, n)
}

metric <- metric_set(roc_auc)



# Recettes ----
# Recette complète : pour modèles gérant la haute dimension (glmnet, RF, XGB)
make_recipe_full <- function(data, outcome_var = "als") {
  recipe(as.formula(paste(outcome_var, "~ .")), data = data) |>
    update_role(match, new_role = "id") |>
    step_mutate(
      smoking_2cat_i = as.numeric(smoking_2cat_i),
      sex            = as.numeric(sex)) |>
    step_zv(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors())
}

# Recette filtrage top-N : pour modèles fragiles (SVM, GAM, bayesGLM, MARS)
# Utilise notre step custom basé sur pROC (100% CRAN)
make_recipe_top <- function(data, outcome_var = "als", n_top = 20) {
  recipe(as.formula(paste(outcome_var, "~ .")), data = data) |>
    update_role(match, new_role = "id") |>
    step_mutate(
      smoking_2cat_i = as.numeric(smoking_2cat_i),
      sex            = as.numeric(sex)) |>
    step_zv(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_select_auc_top(
      all_numeric_predictors(),
      outcome = outcome_var,
      top_n   = n_top)
}



# TEST 1 : Protéines + covariables (sans follow_up) ----
cat("\n========== TEST 1 : protéines + covariables ==========\n")

data_t1 <- bdd_danish |>
  select(als, birth_year, sex, bmi, smoking_2cat_i, match, all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("control", "case")))

cat("Dimensions :", nrow(data_t1), "×", ncol(data_t1) - 2, "prédicteurs\n")
cat("Outcome :", table(data_t1$als), "\n")

folds_t1 <- make_folds_checked(data_t1)

rec_full_t1  <- make_recipe_full(data_t1)
rec_top20_t1 <- make_recipe_top(data_t1, n_top = 20)
rec_top10_t1 <- make_recipe_top(data_t1, n_top = 10)

ncol_preds_t1 <- ncol(data_t1) - 3  # -als -match + les covariables



## TEST 1 – PHASE A : Tuning hyperparamètres ----

### A1. glmnet ----

glmnet_spec_tune <- logistic_reg(penalty = tune(), mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t1 <- workflow() |>
  add_recipe(rec_full_t1) |>
  add_model(glmnet_spec_tune)

set.seed(1996)
tune_glmnet_t1 <- tune_bayes(
  wf_glmnet_t1,
  resamples  = folds_t1,
  param_info = parameters(
    dials::penalty(range = c(-5, 0)),   # log10(lambda) : 1e-5 à 1
    mixture(range = c(0, 1))),   # alpha : Ridge (0) → Lasso (1)
  iter       = 40,               # 40 itérations Bayésiennes >> 140 grille
  initial    = 10,               # 10 pts aléatoires pour amorcer le GP
  metrics    = metric,
  control    = control_bayes(save_pred = TRUE, verbose = TRUE))

autoplot(tune_glmnet_t1, type = "performance") +
  theme_minimal() + ggtitle("Test 1 – glmnet (tune_bayes)")
show_best(tune_glmnet_t1, metric = "roc_auc", n = 5)
best_glmnet_t1 <- select_best(tune_glmnet_t1, metric = "roc_auc")
cat("Best glmnet T1 — mixture (alpha):", round(best_glmnet_t1$mixture, 3),
    "| penalty (lambda):", round(best_glmnet_t1$penalty, 5), "\n")


### A2. Random Forest ----

rf_spec_tune <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
  set_engine("ranger",
             num.threads  = parallel::detectCores() - 1,
             importance   = "permutation") |>  # pour vip() plus tard
  set_mode("classification")

tune_rf_results_t1 <- map(c("gini", "extratrees"), function(sr) {
  
  rf_spec_sr <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
    set_engine("ranger", splitrule = sr,
               num.threads = parallel::detectCores() - 1,
               importance  = "permutation") |>
    set_mode("classification")
  
  wf_rf_sr <- workflow() |> add_recipe(rec_full_t1) |> add_model(rf_spec_sr)
  
  set.seed(1996)
  res <- tune_bayes(
    wf_rf_sr,
    resamples  = folds_t1,
    param_info = parameters(
      mtry(range = c(5, min(30, ncol_preds_t1))),
      min_n(range = c(1, 20))),
    iter       = 25,
    initial    = 8,
    metrics    = metric,
    control    = control_bayes(save_pred = TRUE, verbose = TRUE))
  
  list(splitrule = sr, tuned = res, best = select_best(res, metric = "roc_auc"))
})

# Meilleur splitrule global
best_rf_auc_t1 <- map_dbl(tune_rf_results_t1, ~ show_best(.x$tuned, metric = "roc_auc", n=1)$mean)
best_rf_idx_t1 <- which.max(best_rf_auc_t1)
best_rf_t1     <- tune_rf_results_t1[[best_rf_idx_t1]]
cat("Best RF T1 — splitrule:", best_rf_t1$splitrule,
    "| mtry:", best_rf_t1$best$mtry,
    "| min_n:", best_rf_t1$best$min_n, "\n")


### A3. XGBoost : grid_space_filling + early stopping ----
# Phase 1 : tune_grid sur hyperparamètres structurels (depth, eta, mtry, min_n, gamma)
#           subsample fixé à 0.8 
# Phase 2 : early stopping via xgb.cv pour trouver le nombre optimal d'arbres

xgb_spec_tune <- boost_tree(
  trees          = 1000,
  tree_depth     = tune(),
  learn_rate     = tune(),
  mtry           = tune(),
  min_n          = tune(),
  loss_reduction = tune(),
  sample_size    = 0.8       # ← fixé ici, pas dans set_engine()
) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

# grid_space_filling() = nouveau nom de grid_latin_hypercube() depuis dials 1.3.0
# 50 points LHS couvrent mieux l'espace 5D qu'une grille régulière 3^5 = 243 pts
xgb_grid_t1 <- grid_space_filling(
  tree_depth(range     = c(1, 4)),
  learn_rate(range     = c(-3, -1)),  # 0.001 à 0.1 en log10
  mtry(range           = c(floor(ncol_preds_t1 * 0.1),
                           floor(ncol_preds_t1 * 0.6))),
  min_n(range          = c(1, 20)),
  loss_reduction(range = c(-5, 0)),   # gamma
  size = 50)

wf_xgb_t1 <- workflow() |> add_recipe(rec_full_t1) |> add_model(xgb_spec_tune)

set.seed(1996)
tune_xgb_t1 <- tune_grid(
  wf_xgb_t1,
  resamples = folds_t1,
  grid      = xgb_grid_t1,
  metrics   = metric,
  control   = control_grid(save_pred = TRUE, verbose = TRUE))

autoplot(tune_xgb_t1) + theme_minimal() + ggtitle("Test 1 – XGBoost")
show_best(tune_xgb_t1, metric = "roc_auc", n = 5)
best_xgb_params_t1 <- select_best(tune_xgb_t1, metric = "roc_auc")
cat("Best XGB structural params T1:\n"); print(best_xgb_params_t1)

# Early stopping XGBoost 
# On reprend les meilleurs hyperparamètres structurels et on laisse xgb.cv trouver le nb optimal d'arbres via early stopping

data_t1_model <- data_t1 |>
  mutate(
    als_num        = as.numeric(als == "case"),
    smoking_2cat_i = as.numeric(smoking_2cat_i),
    sex            = as.numeric(sex)) |>
  select(-als, -match)

X_mat_t1 <- data_t1_model |> select(-als_num) |> as.matrix()
Y_t1      <- data_t1$als
Y_num_t1  <- as.numeric(Y_t1 == "case")

dtrain_t1 <- xgb.DMatrix(data = X_mat_t1, label = Y_num_t1)

# Folds xgb.cv (indices 1-based → correspondance avec folds_t1)
xgb_folds_t1 <- lapply(seq_len(nrow(folds_t1)), function(i) {
  # Récupérer les match_ids du fold test
  test_match_ids <- testing(folds_t1$splits[[i]])$match
  which(data_t1$match %in% test_match_ids)
})

set.seed(1996)
xgb_cv_es_t1 <- xgb.cv(
  params = list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = best_xgb_params_t1$tree_depth,
    eta              = best_xgb_params_t1$learn_rate,
    colsample_bytree = best_xgb_params_t1$mtry / ncol_preds_t1,
    min_child_weight = best_xgb_params_t1$min_n,
    gamma            = best_xgb_params_t1$loss_reduction,
    subsample        = 0.8,
    nthread          = parallel::detectCores() - 1),
  data                  = dtrain_t1,
  nrounds               = 5000,
  folds                 = xgb_folds_t1,
  early_stopping_rounds = 50,    # plus conservateur que 30 (données petites)
  maximize              = TRUE,
  verbose               = 1,
  print_every_n         = 100)

best_nrounds_t1 <- which.max(xgb_cv_es_t1$evaluation_log$test_auc_mean)
best_row_t1     <- xgb_cv_es_t1$evaluation_log[best_nrounds_t1, ]
cat("\n✅ Early stopping T1 — meilleur nrounds:", best_nrounds_t1, "\n")
cat("AUC CV (mean ± SD):", round(best_row_t1$test_auc_mean, 3),
    "±", round(best_row_t1$test_auc_std, 3), "\n")

# Courbe d'apprentissage
xgb_cv_es_t1$evaluation_log |>
  ggplot(aes(x = iter)) +
  geom_line(aes(y = train_auc_mean, color = "Train")) +
  geom_line(aes(y = test_auc_mean,  color = "CV Test")) +
  geom_ribbon(aes(ymin = test_auc_mean - test_auc_std,
                  ymax = test_auc_mean + test_auc_std), alpha = 0.15) +
  geom_vline(xintercept = best_nrounds_t1, linetype = "dashed", color = "grey40") +
  labs(title = "Learning curve XGBoost – Test 1",
       x = "Nb d'arbres", y = "AUC", color = "") +
  theme_minimal()


### A4. SVM linéaire + top20 ----

svm_linear_spec_tune <- svm_linear(cost = tune()) |>
  set_engine("kernlab") |>
  set_mode("classification")

wf_svm_t1 <- workflow() |> add_recipe(rec_top20_t1) |> add_model(svm_linear_spec_tune)

set.seed(1996)
tune_svm_t1 <- tune_grid(
  wf_svm_t1,
  resamples = folds_t1,
  grid      = grid_regular(cost(range = c(-2, 2)), levels = 7),
  metrics   = metric,
  control   = control_grid(save_pred = TRUE))

show_best(tune_svm_t1, metric = "roc_auc", n = 5)
best_svm_t1 <- select_best(tune_svm_t1, metric = "roc_auc")
cat("Best SVM T1 — cost:", round(best_svm_t1$cost, 4), "\n")


### A5. MARS + top20 ----
mars_spec_tune <- mars(num_terms = tune(), prod_degree = tune()) |>
  set_engine("earth") |>
  set_mode("classification")

wf_mars_t1 <- workflow() |> add_recipe(rec_top20_t1) |> add_model(mars_spec_tune)

set.seed(1996)
tune_mars_t1 <- tune_grid(
  wf_mars_t1,
  resamples = folds_t1,
  grid      = expand_grid(
    num_terms   = c(5, 10, 20, 50),
    prod_degree = c(1L, 2L)),
  metrics   = metric,
  control   = control_grid(save_pred = TRUE))

show_best(tune_mars_t1, metric = "roc_auc", n = 5)
best_mars_t1 <- select_best(tune_mars_t1, metric = "roc_auc")



# TEST 1 – PHASE B : Comparaison des algorithmes entre eux  ----

glmnet_final_t1 <- logistic_reg(
  penalty = best_glmnet_t1$penalty,
  mixture = best_glmnet_t1$mixture) |>
  set_engine("glmnet") |> set_mode("classification")

rf_final_t1 <- rand_forest(
  mtry  = best_rf_t1$best$mtry,
  min_n = best_rf_t1$best$min_n,
  trees = 1000) |>
  set_engine("ranger",
             splitrule   = best_rf_t1$splitrule,
             importance  = "permutation",
             num.threads = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_final_spec_t1 <- boost_tree(
  trees        = best_nrounds_t1,
  tree_depth   = best_xgb_params_t1$tree_depth,
  learn_rate   = best_xgb_params_t1$learn_rate,
  mtry         = best_xgb_params_t1$mtry,
  min_n        = best_xgb_params_t1$min_n,
  loss_reduction = best_xgb_params_t1$loss_reduction,
  sample_size  = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

svm_final_t1 <- svm_linear(cost = best_svm_t1$cost) |>
  set_engine("kernlab") |> set_mode("classification")

mars_final_t1 <- mars(
  num_terms   = best_mars_t1$num_terms,
  prod_degree = best_mars_t1$prod_degree) |>
  set_engine("earth") |> set_mode("classification")

ridge_bayes_final <- logistic_reg(penalty = 0.05, mixture = 0) |>
  set_engine("glmnet") |> set_mode("classification")

gam_final <- logistic_reg(penalty = 0.01, mixture = 0) |>
  set_engine("glmnet") |>
  set_mode("classification")

null_final <- null_model() |>
  set_engine("parsnip") |> set_mode("classification")

# Construction du workflow_set
wf_set_t1 <- bind_rows(
  workflow_set(preproc = list(full = rec_full_t1),
               models = list(glmnet  = glmnet_final_t1,
                             rf      = rf_final_t1,
                             xgboost = xgb_final_spec_t1,
                             null    = null_final)),
  workflow_set(preproc = list(top20 = rec_top20_t1),
               models = list(svm         = svm_final_t1,
                             mars        = mars_final_t1,
                             ridge_bayes = ridge_bayes_final)),
  workflow_set(preproc = list(top10 = rec_top10_t1),
               models = list(ridge_top10 = gam_final)))


set.seed(1996)
results_t1 <- workflow_map(
  wf_set_t1,
  fn        = "fit_resamples",
  resamples = folds_t1,
  metrics   = metric,
  control   = control_resamples(save_pred = TRUE),
  verbose   = TRUE)

summary_t1 <- summarise_wf_results(results_t1)
cat("\n=== Résultats Test 1 ===\n")
print(summary_t1)

autoplot(results_t1, metric = "roc_auc") +
  theme_minimal() +
  geom_text(aes(label = round(mean, 3)), hjust = -0.1, size = 3) +
  labs(title = "Test 1 – Comparaison des modèles (AUC, 10-fold groupé)",
       x = "AUC moyen (CV)", y = "")

# Meilleur modèle Test 1
best_wf_id_t1 <- summary_t1$wflow_id[1]
cat("✅ Meilleur modèle Test 1 :", best_wf_id_t1, "\n")
cat("   AUC moyen :", round(summary_t1$mean[1], 3),
    "± SE:", round(summary_t1$std_err[1], 3), "\n")



# TEST 2 : Protéines + covariables + follow_up_no_na_y ----

cat("\n========== TEST 2 : + follow_up ==========\n")

data_t2 <- bdd_danish |>
  select(als, follow_up_no_na_y, birth_year, sex, bmi, smoking_2cat_i,
         match, all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("control", "case")))

folds_t2     <- make_folds_checked(data_t2)
ncol_preds_t2 <- ncol(data_t2) - 3

rec_full_t2  <- make_recipe_full(data_t2)
rec_top20_t2 <- make_recipe_top(data_t2, n_top = 20)
rec_top10_t2 <- make_recipe_top(data_t2, n_top = 10)

## TEST 2 – PHASE A : Tuning hyperparamètres ----

### A1. glmnet ----
wf_glmnet_t2 <- workflow() |> add_recipe(rec_full_t2) |> add_model(glmnet_spec_tune)
set.seed(1996)
tune_glmnet_t2 <- tune_bayes(wf_glmnet_t2, resamples = folds_t2,
                             param_info = parameters(dials::penalty(range = c(-5, 0)), mixture(range = c(0, 1))),
                             iter = 40, initial = 10, metrics = metric,
                             control = control_bayes(save_pred = TRUE, verbose = TRUE))
best_glmnet_t2 <- select_best(tune_glmnet_t2, metric = "roc_auc")

### A2. Random Forest ----
tune_rf_results_t2 <- map(c("gini", "extratrees"), function(sr) {
  rf_spec_sr <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
    set_engine("ranger", splitrule = sr, num.threads = parallel::detectCores() - 1,
               importance = "permutation") |>
    set_mode("classification")
  wf_rf_sr <- workflow() |> add_recipe(rec_full_t2) |> add_model(rf_spec_sr)
  set.seed(1996)
  res <- tune_bayes(wf_rf_sr, resamples = folds_t2,
                    param_info = parameters(mtry(range = c(5, min(30, ncol_preds_t2))),
                                            min_n(range = c(1, 20))),
                    iter = 25, initial = 8, metrics = metric,
                    control = control_bayes(save_pred = TRUE, verbose = TRUE))
  list(splitrule = sr, tuned = res, best = select_best(res, metric = "roc_auc"))
})
best_rf_auc_t2 <- map_dbl(tune_rf_results_t2, ~ show_best(.x$tuned, metric = "roc_auc", n=1)$mean)
best_rf_t2     <- tune_rf_results_t2[[which.max(best_rf_auc_t2)]]

### A3. XGBoost : grid_space_filling + early stopping ----
xgb_grid_t2 <- grid_space_filling(
  tree_depth(range     = c(1, 4)), learn_rate(range = c(-3, -1)),
  mtry(range           = c(floor(ncol_preds_t2 * 0.1), floor(ncol_preds_t2 * 0.6))),
  min_n(range          = c(1, 20)), loss_reduction(range = c(-5, 0)), size = 50)

wf_xgb_t2 <- workflow() |> add_recipe(rec_full_t2) |> add_model(xgb_spec_tune)
set.seed(1996)
tune_xgb_t2 <- tune_grid(wf_xgb_t2, resamples = folds_t2, grid = xgb_grid_t2,
                         metrics = metric, control = control_grid(save_pred = TRUE, verbose = TRUE))
best_xgb_params_t2 <- select_best(tune_xgb_t2, metric = "roc_auc")

# Early stopping T2
data_t2_model <- data_t2 |>
  mutate(als_num = as.numeric(als == "case"),
         smoking_2cat_i = as.numeric(smoking_2cat_i), sex = as.numeric(sex)) |>
  select(-als, -match)

X_mat_t2  <- data_t2_model |> select(-als_num) |> as.matrix()
Y_num_t2  <- as.numeric(data_t2$als == "case")
dtrain_t2 <- xgb.DMatrix(data = X_mat_t2, label = Y_num_t2)

xgb_folds_t2 <- lapply(seq_len(nrow(folds_t2)), function(i) {
  test_match_ids <- testing(folds_t2$splits[[i]])$match
  which(data_t2$match %in% test_match_ids)
})

set.seed(1996)
xgb_cv_es_t2 <- xgb.cv(
  params = list(objective = "binary:logistic", eval_metric = "auc",
                max_depth = best_xgb_params_t2$tree_depth, eta = best_xgb_params_t2$learn_rate,
                colsample_bytree = best_xgb_params_t2$mtry / ncol_preds_t2,
                min_child_weight = best_xgb_params_t2$min_n,
                gamma = best_xgb_params_t2$loss_reduction, subsample = 0.8,
                nthread = parallel::detectCores() - 1),
  data = dtrain_t2, nrounds = 5000, folds = xgb_folds_t2,
  early_stopping_rounds = 50, maximize = TRUE, verbose = 1, print_every_n = 100)

best_nrounds_t2 <- which.max(xgb_cv_es_t2$evaluation_log$test_auc_mean)
cat("✅ Early stopping T2 — meilleur nrounds:", best_nrounds_t2, "\n")

### A4. SVM linéaire + top20 ----
wf_svm_t2 <- workflow() |> add_recipe(rec_top20_t2) |> add_model(svm_linear_spec_tune)
set.seed(1996)
tune_svm_t2 <- tune_grid(wf_svm_t2, resamples = folds_t2,
                         grid = grid_regular(cost(range = c(-2, 2)), levels = 7),
                         metrics = metric, control = control_grid(save_pred = TRUE))
best_svm_t2 <- select_best(tune_svm_t2, metric = "roc_auc")

### A5. MARS + top20 ----
wf_mars_t2 <- workflow() |> add_recipe(rec_top20_t2) |> add_model(mars_spec_tune)
set.seed(1996)
tune_mars_t2 <- tune_grid(wf_mars_t2, resamples = folds_t2,
                          grid = expand_grid(num_terms = c(5, 10, 20, 50), prod_degree = c(1L, 2L)),
                          metrics = metric, control = control_grid(save_pred = TRUE))
best_mars_t2 <- select_best(tune_mars_t2, metric = "roc_auc")


## TEST 2 – PHASE B : Comparaison des algorithmes entre eux  ----

glmnet_final_t2 <- logistic_reg(
  penalty = best_glmnet_t2$penalty, mixture = best_glmnet_t2$mixture) |>
  set_engine("glmnet") |> set_mode("classification")

rf_final_t2 <- rand_forest(
  mtry = best_rf_t2$best$mtry, min_n = best_rf_t2$best$min_n, trees = 1000) |>
  set_engine("ranger", splitrule = best_rf_t2$splitrule,
             importance = "permutation", num.threads = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_final_spec_t2 <- boost_tree(
  trees = best_nrounds_t2, tree_depth = best_xgb_params_t2$tree_depth,
  learn_rate = best_xgb_params_t2$learn_rate, mtry = best_xgb_params_t2$mtry,
  min_n = best_xgb_params_t2$min_n, loss_reduction = best_xgb_params_t2$loss_reduction) |>
  set_engine("xgboost", subsample = 0.8, nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

svm_final_t2  <- svm_linear(cost = best_svm_t2$cost) |>
  set_engine("kernlab") |> set_mode("classification")
mars_final_t2 <- mars(num_terms = best_mars_t2$num_terms,
                      prod_degree = best_mars_t2$prod_degree) |>
  set_engine("earth") |> set_mode("classification")

wf_set_t2 <- bind_rows(
  workflow_set(preproc = list(full = rec_full_t2),
               models = list(glmnet = glmnet_final_t2, rf = rf_final_t2,
                             xgboost = xgb_final_spec_t2, null = null_final)),
  workflow_set(preproc = list(top20 = rec_top20_t2),
               models = list(svm = svm_final_t2, mars = mars_final_t2,
                             ridge_bayes = ridge_bayes_final)),
  workflow_set(preproc = list(top10 = rec_top10_t2),
               models = list(gam = gam_final)))

set.seed(1996)
results_t2 <- workflow_map(wf_set_t2, fn = "fit_resamples", resamples = folds_t2,
                           metrics = metric, control = control_resamples(save_pred = TRUE), verbose = TRUE)

summary_t2 <- summarise_wf_results(results_t2)
cat("\n=== Résultats Test 2 ===\n")
print(summary_t2)

autoplot(results_t2, metric = "roc_auc") +
  theme_minimal() +
  labs(title = "Test 2 – Comparaison des modèles (AUC, 10-fold groupé)", x = "AUC moyen")

best_wf_id_t2 <- summary_t2$wflow_id[1]
cat("✅ Meilleur modèle Test 2 :", best_wf_id_t2, "\n")



# Analyse détaillée des meilleurs modèles ----
# On entraîne le meilleur modèle sur TOUTES les données (fit final)
# puis on calcule : SHAP, importance, calibration, courbe ROC, 
# stabilité des prédictions par fold

cat("\n========== PHASE C : Analyse détaillée ==========\n")


# C-helper : entraîner le workflow gagnant sur toutes les données ─────────
fit_best_workflow <- function(wf_set_results, wf_id, data) {
  best_wf <- extract_workflow(wf_set_results, id = wf_id)
  fit(best_wf, data = data)
}


## Analyse du meilleur modèle TEST 1 ----
cat("\n--- Meilleur modèle Test 1 :", best_wf_id_t1, "---\n")

best_fit_t1 <- fit_best_workflow(results_t1, best_wf_id_t1, data_t1)


# Alternative pour récupérer les prédictions OOF selon tidymodels version :
oof_preds_t1 <- results_t1 |>
  extract_workflow_set_result(best_wf_id_t1) |>
  collect_predictions()

### courbe ROC ----
roc_data_t1 <- roc_curve(oof_preds_t1, truth = als, .pred_case)
autoplot(roc_data_t1) +
  labs(title = paste("ROC – Test 1 –", best_wf_id_t1),
       subtitle = paste("AUC OOF:", round(summary_t1$mean[1], 3))) +
  theme_minimal()

### Calibration des probabilités ----
# Est-ce que le modèle prédit 30% de risque → vraiment ~30% de cas ?
cal_plot_windowed(oof_preds_t1, truth = als, estimate = .pred_case,
                  window_size = 0.2, step_size = 0.025) +
  labs(title = paste("Calibration –", best_wf_id_t1)) +
  theme_minimal()

### Variable importance (applicable si RF ou XGBoost est le meilleur) ----
if (grepl("rf|xgboost", best_wf_id_t1)) {
  
  best_model_t1 <- extract_fit_parsnip(best_fit_t1)
  
  # VIP basique
  vip(best_model_t1, num_features = 25, geom = "col") +
    labs(title = paste("Variable Importance –", best_wf_id_t1)) +
    theme_minimal()
  
  # --- SHAP values (XGBoost) 
  if (grepl("xgboost", best_wf_id_t1)) {
    
    xgb_model_t1 <- extract_fit_engine(best_fit_t1)
    
    # Récupérer la matrice X preprocessée (après recette)
    X_prepped_t1 <- rec_full_t1 |>
      prep(training = data_t1) |>
      bake(new_data = data_t1) |>
      select(-als, -match) |>
      as.matrix()
    
    shap_t1 <- shapviz(xgb_model_t1, X_pred = X_prepped_t1)
    
    # Beeswarm : direction + magnitude pour chaque individu
    sv_importance(shap_t1, kind = "beeswarm", max_display = 25) +
      labs(title = "SHAP – Test 1 – Beeswarm (top 25)")
    
    # Bar : moyenne |SHAP| par variable
    sv_importance(shap_t1, kind = "bar", max_display = 25) +
      labs(title = "SHAP – Test 1 – Importance moyenne (top 25)")
    
    # Résumé SHAP tabulaire (protéines uniquement)
    shap_matrix_t1 <- shap_t1$S
    shap_summary_t1 <- tibble(
      Feature       = colnames(shap_matrix_t1),
      mean_abs_shap = colMeans(abs(shap_matrix_t1)),
      mean_shap     = colMeans(shap_matrix_t1)) |>
      filter(Feature %in% proteomic) |>
      arrange(desc(mean_abs_shap)) |>
      mutate(direction = case_when(
        mean_shap > 0  ~ "↑ risque ALS",
        mean_shap < 0  ~ "↓ risque ALS",
        TRUE           ~ "neutre"))
    
    cat("\nTop 20 protéines (SHAP) – Test 1 :\n")
    print(shap_summary_t1, n = 20)
    
    # Dependence plots top 6
    top6_t1 <- shap_summary_t1$Feature[1:6]
    for (prot in top6_t1) {
      print(sv_dependence(shap_t1, v = prot) +
              labs(title = paste("SHAP dependence –", prot)))
    }
    
    # Waterfall pour le 1er cas ALS (diagnostic individuel)
    idx_als_t1 <- which(data_t1$als == "case")[1]
    sv_waterfall(shap_t1, row_id = idx_als_t1) +
      labs(title = paste("Waterfall – 1er cas ALS, Test 1"))
  }
  
  # --- SHAP values (RF via treeshap) 
  if (grepl("^full_rf|^full_rf", best_wf_id_t1)) {
    # Pour Random Forest, shapviz peut utiliser treeshap
    # Nécessite le package treeshap (CRAN)
    # library(treeshap)
    # rf_engine_t1 <- extract_fit_engine(best_fit_t1)
    # unified_t1   <- ranger.unify(rf_engine_t1, X_prepped_t1)
    # shap_rf_t1   <- treeshap(unified_t1, X_prepped_t1)
    # shap_t1      <- shapviz(shap_rf_t1)
    # sv_importance(shap_t1, kind = "beeswarm", max_display = 25)
    cat("Note: pour SHAP RF, décommenter le bloc treeshap ci-dessus\n")
  }
}


##  Analyse du meilleur modèle TEST 2 ----
cat("\n--- Meilleur modèle Test 2 :", best_wf_id_t2, "---\n")

best_fit_t2 <- fit_best_workflow(results_t2, best_wf_id_t2, data_t2)

oof_preds_t2 <- results_t2 |>
  extract_workflow_set_result(best_wf_id_t2) |>
  collect_predictions()

### Courbe ROC ----
roc_data_t2 <- roc_curve(oof_preds_t2, truth = als, .pred_case)
autoplot(roc_data_t2) +
  labs(title = paste("ROC – Test 2 –", best_wf_id_t2),
       subtitle = paste("AUC OOF:", round(summary_t2$mean[1], 3))) +
  theme_minimal()

### Calibration ----
cal_plot_windowed(oof_preds_t2, truth = als, estimate = .pred_case,
                  window_size = 0.2, step_size = 0.025) +
  labs(title = paste("Calibration –", best_wf_id_t2)) +
  theme_minimal()

if (grepl("xgboost", best_wf_id_t2)) {
  xgb_model_t2 <- extract_fit_engine(best_fit_t2)
  X_prepped_t2 <- rec_full_t2 |>
    prep(training = data_t2) |>
    bake(new_data = data_t2) |>
    select(-als, -match) |>
    as.matrix()
  
  shap_t2 <- shapviz(xgb_model_t2, X_pred = X_prepped_t2)
  sv_importance(shap_t2, kind = "beeswarm", max_display = 25) +
    labs(title = "SHAP – Test 2 – Beeswarm")
  sv_importance(shap_t2, kind = "bar", max_display = 25) +
    labs(title = "SHAP – Test 2 – Importance moyenne")
  
  shap_matrix_t2 <- shap_t2$S
  shap_summary_t2 <- tibble(
    Feature       = colnames(shap_matrix_t2),
    mean_abs_shap = colMeans(abs(shap_matrix_t2)),
    mean_shap     = colMeans(shap_matrix_t2)) |>
    filter(Feature %in% proteomic) |>
    arrange(desc(mean_abs_shap)) |>
    mutate(direction = case_when(
      mean_shap > 0 ~ "↑ risque ALS",
      mean_shap < 0 ~ "↓ risque ALS",
      TRUE          ~ "neutre"))
  
  cat("\nTop 20 protéines (SHAP) – Test 2 :\n")
  print(shap_summary_t2, n = 20)
  
  # Analyse spécifique Test 2 : rôle du follow_up dans les prédictions
  # Est-ce que la contribution SHAP de follow_up varie selon le délai ?
  if ("follow_up_no_na_y" %in% colnames(shap_matrix_t2)) {
    sv_dependence(shap_t2, v = "follow_up_no_na_y") +
      labs(title = "SHAP – Contribution du temps de suivi",
           x = "follow_up (années)", y = "Valeur SHAP") +
      theme_minimal()
    
    # Interaction NEFL × follow_up (si NEFL est une variable clé)
    if ("proteomic_neuro_explo_NEFL" %in% colnames(shap_matrix_t2)) {
      sv_dependence(shap_t2, v = "proteomic_neuro_explo_NEFL",
                    color_var = "follow_up_no_na_y") +
        labs(title = "SHAP NEFL × follow_up – Test 2") +
        theme_minimal()
    }
    
    # Interaction globale entre les SHAP values
    sv_interaction(shap_t2, max_display = 10) +
      labs(title = "SHAP interactions – Test 2 (top 10)")
  }
  
  # Waterfall 1er cas ALS
  idx_als_t2 <- which(data_t2$als == "case")[1]
  sv_waterfall(shap_t2, row_id = idx_als_t2) +
    labs(title = "Waterfall – 1er cas ALS, Test 2")
}



# SAUVEGARDE ----

results_proteomic_ALS_occurrence_tidymodels <- list(
  
  test_1 = list(
    data             = data_t1,
    folds            = folds_t1,
    # Tuning Phase A
    tune_glmnet      = tune_glmnet_t1,
    tune_rf          = tune_rf_results_t1,
    tune_xgb         = tune_xgb_t1,
    xgb_cv_early_stop = xgb_cv_es_t1,
    tune_svm         = tune_svm_t1,
    tune_mars        = tune_mars_t1,
    # Hyperparamètres optimaux
    best_glmnet      = best_glmnet_t1,
    best_rf          = best_rf_t1,
    best_xgb_params  = best_xgb_params_t1,
    best_nrounds     = best_nrounds_t1,
    best_svm         = best_svm_t1,
    best_mars        = best_mars_t1,
    # Phase B : comparaison finale
    wf_set_results   = results_t1,
    summary          = summary_t1,
    best_model_id    = best_wf_id_t1,
    best_fit         = best_fit_t1),
  
  test_2 = list(
    data             = data_t2,
    folds            = folds_t2,
    tune_glmnet      = tune_glmnet_t2,
    tune_rf          = tune_rf_results_t2,
    tune_xgb         = tune_xgb_t2,
    xgb_cv_early_stop = xgb_cv_es_t2,
    tune_svm         = tune_svm_t2,
    tune_mars        = tune_mars_t2,
    best_glmnet      = best_glmnet_t2,
    best_rf          = best_rf_t2,
    best_xgb_params  = best_xgb_params_t2,
    best_nrounds     = best_nrounds_t2,
    best_svm         = best_svm_t2,
    best_mars        = best_mars_t2,
    wf_set_results   = results_t2,
    summary          = summary_t2,
    best_model_id    = best_wf_id_t2,
    best_fit         = best_fit_t2)
)

saveRDS(
  results_proteomic_ALS_occurrence_tidymodels,
  file = "~/Documents/POP_ALS_2025_02_03/2_output/2.6.3_results_proteomic_ALS_occurrence_tidymodels.rds")

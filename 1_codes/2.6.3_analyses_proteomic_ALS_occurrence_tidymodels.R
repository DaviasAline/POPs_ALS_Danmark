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
library(vip)

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


bake.step_select_auc_top <- function(object, new_data, ...) {

  non_selected_preds <- setdiff(
    names(new_data)[sapply(new_data, is.numeric)],
    object$selected_vars)

  cols_to_keep <- setdiff(names(new_data), object$rejected_vars)
  new_data |> dplyr::select(dplyr::all_of(cols_to_keep))
}


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


print.step_select_auc_top <- function(x, width = max(20, options()$width - 35), ...) {
  cat("AUC-based top-", x$top_n, " feature selection\n", sep = "")
  invisible(x)
}





# Folds groupés sur `match` + vérification anti-leakage
make_folds_checked <- function(data, v = 10, seed = 1996) {
  data_with_match <- data  # match est déjà dans data
  
  set.seed(seed)
  folds <- group_vfold_cv(data_with_match, group = match, v = v)    # garantit qu'un même triplet se retrouve toujours dans le meme fold et qu'un cas n'est pas separé de ses témoins
  
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
  if (!has_leakage) cat("Vérification leakage OK — aucun triplet splitté\n")
  
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

# C-helper : entraîner le workflow gagnant sur toutes les données 

fit_best_workflow <- function(wf_set_results, wf_id, data, seed = 1996) {
  best_wf <- extract_workflow(wf_set_results, id = wf_id)
  set.seed(seed)
  fit(best_wf, data = data)
}


# Fonction pour extraire les coefficients glmnet d'un fit tidymodels
extract_glmnet_coefs <- function(best_fit) {
  tryCatch({
    glmnet_fit <- best_fit |> extract_fit_engine()
    lambda_opt <- best_fit |>
      extract_fit_parsnip() |>
      pluck("spec", "args", "penalty") |>
      as.numeric()
    
    # Extraire les coefficients de manière plus robuste
    coef_matrix <- coef(glmnet_fit, s = lambda_opt)
    
    # Convertir en tibble correctement
    coefs <- as_tibble(
      data.frame(
        Feature = rownames(coef_matrix),
        coefficient = as.numeric(coef_matrix[, 1]),
        stringsAsFactors = FALSE),
      rownames = NA) |>
      filter(Feature != "(Intercept)" & coefficient != 0) |>
      arrange(desc(abs(coefficient)))
    
    coefs_annotated <- coefs |>
      mutate(
        direction  = case_when(
          coefficient > 0 ~ "↑ ALS risk", 
          coefficient < 0 ~ "↓ ALS risk", 
          TRUE ~ "neutral"),
        abs_coef = abs(coefficient))
    
    return(list(
      success = TRUE, 
      coefs_all = coefs_annotated, 
      coefs_protein = coefs_annotated, 
      lambda_opt = lambda_opt, 
      n_selected = nrow(coefs_annotated)))
    
  }, error = function(e) {
    return(list(success = FALSE, error = as.character(e)))
  })
}

# Fonction pour extraire XGBoost et calculer SHAP
extract_xgb_and_shap <- function(fit_obj, data, test_label) {
  
  # Extraire le modèle XGBoost brut
  xgb_model <- fit_obj |> extract_fit_engine()
  
  # Préparer les données (standardisées via la recette)
  rec_prepped <- fit_obj |> extract_preprocessor() |> prep(data)
  X_mat <- bake(rec_prepped, new_data = data) |>
    select(-als, -match) |>
    as.matrix()
  
  # Calculer SHAP
  shap <- shapviz(xgb_model, X_pred = X_mat)
  shap_matrix <- shap$S
  
  # Résumé SHAP
  shap_summary <- tibble(
    Feature       = colnames(shap_matrix),
    mean_abs_shap = colMeans(abs(shap_matrix)),
    mean_shap     = colMeans(shap_matrix)) |>
    mutate(
      direction  = case_when(
        mean_shap > 0 ~ "↑ ALS risk",
        mean_shap < 0 ~ "↓ ALS risk",
        TRUE          ~ "neutral")) |>
    arrange(desc(mean_abs_shap))
  
  return(list(
    xgb_model      = xgb_model,
    shap           = shap,
    shap_matrix    = shap_matrix,
    shap_summary   = shap_summary,
    X_mat          = X_mat,
    test_label     = test_label))
}


# Recettes ----
# Full recipe (for high-dimensional models: glmnet, RF, XGB)
# Converts sex and smoking to numeric.
# Removes zero-variance predictors (step_zv).
# Standardizes all numeric predictors to mean = 0, SD = 1 (step_normalize).
make_recipe_full <- function(data, outcome_var = "als") {
  recipe(as.formula(paste(outcome_var, "~ .")), data = data) |>
    update_role(match, new_role = "id") |>
    step_mutate(
      smoking_2cat_i = as.numeric(smoking_2cat_i),
      sex            = as.numeric(sex)) |>
    step_zv(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors())
}

# Recette filtrage top-N (for fragile models: SVM, GAM, bayesGLM, MARS)
# Identical to rec_full, plus a custom step step_select_auc_top() that ranks all numeric predictors 
# by their univariate AUC (computed via pROC) 
# and retains only the top N. This is necessary for SVM and MARS, which become numerically unstable when p >> n .
# This step is computed on the training fold only during cross-validation, preventing data leakage.
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

data_t1 <- bdd_danish |>
  select(als, birth_year, sex, bmi, smoking_2cat_i, match, all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("control", "case")))

cat("Dimensions :", nrow(data_t1), "×", ncol(data_t1) - 2, "prédicteurs\n")
cat("Outcome :", table(data_t1$als), "\n")

folds_t1 <- make_folds_checked(data_t1)

rec_full_t1  <- make_recipe_full(data_t1)
rec_top20_t1 <- make_recipe_top(data_t1, n_top = 20)
#rec_top10_t1 <- make_recipe_top(data_t1, n_top = 10)

ncol_preds_t1 <- ncol(data_t1) - 2  # -als -match 



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
best_glmnet_t1 <- select_best(tune_glmnet_t1, metric = "roc_auc")
rm(glmnet_spec_tune, wf_glmnet_t1)


### A2. Random Forest ----
rf_spec_tune <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
  set_engine("ranger",
             num.threads  = parallel::detectCores() - 1,
             importance   = "permutation") |>  
  set_mode("classification")

tune_rf_t1 <- map(c("gini", "extratrees"), function(sr) {
  
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
best_rf_auc_t1 <- map_dbl(tune_rf_t1, ~ show_best(.x$tuned, metric = "roc_auc", n=1)$mean)
best_rf_idx_t1 <- which.max(best_rf_auc_t1)
best_rf_t1     <- tune_rf_t1[[best_rf_idx_t1]]
rm(rf_spec_tune, best_rf_auc_t1, best_rf_idx_t1)


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
  sample_size    = 0.8) |>
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
  loss_reduction(range = c(-5, 0)),   # loss_reduction (gamma) en log10
  size = 50)

wf_xgb_t1 <- workflow() |> add_recipe(rec_full_t1) |> add_model(xgb_spec_tune)

set.seed(1996)
tune_xgb_t1 <- tune_grid(
  wf_xgb_t1,
  resamples = folds_t1,
  grid      = xgb_grid_t1,
  metrics   = metric,
  control   = control_grid(save_pred = TRUE, verbose = TRUE))
best_xgb_params_t1 <- select_best(tune_xgb_t1, metric = "roc_auc")


# Early stopping XGBoost 
# On reprend les meilleurs hyperparamètres structurels et on laisse xgb.cv trouver le nb optimal d'arbres via early stopping

rec_prepped_t1 <- rec_full_t1 |> prep(training = data_t1)
X_mat_t1 <- bake(rec_prepped_t1, new_data = data_t1) |>
  select(-als, -match) |>
  as.matrix()
Y_num_t1 <- as.numeric(data_t1$als == "case")
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
    min_child_weight = 1,  # a revoir, equivalent a min_n de tidyverse mais pas vraiment
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

best_xgb_trees_t1 <- which.max(xgb_cv_es_t1$evaluation_log$test_auc_mean)
rm(xgb_spec_tune, xgb_grid_t1, wf_xgb_t1, rec_prepped_t1, X_mat_t1, Y_num_t1, dtrain_t1, xgb_folds_t1)


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

best_svm_t1 <- select_best(tune_svm_t1, metric = "roc_auc")
rm(svm_linear_spec_tune, wf_svm_t1)


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

best_mars_t1 <- select_best(tune_mars_t1, metric = "roc_auc")
rm(mars_spec_tune, wf_mars_t1)


# TEST 1 – PHASE B : Comparaison des algorithmes entre eux  ----
glmnet_final_t1 <- logistic_reg(
  penalty = best_glmnet_t1$penalty,
  mixture = best_glmnet_t1$mixture) |>
  set_engine("glmnet") |> 
  set_mode("classification")

rf_final_t1 <- rand_forest(
  mtry  = best_rf_t1$best$mtry,
  min_n = best_rf_t1$best$min_n,
  trees = 1000) |>
  set_engine("ranger",
             splitrule   = best_rf_t1$splitrule,
             importance  = "permutation",
             num.threads = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_final_t1 <- boost_tree(
  trees        = best_xgb_trees_t1,
  tree_depth   = best_xgb_params_t1$tree_depth,
  learn_rate   = best_xgb_params_t1$learn_rate,
  mtry         = best_xgb_params_t1$mtry,
  min_n        = best_xgb_params_t1$min_n,
  loss_reduction = best_xgb_params_t1$loss_reduction,
  sample_size  = 0.8) |>
  set_engine("xgboost", 
             nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

svm_final_t1 <- svm_linear(cost = best_svm_t1$cost) |>
  set_engine("kernlab") |> 
  set_mode("classification")

mars_final_t1 <- mars(
  num_terms   = best_mars_t1$num_terms,
  prod_degree = best_mars_t1$prod_degree) |>
  set_engine("earth") |> 
  set_mode("classification")

# Construction du workflow_set
wf_set_t1 <- bind_rows(
  workflow_set(preproc = list(full = rec_full_t1),
               models = list(glmnet  = glmnet_final_t1,
                             rf      = rf_final_t1,
                             xgboost = xgb_final_t1)),
  workflow_set(preproc = list(top20 = rec_top20_t1),
               models = list(svm         = svm_final_t1,
                             mars        = mars_final_t1)))


set.seed(1996)
results_t1 <- workflow_map(
  wf_set_t1,
  fn        = "fit_resamples",
  resamples = folds_t1,
  metrics   = metric,
  control   = control_resamples(save_pred = TRUE),
  verbose   = TRUE)
summary_t1 <- summarise_wf_results(results_t1)

best_wf_id_t1 <- summary_t1$wflow_id[1]      # full_glmnet
best_fit_t1 <- fit_best_workflow(results_t1, best_wf_id_t1, data_t1)

xgb_id_t1 <- summary_t1$wflow_id[5]        # full_xgboost
xgb_fit_t1 <- fit_best_workflow(results_t1, xgb_id_t1, data_t1)

rm(wf_set_t1)


# TEST 2 : Protéines + covariables + follow_up_no_na_y ----

data_t2 <- bdd_danish |>
  select(als, follow_up_no_na_y, birth_year, sex, bmi, smoking_2cat_i,
         match, all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("control", "case")))

folds_t2     <- make_folds_checked(data_t2)
ncol_preds_t2 <- ncol(data_t2) - 3

rec_full_t2  <- make_recipe_full(data_t2)
rec_top20_t2 <- make_recipe_top(data_t2, n_top = 20)
#rec_top10_t2 <- make_recipe_top(data_t2, n_top = 10)

## TEST 2 – PHASE A : Tuning hyperparamètres ----

### A1. glmnet ----
glmnet_spec_tune <- logistic_reg(penalty = tune(), mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t2 <- workflow() |> 
  add_recipe(rec_full_t2) |> 
  add_model(glmnet_spec_tune)

set.seed(1996)
tune_glmnet_t2 <- tune_bayes(wf_glmnet_t2, resamples = folds_t2,
                             param_info = parameters(dials::penalty(range = c(-5, 0)), 
                                                     mixture(range = c(0, 1))),
                             iter = 40, initial = 10, metrics = metric,
                             control = control_bayes(save_pred = TRUE, verbose = TRUE))
best_glmnet_t2 <- select_best(tune_glmnet_t2, metric = "roc_auc")

rm(glmnet_spec_tune, wf_glmnet_t2)


### A2. Random Forest ----
tune_rf_t2 <- map(c("gini", "extratrees"), function(sr) {
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
best_rf_auc_t2 <- map_dbl(tune_rf_t2, ~ show_best(.x$tuned, metric = "roc_auc", n=1)$mean)
best_rf_t2     <- tune_rf_t2[[which.max(best_rf_auc_t2)]]
rm(best_rf_auc_t2)


### A3. XGBoost : grid_space_filling + early stopping ----
xgb_spec_tune <- boost_tree(
  trees          = 1000,
  tree_depth     = tune(),
  learn_rate     = tune(),
  mtry           = tune(),
  min_n          = tune(),
  loss_reduction = tune(),
  sample_size    = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_grid_t2 <- grid_space_filling(
  tree_depth(range     = c(1, 4)), learn_rate(range = c(-3, -1)),
  mtry(range           = c(floor(ncol_preds_t2 * 0.1), floor(ncol_preds_t2 * 0.6))),
  min_n(range          = c(1, 20)), loss_reduction(range = c(-5, 0)), size = 50)

wf_xgb_t2 <- workflow() |> add_recipe(rec_full_t2) |> add_model(xgb_spec_tune)
set.seed(1996)
tune_xgb_t2 <- tune_grid(wf_xgb_t2, 
                         resamples = folds_t2, 
                         grid = xgb_grid_t2,
                         metrics = metric, 
                         control = control_grid(save_pred = TRUE, verbose = TRUE))
best_xgb_params_t2 <- select_best(tune_xgb_t2, metric = "roc_auc")
rm(xgb_spec_tune, xgb_grid_t2, wf_xgb_t2)

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
                min_child_weight = 1,  # A revoir 
                gamma = best_xgb_params_t2$loss_reduction, subsample = 0.8,
                nthread = parallel::detectCores() - 1),
  data = dtrain_t2, nrounds = 5000, folds = xgb_folds_t2,
  early_stopping_rounds = 50, maximize = TRUE, verbose = 1, print_every_n = 100)

best_xgb_trees_t2 <- which.max(xgb_cv_es_t2$evaluation_log$test_auc_mean)
cat("Early stopping T2 — meilleur nrounds:", best_xgb_trees_t2, "\n")

rm(data_t2_model, X_mat_t2, Y_num_t2, dtrain_t2, xgb_folds_t2)


### A4. SVM linéaire + top20 ----
svm_linear_spec_tune <- svm_linear(cost = tune()) |>
  set_engine("kernlab") |>
  set_mode("classification")

wf_svm_t2 <- workflow() |> 
  add_recipe(rec_top20_t2) |> 
  add_model(svm_linear_spec_tune)
set.seed(1996)
tune_svm_t2 <- tune_grid(wf_svm_t2, resamples = folds_t2,
                         grid = grid_regular(cost(range = c(-2, 2)), levels = 7),
                         metrics = metric, control = control_grid(save_pred = TRUE))
best_svm_t2 <- select_best(tune_svm_t2, metric = "roc_auc")
rm(svm_linear_spec_tune, wf_svm_t2)

### A5. MARS + top20 ----
mars_spec_tune <- mars(num_terms = tune(), prod_degree = tune()) |>
  set_engine("earth") |>
  set_mode("classification")

wf_mars_t2 <- workflow() |> 
  add_recipe(rec_top20_t2) |> 
  add_model(mars_spec_tune)
set.seed(1996)
tune_mars_t2 <- tune_grid(wf_mars_t2, 
                          resamples = folds_t2,
                          grid = expand_grid(num_terms = c(5, 10, 20, 50), 
                                             prod_degree = c(1L, 2L)),
                          metrics = metric, 
                          control = control_grid(save_pred = TRUE))
best_mars_t2 <- select_best(tune_mars_t2, metric = "roc_auc")
rm(mars_spec_tune, wf_mars_t2)


## TEST 2 – PHASE B : Comparaison des algorithmes entre eux  ----

glmnet_final_t2 <- logistic_reg(
  penalty = best_glmnet_t2$penalty, 
  mixture = best_glmnet_t2$mixture) |>
  set_engine("glmnet") |> 
  set_mode("classification")

rf_final_t2 <- rand_forest(
  mtry = best_rf_t2$best$mtry, 
  min_n = best_rf_t2$best$min_n, 
  trees = 1000) |>
  set_engine("ranger", 
             splitrule = best_rf_t2$splitrule,
             importance = "permutation", 
             num.threads = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_final_t2 <- boost_tree(
  trees = best_xgb_trees_t2, 
  tree_depth = best_xgb_params_t2$tree_depth,
  learn_rate = best_xgb_params_t2$learn_rate, 
  mtry = best_xgb_params_t2$mtry,
  min_n = best_xgb_params_t2$min_n, 
  loss_reduction = best_xgb_params_t2$loss_reduction, 
  sample_size = 0.8) |>
  set_engine("xgboost", 
             nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

xgb_tree_depth_2_t2 <- boost_tree(
  trees = best_xgb_trees_t2, 
  tree_depth = 2,
  learn_rate = best_xgb_params_t2$learn_rate, 
  mtry = best_xgb_params_t2$mtry,
  min_n = best_xgb_params_t2$min_n, 
  loss_reduction = best_xgb_params_t2$loss_reduction, 
  sample_size = 0.8) |>
  set_engine("xgboost", 
             nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

svm_final_t2  <- svm_linear(cost = best_svm_t2$cost) |>
  set_engine("kernlab") |> 
  set_mode("classification")

mars_final_t2 <- mars(num_terms = best_mars_t2$num_terms,
                      prod_degree = best_mars_t2$prod_degree) |>
  set_engine("earth") |> 
  set_mode("classification")

wf_set_t2 <- bind_rows(
  workflow_set(preproc = list(full = rec_full_t2),
               models = list(glmnet = glmnet_final_t2, 
                             rf = rf_final_t2,
                             xgboost = xgb_final_t2, 
                             xgboost_tree_depth_2 = xgb_tree_depth_2_t2)),
  workflow_set(preproc = list(top20 = rec_top20_t2),
               models = list(svm = svm_final_t2, 
                             mars = mars_final_t2)))

set.seed(1996)
results_t2 <- workflow_map(wf_set_t2, 
                           fn = "fit_resamples", 
                           resamples = folds_t2,
                           metrics = metric, 
                           control = control_resamples(save_pred = TRUE), 
                           verbose = TRUE)

summary_t2 <- summarise_wf_results(results_t2)

best_wf_id_t2 <- summary_t2$wflow_id[1]  #full_glmnet
best_fit_t2 <- fit_best_workflow(results_t2, best_wf_id_t2, data_t2)

xgb_id_t2 <- summary_t2$wflow_id[4] # full_xgboost
xgb_fit_t2 <- fit_best_workflow(results_t2, xgb_id_t2, data_t2)

xgb_tree2_id_t2 <- summary_t2$wflow_id[5]   # full_xgboost_tree_depth_2
xgb_tree2_fit_t2 <- fit_best_workflow(results_t2, xgb_tree2_id_t2, data_t2)

rm(wf_set_t2)


# Sauvegarde ----

results_proteomic_ALS_occurrence_tidymodels <- list(
  
  test_1 = list(
    #data_t1 = data_t1,
    #folds_t1 = folds_t1,
    # Tuning test 1
    tune_glmnet_t1 = tune_glmnet_t1,
    tune_rf_t1 = tune_rf_t1,
    tune_xgb_t1 = tune_xgb_t1,
    xgb_cv_es_t1 = xgb_cv_es_t1,
    tune_svm_t1 = tune_svm_t1,
    tune_mars_t1 = tune_mars_t1,
    # Hyperparamètres optimaux test 1
    best_glmnet_t1 = best_glmnet_t1,
    best_rf_t1 = best_rf_t1,
    best_xgb_params_t1 = best_xgb_params_t1,
    best_xgb_trees_t1 = best_xgb_trees_t1,
    best_svm_t1 = best_svm_t1,
    best_mars_t1  = best_mars_t1,
    # comparaison finale test 1 
    results_t1 = results_t1,
    summary_t1  = summary_t1,
    best_wf_id_t1 = best_wf_id_t1,
    best_fit_t1  = best_fit_t1, 
    xgb_id_t1 = xgb_id_t1, 
    xgb_fit_t1 = xgb_fit_t1, 
    glmnet_final_t1 = glmnet_final_t1, 
    rf_final_t1 = rf_final_t1, 
    xgb_final_t1 = xgb_final_t1, 
    svm_final_t1 = svm_final_t1, 
    mars_final_t1 = mars_final_t1),
  
  test_2 = list(
    #data_t2 = data_t2,
    #folds_t2 = folds_t2,
    # Tuning test 2 
    tune_glmnet_t2 = tune_glmnet_t2,
    tune_rf_t2 = tune_rf_t2,
    tune_xgb_t2 = tune_xgb_t2,
    xgb_cv_es_t2 = xgb_cv_es_t2,
    tune_svm_t2  = tune_svm_t2,
    tune_mars_t2  = tune_mars_t2,
    # Hyperparamètres optimaux
    best_glmnet_t2  = best_glmnet_t2,
    best_rf_t2  = best_rf_t2,
    best_xgb_params_t2  = best_xgb_params_t2,
    best_xgb_trees_t2  = best_xgb_trees_t2,
    best_svm_t2  = best_svm_t2,
    best_mars_t2  = best_mars_t2,
    # comparaison finale test 2
    results_t2 = results_t2,
    summary_t2  = summary_t2,
    best_wf_id_t2 = best_wf_id_t2,
    best_fit_t2  = best_fit_t2, 
    xgb_id_t2 = xgb_id_t2, 
    xgb_fit_t2 = xgb_fit_t2, 
    xgb_tree2_id_t2 = xgb_tree2_id_t2, 
    xgb_tree2_fit_t2 = xgb_tree2_fit_t2, 
    glmnet_final_t2 = glmnet_final_t2, 
    rf_final_t2 = rf_final_t2, 
    xgb_final_t2 = xgb_final_t2, 
    svm_final_t2 = svm_final_t2, 
    mars_final_t2 = mars_final_t2))


rm(tune_glmnet_t1,
   tune_rf_t1,
   tune_xgb_t1,
   xgb_cv_es_t1,
   tune_svm_t1,
   tune_mars_t1,
   best_glmnet_t1,
   best_rf_t1,
   best_xgb_params_t1,
   best_xgb_trees_t1,
   best_svm_t1,
   best_mars_t1,
   glmnet_final_t1, 
   rf_final_t1, 
   xgb_final_t1, 
   svm_final_t1, 
   mars_final_t1, 
   results_t1,
   summary_t1,
   best_wf_id_t1,
   best_fit_t1,
   xgb_id_t1, 
   xgb_fit_t1, 
   
   tune_glmnet_t2,
   tune_rf_t2,
   tune_xgb_t2,
   xgb_cv_es_t2,
   tune_svm_t2,
   tune_mars_t2,
   best_glmnet_t2,
   best_rf_t2,
   best_xgb_params_t2,
   best_xgb_trees_t2,
   best_svm_t2,
   best_mars_t2,
   glmnet_final_t2, 
   rf_final_t2, 
   xgb_final_t2, 
   xgb_tree_depth_2_t2, 
   svm_final_t2, 
   mars_final_t2, 
   results_t2,
   summary_t2,
   best_wf_id_t2,
   best_fit_t2, 
   xgb_id_t2, 
   xgb_fit_t2, 
   xgb_tree2_id_t2, 
   xgb_tree2_fit_t2)
   

# Analyses poussées : glmnet + XGBoost pour Test 1 et Test 2 ----
# + investigation tree_depth=2 pour Test 2

cat("Test 1 best model:", results_proteomic_ALS_occurrence_tidymodels$test_1$best_wf_id_t1, "\n")
cat("Test 2 best model:", results_proteomic_ALS_occurrence_tidymodels$test_2$best_wf_id_t2, "\n\n")


## Test 1 interpretation : glmnet (winner) + XGBoost (comparison) ----

### glmnet test 1 ----
glmnet_t1_result <- extract_glmnet_coefs(
  results_proteomic_ALS_occurrence_tidymodels$test_1$best_fit_t1)

cat("  Lambda optimal:", round(glmnet_t1_result$lambda_opt, 6), "\n")
cat("  Selected variables:", glmnet_t1_result$n_selected, "\n")
cat("  Selected proteins:", nrow(glmnet_t1_result$coefs_protein), "\n\n")
  
# Table Top 89 predictors (by |coefficient|)
results_proteomic_ALS_occurrence_tidymodels$test_1$t_glmnet_t1 <- 
  glmnet_t1_result$coefs_protein |> 
    slice_head(n = 89) |>
    as.data.frame() |>
    mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  select(-abs_coef, -direction) |>
    flextable() |> 
    colformat_double(digits = 3) |>
    flextable::font(fontname = "Calibri", part = "all") |> 
    flextable::fontsize(size = 10, part = "all") |>
    padding(padding.top = 0, padding.bottom = 0, part = "all") |>
    set_table_properties(align = "left") |>
    autofit()
  
# Plot des 89 predicteurs selectionnés
results_proteomic_ALS_occurrence_tidymodels$test_1$f_glmnet_t1 <- glmnet_t1_result$coefs_protein |>
    slice_head(n = 89) |>
    mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_metabolism_|proteomic_immun_res_"), 
           Feature = fct_reorder(Feature, abs_coef)) |>
    ggplot(aes(x = coefficient, y = Feature, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("↑ ALS risk" = "firebrick", "↓ ALS risk" = "steelblue")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Test 1 — glmnet coefficients (top 89 predictors)",
         x = "Coefficient (standardized)", y = NULL) +
    theme_minimal() + 
  theme(legend.position = "bottom")


### xgboost test 1 ----
xgb_t1_shap <- extract_xgb_and_shap(
  results_proteomic_ALS_occurrence_tidymodels$test_1$xgb_fit_t1,
  data = data_t1,
  test_label  = "Test 1 XGBoost")

# Top 20 variables by SHAP importance
results_proteomic_ALS_occurrence_tidymodels$test_1$t_xgboost_t1 <- xgb_t1_shap$shap_summary |> 
  slice_head(n = 20) |> 
  as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  flextable() |> 
  #colformat_double(digits = 5) |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()
  

# SHAP beeswarm plot
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_beeswarm_t1 <- sv_importance(xgb_t1_shap$shap, kind = "beeswarm", max_display = 20) +
        labs(title = "Test 1 XGBoost — SHAP beeswarm (top 20)")

# SHAP bar plot
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_barplot_t1 <- sv_importance(xgb_t1_shap$shap, kind = "bar", max_display = 20) +
        labs(title = "Test 1 XGBoost — SHAP importance (top 20)")

# Dependence plots (top 6)
top6_t1 <- xgb_t1_shap$shap_summary |> slice_head(n = 6) |> pull(Feature)
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_dependanceplot_t1 <- lapply(top6_t1, function(prot) {
  sv_dependence(xgb_t1_shap$shap, v = prot) +
    labs(title = paste("Test 1 —", prot)) +
    theme_minimal()
})
names(results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_dependanceplot_t1) <- top6_t1
rm(top6_t1)


## Test 2 interpreation : glmnet (winner) + XGBoost (comparison) ----

### glmnet test 2 ----
glmnet_t2_result <- extract_glmnet_coefs(
  results_proteomic_ALS_occurrence_tidymodels$test_2$best_fit_t2)

cat("  Lambda optimal:", round(glmnet_t2_result$lambda_opt, 6), "\n")
cat("  Selected variables:", glmnet_t2_result$n_selected, "\n")
cat("  Selected proteins:", nrow(glmnet_t2_result$coefs_protein), "\n\n")

# Top 110 predictors (by |coefficient|)
results_proteomic_ALS_occurrence_tidymodels$test_2$t_glmnet_t2 <- glmnet_t2_result$coefs_protein |> 
    slice_head(n = 110) |>
    as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  select(-abs_coef, -direction) |>
    flextable() |> 
    colformat_double(digits = 2) |>
    flextable::font(fontname = "Calibri", part = "all") |> 
    flextable::fontsize(size = 10, part = "all") |>
    padding(padding.top = 0, padding.bottom = 0, part = "all") |>
    set_table_properties(align = "left") |>
    autofit()
  
# Plot
results_proteomic_ALS_occurrence_tidymodels$test_2$p_glmnet_t2 <- 
  glmnet_t2_result$coefs_protein |>
    slice_head(n = 110) |>
    mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_metabolism_|proteomic_immun_res_"), 
           Feature = fct_reorder(Feature, abs_coef)) |>
    ggplot(aes(x = coefficient, y = Feature, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("↑ ALS risk" = "firebrick", "↓ ALS risk" = "steelblue")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Test 2 — glmnet coefficients (top 110 predictors)",
         x = "Coefficient (standardized)", y = NULL) +
    theme_minimal() + 
  theme(legend.position = "bottom")


### xgboost test 2 ----
xgb_t2_shap <- extract_xgb_and_shap(
  results_proteomic_ALS_occurrence_tidymodels$test_2$xgb_fit_t2,
  data = data_t2,
  test_label  = "Test 2 XGBoost")

# Top 20 variables by SHAP importance
results_proteomic_ALS_occurrence_tidymodels$test_2$t_xgboost_t2 <- xgb_t2_shap$shap_summary |> 
  slice_head(n = 20) |> 
  as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  flextable() |> 
  colformat_double(digits = 2) |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()

# SHAP beeswarm plot
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_beeswarm_t2 <- sv_importance(xgb_t2_shap$shap, kind = "beeswarm", max_display = 20) +
        labs(title = "Test 2 XGBoost — SHAP beeswarm (top 20)")

# SHAP bar plot
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_barplot_t2 <- sv_importance(xgb_t2_shap$shap, kind = "bar", max_display = 20) +
        labs(title = "Test 2 XGBoost — SHAP importance (top 20)")

# Dependence plots (top 6)
top6_t2 <- xgb_t2_shap$shap_summary |> 
  slice_head(n = 6) |> 
  pull(Feature)

results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_dependanceplot_t2 <- lapply(top6_t2, function(prot) {
  sv_dependence(xgb_t2_shap$shap, v = prot) +
    labs(title = paste("Test 2 —", prot)) +
    theme_minimal()
})
names(results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_dependanceplot_t2) <- top6_t2

# Dependence plots by follow-up duration (top 6)
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_dependanceplot_by_followup_t2 <- lapply(top6_t2, function(prot) {
  sv_dependence(xgb_t2_shap$shap, 
                v = prot, 
                color_var = "follow_up_no_na_y") +
    scale_color_viridis_c(name = "Follow-up\n(years)") +
    labs(title = paste("Test 2 —", prot, "× follow_up")) +
    theme_minimal()
})
names(results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_dependanceplot_by_followup_t2) <- top6_t2
rm(top6_t2)

# Interaction follow_up × protéines
fu_breaks <- quantile(data_t2$follow_up_no_na_y, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
fu_group_t2 <- cut(data_t2$follow_up_no_na_y,
                   breaks = fu_breaks,
                   labels = c("Early (≤1/3)", "Intermediate", "Late (≥2/3)"),
                   include.lowest = TRUE)

top15_prot_t2 <- xgb_t2_shap$shap_summary |> slice_head(n = 15) |> pull(Feature)

shap_by_fu <- as_tibble(xgb_t2_shap$shap_matrix[, top15_prot_t2]) |>
  mutate(fu_group = fu_group_t2) |>
  pivot_longer(-fu_group, names_to = "Feature", values_to = "shap") |>
  group_by(fu_group, Feature) |>
  summarise(
    mean_abs_shap = mean(abs(shap)),
    mean_shap     = mean(shap),
    .groups = "drop") |>
  mutate(direction = if_else(mean_shap > 0, "↑ risk", "↓ risk"))

cat("\nSHAP |value| by follow-up tertile:\n")
results_proteomic_ALS_occurrence_tidymodels$test_2$f_fu_shap <- ggplot(shap_by_fu, aes(x = fu_group, y = mean_abs_shap, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ reorder(Feature, -mean_abs_shap), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("↑ risk" = "firebrick", "↓ risk" = "steelblue")) +
  labs(title = "Test 2 — SHAP |value| by follow-up tertile (top 15 proteins)",
       x = "Follow-up tertile", y = "Mean |SHAP|", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Trajectory analysis
cat("\nProtein signal trajectory:\n")
results_proteomic_ALS_occurrence_tidymodels$test_2$shap_fu_trend <- shap_by_fu |>
  select(fu_group, Feature, mean_abs_shap) |>
  pivot_wider(names_from = fu_group, values_from = mean_abs_shap) |>
  rename(early = `Early (≤1/3)`, intermediate = Intermediate, late = `Late (≥2/3)`) |>
  mutate(
    trend = case_when(
      late > early * 1.3  ~ "Late marker (↑ closer to diagnosis)",
      early > late * 1.3  ~ "Early marker (↑ years before)",
      TRUE                ~ "Stable across time"),
    ratio_late_early = round(late / (early + 1e-6), 2)) |>
  arrange(desc(ratio_late_early)) |> 
  as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  flextable() |> 
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()


rm(fu_breaks, fu_group_t2, top15_prot_t2, shap_by_fu)

## Comparison: glmnet vs XGBoost across both tests ----
cat("COMPARISON: glmnet vs XGBoost COEFFICIENTS/SHAP\n")

# Protéines dans les top 20 glmnet (T1 et T2)
if (glmnet_t1_result$success) {
  glmnet_t1_top <- glmnet_t1_result$coefs_protein |> slice_head(n = 20) |> pull(Feature)
} else {
  glmnet_t1_top <- c()
}

if (glmnet_t2_result$success) {
  glmnet_t2_top <- glmnet_t2_result$coefs_protein |> slice_head(n = 20) |> pull(Feature)
} else {
  glmnet_t2_top <- c()
}

# Protéines dans les top 20 XGBoost SHAP (T1 et T2)
xgb_t1_top <- xgb_t1_shap$shap_summary  |> slice_head(n = 20) |> pull(Feature)
xgb_t2_top <- xgb_t2_shap$shap_summary |> slice_head(n = 20) |> pull(Feature)

# Venn-like analysis : overlap
cat("Test 1 — glmnet vs XGBoost overlap:\n")
results_proteomic_ALS_occurrence_tidymodels$test_1$overlap_t1 <- intersect(glmnet_t1_top, xgb_t1_top)
cat("  Common top 20 proteins:", length(results_proteomic_ALS_occurrence_tidymodels$test_1$overlap_t1), "\n")
if (length(results_proteomic_ALS_occurrence_tidymodels$test_1$overlap_t1) > 0) {
  cat("  Proteins:", paste(results_proteomic_ALS_occurrence_tidymodels$test_1$overlap_t1, collapse = ", "), "\n")
}

cat("\nTest 2 — glmnet vs XGBoost overlap:\n")
results_proteomic_ALS_occurrence_tidymodels$test_2$overlap_t2 <- intersect(glmnet_t2_top, xgb_t2_top)
cat("  Common top 20 proteins:", length(results_proteomic_ALS_occurrence_tidymodels$test_2$overlap_t2), "\n")
if (length(results_proteomic_ALS_occurrence_tidymodels$test_2$overlap_t2) > 0) {
  cat("  Proteins:", paste(results_proteomic_ALS_occurrence_tidymodels$test_2$overlap_t2, collapse = ", "), "\n")
}

rm(glmnet_t1_result, glmnet_t1_top, xgb_t1_shap, xgb_t1_top,  
   glmnet_t2_result, glmnet_t2_top, xgb_t2_shap, xgb_t2_top)



## Additional investigation: xgboost tree_depth = 2 for Test 2 ----
xgb_tree2_t2_shap <- extract_xgb_and_shap(
  results_proteomic_ALS_occurrence_tidymodels$test_2$xgb_tree2_fit_t2,
  data = data_t2,
  test_label  = "Test 2 XGBoost - tree depth = 2")

# Top 20 variables by SHAP importance
results_proteomic_ALS_occurrence_tidymodels$test_2$t_xgb_tree2_t2 <- xgb_tree2_t2_shap$shap_summary |> 
  slice_head(n = 20) |> 
  as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  flextable() |> 
  colformat_double(digits = 2) |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()

# SHAP beeswarm plot
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_beeswarm_t2 <- 
  sv_importance(xgb_tree2_t2_shap$shap, kind = "beeswarm", max_display = 20) +
  labs(title = "Test 2 XGBoost - tree depth = 2 — SHAP beeswarm (top 20)")

# SHAP bar plot
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_barplot_t2 <- 
  sv_importance(xgb_tree2_t2_shap$shap, kind = "bar", max_display = 20) +
  labs(title = "Test 2 XGBoost - tree depth = 2 — SHAP importance (top 20)")

# Dependence plots (top 6)
top6_t2_tree2 <- xgb_tree2_t2_shap$shap_summary |> 
  slice_head(n = 6) |> 
  pull(Feature)

results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_dependanceplot_t2 <- 
  lapply(top6_t2_tree2, function(prot) {
  sv_dependence(xgb_tree2_t2_shap$shap, v = prot) +
    labs(title = paste("Test 2 XGBoost - tree depth = 2 —", prot)) +
    theme_minimal()
})
names(results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_dependanceplot_t2) <- top6_t2_tree2

# Dependence plots by follow-up duration (top 6)
results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_dependanceplot_by_followup_t2 <- 
  lapply(top6_t2_tree2, function(prot) {
  sv_dependence(xgb_tree2_t2_shap$shap, 
                v = prot, 
                color_var = "follow_up_no_na_y") +
    scale_color_viridis_c(name = "Follow-up\n(years)") +
    labs(title = paste("Test 2 - tree depth = 2 - ", prot, "× follow_up")) +
    theme_minimal()
})
names(results_proteomic_ALS_occurrence_tidymodels$test_2$f_xgb_tree2_dependanceplot_by_followup_t2) <- top6_t2_tree2
rm(top6_t2_tree2)

# Interaction follow_up × protéines
fu_breaks_tree2 <- quantile(data_t2$follow_up_no_na_y, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
fu_group_t2_tree2 <- cut(data_t2$follow_up_no_na_y,
                         breaks = fu_breaks_tree2,
                         labels = c("Early (≤1/3)", "Intermediate", "Late (≥2/3)"),
                         include.lowest = TRUE)

top15_prot_t2_tree2 <- xgb_tree2_t2_shap$shap_summary |> slice_head(n = 15) |> pull(Feature)

shap_by_fu_tree2 <- as_tibble(xgb_tree2_t2_shap$shap_matrix[, top15_prot_t2_tree2]) |>
  mutate(fu_group = fu_group_t2_tree2) |>
  pivot_longer(-fu_group, names_to = "Feature", values_to = "shap") |>
  group_by(fu_group, Feature) |>
  summarise(
    mean_abs_shap = mean(abs(shap)),
    mean_shap     = mean(shap),
    .groups = "drop") |>
  mutate(direction = if_else(mean_shap > 0, "↑ risk", "↓ risk"))

cat("\nSHAP |value| by follow-up tertile:\n")
results_proteomic_ALS_occurrence_tidymodels$test_2$f_fu_shap_tree2 <- 
  ggplot(shap_by_fu_tree2, aes(x = fu_group, y = mean_abs_shap, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ reorder(Feature, -mean_abs_shap), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("↑ risk" = "firebrick", "↓ risk" = "steelblue")) +
  labs(title = "Test 2 — tree depth = 2 - SHAP |value| by follow-up tertile (top 15 proteins)",
       x = "Follow-up tertile", y = "Mean |SHAP|", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Trajectory analysis
cat("\nProtein signal trajectory:\n")
results_proteomic_ALS_occurrence_tidymodels$test_2$shap_fu_trend_tree2 <- shap_by_fu_tree2 |>
  select(fu_group, Feature, mean_abs_shap) |>
  pivot_wider(names_from = fu_group, values_from = mean_abs_shap) |>
  rename(early = `Early (≤1/3)`, intermediate = Intermediate, late = `Late (≥2/3)`) |>
  mutate(
    trend = case_when(
      late > early * 1.3  ~ "Late marker (↑ closer to diagnosis)",
      early > late * 1.3  ~ "Early marker (↑ years before)",
      TRUE                ~ "Stable across time"),
    ratio_late_early = round(late / (early + 1e-6), 2)) |>
  arrange(desc(ratio_late_early)) |> 
  as.data.frame() |>
  mutate(Feature = str_remove(Feature, "proteomic_neuro_explo_|proteomic_immun_res_|proteomic_metabolism_")) |>
  flextable() |> 
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()


rm(fu_breaks_tree2, fu_group_t2_tree2, top15_prot_t2_tree2, shap_by_fu_tree2, xgb_tree2_t2_shap)

rm(step_select_auc_top, bake.step_select_auc_top, prep.step_select_auc_top, 
   print.step_select_auc_top, make_folds_checked, summarise_wf_results, fit_best_workflow, 
   make_recipe_full, make_recipe_top, 
   data_t1, folds_t1,
   rec_full_t1, rec_top20_t1, ncol_preds_t1, 
   data_t2, folds_t2,
   rec_full_t2, rec_top20_t2, ncol_preds_t2)

saveRDS(
  results_proteomic_ALS_occurrence_tidymodels,
  file = "~/Documents/POP_ALS_2025_02_03/2_output/2.6.3_results_proteomic_ALS_occurrence_tidymodels.rds")



# Aline Davias
# Juin 2026
# Prédiction ALS par protéomique – tidymodels 


# Packages and data ----
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
        coefficient_raw = as.numeric(coef_matrix[, 1]), # On stocke la valeur brute
        stringsAsFactors = FALSE),
      rownames = NA) |>
      filter(Feature != "(Intercept)" & coefficient_raw != 0)
    
    coefs_annotated <- coefs |>
      mutate(
        # CLÉ DU ROBUST-REPAIR : Comme glmnet modélise l'alternative à la référence,
        # on inverse le signe pour obtenir l'effet direct sur le groupe "case" (SLA).
        coefficient = -coefficient_raw, 
        direction  = case_when(
          coefficient > 0 ~ "↑ ALS risk",  # Un coeff positif augmente désormais le risque
          coefficient < 0 ~ "↓ ALS risk",  # Un coeff négatif diminue le risque
          TRUE ~ "neutral"
        ),
        abs_coef = abs(coefficient)) |>
      arrange(desc(abs_coef)) # Tri par importance absolue globale
    
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
  select(als,  match, 
         birth_year, sex, bmi, smoking_2cat_i, all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(1, 0), labels = c("case", "control")))

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


### A3. XGBoost : Tuning complet via tidymodels (sans étape xgb.cv externe) ----

xgb_spec_tune <- boost_tree(
  trees          = tune(),          # Mis en tune() pour remplacer l'early stopping
  tree_depth     = tune(),
  learn_rate     = tune(),
  mtry           = tune(),
  min_n          = tune(),          # Traduction native de min_child_weight
  loss_reduction = tune(),
  sample_size    = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

# Définition de la grille de recherche (on y intègre les arbres)
xgb_grid_t1 <- grid_space_filling(
  trees(range          = c(100, 2000)), # Espace de recherche sur le volume d'arbres
  tree_depth(range     = c(1, 4)),
  learn_rate(range     = c(-3, -1)),    # 0.001 à 0.1 en log10
  mtry(range           = c(floor(ncol_preds_t1 * 0.1), 
                           floor(ncol_preds_t1 * 0.6))),
  min_n(range          = c(1, 20)),
  loss_reduction(range = c(-5, 0)),     # gamma en log10
  size = 50)

wf_xgb_t1 <- workflow() |> 
  add_recipe(rec_full_t1) |> 
  add_model(xgb_spec_tune)

set.seed(1996)
tune_xgb_t1 <- tune_grid(
  wf_xgb_t1,
  resamples = folds_t1,
  grid      = xgb_grid_t1,
  metrics   = metric,
  control   = control_grid(save_pred = TRUE, verbose = TRUE))

# Extraction directe des meilleurs paramètres globaux
best_xgb_params_t1 <- select_best(tune_xgb_t1, metric = "roc_auc")



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

# xgb_final_t1 <- boost_tree(
#   trees        = best_xgb_trees_t1,
#   tree_depth   = best_xgb_params_t1$tree_depth,
#   learn_rate   = best_xgb_params_t1$learn_rate,
#   mtry         = best_xgb_params_t1$mtry,
#   min_n        = best_xgb_params_t1$min_n,
#   loss_reduction = best_xgb_params_t1$loss_reduction,
#   sample_size  = 0.8) |>
#   set_engine("xgboost", 
#              nthread = parallel::detectCores() - 1) |>
#   set_mode("classification")

xgb_final_t1 <- boost_tree(
  trees          = best_xgb_params_t1$trees,
  tree_depth     = best_xgb_params_t1$tree_depth,
  learn_rate     = best_xgb_params_t1$learn_rate,
  mtry           = best_xgb_params_t1$mtry,
  min_n          = best_xgb_params_t1$min_n,
  loss_reduction = best_xgb_params_t1$loss_reduction,
  sample_size    = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
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

xgb_id_t1 <- summary_t1$wflow_id[4]        # full_xgboost
xgb_fit_t1 <- fit_best_workflow(results_t1, xgb_id_t1, data_t1)

rm(wf_set_t1)


# TEST 1 - sensitivity analysis (remove cases and controls with follow-up<5 years) ----
# TEST 1 sensi : Proteins and covariates 
data_t1_sensi <- bdd_danish |>
  filter(follow_up_no_na_y >= 5) |>
  select(als, match,
         birth_year, sex, bmi, smoking_2cat_i,
         all_of(proteomic)) |>
  mutate(als = factor(als, levels = c(1, 0), labels = c("case", "control"))) 

cat("Dimensions :", nrow(data_t1_sensi), "×", ncol(data_t1_sensi) - 2, "prédicteurs\n")
cat("Outcome :", table(data_t1_sensi$als), "\n")

folds_t1_sensi     <- make_folds_checked(data_t1_sensi)
ncol_preds_t1_sensi <- ncol(data_t1_sensi) - 2  # -als -match 

rec_full_t1_sensi  <- make_recipe_full(data_t1_sensi)
rec_top20_t1_sensi <- make_recipe_top(data_t1_sensi, n_top = 20)


## Preparation of the recipe with interactions 
rec_interact_t1_sensi <- make_recipe_full(data_t1_sensi) 

## TEST 1 sensi - PHASE A : Tuning des hyperparamètres ----
glmnet_spec_tune_t1_sensi <- 
  logistic_reg(penalty = tune(), 
               mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t1_sensi <- workflow() |> 
  add_recipe(rec_interact_t1_sensi) |> 
  add_model(glmnet_spec_tune_t1_sensi)

# Lancement de l'optimisation bayésienne
set.seed(1996)
tune_glmnet_t1_sensi <- tune_bayes(
  wf_glmnet_t1_sensi, 
  resamples = folds_t1_sensi,                                                         # on garde toujours les memes folds
  param_info = parameters(dials::penalty(range = c(-5, 0)), 
                          mixture(range = c(0, 1))),
  iter = 40, 
  initial = 10, 
  metrics = metric,
  control = control_bayes(save_pred = TRUE, verbose = TRUE))

best_glmnet_t1_sensi <- select_best(tune_glmnet_t1_sensi, metric = "roc_auc")               # Extraction des meilleurs hyperparamètres 

rm(glmnet_spec_tune_t1_sensi, wf_glmnet_t1_sensi)

## TEST 1 sensi - PHASE B : Entraînement du modèle final sur toutes les données ----
# ici, par rapport à test 1 et 2, on n'utilise pas workflow_map() et fit_best_workflow mais directement fit() car on a qu'un seul model

glmnet_final_t1_sensi <- logistic_reg(
  penalty = best_glmnet_t1_sensi$penalty, 
  mixture = best_glmnet_t1_sensi$mixture) |>
  set_engine("glmnet") |> 
  set_mode("classification")

wf_final_set_t1_sensi <- workflow() |> 
  add_recipe(rec_interact_t1_sensi) |> 
  add_model(glmnet_final_t1_sensi)

set.seed(1996)
best_fit_t1_sensi <- fit(wf_final_set_t1_sensi, data = data_t1_sensi)     

glmnet_t1_result_sensi <- extract_glmnet_coefs(best_fit_t1_sensi)




# TEST 2 : Protéines + covariables + follow_up_no_na_y ----

data_t2 <- bdd_danish |>
  select(als, match,
         follow_up_no_na_y, birth_year, sex, bmi, smoking_2cat_i,
          all_of(proteomic)) |>
 mutate(als = factor(als, levels = c(1, 0), labels = c("case", "control")))

cat("Dimensions :", nrow(data_t2), "×", ncol(data_t2) - 2, "prédicteurs\n")
cat("Outcome :", table(data_t2$als), "\n")

folds_t2     <- make_folds_checked(data_t2)
ncol_preds_t2 <- ncol(data_t2) - 2  # -als -match 

rec_full_t2  <- make_recipe_full(data_t2)
rec_top20_t2 <- make_recipe_top(data_t2, n_top = 20)
#rec_top10_t2 <- make_recipe_top(data_t2, n_top = 10)

## TEST 2 – PHASE A : Tuning hyperparamètres ----

### A1. glmnet ----
glmnet_spec_tune <- 
  logistic_reg(penalty = tune(), 
               mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t2 <- workflow() |> 
  add_recipe(rec_full_t2) |> 
  add_model(glmnet_spec_tune)

set.seed(1996)
tune_glmnet_t2 <- tune_bayes(
  wf_glmnet_t2, 
  resamples = folds_t2,
  param_info = parameters(dials::penalty(range = c(-5, 0)), 
                          mixture(range = c(0, 1))),
  iter = 40, 
  initial = 10, 
  metrics = metric,
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


### A3. XGBoost : Tuning complet via tidymodels (sans étape xgb.cv externe) ----

xgb_spec_tune <- boost_tree(
  trees          = tune(),          # Mis en tune() pour remplacer l'early stopping
  tree_depth     = tune(),
  learn_rate     = tune(),
  mtry           = tune(),
  min_n          = tune(),          # Traduction native de min_child_weight
  loss_reduction = tune(),
  sample_size    = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
  set_mode("classification")

# Définition de la grille de recherche (on y intègre les arbres)
xgb_grid_t2 <- grid_space_filling(
  trees(range          = c(100, 2000)), # Espace de recherche sur le volume d'arbres
  tree_depth(range     = c(1, 4)),
  learn_rate(range     = c(-3, -1)),    # 0.001 à 0.1 en log10
  mtry(range           = c(floor(ncol_preds_t1 * 0.1), 
                           floor(ncol_preds_t1 * 0.6))),
  min_n(range          = c(1, 20)),
  loss_reduction(range = c(-5, 0)),     # gamma en log10
  size = 50)

wf_xgb_t2 <- workflow() |> 
  add_recipe(rec_full_t2) |> 
  add_model(xgb_spec_tune)

set.seed(1996)
tune_xgb_t2 <- tune_grid(
  wf_xgb_t2,
  resamples = folds_t2,
  grid      = xgb_grid_t2,
  metrics   = metric,
  control   = control_grid(save_pred = TRUE, verbose = TRUE))

# Extraction directe des meilleurs paramètres globaux
best_xgb_params_t2 <- select_best(tune_xgb_t2, metric = "roc_auc")



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

# xgb_final_t2 <- boost_tree(
#   trees = best_xgb_trees_t2, 
#   tree_depth = best_xgb_params_t2$tree_depth,
#   learn_rate = best_xgb_params_t2$learn_rate, 
#   mtry = best_xgb_params_t2$mtry,
#   min_n = best_xgb_params_t2$min_n, 
#   loss_reduction = best_xgb_params_t2$loss_reduction, 
#   sample_size = 0.8) |>
#   set_engine("xgboost", 
#              nthread = parallel::detectCores() - 1) |>
#   set_mode("classification")

xgb_final_t2 <- boost_tree(
  trees          = best_xgb_params_t2$trees,
  tree_depth     = best_xgb_params_t2$tree_depth,
  learn_rate     = best_xgb_params_t2$learn_rate,
  mtry           = best_xgb_params_t2$mtry,
  min_n          = best_xgb_params_t2$min_n,
  loss_reduction = best_xgb_params_t2$loss_reduction,
  sample_size    = 0.8) |>
  set_engine("xgboost", nthread = parallel::detectCores() - 1) |>
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
                             xgboost = xgb_final_t2)),
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

rm(wf_set_t2)


# TEST 3 : GLMNET only: Proteins, covariates with two-way interactions with follow_up_no_na_y ----

## Preparation of the recipe with interactions 
rec_interact_t3 <- make_recipe_full(data_t2) |> 
  step_interact(terms = ~ all_of(proteomic):follow_up_no_na_y)

## TEST 3 - PHASE A : Tuning des hyperparamètres ----
glmnet_spec_tune_t3 <- 
  logistic_reg(penalty = tune(), 
               mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t3 <- workflow() |> 
  add_recipe(rec_interact_t3) |> 
  add_model(glmnet_spec_tune_t3)

# Lancement de l'optimisation bayésienne sur les mêmes folds que test 1 et test 2 pour meilleure comparaison
set.seed(1996)
tune_glmnet_t3 <- tune_bayes(
  wf_glmnet_t3, 
  resamples = folds_t2,                                                         # on garde toujours les memes folds
  param_info = parameters(dials::penalty(range = c(-5, 0)), 
                          mixture(range = c(0, 1))),
  iter = 40, 
  initial = 10, 
  metrics = metric,
  control = control_bayes(save_pred = TRUE, verbose = TRUE))

best_glmnet_t3 <- select_best(tune_glmnet_t3, metric = "roc_auc")               # Extraction des meilleurs hyperparamètres 

rm(glmnet_spec_tune_t3, wf_glmnet_t3)

## TEST 3 - PHASE B : Entraînement du modèle final sur toutes les données ----
# ici, par rapport à test 1 et 2, on n'utilise pas workflow_map() et fit_best_workflow mais directement fit() car on a qu'un seul model

glmnet_final_t3 <- logistic_reg(
  penalty = best_glmnet_t3$penalty, 
  mixture = best_glmnet_t3$mixture) |>
  set_engine("glmnet") |> 
  set_mode("classification")

wf_final_set_t3 <- workflow() |> 
  add_recipe(rec_interact_t3) |> 
  add_model(glmnet_final_t3)

set.seed(1996)
best_fit_t3 <- fit(wf_final_set_t3, data = data_t2)     

glmnet_t3_result <- extract_glmnet_coefs(best_fit_t3)


# TEST 3 - sensitivity analysis (remove cases and controls with follow-up<5 years) ----
# TEST 3 sensi : Proteins, covariates with two-way interactions with follow_up_no_na_y

data_t3_sensi <- bdd_danish |>
  select(als, match,
         follow_up_no_na_y, birth_year, sex, bmi, smoking_2cat_i,
         all_of(proteomic)) |>
  filter(follow_up_no_na_y >= 5) |>
  mutate(als = factor(als, levels = c(1, 0), labels = c("case", "control")))

cat("Dimensions :", nrow(data_t3_sensi), "×", ncol(data_t3_sensi) - 2, "prédicteurs\n")
cat("Outcome :", table(data_t3_sensi$als), "\n")

folds_t3_sensi     <- make_folds_checked(data_t3_sensi)
ncol_preds_t3_sensi <- ncol(data_t3_sensi) - 2  # -als -match 

rec_full_t3_sensi  <- make_recipe_full(data_t3_sensi)
rec_top20_t3_sensi <- make_recipe_top(data_t3_sensi, n_top = 20)


## Preparation of the recipe with interactions 
rec_interact_t3_sensi <- make_recipe_full(data_t3_sensi) |> 
  step_interact(terms = ~ all_of(proteomic):follow_up_no_na_y)

## TEST 3 sensi - PHASE A : Tuning des hyperparamètres ----
glmnet_spec_tune_t3_sensi <- 
  logistic_reg(penalty = tune(), 
               mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")

wf_glmnet_t3_sensi <- workflow() |> 
  add_recipe(rec_interact_t3_sensi) |> 
  add_model(glmnet_spec_tune_t3_sensi)

# Lancement de l'optimisation bayésienne sur les mêmes folds que test 1 et test 2 pour meilleure comparaison
set.seed(1996)
tune_glmnet_t3_sensi <- tune_bayes(
  wf_glmnet_t3_sensi, 
  resamples = folds_t3_sensi,                                                         # on garde toujours les memes folds
  param_info = parameters(dials::penalty(range = c(-5, 0)), 
                          mixture(range = c(0, 1))),
  iter = 40, 
  initial = 10, 
  metrics = metric,
  control = control_bayes(save_pred = TRUE, verbose = TRUE))

best_glmnet_t3_sensi <- select_best(tune_glmnet_t3_sensi, metric = "roc_auc")               # Extraction des meilleurs hyperparamètres 

rm(glmnet_spec_tune_t3_sensi, wf_glmnet_t3_sensi)

## TEST 3 sensi - PHASE B : Entraînement du modèle final sur toutes les données ----
# ici, par rapport à test 1 et 2, on n'utilise pas workflow_map() et fit_best_workflow mais directement fit() car on a qu'un seul model

glmnet_final_t3_sensi <- logistic_reg(
  penalty = best_glmnet_t3_sensi$penalty, 
  mixture = best_glmnet_t3_sensi$mixture) |>
  set_engine("glmnet") |> 
  set_mode("classification")

wf_final_set_t3_sensi <- workflow() |> 
  add_recipe(rec_interact_t3_sensi) |> 
  add_model(glmnet_final_t3_sensi)

set.seed(1996)
best_fit_t3_sensi <- fit(wf_final_set_t3_sensi, data = data_t3_sensi)     

glmnet_t3_result_sensi <- extract_glmnet_coefs(best_fit_t3_sensi)


# Sauvegarde ----
results_proteomic_ALS_occurrence_tidymodels <- list(
  
  test_1 = list(
    #data_t1 = data_t1,
    #folds_t1 = folds_t1,
    # Tuning test 1
    tune_glmnet_t1 = tune_glmnet_t1,
    tune_rf_t1 = tune_rf_t1,
    tune_xgb_t1 = tune_xgb_t1,
    #xgb_cv_es_t1 = xgb_cv_es_t1,
    tune_svm_t1 = tune_svm_t1,
    tune_mars_t1 = tune_mars_t1,
    # Hyperparamètres optimaux test 1
    best_glmnet_t1 = best_glmnet_t1,
    best_rf_t1 = best_rf_t1,
    best_xgb_params_t1 = best_xgb_params_t1,
    #best_xgb_trees_t1 = best_xgb_trees_t1,
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
  
  test_1_sensi =
    list( 
      #data_t1_sensi = data_t1_sensi,        
      #folds_t1_sensi = folds_t1_sensi,     
      # Tuning test 1 sensi
      tune_glmnet_t1_sensi = tune_glmnet_t1_sensi,
      # Hyperparamètres optimaux
      best_glmnet_t1_sensi  = best_glmnet_t1_sensi,
      # comparaison finale test 1 sensi
      glmnet_final_t1_sensi = glmnet_final_t1_sensi, 
      wf_final_set_t1_sensi = wf_final_set_t1_sensi, 
      best_fit_t1_sensi  = best_fit_t1_sensi, 
      glmnet_t1_result_sensi = glmnet_t1_result_sensi), 
  
  test_2 = list(
    #data_t2 = data_t2,
    #folds_t2 = folds_t2,
    # Tuning test 2 
    tune_glmnet_t2 = tune_glmnet_t2,
    tune_rf_t2 = tune_rf_t2,
    tune_xgb_t2 = tune_xgb_t2,
    #xgb_cv_es_t2 = xgb_cv_es_t2,
    tune_svm_t2  = tune_svm_t2,
    tune_mars_t2  = tune_mars_t2,
    # Hyperparamètres optimaux
    best_glmnet_t2  = best_glmnet_t2,
    best_rf_t2  = best_rf_t2,
    best_xgb_params_t2  = best_xgb_params_t2,
    #best_xgb_trees_t2  = best_xgb_trees_t2,
    best_svm_t2  = best_svm_t2,
    best_mars_t2  = best_mars_t2,
    # comparaison finale test 2
    results_t2 = results_t2,
    summary_t2  = summary_t2,
    best_wf_id_t2 = best_wf_id_t2,
    best_fit_t2  = best_fit_t2, 
    xgb_id_t2 = xgb_id_t2, 
    xgb_fit_t2 = xgb_fit_t2,  
    glmnet_final_t2 = glmnet_final_t2, 
    rf_final_t2 = rf_final_t2, 
    xgb_final_t2 = xgb_final_t2, 
    svm_final_t2 = svm_final_t2, 
    mars_final_t2 = mars_final_t2), 
  
  test_3 = 
  list( 
    #data_t2 = data_t2,         # we used the same as t2
    #folds_t2 = folds_t2,       # we used the same as t2
    # Tuning test 3 
    tune_glmnet_t3 = tune_glmnet_t3,
    # Hyperparamètres optimaux
    best_glmnet_t3  = best_glmnet_t3,
    # comparaison finale test 3
    glmnet_final_t3 = glmnet_final_t3, 
    wf_final_set_t3 = wf_final_set_t3, 
    best_fit_t3  = best_fit_t3, 
    glmnet_t3_result = glmnet_t3_result), 
  
  test_3_sensi =
  list( 
    #data_t3_sensi = data_t3_sensi,        
    #folds_t3_sensi = folds_t3_sensi,     
    # Tuning test 3 sensi
    tune_glmnet_t3_sensi = tune_glmnet_t3_sensi,
    # Hyperparamètres optimaux
    best_glmnet_t3_sensi  = best_glmnet_t3_sensi,
    # comparaison finale test 3 sensi
    glmnet_final_t3_sensi = glmnet_final_t3_sensi, 
    wf_final_set_t3_sensi = wf_final_set_t3_sensi, 
    best_fit_t3_sensi  = best_fit_t3_sensi, 
    glmnet_t3_result_sensi = glmnet_t3_result_sensi))

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
   xgb_tree2_fit_t2,
   
   tune_glmnet_t3, 
   best_glmnet_t3, 
   glmnet_final_t3, 
   wf_final_set_t3, 
   best_fit_t3, 
   selected_predictors, 
   glmnet_t3_result, 
   rec_interact_t3, 
   
   data_t3_sensi, 
   folds_t3_sensi,   
   rec_full_t3_sensi, 
   rec_top20_t3_sensi, 
   rec_interact_t3_sensi, 
   tune_glmnet_t3_sensi,
   best_glmnet_t3_sensi,
   glmnet_final_t3_sensi,
   wf_final_set_t3_sensi, 
   best_fit_t3_sensi, 
   glmnet_t3_result_sensi)



# Analyses poussées : glmnet + XGBoost pour Test 1 et Test 2 ----
# + investigation tree_depth=2 pour Test 2

cat("Test 1 best model:", results_proteomic_ALS_occurrence_tidymodels$test_1$best_wf_id_t1, "\n")
cat("Test 2 best model:", results_proteomic_ALS_occurrence_tidymodels$test_2$best_wf_id_t2, "\n\n")


## Test 1 interpretation : glmnet (winner) + XGBoost (comparison) ----

### glmnet test 1 ----
results_proteomic_ALS_occurrence_tidymodels$test_1$glmnet_t1_result <- extract_glmnet_coefs(
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

# AUC moyenne glmnet 
results_proteomic_ALS_occurrence_tidymodels$test_1$best_auc_value <- show_best(                                                    # Extraction de la valeur numérique moyenne du meilleur AUC
  results_proteomic_ALS_occurrence_tidymodels$test_1$tune_glmnet_t1, 
  metric = "roc_auc", 
  n = 1)$mean

best_predictions_t1 <-                                                          # Récupération des prédictions hors-échantillon
  results_proteomic_ALS_occurrence_tidymodels$test_1$tune_glmnet_t1 |> 
  collect_predictions(
    parameters = results_proteomic_ALS_occurrence_tidymodels$test_1$best_glmnet_t1)

roc_curve_data_t1 <-                                                            # Calcul des points de la courbe ROC
  best_predictions_t1 |>
  roc_curve(truth = als, .pred_case)

results_proteomic_ALS_occurrence_tidymodels$test_1$roc_curve_data_t1 <-         # plot
  roc_curve_data_t1 |>
  autoplot() +
  labs(
    title = "AUC curve - Test 1 (Elastic Net)",
    subtitle = paste0("CV 10-folds — Mean AUC: ", 
                      round(results_proteomic_ALS_occurrence_tidymodels$test_1$best_auc_value, 3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10))

rm(best_predictions_t1, roc_curve_data_t1)


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
  colformat_double(digits = 5) |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  flextable::fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  set_table_properties(align = "left") |>
  autofit()

# SHAP beeswarm plot
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_beeswarm_t1 <- 
  sv_importance(xgb_t1_shap$shap, 
                kind = "beeswarm", 
                max_display = 20) +
  labs(title = "Test 1 XGBoost — SHAP beeswarm (top 20)")

# SHAP bar plot
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_barplot_t1 <- 
  sv_importance(xgb_t1_shap$shap, 
                kind = "bar", 
                max_display = 20) +
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
results_proteomic_ALS_occurrence_tidymodels$test_1$f_xgb_dependanceplot_t1$proteomic_neuro_explo_NEFL


## Test 2 interpreation : glmnet (winner) + XGBoost (comparison) ----

### glmnet test 2 ----
results_proteomic_ALS_occurrence_tidymodels$test_2$glmnet_t2_result <- extract_glmnet_coefs(
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

# AUC moyenne glmnet 
results_proteomic_ALS_occurrence_tidymodels$test_2$best_auc_value <- show_best(                                                    # Extraction de la valeur numérique moyenne du meilleur AUC
  results_proteomic_ALS_occurrence_tidymodels$test_2$tune_glmnet_t2, 
  metric = "roc_auc", 
  n = 1)$mean

best_predictions_t2 <-                                                          # Récupération des prédictions hors-échantillon
  results_proteomic_ALS_occurrence_tidymodels$test_2$tune_glmnet_t2 |> 
  collect_predictions(
    parameters = results_proteomic_ALS_occurrence_tidymodels$test_2$best_glmnet_t2)

roc_curve_data_t2 <-                                                            # Calcul des points de la courbe ROC
  best_predictions_t2 |>
  roc_curve(truth = als, .pred_case)

results_proteomic_ALS_occurrence_tidymodels$test_2$roc_curve_data_t2 <-         # plot
  roc_curve_data_t2 |>
  autoplot() +
  labs(
    title = "AUC curve - Test 2 (Elastic Net including follow-up as predictor)",
    subtitle = paste0("CV 10-folds — Mean AUC: ", 
                      round(results_proteomic_ALS_occurrence_tidymodels$test_2$best_auc_value, 3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10))

rm(best_predictions_t2, roc_curve_data_t2)


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
  colformat_double(digits = 5) |>
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
   glmnet_t2_result, glmnet_t2_top, xgb_t2_shap, xgb_t2_top, 
   fu_breaks_tree2, fu_group_t2_tree2, top15_prot_t2_tree2, shap_by_fu_tree2, xgb_tree2_t2_shap)

## Test 3 interpretation ----
results_proteomic_ALS_occurrence_tidymodels$test_3$best_auc_value <- show_best(                                                    # Extraction de la valeur numérique moyenne du meilleur AUC
  results_proteomic_ALS_occurrence_tidymodels$test_3$tune_glmnet_t3, 
  metric = "roc_auc", 
  n = 1)$mean

best_predictions_t3 <-                                                          # Récupération des prédictions hors-échantillon
  results_proteomic_ALS_occurrence_tidymodels$test_3$tune_glmnet_t3 |> 
  collect_predictions(
    parameters = results_proteomic_ALS_occurrence_tidymodels$test_3$best_glmnet_t3)

roc_curve_data_t3 <-                                                            # Calcul des points de la courbe ROC
  best_predictions_t3 |>
  roc_curve(truth = als, .pred_case)

results_proteomic_ALS_occurrence_tidymodels$test_3$roc_curve_data_t3 <-         # plot
  roc_curve_data_t3 |>
  autoplot() +
  labs(
    title = "AUC curve - Test 3 (Elastic Net with follow-up interactions)",
    subtitle = paste0("CV 10-folds — Mean AUC: ", 
                      round(results_proteomic_ALS_occurrence_tidymodels$test_3$best_auc_value, 3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10))

rm(best_predictions_t3, roc_curve_data_t3)


rm(step_select_auc_top, bake.step_select_auc_top, prep.step_select_auc_top, 
   print.step_select_auc_top, make_folds_checked, summarise_wf_results, fit_best_workflow, 
   make_recipe_full, make_recipe_top, 
   data_t1, folds_t1,
   rec_full_t1, rec_top20_t1, ncol_preds_t1, 
   data_t2, folds_t2,
   rec_full_t2, rec_top20_t2, ncol_preds_t2)


# Sauvegarde ----
saveRDS(
  results_proteomic_ALS_occurrence_tidymodels,
  file = "~/Documents/POP_ALS_2025_02_03/2_output/2.6.3_results_proteomic_ALS_occurrence_tidymodels.rds")



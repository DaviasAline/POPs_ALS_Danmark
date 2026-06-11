# XGBoost Hyperparameter Tuning avec tidymodels
# A. Davias
# June 4th 2026
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
library(tidymodels)
library(xgboost)
library(shapviz)

# XGBoost ----
## 1. Préparation des données ----
covar_xgboost <- c(proteomic, "birth_year", "sex", "follow_up_no_na_y",
                   "smoking_2cat_i", "bmi")

X_data <- bdd_danish |>
  select(all_of(covar_xgboost), match, als) |>
  mutate(
    sex            = as.numeric(sex),
    smoking_2cat_i = as.numeric(smoking_2cat_i),
    als            = factor(als, levels = c(0, 1), labels = c("No", "Yes")))

id_match <- X_data$match
Y_data   <- X_data$als
X_mat    <- X_data |> select(-match, -als)

cat("Dimensions X_mat:", nrow(X_mat), "x", ncol(X_mat), "\n")
cat("Distribution outcome:", table(Y_data), "\n")

## 2. Folds CV respectant la structure d'appariement ----
# group_vfold_cv garantit que les 2 membres d'une paire matched
# restent toujours dans le même fold → pas de data leakage

# On inclut match_id uniquement pour construire les folds, pas comme prédicteur
data_with_group <- bind_cols(X_mat, als = Y_data, match_id = id_match)
data_model      <- bind_cols(X_mat, als = Y_data)   # sans match_id → pour le modèle

set.seed(1996)
folds_tidy <- group_vfold_cv(
  data  = data_with_group,
  group = match_id,
  v     = 10)

# Vérification : tailles des folds
cat("Nb de folds:", nrow(folds_tidy), "\n")
cat("Taille fold train 1:", nrow(training(folds_tidy$splits[[1]])), "\n")
cat("Taille fold test 1:",  nrow(testing(folds_tidy$splits[[1]])),  "\n")

# Vérification leakage
for (i in seq_len(nrow(folds_tidy))) {
  train_matches <- training(folds_tidy$splits[[i]])$match_id
  test_matches  <- testing(folds_tidy$splits[[i]])$match_id
  leakage <- intersect(train_matches, test_matches)
  if (length(leakage) > 0) cat("Fold", i, ": LEAKAGE détecté !\n")
}
cat("Vérification leakage terminée\n")

## 3. Spécification du modèle ----
ncol_X <- ncol(X_mat)
cat("Nb de prédicteurs:", ncol_X, "\n")

xgb_spec <-
  boost_tree(
    trees      = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    mtry       = tune(),       # nb de colonnes tirées à chaque split
    min_n      = tune()) |>    # min observations par nœud terminal
  set_engine("xgboost", nthread = 1) |>
  set_mode("classification")

## 4. Grille d'hyperparamètres -----
tune_grid_tidy <- grid_regular(
  trees(range      = c(500, 1500)),
  tree_depth(range = c(2, 4)),
  learn_rate(range = c(-2, -1), trans = log10_trans()),  # 0.01 à 0.1
  mtry(range       = c(floor(ncol_X * 0.2),              # ~20% des colonnes
                       floor(ncol_X * 0.8))),             # ~80% des colonnes
  min_n(range      = c(1, 10)),
  levels = 3)   # 3 valeurs par paramètre → 3^5 = 243 combinaisons


cat("Nb de combinaisons à tester:", nrow(tune_grid_tidy), "\n")

## 5. Workflow ----
# add_formula(als ~ .) sur data_model qui ne contient PAS match_id
xgb_wf <- workflow() |>
  add_formula(als ~ .) |>
  add_model(xgb_spec)

## 6. Tuning ----
# set.seed(1996)
# xgb_tuned <- tune_grid(
#   xgb_wf,
#   resamples = folds_tidy,
#   grid      = tune_grid_tidy,
#   metrics   = metric_set(roc_auc),
#   control   = control_grid(
#     save_pred    = TRUE,
#     verbose      = TRUE))    # pour suivre la progression
  

#saveRDS(xgb_tuned, file = "~/Documents/POP_ALS_2025_02_03/2_output/xgb_tuned.rds")
xgb_tuned <- readRDS("~/Documents/POP_ALS_2025_02_03/2_output/xgb_tuned.rds")

## 7. Résultats du tuning ----
autoplot(xgb_tuned)
show_best(xgb_tuned, metric = "roc_auc", n = 10)

best_params <- select_best(xgb_tuned, metric = "roc_auc")
cat("\nMeilleurs hyperparamètres:\n")
print(best_params)

## 8. Early stopping pour affiner nrounds ----
# On reprend les meilleurs hyperparamètres et on laisse xgb.cv
# trouver le nb optimal de trees via early stopping
dtrain <- xgb.DMatrix(
  data  = as.matrix(X_mat),
  label = as.numeric(Y_data == "Yes"))

# Folds au format xgb.cv (liste d'indices 0-based)
xgb_folds <- lapply(seq_len(nrow(folds_tidy)), function(i) {
  test_idx <- which(id_match %in% testing(folds_tidy$splits[[i]])$match_id)
  test_idx
})

# Vérifier
cat("Fold 1 indices (premiers):", head(xgb_folds[[1]]), "\n")
cat("Max index fold 1:", max(xgb_folds[[1]]), "\n")
cat("N total:", nrow(dtrain), "\n")

set.seed(1996)
xgb_cv_es <- xgb.cv(
  params            = list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = best_params$tree_depth,
    eta              = best_params$learn_rate,
    colsample_bytree = best_params$mtry / ncol_X,
    min_child_weight = best_params$min_n,
    subsample        = 0.8),
  data                  = dtrain,
  nrounds               = 5000,
  folds                 = xgb_folds,
  early_stopping_rounds = 30,
  maximize              = TRUE,
  verbose               = 1)

#best_nrounds <- xgb_cv_es$best_iteration #proposé par claude mais ca n'a pas l'air de marcher
best_nrounds <- which.max(xgb_cv_es$evaluation_log$test_auc_mean)
cat("\nMeilleur nrounds avec early stopping:", best_nrounds, "\n")


# 1. Visualiser la courbe d'apprentissage pour comprendre la dynamique
xgb_cv_es$evaluation_log |>
  ggplot(aes(x = iter)) +
  geom_line(aes(y = train_auc_mean, color = "Train")) +
  geom_line(aes(y = test_auc_mean,  color = "Test")) +
  geom_ribbon(aes(ymin = test_auc_mean - test_auc_std,
                  ymax = test_auc_mean + test_auc_std), alpha = 0.2) +
  labs(title = "Learning curve XGBoost",
       x = "Nb de trees", y = "AUC", color = "") +
  theme_minimal()

# 2. AUC CV finale avec IC
best_row <- xgb_cv_es$evaluation_log[best_nrounds, ]
cat("AUC CV (mean ± SD):", round(best_row$test_auc_mean, 3), 
    "±", round(best_row$test_auc_std, 3), "\n")

# 3. Importance normalisée — juste les protéines
importance_matrix |>
  filter(Feature %in% proteomic) |>
  mutate(Gain_pct = Gain / sum(Gain) * 100) |>
  slice_head(n = 20) |>
  ggplot(aes(x = reorder(Feature, Gain_pct), y = Gain_pct)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Gain (%)", title = "Top 20 protéines - Variable importance") +
  theme_minimal()

## 9. Modèle final ----
xgb_final <- xgb.train(
  data    = dtrain,
  nrounds = best_nrounds,
  params  = list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = best_params$tree_depth,
    eta              = best_params$learn_rate,
    colsample_bytree = best_params$mtry / ncol_X,
    min_child_weight = best_params$min_n,
    subsample        = 0.8),
  verbose = 0)

# Importance des variables
importance_matrix <- xgb.importance(model = xgb_final)
xgb.plot.importance(importance_matrix, top_n = 250)

# Top 20 protéines les plus importantes
importance_matrix |>
  filter(Feature %in% proteomic) |>
  slice_head(n = 50) 

## 10. SHAP values plutot que gain ----
# Calcul des SHAP values
shap <- shapviz(xgb_final, X_pred = as.matrix(X_mat))       

# Beeswarm plot: Montre direction + magnitude pour chaque individu
sv_importance(shap, kind = "beeswarm", max_display = 20)         

# Bar plot classique (moyenne |SHAP|)
sv_importance(shap, kind = "bar", max_display = 20)          

# Dependence plot pour NEFL: Relation entre valeur de NEFL et sa contribution SHAP
sv_dependence(shap, v = "proteomic_neuro_explo_NEFL")                           

# Dependence plots top 6 protéines 
top6 <- importance_matrix$Feature[1:6]
for (prot in top6) {
  print(sv_dependence(shap, v = prot))
}

# ── 6. Waterfall pour cas individuel (ex: 1er cas ALS) ──
# Utile pour comprendre ce qui drive la prédiction d'un individu
idx_als <- which(Y_data == "Yes")[1]
sv_waterfall(shap, row_id = idx_als) +
  labs(title = paste("Patient ALS #", idx_als))

# SHAP moyen par protéine → tableau propre ──
shap_matrix <- shap$S   # matrice n x p de SHAP values

shap_summary <- tibble(
  Feature    = colnames(shap_matrix),
  mean_abs_shap = colMeans(abs(shap_matrix)),
  mean_shap     = colMeans(shap_matrix)) |>       # direction moyenne
  filter(Feature %in% proteomic) |>
  arrange(desc(mean_abs_shap)) |>
  mutate(
    direction = case_when(
      mean_shap > 0  ~ "↑ risque ALS",
      mean_shap < 0  ~ "↓ risque ALS",
      TRUE           ~ "neutre"))

print(shap_summary, n = 20)

# ── 1. Comprendre la bimodalité NEFL ──
# Regarder si le pattern SHAP correspond aux cas vs témoins
shap_nefl <- tibble(
  nefl_value = as.matrix(X_mat)[, "proteomic_neuro_explo_NEFL"],
  shap_value = shap$S[, "proteomic_neuro_explo_NEFL"],
  outcome    = Y_data,
  match_id   = id_match)

ggplot(shap_nefl, aes(x = nefl_value, y = shap_value, color = outcome)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "steelblue", "Yes" = "firebrick")) +
  labs(title = "NEFL : SHAP value selon outcome",
       x = "Valeur NEFL", y = "SHAP value") +
  theme_minimal()

# ── 2. Distribution de NEFL par outcome ──
ggplot(shap_nefl, aes(x = nefl_value, fill = outcome)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("No" = "steelblue", "Yes" = "firebrick")) +
  labs(title = "Distribution NEFL par outcome") +
  theme_minimal()

# ── 3. Interaction NEFL × follow-up ──
shap_interaction <- tibble(
  nefl  = as.matrix(X_mat)[, "proteomic_neuro_explo_NEFL"],
  follow_up  = as.matrix(X_mat)[, "follow_up_no_na_y"],
  shap_nefl = shap$S[, "proteomic_neuro_explo_NEFL"],
  outcome   = Y_data)

ggplot(shap_interaction, aes(x = nefl, y = shap_nefl, color = follow_up)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c() +
  facet_wrap(~outcome) +  # séparer cas et témoins
  labs(title = "NEFL SHAP selon outcome et follow-up duration") +
  theme_minimal()

# ── 8. Sauvegarder ──
results_proteomic_ALS_occurrence_y_to_als$logistic_models_including_time$machine_learning_xgboost$shap <- list(
  shapviz_obj  = shap,
  shap_matrix  = shap_matrix,
  shap_summary = shap_summary)

## 11. Sauvegarde ----
# results_proteomic_ALS_occurrence_y_to_als$logistic_models_including_time$machine_learning_xgboost <- list(
#   tuning        = xgb_tuned,
#   best_params   = best_params,
#   early_stop_cv = xgb_cv_es,
#   best_nrounds  = best_nrounds,
#   final_model   = xgb_final,
#   importance    = importance_matrix)

# Nettoyage
rm(covar_xgboost, X_data, X_mat, Y_data, id_match,
   data_with_group, data_model, folds_tidy,
   tune_grid_tidy, xgb_spec, xgb_wf,
   dtrain, xgb_folds, xgb_cv_es,
   best_params, best_nrounds, ncol_X, i)
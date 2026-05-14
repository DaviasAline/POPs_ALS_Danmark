# Test du package caret 
# May 12, 2026
# Aline Davias

# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
library(caret)


# pre selection of the proteins (identified with proteome wide associations) ----
significant_proteins <- 
  proteomic |> 
  str_subset(
    str_flatten(c("DFFA", "HNMT", "PIK3AP1", 
                  "PSIP1", "SH2D1A", "ANGPTL7", 
                  "APEX1", "ARG1", "CD164", 
                  "CLMP", "GLRX", "HDGF", 
                  "LRP11", "MCFD2", "NPTXR", 
                  "QDPR", "AKT1S1", "CD33", 
                  "IL32", "MAD1L1", "NEFL", 
                  "NPM1", "SRP14"), collapse = "|"))


# Sans matching -----
set.seed(1996)
bdd_ml_da <- bdd_danish |>
  haven::zap_labels() |>
  select(als, bmi, smoking_2cat_i, all_of(significant_proteins)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("Control", "Case"))) 

trainIndex <- createDataPartition(bdd_ml_da$als, p = 0.8, list = FALSE)
train_set  <- bdd_ml_da[trainIndex, ]
test_set   <- bdd_ml_da[-trainIndex, ]

fit_control <- trainControl(
  method = "repeatedcv", 
  number = 5, 
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)

# Liste des algorithmes de Chia et al. à tester
algos <- c("rf",        # random forest 
           "glmnet",    # elastic net/lasso/ridge
           "svmRadial", # support vector machine
           "gbm",       # gradient boosting machine
           "pls",       # partial least squares 
           "knn",       # k-nearest neighbors
           "nb",        # naive bayes
           "lda")       # linear discriminant analysis

results_list <- list()

for(a in algos) {
  message("Training model: ", a)
  
  set.seed(1996)
  
  # Entraînement avec centrage et réduction
  mod <- train(
    als ~ ., 
    data = train_set,
    method = a,
    trControl = fit_control,
    metric = "ROC",
    preProcess = c("center", "scale"))
  
  # Prédiction sur le Testing Set (les 20% restants)
  preds_prob <- predict(mod, test_set, type = "prob")
  preds_class <- predict(mod, test_set)
  
  # Calcul des métriques (AUC et Balanced Accuracy)
  roc_obj <- pROC::roc(test_set$als, preds_prob$Case, quiet = TRUE)
  cm <- confusionMatrix(preds_class, test_set$als)
  
  results_list[[a]] <- data.frame(
    Model = a,
    AUC = as.numeric(pROC::auc(roc_obj)),
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"])
}

# Résultat final
final_comparison <- do.call(rbind, results_list) |> arrange(desc(AUC))
print(final_comparison)



# GA
ctrl_ga <- gafsControl(functions = rfGA, # ou caretGA
                       method = "cv", 
                       number = 5,
                       allowParallel = TRUE)

set.seed(1996)
ga_obj <- gafs(x = train_set[, -1],  # On enlève la colonne 'als'
               y = train_set$als,
               iters = 25,           # Comme dans le papier (maxiter=25)
               popSize = 50,         # Réduit par rapport à leur 500 pour tester
               gafsControl = ctrl_ga,
               # Paramètres du modèle interne
               method = "glmnet", 
               trControl = trainControl(method = "cv", number = 3))

# Quelles protéines le GA a-t-il gardées ?
print(ga_obj$optVariables)

rm(trainIndex, train_set, test_set, fit_control, algos, a, mod, preds_prob, preds_class, roc_obj, cm, ctrl_ga)



# Sans matching mais avec ajustement -----
set.seed(1996)
bdd_ml_da_adj <- bdd_danish |>
  haven::zap_labels() |>
  select(als, bmi, smoking_2cat_i, sex, birth_year, all_of(significant_proteins)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("Control", "Case"))) 

trainIndex_adj <- createDataPartition(bdd_ml_da$als, p = 0.8, list = FALSE)
train_set_adj  <- bdd_ml_da[trainIndex_adj, ]
test_set_adj   <- bdd_ml_da[-trainIndex_adj, ]

fit_control_adj <- trainControl(
  method = "repeatedcv", 
  number = 5, 
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)

# Liste des algorithmes de Chia et al. à tester
algos <- c("rf", "glmnet", "svmRadial", "gbm", "pls", "knn", "nb", "lda")

results_list_adj <- list()

for(a in algos) {
  message("Training model: ", a)
  
  set.seed(1996)
  
  # Entraînement avec centrage et réduction
  mod <- train(
    als ~ ., 
    data = train_set_adj,
    method = a,
    trControl = fit_control_adj,
    metric = "ROC",
    preProcess = c("center", "scale"))
  
  # Prédiction sur le Testing Set (les 20% restants)
  preds_prob_adj <- predict(mod, test_set_adj, type = "prob")
  preds_class_adj <- predict(mod, test_set_adj)
  
  # Calcul des métriques (AUC et Balanced Accuracy)
  roc_obj <- pROC::roc(test_set_adj$als, preds_prob_adj$Case, quiet = TRUE)
  cm <- confusionMatrix(preds_class_adj, test_set_adj$als)
  
  results_list_adj[[a]] <- data.frame(
    Model = a,
    AUC = as.numeric(pROC::auc(roc_obj)),
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"]
  )
}

# Résultat final
final_comparison_adj <- do.call(rbind, results_list_adj) |> arrange(desc(AUC))
print(final_comparison_adj)


# GA
ctrl_ga_adj <- gafsControl(functions = rfGA, # ou caretGA
                       method = "cv", 
                       number = 5,
                       allowParallel = TRUE)

set.seed(1996)
ga_obj_adj <- gafs(x = train_set_adj[, -1],  # On enlève la colonne 'als'
               y = train_set_adj$als,
               iters = 25,           # Comme dans le papier (maxiter=25)
               popSize = 50,         # Réduit par rapport à leur 500 pour tester
               gafsControl = ctrl_ga,
               # Paramètres du modèle interne
               method = "glmnet", 
               trControl = trainControl(method = "cv", number = 3))

# Quelles protéines le GA a-t-il gardées ?
print(ga_obj_adj$optVariables)

rm(trainIndex_adj, train_set_adj, test_set_adj, fit_control_adj, algos, a, mod, preds_prob_adj, preds_class_adj, roc_obj, cm, ctrl_ga_adj)




# Avec matching ----
set.seed(1996)
unique_matches <- unique(bdd_danish$match)
train_match_ids <- sample(unique_matches, size = round(0.8 * length(unique_matches)))

train_set_match <- bdd_danish |> 
  haven::zap_labels() |>
  filter(match %in% train_match_ids) |>
  select(als, bmi, smoking_2cat_i, all_of(significant_proteins)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("Control", "Case"))) 

test_set_match <- bdd_danish |> 
  haven::zap_labels() |>
  filter(!(match %in% train_match_ids)) |>
  select(als, bmi, smoking_2cat_i, all_of(significant_proteins)) |>
  mutate(als = factor(als, levels = c(0, 1), labels = c("Control", "Case")))

# Group-based CV
train_match_vector <- bdd_danish$match[bdd_danish$match %in% train_match_ids]
folds <- groupKFold(train_match_vector, k = 5)

fit_control_match <- trainControl(
  index = folds, method = "cv", 
  classProbs = TRUE, summaryFunction = twoClassSummary)

# Liste des algorithmes de Chia et al. à tester
algos <- c("rf", "glmnet", "svmRadial", "gbm", "pls", "knn", "nb", "lda")


# Initialisation du tableau de résultats
results_list_match <- list()

for(a in algos) {
  message("Training model: ", a)
  
  set.seed(1996)
  # Entraînement avec centrage et réduction 
  mod <- train(
    als ~ ., 
    data = train_set_match,
    method = a,
    trControl = fit_control_match,
    metric = "ROC",
    preProcess = c("center", "scale"))
  
  # Prédiction sur le Testing Set (20% restants)
  preds_prob <- predict(mod, test_set_match, type = "prob")
  preds_class <- predict(mod, test_set_match)
  
  # Calcul des métriques (AUC et Balanced Accuracy)
  roc_obj <- pROC::roc(test_set_match$als, preds_prob$Case, quiet = TRUE)
  cm <- confusionMatrix(preds_class, test_set_match$als)
  
  results_list_match[[a]] <- data.frame(
    Model = a,
    AUC = as.numeric(pROC::auc(roc_obj)),
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"]
  )
}

# Résultat final
final_comparison_match <- do.call(rbind, results_list_match) |> arrange(desc(AUC))
print(final_comparison_match)

# genetic algorithm
ctrl_ga_match <- gafsControl(functions = rfGA, # ou caretGA
                           method = "cv", 
                           number = 5,
                           allowParallel = TRUE)

set.seed(1996)
ga_obj_match <- gafs(x = train_set_match[, -1],  # On enlève la colonne 'als'
                   y = train_set_match$als,
                   iters = 25,           # Comme dans Chia et al  (maxiter=25)
                   popSize = 50,         # Réduit par rapport à leur 500 pour tester
                   gafsControl = ctrl_ga,
                   # Paramètres du modèle interne
                   method = "glmnet", 
                   trControl = trainControl(method = "cv", number = 3))

# Quelles protéines le GA a-t-il gardées ?
print(ga_obj_match$optVariables)


rm(unique_matches, train_match_ids, train_set_match, test_set_match, train_match_vector, folds, fit_control_match, algos, a, mod, preds_prob, preds_class, roc_obj, cm, ctrl_ga_match)


# work on machine learning 
# 04/12/2025

# supervised machine learning : random forest method ----
## data preparation ----
proteomic_selected <-                                                           # selection des prot avec p<0.05
  results_proteomic_ALS_occurrence$main$main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         p_value_raw<0.05) |> 
  pull(explanatory) 

proteomic_selected <-                                                           # selection des noms des var complets
  proteomic[
  sapply(proteomic_selected, \(g) grepl(g, proteomic)) |> rowSums() > 0
]


bdd_danish |>                                                                   # visu data
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, 
         all_of(proteomic_selected)) |>
  str()

bdd_rf <-                                                                       # creation df
  bdd_danish |>  
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, 
         all_of(proteomic_selected)) |>
  mutate(
    als = factor(as.character(als)),        
    match = factor(as.character(match)))



## unconditional random forest -----
### chargement des packages ----
library(randomForest)
library(caret)

### data preparation ----
set.seed(1996)
matches <- unique(bdd_rf$match)                                                 # 70% for training, 30% for testing
train_matches <- sample(matches, size = floor(0.7 * length(matches)))

train_data <- bdd_rf |> filter(match %in% train_matches) |> select(-match)      # we remove match because we are already adjusting on baseline age and sex
test_data  <- bdd_rf |> filter(!match %in% train_matches) |> select(-match)     # we remove match because we are already adjusting on baseline age and sex

rm(matches, train_matches)

### random forest with randomForest() ----
set.seed(1996)

rf_model <- randomForest(
  als ~ .,
  data = train_data,
  ntree = 500,           # nombre d’arbres
  mtry = floor(sqrt(ncol(train_data)-1)), # nombre de variables testées à chaque split
  importance = TRUE)

print(rf_model)

### importance des variables ----
importance(rf_model)
varImpPlot(rf_model)

### performance predictive ----
pred <- predict(rf_model, newdata = test_data, type = "response")
confusionMatrix(pred, test_data$als)



## conditionnal random forest ----
### chargement des packages ----
library(party)
library(pROC)
library(groupdata2)   # pour CV / split groupé

### data preparation ----
set.seed(1996)
matches <- unique(bdd_rf$match)                                                 # 70% for training, 30% for testing
train_matches <- sample(matches, size = floor(0.7 * length(matches)))

train_data <- 
  bdd_rf |> filter(match %in% train_matches) |> select(-baseline_age, -sex)     # we remove baseline age and sex because we are already matching on it with variable match
test_data  <- 
  bdd_rf |> filter(!match %in% train_matches) |> select(-baseline_age, -sex)    # we remove baseline age and sex because we are already matching on it with variable match

rm(matches, train_matches)


### conditional random forest with cforest() ----
set.seed(1996)

cforest_ctrl <- cforest_unbiased(
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data) - 2)))                                     # nb de variables testées à chaque split

rf_model <- cforest(
  als ~ .,                                                                      # car train_data n’a que les variables utiles
  data = train_data,
  controls = cforest_ctrl)

### importance des variables ----
varimp_cond <- varimp(rf_model, conditional = TRUE)

importance_sorted <- sort(varimp_cond, decreasing = TRUE)                       # Classement décroissant
importance_sorted

### performance predictive ----
pred_prob <- predict(rf_model, newdata = test_data, type = "prob")              # probabilités
prob_als1 <- pred_prob[, "1"]

pred_class <- ifelse(prob_als1 > 0.5, 1, 0)                                     # classes predites
pred_class <- factor(pred_class)

roc_obj <- roc(test_data$als, prob_als1)                                        # AUC
auc_value <- auc(roc_obj)
auc_value

plot(roc_obj, main = sprintf("ROC Curve - AUC = %.3f", auc_value))              # Courbe ROC

accuracy <- mean(pred_class == test_data$als)                                   # Accuracy
accuracy



## conditional random forest avec 10 fold cross validation ----
### chargement des packages ----
library(party)
library(pROC)
library(groupdata2)   # group K-fold CV

### data preparation ----
bdd_rf <- 
  bdd_danish |>  
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, 
         all_of(proteomic_selected)) |>
  mutate(
    als   = factor(as.character(als)),        
    match = factor(as.character(match)))





### 10-fold group cross-validation (group = match) ----
set.seed(1996)
bdd_cv <- fold(                                     # groupdata2 construit une colonne .folds qui respecte les groupes
  data   = bdd_rf |> select(-baseline_age, -sex),   # on retire ces deux var car on match avec var match
  k      = 10,
  id_col = "match",    # groupement par match
  method = "n_dist")


### Boucle CV sur les 10 folds ----
auc_values <- c()                                  # Stockage des métriques
accuracy_values <- c()

for (fold in 1:10) {
  
  cat("Running fold", fold, "...\n")
  
  train_data <- bdd_cv |> filter(.folds != fold)
  test_data  <- bdd_cv |> filter(.folds == fold)
  
  cforest_ctrl <- cforest_unbiased(            # Construction du modèle
    ntree = 500,
    mtry = floor(sqrt(ncol(train_data) - 2)))   # -2 car als + match
  
  rf_model <- cforest(
    als ~ .,
    data = train_data |> select(-.folds),
    controls = cforest_ctrl)
  
  pred_prob <- predict(rf_model, newdata = test_data, type = "prob")            # Prédictions probabilistes
  prob_als1 <- pred_prob[, "1"]
  
  roc_obj <- roc(test_data$als, prob_als1)                                      # AUC
  auc_values <- c(auc_values, auc(roc_obj))
  
  pred_class <- factor(ifelse(prob_als1 > 0.5, 1, 0))                           # Accuracy
  accuracy_values <- c(accuracy_values, mean(pred_class == test_data$als))
}


### Results ----
mean_auc <- mean(auc_values)
sd_auc   <- sd(auc_values)

mean_acc <- mean(accuracy_values)
sd_acc   <- sd(accuracy_values)

cat("\n=== Cross-Validated Performance ===\n")
cat(sprintf("AUC (mean ± SD): %.3f ± %.3f\n", mean_auc, sd_auc))
cat(sprintf("Accuracy (mean ± SD): %.3f ± %.3f\n", mean_acc, sd_acc))

### Importance des variables sur modèle final (entraînement complet) -----
cforest_ctrl_final <- cforest_unbiased(
  ntree = 500,
  mtry = floor(sqrt(ncol(bdd_rf) - 2)))

rf_model_final <- cforest(
  als ~ .,
  data = bdd_rf |> select(-baseline_age, -sex),
  controls = cforest_ctrl_final)

# Importance conditionnelle (permutation)
varimp_cond <- varimp(rf_model_final, conditional = TRUE)
importance_sorted <- sort(varimp_cond, decreasing = TRUE)

importance_sorted


# unsupervised machine learning : K-mean clustering method ----

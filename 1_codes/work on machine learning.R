# travail de machine learning pour predire le risque de ALS ----
# 04/12/2025

proteomic_selected <- results_proteomic_ALS_occurrence$main$main_results |> 
  filter(analysis == "main", 
         model == "adjusted", 
         term == "Continuous", 
         p_value_raw<0.05) |> 
  pull(explanatory) 

proteomic_selected <- proteomic[
  sapply(proteomic_selected, \(g) grepl(g, proteomic)) |> rowSums() > 0
]


bdd_danish |> 
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, 
         all_of(proteomic_selected)) |>
  str()

bdd_rf <- bdd_danish |>
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, all_of(proteomic_selected)) |>
  mutate(
    als = factor(as.character(als)),        
    match = factor(as.character(match)))



## unconditional random forest -----
library(randomForest)
library(caret)


set.seed(1996)

# 80% pour training, 20% pour testing
matches <- unique(bdd_rf$match)
train_matches <- sample(matches, size = floor(0.8 * length(matches)))

train_data <- bdd_rf |> filter(match %in% train_matches) |> select(-match)
test_data  <- bdd_rf |> filter(!match %in% train_matches) |> select(-match)


set.seed(1996)

rf_model <- randomForest(
  als ~ .,
  data = train_data,
  ntree = 500,           # nombre d’arbres
  mtry = floor(sqrt(ncol(train_data)-1)), # nombre de variables testées à chaque split
  importance = TRUE
)

print(rf_model)

pred <- predict(rf_model, newdata = test_data, type = "response")
confusionMatrix(pred, test_data$als)

# Importance globale
importance(rf_model)
varImpPlot(rf_model)




## conditional random forest ----
library(party)

set.seed(1996)

cforest_ctrl <- cforest_unbiased(
  ntree = 500,               # nombre d’arbres
  mtry = floor(sqrt(ncol(bdd_rf) - 2))  # nombre de variables testées à chaque split
)

rf_model <- cforest(
  als ~ baseline_age + sex + smoking_2cat_i + bmi + .,
  data = bdd_rf,
  controls = cforest_ctrl,
  # on inclut le facteur 'match' pour l'appariement
  weights = NULL
)

pred_prob <- predict(rf_model, type = "prob")
# probabilité que als = 1
pred_class <- factor(ifelse(pred_prob[, "1"] > 0.5, 1, 0))

var_importance <- varimp(rf_model, conditional = TRUE)
var_importance_sorted <- sort(var_importance, decreasing = TRUE)
print(var_importance_sorted)



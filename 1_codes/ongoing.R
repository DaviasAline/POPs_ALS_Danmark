
# Loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
library(survival)
library(glmnet)
library(penalized)
library(survminer)
library(xgboost)

# Cox x Xgboost ----
X_matrix <- bdd_danish |> 
  select(birth_year, sex, bmi, smoking_2cat_i, all_of(proteomic)) |>
  mutate(
    sex_male = ifelse(sex == "Male", 1, 0),
    smoking_ever = ifelse(smoking_2cat_i == "Ever", 1, 0)) |>
  select(-sex, -smoking_2cat_i) |> 
  as.matrix()

y_survival <- ifelse(bdd_danish$als == 1, bdd_danish$follow_up_bis, -bdd_danish$follow_up_bis)
dtrain <- xgb.DMatrix(data = X_matrix, label = y_survival)

# Hyperparamètres ultra-régularisés pour éviter le surapprentissage immédiat
params <- list(
  objective        = "survival:cox",
  eval_metric      = "cox-nloglik",
  eta              = 0.01,           # On ralentit l'apprentissage par 5
  max_depth        = 2,              # Des arbres très simples (souches) pour commencer
  subsample        = 0.7,            # Sous-échantillonnage des lignes
  colsample_bytree = 0.3,            # On ne teste que 30% des protéines par arbre (évite qu'une protéine écrase tout)
  alpha            = 1,              # Régularisation L1 (Lasso-like pour XGBoost)
  lambda           = 1)              # Régularisation L2 (Ridge-like pour XGBoost)

set.seed(1996)
cv_model <- xgb.cv(
  params                = params,
  data                  = dtrain,
  nrounds               = 1000,          # On augmente le max car eta est plus petit
  nfold                 = 5,
  early_stopping_rounds = 50,          # Plus de patience pour laisser le modèle converger
  verbose               = TRUE,
  print_every_n         = 20)

# Si le modèle préfère toujours l'itération 0, on force au moins 15 arbres
best_nrounds <- cv_model$best_iteration
if(is.null(best_nrounds) || best_nrounds == 0) {
  best_nrounds <- 15
  cat("\n[Info] Le modèle a encore du mal à converger en CV. On force à 15 arbres pour l'extraction.\n")
} else {
  cat("\n---> Le nombre optimal d'arbres retenu est de :", best_nrounds, "\n\n")
}

# Entraînement final 
final_xgb_cox <- xgb.train(
  params  = params,
  data    = dtrain,
  nrounds = best_nrounds)

# Extraction de l'importance
importance_matrix <- xgb.importance(feature_names = colnames(X_matrix), model = final_xgb_cox)
print(head(importance_matrix, 20))

# Graphique d'importance 
xgb.plot.importance(
  importance_matrix = importance_matrix[1:20, ], 
  main = "Top 20 de l'importance des variables (Cox XGBoost)")

rm(X_matrix, y_survival, dtrain, params, best_nrounds)

# Logistic x Xgboost ----
X_matrix_binary <- bdd_danish |> 
  select(follow_up_bis, birth_year, sex, bmi, smoking_2cat_i, all_of(proteomic)) |>
  mutate(
    sex_male = ifelse(sex == "Male", 1, 0),
    smoking_ever = ifelse(smoking_2cat_i == "Ever", 1, 0)) |>
  select(-sex, -smoking_2cat_i) |> 
  as.matrix()

y_binary <- bdd_danish$als
dtrain_binary <- xgb.DMatrix(data = X_matrix_binary, label = y_binary)          # Encapsulation dans le format optimisé d'XGBoost

# Paramètres logistiques stabilisés pour la haute dimension
params_binary <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",            # Notre métrique cible est l'Aire Sous la Courbe
  eta              = 0.01,             # Apprentissage très progressif
  max_depth        = 3,                # Profondeur de 3 pour capter les interactions complexes
  subsample        = 0.8,              # 80% des lignes par arbre
  colsample_bytree = 0.4,              # 40% des protéines tirées au sort par arbre
  alpha            = 1,                # Pénalisation L1 (Lasso)
  lambda           = 1)                # Pénalisation L2 (Ridge)


set.seed(1996)
cv_binary <- xgb.cv(
  params                = params_binary,
  data                  = dtrain_binary,
  nrounds               = 1000, 
  nfold                 = 5,           # 5 folds
  early_stopping_rounds = 50,          # Arrêt dès que l'AUC de test stagne pendant 50 arbres
  verbose               = TRUE,
  print_every_n         = 20)

best_nrounds_binary <- cv_binary$best_iteration

# Sécurité mathématique : si l'index renvoie 0 ou NULL à cause du codage interne d'XGBoost,
# on extrait manuellement l'index exact où l'AUC de test a atteint son maximum historique.
if (is.null(best_nrounds_binary) || length(best_nrounds_binary) == 0 || best_nrounds_binary == 0) {
  best_nrounds_binary <- which.max(cv_binary$evaluation_log$test_auc_mean)
}

cat("\n---> Nombre optimal d'arbres retenu pour le modèle binaire :", best_nrounds_binary, "\n\n")


# Entraînement final du classifieur sur l'ensemble du dataset & IMPORTANCE DES VARIABLES
final_xgb_binary <- xgb.train(
  params  = params_binary,
  data    = dtrain_binary,
  nrounds = best_nrounds_binary)

# Extraction et calcul des gains d'importance pour chaque variable
importance_binary <- xgb.importance(feature_names = colnames(X_matrix_binary), model = final_xgb_binary)

# Top 20
cat("Top 20 des variables les plus prédictives du statut ALS (modèle binaire) :\n")
print(head(importance_binary, 20))

# Graphique d'importance 
xgb.plot.importance(
  importance_matrix = importance_binary[1:20, ], 
  main = "Top 20 de l'importance des variables (XGBoost Logistique)")

rm(X_matrix_binary, y_binary, dtrain_binary, params_binary, best_nrounds_binary)

# Cox x lasso ----
## Model 1 with matching ----
fit_penalized <- penalized(
  response = surv_obj,                         
  penalized = bdd_danish[, proteomic_sd],                      
  unpenalized = ~ bmi + smoking_2cat_i + strata(match), 
  lambda1 = 1, lambda2 = 1,                         # Équivalent Elastic Net (Lasso + Ridge)
  data = bdd_danish)


## Model 2 with adjustment ----
### model ----
X_clinical <- model.matrix(~ bmi + smoking_2cat_i + birth_year + sex, data = bdd_danish)[, -1]
X_proteins <- as.matrix(bdd_danish[, proteomic_sd])
storage.mode(X_proteins) <- "double" 
X_matrix_adj <- cbind(X_clinical, X_proteins)
n_clinical_vars <- ncol(X_clinical) # Dynamique selon le nombre de niveaux de 'sex'
n_protein_vars  <- ncol(X_proteins)  # Égal à 273
p.fac <- c(rep(0, n_clinical_vars), rep(1, n_protein_vars))

set.seed(1996) 
cv_fit_adj <- cv.glmnet(
  x = X_matrix_adj, 
  y = surv_obj, 
  family = "cox",       
  alpha = 0.5,             # 0.5 = Elastic Net (Équilibre entre sélection Lasso et corrélation Ridge)
  penalty.factor = p.fac,  # Force les covariables à etre retenues
  nfolds = 10, 
  type.measure = "C")            

## Graphique de la validation croisée (Évolution de l'erreur selon lambda) ----
plot(cv_fit_adj)

## Performance globale (C-index de Harrell) ----
c_index <- max(cv_fit_adj$cvm)
cat("Performance globale du modèle (C-index cross-validé) :", round(c_index, 4), "\n")

## Extraction des coefficients non nuls au Lambda Optimal ----
coefficients <- coef(cv_fit_adj, s = "lambda.min") |> as.matrix()
selected_features <- coefficients[coefficients[, 1] != 0, , drop = FALSE]

print("Variables et protéines retenues dans le modèle prédictif final :")
print(selected_features)

df_plot <- as.data.frame(selected_features)
df_plot$Variable <- rownames(df_plot)
colnames(df_plot)[1] <- "Coefficient"

ggplot(df_plot, aes(x = reorder(Variable, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = ifelse(df_plot$Coefficient > 0, "firebrick", "steelblue")) +
  coord_flip() +
  labs(title = "Signature Protéique Prédictive de la SLA", x = "Protéines / Covariables", y = "Poids dans le modèle (Coefficient)") +
  theme_minimal()

### Courbes de survie ----
# On crée deux patients types (lignes) avec les mêmes colonnes que X_matrix_adj
new_profiles <- matrix(0, nrow = 2, ncol = ncol(X_matrix_adj))
colnames(new_profiles) <- colnames(X_matrix_adj)
rownames(new_profiles) <- c("Profil Moyen", "Profil Haut Risque")

# Remplissage du profil Moyen (Standard) avec les moyennes cliniques
new_profiles[, "bmi"] <- mean(bdd_danish$bmi)
new_profiles[, "birth_year"] <- mean(bdd_danish$birth_year)

# Remplissage du profil Haut Risque (on force les biomarqueurs clés trouvés)
# Par exemple, on simule un Neurofilament (NEFL) à +2 écart-types
if("proteomic_neuro_explo_NEFL_sd" %in% colnames(X_matrix_adj)) {
  new_profiles["Profil Haut Risque", "proteomic_neuro_explo_NEFL_sd"] <- 2
}
# On peut ajouter d'autres protéines ici (ex: CD33 = 1.5, IL32 = -1.5)

# Calcul des trajectoires de survie pour ces profils
fit_plot <- survfit(cv_fit_adj, s = "lambda.min", x = X_matrix_adj, y = surv_obj, newx = new_profiles)

# Graphique final avec survminer
ggsurvplot(fit_plot, 
           data = bdd_danish,
           palette = c("#2E9FDF", "#E7B800"),
           legend.labs = c("Individu Moyen", "Individu à Risque (NEFL+)"),
           xlab = "Années depuis la Baseline", 
           ylab = "Probabilité d'être sans diagnostic de SLA",
           ggtheme = theme_minimal(),
           title = "Prédiction de l'incidence de la SLA (Elastic Net Ajusté)")




# Prise en compte differences de proteins EPIC vs UK biobank ----
library(readxl)
library(tidyverse)
library(questionr)
library(gtsummary)
Olink_explore <- read_excel("~/Documents/Appels à projets/9. UK biobank access 2026/Olink Explore protein list.xltx")
Olink_target <- read_excel("~/Documents/Appels à projets/9. UK biobank access 2026/Olink Target 96 protein list.xltx")


test <- left_join(Olink_target, Olink_explore, by = c("UniProt ID", "Gene", "Protein name"))


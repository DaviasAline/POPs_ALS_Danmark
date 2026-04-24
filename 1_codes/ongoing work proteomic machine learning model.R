# Charger les bibliothèques nécessaires
library(glmnet)
library(survival)
library(dplyr)
library(pROC)


# Lasso ----
## Preparation des données ----
# Création de la matrice X (prédicteurs)
# Le "- 1" supprime l'intercept car glmnet le gère automatiquement
X_full <- model.matrix(~ bmi + smoking_2cat_i + . - 1, 
                       data = bdd_danish[, c("bmi", "smoking_2cat_i", proteomic_sd)])

# Création de l'objet de réponse (Surv)
# On met un temps fictif de 1 pour tout le monde (standard pour le logit conditionnel)
y_cox <- Surv(time = rep(1, nrow(bdd_danish)), event = bdd_danish$als)

## Base model lasso ----
# Extraction uniquement des colonnes de protéines
X_prot <- as.matrix(bdd_danish[, proteomic_sd])

set.seed(123) # Pour la reproductibilité de la cross-validation
cv_m1 <- cv.glmnet(
  x = X_prot, 
  y = y_cox, 
  family = "cox", 
  strata = bdd_danish$match, # Indispensable pour respecter l'appariement
  alpha = 1                  # Lasso
)

plot(cv_m1)
title("Modèle 1 : Protéines seules", line = 3)


# Extraire les protéines sélectionnées par le Lasso (Modèle 1)
# On utilise lambda.min (le plus précis) ou lambda.1se (le plus parcimonieux)
coefs <- coef(cv_m1, s = "lambda.min")
selected_vars <- row.names(coefs)[as.matrix(coefs) != 0]

# Retirer les covariables de la liste pour n'avoir que les protéines
selected_proteins_1 <- setdiff(selected_vars, colnames(X_full)[idx_covar])

cat("Nombre de protéines sélectionnées :", length(selected_proteins_1), "\n")
print(selected_proteins_1)


# Une fois le Lasso terminé, on repasse en clogit classique pour les stats
if(length(selected_proteins) > 0) {
  
  formula_final <- as.formula(paste(
    "als ~ bmi + smoking_2cat_i +", 
    paste(selected_proteins, collapse = " + "), 
    "+ strata(match)"
  ))
  
  fit_final <- clogit(formula_final, data = bdd_danish)
  
  print(summary(fit_final))
  
} else {
  message("Aucune protéine n'a été sélectionnée par le Lasso.")
}

# 1. Calculer le score de risque pour le Modèle 1 (X_prot uniquement)
# On applique directement le signe (-) pour redresser la prédiction
risk_scores_m1 <- -as.numeric(predict(cv_m1, newx = X_prot, s = "lambda.min", type = "link"))

# 2. Calculer la concordance (C-index) avec l'appariement
perf_m1 <- concordance(y_cox ~ risk_scores_m1 + strata(bdd_danish$match))

cat("--- RÉSULTATS MODÈLE 1 (PROTÉINES SEULES) ---\n")
print(perf_m1)

# 3. Générer la courbe ROC pour le Modèle 1
roc_m1 <- roc(bdd_danish$als, risk_scores_m1, quiet = TRUE)

plot(roc_m1, col = "blue", lwd = 2, 
     main = paste("Comparaison des Modèles\nAUC M1 (Bleu):", round(auc(roc_m1), 3)))



## Adjusted model lasso ----
# On définit le vecteur de pénalité dynamiquement
p.fac <- rep(1, ncol(X_full))

# On identifie les colonnes à NE PAS pénaliser (BMI et Smoking)
# On cherche les colonnes créées par model.matrix qui contiennent ces noms
idx_covar <- grep("bmi|smoking", colnames(X_full))
p.fac[idx_covar] <- 0

set.seed(123)
cv_m2 <- cv.glmnet(
  x = X_full, 
  y = y_cox, 
  family = "cox", 
  strata = bdd_danish$match,
  penalty.factor = p.fac,    # 0 pour les covariables, 1 pour les protéines
  alpha = 1
)

plot(cv_m2)
title("Modèle 2 : Ajusté (BMI + Smoking)", line = 3)


# Extraire les protéines sélectionnées par le Lasso (Modèle 2)
# On utilise lambda.min (le plus précis) ou lambda.1se (le plus parcimonieux)
coefs <- coef(cv_m2, s = "lambda.min")
selected_vars <- row.names(coefs)[as.matrix(coefs) != 0]

# Retirer les covariables de la liste pour n'avoir que les protéines
selected_proteins_2 <- setdiff(selected_vars, colnames(X_full)[idx_covar])

cat("Nombre de protéines sélectionnées :", length(selected_proteins_2), "\n")
print(selected_proteins_2)


# Une fois le Lasso terminé, on repasse en clogit classique pour les stats
if(length(selected_proteins) > 0) {
  
  formula_final <- as.formula(paste(
    "als ~ bmi + smoking_2cat_i +", 
    paste(selected_proteins, collapse = " + "), 
    "+ strata(match)"
  ))
  
  fit_final <- clogit(formula_final, data = bdd_danish)
  
  print(summary(fit_final))
  
} else {
  message("Aucune protéine n'a été sélectionnée par le Lasso.")
}

# 1. Calculer le score de risque (linear predictor) pour chaque individu
# On utilise le modèle cv_m2 avec le meilleur lambda
risk_scores <- predict(cv_m2, newx = X_full, s = "lambda.min", type = "link")

# 2. Calculer la concordance en tenant compte de l'appariement

test_perf <- concordance(y_cox ~ risk_scores + strata(bdd_danish$match))
print(test_perf)


# Inversion du signe pour que le score suive le risque ALS
real_risk_scores <- -as.numeric(predict(cv_m2, newx = X_full, s = "lambda.min", type = "link"))

# Recalcul de la concordance
final_perf <- concordance(y_cox ~ real_risk_scores + strata(bdd_danish$match))
print(final_perf)


# On crée un objet ROC classique sur les scores de risque
roc_obj <- roc(bdd_danish$als, as.numeric(real_risk_scores))

plot(roc_obj, main = paste("Courbe ROC - AUC:", round(auc(roc_obj), 3)))

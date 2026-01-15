

# Grouped (quartiles) elastic net regression (for ungroup POPs) -----
## data preparation ----
bdd_cases_danish_bis <- 
  bdd_danish |>
  filter (als == 1) |>
  filter(follow_up_death>0) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_group)) |>
  mutate(across(all_of(POPs_group), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))

POPs_group_quart_bis <- setdiff(POPs_group_quart, "PCB_4_quart")    
covariates_danish <- c("sex", "diagnosis_age", "smoking_2cat_i", "bmi", "marital_status_2cat_i")



# Créer les dummies tidyverse-style pour chaque variable
#    - Les expositions (quartiles) sont encodées en one-hot 4 colonnes
#    - Les covariables forcées sont gardées telles quelles
#    - On obtient X_expo et X_cov forcé

# Expositions (chaque quartile -> one-hot)
X_expo <- bdd_cases_danish_bis |>                                               # creer des dummies pour avoir trois varriables par expo au lieu de une variable à 4 catégories                           
  select(all_of(POPs_group_quart_bis)) |>
  fastDummies::dummy_cols(remove_first_dummy = TRUE, 
                          remove_selected_columns = TRUE) |>
  as_tibble()

# Covariables forcées
X_cov <- bdd_cases_danish_bis |>
  select(all_of(covariates_danish)) |>
  mutate(across(where(is.factor), ~ as.numeric(. != levels(.)[1]))) |>          # convertit 2-level factors en 0/1
  as_tibble()

# Fusionner matrices
X <- bind_cols(X_expo, X_cov) |> as.matrix()

## Création des groupes de variables ----
# - Chaque exposition = un groupe (4 variables d'u même facteur'une meme expo)
# - Covariables = groupes de variables non pénalisés

# Nombre de colonnes par exposition
ncol_expo_each <-                                                               # 3 par variable
  X_expo |>                        
  select(starts_with(POPs_group_quart_bis[1])) |>
  ncol() 

group_sizes <-                                                                  # les covariables en bloc à la fin
  c(rep(ncol_expo_each, length(POPs_group_quart_bis)),
    ncol(X_cov))  

group_starts <- cumsum(c(1, head(group_sizes, -1)))

# Penalité par groupe : 1 pour expos, 0 pour covariables
penalty_vec <- c(rep(1, length(POPs_group_quart_bis)), rep(0, 1))

## Modele de cox -----
glm_family <- 
  glm.cox(stop = bdd_cases_danish_bis$follow_up_death, status = bdd_cases_danish_bis$als)

set.seed(1996)
cvfit <- cv.grpnet(
  X,
  glm_family,
  groups = group_starts,
  penalty = penalty_vec,
  alpha = 0.1,                                                                  # ranging from 0 (ridge) to 1 (lasso)
  n_folds = 5,
  standardize = FALSE)

## Résultats -----
print(cvfit)
plot(cvfit)
coef(cvfit)


# Extraire les coefficients au lambda optimal
coef_min <- coef(cvfit, lambda = "lambda.min")

# Transformer en vecteur nommé
beta_vec <- as.numeric(coef_min$betas)
names(beta_vec) <- colnames(X) 

# Créer un tibble propre
coef_tbl <- tibble(
  variable = names(beta_vec),
  beta = beta_vec)

# Filtrer les variables sélectionnées (non nulles)
selected_vars <- 
  coef_tbl |>
  filter(beta != 0) |>
  arrange(desc(abs(beta))) |>
  mutate(
    HR = exp(beta))

selected_vars

rm(selected_vars, coef_min, cvfit, coef_tbl, beta_vec, glm_family, penalty_vec, group_starts, group_sizes, ncol_expo_each, X, X_cov, X_expo, dat, covariates_danish, POPs_group_quart_bis, bdd_cases_danish_bis)



# Grouped (quartiles) elastic net regression (for ungrouped POPs) -----
## data preparation ----
bdd_cases_danish_bis <- 
  bdd_danish |>
  filter (als == 1) |>
  filter(follow_up_death>0) |>
  filter(study == "Danish") |>
  select(study, als, follow_up_death, status_death, sex, baseline_age, diagnosis_age, death_age,
         bmi, marital_status_2cat_i, smoking_i, smoking_2cat_i, education_i, cholesterol_i, 
         all_of(POPs_included)) |>
  mutate(across(all_of(POPs_included), ~ factor(ntile(.x, 4),                      # creation of POPs quartiles (cohort and cases specific)                        
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat_i = fct_relevel(smoking_2cat_i, "Ever", "Never"), 
         marital_status_2cat_i = fct_relevel(marital_status_2cat_i, "Married/cohabit", "Other"))

  
covariates_danish <- c("sex", "diagnosis_age", "smoking_2cat_i", "bmi", "marital_status_2cat_i")


# Créer les dummies tidyverse-style pour chaque variable
#    - Les expositions (quartiles) sont encodées en one-hot 4 colonnes
#    - Les covariables forcées sont gardées telles quelles
#    - On obtient X_expo et X_cov forcé

# Expositions (chaque quartile -> one-hot)
X_expo <- bdd_cases_danish_bis |>                                               # creer des dummies pour avoir trois varriables par expo au lieu de une variable à 4 catégories                           
  select(all_of(POPs_included_quart)) |>
  fastDummies::dummy_cols(remove_first_dummy = TRUE, 
                          remove_selected_columns = TRUE) |>
  as_tibble()

# Covariables forcées
X_cov <- bdd_cases_danish_bis |>
  select(all_of(covariates_danish)) |>
  mutate(across(where(is.factor), ~ as.numeric(. != levels(.)[1]))) |>          # convertit 2-level factors en 0/1
  as_tibble()

# Fusionner matrices
X <- bind_cols(X_expo, X_cov) |> as.matrix()

## Création des groupes de variables ----
# - Chaque exposition = un groupe (4 variables d'u même facteur'une meme expo)
# - Covariables = groupes de variables non pénalisés

# Nombre de colonnes par exposition
ncol_expo_each <-                                                               # 3 par variable
  X_expo |>                        
  select(starts_with(POPs_included_quart[1])) |>
  ncol() 

group_sizes <-                                                                  # les covariables en bloc à la fin
  c(rep(ncol_expo_each, length(POPs_included_quart)),
    ncol(X_cov))  

group_starts <- cumsum(c(1, head(group_sizes, -1)))

# Penalité par groupe : 1 pour expos, 0 pour covariables
penalty_vec <- c(rep(1, length(POPs_included_quart)), rep(0, 1))

## Modele de cox -----
glm_family <- 
  glm.cox(stop = bdd_cases_danish_bis$follow_up_death, status = bdd_cases_danish_bis$als)

set.seed(1996)
cvfit <- cv.grpnet(
  X,
  glm_family,
  groups = group_starts,
  penalty = penalty_vec,
  alpha = 0.1,                                                                  # ranging from 0 (ridge) to 1 (lasso)
  n_folds = 5,
  standardize = FALSE)

## Résultats -----
print(cvfit)
plot(cvfit)
coef(cvfit)


# Extraire les coefficients au lambda optimal
coef_min <- coef(cvfit, lambda = "lambda.min")

# Transformer en vecteur nommé
beta_vec <- as.numeric(coef_min$betas)
names(beta_vec) <- colnames(X) 

# Créer un tibble propre
coef_tbl <- tibble(
  variable = names(beta_vec),
  beta = beta_vec)

# Filtrer les variables sélectionnées (non nulles)
selected_vars <- 
  coef_tbl |>
  filter(beta != 0) |>
  arrange(desc(abs(beta))) |>
  mutate(
    HR = exp(beta))

selected_vars

rm(selected_vars, coef_min, cvfit, coef_tbl, beta_vec, glm_family, penalty_vec, group_starts, group_sizes, ncol_expo_each, X, X_cov, X_expo, dat, covariates_danish, bdd_cases_danish_bis)

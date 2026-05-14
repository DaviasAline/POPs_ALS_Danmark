# Aline Davias
# May 1, 2026 
# Super learner


# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# Just proteins + covariates -----
## Data preparation ----
listWrappers()     # Check available algorithms

data_prep <- bdd_danish |>
  select(als, bmi, smoking_2cat_i, match, all_of(proteomic)) |>
  mutate(smoking_2cat_i = as.numeric(smoking_2cat_i) - 1) |>
  mutate(across(-c(als, match), ~ as.numeric(scale(.x))))

Y_sl <- data_prep$als
X_sl <- data_prep |> select(-als, -match) 
match_ids <- data_prep$match


## Test 1 a : method AUC, many algorithms, no screening, default hyperparameters ----

list_lib_1_a <- c(
  # 1. Linear & Penalized Models
  "SL.leekasso",    # Performs internal screening; expected to be robust
  "SL.glmnet",      # Uses L1/L2 regularization to handle high-dimensional data
  "SL.bayesglm",    # Bayesian approach; may struggle or fail if p > n
  "SL.stepAIC",     # Likelihood-based selection; likely to fail without screening (p > n)
  
  # 2. Tree-Based Models
  "SL.ranger",      # Random Forest; handles many predictors via random feature subsets
  "SL.xgboost",     # Gradient Boosting; capable of handling high dimensionality
  
  # 3. Non-Linear Models
  "SL.ksvm",        # Support Vector Machines; high risk of overfitting without screening
  "SL.gam",         # Spline-based; expected to fail or become unstable with 276 predictors
  
  # 4. Adaptive Splines & Thresholds
  "SL.polymars",    # Logic-based splines; may crash if dimensionality is too high
  "SL.earth",       # MARS algorithm; built-in selection but prone to noise when p is large
  
  # 5. Baseline
  "SL.mean")         # Baseline model; predicts the prevalence (AUC should be 0.50)


set.seed(1996)
sl_fit_CV_1_a <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_a,                                                    # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_a)
plot(sl_fit_CV_1_a) + theme_minimal() 

rm(list_lib_1_a)


## Test des hyper parameters specifique pour chaque algorythm----
### Super learner special SL.ranger ----
screen.top20 <- function(..., ntokeep = 20) {
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

ranger_grid <- expand.grid(
  mtry = c(5, 15, 30),                # Exploration de la dimensionnalité
  min.node.size = c(1, 10, 20),       # Contrôle du lissage/overfitting
  splitrule = c("gini", "extratrees")) # Gini standard vs Variantes robustes


ranger_learners <- unlist(apply(ranger_grid, MARGIN = 1, function(params) {
  
  m_val <- as.numeric(params["mtry"])
  s_val <- as.numeric(params["min.node.size"])
  r_val <- as.character(params["splitrule"])
  
  base_name <- paste0("SL.ranger.m", m_val, ".s", s_val, ".", r_val)
  
  fn <- function(Y, X, newX, family, ...) {
    n_features <- ncol(X)
    mtry_adj <- min(m_val, n_features) 
    SL.ranger(Y = Y, X = X, newX = newX, family = family, 
              mtry = mtry_adj, min.node.size = s_val, 
              splitrule = r_val, num.trees = 1000, ...)
  }
  assign(base_name, fn, envir = .GlobalEnv)
  
  return(list(
    base_name,               
    c(base_name, "screen.top20")))
}), recursive = FALSE)


set.seed(1996)
sl_fit_CV_ranger_1 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = ranger_learners,                                                 # list of hyper parameter combinaisons to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_ranger_1)
plot(sl_fit_CV_ranger_1) + theme_minimal() 

rm(screen.top20, ranger_grid, ranger_learners)

### Super learner special SL.xgboost ----
xgb_grid <- expand.grid(
  max_depth = c(1, 2, 4),      # Arbres très simples à modérés
  eta = c(0.001, 0.01, 0.1),   # Vitesse d'apprentissage
  nrounds = c(500, 1000))       # Nombre d'itérations


xgboost_learners <- unlist(apply(xgb_grid, MARGIN = 1, function(params) {
  
  d_val <- as.numeric(params["max_depth"])
  e_val <- as.numeric(params["eta"])
  n_val <- as.numeric(params["nrounds"])
  
  base_name <- paste0("SL.xgb.d", d_val, ".e", e_val, ".n", n_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.xgboost(Y = Y, X = X, newX = newX, family = family,
               max_depth = d_val, eta = e_val, nrounds = n_val,
               verbose = 0, ...)
  }
  assign(base_name, fn, envir = .GlobalEnv)
  
  return(list(
    base_name,                
    c(base_name, "screen.top20") 
  ))
}), recursive = FALSE)

set.seed(1996)
sl_fit_CV_xgboost_1 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = xgboost_lerners,                                                 # list of hyper parameter combinaisons to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_xgboost_1)
plot(sl_fit_CV_xgboost_1) + theme_minimal() 

rm(xgb_grid, xgboost_learners)

### Super learner special glmnet ----
glmnet_grid <- expand.grid(
  alpha = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)) # Du Ridge pur au Lasso pur


glmnet_learners <- apply(glmnet_grid, MARGIN = 1, function(params) {
  
  a_val <- as.numeric(params["alpha"])
  name <- paste0("SL.glmnet.a", a_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.glmnet(Y = Y, X = X, newX = newX, family = family, 
              alpha = a_val, ...)
  }
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

set.seed(1996)
sl_fit_CV_glmnet_1 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = glmnet_learners, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_glmnet_1)
plot(sl_fit_CV_glmnet_1) + theme_minimal() 

rm(glmnet_grid, glmnet_learners)


### Super Learner special SVM (ksvm) ----
screen.top20 <- function(..., ntokeep = 20) {      
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

ksvm_grid <- expand.grid(
  kernel = c("rbfdot", "polydot", "vanilladot"), # Radial, Polynomial, Linéaire
  C = c(0.1, 1, 10)                              # Coût (régularisation)
)

ksvm_learners <- apply(ksvm_grid, MARGIN = 1, function(params) {
  k_val <- as.character(params["kernel"])
  c_val <- as.numeric(params["C"])
  name <- paste0("SL.ksvm.k.", k_val, ".c", c_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.ksvm(Y = Y, X = X, newX = newX, family = family, 
            kernel = k_val, C = c_val, ...)
  }
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

ksvm_learners_screened <- apply(ksvm_grid, MARGIN = 1, function(params) {
  k_val <- as.character(params["kernel"])
  c_val <- as.numeric(params["C"])
  name <- paste0("SL.ksvm.k.", k_val, ".c", c_val, ".screen20")
  
  fn <- function(Y, X, newX, family, ...) {
    SL.ksvm(Y = Y, X = X, newX = newX, family = family, 
            kernel = k_val, C = c_val, ...)
  }
  assign(name, fn, envir = .GlobalEnv)
  
  return(c(name, "screen.top20"))
})

ksvm_learners_list <- lapply(ksvm_learners, function(x) x)
ksvm_learners_screened_list <- as.list(as.data.frame(ksvm_learners_screened))
ksvm_library_combined <- c(ksvm_learners_list, ksvm_learners_screened_list)

set.seed(1996)
sl_fit_CV_ksvm_1 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = ksvm_library_combined, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_ksvm_1)
plot(sl_fit_CV_ksvm_1) + theme_minimal() 

rm(screen.top20, ksvm_grid, ksvm_learners, ksvm_learners_list, ksvm_learners_screened_list, ksvm_library_combined)


### Super Learner special SL.earth ----
screen.top20 <- function(..., ntokeep = 20) {      
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

earth_grid <- expand.grid(
  degree = c(1, 2),          # 1 = additif, 2 = interactions
  nprune = c(5, 10, 20, 50), # limitation de la complexité
  stringsAsFactors = FALSE)

earth_names <- apply(earth_grid, MARGIN = 1, function(params) {
  d_val <- as.numeric(params["degree"])
  n_val <- as.numeric(params["nprune"])
  
  name <- paste0("SL.earth.d", d_val, ".n", n_val)
  
  fn <- eval(substitute(
    function(Y, X, newX, family, ...) {
      SL.earth(Y = Y, X = X, newX = newX, family = family, 
               degree = D, nprune = N, ...)
    }, list(D = d_val, N = n_val)
  ))
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

earth_library_screened <- lapply(earth_names, function(algo) {
  c(algo, "screen.top20")
})

earth_library <- c(as.list(earth_names), earth_library_screened)

set.seed(1996)
sl_fit_CV_earth_1 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = earth_library, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_earth_1)
plot(sl_fit_CV_earth_1) + theme_minimal() 

rm(screen.top20, earth_grid, earth_names, earth_library_screened, earth_library)


## Test 1 b : method AUC, many algorithms, screening when necessary, optimized hyperparameters ----
### top 20 ----
screen.top20 <- function(..., ntokeep = 20) {      
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

# GLMNet optimal (Alpha 0.9)
SL.glmnet.opt <- function(..., alpha = 0.9) { SL.glmnet(..., alpha = alpha) }

# KSVM optimal (Vanilladot C=10)
SL.ksvm.opt <- function(..., kernel = "vanilladot", C = 10) { 
  SL.ksvm(..., kernel = kernel, C = C) 
}

# Ranger optimal (Extratrees, mtry 15, node size 1)
SL.ranger.opt <- function(Y, X, newX, family, mtry = 15, min.node.size = 1, splitrule = "extratrees", ...) {
  mtry_safe <- min(mtry, ncol(X))
  
  SL.ranger(Y = Y, X = X, newX = newX, family = family, 
            mtry = mtry_safe, 
            min.node.size = min.node.size, 
            splitrule = splitrule, 
            num.trees = 1000, ...)
}

# XGBoost optimal (Prof 2, Learning rate lent)
SL.xgb.opt <- function(..., max_depth = 2, eta = 0.001, nrounds = 500) {
  SL.xgboost(..., max_depth = max_depth, eta = eta, nrounds = nrounds)
}

# Earth optimal (Degree 2 pour les interactions)
SL.earth.opt <- function(..., degree = 2, nfold = 5) {
  SL.earth(..., degree = degree, nfold = nfold)
}


list_lib_1_b_optimized <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top20"),
  c("SL.stepAIC", "screen.top20"),
  
  # 2. Support Vector Machines
  c("SL.ksvm.opt", "screen.top20"),        # Vanilladot C=10 (Robuste à 0.58)
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top20"),      # Extratrees (Le meilleur des arbres)
  c("SL.xgb.opt", "screen.top20"),
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top20"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top20"),
  c("SL.polymars", "screen.top20"),
  
  # 5. Baseline
  "SL.mean")


set.seed(1996)
sl_fit_CV_1_b_opt_20 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_b_optimized,                                          # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_b_opt_20)
plot(sl_fit_CV_1_b_opt_20) + theme_minimal() 

rm(SL.glmnet.opt, SL.glmnet.opt, SL.ksvm.opt, SL.ranger.opt, SL.xgb.opt)


### top 10 ----
screen.top10 <- function(..., ntokeep = 10) {       # on selectionnera les 10 prot les + significatives par t tetst en CV par SL 
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}


list_lib_1_b_optimized_10 <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top10"),
  c("SL.stepAIC", "screen.top10"),
  
  # 2. Support Vector Machines
  c("SL.ksvm.opt", "screen.top10"),        # Vanilladot C=10 (Robuste à 0.58)
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top10"),      # Extratrees (Le meilleur des arbres)
  c("SL.xgb.opt", "screen.top10"),
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top10"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top10"),
  c("SL.polymars", "screen.top10"),
  
  # 5. Baseline
  "SL.mean")


set.seed(1996)
sl_fit_CV_1_b_opt_10 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_b_optimized_10,                                       # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_b_opt_10)
plot(sl_fit_CV_1_b_opt_10) + theme_minimal() 


## Test 1 c : method AUC, selected algorithms, screening when necessary, optimized hyperparameters ----
### top 20 ----
list_lib_1_c_optimized_reduit <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top20"),
  c("SL.stepAIC", "screen.top20"),
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top20"),      # Extratrees (Le meilleur des arbres)
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top20"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top20"),
  
  # 5. Baseline
  "SL.mean")

sl_fit_CV_1_c_opt_20 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_c_optimized_reduit,                                   # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_c_opt_20)
plot(sl_fit_CV_1_c_opt_20) + theme_minimal()

### top 10 ----
list_lib_1_c_optimized_reduit_10 <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top10"),
  c("SL.stepAIC", "screen.top10"),
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top10"),      # Extratrees (Le meilleur des arbres)
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top10"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top10"),
  
  # 5. Baseline
  "SL.mean")

sl_fit_CV_1_c_opt_10 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_c_optimized_reduit_10,                                # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_c_opt_10)
plot(sl_fit_CV_1_c_opt_10) + theme_minimal()

rm(data_prep, X_sl, Y_sl, match_ids)

# Proteins + covariates + time from baseline to diagnosis -----
listWrappers()     # Check available algorithms

data_prep <- bdd_danish |>
  select(als, bmi, smoking_2cat_i, follow_up, match, all_of(proteomic)) |>
  mutate(smoking_2cat_i = as.numeric(smoking_2cat_i) - 1) |>
  group_by(match) |>
  mutate(follow_up = max(follow_up, na.rm = TRUE)) |>
  ungroup() |>
  mutate(across(-c(als, match), ~ as.numeric(scale(.x))))

Y_sl <- data_prep$als
X_sl <- data_prep |> select(-als, -match) 
match_ids <- data_prep$match


## Test 2 a : method AUC, many algorithms, no screening, default hyperparameters ----

list_lib_2_a <- c(
  # 1. Linear & Penalized Models
  "SL.leekasso",     # Performs internal screening; expected to be robust
  "SL.glmnet",      # Uses L1/L2 regularization to handle high-dimensional data
  "SL.bayesglm",    # Bayesian approach; may struggle or fail if p > n
  "SL.stepAIC",     # Likelihood-based selection; likely to fail without screening (p > n)
  
  # 2. Tree-Based Models
  "SL.ranger",      # Random Forest; handles many predictors via random feature subsets
  "SL.xgboost",     # Gradient Boosting; capable of handling high dimensionality
  
  # 3. Non-Linear Models
  "SL.ksvm",        # Support Vector Machines; high risk of overfitting without screening
  "SL.gam",         # Spline-based; expected to fail or become unstable with 276 predictors
  
  # 4. Adaptive Splines & Thresholds
  "SL.polymars",     # Logic-based splines; may crash if dimensionality is too high
  "SL.earth",       # MARS algorithm; built-in selection but prone to noise when p is large
  
  # 5. Baseline
  "SL.mean")         # Baseline model; predicts the prevalence (AUC should be 0.50)


set.seed(1996)
sl_fit_CV_2_a <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_2_a,                                                    # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_2_a)
plot(sl_fit_CV_2_a) + theme_minimal() 


## Test des hyper parameters specifique pour chaque algorythm----
### Super learner special SL.ranger ----
screen.top20 <- function(..., ntokeep = 20) {
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

ranger_grid <- expand.grid(
  mtry = c(5, 15, 30),                # Exploration de la dimensionnalité
  min.node.size = c(1, 10, 20),       # Contrôle du lissage/overfitting
  splitrule = c("gini", "extratrees") # Gini standard vs Variantes robustes
)

ranger_learners <- unlist(apply(ranger_grid, MARGIN = 1, function(params) {
  
  m_val <- as.numeric(params["mtry"])
  s_val <- as.numeric(params["min.node.size"])
  r_val <- as.character(params["splitrule"])
  
  base_name <- paste0("SL.ranger.m", m_val, ".s", s_val, ".", r_val)
  
  fn <- function(Y, X, newX, family, ...) {
    n_features <- ncol(X)
    mtry_adj <- min(m_val, n_features) 
    SL.ranger(Y = Y, X = X, newX = newX, family = family, 
              mtry = mtry_adj, min.node.size = s_val, 
              splitrule = r_val, num.trees = 1000, ...)
  }
  assign(base_name, fn, envir = .GlobalEnv)
  
  return(list(
    base_name,               
    c(base_name, "screen.top20") 
  ))
}), recursive = FALSE)


set.seed(1996)
sl_fit_CV_ranger_2 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = ranger_learners,                                                 # list of hyper parameter combinaisons to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_ranger_2)
plot(sl_fit_CV_ranger_2) + theme_minimal() 



### Super learner special SL.xgboost ----
xgb_grid <- expand.grid(
  max_depth = c(1, 2, 4),      # Arbres très simples à modérés
  eta = c(0.001, 0.01, 0.1),   # Vitesse d'apprentissage
  nrounds = c(500, 1000)       # Nombre d'itérations
)

xgboost_learners <- unlist(apply(xgb_grid, MARGIN = 1, function(params) {
  
  d_val <- as.numeric(params["max_depth"])
  e_val <- as.numeric(params["eta"])
  n_val <- as.numeric(params["nrounds"])
  
  base_name <- paste0("SL.xgb.d", d_val, ".e", e_val, ".n", n_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.xgboost(Y = Y, X = X, newX = newX, family = family,
               max_depth = d_val, eta = e_val, nrounds = n_val,
               verbose = 0, ...)
  }
  assign(base_name, fn, envir = .GlobalEnv)
  
  return(list(
    base_name,                
    c(base_name, "screen.top20") 
  ))
}), recursive = FALSE)

set.seed(1996)
sl_fit_CV_xgboost_2 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = xgboost_learners,                                                # list of hyper parameter combinaisons to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_xgboost_2)
plot(sl_fit_CV_xgboost_2) + theme_minimal() 

### Super learner special glmnet ----
glmnet_grid <- expand.grid(
  alpha = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # Du Ridge pur au Lasso pur
)

glmnet_learners <- apply(glmnet_grid, MARGIN = 1, function(params) {
  
  a_val <- as.numeric(params["alpha"])
  name <- paste0("SL.glmnet.a", a_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.glmnet(Y = Y, X = X, newX = newX, family = family, 
              alpha = a_val, ...)
  }
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

set.seed(1996)
sl_fit_CV_glmnet_2 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = glmnet_learners, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_glmnet_2)
plot(sl_fit_CV_glmnet_2) + theme_minimal() 


### Super Learner special SVM (ksvm) ----
screen.top20 <- function(..., ntokeep = 20) {      
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

ksvm_grid <- expand.grid(
  kernel = c("rbfdot", "polydot", "vanilladot"), # Radial, Polynomial, Linéaire
  C = c(0.1, 1, 10)                              # Coût (régularisation)
)

ksvm_learners <- apply(ksvm_grid, MARGIN = 1, function(params) {
  k_val <- as.character(params["kernel"])
  c_val <- as.numeric(params["C"])
  name <- paste0("SL.ksvm.k.", k_val, ".c", c_val)
  
  fn <- function(Y, X, newX, family, ...) {
    SL.ksvm(Y = Y, X = X, newX = newX, family = family, 
            kernel = k_val, C = c_val, ...)
  }
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

ksvm_learners_screened <- apply(ksvm_grid, MARGIN = 1, function(params) {
  k_val <- as.character(params["kernel"])
  c_val <- as.numeric(params["C"])
  name <- paste0("SL.ksvm.k.", k_val, ".c", c_val, ".screen20")
  
  fn <- function(Y, X, newX, family, ...) {
    SL.ksvm(Y = Y, X = X, newX = newX, family = family, 
            kernel = k_val, C = c_val, ...)
  }
  assign(name, fn, envir = .GlobalEnv)
  
  return(c(name, "screen.top20"))
})

ksvm_learners_list <- lapply(ksvm_learners, function(x) x)
ksvm_learners_screened_list <- as.list(as.data.frame(ksvm_learners_screened))
ksvm_library_combined <- c(ksvm_learners_list, ksvm_learners_screened_list)

set.seed(1996)
sl_fit_CV_ksvm_2 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = ksvm_library_combined, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_ksvm_2)
plot(sl_fit_CV_ksvm_2) + theme_minimal()


### Super Learner special SL.earth ----
earth_grid <- expand.grid(
  degree = c(1, 2),          # 1 = additif, 2 = interactions
  nprune = c(5, 10, 20, 50), # limitation de la complexité
  stringsAsFactors = FALSE)

earth_names <- apply(earth_grid, MARGIN = 1, function(params) {
  d_val <- as.numeric(params["degree"])
  n_val <- as.numeric(params["nprune"])
  
  name <- paste0("SL.earth.d", d_val, ".n", n_val)
  
  fn <- eval(substitute(
    function(Y, X, newX, family, ...) {
      SL.earth(Y = Y, X = X, newX = newX, family = family, 
               degree = D, nprune = N, ...)
    }, list(D = d_val, N = n_val)
  ))
  
  assign(name, fn, envir = .GlobalEnv)
  return(name)
})

earth_library_screened <- lapply(earth_names, function(algo) {
  c(algo, "screen.top20")
})

earth_library <- c(as.list(earth_names), earth_library_screened)

set.seed(1996)
sl_fit_CV_earth_2 <- CV.SuperLearner(
  Y = Y_sl, 
  X = X_sl, 
  family = binomial(), 
  SL.library = earth_library, 
  id = match_ids,
  method = "method.AUC", 
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_earth_2)
plot(sl_fit_CV_earth_2) + theme_minimal()

## Test 2 b : method AUC, many algorithms, screening when necessary, optimized hyperparameters ----
### top 20 ----
screen.top20 <- function(..., ntokeep = 20) {      
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

# GLMNet optimal (Alpha 0.9)
SL.glmnet.opt <- function(..., alpha = 0.9) { SL.glmnet(..., alpha = alpha) }

# KSVM optimal (Vanilladot C=10)
SL.ksvm.opt <- function(..., kernel = "vanilladot", C = 10) { 
  SL.ksvm(..., kernel = kernel, C = C) 
}

# Ranger optimal (Extratrees, mtry 15, node size 1)
SL.ranger.opt <- function(Y, X, newX, family, mtry = 15, min.node.size = 1, splitrule = "extratrees", ...) {
  mtry_safe <- min(mtry, ncol(X))
  
  SL.ranger(Y = Y, X = X, newX = newX, family = family, 
            mtry = mtry_safe, 
            min.node.size = min.node.size, 
            splitrule = splitrule, 
            num.trees = 1000, ...)
}

# XGBoost optimal (Prof 2, Learning rate lent)
SL.xgb.opt <- function(..., max_depth = 2, eta = 0.001, nrounds = 500) {
  SL.xgboost(..., max_depth = max_depth, eta = eta, nrounds = nrounds)
}

# Earth optimal (Degree 2 pour les interactions)
SL.earth.opt <- function(..., degree = 2, nfold = 5) {
  SL.earth(..., degree = degree, nfold = nfold)
}


list_lib_2_b_optimized <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top20"),
  c("SL.stepAIC", "screen.top20"),
  
  # 2. Support Vector Machines
  c("SL.ksvm.opt", "screen.top20"),        # Vanilladot C=10 (Robuste à 0.58)
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top20"),      # Extratrees (Le meilleur des arbres)
  c("SL.xgb.opt", "screen.top20"),
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top20"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top20"),
  c("SL.polymars", "screen.top20"),
  
  # 5. Baseline
  "SL.mean")


set.seed(1996)
sl_fit_CV_2_b_opt_20 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_b_optimized,                                          # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_2_b_opt_20)
plot(sl_fit_CV_2_b_opt_20) + theme_minimal()


### top 10 ----
screen.top10 <- function(..., ntokeep = 10) {       # on selectionnera les 10 prot les + significatives par t tetst en CV par SL 
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}


list_lib_2_b_optimized_10 <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top10"),
  c("SL.stepAIC", "screen.top10"),
  
  # 2. Support Vector Machines
  c("SL.ksvm.opt", "screen.top10"),        # Vanilladot C=10 (Robuste à 0.58)
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top10"),      # Extratrees (Le meilleur des arbres)
  c("SL.xgb.opt", "screen.top10"),
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top10"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top10"),
  c("SL.polymars", "screen.top10"),
  
  # 5. Baseline
  "SL.mean")


set.seed(1996)
sl_fit_CV_2_b_opt_10 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_b_optimized_10,                                       # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_2_b_opt_10)
plot(sl_fit_CV_2_b_opt_10) + theme_minimal()


## Test 2 c : method AUC, selected algorithms, screening when necessary, optimized hyperparameters ----
### top 20 ----
list_lib_2_c_optimized_reduit <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top20"),
  c("SL.stepAIC", "screen.top20"),
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top20"),      # Extratrees (Le meilleur des arbres)
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top20"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top20"),
  
  # 5. Baseline
  "SL.mean")

sl_fit_CV_2_c_opt_20 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_c_optimized_reduit,                                   # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_2_c_opt_20)
plot(sl_fit_CV_2_c_opt_20) + theme_minimal()

### top 10 ----
list_lib_2_c_optimized_reduit_10 <- list(
  # 1. Modèles Linéaires & Pénalisés
  "SL.leekasso",
  "SL.glmnet.opt",                         # Alpha 0.9 (Le champion à 0.618)
  c("SL.bayesglm", "screen.top10"),
  c("SL.stepAIC", "screen.top10"),
  
  # 3. Tree-Based Ensemble
  c("SL.ranger.opt", "screen.top10"),      # Extratrees (Le meilleur des arbres)
  
  # 4. Adaptive Splines (Version optimisée)
  c("SL.earth.opt", "screen.top10"),       # Degré 2 (0.547 vs 0.520)
  c("SL.gam", "screen.top10"),
  
  # 5. Baseline
  "SL.mean")

sl_fit_CV_2_c_opt_10 <- CV.SuperLearner(
  Y = Y_sl,                                                                     # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_c_optimized_reduit_10,                                # list of algorithms to try
  id = match_ids,                                                               # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_2_c_opt_10)
plot(sl_fit_CV_2_c_opt_10) + theme_minimal()

# assemblage et export ----
results_proteomic_ALS_occurrence_SL <-
  list(
    test_1 = list(sl_fit_CV_1_a = sl_fit_CV_1_a, 
       sl_fit_CV_ranger_1 = sl_fit_CV_ranger_1, 
       sl_fit_CV_xgboost_1 = sl_fit_CV_xgboost_1, 
       sl_fit_CV_glmnet_1 = sl_fit_CV_glmnet_1, 
       sl_fit_CV_ksvm_1 = sl_fit_CV_ksvm_1, 
       sl_fit_CV_earth_1 = sl_fit_CV_earth_1, 
       sl_fit_CV_1_b_opt_20 = sl_fit_CV_1_b_opt_20, 
       sl_fit_CV_1_b_opt_10 = sl_fit_CV_1_b_opt_10, 
       sl_fit_CV_1_c_opt_20 = sl_fit_CV_1_c_opt_20, 
       sl_fit_CV_1_c_opt_10 = sl_fit_CV_1_c_opt_10), 
       
       sl_fit_CV_2_a = sl_fit_CV_2_a, 
       sl_fit_CV_ranger_2 = sl_fit_CV_ranger_2, 
       sl_fit_CV_xgboost_2 = sl_fit_CV_xgboost_2, 
       sl_fit_CV_glmnet_2 = sl_fit_CV_glmnet_2, 
       sl_fit_CV_ksvm_2 = sl_fit_CV_ksvm_2, 
       sl_fit_CV_earth_2 = sl_fit_CV_earth_2, 
       sl_fit_CV_2_b_opt_20 = sl_fit_CV_2_b_opt_20, 
       sl_fit_CV_2_b_opt_10 = sl_fit_CV_2_b_opt_10, 
       sl_fit_CV_2_c_opt_20 = sl_fit_CV_2_c_opt_20, 
       sl_fit_CV_2_c_opt_10 = sl_fit_CV_2_c_opt_10)



saveRDS(results_proteomic_ALS_occurrence_SL, file = "~/Documents/POP_ALS_2025_02_03/2_output/results_proteomic_ALS_occurrence_SL.rds")


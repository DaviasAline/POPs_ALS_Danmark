# Aline Davias
# May 1, 2026 
# Super learner


# Data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")




listWrappers()     # Check available algorithms
X_sl <- bdd_danish |> 
  select(bmi, smoking_2cat_i, all_of(proteomic_sd)) |>                         
  mutate(smoking_2cat_i = as.numeric(smoking_2cat_i) - 1)                       # SuperLearner demande des variables numériques


## test 1 ----
### test 1 a : method AUC, many algorithms, no screening, default hyperparameters ----

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
  "SL.mean"         # Baseline model; predicts the prevalence (AUC should be 0.50)
)

set.seed(1996)
sl_fit_CV_1_a <- CV.SuperLearner(
  Y = bdd_danish$als,                                                           # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_a,                                                    # list of algorithms to try
  id = bdd_danish$match,                                                        # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_a)
plot(sl_fit_CV_1_a) + theme_minimal() 


### test 1 b : method AUC, many algorithms, screening when necessary, default hyperparameters ----
screen.top20 <- function(..., ntokeep = 20) {       # on selectionnera les 20 prot les + significatives par t tetst en CV par SL 
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

list_lib_1_b <- list(
  # 1. Linear & Penalized Models
  "SL.leekasso",                      # Built-in screening; robust for high-dimensional proteomics
  "SL.glmnet",                        # Lasso/Elastic Net; natively handles p > n via regularization
  c("SL.bayesglm", "screen.top20"),   # Bayesian GLM; requires p < n for stability/convergence
  c("SL.stepAIC", "screen.top20"),    # Stepwise selection; strictly requires p < n to function
  
  # 2. Tree-Based Ensemble Methods
  "SL.ranger",                        # Random Forest; handles high-dimensional data via feature sampling
  "SL.xgboost",                       # Gradient Boosting; robust to p > n through internal pruning
  
  # 3. Non-Linear & Spline-Based Models
  c("SL.ksvm", "screen.top20"),       # Support Vector Machines; screening prevents overfitting in high dimensions
  c("SL.gam", "screen.top20"),        # Generalized Additive Models; requires p < n to estimate splines
  
  # 4. Adaptive Splines & Threshold Models
  c("SL.polymars", "screen.top20"),   # Polyclass MARS; screening improves knot selection stability
  c("SL.earth", "screen.top20"),      # Multivariate Adaptive Regression Splines; screening reduces noise
  
  # 5. Baseline Reference
  "SL.mean"                           # Simple average; used as a benchmark for "no predictive power"
)

set.seed(1996)
sl_fit_CV_1_b <- CV.SuperLearner(
  Y = bdd_danish$als,                                                           # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_b,                                                    # list of algorithms to try
  id = bdd_danish$match,                                                        # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_b)
plot(sl_fit_CV_1_a) + theme_minimal() 

### test 1 c : method AUC, many algorithms, screening when necessary, tunned hyperparameters ----
screen.top20 <- function(..., ntokeep = 20) {       # on selectionnera les 20 prot les + significatives par t tetst en CV par SL 
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

if("package:gam" %in% search()) detach("package:gam", unload=TRUE)
library(mgcv) # On garde uniquement mgcv qui est plus moderne

# Elastic Net (Balance between Lasso and Ridge)
SL.glmnet.alpha5 <- function(..., alpha = 0.5) {  #  Elastic Net mixing parameter
  SL.glmnet(..., alpha = alpha)
}

# Random Forest: Reducing complexity to avoid overfitting
SL.ranger.tuned <- function(..., 
                            num.trees = 1000,          # number of trees to be grown in the forest.
                            min.node.size = 10) {      # minimum number of observations required in a terminal node
  SL.ranger(..., num.trees = num.trees, min.node.size = min.node.size)
}

# XGBoost: Slower learning rate and shallower trees
SL.xgb.tuned <- function(..., 
                         nrounds = 500,     # The number of boosting iterations (trees) to be performed.
                         max_depth = 2,     # The maximum depth of a tree
                         eta = 0.01) {      # earning rate (or "shrinkage"). It scales the contribution of each new tree. A very low value (0.01) means the model learns slowly, which is the most effective way to prevent overfitting in boosting.
  SL.xgboost(..., nrounds = nrounds, max_depth = max_depth, eta = eta, verbose = 0)
}

# Earth: Testing interactions (degree = 2)
SL.earth.tuned <- function(..., 
                           degree = 2) {    # maximum degree of interaction.
  SL.earth(..., degree = degree)
}

list_lib_1_c <- list(
  "SL.leekasso",                                 # Internal screening
  "SL.glmnet",                                   # Default Lasso
  c("SL.glmnet.alpha5", "screen.top20"),         # Tuned Elastic Net
  c("SL.bayesglm", "screen.top20"),              # Standard Bayes
  c("SL.stepAIC", "screen.top20"),               # Standard Stepwise
  c("SL.ranger.tuned", "screen.top20"),          # Tuned Random Forest
  c("SL.xgb.tuned", "screen.top20"),             # Tuned XGBoost
  c("SL.ksvm", "screen.top20"),                  # Standard SVM (RBF)
  c("SL.gam", "screen.top20"),                   # Standard GAM
  c("SL.polymars", "screen.top20"),              # Standard Polymars
  c("SL.earth.tuned", "screen.top20"),           # Tuned Earth (Interactions)
  "SL.mean"                                      # Baseline
)

set.seed(1996)
sl_fit_CV_1_c <- CV.SuperLearner(
  Y = bdd_danish$als,                                                           # outcome
  X = X_sl,                                                                     # predictors (proteins + covariates)
  family = binomial(),                                                          # binary outcome (0/1)
  SL.library = list_lib_1_c,                                                    # list of algorithms to try
  id = bdd_danish$match,                                                        # matching
  method = "method.AUC",                                                        # performance method (AUC), default would be Non-Negative Least Squares
  cvControl = list(V = 10),                                                     # external CV: 10 folds to evaluate the overall performance (generalizability)
  innerCvControl = list(list(V = 10)),                                          # internal CV: 10 folds to estimate the optimal weights for each algorithm
  control = list(
    saveFitLibrary = TRUE,                                                      # keeps individual model fits in memory (required for variable importance later)
    trimLogit = 0.001))                                                         # Prevents numerical instability for probabilities near 0 or 1

summary(sl_fit_CV_1_c)
plot(sl_fit_CV_1_c) + theme_minimal() 


### test 1 d : method AUC, selected algorithms, screening when necessary, tunned hyperparameters ----
# Relaxed Lasso: Better coefficient estimation
SL.glmnet.relaxed <- function(..., relax = TRUE) {
  SL.glmnet(..., relax = relax)
}

# Screener Top 50: Giving more "room" to non-linear models
screen.top50 <- function(..., ntokeep = 50) {
  screen.ttest(..., pvp = 1, minstep = ntokeep)
}

# --- 2. Elite Library 1d ---

list_lib_1_d <- list(
  "SL.leekasso",                                  # The current champion (0.614)
  "SL.glmnet.relaxed",                            # Improved version of your 2nd best
  "SL.glmnet",                                    # Standard Lasso (All proteins)
  c("SL.earth.tuned", "screen.top50"),            # Interactions with more variables
  c("SL.gam", "screen.top50"),                    # Flexibility with more variables
  c("SL.bayesglm", "screen.top20"),               # Robustness baseline
  "SL.mean"                                       # Benchmarking
)

# --- 3. Execution ---

set.seed(1996)
sl_fit_CV_1_d <- CV.SuperLearner(
  Y = bdd_danish$als,
  X = X_sl,
  family = binomial(),
  SL.library = list_lib_1_d,
  id = bdd_danish$match,
  method = "method.AUC",
  cvControl = list(V = 10),
  innerCvControl = list(list(V = 10)),
  control = list(saveFitLibrary = TRUE, trimLogit = 0.001))

summary(sl_fit_CV_1_d)

### test 1 e ----
#### nested cross validation ----
list_lib_1_e <- list(
  "SL.leekasso",                         
  "SL.glmnet",                           
  c("SL.earth.tuned", "screen.top20"),  
  c("SL.gam", "screen.top20"),           
  "SL.mean")

set.seed(1996)
sl_fit_CV_1_e <- CV.SuperLearner(
  Y = bdd_danish$als,
  X = X_sl,
  family = binomial(),
  SL.library = list_lib_1_e,
  id = bdd_danish$match,
  method = "method.AUC",
  cvControl = list(V = 10),
  innerCvControl = list(list(V = 10)),
  control = list(
    saveFitLibrary = TRUE, 
    trimLogit = 0.001))

#### weights of the algorithms ----
set.seed(100)
sl_fit_1_e <- SuperLearner(
  Y = bdd_danish$als, 
  X = X_sl, 
  family = binomial(), 
  SL.library = list_lib_1_e, 
  method = "method.AUC", 
  id = bdd_danish$match)

weights_1_e <- data.frame(
  Algorithm = names(sl_fit_1_e$coef),
  Weight = sl_fit_1_e$coef)

weights_1_e[order(-weights$Weight), ]

#### extraction of proteins selected by lasso (SL.glmnet_All) ----
lasso_model <- sl_fit_1_e$fitLibrary$SL.glmnet_All$object
lasso_coefs <- as.matrix(coef(lasso_model, s = "lambda.min"))
selected_proteins_lasso_1_e <- data.frame(
  Protein = rownames(lasso_coefs),
  Coefficient = as.numeric(lasso_coefs)) |> 
  filter(Coefficient != 0 & Protein != "(Intercept)") |>
  arrange(desc(abs(Coefficient)))

rm(lasso_model, lasso_coefs)

# assemblage et export ----
results_proteomic_ALS_occurrence_SL <-
  list(sl_fit_CV_1_a = sl_fit_CV_1_a, 
       sl_fit_CV_1_b = sl_fit_CV_1_b, 
       sl_fit_CV_1_c =sl_fit_CV_1_c, 
       sl_fit_CV_1_d = sl_fit_CV_1_d, 
       sl_fit_CV_1_e = sl_fit_CV_1_e, 
       weights_1_e = weights_1_e, 
       selected_proteins_lasso_1_e = selected_proteins_lasso_1_e)

saveRDS(results_proteomic_ALS_occurrence_SL, file = "~/Documents/POP_ALS_2025_02_03/2_output/results_proteomic_ALS_occurrence_SL.rds")


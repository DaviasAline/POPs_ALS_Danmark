# Work on machine learning 04/12/2025 ----

## supervised machine learning : random forest method ----
### data preparation ----
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



### unconditional random forest -----
#### chargement des packages ----
library(randomForest)
library(caret)

#### data preparation ----
set.seed(1996)
matches <- unique(bdd_rf$match)                                                 # 70% for training, 30% for testing
train_matches <- sample(matches, size = floor(0.7 * length(matches)))

train_data <- bdd_rf |> filter(match %in% train_matches) |> select(-match)      # we remove match because we are already adjusting on baseline age and sex
test_data  <- bdd_rf |> filter(!match %in% train_matches) |> select(-match)     # we remove match because we are already adjusting on baseline age and sex

rm(matches, train_matches)

#### random forest with randomForest() ----
set.seed(1996)

rf_model <- randomForest(
  als ~ .,
  data = train_data,
  ntree = 500,           # nombre d’arbres
  mtry = floor(sqrt(ncol(train_data)-1)), # nombre de variables testées à chaque split
  importance = TRUE)

print(rf_model)

#### importance des variables ----
importance(rf_model)
varImpPlot(rf_model)

#### performance predictive ----
pred <- predict(rf_model, newdata = test_data, type = "response")
confusionMatrix(pred, test_data$als)



### conditionnal random forest ----
#### chargement des packages ----
library(party)
library(pROC)
library(groupdata2)   # pour CV / split groupé

#### data preparation ----
set.seed(1996)
matches <- unique(bdd_rf$match)                                                 # 70% for training, 30% for testing
train_matches <- sample(matches, size = floor(0.7 * length(matches)))

train_data <- 
  bdd_rf |> filter(match %in% train_matches) |> select(-baseline_age, -sex)     # we remove baseline age and sex because we are already matching on it with variable match
test_data  <- 
  bdd_rf |> filter(!match %in% train_matches) |> select(-baseline_age, -sex)    # we remove baseline age and sex because we are already matching on it with variable match

rm(matches, train_matches)


#### conditional random forest with cforest() ----
set.seed(1996)

cforest_ctrl <- cforest_unbiased(
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data) - 2)))                                     # nb de variables testées à chaque split

rf_model <- cforest(
  als ~ .,                                                                      # car train_data n’a que les variables utiles
  data = train_data,
  controls = cforest_ctrl)

#### importance des variables ----
varimp_cond <- varimp(rf_model, conditional = TRUE)

importance_sorted <- sort(varimp_cond, decreasing = TRUE)                       # Classement décroissant
importance_sorted

#### performance predictive ----
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



### conditional random forest avec 10 fold cross validation ----
#### chargement des packages ----
library(party)
library(pROC)
library(groupdata2)   # group K-fold CV

#### data preparation ----
bdd_rf <- 
  bdd_danish |>  
  select(als, match, baseline_age, sex, smoking_2cat_i, bmi, 
         all_of(proteomic_selected)) |>
  mutate(
    als   = factor(as.character(als)),        
    match = factor(as.character(match)))





#### 10-fold group cross-validation (group = match) ----
set.seed(1996)
bdd_cv <- fold(                                     # groupdata2 construit une colonne .folds qui respecte les groupes
  data   = bdd_rf |> select(-baseline_age, -sex),   # on retire ces deux var car on match avec var match
  k      = 10,
  id_col = "match",    # groupement par match
  method = "n_dist")


#### Boucle CV sur les 10 folds ----
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


#### Results ----
mean_auc <- mean(auc_values)
sd_auc   <- sd(auc_values)

mean_acc <- mean(accuracy_values)
sd_acc   <- sd(accuracy_values)

cat("\n=== Cross-Validated Performance ===\n")
cat(sprintf("AUC (mean ± SD): %.3f ± %.3f\n", mean_auc, sd_auc))
cat(sprintf("Accuracy (mean ± SD): %.3f ± %.3f\n", mean_acc, sd_acc))

#### Importance des variables sur modèle final (entraînement complet) -----
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


## unsupervised machine learning : K-mean clustering method ----






# Prise en compte differences de proteins EPIC vs UK biobank ----
library(readxl)
library(tidyverse)
library(questionr)
library(gtsummary)
Olink_explore <- read_excel("~/Documents/Appels à projets/9. UK biobank access 2026/Olink Explore protein list.xltx")
Olink_target <- read_excel("~/Documents/Appels à projets/9. UK biobank access 2026/Olink Target 96 protein list.xltx")


test <- left_join(Olink_target, Olink_explore, by = c("UniProt ID", "Gene", "Protein name"))


# investigation des autres maladies neurologiques -----

library(haven)
lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long <- 
  read_dta("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/lpr_DKEPIC_ALS_sum_3k(johnni@cancer.dk).dta")

lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_ <- 
  lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |> 
  rename(sample = code, serial_no = loebenr, match = saet) |>
  group_by(sample, serial_no, match) |>
  mutate(id_temp = row_number()) |>
  ungroup() |>
  pivot_wider(
    names_from = id_temp, 
    values_from = -c(sample, serial_no, match, id_temp),
    names_sep = "_",
    names_vary = "slowest")

lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long <-   
  lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |> 
  rename(sample = code, serial_no = loebenr, match = saet)  |> 
  mutate(c_diag = as.factor(c_diag))

als_data <- 
  bdd_danish |> 
  select(sample, serial_no, match, als, proteomic_neuro_explo_NEFL) |>
  mutate(sample = as.numeric(sample))

lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long <- 
  lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |> 
  left_join(als_data, by = c("sample", "serial_no", "match")) 
rm(als_data)

lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long <- 
  lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |>
  mutate(
    c_diag_name = 
      case_when(
        str_detect(c_diag, "G35|34009") ~ "multiple sclerosis (MS)", 
        str_detect(c_diag, "G30|29010") ~ "Alzheimer's disease", 
        str_detect(c_diag, "G20|33200") ~ "Parkinson's disease", 
        str_detect(c_diag, "G310|F020|29011") ~ "frontotemporal dementia (FTD)", 
        str_detect(c_diag, "G318|F028|29019") ~ "dementia with Lewy bodies (DLB)", 
        str_detect(c_diag, "E752|27280") ~ "Niemann-Pick disease type C1 (NPC1)",
        str_detect(c_diag, "G10|33100") ~ "Huntington disease (HD)", 
        str_detect(c_diag, "G231|33310") ~ "progressive supranuclear palsy (PSP)", 
        str_detect(c_diag, "G3185|33199") ~ "corticobasal degeneration (CBD)", 
        str_detect(c_diag, "G232|33311") ~ "multiple system atrophy (MSA)", 
        str_detect(c_diag, "A810|06500") ~ "prion diseases (CJD)",
        TRUE ~ NA_character_)) |>
  select(sample, serial_no, als, match, c_diag, c_diag_name, everything()) |>
  arrange(sample, match, serial_no, als)


lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |>
  filter(
    str_starts(c_diag, "DG35") | str_starts(c_diag, "34009") |                  # ICD-10 and ICD-8 for multiple sclerosis (MS)
      str_starts(c_diag, "DG30") | str_starts(c_diag, "29010") |                  # ICD-10 and ICD-8 for Alzheimer's disease
      str_starts(c_diag, "DG20") | str_starts(c_diag, "33200") |                  # ICD-10 and ICD-8 for Parkinson's disease
      str_starts(c_diag, "DG310") | str_starts(c_diag, "DF020") | str_starts(c_diag, "29011") | # ICD-10 and ICD-8 for frontotemporal dementia (FTD)
      str_starts(c_diag, "DG318") | str_starts(c_diag, "DF028") | str_starts(c_diag, "29019") | # ICD-10 and ICD-8 for dementia with Lewy bodies (DLB)
      str_starts(c_diag, "DE752") | str_starts(c_diag, "27280") |                 # ICD-10 and ICD-8 for Niemann-Pick disease type C1 (NPC1)
      str_starts(c_diag, "DG10") | str_starts(c_diag, "33100") |                  # ICD-10 and ICD-8 for Huntington disease (HD)
      str_starts(c_diag, "DG231") | str_starts(c_diag, "33310") |                 # ICD-10 and ICD-8 for progressive supranuclear palsy (PSP)
      str_starts(c_diag, "DG3185") | str_starts(c_diag, "33199") |                # ICD-10 and ICD-8 for corticobasal degeneration (CBD)
      str_starts(c_diag, "DG232") | str_starts(c_diag, "33311") |                 # ICD-10 and ICD-8 for multiple system atrophy (MSA)
      str_starts(c_diag, "DA810") | str_starts(c_diag, "6500")) |>                # ICD-10 and ICD-8 for prion diseases (CJD)
  distinct(sample, c_diag_name, .keep_all = TRUE) |>
  arrange(als) |>
  select(c_diag_name, als) |>
  tbl_summary(
    by = als,
    label = list(c_diag_name ~ "Other neurological diseases"),
    digits = everything() ~ 0 ) |>
  modify_header(
    update = list(
      stat_1 ~ "**Controls**  \n(N = {n})",
      stat_2 ~ "**Cases**  \n(N = {n})")) |>
  modify_spanning_header(all_stat_cols() ~ "**ALS**") |>
  bold_labels()


ident_other_neuro <- lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |>
  #filter(als == 0) |>
  filter(
    str_starts(c_diag, "DG35") | str_starts(c_diag, "34009") |                  # ICD-10 and ICD-8 for multiple sclerosis (MS)
      str_starts(c_diag, "DG30") | str_starts(c_diag, "29010") |                  # ICD-10 and ICD-8 for Alzheimer's disease
      str_starts(c_diag, "DG20") | str_starts(c_diag, "33200") |                  # ICD-10 and ICD-8 for Parkinson's disease
      str_starts(c_diag, "DG310") | str_starts(c_diag, "DF020") | str_starts(c_diag, "29011") | # ICD-10 and ICD-8 for frontotemporal dementia (FTD)
      str_starts(c_diag, "DG318") | str_starts(c_diag, "DF028") | str_starts(c_diag, "29019") | # ICD-10 and ICD-8 for dementia with Lewy bodies (DLB)
      str_starts(c_diag, "DE752") | str_starts(c_diag, "27280") |                 # ICD-10 and ICD-8 for Niemann-Pick disease type C1 (NPC1)
      str_starts(c_diag, "DG10") | str_starts(c_diag, "33100") |                  # ICD-10 and ICD-8 for Huntington disease (HD)
      str_starts(c_diag, "DG231") | str_starts(c_diag, "33310") |                 # ICD-10 and ICD-8 for progressive supranuclear palsy (PSP)
      str_starts(c_diag, "DG3185") | str_starts(c_diag, "33199") |                # ICD-10 and ICD-8 for corticobasal degeneration (CBD)
      str_starts(c_diag, "DG232") | str_starts(c_diag, "33311") |                 # ICD-10 and ICD-8 for multiple system atrophy (MSA)
      str_starts(c_diag, "DA810") | str_starts(c_diag, "6500")) |>                # ICD-10 and ICD-8 for prion diseases (CJD)
  distinct(sample, c_diag_name, .keep_all = TRUE) |>
  pull(sample)

diagnostics_lookup <- lpr_DKEPIC_ALS_sum_3k_johnni_cancer_dk_long |>
  filter(!is.na(c_diag_name)) |>
  distinct(sample, c_diag_name) |>
  mutate(sample = as.character(sample))

bdd_plot <- bdd_danish |>
  filter(match != 159) |>
  left_join(diagnostics_lookup, by = "sample") |>
  mutate(
    highlight = if_else(is.na(c_diag_name), "None", c_diag_name),
    highlight = factor(highlight),
    als_label = factor(als, levels = c(0, 1), labels = c("Controls", "Cases"))) |>
  select(sample, als, highlight, als_label, proteomic_neuro_explo_NEFL)

descrip_num(data = bdd_plot, vars = "proteomic_neuro_explo_NEFL")

library(RColorBrewer)


n_diseases <- length(unique(bdd_plot$highlight)) - 1 

ggplot(bdd_plot, aes(x = als_label, y = proteomic_neuro_explo_NEFL)) +
  geom_violin(fill = "grey95", color = "grey50", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "grey30", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = highlight, alpha = highlight), width = 0.15, size = 2) +
  coord_cartesian(ylim = c(1, 5.1)) + 
  scale_color_manual(values = c("None" = "grey30", 
                                setNames(brewer.pal(min(n_diseases, 9), "Set1"), 
                                         setdiff(unique(bdd_plot$highlight), "None")))) +
  scale_alpha_manual(
    values = c("None" = 0.2, setNames(rep(1, n_diseases), setdiff(unique(bdd_plot$highlight), "None"))), 
    guide = "none") +
  theme_minimal() +
  labs(
    x = "ALS",
    y = "NfL distribution",
    color = "Other neurologic\ndiagnostic")




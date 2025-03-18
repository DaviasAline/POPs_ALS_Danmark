# Aline Davias
# 03/02/2025

# package loading ----
library(haven)
library(tidyverse)
library(gtsummary)
library(readxl)
library(splines)
library(questionr)
library(esquisse)
library(survival)
library(splines)
library(writexl)
library(mice)
library(see)
source("~/Documents/POP_ALS_2025_02_03/1_codes/0_functions.R")

# data loading ----
bdd_questionaire <- read_sas("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/pop_als_1ak.sas7bdat", 
                        NULL) %>%
  mutate(code = as.character(code)) %>%
  select(sample = code, serial_no = loebenr, match = saet,                      # identification variables 
         birth_date = fsdato, baseline_date = mdato,                            # important dates
         als_date = als_ddat, death_date = DEATH_DATE,                          
         sex, marital_status = civil,                                           # important metadata
         high_educ = HQJ_SKO, medium_educ = MEL_SKO, low_educ = LAV_SKO,        # education                          
         alcohol = alko, smoking = rygning,                                     # alcohol and smoking 
         bmi, heart_blood_clot = S25A01N, brain_blood_clot = S25C01N,           # important health information 
         hypertension = S25D01N, hypertension_medicine = S25D03N, diabetes = S25F01N,                      
         cholesterol = blodprqv, blod_dias = blodtdia, blod_sys = blodtsys, 
         diy_work = gqrselv, gardening = havearb, house_work = husarb,          # domestic work
         house_work__hours = thusarb, 
         everything())                                                          # occupations

bdd_POPs <- read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/2024_ALS-Denmark-POP-Final-Results.xlsx", 
                          sheet = "Results-blank-subtracted", 
                          range = "A3:AC505") %>% 
  slice(-1) %>%
  rename(sample = Sample) %>%
  mutate(sample = as.character(sample)) 

bdd_danish <- left_join(bdd_questionaire, bdd_POPs, by = "sample")                     # Delated 40 "T" idents (check with Marc)

bdd_loq <- 
  read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/2024_ALS-Denmark-POP-Final-Results.xlsx", 
                      sheet = "Results-blank-subtracted", 
                      range = "D1:AB3") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("LOQ") %>%
  select(variable = V2, "% > LOQ" = V1, LOQ) %>%
  mutate(
    LOQ = gsub("\\.\\.\\..*", "", LOQ), 
    LOQ =  as.numeric(LOQ), 
    `% > LOQ` = as.numeric(`% > LOQ`), 
    `% > LOQ` = round(`% > LOQ`, 1), 
    variable = gsub(",", "", variable),  
    variable = gsub("'", "", variable),
    variable = gsub("-", "_", variable))  

# data cleaning ----
POPs <- bdd_POPs %>% select(-c(n, Batch, sample, Comment)) %>% colnames()
bdd_danish <- bdd_danish %>%  mutate(across(all_of(POPs), as.numeric)) 
names(bdd_danish) <- ifelse(names(bdd_danish) %in% POPs, gsub("-", "_", names(bdd_danish)), names(bdd_danish))
names(bdd_POPs) <- ifelse(names(bdd_POPs) %in% POPs, gsub("-", "_", names(bdd_POPs)), names(bdd_POPs))

bdd_POPs <- bdd_POPs %>% 
  rename(pp_DDT = "p,p'_DDT", 
         pp_DDE = "p,p'_DDE")

bdd_danish <- bdd_danish %>% 
  rename(pp_DDT = "p,p'_DDT", 
         pp_DDE = "p,p'_DDE") %>%
  mutate(
    smoking = as.character(smoking), 
    smoking = fct_recode(smoking, 
                         "Never" = "1",
                         "Previous" = "2",
                         "Current" = "3"), 
    marital_status  = as.character(marital_status),
    marital_status = fct_recode(marital_status, 
                       "Widowed" = "1",
                       "Divorced" = "2",
                       "Married/cohabit" = "3",
                       "Unmarried" = "4",
                       NULL = "9"), 
    sex = fct_recode(sex, 
                     "Female" = "F",
                     "Male" = "M"))


# vectors creation ----
POPs <- bdd_POPs %>% select(-c(n, Batch, sample, Comment)) %>% colnames()
POPs_quart <- paste0(POPs, "_quart")
POPs_ter <- paste0(POPs, "_ter")

PCBs <- bdd_POPs %>% select(all_of(POPs)) %>% select(starts_with("PCB")) %>% colnames()
PCBs_quart <- paste0(PCBs, "_quart")
PCBs_ter <- paste0(PCBs, "_ter")

PBDEs<- bdd_POPs %>% select(all_of(POPs)) %>% select(starts_with("BDE")) %>% colnames()
PBDEs_quart <- paste0(PBDEs, "_quart")
PBDEs_ter <- paste0(PBDEs, "_ter")

OCPs <- bdd_POPs %>% select("PeCB", "HCB", "α_HCH", "β_HCH", "γ_HCH", 
                            "Oxychlordane", "Transnonachlor", "pp_DDT", "pp_DDE") %>% 
  colnames()
OCPs_quart <- paste0(OCPs, "_quart")
OCPs_ter <- paste0(OCPs, "_ter")

POPs_group <- c("PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "β_HCH", "Σchlordane", "ΣPBDE")
POPs_group_quart <- paste0(POPs_group, "_quart")
POPs_group_ter <- paste0(POPs_group, "_ter")
POPs_group_95 <- paste0(POPs_group, "_95")
POPs_group_outlier <- paste0(POPs_group, "_outlier")

# 0 values vizualization ----
# bdd_danish %>% 
#   select(all_of(POPs)) %>% 
#   filter(if_any(everything(), ~ . == 0)) %>% 
#   select(where(~ any(. == 0))) %>%
#   View()
# 
# bdd_danish %>% select(PeCB) %>% filter(PeCB ==0) %>% View() # 5 valeurs 0
# bdd_danish %>% select(`α-HCH`) %>% filter(`α-HCH` ==0) %>% View()  # 372 valeurs 0
# bdd_danish %>% select(`γ-HCH`) %>% filter(`γ-HCH` ==0) %>% View()  # 385 valeurs 0
# bdd_danish %>% select(`BDE-99`) %>% filter(`BDE-99` ==0) %>% View()  # 1 valeur 0


# variables creation ----
## exposures ----
bdd_danish <- bdd_danish %>%
  mutate(across(all_of(POPs), ~ factor(ntile(.x, 4), labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) %>%
  mutate(across(all_of(POPs), ~ factor(ntile(.x, 3), labels = c("T1", "T2", "T3")),
                .names = "{.col}_ter")) %>% 
  mutate(
    PCB_DL = PCB_118 + PCB_156, 
    PCB_NDL = PCB_28 + PCB_52 + PCB_74 + PCB_99 + PCB_101 + PCB_138 + PCB_153 + 
      PCB_170 + PCB_180 + PCB_183 + PCB_187, 
    PCB_4 = PCB_118 + PCB_138 + PCB_153 + PCB_180, 
    ΣPBDE = BDE_47 + BDE_99 + BDE_153, 
    ΣDDT = pp_DDT + pp_DDE,  
    ΣHCH = β_HCH + γ_HCH, 
    Σchlordane = Transnonachlor + Oxychlordane) %>% # alpha HCH, gamma HCH exluded, beta HCH studied alone
  mutate(
    across(c("PCB_DL", 
             "PCB_NDL", 
             "PCB_4", 
             "ΣPBDE", 
             "ΣDDT", 
             "ΣHCH", 
             "Σchlordane"), ~ factor(ntile(.x, 4), labels = c("Q1", "Q2", "Q3", "Q4")),
           .names = "{.col}_quart")) %>%
  mutate(
    across(c("PCB_DL", 
             "PCB_NDL", 
             "PCB_4", 
             "ΣPBDE", 
             "ΣDDT", 
             "ΣHCH", 
             "Σchlordane"), ~ factor(ntile(.x, 3), labels = c("T1", "T2", "T3")),
           .names = "{.col}_ter")) %>%
  mutate(across(all_of(POPs_group), 
                ~ if_else(. > quantile(., 0.95, na.rm = TRUE), NA_real_, .),
                .names = "{.col}_95")) %>%
  mutate(across(all_of(POPs_group), ~ {
    Q1 <- quantile(., 0.25, na.rm = TRUE)
    Q3 <- quantile(., 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    ifelse(. < lower_bound | . > upper_bound, NA, .)
    }, .names = "{.col}_outlier")) 


replace_with_median <- function(data, var, quartile_var) {                       
  new_var_name <- paste0(rlang::as_name(ensym(var)), "_quart_med")  
  data %>%
    group_by({{quartile_var}}) %>%
    mutate(!!new_var_name := median({{var}}, na.rm = TRUE)) %>%
    ungroup()
}


bdd_danish <- bdd_danish %>%
  replace_with_median(HCB, HCB_quart) %>%
  replace_with_median(PCB_DL, PCB_DL_quart) %>%
  replace_with_median(PCB_NDL, PCB_NDL_quart) %>%
  replace_with_median(PCB_4, PCB_4_quart) %>%
  replace_with_median(ΣPBDE, ΣPBDE_quart) %>%
  replace_with_median(ΣDDT, ΣDDT_quart) %>%
  replace_with_median(ΣHCH, ΣHCH_quart) %>%
  replace_with_median(Σchlordane, Σchlordane_quart) %>%
  replace_with_median(PeCB, PeCB_quart) %>%
  replace_with_median(β_HCH, β_HCH_quart)
rm(replace_with_median)

POPs_group_quart_med <- bdd_danish %>% select(contains("quart_med")) %>% colnames()
descrip_num(data = bdd_danish, vars = POPs_group_quart_med)
bdd_danish %>% select(all_of(POPs_group_quart_med)) %>% tbl_summary()

## metadata ----
### education, age variables, follow up ----
bdd_danish <- bdd_danish %>% 
  mutate(
    als = ifelse(!is.na(als_date), 1, 0), 
    education = case_when(
      high_educ == 1 & medium_educ == 0 & low_educ == 0 ~ "≥7 years of primary school", 
      high_educ == 0 & medium_educ == 1 & low_educ == 0 ~ "7-10 years of primary school", 
      high_educ == 0 & medium_educ == 0 & low_educ == 1 ~ "≤7 years of primary school", 
      .default = NA),
    education = as.factor(education),
    education =  fct_relevel(education, 
                             "≤7 years of primary school", 
                             "7-10 years of primary school",
                             "≥7 years of primary school"),
    baseline_age = as.numeric(difftime(baseline_date, birth_date, units = "days")) / 365.25, 
    diagnosis_age = as.numeric(difftime(als_date, birth_date, units = "days")) / 365.25, 
    death_age = as.numeric(difftime(death_date, birth_date, units = "days")) / 365.25, 
    follow_up = as.numeric(difftime(als_date, baseline_date, units = "days"))/365.25) %>%
  select(
    "sample", "serial_no", "match", 
    "birth_date", "baseline_date", "als_date", "death_date",
    als, baseline_age, diagnosis_age, death_age, follow_up, 
    "sex", "marital_status", 
    "high_educ", "medium_educ", "low_educ", education, 
    everything())

### occupation ----
#bdd_questionaire %>% select(sample, starts_with("S59")) %>% View()
#bdd_questionaire %>% select(starts_with("S59")) %>% tbl_summary()

# First category should be "White-collar workers" (eg technical, medical, scientific, administrative, clerical etc. 
## professional occupations, including nurses, artistic occupations and students)
## I've put "S59X44N" # hair dresser, "S59X47N" # waiter/waitress, "S59X48N" # chef (cook chef?), "S59X49N" # Health care

# Second category should be "farmers and related agricultural workers"
## I've put "S59X01N" # livestock, "S59X02N" # agriculture, "S59X03N" # agriculture - insects / weeds

# Third category should be "Industrial workers" (mainly metal, machine, paper and construction workers)
## I've put "S59X04N" # Mining, stone and limestone quarrying, "S59X05N" # galvanization industry, 
## "S59X06N" # chemical industry, "S59X07N" # chemistry - refinery, "S59X08N" # chemistry - dyeing, 
## "S59X09N" # Chemistry - Chemical laboratory, "S59X1ON" # rubber industry, "S59X11N" # textile industry, 
## "S59X12N" # textile dyeing, "S59X13N" # textile - spinning, weaving, "S59X14N" # leather tanning and dyeing, 
## "S59X15N" # leather products (shoes), "S59X16N" # wood industry, "S59X17N" # wood industry - furniture manufacturing, 
## "S59X18N" # metal processing, "S59X19N" # metal turning etc, "S59X20N" # metal welding, "S59X21N" # metal painting, 
## "S59X22N" # foundry, "S59X23N" # foundry - steel rolling mill, "S59X24N" # foundry - alloy rolling mill, 
## "S59X25N" # shipbuiding industry, "S59X26N" # electronics industry, "S59X27N" # glass industry,
## "S59X28N" # graphic industry, "S59X29N" # construction industry, "S59X30N" # construction - roofer,
## "S59X31N" # construction - asphalt worker, "S59X32N" # construction - demolition worker,
## "S59X36N" # paper of paper product manufacturing, "S59X37N" # asbestos or cement, 
## "S59X38N" # asbestos insulation, "S59X39N" # cement product industry, "S59X40N" # porcelain/faience/earthenware
## "S59X42N" # painter, "S59X43N" # welder (soudeur)

# Third category should be Other blue-collar workers, or occupation unknown (eg cleaners, caretakers, other service workers, housewives)
## I've put "S59X33N" # transportation, "S59X34N" # transportation - truck driver, 
## "S59X35N" # transportation - bus or taxi, "S59X41N" # butcher, 
## "S59X45N" # tank attendant (pompiste), "S59X46N" # automechanic

bdd_occupation <- bdd_danish %>% 
  select(sample, starts_with("S59X")) 

bdd_occupation <- bdd_occupation %>%
  mutate(
    occupation_verif = rowSums(bdd_occupation[2:50] == 1, na.rm = TRUE),
    occupation_verif_white_collar = rowSums(bdd_occupation[c("S59X44N", # hair dresser 
                                                             "S59X47N", # waiter/waitress
                                                             "S59X48N", # chef (cook chef?) 
                                                             "S59X49N"  # Health care
    )] == 1, na.rm = TRUE),
    occupation_verif_white_collar_percent = occupation_verif_white_collar/4*100,
    occupation_verif_agricultural = rowSums(bdd_occupation[c("S59X01N", # livestock
                                                             "S59X02N", # agriculture
                                                             "S59X03N"  # agriculture - insects / weeds
    )] == 1, na.rm = TRUE),
    occupation_verif_agricultural_percent = occupation_verif_agricultural/3*100,
    occupation_verif_indus = rowSums(bdd_occupation[c("S59X04N", # Mining, stone and limestone quarrying
                                                      "S59X05N", # galvanization industry
                                                      "S59X06N", # chemical industry
                                                      "S59X07N", # chemistry - refinery
                                                      "S59X08N", # chemistry - dyeing
                                                      "S59X09N", # Chemistry - Chemical laboratory
                                                      "S59X10N", # rubber industry
                                                      "S59X11N", # textile industry
                                                      "S59X12N", # texttile dyeing
                                                      "S59X13N", # textile - spinning, weaving
                                                      "S59X14N", # leather tanning and dyeing
                                                      "S59X15N", # leather products (shoes)
                                                      "S59X16N", # wood industry
                                                      "S59X17N", # wood industry - furniture manufacturing 
                                                      "S59X18N", # metal processing
                                                      "S59X19N", # metal turning etc
                                                      "S59X20N", # metal welding 
                                                      "S59X21N", # metal painting 
                                                      "S59X22N", # foundry
                                                      "S59X23N", # foundry - steel rolling mill
                                                      "S59X24N", # foundry - alloy rolling mill
                                                      "S59X25N", # shipbuiding industry
                                                      "S59X26N", # electronics industry
                                                      "S59X27N", # glass industry
                                                      "S59X28N", # graphic industry
                                                      "S59X29N", # constrution industry 
                                                      "S59X30N", # construciton - roofer
                                                      "S59X31N", # construction - asphalt worker
                                                      "S59X32N", # construction - demolition worker
                                                      "S59X36N", # paper of paper product manufacturing 
                                                      "S59X37N", # asbestos or cement 
                                                      "S59X38N", # asbestos insulation
                                                      "S59X39N", # cement product industry
                                                      "S59X40N", # porcelain/faience/earthenware
                                                      "S59X42N", # painter
                                                      "S59X43N"  # welder (soudeur)
    )] == 1, na.rm = TRUE),
    occupation_verif_indus_percent = occupation_verif_indus/36*100,
    occupation_verif_other_blue = rowSums(bdd_occupation[c("S59X33N",  # transportation
                                                           "S59X34N", # transportation - truck driver
                                                           "S59X35N", # transportation - bus or taxi
                                                           "S59X41N", # butcher
                                                           "S59X45N", # tank attendant (pompiste)
                                                           "S59X46N"  # automechanic
    )] == 1, na.rm = TRUE),
    occupation_verif_other_blue_percent = occupation_verif_other_blue/6*100, 
    occupation = case_when(
      occupation_verif_white_collar_percent == 0 & 
        occupation_verif_agricultural_percent == 0 & 
        occupation_verif_indus_percent == 0 & 
        occupation_verif_other_blue_percent == 0 ~ NA_character_,
      occupation_verif_white_collar_percent == pmax(occupation_verif_white_collar_percent, occupation_verif_agricultural_percent, occupation_verif_indus_percent, occupation_verif_other_blue_percent) ~ "White collar",
      occupation_verif_agricultural_percent == pmax(occupation_verif_white_collar_percent, occupation_verif_agricultural_percent, occupation_verif_indus_percent, occupation_verif_other_blue_percent) ~ "Agricultural related work",
      occupation_verif_indus_percent == pmax(occupation_verif_white_collar_percent, occupation_verif_agricultural_percent, occupation_verif_indus_percent, occupation_verif_other_blue_percent) ~ "Industrial work",
      occupation_verif_other_blue_percent == pmax(occupation_verif_white_collar_percent, occupation_verif_agricultural_percent, occupation_verif_indus_percent, occupation_verif_other_blue_percent) ~ "Other blue collar"))


#bdd_occupation %>% filter(is.na(occupation)) %>% View()
# 191 persons (38%) did not answer or anwsered no (2)

#bdd_occupation %>% filter(occupation_verif>1) %>% select(sample, starts_with("occupation")) %>% View()  
# 178 persons (36%) responded to more than one occupation 

# bdd_occupation %>% 
#   filter(rowSums(across(c(occupation_verif_white_collar, 
#                           occupation_verif_agricultural,
#                           occupation_verif_indus, 
#                           occupation_verif_other_blue), ~ .x > 0)) > 1) %>% 
#   select(sample, starts_with("occupation")) %>% View()  
# 80 persons (16%) have occupations in more than only one big category

# bdd_occupation %>% 
#   filter(rowSums(across(c(occupation_verif_white_collar, 
#                           occupation_verif_agricultural,
#                           occupation_verif_indus, 
#                           occupation_verif_other_blue), ~ .x > 0)) > 1) %>% 
#   select(sample, 
#          "occupation_verif", 
#          "occupation_verif_white_collar_percent",
#          "occupation_verif_agricultural_percent",              
#          "occupation_verif_indus_percent",          
#          "occupation_verif_other_blue_percent",  
#          "occupation" ) %>% View()  
# For this 80 persons, we put then in the category that has the most "yes" answers


# bdd_occupation %>% 
#   select(sample, 
#          "occupation_verif", 
#          "occupation_verif_white_collar_percent",
#          "occupation_verif_agricultural_percent",              
#          "occupation_verif_indus_percent",          
#          "occupation_verif_other_blue_percent",  
#          "occupation" ) %>% View()  

bdd_occupation %>% select(occupation) %>% tbl_summary()
# Okay, done, but still 38% of missing data 

bdd_danish <- left_join(bdd_danish, bdd_occupation[c("sample", "occupation")], by = "sample")
bdd_danish <- bdd_danish %>% relocate(occupation, .after = education)

rm(bdd_occupation)

## data verification ----
### check birth date matching ----
bdd_danish <- bdd_danish %>%
  group_by(match) %>%
  mutate(case_date = birth_date[serial_no == "0"]) %>%
  mutate(kfun_age = as.numeric(difftime(birth_date, case_date, units = "days")) / 365.25, 
         kfun_age = as.integer(kfun_age)) %>%
  ungroup() %>%
  select("sample","serial_no","match", kfun_age, everything()) %>%
  select(-case_date)

 # bdd_danish %>% 
 #   select(sample, serial_no, match, birth_date, kfun_age) %>%
 #   filter(!kfun_age == 0) %>%
 #   View()
# 
# bdd_danish %>%
#   select(sample, serial_no, match, birth_date, sex, kfun_age) %>%
#   filter(match == 72) %>%
#   View()

### check sex matching ----
bdd_danish <- bdd_danish %>%
  group_by(match) %>%
  mutate(kfun_sex = ifelse(n_distinct(sex) == 1, 1, 2)) %>%
  ungroup()
 bdd_danish <- bdd_danish %>% relocate(kfun_sex, .after = kfun_age)

# bdd_danish %>% 
#   select(sample, serial_no, match, als, sex, kfun_sex) %>%
#   filter(!kfun_sex == 1) %>%
#   View()

### check chronology ----
descrip_num(data = bdd_danish, vars = c("baseline_age", "diagnosis_age", "death_age")) 
# bdd_danish %>%
#   mutate(
#     check_baseline_diagnosis = baseline_age < diagnosis_age,
#     check_baseline_death = ifelse(!is.na(death_date), baseline_age < death_date, TRUE),
#     check_diagnosis_death = ifelse(!is.na(death_date), diagnosis_age < death_date, TRUE)) %>%
#   filter(!check_baseline_diagnosis | !check_baseline_death | !check_diagnosis_death) %>%
#   View()

### decisions ----
# we decided to remove match 72 because the controls doesn't match in term of sex and age with the case
bdd_danish <- bdd_danish %>% filter(match != 72) 
bdd_POPs <- bdd_POPs %>% filter(!sample %in% c("208", "209", "210")) 
bdd_questionaire <- bdd_questionaire %>% filter(!sample %in% c("208", "209", "210")) 

# wrong als diagnosis dates for samples 283 and 300, update the data with the correct dates 
bdd_danish <- bdd_danish %>%
  mutate(
    als_date = as.character(als_date), 
    als_date = case_when(sample == "283" ~ "2006-10-24", 
                         sample == "300" ~ "1999-07-14",
                        TRUE ~ als_date), 
    als_date = as.Date(als_date), 
    diagnosis_age = as.numeric(difftime(als_date, birth_date, units = "days")) / 365.25)

# be careful, serial_no variable is not always 0, 1, 2 


# variable coding ----
bdd_danish <- bdd_danish %>%
  mutate(
    marital_status_2cat = 
      fct_recode(marital_status, 
                 "Other" = "Widowed",
                 "Other" = "Divorced",
                 "Other" = "Unmarried"), 
    smoking_2cat = fct_recode(smoking, 
                              "Ever" = "Current", 
                              "Ever" = "Previous"))

# missing values imputation ----
covar_a_imputer <- bdd_danish %>% 
  select(baseline_age, sex, marital_status_2cat, smoking, smoking_2cat, bmi, cholesterol, education) %>%
  colnames()

visu_na <- 
      bdd_danish %>% 
      select(all_of(covar_a_imputer)) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels()

bdd_danish_i <- bdd_danish %>%
  select(-starts_with("S59X"), -Comment) 
    # %>%
    #   mutate(
    #     marital_status = 
    #       fct_recode(marital_status, 
    #                  "1" = "Widowed", 
    #                  "2" = "Divorced",
    #                  "3" = "Married/cohabit", 
    #                  "4" = "Unmarried"))


bdd_danish_ii<- mice(bdd_danish_i, seed = 1996, maxit = 0)

method <- bdd_danish_ii$method
method[!names(method) %in% covar_a_imputer] <- ""
method 

pred <- quickpred(bdd_danish_i, mincor = 0.2, minpuc = 0.4)

bdd_danish_ii <- mice(bdd_danish_i, m=1, meth = method, pred = pred, seed = 11111)

bdd_danish_ii <- 
  complete(bdd_danish_ii) %>% 
  as.data.frame() 
# %>%
#   mutate(
#     marital_status = 
#       fct_recode(marital_status, 
#                  "Widowed" = "1", 
#                  "Divorced" = "2",
#                  "Married/cohabit" = "3", 
#                  "Unmarried" = "4"))

visu_na_i <- 
  bdd_danish_ii %>% 
  select(all_of(covar_a_imputer)) %>%
  is.na() %>%
  as.data.frame() %>%
  tbl_summary(
    missing = "no",
    statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
  bold_labels()

tbl_merge(tbls = list(visu_na, visu_na_i), 
          tab_spanner = c("**Not imputed**", "**Imputed**"))

tbl_merge(tbls = list(
  tbl1 = bdd_danish %>% select(marital_status_2cat) %>% tbl_summary(), 
  tbl2 = bdd_danish_ii %>% select(marital_status_2cat) %>% tbl_summary()), 
  tab_spanner = c("**Not imputed**", "**Imputed**"))

bdd_danish_ii <- 
  bdd_danish_ii %>% 
  select(sample, 
         marital_status_2cat_i = marital_status_2cat, 
         education_i = education, 
         smoking_i = smoking, 
         smoking_2cat_i = smoking_2cat, 
         cholesterol_i = cholesterol)

bdd_danish <- bdd_danish_ii %>%
  select("sample",
         "marital_status_2cat_i", 
         "education_i", 
         "smoking_i", 
         "smoking_2cat_i", 
         "cholesterol_i") %>%
  left_join(bdd_danish, by = "sample")

rm(pred, visu_na, visu_na_i, method, bdd_danish_i, bdd_danish_ii, covar_a_imputer)

# finnish data ----
bdd_finnish  <- 
  read_csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Ian/Excel Data/als_101322.csv", 
           show_col_types = FALSE) %>%
  rename(sample = BARCODE, 
         study = STUDY,
         match = STRATUM, 
         als = CASE,
         #kfun = KFUN, # not found in Ian's data 
         marital_status = MARITAL_STATUS,
         sex = SEX, 
         occupation = OCCUPATION, 
         smoking = SMOKING, 
         alcohol = ALCOHOL, 
         education = EDUCATION, 
         bmi = BMI, 
         #diabetes = DIABETES,  # not found in Ian's data 
         blod_sys = SBP, 
         blod_dias = DBP,
         baseline_date = DATE_BASELINE, 
         baseline_age = AGE_BASELINE, 
         death_age = age_death, 
         diagnosis_age = AGE_ALS, 
         follow_up = FOLLOW_UP_ALS, 
         "α_HCH" = "alfa_HCH", 
         "β_HCH" = "beta_HCH", 
         "γ_HCH" = "gamma_HCH", 
         "Transnonachlor" = "Trans_nonachlor") %>%
  mutate(sample = as.character(sample), 
         sex = as.character(sex), 
         sex = fct_recode(sex, "Male" = "1", "Female" = "2"),
         smoking = as.character(smoking), 
         smoking = fct_recode(smoking, 
                              "Never" = "1",
                              "Previous" = "2",
                              "Current" = "3",
                              "Current" = "4",
                              "Current" = "5"), 
         marital_status = as.character(marital_status), 
         marital_status = fct_recode(marital_status, 
                                     "Unmarried" = "1",
                                     "Married/cohabit" = "2",
                                     "Widowed" = "3",
                                     "Divorced" = "4"), 
         education = as.character(education), 
         education = fct_recode(education, 
                                "<7 years" = "1",
                                "7-12 years" = "2",
                                ">12 years" = "3"), 
         study = as.character(study), 
         study = fct_recode(study, 
                            "Finnish_1" = "1",
                            "Finnish_2" = "2",
                            "Finnish_3" = "3")) %>% 
  select(als, sample, study, match, sex, marital_status, occupation, smoking, alcohol, education, 
         bmi, chol, blod_sys, blod_dias, 
         baseline_age, diagnosis_age, 
         all_of(POPs), 
         -PeCB,  -γ_HCH, -Oxychlordane, -PCB_28, -PCB_52, 
         PeCB_all, Oxychlordane_all, PCB.28_all, PCB.52_all, 
         gamma.HCH_all) %>%
  rename(cholesterol = chol, PeCB = PeCB_all, Oxychlordane = Oxychlordane_all, PCB_28 = PCB.28_all, PCB_52 = PCB.52_all, 
         γ_HCH = gamma.HCH_all) %>%
  relocate(PeCB, .after = "diagnosis_age") %>%
  relocate(γ_HCH, .after = "β_HCH")

bdd_danish_red <- bdd_danish %>% select(als, sample, match, sex, marital_status, occupation, smoking, alcohol, education, 
                             bmi, cholesterol, blod_sys, blod_dias, 
                             baseline_age, diagnosis_age, 
                             all_of(POPs)) %>%
  mutate(study =  "Danish", 
         alcohol = alcohol*7) %>%
  select(sample, study, everything())

bdd_tot <- rbind(bdd_finnish, bdd_danish_red)

bdd_tot <- bdd_tot %>%
  mutate(PCB_DL = PCB_118 + PCB_156, 
         PCB_NDL = PCB_28 + PCB_52 + PCB_74 + PCB_99 + PCB_101 + PCB_138 + PCB_153 + 
           PCB_170 + PCB_180 + PCB_183 + PCB_187, 
         PCB_4 = PCB_118 + PCB_138 + PCB_153 + PCB_180, 
         ΣDDT = pp_DDT + pp_DDE,  
         ΣHCH = β_HCH + γ_HCH, 
         Σchlordane = Transnonachlor + Oxychlordane)


bdd_tot_long <- bdd_tot %>%
  select(sample, study, all_of(POPs), 
         "PCB_DL", "PCB_NDL", "PCB_4", "HCB", "ΣDDT", "ΣHCH", "β_HCH", "Σchlordane") %>%
  select(-α_HCH, -BDE_47, -BDE_99, -BDE_153) %>%
  pivot_longer(cols = c(-sample, -study), names_to = "POPS", values_to = "values")


# variable labels ---- 
var_label(bdd_danish) <- list(
  sample = "Identifcation", 
  match = "match", 
  sex = "Sex", 
  marital_status = "Marital status",
  marital_status_2cat = "Marital status", 
  marital_status_2cat_i = "Marital status", 
  occupation = "Occupation",
  smoking = "Smoking status", 
  smoking_i = "Smoking status", 
  smoking_2cat = "Smoking status", 
  smoking_2cat_i = "Smoking status", 
  alcohol = "Alcohol consumption (g/week)", 
  education = "Education", 
  education_i = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  cholesterol_i = "Serum cholesterol (mmol/L)",
  blod_sys = "Systolic blood presure (mmHg)",
  blod_dias = "Diastolic blood presure (mmHg)",
  baseline_age = "Age at baseline",
  diagnosis_age = "Age at ALS diagnosis", 
  death_age = "Age at death",
  PeCB = "Pentachlorobenzene (PeCB)",            
  HCB = "HCB",            
  α_HCH = "α-HCH",
  β_HCH = "β-HCH", 
  γ_HCH = "γ-HCH", 
  Oxychlordane = "Oxychlordane", 
  Transnonachlor = "Trans-nonachlor", 
  pp_DDT = "p,p'-DDT", 
  pp_DDE = "p,p'-DDE", 
  PCB_28 = "PCB-28", 
  PCB_52 = "PCB-52", 
  PCB_74 = "PCB-74", 
  PCB_99 = "PCB-99", 
  PCB_101 = "PCB-101", 
  PCB_118 = "PCB-118", 
  PCB_138 = "PCB-138", 
  PCB_153 = "PCB-153", 
  PCB_156 = "PCB-156", 
  PCB_170 = "PCB-170", 
  PCB_180 = "PCB-180", 
  PCB_183 = "PCB-183", 
  PCB_187 = "PCB-187", 
  BDE_47 = "BDE-47", 
  BDE_99 = "BDE-99", 
  BDE_153 = "BDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "PCB-118,138,153,180"
  )

var_label(bdd_tot) <- list(
  sample = "Identifcation", 
  match = "match", 
  sex = "Sex", 
  marital_status = "Marital status",
  occupation = "Occupation",
  smoking = "Smoking status (Never/Previous/Current)", 
  alcohol = "Alcohol consumption (g/week)", 
  education = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  blod_sys = "Systolic blood presure (mmHg)",
  blod_dias = "Diastolic blood presure (mmHg)",
  baseline_age = "Age at baseline",
  diagnosis_age = "Age at ALS diagnosis")

rm(bdd_danish_red, bdd_questionaire, bdd_POPs)


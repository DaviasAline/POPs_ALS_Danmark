# A. Davias 
# lipid work on the finish data 

# packages loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/0_functions.R")

# data loading ----
## danish data ----
bdd_danish <- 
  read_sas("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/pop_als_1ak.sas7bdat", NULL) |>
  mutate(code = as.character(code))
bdd_danish_POPs <- read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/2024_ALS-Denmark-POP-Final-Results.xlsx", 
                              sheet = "Results-blank-subtracted", 
                              range = "A3:AC505") |> 
  slice(-1) |>                                                                  # removing of the additional rows
  rename(code = Sample) 
bdd_danish_lipids <- 
  read.csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Danish Lipids/EPIC_final_results.csv", sep=";") |>
  filter(!grepl("T", X12))|>                                                    # removing of the additional rows
  pivot_wider(names_from = "ANALYSEX",                                          # pivotwider to get only one row per individual  
              values_from = "MEANVAL") |>
  rename(Barcode = BATCHNO)
bdd_danish_fattyacids <- 
  read.csv2("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Danish Lipids/ALS_plasma_FattyAcids_results.csv") |>
  rename(code = CODE) |>
  mutate(code = as.character(code))

bdd_danish <- left_join(bdd_danish, bdd_danish_POPs, by = "code")               # merging all the different datasets for the danish data
bdd_danish <- left_join(bdd_danish, bdd_danish_fattyacids, by = "code")
bdd_danish <- left_join(bdd_danish, bdd_danish_lipids, by = "Barcode")
rm(bdd_danish_POPs, bdd_danish_lipids, bdd_danish_fattyacids)

## finnish data ----
bdd_finnish <- 
  read.delim("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/THLBB2022_24/THLBB2020_21_H2011_phenotypes_14092022.txt") |>
  rename(PSEUDO_ID = BARCODE)
bdd_finnish_POPs_raw <- 
  read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/ALS_Autoklinikka_Mikro-POP_results.xlsx", 
             sheet = "Raw Data") |>
  slice(-1, -2) |>                                                              # removing of the additional rows
  rename_with(~ paste0(.x, "_raw"), .cols = !c("Batch", "Code")) |>             # adding "_raw" because the POPs variables are already present in the metadata dataset
  rename(PSEUDO_ID = Code) |>
  mutate(across(
    .cols = -c(Batch, PSEUDO_ID),
    .fns = ~ as.numeric(.)))
# bdd_finnish_POPs_not_raw <- 
#   read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/ALS_Autoklinikka_Mikro-POP_results.xlsx", 
#              sheet = "Results ") |>
#   slice(-1, -2) |>                                                              # removing of the additional rows
#   rename_with(~ paste0(.x, "_not_raw"), .cols = !c("Batch", "Code")) |>         # adding "_raw" because the POPs variables are already present in the metadata dataset
#   rename(PSEUDO_ID = Code)
bdd_finnish_lipids <- 
  read.csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_Lipid.csv", sep=";") |>
  filter(!is.na(PSEUDO_ID)) |>                                                  # removing of the additional rows
  pivot_wider(names_from = "ANALYSEX",                                          # pivotwider to get only one row per individual
              values_from = "MEANVAL") |>
  rename(BARCODE = Barcode)
bdd_finnish_fattyacids <- 
  read.csv2("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_FattyAcid.csv") |>
  filter(!is.na(PSEUDO_ID)) |>                                                  # removing of the additional rows
  rename(BARCODE = Barcode)                                               

bdd_finnish <- left_join(bdd_finnish, bdd_finnish_POPs_raw, by = c("PSEUDO_ID", "Batch"))     # merging all the different datasets for the finnish data
# bdd_finnish <- left_join(bdd_finnish, bdd_finnish_POPs_not_raw, by = c("PSEUDO_ID", "Batch"))
bdd_finnish <- full_join(bdd_finnish, bdd_finnish_lipids, by = "PSEUDO_ID")
bdd_finnish <- full_join(bdd_finnish, bdd_finnish_fattyacids, by = c("BARCODE", "PSEUDO_ID"))
bdd_finnish <- bdd_finnish |> filter(!PSEUDO_ID == "9047039803")                  # removing the personn that wants to be removed from the analyses 
rm(bdd_finnish_POPs_raw, bdd_finnish_lipids, bdd_finnish_fattyacids)

# data cleaning ----
## danish data ----
bdd_danish <- bdd_danish |>
  rename_with(~ gsub("PCB-", "PCB_", .x), .cols = contains("PCB-")) |>
  rename_with(~ gsub("BDE-", "PBDE_", .x), .cols = contains("BDE-")) |>
  select(Barcode, sample = code, serial_no = loebenr, match = saet,             # identification variables 
         Batch, n, X12, 
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
         contains("PCB_"),                                                      # POPs
         OCP_PeCB = "PeCB", 
         OCP_HCB = "HCB", 
         OCP_α_HCH = "α-HCH", 
         OCP_β_HCH = "β-HCH", 
         OCP_γ_HCH = "γ-HCH",                   
         OCP_oxychlordane = "Oxychlordane", 
         OCP_transnonachlor = "Transnonachlor",
         OCP_pp_DDT = "p,p'-DDT", 
         OCP_pp_DDE = "p,p'-DDE", 
         contains("PBDE_"),
         "myristic_acid_sat" = "MYR",                                           # saturated acids 
         "pentadecylic_acid_sat" = "PENT", 
         "palmitic_acid_sat" = "PAL",
         "margaric_acid_sat" = "HEPT", 
         "stearic_acid_sat" =  "STEA", 
         "arachidic_acid_sat" = "EIKO",
         "dehenic_acid_sat" = "DOKO", 
         "lignoceric_acid_sat" = "TCOSA", 
         "dma16x_acid_ω9" = "DMA16x",                                           # omega 9 acids
         "ptol2_acid_ω9" = "PTOL2", 
         "dma18x_acid_ω9" = "DMA18x", 
         "oleic_acid_ω9" = "OLE", 
         "gadoleic_acid_ω9" = "EIKE", 
         "nervonic_acid_ω9" = "NERVO", 
         "palmitoleic_acid_ω7" = "PTOL",                                        # omega 7 acids
         "i17_acid_ω7" = "i17", 
         "ai17_acid_ω7" = "ai17", 
         "vaccenic_acid_ω7" = "cVAKS",
         "rumenic_acid_ω6" = "CLA1",                                            # omega 6 acids
         "linoleic_acid_ω6" = "LA", 
         "dihomo_γ_linolenic_acid_ω6" = "DGLA", 
         "arachidonic_acid_ω6" = "ARA", 
         "adrenic_acid_ω6" = "DTETR",
         "α_linolenic_acid_ω3" = "ALA",                                         # omega 3 acids
         "timnodonic_acid_ω3" = "EPA", 
         "clupanodonic_acid_ω3" = "DPA", 
         "cervonic_acid_ω3" = "DHA",
         "fS_Kol" = "fS-Kol",                                                   # cholesterol            
         "fS_Trigly"= "fS-Trigly",                                              # triglycerides
         unit_lipids = UNIT, 
         unit_fatty_acids = Unit,
         "COMMISNO_bdd_lipids" = "COMMISNO",                                    # other variables 
         "ACCOUNT_bdd_lipids" = "ACCOUNT", 
         "RESULDAT_bdd_lipids" = "RESULDAT", 
         "TESTDAT_bdd_lipids" = "TESTDAT",
         "APPRIVID_bdd_lipids" = "APPRIVID", 
         "TESTNO_bdd_lipids" =  "TESTNO", 
         "SERIALNO_bdd_lipids" =  "SERIALNO", 
         "TEXT_bdd_lipids" = "TEXT",   
         "SampleStatusText_bdd_fattyacids" = "SampleStatusText", 
         "Comment_bdd_POPs" = "Comment", 
         everything())  |>                                                      # occupations
  select(-X1, -X7, -X3, -XDATE2, -LIM1, -LIM2, -LIM3, -LIM4, -X, -XInt3, -XINT1)# empty variables to delate


bdd_danish <- bdd_danish |>
  mutate(
    smoking = as.character(smoking),                                            # recoding some variables 
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
                     "Male" = "M")) |>                       
  mutate(across(c( "fS_Trigly", "fS_Kol"),                                      # adjusting the class of some variables
                ~as.numeric(gsub(",", ".", ., fixed = TRUE)))) |>
  mutate(across(c(contains("PCB_"), contains("OCP_"), contains("PBDE_")), as.numeric)) 


## finnish data ----
bdd_finnish  <-
  bdd_finnish |>
  rename_with(~ gsub("BDE_", "PBDE_", .x), .cols = contains("BDE_")) |>
  rename_with(~ gsub("BDE-", "PBDE_", .x), .cols = contains("BDE-")) |>
  rename_with(~ gsub("PCB-", "PCB_", .x), .cols = contains("PCB-")) |>
  select(sample = "PSEUDO_ID",
         study = STUDY,
         "SUBJECT_ID" , 
         "BARCODE", 
         "status_als", 
         match = STRATUM,
         als = CASE,
         als_class =  ALS_CLASS, 
         marital_status = MARITAL_STATUS,
         sex = SEX,
         occupation = OCCUPATION,
         occupation_cat = OCCUPATION_CAT,
         smoking = SMOKING,
         alcohol = ALCOHOL,
         education = EDUCATION,
         bmi = BMI,
         blod_sys = SBP,
         blod_dias = DBP,
         baseline_date = DATE_BASELINE,
         baseline_age = AGE_BASELINE,
         death_age = age_death,
         diagnosis_age = AGE_ALS,
         follow_up = FOLLOW_UP_ALS,
         "THAWED",  "MUNICIPALITY",  "LEVEL_URBANISATION",  
         "LEISURE_EXERCISE", "chol", "cholesterol", "GHEALTH", "HEALTH", "Batch", 
         "FOLLOW_UP_Death", "status_death",
         contains("PCB_"),
         OCP_PeCB = "PeCB", 
         OCP_PeCB_raw = "PeCB_raw", 
         # OCP_PeCB_not_raw = "PeCB_not_raw", 
         OCP_HCB = "HCB", 
         OCP_HCB_raw = "HCB_raw", 
         # OCP_HCB_not_raw = "HCB_not_raw", 
         "OCP_α_HCH" = "alfa_HCH",
         "OCP_α_HCH_raw" = "α-HCH_raw", 
         # "OCP_α_HCH_not_raw" = "α-HCH_not_raw", 
         "OCP_β_HCH" = "beta_HCH",
         "OCP_β_HCH_raw" = "β-HCH_raw", 
         # "OCP_β_HCH_not_raw" = "β-HCH_not_raw", 
         "OCP_γ_HCH" = "gamma_HCH",
         "OCP_γ_HCH_raw" = "γ-HCH_raw", 
         # "OCP_γ_HCH_not_raw" = "γ-HCH_not_raw", 
         OCP_oxychlordane = "Oxychlordane", 
         "OCP_oxychlordane_raw" =  "Oxy\r\nchlordane_raw", 
         # "OCP_oxychlordane_not_raw"  = "Oxy\r\nchlordane_not_raw", 
         "OCP_transnonachlor" = "Trans_nonachlor", 
         "OCP_transnonachlor_raw" = "Trans-nona\r\nchlor_raw", 
         # "OCP_transnonachlor_not_raw" = "Trans-nona\r\nchlor_not_raw", 
         OCP_pp_DDT = "pp_DDT",  
         "OCP_pp_DDT_raw" = "p,p'-DDT_raw",
         # "OCP_pp_DDT_not_raw" = "p,p'-DDT_not_raw", 
         OCP_pp_DDE = "pp_DDE",
         "OCP_pp_DDE_raw" = "p,p'-DDE_raw", 
         # "OCP_pp_DDE_not_raw" = "p,p'-DDE_not_raw", 
         contains("PBDE"),
         "myristic_acid_sat" = "MYR",                                           # saturated acids 
         "pentadecylic_acid_sat" = "PENT", 
         "palmitic_acid_sat" = "PAL",
         "margaric_acid_sat" = "HEPT", 
         "stearic_acid_sat" =  "STEA", 
         "arachidic_acid_sat" = "EIKO",
         "dehenic_acid_sat" = "DOKO", 
         "lignoceric_acid_sat" = "TCOSA", 
         "dma16x_acid_ω9" = "DMA16x",                                           # omega 9 acids
         "ptol2_acid_ω9" = "PTOL2", 
         "dma18x_acid_ω9" = "DMA18x", 
         "oleic_acid_ω9" = "OLE", 
         "gadoleic_acid_ω9" = "EIKE", 
         "nervonic_acid_ω9" = "NERVO", 
         "palmitoleic_acid_ω7" = "PTOL",                                        # omega 7 acids
         "i17_acid_ω7" = "i17", 
         "ai17_acid_ω7" = "ai17", 
         "vaccenic_acid_ω7" = "cVAKS",
         "rumenic_acid_ω6" = "CLA1",                                            # omega 6 acids
         "linoleic_acid_ω6" = "LA", 
         "dihomo_γ_linolenic_acid_ω6" = "DGLA", 
         "arachidonic_acid_ω6" = "ARA", 
         "adrenic_acid_ω6" = "DTETR",
         "α_linolenic_acid_ω3" = "ALA",                                         # omega 3 acids
         "timnodonic_acid_ω3" = "EPA", 
         "clupanodonic_acid_ω3" = "DPA", 
         "cervonic_acid_ω3" = "DHA",
         "fS_Kol" = "fS-Kol",                                                   # cholesterol            
         "fS_Trigly"= "fS-Trigly",                                              # triglycerides
         "S_Ca"  = "S-Ca",
         unit_lipids = UNIT, 
         unit_fatty_acids = Unit,
         "COMMISNO_bdd_lipids" = "COMMISNO",                                    # other variables 
         "ACCOUNT_bdd_lipids" = "ACCOUNT", 
         "RESULDAT_bdd_lipids" = "RESULDAT", 
         "TESTDAT_bdd_lipids" = "TESTDAT",
         "APPRIVID_bdd_lipids" = "APPRIVID", 
         "TESTNO_bdd_lipids" =  "TESTNO", 
         "SERIALNO_bdd_lipids" =  "SERIALNO", 
         "TEXT_bdd_lipids" = "TEXT",   
         everything()) |>
  mutate(sex = as.character(sex),
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
                            "Finnish_3" = "3")) |> 
         mutate(across(c( "fS_Trigly", "fS_Kol","S_Ca"),                        # adjusting the class of some variables
                       ~as.numeric(gsub(",", ".", ., fixed = TRUE)))) |>
  select(-X1, -X7, -X3, -XDATE2, -LIM1, -LIM2, -LIM3, -LIM4, -X, -XInt3, -XINT1, - X12) # empty variables to delate


# vector creation ----
POPs <- bdd_danish |> select(contains("PCB_"), contains("OCP_"), contains("PBDE_")) |> colnames()
POPs_quart <- paste0(POPs, "_quart")

POPs_group <- c("PCB_DL", "PCB_NDL", "PCB_4", "OCP_HCB", "ΣDDT", "OCP_β_HCH", "Σchlordane", "ΣPBDE")
POPs_group_quart <- paste0(POPs_group, "_quart")
POPs_group_outlier <- paste0(POPs_group, "_outlier")

POPs_group_quart_med <- paste0(POPs_group, "_quart_med")

POPs_included <- bdd_danish |> select(all_of(POPs)) |> select(-OCP_PeCB, - OCP_α_HCH, -OCP_γ_HCH) |> colnames()
POPs_included_quart <- paste0(POPs_included, "_quart")
POPs_included_outlier <- paste0(POPs_included, "_outlier")

fattyacids <- bdd_danish |> select(contains("_sat"), contains("_ω9"), contains("_ω7"), contains("_ω6"), contains("_ω3")) |> colnames()
explanatory <- c("pufas", "pufas_ω9", "pufas_ω3", "pufas_ω6", 
                 "rumenic_acid_ω6", "linoleic_acid_ω6", "dihomo_γ_linolenic_acid_ω6", "arachidonic_acid_ω6", "adrenic_acid_ω6",
                 "α_linolenic_acid_ω3", "timnodonic_acid_ω3", "clupanodonic_acid_ω3", "cervonic_acid_ω3") 

covariates_danish <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'cholesterol_i', 'marital_status_2cat_i', 'education_i')


POPs_finnish <- ifelse(POPs %in% c("PCB_28", "PCB_52", "OCP_PeCB", "OCP_α_HCH", "OCP_γ_HCH", 
                           "OCP_oxychlordane", "PBDE_47", "PBDE_99", "PBDE_153"), paste0(POPs, "_raw"), POPs)
covariates_finnish <- c('sex', 'baseline_age', 'smoking', 'bmi', 'cholesterol', 'marital_status')                                 # education removed because missing in one finnish cohort 

# variable creation ----
## danish data ----
bdd_danish <- bdd_danish |>
  mutate(
    study = "Danish", 
    als = ifelse(!is.na(als_date), 1, 0),                                       # metadata variables creation 
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
    follow_up = as.numeric(difftime(als_date, baseline_date, units = "days"))/365.25,
    marital_status_2cat = 
      fct_recode(marital_status, 
                 "Other" = "Widowed",
                 "Other" = "Divorced",
                 "Other" = "Unmarried"), 
    smoking_2cat = fct_recode(smoking, 
                              "Ever" = "Current", 
                              "Ever" = "Previous")) |>
  
  mutate(across(all_of(POPs), ~ factor(ntile(.x, 4),                            # POP variables creation
                                       labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart")) |>
  mutate(
    PCB_DL = PCB_118 + PCB_156, 
    PCB_NDL = PCB_28 + PCB_52 + PCB_74 + PCB_99 + PCB_101 + PCB_138 + PCB_153 + 
      PCB_170 + PCB_180 + PCB_183 + PCB_187, 
    PCB_4 = PCB_118 + PCB_138 + PCB_153 + PCB_180, 
    ΣPBDE = PBDE_47 + PBDE_99 + PBDE_153, 
    ΣDDT = OCP_pp_DDT + OCP_pp_DDE,  
    ΣHCH = OCP_β_HCH + OCP_γ_HCH, 
    Σchlordane = OCP_transnonachlor + OCP_oxychlordane) %>% # alpha HCH, gamma HCH exluded, beta HCH studied alone
  mutate(
    across(c("PCB_DL", 
             "PCB_NDL", 
             "PCB_4", 
             "ΣPBDE", 
             "ΣDDT", 
             "ΣHCH", 
             "Σchlordane"), ~ factor(ntile(.x, 4), labels = c("Q1", "Q2", "Q3", "Q4")),
           .names = "{.col}_quart")) %>%
  mutate(across(all_of(POPs_group), ~ {
    Q1 <- quantile(., 0.25, na.rm = TRUE)
    Q3 <- quantile(., 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    ifelse(. < lower_bound | . > upper_bound, NA, .)
  }, .names = "{.col}_outlier")) 


bdd_danish <- bdd_danish %>%
  replace_with_median(OCP_HCB, OCP_HCB_quart) %>%
  replace_with_median(PCB_DL, PCB_DL_quart) %>%
  replace_with_median(PCB_NDL, PCB_NDL_quart) %>%
  replace_with_median(PCB_4, PCB_4_quart) %>%
  replace_with_median(ΣPBDE, ΣPBDE_quart) %>%
  replace_with_median(ΣDDT, ΣDDT_quart) %>%
  replace_with_median(ΣHCH, ΣHCH_quart) %>%
  replace_with_median(Σchlordane, Σchlordane_quart) %>%
  replace_with_median(OCP_PeCB, OCP_PeCB_quart) %>%
  replace_with_median(OCP_β_HCH, OCP_β_HCH_quart)


bdd_danish <- bdd_danish |>
  mutate(
    pufas_ω9 = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9, 
    pufas_ω3 = α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3, 
    pufas_ω6 = rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6, 
    pufas = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9 +
      palmitoleic_acid_ω7 + i17_acid_ω7 + ai17_acid_ω7 + vaccenic_acid_ω7 +
      rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6 +
      α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3) 

## finnish data ----
bdd_finnish <- bdd_finnish |>
    mutate(across(all_of(POPs_finnish), ~ factor(ntile(.x, 4),                  # POP variables creation
                                         labels = c("Q1", "Q2", "Q3", "Q4")),
                  .names = "{.col}_quart")) |>
      mutate(
        PCB_DL = PCB_118 + PCB_156, 
        PCB_NDL = PCB_28_raw + PCB_52_raw + PCB_74 + PCB_99 + PCB_101 + PCB_138 + PCB_153 + 
          PCB_170 + PCB_180 + PCB_183 + PCB_187, 
        PCB_4 = PCB_118 + PCB_138 + PCB_153 + PCB_180, 
        ΣPBDE = PBDE_47_raw + PBDE_99_raw + PBDE_153_raw, 
        ΣDDT = OCP_pp_DDT + OCP_pp_DDE,  
        ΣHCH = OCP_β_HCH + OCP_γ_HCH_raw, 
        Σchlordane = OCP_transnonachlor + OCP_oxychlordane_raw) %>% # alpha HCH, gamma HCH exluded, beta HCH studied alone
      mutate(
        across(c("PCB_DL", 
                 "PCB_NDL", 
                 "PCB_4", 
                 "ΣPBDE", 
                 "ΣDDT", 
                 "ΣHCH", 
                 "Σchlordane"), ~ factor(ntile(.x, 4), labels = c("Q1", "Q2", "Q3", "Q4")),
               .names = "{.col}_quart")) %>%
      mutate(across(all_of(POPs_group), ~ {
        Q1 <- quantile(., 0.25, na.rm = TRUE)
        Q3 <- quantile(., 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        lower_bound <- Q1 - 1.5 * IQR
        upper_bound <- Q3 + 1.5 * IQR
        ifelse(. < lower_bound | . > upper_bound, NA, .)
      }, .names = "{.col}_outlier")) |>
  mutate(
    smoking_2cat = fct_recode(smoking, 
                              "Ever" = "Current", 
                              "Ever" = "Previous"), 
    marital_status_2cat = 
      fct_recode(marital_status, 
                 "Other" = "Widowed",
                 "Other" = "Divorced",
                 "Other" = "Unmarried"))
    

bdd_finnish <- bdd_finnish %>%
      replace_with_median(OCP_HCB, OCP_HCB_quart) %>%
      replace_with_median(PCB_DL, PCB_DL_quart) %>%
      replace_with_median(PCB_NDL, PCB_NDL_quart) %>%
      replace_with_median(PCB_4, PCB_4_quart) %>%
      replace_with_median(ΣPBDE, ΣPBDE_quart) %>%
      replace_with_median(ΣDDT, ΣDDT_quart) %>%
      replace_with_median(ΣHCH, ΣHCH_quart) %>%
      replace_with_median(Σchlordane, Σchlordane_quart) %>%
      replace_with_median(OCP_PeCB_raw, OCP_PeCB_raw_quart) %>%
      replace_with_median(OCP_β_HCH, OCP_β_HCH_quart) %>%
      replace_with_median(OCP_γ_HCH_raw, OCP_γ_HCH_raw_quart)

    
bdd_finnish <- bdd_finnish |>
      mutate(
        pufas_ω9 = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9, 
        pufas_ω3 = α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3, 
        pufas_ω6 = rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6, 
        pufas = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9 +
          palmitoleic_acid_ω7 + i17_acid_ω7 + ai17_acid_ω7 + vaccenic_acid_ω7 +
          rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6 +
          α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3)
    

# missing values imputation (covariates) ----
## danish data ----
covar_a_imputer <- bdd_danish %>% 
  select(marital_status_2cat, smoking, smoking_2cat, cholesterol, education) %>%
  colnames()

bdd_danish_i <- bdd_danish %>%
  select(sample, all_of(POPs), all_of(fattyacids), 
         "sex", "baseline_age", "smoking", "smoking_2cat", "bmi", "cholesterol", 
         "marital_status_2cat", "education", als) 

bdd_danish_ii <- mice(bdd_danish_i, seed = 1996, maxit = 0)

method <- bdd_danish_ii$method
method[!names(method) %in% covar_a_imputer] <- ""
method 

pred <- quickpred(bdd_danish_i, mincor = 0.2, minpuc = 0.4)

bdd_danish_ii <- mice(bdd_danish_i, m=1, meth = method, pred = pred, seed = 11111)

bdd_danish_ii <- 
  complete(bdd_danish_ii) %>% 
  as.data.frame() 

tbl_merge(tbls = list(
  visu_na_i = 
    bdd_danish_ii %>% 
    select(all_of(covar_a_imputer)) %>%
    is.na() %>%
    as.data.frame() %>%
    tbl_summary(
      missing = "no",
      statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
    bold_labels(), 
  visu_na = 
    bdd_danish %>% 
    select(all_of(covar_a_imputer)) %>%
    is.na() %>%
    as.data.frame() %>%
    tbl_summary(
      missing = "no",
      statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
    bold_labels()), 
          tab_spanner = c("**Imputed**", "**Not imputed**"))

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

rm(pred, method, bdd_danish_i, bdd_danish_ii, covar_a_imputer)

## finnish data ----
# no missing values on the covariates of the finnish cohorts except education that is missing for one entire cohort

# variable labels ----
var_label(bdd_danish) <- list(
  sample = "Identifcation", 
  match = "match", 
  sex = "Sex", 
  marital_status = "Marital status",
  marital_status_2cat = "Marital status", 
  marital_status_2cat_i = "Marital status", 
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
  OCP_PeCB = "Pentachlorobenzene (PeCB)",            
  OCP_HCB = "HCB",            
  OCP_α_HCH = "α-HCH",
  OCP_β_HCH = "β-HCH", 
  OCP_γ_HCH = "γ-HCH", 
  OCP_oxychlordane = "Oxychlordane", 
  OCP_transnonachlor = "Trans-nonachlor", 
  OCP_pp_DDT = "p,p'-DDT", 
  OCP_pp_DDE = "p,p'-DDE", 
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
  PBDE_47 = "BDE-47", 
  PBDE_99 = "BDE-99", 
  PBDE_153 = "BDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "PCB-118,138,153,180",
  "myristic_acid_sat" =  "Myristic acid", 
  "pentadecylic_acid_sat" = "Pentadecylic acid", 
  "palmitic_acid_sat" =  "Palmitic acid", 
  "margaric_acid_sat" =  "Margaric acid", 
  "stearic_acid_sat" = "Stearic acid",           
  "arachidic_acid_sat" = "Arachidic acid",
  "dehenic_acid_sat" = "Dehenic acid", 
  "lignoceric_acid_sat" = "Lignoceric acid", 
  "dma16x_acid_ω9" = "dma16x acid ω9", 
  "ptol2_acid_ω9"  = "ptol2 acid ω9", 
  "dma18x_acid_ω9" = "dma18x acid ω9", 
  "oleic_acid_ω9" = "Oleic acid ω9", 
  "gadoleic_acid_ω9" = "Gadoleic acid ω9", 
  "nervonic_acid_ω9" = "Nervonic acid ω9",
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω7", 
  "i17_acid_ω7" = "i17 acid ω7", 
  "ai17_acid_ω7" = "ai17 acid ω7", 
  "vaccenic_acid_ω7" = "Vaccenic acid ω7", 
  "rumenic_acid_ω6" = "Rumenic acid ω6", 
  "linoleic_acid_ω6" = "Linoleic acid ω6",          
  "dihomo_γ_linolenic_acid_ω6" = "Dihomo-γ-linolenic acid ω6", 
  "arachidonic_acid_ω6" = "Arachidonic acid ω6", 
  "adrenic_acid_ω6" =  "Adrenic acid ω6", 
  "α_linolenic_acid_ω3" = "α-linolenic acid (ALA) ω3", 
  "timnodonic_acid_ω3" = "Timnodonic acid (EPA) ω3", 
  "clupanodonic_acid_ω3"  ="Clupanodonic acid (DPA) ω3", 
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3", 
  "pufas_ω9" = "ω9 instaturated acids", 
  "pufas_ω6" = "ω6 insaturated acids", 
  "pufas_ω3" = "ω3 instaurated acids", 
  "pufas" = "Insaturared acids")

var_label(bdd_finnish) <- list(
  sample = "Identifcation", 
  match = "match", 
  sex = "Sex", 
  marital_status = "Marital status",
  smoking = "Smoking status", 
  alcohol = "Alcohol consumption (g/week)", 
  education = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  baseline_age = "Age at baseline",
  diagnosis_age = "Age at ALS diagnosis", 
  death_age = "Age at death",
  OCP_PeCB = "Pentachlorobenzene (PeCB)",            
  OCP_HCB = "HCB",            
  OCP_α_HCH = "α-HCH",
  OCP_β_HCH = "β-HCH", 
  OCP_γ_HCH = "γ-HCH", 
  OCP_oxychlordane = "Oxychlordane", 
  OCP_transnonachlor = "Trans-nonachlor", 
  OCP_pp_DDT = "p,p'-DDT", 
  OCP_pp_DDE = "p,p'-DDE", 
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
  PBDE_47 = "BDE-47", 
  PBDE_99 = "BDE-99", 
  PBDE_153 = "BDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "PCB-118,138,153,180",
  "myristic_acid_sat" =  "Myristic acid", 
  "pentadecylic_acid_sat" = "Pentadecylic acid", 
  "palmitic_acid_sat" =  "Palmitic acid", 
  "margaric_acid_sat" =  "Margaric acid", 
  "stearic_acid_sat" = "Stearic acid",           
  "arachidic_acid_sat" = "Arachidic acid",
  "dehenic_acid_sat" = "Dehenic acid", 
  "lignoceric_acid_sat" = "Lignoceric acid", 
  "dma16x_acid_ω9" = "dma16x acid ω9", 
  "ptol2_acid_ω9"  = "ptol2 acid ω9", 
  "dma18x_acid_ω9" = "dma18x acid ω9", 
  "oleic_acid_ω9" = "Oleic acid ω9", 
  "gadoleic_acid_ω9" = "Gadoleic acid ω9", 
  "nervonic_acid_ω9" = "Nervonic acid ω9",
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω7", 
  "i17_acid_ω7" = "i17 acid ω7", 
  "ai17_acid_ω7" = "ai17 acid ω7", 
  "vaccenic_acid_ω7" = "Vaccenic acid ω7", 
  "rumenic_acid_ω6" = "Rumenic acid ω6", 
  "linoleic_acid_ω6" = "Linoleic acid ω6",          
  "dihomo_γ_linolenic_acid_ω6" = "Dihomo-γ-linolenic acid ω6", 
  "arachidonic_acid_ω6" = "Arachidonic acid ω6", 
  "adrenic_acid_ω6" =  "Adrenic acid ω6", 
  "α_linolenic_acid_ω3" = "α-linolenic acid (ALA) ω3", 
  "timnodonic_acid_ω3" = "Timnodonic acid (EPA) ω3", 
  "clupanodonic_acid_ω3"  ="Clupanodonic acid (DPA) ω3", 
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3", 
  "pufas_ω9" = "ω9 instaturated acids", 
  "pufas_ω6" = "ω6 insaturated acids", 
  "pufas_ω3" = "ω3 instaurated acids", 
  "pufas" = "Insaturared acids")

# merged dataset ----
bdd_danish_red <- bdd_danish |> 
  select(sample, als, study, 
         all_of(covariates_danish), 
         "baseline_age", "death_age", "diagnosis_age",  alcohol, smoking, marital_status_2cat, 
         all_of(POPs), 
         all_of(POPs_group), 
         all_of(POPs_group_quart), 
         all_of(POPs_group_quart_med),
         all_of(fattyacids), 
         all_of(explanatory)) |>
  rename(smoking_2cat = smoking_2cat_i, 
         cholesterol = cholesterol_i, 
         education = education_i)

bdd_finnish_red <- bdd_finnish |> 
  select(sample, als, study, 
         all_of(covariates_finnish), education, alcohol, smoking_2cat, marital_status_2cat, 
         "baseline_age", "death_age", "diagnosis_age",                   
         all_of(POPs_finnish), 
         all_of(POPs_group), 
         all_of(POPs_group_quart), 
         all_of(POPs_group_quart_med),
         all_of(fattyacids), 
         all_of(explanatory)) |>
  rename_with(~ gsub("_raw", "", .x)) |>
  mutate(sample = as.character(sample))

bdd <- bind_rows(bdd_danish_red, bdd_finnish_red)
rm(bdd_danish_red, bdd_finnish_red)

# analyses ----




## main analysis ----
model1 <- map_dfr(explanatory, function(var) {                                   # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  
  tibble(
    model = "model 1",
    variable = var,
    df = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

model2 <- map_dfr(explanatory, function(var) {                                   # map_dfr() met tout dans un seul dataframe par rapport a map() qui renvoit une liste
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model)
  
  tibble(
    model = "model 2",
    variable = var,
    df = model_summary$term,
    OR = exp(model_summary$estimate),
    lower_CI = exp(model_summary$estimate - 1.96 * model_summary$std.error),
    upper_CI = exp(model_summary$estimate + 1.96 * model_summary$std.error),
    `p-value` = model_summary$p.value)
})

main_results <- bind_rows(model1, model2) |>
  filter(variable == df) |>
  mutate(
    OR = as.numeric(sprintf("%.1f", OR)),
    lower_CI = as.numeric(sprintf("%.1f", lower_CI)),
    upper_CI = as.numeric(sprintf("%.1f", upper_CI)),, 
    "95%CI" = paste(lower_CI, ", ", upper_CI, sep = ''),
    `p-value_raw`= `p-value`, 
    `p-value` = ifelse(`p-value` < 0.01, "<0.01", number(`p-value`, accuracy = 0.01, decimal.mark = "."))) 

rm(model1, model2)

## results presentation ----
### table 1 ----
table_1 <- bdd_danish |>
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) |>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status_2cat, education, alcohol, smoking_2cat, bmi, cholesterol) |>
  tbl_summary(by = als, 
              missing = 'no', 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1)) |>
  bold_labels() |>
  add_p(include = -diagnosis_age) |>
  add_overall() |>
  add_n() |>
  as_flex_table() |>
  # font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  


### table 2 ----
table_2 <- tbl_merge(
  tbls = list(
    tbl_1 = bdd_danish |>
      select(als, all_of(covariates_danish)) |>
      tbl_uvregression(
        y = als,
        method = glm,
        method.args = list(family = binomial), 
        exponentiate = TRUE, 
        estimate_fun = label_number(accuracy = 0.1, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun)|>
      bold_labels(), 
    tbl_2 = clogit(als ~ sex + baseline_age +
                     smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
                   data = bdd_danish) |>
      tbl_regression(exponentiate = TRUE, 
                     estimate_fun = label_number(accuracy = .1, decimal.mark = "."),
                     pvalue_fun = custom_pvalue_fun) |>
      bold_labels()), 
  tab_spanner = c("**Univariate**", "**Adjusted**")) |> 
  as_flex_table()|>
  # font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all")  |>
  add_footer_lines(
    "1Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category.
  2CI: Confidence interval.
  3Estimated risk of ALS when the characteristic is increasing by one unit, or compared to the reference category; adjusted for all the variables in the table.")

### figure 1 ----
figure_1 <- bdd_danish |>
  select(als, all_of(explanatory)) |>
  pivot_longer(cols = -als, names_to = "PUFAs", values_to = "Values") |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1"), 
    PUFAs =   fct_relevel(PUFAs, 
      "pufas", "pufas_ω9", "pufas_ω6", "pufas_ω3", "linoleic_acid_ω6",
      "dihomo_γ_linolenic_acid_ω6", "arachidonic_acid_ω6", "adrenic_acid_ω6",
      "rumenic_acid_ω6", "α_linolenic_acid_ω3", "timnodonic_acid_ω3",
      "clupanodonic_acid_ω3", "cervonic_acid_ω3"), 
    PUFAs = fct_rev(PUFAs)) |> 
  arrange(PUFAs) |>
  ggplot() +
  aes(x = Values, y = PUFAs, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_x_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(fill = "ALS ", x = "Values (pg/ml, log transformed)", y = "PUFAs") +
  theme_lucid() + 
  theme(legend.position = "bottom")


### figure 2 ----
figure_2 <- main_results %>% 
  filter(!variable == "adrenic_acid_ω6") |>
  mutate(p.value_shape = ifelse('p-value_raw'<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Base model" = "model 1",
                            "Adjusted model" = "model 2"),
         model = fct_relevel(model, 'Base model', 'Adjusted model')) %>%
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(cols = dplyr::vars(model), switch = "y", scales = "free_x") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "PUFAs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()


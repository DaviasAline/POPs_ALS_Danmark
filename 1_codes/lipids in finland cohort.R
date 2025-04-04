# A. Davias 
# lipid work on the finish data 

# packages loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/0_functions.R")
# source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
# rm(OCPs, OCPs_quart, PBDEs, PBDEs_quart, PCBs, PCBs_quart, POPs, POPs_group, 
#    POPs_group_outlier, POPs_group_quart, POPs_group_quart_med, POPs_included, 
#    POPs_included_quart, POPs_quart, bdd_loq)

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
  read.delim("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/THLBB2022_24/THLBB2020_21_H2011_phenotypes_14092022.txt")
bdd_finnish_POPs_raw <- 
  read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/ALS_Autoklinikka_Mikro-POP_results.xlsx", 
             sheet = "Raw Data") |>
  slice(-1, -2) |>                                                              # removing of the additional rows
  rename_with(~ paste0(.x, "_raw"), .cols = !c("Batch", "Code")) |>             # adding "_raw" because the POPs variables are already present in the metadata dataset
  rename(BARCODE = Code)
bdd_finnish_POPs_not_raw <- 
  read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Finnish data/ALS_Autoklinikka_Mikro-POP_results.xlsx", 
             sheet = "Results ") |>
  slice(-1, -2) |>                                                              # removing of the additional rows
  rename_with(~ paste0(.x, "_not_raw"), .cols = !c("Batch", "Code")) |>         # adding "_raw" because the POPs variables are already present in the metadata dataset
  rename(BARCODE = Code)
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

bdd_finnish <- left_join(bdd_finnish, bdd_finnish_POPs_raw, by = "BARCODE")     # merging all the different datasets for the finnish data
bdd_finnish <- left_join(bdd_finnish, bdd_finnish_POPs_not_raw, by = "BARCODE")
bdd_finnish <- left_join(bdd_finnish, bdd_finnish_lipids, by = "BARCODE")
bdd_finnish <- left_join(bdd_finnish, bdd_finnish_fattyacids, by = "BARCODE")
bdd_finnish <- bdd_finnish |> filter(!BARCODE == "9047039803")                  # removing the personn that wants to be removed from the analyses 
rm(bdd_finnish_POPs_raw, bdd_finnish_POPs_not_raw, bdd_finnish_lipids, bdd_finnish_fattyacids)


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
  read_csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Ian/Excel Data/als_101322.csv", 
           show_col_types = FALSE) %>%
  rename(sample = BARCODE, 
         study = STUDY,
         match = STRATUM, 
         als = CASE,
         marital_status = MARITAL_STATUS,
         sex = SEX, 
         occupation = OCCUPATION, 
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



# vector creation ----
POPs <- bdd_danish |> select(contains("PCB_"), contains("OCP_"), contains("PBDE_")) |> colnames()
POPs_quart <- paste0(POPs, "_quart")

POPs_group <- c("PCB_DL", "PCB_NDL", "PCB_4", "OCP_HCB", "ΣDDT", "OCP_β_HCH", "Σchlordane", "ΣPBDE")
POPs_group_quart <- paste0(POPs_group, "_quart")
POPs_group_outlier <- paste0(POPs_group, "_outlier")

POPs_group_quart_med <- paste0(POPs, "quart_med")

POPs_included <- bdd_danish |> select(all_of(POPs)) |> select(-OCP_PeCB, - OCP_α_HCH, -OCP_γ_HCH) |> colnames()
POPs_included_quart <- paste0(POPs_included, "_quart")
POPs_included_outlier <- paste0(POPs_included, "_outlier")

fattyacids <- bdd_danish |> select(contains("_sat"), contains("_ω9"), contains("_ω7"), contains("_ω6"), contains("_ω3")) |> colnames()

# variable creation ----
## danish data ----
bdd_danish <- bdd_danish |>
  mutate(
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
rm(replace_with_median)

bdd_danish <- bdd_danish |>
  mutate("pufas_ω3" = "α_linolenic_acid_ω3" + "timnodonic_acid_ω3" + "clupanodonic_acid_ω3" + "cervonic_acid_ω3", 
         "pufas_ω6" = "rumenic_acid_ω6" +  "linoleic_acid_ω6" + "dihomo_γ_linolenic_acid_ω6" + "arachidonic_acid_ω6" + "adrenic_acid_ω6")

## finnish data ----

 
# variable labels ----
var_label(bdd_finnish) <- list(
  sample = "Identifcation", 
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
  "palmitoleic_acid_ω9" = "Palmitoleic acid ω9", 
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
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3")






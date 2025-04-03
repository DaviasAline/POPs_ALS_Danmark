# A. Davias 
# lipid work on the finish data 

# packages loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/0_functions.R")
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
rm(OCPs, OCPs_quart, PBDEs, PBDEs_quart, PCBs, PCBs_quart, POPs, POPs_group, 
   POPs_group_outlier, POPs_group_quart, POPs_group_quart_med, POPs_included, 
   POPs_included_quart, POPs_quart, bdd_loq)

# data loading ----
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


bdd_fattyacids_finnish <- 
  read.csv2("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_FattyAcid.csv") |>
  mutate(PSEUDO_ID = as.character(PSEUDO_ID)) |>
  filter(!is.na(PSEUDO_ID)) |>     # on enleve les lignes supplementaires
  select("sample" = "PSEUDO_ID", 
         "Barcode", "Unit", 
         "myristic_acid_sat" = "MYR", 
         "pentadecylic_acid_sat" = "PENT", 
         "palmitic_acid_sat" = "PAL",
         "margaric_acid_sat" = "HEPT", 
         "stearic_acid_sat" =  "STEA", 
         "arachidic_acid_sat" = "EIKO",
         "dehenic_acid_sat" = "DOKO", 
         "lignoceric_acid_sat" = "TCOSA", 
         "dma16x_acid_ω9" = "DMA16x", 
         "ptol2_acid_ω9" = "PTOL2", 
         "dma18x_acid_ω9" = "DMA18x", 
         "oleic_acid_ω9" = "OLE", 
         "gadoleic_acid_ω9" = "EIKE", 
         "nervonic_acid_ω9" = "NERVO", 
         "palmitoleic_acid_ω7" = "PTOL", 
         "i17_acid_ω7" = "i17", 
         "ai17_acid_ω7" = "ai17", 
         "vaccenic_acid_ω7" = "cVAKS",
         "rumenic_acid_ω6" = "CLA1", 
         "linoleic_acid_ω6" = "LA", 
         "dihomo_γ_linolenic_acid_ω6" = "DGLA", 
         "arachidonic_acid_ω6" = "ARA", 
         "adrenic_acid_ω6" = "DTETR",
         "α_linolenic_acid_ω3" = "ALA", 
         "timnodonic_acid_ω3" = "EPA", 
         "clupanodonic_acid_ω3" = "DPA", 
         "cervonic_acid_ω3" = "DHA")


bdd_lipids_finnish <- 
  read.csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_Lipid.csv", sep=";") |>
  rename(sample = PSEUDO_ID) |>
  mutate(sample = as.character(sample)) |>
  filter(!is.na(sample)) |>   # on enleve les lignes supplementaires
  pivot_wider(names_from = "ANALYSEX", 
              values_from = "MEANVAL") |>
  select(-X12, -X1, -X7, -X3, -XDATE2, -LIM1, -LIM2, -LIM3, -LIM4, -X, -XInt3, -"XINT1") |>
  mutate(across(c("S-Ca", "fS-Trigly", "fS-Kol"), ~as.numeric(gsub(",", ".", ., fixed = TRUE))))
    

## danish data ----
bdd_lipids_danish <- 
  read.csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Danish Lipids/EPIC_final_results.csv", sep=";") |>
  pivot_wider(names_from = "ANALYSEX", 
              values_from = "MEANVAL")

bdd_fattyacids_danish <- 
  read.csv2("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Danish Lipids/ALS_plasma_FattyAcids_results.csv")



# vector creation ----
POPs <- c("PCB_118", "PCB_156", "PeCB", 
          "PCB_28", "PCB_52", "PCB_74", "PCB_99", "PCB_101",  "PCB_138", "PCB_153",  "PCB_170", "PCB_180", 
          "PCB_183", "PCB_187",
          "HCB", "β_HCH", "γ_HCH", "Oxychlordane", "Transnonachlor", "pp_DDT", "pp_DDE")

fatty_acids <- bdd_fattyacids_finnish |> select(-sample, -Barcode, -Unit) |> colnames()

 
# merging datasets ----
bdd_finnish <- full_join(bdd_finnish, bdd_fattyacids_finnish |> select(-'Unit'), by = 'sample')
bdd_finnish <- full_join(bdd_finnish, bdd_lipids_finnish |> select(sample, 'Barcode', "S-Ca", "fS-Trigly", "fS-Kol"), by = c('sample', 'Barcode'))
bdd_finnish <- bdd_finnish |> filter(!sample == "9047039803") # on enleve la personne qui veut etre enlever des analyses

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
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω9", 
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

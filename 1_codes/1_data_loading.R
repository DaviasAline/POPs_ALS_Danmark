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

bdd_danish_proteomic <- read_excel("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Danish EPIC data/Weisskopf_NEX_IR_MET_NPX.xlsx") |>
  mutate(UniProt = as.factor(UniProt),                                          # what is index variable ? 90 different values 
         Assay = as.factor(Assay), 
         Panel = as.factor(Panel), 
         Panel_Version = as.factor(Panel_Version), 
         PlateID = as.factor(PlateID), 
         QC_Warning = as.factor(QC_Warning), 
         Normalization = as.factor(Normalization), 
         "Olink NPX Signature Version" = as.factor("Olink NPX Signature Version"))

bdd_danish_proteomic_wide <- bdd_danish_proteomic %>%
  select(SampleID, Assay, NPX) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Assay,
    values_from = NPX) |>
  rename(code = SampleID)

bdd_danish <- left_join(bdd_danish, bdd_danish_POPs, by = "code")               # merging all the different datasets for the danish data
bdd_danish <- left_join(bdd_danish, bdd_danish_fattyacids, by = "code")
bdd_danish <- left_join(bdd_danish, bdd_danish_lipids, by = "Barcode")
bdd_danish <- left_join(bdd_danish, bdd_danish_proteomic_wide, by = "code")
rm(bdd_danish_POPs, bdd_danish_lipids, bdd_danish_fattyacids, bdd_danish_proteomic_wide)

# bdd_danish |> filter(saet == 72) |> View()                           
bdd_danish <- bdd_danish %>% filter(!code %in% c("208", "209", "210"))        # we decided to remove match 72 because the controls doesn't match in term of sex and age with the case

bdd_danish_loq <- 
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
    variable = gsub("-", "_", variable), 
    variable = fct_recode(variable, 
      "PBDE_153" = "BDE_153",
      "PBDE_47" = "BDE_47",
      "PBDE_99" = "BDE_99",
      "OCP_HCB" = "HCB",
      "OCP_xychlordane" = "Oxychlordane",
      "OCP_PeCB" = "PeCB",
      "OCP_pp_DDE" = "pp_DDE",
      "OCP_pp_DDT" = "pp_DDT",
      "OCP_transnonachlor" = "Transnonachlor",
      "OCP_α_HCH" = "α_HCH",
      "OCP_β_HCH" = "β_HCH",
      "OCP_γ_HCH" = "γ_HCH"))  

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
bdd_finnish <- bdd_finnish |> filter(!PSEUDO_ID == "9047039803")                # removing the personn that wants to be removed from the analyses 
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
         municipality = MUNICIPALITY, 
         level_urbanization = "LEVEL_URBANISATION",  
         thawed = THAWED, 
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
         follow_up_als = FOLLOW_UP_ALS,
         "LEISURE_EXERCISE", "chol", "cholesterol", "GHEALTH", "HEALTH", "Batch", 
         follow_up_death ="FOLLOW_UP_Death", 
         "status_death",
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
                            "FMC" = "1",
                            "FMCF" = "2",
                            "MFH" = "3"), 
         municipality = as.factor(as.character(municipality)), 
         level_urbanization = as.factor(as.character(level_urbanization)), 
         # status_death = as.factor(as.character(status_death)),                # do not factor the status_death variable (if yes, there is an issue in the cox models)
         follow_up_death = follow_up_death/ 30.44,
         thawed = as.factor(as.character(thawed))) |> 
         mutate(across(c( "fS_Trigly", "fS_Kol","S_Ca"),                        # adjusting the class of some variables
                       ~as.numeric(gsub(",", ".", ., fixed = TRUE)))) |>
  select(-X1, -X7, -X3, -XDATE2, -LIM1, -LIM2, -LIM3, -LIM4, -X, -XInt3, -XINT1, - X12) # empty variables to delate


# vector creation ----
POPs <- bdd_danish |> select(contains("PCB_"), contains("OCP_"), contains("PBDE_")) |> colnames()
POPs_quart <- paste0(POPs, "_quart")

POPs_tot <- c("PCB_4",  
              "PCB_DL", "PCB_118", "PCB_156", 
              "PCB_NDL", "PCB_28", "PCB_52", "PCB_74", "PCB_99", "PCB_101", 
              "PCB_138", "PCB_153", "PCB_170", "PCB_180", "PCB_183", "PCB_187", 
              "OCP_HCB",  
              "ΣDDT", "OCP_pp_DDT",  "OCP_pp_DDE",
              "OCP_α_HCH", "OCP_β_HCH", "OCP_γ_HCH", 
              "Σchlordane", "OCP_oxychlordane", "OCP_transnonachlor", 
              "OCP_PeCB",   
              "ΣPBDE", "PBDE_47", "PBDE_99", "PBDE_153")

POPs_group <- c("PCB_4", "PCB_DL", "PCB_NDL", "OCP_HCB", "ΣDDT", "OCP_β_HCH", "Σchlordane", "ΣPBDE")
POPs_group_quart <- paste0(POPs_group, "_quart")
POPs_group_outlier <- paste0(POPs_group, "_outlier")
POPs_group_sd <- paste0(POPs_group, "_sd")
POPs_group_quart_med <- paste0(POPs_group, "_quart_med")

POPs_group_finnish <- c("PCB_4", "PCB_DL", "PCB_NDL", "OCP_HCB", "ΣDDT", "ΣHCH", "OCP_β_HCH", "OCP_γ_HCH", "Σchlordane", "OCP_PeCB")
POPs_group_quart_finnish<- paste0(POPs_group_finnish, "_quart")
POPs_group_outlier_finnish <- paste0(POPs_group_finnish, "_outlier")
POPs_group_sd_finnish <- paste0(POPs_group_finnish, "_sd")
POPs_group_quart_med_finnish <- paste0(POPs_group_finnish, "_quart_med")

POPs_included <- bdd_danish |> select(all_of(POPs)) |> select(-OCP_PeCB, - OCP_α_HCH, -OCP_γ_HCH) |> colnames()
POPs_included_quart <- paste0(POPs_included, "_quart")
POPs_included_outlier <- paste0(POPs_included, "_outlier")

fattyacids <- bdd_danish |> 
  select(contains("_sat"), contains("_ω9"), contains("_ω7"), contains("_ω6"), contains("_ω3")) |> 
  colnames()
fattyacids <- c("pufas", "pufas_ω9", "pufas_ω7", "pufas_ω6", "pufas_ω3", "ratio_ω6_ω3", fattyacids)
explanatory_raw <- c("pufas", "pufas_ω9", "pufas_ω7", "pufas_ω6", "pufas_ω3", "ratio_ω6_ω3",
                 "rumenic_acid_ω6", "linoleic_acid_ω6", "dihomo_γ_linolenic_acid_ω6", "arachidonic_acid_ω6", "adrenic_acid_ω6",
                 "α_linolenic_acid_ω3", "timnodonic_acid_ω3", "clupanodonic_acid_ω3", "cervonic_acid_ω3") 
explanatory <- c("pufas_sd", "pufas_ω9_sd", "pufas_ω7_sd", "pufas_ω6_sd", "pufas_ω3_sd", "ratio_ω6_ω3_sd",
                 "rumenic_acid_ω6_sd", "linoleic_acid_ω6_sd", "dihomo_γ_linolenic_acid_ω6_sd", "arachidonic_acid_ω6_sd", "adrenic_acid_ω6_sd",
                 "α_linolenic_acid_ω3_sd", "timnodonic_acid_ω3_sd", "clupanodonic_acid_ω3_sd", "cervonic_acid_ω3_sd") 
explanatory_quart <- c("pufas_quart", "pufas_ω9_quart", "pufas_ω7_quart", "pufas_ω6_quart", "pufas_ω3_quart", "ratio_ω6_ω3_quart",
                 "rumenic_acid_ω6_quart", "linoleic_acid_ω6_quart", "dihomo_γ_linolenic_acid_ω6_quart", "arachidonic_acid_ω6_quart", "adrenic_acid_ω6_quart",
                 "α_linolenic_acid_ω3_quart", "timnodonic_acid_ω3_quart", "clupanodonic_acid_ω3_quart", "cervonic_acid_ω3_quart") 

covariates_danish <- c('sex', 'baseline_age', 'smoking_2cat_i', 'bmi', 'cholesterol_i', 'marital_status_2cat_i', 'education_i')
covariates_finnish <- c("marital_status_2cat", 'smoking_2cat', 'bmi', 'cholesterol')     # education removed because missing in one finnish cohort 

POPs_finnish <- ifelse(POPs %in% c("PCB_28", "PCB_52", "OCP_PeCB", "OCP_α_HCH", "OCP_γ_HCH", 
                           "OCP_oxychlordane", "PBDE_47", "PBDE_99", "PBDE_153"), paste0(POPs, "_raw"), POPs)
proteomic <- unique(bdd_danish_proteomic$Assay) |> as.character()

# variable creation ----
## danish data ----
bdd_danish <- bdd_danish |>
  mutate(
    study = "Danish", 
    als = ifelse(!is.na(als_date), 1, 0),                                       # metadata variables creation 
    education = case_when(
      high_educ == 1 & medium_educ == 0 & low_educ == 0 ~ ">10 years of primary school", 
      high_educ == 0 & medium_educ == 1 & low_educ == 0 ~ "7-10 years of primary school", 
      high_educ == 0 & medium_educ == 0 & low_educ == 1 ~ "<7 years of primary school", 
      .default = NA),
    education = as.factor(education),
    education =  fct_relevel(education, 
                             "<7 years of primary school", 
                             "7-10 years of primary school",
                             ">10 years of primary school"),
    education_merged = fct_recode(education, 
                                  "Low" = "<7 years of primary school",
                                  "Medium" = "7-10 years of primary school",
                                  "High" = ">10 years of primary school"), 
    als_date = as.character(als_date), 
    als_date = case_when(sample == "283" ~ "2006-10-24", 
                         sample == "300" ~ "1999-07-14",
                         TRUE ~ als_date), 
    als_date = as.Date(als_date), 
    baseline_age = as.numeric(difftime(baseline_date, birth_date, units = "days")) / 365.25, 
    diagnosis_age = as.numeric(difftime(als_date, birth_date, units = "days")) / 365.25, 
    death_age = as.numeric(difftime(death_date, birth_date, units = "days")) / 365.25, 
    follow_up = as.numeric(difftime(als_date, baseline_date, units = "days"))/30.44,
    follow_up_death = case_when(als == 1 & !is.na(death_date) ~ as.numeric(difftime(death_date, als_date, units = "days")), 
                                als == 1 & is.na(death_date) ~ as.numeric(difftime("2017-09-19", als_date, units = "days")), 
                                als == 0 ~ NA), 
    follow_up_death = follow_up_death/ 30.44,
    status_death = case_when(als == 1 & !is.na(death_date) ~ 1, 
                             als == 1 & is.na(death_date) ~ 0, 
                             als == 0 ~ NA),
    # status_death = as.factor(as.character(status_death)),                     # do not factor the status_death variable (if yes, there is an issue in the cox models)
    time_baseline_diagnosis = diagnosis_age - baseline_age, 
    time_baseline_death = death_age - baseline_age, 
    time_diagnosis_death = death_age - diagnosis_age, 
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
    PCB_NDL = PCB_28 + PCB_52 + PCB_74 + PCB_99 + PCB_101 + PCB_138 + PCB_153 + PCB_170 + PCB_180 + PCB_183 + PCB_187, 
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
  }, .names = "{.col}_outlier")) |>
  mutate(
    acids_sat = myristic_acid_sat + pentadecylic_acid_sat + palmitic_acid_sat + margaric_acid_sat + stearic_acid_sat +  
      arachidic_acid_sat + dehenic_acid_sat + lignoceric_acid_sat,  
    pufas_ω9 = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9, 
    pufas_ω7 = palmitoleic_acid_ω7 + i17_acid_ω7 + ai17_acid_ω7 + vaccenic_acid_ω7,   
    pufas_ω6 = rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6, 
    pufas_ω3 = α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3, 
    pufas = pufas_ω9 + pufas_ω7 + pufas_ω6 + pufas_ω3, 
    fatty_acids = acids_sat + pufas, 
    ratio_ω6_ω3 = pufas_ω6/pufas_ω3) |>
  mutate(across(all_of(fattyacids),
                ~as.numeric(scale(.x)),
                .names = "{.col}_sd")) |>
  mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart"))

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


## finnish data ----
bdd_finnish <- bdd_finnish |>
  mutate(
    smoking_2cat = fct_recode(smoking,                                          # metadata creation
                              "Ever" = "Current", 
                              "Ever" = "Previous"), 
    education_merged = fct_recode(education, 
                                  "Low" = "<7 years",
                                  "Medium" = "7-12 years",
                                  "High" = ">12 years"), 
    marital_status_2cat = 
      fct_recode(marital_status, 
                 "Other" = "Widowed",
                 "Other" = "Divorced",
                 "Other" = "Unmarried"), 
    time_baseline_diagnosis = diagnosis_age - baseline_age, 
    follow_up = (diagnosis_age - baseline_age)*12, 
    time_baseline_death = death_age - baseline_age, 
    time_diagnosis_death = death_age - diagnosis_age) |>
  
  mutate(across(all_of(POPs_finnish), ~ factor(ntile(.x, 4),                    # POP variables creation
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
           .names = "{.col}_quart")) |>
  mutate(across(all_of(POPs_group), ~ {
    Q1 <- quantile(., 0.25, na.rm = TRUE)
    Q3 <- quantile(., 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    ifelse(. < lower_bound | . > upper_bound, NA, .)
  }, .names = "{.col}_outlier")) |>
  
  mutate(                                                                       # fatty acids variables 
    acids_sat = myristic_acid_sat + pentadecylic_acid_sat + palmitic_acid_sat + margaric_acid_sat + stearic_acid_sat +  
      arachidic_acid_sat + dehenic_acid_sat + lignoceric_acid_sat,  
    pufas_ω9 = dma16x_acid_ω9 + ptol2_acid_ω9 + dma18x_acid_ω9 + oleic_acid_ω9 + gadoleic_acid_ω9 + nervonic_acid_ω9, 
    pufas_ω7 = palmitoleic_acid_ω7 + i17_acid_ω7 + ai17_acid_ω7 + vaccenic_acid_ω7,   
    pufas_ω6 = rumenic_acid_ω6 + linoleic_acid_ω6 + dihomo_γ_linolenic_acid_ω6 + arachidonic_acid_ω6 + adrenic_acid_ω6, 
    pufas_ω3 = α_linolenic_acid_ω3 + timnodonic_acid_ω3 + clupanodonic_acid_ω3 + cervonic_acid_ω3, 
    pufas = pufas_ω9 + pufas_ω7 + pufas_ω6 + pufas_ω3, 
    fatty_acids = acids_sat + pufas, 
    ratio_ω6_ω3 = pufas_ω6/pufas_ω3) |>
  mutate(across(
    all_of(fattyacids),
    ~as.numeric(scale(.x)),
    .names = "{.col}_sd"))|>
  mutate(across(all_of(fattyacids), ~ factor(ntile(.x, 4),                           
                                             labels = c("Q1", "Q2", "Q3", "Q4")),
                .names = "{.col}_quart"))
    
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


# merged dataset ----
bdd_danish_red <- bdd_danish |> 
  select(sample, als, study, match, 
         all_of(covariates_danish), education_merged, 
         baseline_age, death_age, diagnosis_age,  
         time_baseline_diagnosis, time_baseline_death, time_diagnosis_death,
        follow_up,  follow_up_death, status_death, 
         alcohol, smoking, marital_status_2cat, blod_sys, blod_dias, 
         all_of(POPs), 
         all_of(POPs_group), ΣHCH, 
         all_of(POPs_group_quart), 
         all_of(POPs_group_quart_med),
         all_of(fattyacids), 
         all_of(explanatory), 
         all_of(explanatory_quart)) |>
  rename(smoking_2cat = smoking_2cat_i, 
         cholesterol = cholesterol_i, 
         education = education_i)

bdd_finnish_red <- bdd_finnish |> 
  select(sample, als, study, match, 
        baseline_age, sex,  smoking, bmi, cholesterol, marital_status, education, education_merged, alcohol, smoking_2cat, marital_status_2cat, blod_sys, blod_dias, 
         baseline_age, death_age, diagnosis_age, S_Ca,  
        follow_up, follow_up_death, status_death, 
        municipality, level_urbanization, thawed, 
        time_baseline_diagnosis, time_baseline_death, time_diagnosis_death,
         all_of(POPs_finnish), 
         all_of(POPs_group), ΣHCH, 
         all_of(POPs_group_quart), 
         all_of(POPs_group_quart_med),
         all_of(fattyacids), 
         all_of(explanatory), 
        all_of(explanatory_quart)) |>
  rename_with(~ gsub("_raw", "", .x)) |>
  mutate(sample = as.character(sample))

bdd <- bind_rows(bdd_danish_red, bdd_finnish_red)
rm(bdd_danish_red, bdd_finnish_red)

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
  education_merged = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  cholesterol_i = "Serum cholesterol (mmol/L)",
  blod_sys = "Systolic blood presure (mmHg)",
  blod_dias = "Diastolic blood presure (mmHg)",
  baseline_age = "Age at baseline (years)",
  diagnosis_age = "Age at ALS diagnosis (years)", 
  death_age = "Age at death (years)",
  status_death = "Status at the end of follow-up",
  time_baseline_diagnosis = "Duration between baseline and ALS diagnosis (years)", 
  time_baseline_death = "Duration between baseline and death (years)", 
  time_diagnosis_death = "Duration between diagnosis and death (years)",
  follow_up	= "Length of follow-up from baseline to ALS diagnosis (months)", 
  follow_up_death	= "Length of follow-up from ALS diagnosis (months)", 
  status_death = "Status at end of the follow-up",
  OCP_PeCB = "Pentachlorobenzene (PeCB)",            
  OCP_HCB = "HCB",            
  OCP_α_HCH = "α-HCH",
  OCP_β_HCH = "β-HCH", 
  OCP_γ_HCH = "γ-HCH", 
  OCP_oxychlordane = "Oxychlordane", 
  OCP_transnonachlor = "Transnonachlor", 
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
  PBDE_47 = "PBDE-47", 
  PBDE_99 = "PBDE-99", 
  PBDE_153 = "PBDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "Most prevalent PCBs",
  "myristic_acid_sat" =  "Myristic acid (%)", 
  "pentadecylic_acid_sat" = "Pentadecylic acid (%)", 
  "palmitic_acid_sat" =  "Palmitic acid (%)", 
  "margaric_acid_sat" =  "Margaric acid (%)", 
  "stearic_acid_sat" = "Stearic acid (%)",           
  "arachidic_acid_sat" = "Arachidic acid (%)",
  "dehenic_acid_sat" = "Dehenic acid (%)", 
  "lignoceric_acid_sat" = "Lignoceric acid (%)", 
  "dma16x_acid_ω9" = "dma16x acid ω9 (%)", 
  "ptol2_acid_ω9"  = "ptol2 acid ω9 (%)", 
  "dma18x_acid_ω9" = "dma18x acid ω9 (%)", 
  "oleic_acid_ω9" = "Oleic acid ω9 (%)", 
  "gadoleic_acid_ω9" = "Gadoleic acid ω9 (%)", 
  "nervonic_acid_ω9" = "Nervonic acid ω9 (%)",
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω7 (%)", 
  "i17_acid_ω7" = "i17 acid ω7 (%)", 
  "ai17_acid_ω7" = "ai17 acid ω7 (%)", 
  "vaccenic_acid_ω7" = "Vaccenic acid ω7 (%)", 
  "rumenic_acid_ω6" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9" = "ω9 unsaturated acids (%)", 
  "pufas_ω7" = "ω7 unsaturated acids (%)", 
  "pufas_ω6" = "ω6 unsaturated acids (%)", 
  "pufas_ω3" = "ω3 unsaturated acids (%)", 
  "pufas" = "Unsaturated acids (%)", 
  "ratio_ω6_ω3" = "ω6/ω3 ratio (%)", 
  
  "rumenic_acid_ω6_sd" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6_sd" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6_sd" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6_sd" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6_sd" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3_sd" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3_sd" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3_sd"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3_sd" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9_sd" = "ω9 unsaturated acids (%)", 
  "pufas_ω7_sd" = "ω7 unsaturated acids (%)", 
  "pufas_ω6_sd" = "ω6 unsaturated acids (%)", 
  "pufas_ω3_sd" = "ω3 unsaturated acids (%)", 
  "pufas_sd" = "Unsaturated acids (%)"
  )

var_label(bdd_finnish) <- list(
  sample = "Identifcation", 
  match = "match", 
  sex = "Sex", 
  marital_status = "Marital status",
  marital_status_2cat = "Marital status",
  municipality = "Municipality",
  level_urbanization = "Level of urbanization",
  smoking = "Smoking status", 
  smoking_2cat = "Smoking status", 
  alcohol = "Alcohol consumption (g/week)", 
  education = "Education", 
  education_merged = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  baseline_age = "Age at baseline",
  diagnosis_age = "Age at ALS diagnosis", 
  death_age = "Age at death",
  time_baseline_diagnosis = "Duration between baseline and ALS diagnosis (years)", 
  time_baseline_death = "Duration between baseline and death (years)", 
  time_diagnosis_death = "Duration between diagnosis and death (years)",
  follow_up	= "Length of follow-up from baseline to ALS diagnosis (months)", 
  follow_up_death	= "Length of follow-up from ALS diagnosis (months)", 
  status_death = "Status at end of the follow-up",
  OCP_PeCB = "Pentachlorobenzene (PeCB)",            
  OCP_HCB = "HCB",            
  OCP_α_HCH = "α-HCH",
  OCP_β_HCH = "β-HCH", 
  OCP_γ_HCH = "γ-HCH", 
  OCP_oxychlordane = "Oxychlordane", 
  OCP_transnonachlor = "Transnonachlor", 
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
  PBDE_47 = "PBDE-47", 
  PBDE_99 = "PBDE-99", 
  PBDE_153 = "PBDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "Most prevalent PCBs",
  "myristic_acid_sat" =  "Myristic acid (%)", 
  "pentadecylic_acid_sat" = "Pentadecylic acid (%)", 
  "palmitic_acid_sat" =  "Palmitic acid (%)", 
  "margaric_acid_sat" =  "Margaric acid (%)", 
  "stearic_acid_sat" = "Stearic acid (%)",           
  "arachidic_acid_sat" = "Arachidic acid (%)",
  "dehenic_acid_sat" = "Dehenic acid (%)", 
  "lignoceric_acid_sat" = "Lignoceric acid (%)", 
  "dma16x_acid_ω9" = "dma16x acid ω9 (%)", 
  "ptol2_acid_ω9"  = "ptol2 acid ω9 (%)", 
  "dma18x_acid_ω9" = "dma18x acid ω9 (%)", 
  "oleic_acid_ω9" = "Oleic acid ω9 (%)", 
  "gadoleic_acid_ω9" = "Gadoleic acid ω9 (%)", 
  "nervonic_acid_ω9" = "Nervonic acid ω9 (%)",
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω7 (%)", 
  "i17_acid_ω7" = "i17 acid ω7 (%)", 
  "ai17_acid_ω7" = "ai17 acid ω7 (%)", 
  "vaccenic_acid_ω7" = "Vaccenic acid ω7 (%)", 
  "rumenic_acid_ω6" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9" = "ω9 unsaturated acids (%)", 
  "pufas_ω7" = "ω7 unsaturated acids (%)", 
  "pufas_ω6" = "ω6 unsaturated acids (%)", 
  "pufas_ω3" = "ω3 unsaturated acids (%)", 
  "pufas" = "Unsaturated acids (%)", 
  "ratio_ω6_ω3" = "ω6/ω3 ratio (%)", 
  
  "rumenic_acid_ω6_sd" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6_sd" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6_sd" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6_sd" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6_sd" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3_sd" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3_sd" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3_sd"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3_sd" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9_sd" = "ω9 unsaturated acids (%)", 
  "pufas_ω7_sd" = "ω7 unsaturated acids (%)", 
  "pufas_ω6_sd" = "ω6 unsaturated acids (%)", 
  "pufas_ω3_sd" = "ω3 unsaturated acids (%)", 
  "pufas_sd" = "Unsaturated acids (%)", 
  "ratio_ω6_ω3_sd" = "ω6/ω3 ratio (%)")


var_label(bdd) <- list(
  sample = "Identifcation", 
  als = "Amyotrophic lateral sclerosis",
  study = "Cohort",
  match = "match", 
  sex = "Sex", 
  baseline_age = "Age at baseline",
  marital_status = "Marital status",
  marital_status_2cat = "Marital status", 
  municipality = "Municipality",
  level_urbanization = "Level of urbanization",
  smoking = "Smoking status", 
  smoking_2cat = "Smoking status", 
  alcohol = "Alcohol consumption (g/week)", 
  education = "Education", 
  education_merged = "Education", 
  bmi = "Boby mass index (kg/m²)",
  cholesterol = "Serum cholesterol (mmol/L)",
  blod_sys = "Systolic blood presure (mmHg)",
  blod_dias = "Diastolic blood presure (mmHg)",
  diagnosis_age = "Age at ALS diagnosis", 
  death_age = "Age at death",
  time_baseline_diagnosis = "Duration between baseline and ALS diagnosis (years)", 
  time_baseline_death = "Duration between baseline and death (years)", 
  time_diagnosis_death = "Duration between diagnosis and death (years)",
  follow_up	= "Length of follow-up from baseline to ALS diagnosis (months)", 
  follow_up_death	= "Length of follow-up from ALS diagnosis (months)", 
  status_death = "Status at end of the follow-up",
  OCP_PeCB = "Pentachlorobenzene (PeCB)",            
  OCP_HCB = "HCB",            
  OCP_α_HCH = "α-HCH",
  OCP_β_HCH = "β-HCH", 
  OCP_γ_HCH = "γ-HCH", 
  OCP_oxychlordane = "Oxychlordane", 
  OCP_transnonachlor = "Transnonachlor", 
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
  PBDE_47 = "PBDE-47", 
  PBDE_99 = "PBDE-99", 
  PBDE_153 = "PBDE-153", 
  PCB_DL =  "Dioxin-like PCBs", 
  PCB_NDL = "Non dioxin-like PCBs", 
  PCB_4 = "Most prevalent PCBs",
  "myristic_acid_sat" =  "Myristic acid (%)", 
  "pentadecylic_acid_sat" = "Pentadecylic acid (%)", 
  "palmitic_acid_sat" =  "Palmitic acid (%)", 
  "margaric_acid_sat" =  "Margaric acid (%)", 
  "stearic_acid_sat" = "Stearic acid (%)",           
  "arachidic_acid_sat" = "Arachidic acid (%)",
  "dehenic_acid_sat" = "Dehenic acid (%)", 
  "lignoceric_acid_sat" = "Lignoceric acid (%)", 
  "dma16x_acid_ω9" = "dma16x acid ω9 (%)", 
  "ptol2_acid_ω9"  = "ptol2 acid ω9 (%)", 
  "dma18x_acid_ω9" = "dma18x acid ω9 (%)", 
  "oleic_acid_ω9" = "Oleic acid ω9 (%)", 
  "gadoleic_acid_ω9" = "Gadoleic acid ω9 (%)", 
  "nervonic_acid_ω9" = "Nervonic acid ω9 (%)",
  "palmitoleic_acid_ω7" = "Palmitoleic acid ω7 (%)", 
  "i17_acid_ω7" = "i17 acid ω7 (%)", 
  "ai17_acid_ω7" = "ai17 acid ω7 (%)", 
  "vaccenic_acid_ω7" = "Vaccenic acid ω7 (%)", 
  "rumenic_acid_ω6" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9" = "ω9 unsaturated acids (%)", 
  "pufas_ω7" = "ω7 unsaturated acids (%)", 
  "pufas_ω6" = "ω6 unsaturated acids (%)", 
  "pufas_ω3" = "ω3 unsaturated acids (%)", 
  "pufas" = "Unsaturated acids (%)", 
  "ratio_ω6_ω3" = "ω6/ω3 ratio (%)", 
  
  "rumenic_acid_ω6_sd" = "Rumenic acid ω6 (%)", 
  "linoleic_acid_ω6_sd" = "Linoleic acid ω6 (%)",          
  "dihomo_γ_linolenic_acid_ω6_sd" = "Dihomo-γ-linolenic acid ω6 (%)", 
  "arachidonic_acid_ω6_sd" = "Arachidonic acid ω6 (%)", 
  "adrenic_acid_ω6_sd" =  "Adrenic acid ω6 (%)", 
  "α_linolenic_acid_ω3_sd" = "α-linolenic acid (ALA) ω3 (%)", 
  "timnodonic_acid_ω3_sd" = "Timnodonic acid (EPA) ω3 (%)", 
  "clupanodonic_acid_ω3_sd"  ="Clupanodonic acid (DPA) ω3 (%)", 
  "cervonic_acid_ω3_sd" = "Cervonic acid (DHA) ω3 (%)", 
  "pufas_ω9_sd" = "ω9 unsaturated acids (%)", 
  "pufas_ω7_sd" = "ω7 unsaturated acids (%)", 
  "pufas_ω6_sd" = "ω6 unsaturated acids (%)", 
  "pufas_ω3_sd" = "ω3 unsaturated acids (%)", 
  "pufas_sd" = "Unsaturated acids (%)", 
  "ratio_ω6_ω3_sd" = "ω6/ω3 ratio (%)")


fattyacids_labels <- c(
  "Unsaturated acids (%)" = "pufas", 
  "ω9 unsaturated acids (%)" = "pufas_ω9", 
  "ω7 unsaturated acids (%)" = "pufas_ω7", 
  "ω6 unsaturated acids (%)" = "pufas_ω6", 
  "ω3 unsaturated acids (%)" = "pufas_ω3", 
  "ω6/ω3 ratio" = "ratio_ω6_ω3",
  "Myristic acid (%)" = "myristic_acid_sat", 
  "Pentadecylic acid (%)" = "pentadecylic_acid_sat", 
  "Palmitic acid (%)" = "palmitic_acid_sat", 
  "Margaric acid (%)" = "margaric_acid_sat", 
  "Stearic acid (%)" = "stearic_acid_sat",           
  "Arachidic acid (%)" = "arachidic_acid_sat",
  "Dehenic acid (%)" = "dehenic_acid_sat", 
  "Lignoceric acid (%)" = "lignoceric_acid_sat", 
  "dma16x acid ω9 (%)" = "dma16x_acid_ω9", 
  "ptol2 acid ω9 (%)"  = "ptol2_acid_ω9", 
  "dma18x acid ω9 (%)" = "dma18x_acid_ω9", 
  "Oleic acid ω9 (%)" = "oleic_acid_ω9", 
  "Gadoleic acid ω9 (%)" = "gadoleic_acid_ω9", 
  "Nervonic acid ω9 (%)" = "nervonic_acid_ω9",
  "Palmitoleic acid ω7 (%)" = "palmitoleic_acid_ω7", 
  "i17 acid ω7 (%)" = "i17_acid_ω7", 
  "ai17 acid ω7 (%)" = "ai17_acid_ω7", 
  "Vaccenic acid ω7 (%)" = "vaccenic_acid_ω7", 
  "Rumenic acid ω6 (%)" = "rumenic_acid_ω6", 
  "Linoleic acid ω6 (%)" = "linoleic_acid_ω6",          
  "Dihomo-γ-linolenic acid ω6 (%)" = "dihomo_γ_linolenic_acid_ω6", 
  "Arachidonic acid ω6 (%)" = "arachidonic_acid_ω6", 
  "Adrenic acid ω6 (%)" = "adrenic_acid_ω6", 
  "α-linolenic acid (ALA) ω3 (%)" = "α_linolenic_acid_ω3", 
  "Timnodonic acid (EPA) ω3 (%)" = "timnodonic_acid_ω3", 
  "Clupanodonic acid (DPA) ω3 (%)" = "clupanodonic_acid_ω3", 
  "Cervonic acid (DHA) ω3 (%)" = "cervonic_acid_ω3") 

explanatory_labels <- c(  
  "Unsaturated acids (%)" = "pufas", 
  "ω9 unsaturated acids (%)" = "pufas_ω9", 
  "ω7 unsaturated acids (%)" = "pufas_ω7", 
  "ω6 unsaturated acids (%)" = "pufas_ω6", 
  "ω3 unsaturated acids (%)" = "pufas_ω3", 
  "ω6/ω3 ratio" = "ratio_ω6_ω3",
  "Rumenic acid ω6 (%)" = "rumenic_acid_ω6", 
  "Linoleic acid ω6 (%)" = "linoleic_acid_ω6",          
  "Dihomo-γ-linolenic acid ω6 (%)" = "dihomo_γ_linolenic_acid_ω6", 
  "Arachidonic acid ω6 (%)" = "arachidonic_acid_ω6", 
  "Adrenic acid ω6 (%)" = "adrenic_acid_ω6", 
  "α-linolenic acid (ALA) ω3 (%)" = "α_linolenic_acid_ω3", 
  "Timnodonic acid (EPA) ω3 (%)" = "timnodonic_acid_ω3", 
  "Clupanodonic acid (DPA) ω3 (%)" = "clupanodonic_acid_ω3", 
  "Cervonic acid (DHA) ω3 (%)" = "cervonic_acid_ω3")

explanatory_sd_labels <- c(  
  "Unsaturated acids (%)" = "pufas_sd", 
  "ω9 unsaturated acids (%)" = "pufas_ω9_sd", 
  "ω7 unsaturated acids (%)" = "pufas_ω7_sd", 
  "ω6 unsaturated acids (%)" = "pufas_ω6_sd", 
  "ω3 unsaturated acids (%)" = "pufas_ω3_sd", 
  "ω6/ω3 ratio" = "ratio_ω6_ω3_sd",
  "Rumenic acid ω6 (%)" = "rumenic_acid_ω6_sd", 
  "Linoleic acid ω6 (%)" = "linoleic_acid_ω6_sd",          
  "Dihomo-γ-linolenic acid ω6 (%)" = "dihomo_γ_linolenic_acid_ω6_sd", 
  "Arachidonic acid ω6 (%)" = "arachidonic_acid_ω6_sd", 
  "Adrenic acid ω6 (%)" = "adrenic_acid_ω6_sd", 
  "α-linolenic acid (ALA) ω3 (%)" = "α_linolenic_acid_ω3_sd", 
  "Timnodonic acid (EPA) ω3 (%)" = "timnodonic_acid_ω3_sd", 
  "Clupanodonic acid (DPA) ω3 (%)" = "clupanodonic_acid_ω3_sd", 
  "Cervonic acid (DHA) ω3 (%)" = "cervonic_acid_ω3_sd")


explanatory_quart_labels <- c(
  pufas_quart = "unsaturated acids",
  pufas_ω9_quart = "ω9 unsaturated acids",
  pufas_ω7_quart = "ω7 unsaturated acids",
  pufas_ω6_quart = "ω6 unsaturated acids",
  pufas_ω3_quart = "ω3 unsaturated acids",
  ratio_ω6_ω3_quart = "ω6/ω3 ratio",
  rumenic_acid_ω6_quart = "rumenic acid ω6",
  linoleic_acid_ω6_quart = "linoleic acid ω6",
  dihomo_γ_linolenic_acid_ω6_quart = "dihomo-γ-linolenic acid ω6",
  arachidonic_acid_ω6_quart = "arachidonic acid ω6",
  adrenic_acid_ω6_quart = "adrenic acid ω6",
  α_linolenic_acid_ω3_quart = "α-linolenic acid (ALA) ω3",
  timnodonic_acid_ω3_quart = "timnodonic acid (EPA) ω3",
  clupanodonic_acid_ω3_quart = "clupanodonic acid (DPA) ω3",
  cervonic_acid_ω3_quart = "cervonic acid (DHA) ω3"
)

POPs_labels <- c(
  "Most prevalent PCBs" = "PCB_4",
  "Dioxin-like PCBs" = "PCB_DL",
  "PCB-118" = "PCB_118",
  "PCB-156" = "PCB_156",
  "Non dioxin-like PCBs" = "PCB_NDL",
  "PCB-28" = "PCB_28",
  "PCB-52" = "PCB_52",
  "PCB-74" = "PCB_74",
  "PCB-99" = "PCB_99",
  "PCB-101" = "PCB_101",
  "PCB-138" = "PCB_138",
  "PCB-153" = "PCB_153",
  "PCB-170" = "PCB_170",
  "PCB-180" = "PCB_180",
  "PCB-183" = "PCB_183",
  "PCB-187" = "PCB_187",
  "HCB" = "OCP_HCB",
  "ΣDDT" = "ΣDDT",
  "p,p'-DDT" = "OCP_pp_DDT",
  "p,p'-DDE" = "OCP_pp_DDE",
  "α-HCH" = "OCP_α_HCH",
  "β-HCH" = "OCP_β_HCH",
  "γ-HCH" = "OCP_γ_HCH",
  "Σchlordane" = "Σchlordane",
  "Oxychlordane" = "OCP_oxychlordane",
  "Transnonachlor" = "OCP_transnonachlor",
  "Pentachlorobenzene (PeCB)" = "OCP_PeCB",
  "ΣPBDE" = "ΣPBDE",
  "PBDE-47" = "PBDE_47",
  "PBDE-99" = "PBDE_99",
  "PBDE-153" = "PBDE_153")


POPs_group_labels <- c(
  "Most prevalent PCBs" = "PCB_4",
  "Dioxin-like PCBs" = "PCB_DL",
  "Non dioxin-like PCBs" = "PCB_NDL",
  "HCB" = "OCP_HCB",
  "ΣDDT" = "ΣDDT",
  "β-HCH" = "OCP_β_HCH",
  "Σchlordane" = "Σchlordane",
  "ΣPBDE" = "ΣPBDE")

POPs_group_labels_finnish <- c(
  "Most prevalent PCBs" = "PCB_4",
  "Dioxin-like PCBs" = "PCB_DL",
  "Non dioxin-like PCBs" = "PCB_NDL",
  "HCB" = "OCP_HCB",
  "ΣDDT" = "ΣDDT",
  "ΣHCH" = "ΣHCH", 
  "β-HCH" = "OCP_β_HCH",
  "γ-HCH" = "OCP_γ_HCH",
  "Σchlordane" = "Σchlordane",
  "PeCB" = "OCP_PeCB")

POPs_group_sd_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group_sd)

POPs_group_quart_labels <- c(
  "Most prevalent PCBs" = "PCB_4_quart",
  "Dioxin-like PCBs" = "PCB_DL_quart",
  "Non dioxin-like PCBs" = "PCB_NDL_quart",
  "HCB" = "OCP_HCB_quart",
  "ΣDDT" = "ΣDDT_quart",
  "β-HCH" = "OCP_β_HCH_quart",
  "Σchlordane" = "Σchlordane_quart",
  "ΣPBDE" = "ΣPBDE_quart")
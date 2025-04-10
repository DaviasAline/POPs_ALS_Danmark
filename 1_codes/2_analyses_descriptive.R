# Aline Davias 
# April 9, 2025

source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# danish data ----
## metadata ----
descrip_covar_danish <- bdd_danish|> 
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"))|>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status, education, alcohol, smoking, bmi, cholesterol)|>
  tbl_summary(by = als)|>
  bold_labels()|>
  add_overall() 

## POPs ----
descrip_expo_danish <- bdd_danish |>
  select(all_of(POPs), all_of(POPs_group)) |>
  pivot_longer(cols = everything(), names_to = "POPs", values_to = "values") |>
  mutate(
    POPs = fct_recode(POPs, 
      "PCB-118,138,153,180" = "PCB_4", 
      "Dioxin-like PCBs" = "PCB_DL",
      "Non-dioxin-like PCBs" = "PCB_NDL", 
      "p,p’-DDT"  = "OCP_pp_DDT", 
      "p,p’-DDE" = "OCP_pp_DDE", 
      "Oxychlordane" = "OCP_oxychlordane", 
      "Transnonachlor" = "OCP_transnonachlor"), 
    POPs = gsub("OCP_", "", POPs), 
    POPs = gsub("_", "-", POPs), 
  POPs = fct_relevel(POPs, 
                       "PCB-118,138,153,180", 
                     "Dioxin-like PCBs", "PCB-118", "PCB-156",
                     "Non-dioxin-like PCBs",  "PCB-28",  "PCB-52", "PCB-74", "PCB-99", 
                     "PCB-101",   "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183", "PCB-187", 
                     "HCB","ΣDDT","p,p’-DDT", "p,p’-DDE", 
                     "α-HCH", "β-HCH", "γ-HCH", 
                     "Σchlordane", "Oxychlordane", "Transnonachlor", 
                     "PeCB",   "ΣPBDE", "PBDE-47","PBDE-99","PBDE-153")) |>
  ggplot() +
  aes(x = POPs, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)") +
  coord_flip() +
  theme_lucid()

descrip_expo_danish_by_als <- bdd_danish |>
  select(als, all_of(POPs), all_of(POPs_group)) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(
    POPs = fct_recode(POPs, 
                      "PCB-118,138,153,180" = "PCB_4", 
                      "Dioxin-like PCBs" = "PCB_DL",
                      "Non-dioxin-like PCBs" = "PCB_NDL", 
                      "p,p’-DDT"  = "OCP_pp_DDT", 
                      "p,p’-DDE" = "OCP_pp_DDE", 
                      "Oxychlordane" = "OCP_oxychlordane", 
                      "Transnonachlor" = "OCP_transnonachlor"), 
    POPs = gsub("OCP_", "", POPs), 
    POPs = gsub("_", "-", POPs), 
    POPs = fct_relevel(POPs, 
                       "PCB-118,138,153,180", 
                       "Dioxin-like PCBs", "PCB-118", "PCB-156",
                       "Non-dioxin-like PCBs",  "PCB-28",  "PCB-52", "PCB-74", "PCB-99", 
                       "PCB-101",   "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183", "PCB-187", 
                       "HCB","ΣDDT","p,p’-DDT", "p,p’-DDE", 
                       "α-HCH", "β-HCH", "γ-HCH", 
                       "Σchlordane", "Oxychlordane", "Transnonachlor", 
                       "PeCB",   "ΣPBDE", "PBDE-47","PBDE-99","PBDE-153"), 
    POPs = fct_rev(POPs), 
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)", fill = "ALS") +
  coord_flip() +
  theme_lucid()


## fatty acids ----
descrip_fattyacids_danish <- bdd_danish |>
  select(all_of(explanatory), all_of(fattyacids)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fatty_acid_labels)) |>
  ggplot() +
  aes(x = fattyacids, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fattu acids", y = "Values (log transformed)") +
  coord_flip() +
  theme_lucid()

descrip_fattyacids_danish_by_als <- bdd_danish |>
  select(als, all_of(explanatory), all_of(fattyacids)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fatty_acid_labels)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  ggplot() +
  aes(x = fattyacids, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fatty acids", y = "Values (log transformed)", fill = "ALS") +
  coord_flip() +
  theme_lucid()


# finnish data ----
## metadata ----
descrip_covar_finnish <- bdd_finnish |> 
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"))|>
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status, education, alcohol, smoking, bmi, cholesterol)|>
  tbl_summary(by = als)|>
  bold_labels()|>
  add_overall() 

## POPs ----
descrip_expo_finnish <- bdd |>
  filter(study %in% c("Finnish_1", "Finnish_2", "Finnish_3")) |>
  select(all_of(POPs), all_of(POPs_group)) |>
  pivot_longer(cols = everything(), names_to = "POPs", values_to = "values") |>
  mutate(
    POPs =   fct_relevel(POPs, 
                         "PCB_DL", "PCB_118", "PCB_156", "PCB_NDL", "PCB_4", "PCB_28",
                         "PCB_52", "PCB_74", "PCB_99", "PCB_101", "PCB_138", "PCB_153",
                         "PCB_170", "PCB_180", "PCB_183", "PCB_187", "OCP_HCB", "OCP_PeCB",
                         "ΣDDT", "OCP_pp_DDE", "OCP_pp_DDT", "Σchlordane", "OCP_oxychlordane",
                         "OCP_transnonachlor", "OCP_α_HCH", "OCP_β_HCH", "OCP_γ_HCH",
                         "ΣPBDE", "PBDE_47", "PBDE_99", "PBDE_153"), 
    POPs = fct_rev(POPs)) |>
  ggplot() +
  aes(x = POPs, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)") +
  coord_flip() +
  theme_lucid()

descrip_expo_finnish_by_als <- bdd |>
  filter(study %in% c("Finnish_1", "Finnish_2", "Finnish_3")) |>
  select(als, all_of(POPs), all_of(POPs_group)) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(
    POPs = fct_recode(POPs, 
                      "PCB-118,138,153,180" = "PCB_4", 
                      "Dioxin-like PCBs" = "PCB_DL",
                      "Non-dioxin-like PCBs" = "PCB_NDL", 
                      "p,p’-DDT"  = "OCP_pp_DDT", 
                      "p,p’-DDE" = "OCP_pp_DDE", 
                      "Oxychlordane" = "OCP_oxychlordane", 
                      "Transnonachlor" = "OCP_transnonachlor"), 
    POPs = gsub("OCP_", "", POPs), 
    POPs = gsub("_", "-", POPs), 
    POPs = fct_relevel(POPs, 
                       "PCB-118,138,153,180", 
                       "Dioxin-like PCBs", "PCB-118", "PCB-156",
                       "Non-dioxin-like PCBs",  "PCB-28",  "PCB-52", "PCB-74", "PCB-99", 
                       "PCB-101",   "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183", "PCB-187", 
                       "HCB","ΣDDT","p,p’-DDT", "p,p’-DDE", 
                       "α-HCH", "β-HCH", "γ-HCH", 
                       "Σchlordane", "Oxychlordane", "Transnonachlor", 
                       "PeCB",   "ΣPBDE", "PBDE-47","PBDE-99","PBDE-153"), 
    POPs = fct_rev(POPs), 
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

## fatty acids ----
descrip_fattyacids_finnish <- bdd_finnish |>
  select(all_of(explanatory), all_of(fattyacids)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fatty_acid_labels)) |>
  ggplot() +
  aes(x = fattyacids, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fattu acids", y = "Values (log transformed)") +
  coord_flip() +
  theme_lucid()

descrip_fattyacids_finnish_by_als <- bdd_finnish |>
  select(als, all_of(explanatory), all_of(fattyacids)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fatty_acid_labels)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  ggplot() +
  aes(x = fattyacids, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fatty acids", y = "Values (log transformed)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

# comparison danish / finnish ----
## metadata ----
## POPs ----
## fatty acids ----
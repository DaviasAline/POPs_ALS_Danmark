# Aline Davias 
# April 9, 2025

source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# danish data ----
## metadata ----
covar_danish <- bdd_danish|> 
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
POPs_table_danish <- descrip_num(data = bdd_danish, vars = POPs_tot)
POPs_table_danish <- left_join(POPs_table_danish, bdd_danish_loq, by = "variable") |>
  mutate(variable = factor(variable, levels = POPs_labels), 
         variable = fct_recode(variable, !!!POPs_labels)) |>
  arrange(variable) 

POPs_table_danish_by_als <- bdd_danish |>
  select(als, all_of(POPs_tot)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall() |>
  add_p()

POPs_boxplot_danish <- bdd_danish |>
  select(all_of(POPs_tot)) |>
  pivot_longer(cols = everything(), names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease plasma concentrations (pg/ml)") +
  coord_flip() +
  theme_lucid()

POPs_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(POPs_tot)) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs),
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease plasma concentrations (pg/ml)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

POPs_group_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(POPs_group)) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_group_labels), 
         POPs = fct_recode(POPs, !!!POPs_group_labels), 
         POPs = fct_rev(POPs),
         als = as.character(als), 
         als = fct_recode(als, 
                          "Controls" = "0",
                          "Cases" = "1")) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease plasma concentrations (pg/ml)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

POPs_group_boxplot_danish_by_death <- bdd_danish |>
  filter(als == 1) |>
  select(status_death, all_of(POPs_group)) |>
  pivot_longer(cols = -status_death, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_group_labels), 
         POPs = fct_recode(POPs, !!!POPs_group_labels), 
         POPs = fct_rev(POPs),
         status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1")) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = status_death) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease plasma concentrations (pg/ml)", fill = "") +
  coord_flip() +
  theme_lucid()


POPs_heatmap_danish <- bdd_danish |> 
  select(all_of(POPs_tot)) |>
  rename(!!!POPs_labels) 

POPs_heatmap_danish <- cor(POPs_heatmap_danish, 
                           use = "pairwise.complete.obs", 
                           method = "pearson")

## fatty acids ----
fattyacids_table_danish <- 
  descrip_num(data = bdd_danish, vars = c(fattyacids)) |> 
  mutate(variable = factor(variable, levels = fattyacids_labels), 
         variable = fct_recode(variable, !!!fattyacids_labels)) |>
  arrange(variable) 

fattyacids_table_danish_by_als <- bdd_danish |>
  select(als, all_of(fattyacids)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall() |>
  add_p() 

fattyacids_boxplot_danish <- bdd_danish |>
  select(all_of(fattyacids)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_labels), 
         fattyacids = fct_rev(fattyacids)) |>
  ggplot() +
  aes(x = fattyacids, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fatty acids", y = "Pre-disease serum proportions (%)") +
  coord_flip() +
  theme_lucid()

fattyacids_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(fattyacids)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_labels), 
         fattyacids = fct_rev(fattyacids)) |>
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
  labs(x = "Fatty acids", y = "Pre-disease serum proportions (%)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

fattyacids_heatmap_danish <- bdd_danish |> select(all_of(fattyacids)) |>
  rename(!!!fattyacids_labels) 
fattyacids_heatmap_danish <- cor(fattyacids_heatmap_danish, 
                                 use = "pairwise.complete.obs", 
                                 method = "pearson")

POPs_fattyacids_heatmap_danish <- 
  heatmap_cor_pairwise(fattyacids, POPs, decimal = 1, data = bdd_danish)

## proteomic ----
proteomic_table_danish <- descrip_num(data = bdd_danish, vars = proteomic)

proteomic_table_danish_by_als <- bdd_danish |>
  select(als, all_of(proteomic)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall()

proteomic_boxplot_danish <- bdd_danish |>
  select(all_of(proteomic)) |>
  pivot_longer(cols = everything(), names_to = "Proteomic", values_to = "values") |>
  # mutate(POPs = factor(POPs, levels = POPs_labels), 
  #        POPs = fct_recode(POPs, !!!POPs_labels), 
  #        POPs = fct_rev(POPs)) |>
  # arrange(POPs) |>
  ggplot() +
  aes(x = Proteomic, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  # scale_y_continuous(trans = "log", 
  #                    labels = number_format(accuracy = 1)) +
  labs(x = "Proteomic", y = "Pre-disease plasma concentrations (unit?)") +
  coord_flip() +
  theme_lucid()

proteomic_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(proteomic)) |>
  pivot_longer(cols = -als, names_to = "Proteomic", values_to = "values") |>
  mutate(
    # POPs = factor(POPs, levels = POPs_labels), 
    # POPs = fct_recode(POPs, !!!POPs_labels), 
    # POPs = fct_rev(POPs),
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  # arrange(POPs) |>
  ggplot() +
  aes(x = Proteomic, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  # scale_y_continuous(trans = "log", 
  #                    labels = number_format(accuracy = 1)) +
  labs(x = "Proteomic", y = "Pre-disease plasma concentrations (unit?)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

# proteomic_group_boxplot_danish_by_als <- bdd_danish |>
#   select(als, all_of(proteomic)) |>
#   pivot_longer(cols = -als, names_to = "Proteomic", values_to = "values") |>
#   mutate(
#     # POPs = factor(POPs, levels = POPs_group_labels), 
#     # POPs = fct_recode(POPs, !!!POPs_group_labels), 
#     # POPs = fct_rev(POPs),
#     als = as.character(als), 
#     als = fct_recode(als, 
#                      "Controls" = "0",
#                      "Cases" = "1")) |>
#   # arrange(POPs) |>
#   ggplot() +
#   aes(x = Proteomic, y = values, fill = als) +
#   geom_boxplot() +
#   # scale_fill_hue(direction = 1, 
#   #                guide = guide_legend(reverse = TRUE)) +
#   scale_y_continuous(trans = "log", 
#                      labels = number_format(accuracy = 1)) +
#   labs(x = "Proteomic", y = "Pre-disease plasma concentrations (unit?)", fill = "ALS") +
#   coord_flip() +
#   theme_lucid()

proteomic_boxplot_danish_by_death <- bdd_danish |>
  filter(als == 1) |>
  select(status_death, all_of(proteomic)) |>
  pivot_longer(cols = -status_death, names_to = "Proteomic", values_to = "values") |>
  mutate(
    # POPs = factor(POPs, levels = POPs_group_labels),
    # POPs = fct_recode(POPs, !!!POPs_group_labels),
    # POPs = fct_rev(POPs),
    status_death = as.character(status_death),
    status_death = fct_recode(status_death,
                              "Alive" = "0",
                              "Deceased" = "1")) |>
  # arrange(POPs) |>
  ggplot() +
  aes(x = Proteomic, y = values, fill = status_death) +
  geom_boxplot() +
  scale_fill_hue(direction = 1,
                 guide = guide_legend(reverse = TRUE)) +
  # scale_y_continuous(trans = "log",
  #                    labels = number_format(accuracy = 1)) +
  labs(x = "Proteomic", y = "Pre-disease plasma concentrations (unit?)", fill = "") +
  coord_flip() +
  theme_lucid()


proteomic_heatmap_danish <- bdd_danish |> 
  select(all_of(proteomic)) 
# |> rename(!!!POPs_labels) 

proteomic_heatmap_danish <- cor(proteomic_heatmap_danish, 
                                use = "pairwise.complete.obs", 
                                method = "pearson")


# finnish data ----
## metadata ----
covar_finnish <- bdd_finnish |> 
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
POPs_table_finnish <- bdd |> filter(!study %in% "Danish") 
POPs_table_finnish <- descrip_num(data = POPs_table_finnish, vars = POPs_tot) |> 
  mutate(variable = factor(variable, levels = POPs_labels), 
         variable = fct_recode(variable, !!!POPs_labels)) |>
  arrange(variable) 

POPs_table_finnish_by_als <- bdd |>
  filter(!study %in% c("Danish")) |>
  select(als, all_of(POPs_tot)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall() |>
  add_p(include = -`OCP_α_HCH`) 

POPs_boxplot_finnish <- bdd |>
  filter(study %in% c("FMC", "FMCF", "MFH")) |>
  select(all_of(POPs_tot)) |>
  pivot_longer(cols = everything(), names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease serum concentrations (pg/ml)") +
  coord_flip() +
  theme_lucid()

POPs_boxplot_finnish_by_als <- bdd |>
  filter(study %in% c("FMC", "FMCF", "MFH")) |>
  select(als, all_of(POPs), all_of(POPs_group)) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs),
    als = as.character(als), 
    als = fct_recode(als, 
                     "Controls" = "0",
                     "Cases" = "1")) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease serum concentrations (pg/ml)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

POPs_group_boxplot_finnish_by_als <- bdd |>
  filter(study %in% c("FMC", "FMCF", "MFH")) |>
  select(als, all_of(POPs_group)) |>
  select(!contains("PBDE")) |>
  pivot_longer(cols = -als, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_group_labels), 
         POPs = fct_recode(POPs, !!!POPs_group_labels), 
         POPs = fct_rev(POPs),
         als = as.character(als), 
         als = fct_recode(als, 
                          "Controls" = "0",
                          "Cases" = "1")) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease serum concentrations (pg/ml)", fill = "ALS") +
  coord_flip() +
  theme_lucid()



POPs_heatmap_finnish <- bdd |>
  filter(!study %in% c("Danish")) |>
  select(all_of(POPs_tot)) |>
  rename(!!!POPs_labels) |>
  select(-"α-HCH")

POPs_heatmap_finnish <- cor(POPs_heatmap_finnish, 
                            use = "pairwise.complete.obs", 
                            method = "pearson")

## fatty acids ----
fattyacids_table_finnish <-
  descrip_num(data = bdd_finnish, vars = fattyacids) |> 
  mutate(variable = fct_recode(variable, !!!fattyacids_labels)) 

fattyacids_table_finnish_by_als <- 
  bdd_finnish |>
  select(als, all_of(fattyacids)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall() |>
  add_p() 

fattyacids_boxplot_finnish <- bdd_finnish |>
  select(all_of(fattyacids)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_labels), 
         fattyacids = fct_rev(fattyacids)) |>
  ggplot() +
  aes(x = fattyacids, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fatty acids", y = "Pre-disease serum proportions (%)") +
  coord_flip() +
  theme_lucid()

fattyacids_boxplot_finnish_by_als <- bdd_finnish |>
  select(als, all_of(fattyacids)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_labels), 
         fattyacids = fct_rev(fattyacids), 
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
  labs(x = "Fatty acids", y = "Pre-disease serum proportions (%)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

fattyacids_heatmap_finnish <- bdd_finnish |> 
  select(all_of(fattyacids)) |>
  rename(!!!fattyacids_labels) 
fattyacids_heatmap_finnish <- cor(fattyacids_heatmap_finnish, 
                                  use = "pairwise.complete.obs", 
                                  method = "pearson")

POPs_fattyacids_heatmap_finnish <- heatmap_cor_pairwise(fattyacids, POPs_finnish, decimal = 1, data = bdd_finnish)

# comparison danish / finnish on the total population ----
## metadata ----
covar_comp <- bdd |>
  select("study", "sex", "marital_status", "smoking", "alcohol", "education", 
         "bmi", "cholesterol", "blod_sys", "blod_dias", "baseline_age", "diagnosis_age") |>
  tbl_summary(by = "study") |>
  bold_labels() |>
  add_p(include = -education)

## POPs ----
POPs_table_comp <- bdd |>
  filter(als == 0) |>
  select(study, all_of(POPs_tot)) |>
  tbl_summary(by = "study") |>
  bold_labels() |>
  add_n() 

POPs_boxplot_comp <- bdd |>
  select(study, all_of(POPs_tot)) |>
  pivot_longer(cols = -study, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_labels), 
         POPs = fct_recode(POPs, !!!POPs_labels), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = study) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease concentrations (pg/ml)", fill = "Study") +
  coord_flip() +
  theme_lucid()

# bdd |>
#   select(study, all_of(POPs), all_of(POPs_group)) |>
#   pivot_longer(cols = -study, names_to = "POPs", values_to = "values") |>
#     mutate(
#     POPs =   fct_relevel(POPs, 
#     "PCB_DL", "PCB_118", "PCB_156", "PCB_NDL", "PCB_4", "PCB_28",
#     "PCB_52", "PCB_74", "PCB_99", "PCB_101", "PCB_138", "PCB_153",
#     "PCB_170", "PCB_180", "PCB_183", "PCB_187", "OCP_HCB", "OCP_PeCB",
#     "ΣDDT", "OCP_pp_DDE", "OCP_pp_DDT", "Σchlordane", "OCP_oxychlordane",
#     "OCP_transnonachlor", "OCP_α_HCH", "OCP_β_HCH", "OCP_γ_HCH",
#     "ΣPBDE", "PBDE_47", "PBDE_99", "PBDE_153"), 
#     POPs = fct_rev(POPs)) |>
#   filter(POPs %in% c("OCP_HCB", "PCB_4", "PCB_DL", "PCB_NDL", "OCP_PeCB", "OCP_β_HCH", "Σchlordane", "ΣDDT", "ΣHCH", "ΣPBDE")) |>
#   ggplot() +
#   aes(x = POPs, y = values, fill = study) +
#   geom_boxplot() +
#   scale_fill_hue(direction = 1) +
#   scale_y_continuous(trans = "log", 
#                      labels = number_format(accuracy = 1)) +
#   labs(x = "POPs", y = "Values (pg/ml, log transformed)", fill = "Study") +
#   coord_flip() +
#   theme_lucid() + 
#   scale_fill_discrete()  


## fatty acids ----
fattyacids_table_comp <- bdd |>
  filter(als == 0) |>
  select(study, all_of(fattyacids)) |>
  tbl_summary(by = "study") |>
  bold_labels() |>
  add_n()

fattyacids_boxplot_comp <- bdd |>
  select(study, all_of(fattyacids)) |>
  pivot_longer(cols = -study, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fattyacids_labels)) |>
  ggplot() +
  aes(x = fattyacids, y = values, fill = study) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "Fatty acids", y = "Pre-disease proportions (%)", fill = "Cohort") +
  coord_flip() +
  theme_lucid()

# comparison danish / finnish on the cases only ----
## metadata ----
covar_comp_cases <- bdd |> 
  filter(als == 1) |>
  mutate(status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1")) |>
  select(
    study, baseline_age, diagnosis_age, death_age, 
    follow_up_death, status_death,
    sex, marital_status_2cat, education_merged, alcohol, smoking, bmi, cholesterol)|>
  tbl_summary(by = study, 
              missing = "no", 
              label = list(status_death ~ "Status at end of the follow-up"), 
              digits = list(baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1, 
                            alcohol ~ 1))|>
  bold_labels()|>
  add_overall() 

## POPs ----
POPs_table_comp_cases <- bdd |>
  filter(als == 1) |>
  select(study, all_of(POPs_tot), -OCP_α_HCH) |>
  tbl_summary(by = "study") |>
  bold_labels() 

POPs_group_bis <- c("PCB_4", "PCB_DL", "PCB_NDL", "OCP_HCB", "ΣDDT", "OCP_β_HCH", "OCP_γ_HCH", "Σchlordane", "OCP_PeCB", "ΣPBDE")
POPs_group_labels_bis <- c(
  "Most prevalent PCBs" = "PCB_4",
  "Dioxin-like PCBs" = "PCB_DL",
  "Non dioxin-like PCBs" = "PCB_NDL",
  "HCB" = "OCP_HCB",
  "ΣDDT" = "ΣDDT",
  "β-HCH" = "OCP_β_HCH",
  "γ-HCH" = "OCP_γ_HCH",
  "Σchlordane" = "Σchlordane",
  "PeCB" = "OCP_PeCB",
  "ΣPBDE" = "ΣPBDE")

POPs_boxplot_comp_cases <- bdd |>
  filter(als == 1) |>
  select(study, all_of(POPs_group_bis)) |>
  pivot_longer(cols = -study, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_group_labels_bis), 
         POPs = fct_recode(POPs, !!!POPs_group_labels_bis), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = study) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease concentrations (pg/ml)", fill = "Study") +
  coord_flip() +
  theme_lucid()

rm(POPs_group_bis, POPs_group_labels_bis)

POPs_heatmap_cases <- bdd |>
  filter(als == 1) |>
  select(all_of(POPs_tot)) |>
  rename(!!!POPs_labels) 

POPs_heatmap_cases <- cor(POPs_heatmap_cases, 
                           use = "pairwise.complete.obs", 
                           method = "pearson")

# Assemblage ----
results_descriptive <- list(
  danish = list(
    covar_danish = covar_danish, 
    
    POPs_table_danish = POPs_table_danish,
    POPs_table_danish_by_als = POPs_table_danish_by_als,
    POPs_boxplot_danish = POPs_boxplot_danish,
    POPs_boxplot_danish_by_als = POPs_boxplot_danish_by_als,
    POPs_group_boxplot_danish_by_als = POPs_group_boxplot_danish_by_als,
    POPs_group_boxplot_danish_by_death = POPs_group_boxplot_danish_by_death,    # among cases
    POPs_heatmap_danish = POPs_heatmap_danish,
    
    fattyacids_table_danish = fattyacids_table_danish,
    fattyacids_table_danish_by_als = fattyacids_table_danish_by_als,
    fattyacids_boxplot_danish = fattyacids_boxplot_danish,
    fattyacids_boxplot_danish_by_als = fattyacids_boxplot_danish_by_als, 
    fattyacids_heatmap_danish = fattyacids_heatmap_danish, 
    POPs_fattyacids_heatmap_danish = POPs_fattyacids_heatmap_danish, 
    
    proteomic_table_danish = proteomic_table_danish, 
    proteomic_table_danish_by_als = proteomic_table_danish_by_als, 
    proteomic_boxplot_danish = proteomic_boxplot_danish, 
    proteomic_boxplot_danish_by_als = proteomic_boxplot_danish_by_als, 
    proteomic_boxplot_danish_by_death = proteomic_boxplot_danish_by_death, 
    proteomic_heatmap_danish = proteomic_heatmap_danish),
  
  finnish = list(
    covar_finnish = covar_finnish,
    POPs_table_finnish = POPs_table_finnish,
    POPs_table_finnish_by_als = POPs_table_finnish_by_als,
    POPs_boxplot_finnish = POPs_boxplot_finnish,
    POPs_boxplot_finnish_by_als = POPs_boxplot_finnish_by_als,
    POPs_group_boxplot_finnish_by_als = POPs_group_boxplot_finnish_by_als,
    POPs_heatmap_finnish = POPs_heatmap_finnish,
    fattyacids_table_finnish = fattyacids_table_finnish,
    fattyacids_table_finnish_by_als = fattyacids_table_finnish_by_als,
    fattyacids_boxplot_finnish = fattyacids_boxplot_finnish,
    fattyacids_boxplot_finnish_by_als = fattyacids_boxplot_finnish_by_als, 
    fattyacids_heatmap_finnish = fattyacids_heatmap_finnish, 
    POPs_fattyacids_heatmap_finnish = POPs_fattyacids_heatmap_finnish),
  comp = list(
    covar_comp = covar_comp,
    POPs_table_comp = POPs_table_comp,
    POPs_boxplot_comp = POPs_boxplot_comp,
    fattyacids_table_comp = fattyacids_table_comp,
    fattyacids_boxplot_comp = fattyacids_boxplot_comp, 
    covar_comp_cases = covar_comp_cases, 
    POPs_table_comp_cases = POPs_table_comp_cases, 
    POPs_boxplot_comp_cases = POPs_boxplot_comp_cases, 
    POPs_heatmap_cases = POPs_heatmap_cases))


rm(
  covar_danish,
  POPs_table_danish,
  POPs_table_danish_by_als,
  POPs_boxplot_danish,
  POPs_boxplot_danish_by_als,
  POPs_group_boxplot_danish_by_als, 
  POPs_group_boxplot_danish_by_death, 
  POPs_heatmap_danish, 
  fattyacids_table_danish,
  fattyacids_table_danish_by_als,
  fattyacids_boxplot_danish,
  fattyacids_boxplot_danish_by_als,
  fattyacids_heatmap_danish,
  POPs_fattyacids_heatmap_danish, 
  proteomic_table_danish, 
  proteomic_table_danish_by_als, 
  proteomic_boxplot_danish, 
  proteomic_boxplot_danish_by_als, 
  proteomic_boxplot_danish_by_death,
  proteomic_heatmap_danish, 
  
  covar_finnish,
  POPs_table_finnish,
  POPs_table_finnish_by_als,
  POPs_boxplot_finnish,
  POPs_boxplot_finnish_by_als,
  POPs_group_boxplot_finnish_by_als, 
  POPs_heatmap_finnish, 
  fattyacids_table_finnish,
  fattyacids_table_finnish_by_als,
  fattyacids_boxplot_finnish,
  fattyacids_boxplot_finnish_by_als,
  fattyacids_heatmap_finnish, 
  POPs_fattyacids_heatmap_finnish, 
  
  covar_comp,
  POPs_table_comp,
  POPs_boxplot_comp,
  fattyacids_table_comp,
  fattyacids_boxplot_comp, 
  
  covar_comp_cases, 
  POPs_table_comp_cases, 
  POPs_boxplot_comp_cases, 
  POPs_heatmap_cases)

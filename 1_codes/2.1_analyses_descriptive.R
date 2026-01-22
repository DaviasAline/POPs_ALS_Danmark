# Aline Davias 
# April 9, 2025

source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")

# danish data ----
## metadata ----
covar_danish <- bdd_danish |> 
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"))|>
  select(
    als, baseline_age, diagnosis_age, death_age, birth_year, 
    sex, marital_status, education, alcohol, smoking, bmi, cholesterol)|>
  tbl_summary(by = als)|>
  bold_labels()|>
  add_overall() 

covar_danish_cases <- 
  bdd |> 
  filter(study == "Danish") |>
  filter(als == 1) |>
  mutate(status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1"), 
         sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other")) |>
  select(
    baseline_age, diagnosis_age, death_age, birth_year, 
    follow_up_death, status_death,
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, cholesterol)|>
  tbl_summary(missing = "no", 
              label = list(status_death ~ "Status at end of the follow-up"), 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1, 
                            alcohol ~ 1))|>
  bold_labels() |>
  add_n(statistic = "{N_miss} ({p_miss}%)", col_label = "**N missing**")


covar_danish_by_death <- 
  bdd |> 
  filter(study == "Danish") |>
  filter(als == 1) |>
  mutate(status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1"), 
         sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other")) |>
  select(
    status_death, baseline_age, diagnosis_age, death_age, birth_year, 
    follow_up_death, status_death,
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, cholesterol)|>
  tbl_summary(by = status_death, 
              missing = "no", 
              digits = list(birth_year ~ 0, 
                            baseline_age ~ 0, 
                            diagnosis_age ~ 0, 
                            death_age ~ 0, 
                            bmi ~ 1, 
                            cholesterol ~ 1, 
                            alcohol ~ 1))|>
  bold_labels()


## POPs ----
POPs_table_danish <- descrip_num(data = bdd_danish, vars = POPs_tot)
POPs_table_danish <- left_join(POPs_table_danish, bdd_danish_POPs_loq, by = "variable") |>
  mutate(variable = factor(variable, levels = POPs_tot_labels), 
         variable = fct_recode(variable, !!!POPs_tot_labels)) |>
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

POPs_table_danish_by_death <- bdd_danish |>
  filter(als == 1) |>                         # among cases 
  select(status_death, all_of(POPs_tot)) |>
  mutate(
    status_death = as.character(status_death), 
    status_death = fct_recode(status_death, 
                              "Alive" = "0",
                              "Deceased" = "1")) |>
  tbl_summary(by = status_death, 
              digits = all_continuous() ~1) |>
  bold_labels() 

POPs_boxplot_danish <- bdd_danish |>
  select(all_of(POPs_tot)) |>
  pivot_longer(cols = everything(), names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_tot_labels), 
         POPs = fct_recode(POPs, !!!POPs_tot_labels), 
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
  mutate(POPs = factor(POPs, levels = POPs_tot_labels), 
         POPs = fct_recode(POPs, !!!POPs_tot_labels), 
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
  labs(x = "POPs", y = "Pre-disease plasma concentrations (log10 pg/ml)", fill = "") +
  coord_flip() +
  theme_lucid()

POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  POPs_group_bis,
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"))

POPs_heatmap_danish_group <- 
  bdd_danish |> 
  select(all_of(POPs_group_bis)) |>
  rename(!!!pollutant_labels_bis) 
POPs_heatmap_danish_group <- 
  cor(POPs_heatmap_danish_group, 
      use = "pairwise.complete.obs", 
      method = "pearson")

POPs_heatmap_danish <- 
  bdd_danish |> 
  select(all_of(POPs)) |>
  rename(!!!POPs_labels) 
POPs_heatmap_danish <- 
  cor(POPs_heatmap_danish, 
      use = "pairwise.complete.obs", 
      method = "pearson")


POPs_heatmap_danish_cases_group <- 
  bdd_danish |> 
  filter(als == 1) |>
  select(all_of(POPs_group_bis)) |>
  rename(!!!pollutant_labels_bis) 
POPs_heatmap_danish_cases_group <- 
  cor(POPs_heatmap_danish_cases_group, 
      use = "pairwise.complete.obs", 
      method = "pearson")

POPs_heatmap_danish_cases <- 
  bdd_danish |> 
  filter(als == 1) |>
  select(all_of(POPs)) |>
  rename(!!!POPs_labels) 
POPs_heatmap_danish_cases <- 
  cor(POPs_heatmap_danish_cases, 
      use = "pairwise.complete.obs", 
      method = "pearson")

## fatty acids ----
fattyacids_table_danish <- 
  descrip_num(data = bdd_danish, vars = c(fattyacids_tot)) |> 
  mutate(variable = factor(variable, levels = fattyacids_tot_labels), 
         variable = fct_recode(variable, !!!fattyacids_tot_labels)) |>
  arrange(variable) 

fattyacids_table_danish_by_als <- bdd_danish |>
  select(als, all_of(fattyacids_tot)) |>
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
  select(all_of(fattyacids_tot)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_tot_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_tot_labels), 
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
  select(als, all_of(fattyacids_tot)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_tot_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_tot_labels), 
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

fattyacids_heatmap_danish <- bdd_danish |> select(all_of(fattyacids_tot)) |>
  rename(!!!fattyacids_tot_labels) 
fattyacids_heatmap_danish <- cor(fattyacids_heatmap_danish, 
                                 use = "pairwise.complete.obs", 
                                 method = "pearson")

POPs_fattyacids_heatmap_danish <- 
  heatmap_cor_pairwise(fattyacids_tot, POPs, decimal = 1, data = bdd_danish)

## proteomic ----
### all ----
proteomic_table_danish <- descrip_num(data = bdd_danish, vars = proteomic) |>
  mutate(variable = gsub("^proteomic_(immun_res|neuro_explo|metabolism)_", "", variable))

# proteomic_table_danish_by_als <- bdd_danish |>
#   select(als, all_of(proteomic)) |>
#   mutate(
#     als = as.character(als), 
#     als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
#     als = fct_relevel(als, "Cases", "Controls")) |>
#   tbl_summary(by = als, 
#               digits = all_continuous() ~1) |>
#   bold_labels() |>
#   add_overall() |>
#   add_p()

proteomic_boxplot_danish <- bdd_danish |>
  select(all_of(proteomic)) |>
  remove_var_label() |>
  pivot_longer(cols = everything(), names_to = "proteomic", values_to = "values") |>
  mutate(proteomic = factor(proteomic, levels = proteomic_labels),
         proteomic = fct_recode(proteomic, !!!proteomic_labels),
         proteomic = fct_rev(proteomic)) |>
  arrange(proteomic) |>
  ggplot() +
  aes(x = proteomic, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "proteomic", y = "Pre-disease serum concentrations (unit?)") +
  coord_flip() +
  theme_lucid()

proteomic_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(proteomic)) |>
  remove_var_label() |>
  pivot_longer(cols = -als, names_to = "proteomic", values_to = "values") |>
  mutate(proteomic = factor(proteomic, levels = proteomic_labels), 
         proteomic = fct_recode(proteomic, !!!proteomic_labels), 
         proteomic = fct_rev(proteomic),
         als = as.character(als), 
         als = fct_recode(als, 
                          "Controls" = "0",
                          "Cases" = "1")) |>
  arrange(proteomic) |>
  ggplot() +
  aes(x = proteomic, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "proteomic", y = "Pre-disease serum concentrations (unit?)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

proteomic_heatmap_danish <- 
  bdd_danish |> 
  select(all_of(proteomic)) |>
  rename(!!!proteomic_labels) 
proteomic_heatmap_danish <- 
  cor(proteomic_heatmap_danish, 
      use = "pairwise.complete.obs", 
      method = "pearson")

proteomic_heatmap_danish_neuro_explo <- 
  bdd_danish |> 
  select(all_of(proteomic_neuro_explo)) |>
  rename(!!!proteomic_neuro_explo_labels) 
proteomic_heatmap_danish_neuro_explo <- 
  cor(proteomic_heatmap_danish_neuro_explo, 
      use = "pairwise.complete.obs", 
      method = "pearson")

proteomic_heatmap_danish_immun_res <- 
  bdd_danish |> 
  select(all_of(proteomic_immun_res)) |>
  rename(!!!proteomic_immun_res_labels) 
proteomic_heatmap_danish_immun_res <- 
  cor(proteomic_heatmap_danish_immun_res, 
      use = "pairwise.complete.obs", 
      method = "pearson")

proteomic_heatmap_danish_metabolism <- 
  bdd_danish |> 
  select(all_of(proteomic_metabolism)) |>
  rename(!!!proteomic_metabolism_labels) 
proteomic_heatmap_danish_metabolism <- 
  cor(proteomic_heatmap_danish_metabolism, 
      use = "pairwise.complete.obs", 
      method = "pearson")



proteomic_boxplot_danish_by_death <- bdd_danish |>
  filter(als == 1) |>
  select(status_death, all_of(proteomic)) |>
  remove_var_label() |>
  pivot_longer(cols = -status_death, names_to = "proteomic", values_to = "values") |>
  mutate(proteomic = factor(proteomic, levels = proteomic_labels), 
         proteomic = fct_recode(proteomic, !!!proteomic_labels), 
         proteomic = fct_rev(proteomic),
         status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1")) |>
  arrange(proteomic) |>
  ggplot() +
  aes(x = proteomic, y = values, fill = status_death) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "proteomic", y = "Pre-disease serum concentrations (unit?)", fill = "") +
  coord_flip() +
  theme_lucid()

### NEFL ----
bdd_danish_sensi_3 <- 
  bdd_danish |>                                                                 # densityplot all prot
  filter(match != 159) |>
  group_by(match) |>
  mutate(follow_up_ter = follow_up[als == 1]) |>
  ungroup() |>
  mutate(
    follow_up_ter = cut(
      follow_up_ter,
      breaks = quantile(follow_up_ter, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE, 
      labels = c("Filtered to follow-up < 11.2\nyears (tertile 1)", "Filtered to follow-up 11.2-16.2\nyears (tertile 2)", "Filtered to follow-up > 16.2\nyears (tertile 3)"))) |>
  select(proteomic_neuro_explo_NEFL, als, follow_up_ter) |>
  mutate(als = factor(as.character(als)), 
         als = fct_recode(als, "Controls" = "0", "Cases" = "1")) 


figure_NEFL <- 
  
  bdd_danish |>                                                                 # densityplot all prot
  filter(match != 159) |>
  select(proteomic_neuro_explo_NEFL, als) |>
  mutate(als = factor(as.character(als)), 
         als = fct_recode(als, "Controls" = "0", "Cases" = "1")) |>
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_density(alpha = 0.4) +
  scale_fill_hue(direction = -1) + 
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Density", title = "All", fill = "ALS") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        title = element_text(size = 14))  +


  bdd_danish |>                                                                 # boxplot all prot
  filter(match != 159) |>
  select(proteomic_neuro_explo_NEFL, als) |>
  mutate(als = factor(as.character(als)), 
         als = fct_recode(als, "Controls" = "0", "Cases" = "1")) |>
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_hue(direction = -1) +
  scale_x_continuous(limits = c(1, 5)) +
  labs(x = "Neurofilament light polypeptide (NPX)", fill = "ALS") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),  
        axis.text.y = element_blank())  +


bdd_danish_sensi_3 |>                                                           # densityplot par tertile
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_density(alpha = 0.4) +
  scale_fill_hue(direction = -1) +
  labs(
    x = "Neurofilament light polypeptide (NPX)",
    y = "Density",
    fill = "ALS") +
  theme_minimal() +
  facet_wrap(vars(follow_up_ter)) +
  xlim(1, 5) +
  ylim(0, 1)  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),  
        legend.position = "none", 
        strip.text = element_text(hjust = 0, size = 14), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 8)) +
  
  
  bdd_danish_sensi_3 |>
  ggplot() +
  aes(x = proteomic_neuro_explo_NEFL, fill = als) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_hue(direction = -1) +
  labs(
    x = "Neurofilament light polypeptide (NPX)",
    fill = "ALS") +
  theme_minimal() +
  facet_wrap(vars(follow_up_ter)) +
  xlim(1, 5)  + 
  theme(strip.text = element_blank(), 
        legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  
  plot_layout(heights = c(5, 1, 5, 1 ), ncol = 1)

rm(bdd_danish_sensi_3)




## EVs ----
EVs_table_danish <- descrip_num(data = bdd_danish, vars = EVs)

EVs_table_danish_by_als <- bdd_danish |>
  select(als, all_of(EVs)) |>
  mutate(
    als = as.character(als), 
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"), 
    als = fct_relevel(als, "Cases", "Controls")) |>
  tbl_summary(by = als, 
              digits = all_continuous() ~1) |>
  bold_labels() |>
  add_overall() |>
  add_p()

EVs_boxplot_danish <- bdd_danish |>
  select(all_of(EVs)) |>
  remove_var_label() |>
  pivot_longer(cols = everything(), names_to = "EVs", values_to = "values") |>
  mutate(EVs = factor(EVs, levels = EVs_labels),
         EVs = fct_recode(EVs, !!!EVs_labels),
         EVs = fct_rev(EVs)) |>
  arrange(EVs) |>
  ggplot() +
  aes(x = EVs, y = values) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "EVs", y = "Pre-disease serum concentrations (unit?)") +
  coord_flip() +
  theme_lucid()

EVs_boxplot_danish_by_als <- bdd_danish |>
  select(als, all_of(EVs)) |>
  remove_var_label() |>
  pivot_longer(cols = -als, names_to = "EVs", values_to = "values") |>
  mutate(EVs = factor(EVs, levels = EVs_labels), 
         EVs = fct_recode(EVs, !!!EVs_labels), 
         EVs = fct_rev(EVs),
         als = as.character(als), 
         als = fct_recode(als, 
                          "Controls" = "0",
                          "Cases" = "1")) |>
  arrange(EVs) |>
  ggplot() +
  aes(x = EVs, y = values, fill = als) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "EVs", y = "Pre-disease serum concentrations (unit?)", fill = "ALS") +
  coord_flip() +
  theme_lucid()

EVs_heatmap_danish <- 
  bdd_danish |> 
  select(all_of(EVs)) |>
  rename(!!!EVs_labels) 
EVs_heatmap_danish <- 
  cor(EVs_heatmap_danish, 
      use = "pairwise.complete.obs", 
      method = "pearson")

EVs_boxplot_danish_by_death <- bdd_danish |>
  filter(als == 1) |>
  select(status_death, all_of(EVs)) |>
  remove_var_label() |>
  pivot_longer(cols = -status_death, names_to = "EVs", values_to = "values") |>
  mutate(EVs = factor(EVs, levels = EVs_labels), 
         EVs = fct_recode(EVs, !!!EVs_labels), 
         EVs = fct_rev(EVs),
         status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1")) |>
  arrange(EVs) |>
  ggplot() +
  aes(x = EVs, y = values, fill = status_death) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "EVs", y = "Pre-disease serum concentrations (unit?)", fill = "") +
  coord_flip() +
  theme_lucid()


# finnish data ----
## metadata ----
covar_finnish <- bdd_finnish |> 
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls"))|>
  select(
    als, baseline_age, diagnosis_age, death_age, birth_year, 
    sex, marital_status, education, alcohol, smoking, bmi, cholesterol)|>
  tbl_summary(by = als)|>
  bold_labels()|>
  add_overall() 

## POPs ----
POPs_table_finnish <- bdd |> filter(!study %in% "Danish") 
POPs_table_finnish <- descrip_num(data = POPs_table_finnish, vars = POPs_tot) |> 
  mutate(variable = factor(variable, levels = POPs_tot_labels), 
         variable = fct_recode(variable, !!!POPs_tot_labels)) |>
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
  mutate(POPs = factor(POPs, levels = POPs_tot_labels), 
         POPs = fct_recode(POPs, !!!POPs_tot_labels), 
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
  mutate(POPs = factor(POPs, levels = POPs_tot_labels), 
         POPs = fct_recode(POPs, !!!POPs_tot_labels), 
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
  rename(!!!POPs_tot_labels) |>
  select(-"α-HCH")

POPs_heatmap_finnish <- cor(POPs_heatmap_finnish, 
                            use = "pairwise.complete.obs", 
                            method = "pearson")

## fatty acids ----
fattyacids_table_finnish <-
  descrip_num(data = bdd_finnish, vars = fattyacids_tot) |> 
  mutate(variable = fct_recode(variable, !!!fattyacids_tot_labels)) 

fattyacids_table_finnish_by_als <- 
  bdd_finnish |>
  select(als, all_of(fattyacids_tot)) |>
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
  select(all_of(fattyacids_tot)) |>
  pivot_longer(cols = everything(), names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_tot_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_tot_labels), 
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
  select(als, all_of(fattyacids_tot)) |>
  pivot_longer(cols = -als, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = factor(fattyacids, levels = fattyacids_tot_labels), 
         fattyacids = fct_recode(fattyacids, !!!fattyacids_tot_labels), 
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
  select(all_of(fattyacids_tot)) |>
  rename(!!!fattyacids_tot_labels) 
fattyacids_heatmap_finnish <- cor(fattyacids_heatmap_finnish, 
                                  use = "pairwise.complete.obs", 
                                  method = "pearson")

POPs_fattyacids_heatmap_finnish <- heatmap_cor_pairwise(fattyacids_tot, POPs_finnish, decimal = 1, data = bdd_finnish)

# comparison danish / finnish on the total population ----
## metadata ----
covar_comp <- bdd |>
  select("study", "sex", "marital_status", "smoking", "alcohol", "education", 
         "bmi", "cholesterol", "blod_sys", "blod_dias", 
         "baseline_age", "diagnosis_age", "birth_year") |>
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
  mutate(POPs = factor(POPs, levels = POPs_tot_labels), 
         POPs = fct_recode(POPs, !!!POPs_tot_labels), 
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
  select(study, all_of(fattyacids_tot)) |>
  tbl_summary(by = "study") |>
  bold_labels() |>
  add_n()

fattyacids_boxplot_comp <- bdd |>
  select(study, all_of(fattyacids_tot)) |>
  pivot_longer(cols = -study, names_to = "fattyacids", values_to = "values") |>
  mutate(fattyacids = fct_recode(fattyacids, !!!fattyacids_tot_labels)) |>
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
                                   "Deceased" = "1"), 
         sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other")) |>
  select(
    study, baseline_age, diagnosis_age, death_age, birth_year, 
    follow_up_death, status_death,
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, cholesterol)|>
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

covar_comp_cases_2cat <- bdd |> 
  filter(als == 1) |>
  mutate(status_death = as.character(status_death), 
         status_death = fct_recode(status_death, 
                                   "Alive" = "0",
                                   "Deceased" = "1"), 
         sex = fct_relevel(sex, "Male", "Female"), 
         smoking_2cat = fct_relevel(smoking_2cat, "Ever", "Never"), 
         marital_status_2cat = fct_relevel(marital_status_2cat, "Married/cohabit", "Other"), 
         study_2cat = fct_recode(
           study_2cat, 
           "Danish EPIC" = "Danish",
           "Finnish Health Surveys" = "Finnish")) |>
  select(
    study_2cat, baseline_age, diagnosis_age, death_age, birth_year, 
    follow_up, follow_up_death, status_death,
    sex, marital_status_2cat, education_merged, alcohol, smoking_2cat, bmi, cholesterol)|>
  tbl_summary(by = study_2cat, 
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
### Table 
POPs_table_comp_cases <- bdd |>
  filter(als == 1) |>
  select(study, all_of(POPs_tot), -OCP_α_HCH) |>
  tbl_summary(by = "study") |>
  bold_labels() 

POPs_table_comp_cases_2cat <- bdd |>
  filter(als == 1) |>
  select(study_2cat, all_of(POPs_tot), -OCP_α_HCH) |>
  mutate(
    study_2cat = fct_recode(
      study_2cat, 
      "Danish EPIC" = "Danish",
      "Finnish Health Surveys" = "Finnish")) |>
  tbl_summary(by = "study_2cat") |>
  bold_labels() 

### Boxplots
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

POPs_boxplot_comp_cases_2cat <- bdd |>
  filter(als == 1) |>
  select(study_2cat, all_of(POPs_group_bis)) |>
  mutate(
    study_2cat = fct_recode(
      study_2cat, 
      "Danish EPIC" = "Danish",
      "Finnish Health Surveys" = "Finnish")) |>
  pivot_longer(cols = -study_2cat, names_to = "POPs", values_to = "values") |>
  mutate(POPs = factor(POPs, levels = POPs_group_labels_bis), 
         POPs = fct_recode(POPs, !!!POPs_group_labels_bis), 
         POPs = fct_rev(POPs)) |>
  arrange(POPs) |>
  ggplot() +
  aes(x = POPs, y = values, fill = study_2cat) +
  geom_boxplot() +
  scale_fill_hue(direction = 1, 
                 guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Pre-disease concentrations (pg/ml)", fill = "Study") +
  coord_flip() +
  theme_lucid()

rm(POPs_group_bis, POPs_group_labels_bis)

### Heatmaps
POPs_group_bis <- c(setdiff(POPs_group, "PCB_4"), "OCP_PeCB")
pollutant_labels_bis <- set_names(
  POPs_group_bis,
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE", "Pentachlorobenzene"))

POPs_heatmap_cases_group <- bdd |>
  filter(als == 1) |>
  select(all_of(POPs_group_bis)) |>
  rename(!!!pollutant_labels_bis) 

POPs_heatmap_cases_group <- cor(POPs_heatmap_cases_group, 
                          use = "pairwise.complete.obs", 
                          method = "pearson")

POPs_heatmap_cases <- bdd |>
  filter(als == 1) |>
  select(all_of(POPs)) |>
  rename(!!!POPs_labels) |>
  select(-"α-HCH")

POPs_heatmap_cases <- cor(POPs_heatmap_cases, 
                           use = "pairwise.complete.obs", 
                           method = "pearson")

# Assemblage ----
results_descriptive <- list(
  danish = list(
    covar_danish = covar_danish, 
    covar_danish_cases = covar_danish_cases, 
    covar_danish_by_death = covar_danish_by_death,                              # among cases
    
    POPs_table_danish = POPs_table_danish,
    POPs_table_danish_by_als = POPs_table_danish_by_als,
    POPs_table_danish_by_death = POPs_table_danish_by_death,                    # among cases 
    POPs_boxplot_danish = POPs_boxplot_danish,
    POPs_boxplot_danish_by_als = POPs_boxplot_danish_by_als,
    POPs_group_boxplot_danish_by_als = POPs_group_boxplot_danish_by_als,
    POPs_group_boxplot_danish_by_death = POPs_group_boxplot_danish_by_death,    # among cases
    POPs_heatmap_danish = POPs_heatmap_danish,
    POPs_heatmap_danish_group = POPs_heatmap_danish_group,
    POPs_heatmap_danish_cases_group = POPs_heatmap_danish_cases_group, 
    POPs_heatmap_danish_cases = POPs_heatmap_danish_cases, 
    
    fattyacids_table_danish = fattyacids_table_danish,
    fattyacids_table_danish_by_als = fattyacids_table_danish_by_als,
    fattyacids_boxplot_danish = fattyacids_boxplot_danish,
    fattyacids_boxplot_danish_by_als = fattyacids_boxplot_danish_by_als, 
    fattyacids_heatmap_danish = fattyacids_heatmap_danish, 
    POPs_fattyacids_heatmap_danish = POPs_fattyacids_heatmap_danish, 
    
    proteomic_table_danish = proteomic_table_danish, 
    # proteomic_table_danish_by_als = proteomic_table_danish_by_als, 
    proteomic_boxplot_danish = proteomic_boxplot_danish, 
    proteomic_boxplot_danish_by_als = proteomic_boxplot_danish_by_als, 
    proteomic_heatmap_danish = proteomic_heatmap_danish, 
    proteomic_heatmap_danish_neuro_explo = proteomic_heatmap_danish_neuro_explo, 
    proteomic_heatmap_danish_immun_res = proteomic_heatmap_danish_immun_res, 
    proteomic_heatmap_danish_metabolism = proteomic_heatmap_danish_metabolism, 
    proteomic_boxplot_danish_by_death = proteomic_boxplot_danish_by_death, 
    figure_NEFL = figure_NEFL, 
    
    EVs_table_danish = EVs_table_danish, 
    EVs_table_danish_by_als = EVs_table_danish_by_als, 
    EVs_boxplot_danish = EVs_boxplot_danish, 
    EVs_boxplot_danish_by_als = EVs_boxplot_danish_by_als, 
    EVs_heatmap_danish = EVs_heatmap_danish, 
    EVs_boxplot_danish_by_death = EVs_boxplot_danish_by_death),
  
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
    covar_comp_cases_2cat = covar_comp_cases_2cat, 
    POPs_table_comp_cases = POPs_table_comp_cases, 
    POPs_table_comp_cases_2cat = POPs_table_comp_cases_2cat, 
    POPs_boxplot_comp_cases = POPs_boxplot_comp_cases, 
    POPs_boxplot_comp_cases_2cat = POPs_boxplot_comp_cases_2cat, 
    POPs_heatmap_cases = POPs_heatmap_cases, 
    POPs_heatmap_cases_group = POPs_heatmap_cases_group))


rm(
  bdd_danish_POPs_loq, 
  covar_danish,
  covar_danish_cases, 
  covar_danish_by_death, 
  POPs_table_danish,
  POPs_table_danish_by_als,
  POPs_table_danish_by_death, 
  POPs_boxplot_danish,
  POPs_boxplot_danish_by_als,
  POPs_group_boxplot_danish_by_als, 
  POPs_group_boxplot_danish_by_death, 
  POPs_heatmap_danish, 
  POPs_heatmap_danish_group, 
  POPs_heatmap_danish_cases, 
  POPs_heatmap_danish_cases_group, 
  fattyacids_table_danish,
  fattyacids_table_danish_by_als,
  fattyacids_boxplot_danish,
  fattyacids_boxplot_danish_by_als,
  fattyacids_heatmap_danish,
  POPs_fattyacids_heatmap_danish, 
  proteomic_table_danish, 
  # proteomic_table_danish_by_als, 
  proteomic_boxplot_danish, 
  proteomic_boxplot_danish_by_als, 
  proteomic_heatmap_danish, 
  proteomic_heatmap_danish_neuro_explo, 
  proteomic_heatmap_danish_immun_res, 
  proteomic_heatmap_danish_metabolism,
  proteomic_boxplot_danish_by_death, 
  figure_NEFL, 
  EVs_table_danish, 
  EVs_table_danish_by_als, 
  EVs_boxplot_danish, 
  EVs_boxplot_danish_by_als, 
  EVs_heatmap_danish, 
  EVs_boxplot_danish_by_death, 
  
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
  covar_comp_cases_2cat, 
  POPs_table_comp_cases, 
  POPs_table_comp_cases_2cat, 
  POPs_boxplot_comp_cases, 
  POPs_boxplot_comp_cases_2cat, 
  POPs_heatmap_cases, 
  POPs_heatmap_cases_group)

# Aline Davias
# 03/02/2025



# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/0_functions.R")
source("~/Documents/POP_ALS_2025_02_03/1_codes/1_data_loading.R")
library("see")
library("scales")
library('plotly')
library('gridExtra')  
library('ggeffects')
library('mgcv')
library('purrr')
library('tidyr')
library('ggstance')
library('gratia')
library('glue')

# descriptif ----
## covariates table ----
descrip_covar <- bdd_danish %>% 
  mutate(
    als = as.character(als),
    als = fct_recode(als, "Controls" = "0", "Cases" = "1"),
    als = fct_relevel(als, "Cases", "Controls")) %>%
  select(
    als, baseline_age, diagnosis_age, death_age, 
    sex, marital_status, education, alcohol, smoking, bmi, cholesterol) %>%
  tbl_summary(by = als) %>%
  bold_labels() %>%
  add_overall() 

## exposures ----
bdd_danish_long <- bdd_danish %>% 
  select(sample, all_of(POPs)) %>%
  pivot_longer(cols = -sample, names_to = "POPs", values_to = "Values") %>%
  filter(!POPs %in% c("PeCB", "α_HCH")) %>%
  mutate(POPs =  gsub("_", "-", POPs), 
         POPs = fct_recode(POPs, 
             "Trans-nonachlor" = "Transnonachlor",
             "p,p'-DDE" = "pp-DDE",
             "p,p'-DDT" = "pp-DDT", 
             "PBDE-153" = "BDE-153",
             "PBDE-47" = "BDE-47",
             "PBDE-99" = "BDE-99",), 
         POPs_group =  fct_recode(POPs, 
             "ΣPBDE" = "PBDE-153",
             "ΣPBDE" = "PBDE-47",
             "ΣPBDE" = "PBDE-99",
             "Σchlordane" = "Oxychlordane",
             "PCB_NDL" = "PCB-101",
             "PCB_DL" = "PCB-118",
             "PCB_NDL" = "PCB-138",
             "PCB_NDL" = "PCB-153",
             "PCB_DL" = "PCB-156",
             "PCB_NDL" = "PCB-170",
             "PCB_NDL" = "PCB-180",
             "PCB_NDL" = "PCB-183",
             "PCB_NDL" = "PCB-187",
             "PCB_NDL" = "PCB-28",
             "PCB_NDL" = "PCB-52",
             "PCB_NDL" = "PCB-74",
             "PCB_NDL" = "PCB-99",
             "ΣDDT" = "p,p'-DDE",
             "ΣDDT" = "p,p'-DDT",
             "Σchlordane" = "Trans-nonachlor",
             "ΣHCH" = "β-HCH",
             "ΣHCH" = "γ-HCH"
           ), 
         POPs_group_2 =  fct_recode(POPs, 
                                    "PCBs" = "PCB-101",
                                    "PCBs" = "PCB-118",
                                    "PCBs" = "PCB-138",
                                    "PCBs" = "PCB-153",
                                    "PCBs" = "PCB-156",
                                    "PCBs" = "PCB-170",
                                    "PCBs" = "PCB-180",
                                    "PCBs" = "PCB-183",
                                    "PCBs" = "PCB-187",
                                    "PCBs" = "PCB-28",
                                    "PCBs" = "PCB-52",
                                    "PCBs" = "PCB-74",
                                    "PCBs" = "PCB-99",
                                    
                                    "PBDEs" = "PBDE-153",
                                    "PBDEs" = "PBDE-47",
                                    "PBDEs" = "PBDE-99",
                                    
                                    "OCPs" = "HCB",
                                    #"OCPs" = "α-HCH", 
                                    "OCPs" = "β-HCH",
                                    "OCPs" = "γ-HCH",
                                    "OCPs" = "p,p'-DDE",
                                    "OCPs" = "p,p'-DDT",
                                    "OCPs" = "Oxychlordane",
                                    "OCPs" = "Trans-nonachlor" 
                                    #"OCPs" = "PeCB"
                                    ), 
         POPs_group_2 = fct_relevel(POPs_group_2, 
                                    "PCBs", "PBDEs", "OCPs"),
         POPs = fct_relevel(POPs,      
                            "PCB-28", 
                            "PCB-52", 
                            "PCB-74", 
                            "PCB-99",
                            "PCB-101", 
                            "PCB-118", 
                            "PCB-156",  
                            "PCB-138", 
                            "PCB-153", 
                            "PCB-170", 
                            "PCB-180", 
                            "PCB-183",
                            "PCB-187", 
                            
                            "PBDE-47", 
                            "PBDE-99", 
                            "PBDE-153", 
                            
                            "p,p'-DDT", 
                            "p,p'-DDE",
                            "β-HCH", 
                            "γ-HCH", 
                            "Oxychlordane", 
                            "Trans-nonachlor", 
                            "HCB"))
descrip_expo <- bdd_danish_long %>%
  mutate(
    POPs = fct_relevel(POPs, 
                       "PCB-118", "PCB-156", "PCB-28", "PCB-52", "PCB-74", "PCB-99",
                       "PCB-101", "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183",
                       "PCB-187", "PBDE-47", "PBDE-99", "PBDE-153", "p,p'-DDT", "p,p'-DDE",
                       "β-HCH", "γ-HCH", "Oxychlordane", "Trans-nonachlor", "HCB"), 
    POPs_group = fct_relevel(POPs_group, 
                             "PCB_DL", "PCB_NDL", "ΣPBDE", "ΣDDT", "ΣHCH", "Σchlordane", "HCB"), 
    POPs = factor(POPs, levels = rev(levels(POPs)))) %>%
  arrange(desc(POPs_group)) %>%
  ggplot() +
  aes(x = POPs, y = Values, fill = POPs_group) +
  geom_boxplot() +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)", fill = "POP groups") +
  coord_flip() +
  theme_lucid() +
  scale_fill_brewer(palette = "Set2", direction = 1) 


descrip_expo_group <- bdd_danish_long %>%
  mutate(
    POPs = fct_relevel(POPs, 
                       "PCB-118", "PCB-156", "PCB-28", "PCB-52", "PCB-74", "PCB-99",
                       "PCB-101", "PCB-138", "PCB-153", "PCB-170", "PCB-180", "PCB-183",
                       "PCB-187", "PBDE-47", "PBDE-99", "PBDE-153", "p,p'-DDT", "p,p'-DDE",
                       "β-HCH", "γ-HCH", "Oxychlordane", "Trans-nonachlor", "HCB"), 
    POPs_group = fct_relevel(POPs_group, 
                             "PCB_DL", "PCB_NDL", "ΣPBDE", "ΣDDT", "ΣHCH", "Σchlordane", "HCB"), 
    POPs = factor(POPs, levels = rev(levels(POPs)))) %>%
  arrange(desc(POPs_group)) %>%
  ggplot() +
  aes(x = POPs_group, y = Values, fill = POPs_group) +
  geom_boxplot() +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  labs(x = "POPs", y = "Values (pg/ml, log transformed)", fill = "POP groups") +
  coord_flip() +
  theme_lucid() +
  scale_fill_brewer(palette = "Set2", direction = 1) 

descrip_expo_group <- ggplot(bdd_danish_long) +
  aes(x = POPs_group, y = Values, fill = POPs_group) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log", 
                     labels = number_format(accuracy = 1)) +
  coord_flip() +
  theme_lucid() +
  labs(x = "", y = "POP group concentrations (ng/ml)", fill = "POP groups")

# ggsave("~/Documents/POP_ALS_2025_02_03/2_output/descrip_expo (Figure 1).tiff", 
#        descrip_expo)

# statistics ----
## effects of the covariates on ALS ----
covar <- 
  clogit(als ~ strata(match) + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tbl_regression(exponentiate = TRUE, 
                 estimate_fun = label_number(accuracy = .01, decimal.mark = "."),
                 pvalue_fun = custom_pvalue_fun) |>
  bold_labels()

## main analysis ----
### model 1 ----
#### spline transformation ----
model1_spline <- data.frame(variable = character(),
                            df = integer(),
                            OR = numeric(),
                            lower_CI = numeric(),
                            upper_CI = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_spline <- rbind(model1_spline, data.frame(variable = var,
                                                   df = df_value, 
                                                   OR = OR,
                                                   lower_CI = lower_CI,
                                                   upper_CI = upper_CI,
                                                   "p-value" = p_value))
}

options(scipen = 999)
model1_spline[] <- lapply(model1_spline, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_spline <- model1_spline %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quartile transformation ----
model1_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  model_summary <- tidy(model)
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quart <- rbind(model1_quart, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model1_quart[] <- lapply(model1_quart, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_quart <- model1_quart %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quart") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model1_quadratic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quadratic <- rbind(model1_quadratic, data.frame(variable = var,
                                                         df = df_value, 
                                                         OR = OR,
                                                         lower_CI = lower_CI,
                                                         upper_CI = upper_CI,
                                                         "p-value" = p_value))
}

options(scipen = 999)
model1_quadratic[] <- lapply(model1_quadratic, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_quadratic <- model1_quadratic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)


#### cubic transformation ----
model1_cubic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_cubic <- rbind(model1_cubic, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

options(scipen = 999)
model1_cubic[] <- lapply(model1_cubic, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_cubic <- model1_cubic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### gamm ----
model1_gamm <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model1_gamm[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### model 2 ----
# matched on sex and age, adjusted on for smoking_2cat_i, BMI, serum total cholesterol_i, marital status and education

#### spline transformation ----
model2_spline <- data.frame(variable = character(),
                            df = integer(),
                            OR = numeric(),
                            lower_CI = numeric(),
                            upper_CI = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_spline <- rbind(model2_spline, data.frame(variable = var,
                                                   df = df_value, 
                                                   OR = OR,
                                                   lower_CI = lower_CI,
                                                   upper_CI = upper_CI,
                                                   "p-value" = p_value))
}

options(scipen = 999)
model2_spline[] <- lapply(model2_spline, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_spline <- model2_spline %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quartile transformation ----
model2_quart <- data.frame(variable = character(),
                           df = integer(),
                           OR = numeric(),
                           lower_CI = numeric(),
                           upper_CI = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- clogit(formula, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl(paste0("^", var), term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quart <- rbind(model2_quart, data.frame(variable = var,
                                                 df = df_value, 
                                                 OR = OR,
                                                 lower_CI = lower_CI,
                                                 upper_CI = upper_CI,
                                                 "p-value" = p_value))
}

model2_quart[] <- lapply(model2_quart, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_quart <- model2_quart %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quart") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model2_quadratic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quadratic <- rbind(model2_quadratic, 
                            data.frame(variable = var,
                                       df = df_value, 
                                       OR = OR,
                                       lower_CI = lower_CI,
                                       upper_CI = upper_CI,
                                       "p-value" = p_value))
}

options(scipen = 999)
model2_quadratic[] <- lapply(model2_quadratic, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_quadratic <- model2_quadratic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### cubic transformation ----
model2_cubic <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_cubic <- rbind(model2_cubic, 
                        data.frame(variable = var,
                                   df = df_value, 
                                   OR = OR,
                                   lower_CI = lower_CI,
                                   upper_CI = upper_CI,
                                   "p-value" = p_value))
}

options(scipen = 999)
model2_cubic[] <- lapply(model2_cubic, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_cubic <- model2_cubic %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### gamm ----
model2_gamm <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### models 3 ----
#### spline transformation ----
model3_spline_HCB <- 
  clogit(als ~ 
           ns(HCB, df = 4) + 
           strata(match) + 
           PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("HCB", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_PCB_DL <- 
  clogit(als ~ 
           ns(PCB_DL, df = 4) + 
           strata(match) + 
           HCB + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_PCB_NDL <- 
  clogit(als ~ 
           ns(PCB_NDL, df = 4) + 
           strata(match) + 
           HCB + PCB_DL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_ΣPBDE <- 
  clogit(als ~ 
           ns(ΣPBDE, df = 4) + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_ΣDDT <- 
  clogit(als ~ 
           ns(ΣDDT, df = 4) + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + β_HCH + Σchlordane +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_β_HCH <- 
  clogit(als ~ 
           ns(β_HCH, df = 4) + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("β_HCH", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline_Σchlordane <- 
  clogit(als ~ 
           ns(Σchlordane, df = 4) + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH +
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term)) |>
  mutate(df = c(1, 2, 3, 4))

model3_spline <- bind_rows(
  model3_spline_HCB, model3_spline_PCB_DL, model3_spline_PCB_NDL, model3_spline_ΣPBDE, model3_spline_ΣDDT, model3_spline_β_HCH, model3_spline_Σchlordane) |>
  mutate(
    variable = gsub('ns\\(', '', term), 
    variable = gsub(', df = 4)1', '', variable),
    variable = gsub(', df = 4)2', '', variable),
    variable = gsub(', df = 4)3', '', variable),
    variable = gsub(', df = 4)4', '', variable),
    model = "copollutant_spline", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)
model3_spline[] <- lapply(model3_spline, function(x) if(is.numeric(x)) round(x, 2) else x)

rm(model3_spline_HCB, model3_spline_PCB_DL, model3_spline_PCB_NDL, model3_spline_ΣPBDE, model3_spline_ΣDDT, model3_spline_β_HCH, model3_spline_Σchlordane)

#### quartile tranformation ----
model3_quart_PCB_DL <- 
  clogit(als ~ 
           PCB_DL_quart + 
           strata(match) + 
           HCB + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_PCB_NDL <- 
  clogit(als ~ 
           PCB_NDL_quart + 
           strata(match) + 
           HCB + PCB_DL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_HCB <- 
  clogit(als ~ 
           HCB_quart + 
           strata(match) + 
           PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("HCB_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_ΣDDT <- 
  clogit(als ~ 
           ΣDDT_quart + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_β_HCH <- 
  clogit(als ~ 
           β_HCH_quart + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("β_HCH_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_Σchlordane <- 
  clogit(als ~ 
           Σchlordane_quart + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart_ΣPBDE <- 
  clogit(als ~ 
           ΣPBDE_quart + 
           strata(match) + 
           HCB + PCB_DL + PCB_NDL + ΣDDT + β_HCH + Σchlordane + 
           smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
         data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE_quart", term)) |>
  mutate(df = c(2, 3, 4))

model3_quart <- bind_rows(
  model3_quart_HCB, model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_ΣPBDE, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane) |>
  mutate(
    variable = gsub("_quartQ2", "", term), 
    variable = gsub("_quartQ3", "", variable), 
    variable = gsub("_quartQ4", "", variable), 
    model = "copollutant_quart", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)
model3_quart[] <- lapply(model3_quart, function(x) if(is.numeric(x)) round(x, 2) else x)

rm(model3_quart_HCB, model3_quart_PCB_DL, model3_quart_PCB_NDL, model3_quart_ΣPBDE, model3_quart_ΣDDT, model3_quart_β_HCH, model3_quart_Σchlordane)

#### quadratic tranformation ----
model3_quadra_HCB <- 
  glm(als ~ 
        poly(HCB, 2) + 
        baseline_age + sex + 
        PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("HCB", term)) |>
  mutate(df = c(1, 2))

model3_quadra_PCB_DL <- 
  glm(als ~ 
        poly(PCB_DL, 2) + 
        baseline_age + sex + 
        HCB + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term)) |>
  mutate(df = c(1, 2))

model3_quadra_PCB_NDL <- 
  glm(als ~ 
        poly(PCB_NDL, 2) + 
        baseline_age + sex + 
        HCB + PCB_DL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term)) |>
  mutate(df = c(1, 2))

model3_quadra_ΣPBDE <- 
  glm(als ~ 
        poly(ΣPBDE, 2) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2))

model3_quadra_ΣDDT <- 
  glm(als ~ 
        poly(ΣDDT, 2) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term)) |>
  mutate(df = c(1, 2))

model3_quadra_β_HCH <- 
  glm(als ~ 
        poly(β_HCH, 2) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("β_HCH", term)) |>
  mutate(df = c(1, 2))

model3_quadra_Σchlordane <- 
  glm(als ~ 
        poly(Σchlordane, 2) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term)) |>
  mutate(df = c(1, 2))

model3_quadratic <- bind_rows(
  model3_quadra_HCB, model3_quadra_PCB_DL, model3_quadra_PCB_NDL, model3_quadra_ΣPBDE, model3_quadra_ΣDDT, model3_quadra_β_HCH, model3_quadra_Σchlordane) |>
  mutate(
    variable = gsub("poly\\(", "", term),
    variable = gsub(", 2)1", "", variable),
    variable = gsub(", 2)2", "", variable),
    model = "copollutant_quadra", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)
model3_quadratic[] <- lapply(model3_quadratic, function(x) if(is.numeric(x)) round(x, 2) else x)

rm(model3_quadra_HCB, model3_quadra_PCB_DL, model3_quadra_PCB_NDL, model3_quadra_ΣPBDE, model3_quadra_ΣDDT, model3_quadra_β_HCH, model3_quadra_Σchlordane)

#### cubic tranformation ----
model3_cubic_HCB <- 
  glm(als ~ 
        poly(HCB, 3) + 
        baseline_age + sex + 
        PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("HCB", term)) |>
  mutate(df = c(1, 2, 3))

model3_cubic_PCB_DL <- 
  glm(als ~ 
        poly(PCB_DL, 3) + 
        baseline_age + sex + 
        HCB + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_DL", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_PCB_NDL <- 
  glm(als ~ 
        poly(PCB_NDL, 3) + 
        baseline_age + sex + 
        HCB + PCB_DL + ΣPBDE + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("PCB_NDL", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_ΣPBDE <- 
  glm(als ~ 
        poly(ΣPBDE, 3) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣDDT + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣPBDE", term)) |>
  mutate(df = c(1, 2, 3))

model3_cubic_ΣDDT <- 
  glm(als ~ 
        poly(ΣDDT, 3) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + β_HCH + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("ΣDDT", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_β_HCH <- 
  glm(als ~ 
        poly(β_HCH, 3) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + Σchlordane + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("β_HCH", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic_Σchlordane <- 
  glm(als ~ 
        poly(Σchlordane, 3) + 
        baseline_age + sex + 
        HCB + PCB_DL + PCB_NDL + ΣPBDE + ΣDDT + β_HCH + 
        smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
      family = binomial, data = bdd_danish) |>
  tidy()  |>
  filter(grepl("Σchlordane", term))  |>
  mutate(df = c(1, 2, 3))

model3_cubic <- bind_rows(
  model3_cubic_HCB, model3_cubic_PCB_DL, model3_cubic_PCB_NDL, model3_cubic_ΣPBDE, model3_cubic_ΣDDT, model3_cubic_β_HCH, model3_cubic_Σchlordane) |>
  mutate(
    variable = gsub("poly\\(", "", term),
    variable = gsub(", 3)1", "", variable),
    variable = gsub(", 3)2", "", variable),
    variable = gsub(", 3)3", "", variable),
    model = "copollutant_cubic", 
    df = as.character(df),
    OR = exp(estimate), 
    lower_CI = exp(estimate - 1.96 * std.error), 
    upper_CI = exp(estimate + 1.96 * std.error)) |> 
  select(variable, model, df, OR, lower_CI, upper_CI, p.value)
model3_cubic[] <- lapply(model3_cubic, function(x) if(is.numeric(x)) round(x, 2) else x)

rm(model3_cubic_HCB, model3_cubic_PCB_DL, model3_cubic_PCB_NDL, model3_cubic_ΣPBDE, model3_cubic_ΣDDT, model3_cubic_β_HCH, model3_cubic_Σchlordane)

### heterogeneity tests ----
#### model 1 spline ----
heterogeneity_base_spline <- data.frame(variable = character(),
                                        model = factor(), 
                                        p.value_heterogeneity = numeric(), 
                                        stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_spline <- rbind(heterogeneity_base_spline, 
                                     data.frame(variable = var,
                                                model = "base_spline",
                                                p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

#### model 1 quartile ----
heterogeneity_base_quart <- data.frame(variable = character(),
                                       model = factor(),
                                       p.value_heterogeneity = numeric(), 
                                       stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_base_quart <- rbind(heterogeneity_base_quart, 
                                    data.frame(variable = var,
                                               model = "base_quart",
                                               p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

#### model 2 spline ----
heterogeneity_adjusted_spline <- data.frame(variable = character(),
                                            model = factor(), 
                                            p.value_heterogeneity = numeric(), 
                                            stringsAsFactors = FALSE)

for (var in POPs_group) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_spline <- rbind(heterogeneity_adjusted_spline, 
                                         data.frame(variable = var,
                                                    model = "adjusted_spline",
                                                    p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

#### model 2 quartile ----
heterogeneity_adjusted_quart <- data.frame(variable = character(),
                                           model = factor(), 
                                           p.value_heterogeneity = numeric(), 
                                           stringsAsFactors = FALSE)

for (var in POPs_group_quart) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_heterogeneity <- anova$`Pr(>|Chi|)`[2]
  
  heterogeneity_adjusted_quart <- rbind(heterogeneity_adjusted_quart, 
                                        data.frame(variable = var,
                                                   model = "adjusted_quart",
                                                   p.value_heterogeneity = p.value_heterogeneity))
}
rm(var, test_1, test_2, formula, anova, p.value_heterogeneity)

heterogeneity_tests <- 
  bind_rows(heterogeneity_base_spline, heterogeneity_base_quart, heterogeneity_adjusted_spline, heterogeneity_adjusted_quart) %>%
  mutate(variable = gsub("_quart", "", variable))

### trend tests ----
#### model 1 quartile ----
trend_base <- data.frame(variable = character(),
                         model = factor(), 
                         p.value_trend = numeric(), 
                         stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  test_1 <- clogit(als ~ strata(match), data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match)"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_trend <- anova$`Pr(>|Chi|)`[2]
  
  trend_base <- rbind(trend_base, 
                      data.frame(variable = var,
                                 model = "base_quart",
                                 p.value_trend = p.value_trend))
}
rm(var, test_1, test_2, formula, anova, p.value_trend)

#### model 2 quartile ----
trend_adjusted <- data.frame(variable = character(),
                             model = factor(), 
                             p.value_trend = numeric(), 
                             stringsAsFactors = FALSE)

for (var in POPs_group_quart_med) {
  
  test_1 <- clogit(als ~ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, data = bdd_danish)
  
  formula <- as.formula(paste("als ~", var, "+ strata(match) + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  test_2 <- clogit(formula, data = bdd_danish)
  
  anova <- anova(test_1, test_2, test = "LR")
  p.value_trend <- anova$`Pr(>|Chi|)`[2]
  
  trend_adjusted <- rbind(trend_adjusted, 
                          data.frame(variable = var,
                                     model = "adjusted_quart",
                                     p.value_trend = p.value_trend))
}
rm(var, test_1, test_2, formula, anova, p.value_trend)

trend_tests <- 
  bind_rows(trend_base, trend_adjusted) %>%
  mutate(variable = gsub("_quart_med", "", variable))

### merging the main results ----
main_results <- bind_rows(model1_spline, 
                          model1_quart, 
                          model1_quadratic,
                          model1_cubic, 
                          model2_spline, 
                          model2_quart,
                          model2_quadratic,
                          model2_cubic,
                          model3_spline, 
                          model3_quart, 
                          model3_quadratic, 
                          model3_cubic) %>% 
  mutate(variable = gsub("_quart", "", variable)) %>%
  arrange(variable) %>%
  mutate("95%CI" = paste(lower_CI, ", ", upper_CI)) %>%
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value"), 
         lower_CI, upper_CI) 

main_results <- left_join(main_results, heterogeneity_tests, by = c("variable", "model"))
main_results <- left_join(main_results, trend_tests, by = c("variable", "model"))
main_results[] <- lapply(main_results, function(x) if(is.numeric(x)) round(x, 2) else x)

results_spline <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI) |>
  filter(model %in% c('base_spline', 'adjusted_spline', 'copollutant_spline')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_spline"), 
         contains("adjusted_spline"), 
         contains("copollutant_spline")) 
colnames(results_spline) <- gsub('_spline', '', colnames(results_spline))

results_quart <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI) |>
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         quartiles = df,
         contains("base_quart"), 
         contains("adjusted_quart"), 
         contains("copollutant_quart")) 
colnames(results_quart) <- gsub('_quart', '', colnames(results_quart))

results_quadratic <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI) |>
  filter(model %in% c('base_quadra', 'adjusted_quadra', 'copollutant_quadra')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         degree = df,
         contains("base_quadra"), 
         contains("adjusted_quadra"), 
         contains("copollutant_quadra")) 
colnames(results_quadratic) <- gsub('_quadra', '', colnames(results_quadratic))

results_cubic <- 
  main_results |>
  select(-p.value_heterogeneity, -p.value_trend, -lower_CI, -upper_CI) |>
  filter(model %in% c('base_cubic', 'adjusted_cubic', 'copollutant_cubic')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         degree = df,
         contains("base_cubic"), 
         contains("adjusted_cubic"), 
         contains("copollutant_cubic")) 
colnames(results_cubic) <- gsub('_cubic', '', colnames(results_cubic))


rm(model1_quart, model1_spline, model1_quadratic, model1_cubic,
   model2_quart, model2_spline, model2_quadratic, model2_cubic,
   model3_quart, model3_spline, model3_quadratic, model3_cubic,
   heterogeneity_base_quart, heterogeneity_base_spline,
   heterogeneity_adjusted_spline, heterogeneity_adjusted_quart, 
   trend_base, trend_adjusted, 
   heterogeneity_tests, trend_tests)


## sensitivity analyses 95% ----
### model 1 ----
# adjusted for sex and age

#### spline transformation ----
model1_spline_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age"))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_spline_95 <- rbind(model1_spline_95, data.frame(variable = var,
                                                         df = df_value, 
                                                         OR = OR,
                                                         lower_CI = lower_CI,
                                                         upper_CI = upper_CI,
                                                         "p-value" = p_value))
}

options(scipen = 999)
model1_spline_95[] <- lapply(model1_spline_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_spline_95 <- model1_spline_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model1_quadratic_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age"))
  
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quadratic_95 <- rbind(model1_quadratic_95, data.frame(variable = var,
                                                               df = df_value, 
                                                               OR = OR,
                                                               lower_CI = lower_CI,
                                                               upper_CI = upper_CI,
                                                               "p-value" = p_value))
}

options(scipen = 999)
model1_quadratic_95[] <- lapply(model1_quadratic_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_quadratic_95 <- model1_quadratic_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model1_cubic_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age"))
  
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_cubic_95 <- rbind(model1_cubic_95, data.frame(variable = var,
                                                       df = df_value, 
                                                       OR = OR,
                                                       lower_CI = lower_CI,
                                                       upper_CI = upper_CI,
                                                       "p-value" = p_value))
}

options(scipen = 999)
model1_cubic_95[] <- lapply(model1_cubic_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_cubic_95 <- model1_cubic_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)


### model 2 ----
# adjusted for sex, age, smoking_2cat_i, BMI, serum total cholesterol_i, marital status and education

#### spline transformation ----
model2_spline_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_spline_95 <- rbind(model2_spline_95, data.frame(variable = var,
                                                         df = df_value, 
                                                         OR = OR,
                                                         lower_CI = lower_CI,
                                                         upper_CI = upper_CI,
                                                         "p-value" = p_value))
}

options(scipen = 999)
model2_spline_95[] <- lapply(model2_spline_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_spline_95 <- model2_spline_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model2_quadratic_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quadratic_95 <- rbind(model2_quadratic_95, 
                                data.frame(variable = var,
                                           df = df_value, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           "p-value" = p_value))
}

options(scipen = 999)
model2_quadratic_95[] <- lapply(model2_quadratic_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_quadratic_95 <- model2_quadratic_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model2_cubic_95 <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_95) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_cubic_95 <- rbind(model2_cubic_95, 
                           data.frame(variable = var,
                                      df = df_value, 
                                      OR = OR,
                                      lower_CI = lower_CI,
                                      upper_CI = upper_CI,
                                      "p-value" = p_value))
}

options(scipen = 999)
model2_cubic_95[] <- lapply(model2_cubic_95, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_cubic_95 <- model2_cubic_95 %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

## sensitivity analyses without outliers ----
### model 1  ----
# adjusted for sex and age

#### spline transformation ----
model1_spline_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_spline_outlier <- rbind(model1_spline_outlier, data.frame(variable = var,
                                                                   df = df_value, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   "p-value" = p_value))
}

options(scipen = 999)
model1_spline_outlier[] <- lapply(model1_spline_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_spline_outlier <- model1_spline_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### quadratic transformation ----
model1_quadratic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_quadratic_outlier <- rbind(model1_quadratic_outlier, data.frame(variable = var,
                                                                         df = df_value, 
                                                                         OR = OR,
                                                                         lower_CI = lower_CI,
                                                                         upper_CI = upper_CI,
                                                                         "p-value" = p_value))
}

options(scipen = 999)
model1_quadratic_outlier[] <- lapply(model1_quadratic_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_quadratic_outlier <- model1_quadratic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model1_cubic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  
  model1_cubic_outlier <- rbind(model1_cubic_outlier, data.frame(variable = var,
                                                                 df = df_value, 
                                                                 OR = OR,
                                                                 lower_CI = lower_CI,
                                                                 upper_CI = upper_CI,
                                                                 "p-value" = p_value))
}

options(scipen = 999)
model1_cubic_outlier[] <- lapply(model1_cubic_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model1_cubic_outlier <- model1_cubic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "base_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)


### model 2 ----
# adjusted for sex, age, smoking_2cat_i, BMI, serum total cholesterol_i, marital status and education

#### spline transformation ----
model2_spline_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ ns(", var, ", df=4) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- glm(formula, family = binomial, data = bdd_danish)
  model_summary <- tidy(model) %>% filter(grepl("ns\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_spline_outlier <- rbind(model2_spline_outlier, data.frame(variable = var,
                                                                   df = df_value, 
                                                                   OR = OR,
                                                                   lower_CI = lower_CI,
                                                                   upper_CI = upper_CI,
                                                                   "p-value" = p_value))
}

options(scipen = 999)
model2_spline_outlier[] <- lapply(model2_spline_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_spline_outlier <- model2_spline_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_spline") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var)

#### quadratic transformation ----
model2_quadratic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=2) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_quadratic_outlier <- rbind(model2_quadratic_outlier, 
                                    data.frame(variable = var,
                                               df = df_value, 
                                               OR = OR,
                                               lower_CI = lower_CI,
                                               upper_CI = upper_CI,
                                               "p-value" = p_value))
}

options(scipen = 999)
model2_quadratic_outlier[] <- lapply(model2_quadratic_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_quadratic_outlier <- model2_quadratic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_quadra") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### cubic transformation ----
model2_cubic_outlier <- data.frame(
  variable = character(),
  df = integer(),
  OR = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ poly(", var, ", degree=3) + sex + baseline_age + smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  bdd_danish_red <- bdd_danish %>% filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  model_summary <- tidy(model) %>% filter(grepl("poly\\(", term))
  OR <- exp(model_summary$estimate)
  lower_CI <- exp(model_summary$estimate - 1.96 * model_summary$std.error)
  upper_CI <- exp(model_summary$estimate + 1.96 * model_summary$std.error)
  p_value <- model_summary$p.value
  df_value <- model_summary$term
  model2_cubic_outlier <- rbind(model2_cubic_outlier, 
                                data.frame(variable = var,
                                           df = df_value, 
                                           OR = OR,
                                           lower_CI = lower_CI,
                                           upper_CI = upper_CI,
                                           "p-value" = p_value))
}

options(scipen = 999)
model2_cubic_outlier[] <- lapply(model2_cubic_outlier, function(x) if(is.numeric(x)) round(x, 2) else x)
model2_cubic_outlier <- model2_cubic_outlier %>% 
  mutate(
    df = str_sub(df, start = -1), 
    model = "adjusted_cubic") %>%
  select(variable, model, everything())
rm(model, lower_CI, upper_CI, df_value, formula, p_value, OR, model_summary, var, bdd_danish_red)

#### gamm ----
model2_gamm_outliers <- list()

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste("als ~ s(", var, ") + sex + baseline_age + 
                              smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i"))
  model <- gam(formula, family = binomial, method = 'REML', data = bdd_danish)
  model_summary <- summary(model)
  model2_gamm_outliers[[var]] <- model_summary
}

rm(var, formula, model, model_summary)

### merging the sensitivity results ----
sensitivity_results <- bind_rows(model1_spline_outlier,
                                  model1_quadratic_outlier,
                                 model1_cubic_outlier, 
                                 model2_spline_outlier,
                                model2_quadratic_outlier,
                                 model2_cubic_outlier) %>% 
  arrange(variable) %>%
  mutate(
    "95%CI" = paste(lower_CI, ", ", upper_CI)) %>%
  select(-starts_with("lower"), -starts_with("upper")) %>%
  select(variable, 
         model,
         df,
         starts_with("OR"), 
         starts_with("95%"), 
         starts_with("p.value")) 

results_spline_outliers <- 
  sensitivity_results |>
  filter(model %in% c('base_spline', 'adjusted_spline', 'copollutant_spline')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_spline"), 
         contains("adjusted_spline"), 
         contains("copollutant_spline")) 
colnames(results_spline_outliers) <- gsub('_spline', '', colnames(results_spline_outliers))

results_quadratic_outliers <- 
  sensitivity_results |>
  filter(model %in% c('base_quadra', 'adjusted_quadra', 'copollutant_quadra')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_quadra"), 
         contains("adjusted_quadra"), 
         contains("copollutant_quadra")) 
colnames(results_quadratic_outliers) <- gsub('_quadra', '', colnames(results_quadratic_outliers))

results_cubic_outliers <- 
  sensitivity_results |>
  filter(model %in% c('base_cubic', 'adjusted_cubic', 'copollutant_cubic')) |>
  pivot_wider(
    names_from = model,  
    values_from = c(OR, `95%CI`, p.value)) |>
  select(variable, 
         df,
         contains("base_cubic"), 
         contains("adjusted_cubic"), 
         contains("copollutant_cubic")) 
colnames(results_cubic_outliers) <- gsub('_cubic', '', colnames(results_cubic_outliers))

rm(model1_spline_95, model1_spline_outlier,
   model1_quadratic_95, model1_quadratic_outlier,
   model1_cubic_95, model1_cubic_outlier, 
   model2_spline_95, model2_spline_outlier,
   model2_quadratic_95, model2_quadratic_outlier,
   model2_cubic_95, model2_cubic_outlier)

# figures ----
## quartiles ----
plot_quart <- main_results %>% 
  filter(model %in% c('base_quart', 'adjusted_quart', 'copollutant_quart')) %>%
  mutate(df = fct_recode(df, "Quartile 2" = "2", "Quartile 3" = "3", "Quartile 4" = "4"), 
         p.value_shape = ifelse(p.value<0.05, "p-value<0.05", "p-value≥0.05"), 
         model = fct_recode(model, 
                            "Adjusted model" = "adjusted_quart",
                            "Base model" = "base_quart",
                            "Copollutant model" = "copollutant_quart"),
         model = fct_relevel(model, 'Base model', 'Adjusted model', 'Copollutant model'), 
         df = fct_relevel(df, "Quartile 4", "Quartile 3", "Quartile 2" ), 
         variable = fct_recode(variable, 
                               "Most prevalent\nPCBs" = "PCB_4",
                               "Dioxin-like\nPCBs" = "PCB_DL",
                               "Non-dioxin-like\nPCBs" = "PCB_NDL",
                               "β-HCH" = "β_HCH")) %>%
  ggplot(aes(x = df, y = OR, ymin = lower_CI, ymax = upper_CI, color = p.value_shape)) +
  geom_pointrange(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_grid(rows = dplyr::vars(variable), cols = dplyr::vars(model), switch = "y") +  
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black")) +
  labs(x = "POPs", y = "Odds Ratio (OR)", color = "p-value") +
  theme_lucid() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = "bottom", 
        strip.text.y = element_text(hjust = 0.5)) +
  coord_flip()

## model 1 ----
### spline ----
plot_base_spline <- list()

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + strata(match)"))
  model <- clogit(formula, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  pred <- predict(model, newdata = new_data, type = "lp", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) + 
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = "spline") +
    theme_minimal()
  
  plot_base_spline[[var]] <- plot
}
rm(var, formula, model, new_data, pred, plot)

### spline 95 ----
plot_base_spline_95 <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_base_spline_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### spline outliers ----
plot_base_spline_outlier <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_base_spline_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### quadratic ----
plot_base_quadratic <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic') +
    theme_minimal()
  
  plot_base_quadratic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)


### quadratic 95 ----
plot_base_quadratic_95 <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_base_quadratic_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### quadratic outliers ----
plot_base_quadratic_outlier <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_base_quadratic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic ----
plot_base_cubic <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic') +
    theme_minimal()
  
  plot_base_cubic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### cubic 95 ----
plot_base_cubic_95 <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_base_cubic_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic outliers ----
plot_base_cubic_outlier <- list()

covariates <- c("sex", "baseline_age")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_base_cubic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)



### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

covariates <- c("sex", "baseline_age")

plot_base_gamm <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", round(p_value, 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels, covariates)

### gamm outliers ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group_outlier)

covariates <- c("sex", "baseline_age")

plot_base_gamm_outliers <- map(POPs_group_outlier, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", round(p_value, 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels, covariates)



## model 2 ----
### spline ----

plot_adjusted_spline <- list()

covariates <- c("smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + "), 
                               " + strata(match)"))
  
  model <- clogit(formula, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "lp", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = "spline") +
    theme_minimal()
  
  plot_adjusted_spline[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)


### spline 95 ----
plot_adjusted_spline_95 <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_adjusted_spline_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### spline outliers ----
plot_adjusted_spline_outlier <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ ns(", var, ", df = 4) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'spline without outliers') +
    theme_minimal()
  
  plot_adjusted_spline_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)


### quadratic  ----
plot_adjusted_quadratic <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic') +
    theme_minimal()
  
  plot_adjusted_quadratic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### quadratic 95 ----
plot_adjusted_quadratic_95 <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_adjusted_quadratic_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### quadratic outliers ----
plot_adjusted_quadratic_outlier <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 2) + ", 
                               paste(covariates, collapse = " + ")))
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  # Plot the estimated probability
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +  # Add CI ribbon
    labs(x = var, y = "Predicted probability of ALS", title = 'linear quadratic without outliers') +
    theme_minimal()
  
  plot_adjusted_quadratic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic ----
plot_adjusted_cubic <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  model <- glm(formula, family = binomial, data = bdd_danish)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic') +
    theme_minimal()
  
  plot_adjusted_cubic[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov)

### cubic 95 ----
plot_adjusted_cubic_95 <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_95) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_adjusted_cubic_95[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)

### cubic outliers -----
plot_adjusted_cubic_outlier <- list()

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", "cholesterol_i", "marital_status_2cat_i", "education_i")

for (var in POPs_group_outlier) {
  
  formula <- as.formula(paste0("als ~ poly(", var, ", degree = 3) + ", 
                               paste(covariates, collapse = " + ")))
  
  bdd_danish_red <- bdd_danish |> filter(!is.na(.data[[var]]))
  model <- glm(formula, family = binomial, data = bdd_danish_red)
  
  new_data <- data.frame(seq(min(bdd_danish[[var]], na.rm = TRUE), 
                             max(bdd_danish[[var]], na.rm = TRUE), 
                             length.out = 498))
  colnames(new_data) <- var
  
  new_data$match <- bdd_danish$match[1]
  
  for (cov in covariates) {
    if (is.numeric(bdd_danish[[cov]])) {
      new_data[[cov]] <- mean(bdd_danish[[cov]], na.rm = TRUE)  
    } else {
      new_data[[cov]] <- as.factor(names(which.max(table(bdd_danish[[cov]])))) 
    }
  }
  
  pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$upper <- pred$fit + 1.96 * pred$se.fit
  new_data$lower <- pred$fit - 1.96 * pred$se.fit
  
  new_data$prob <- exp(pred$fit) / (1 + exp(pred$fit))
  new_data$prob_lower <- exp(new_data$lower) / (1 + exp(new_data$lower))
  new_data$prob_upper <- exp(new_data$upper) / (1 + exp(new_data$upper))
  
  plot <- ggplot(new_data, aes_string(x = var, y = "prob")) +
    geom_line(color = "blue", size = 1) +  
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) + 
    labs(x = var, y = "Predicted probability of ALS", title = 'cubic without outliers') +
    theme_minimal()
  
  plot_adjusted_cubic_outlier[[var]] <- plot
}

rm(var, formula, model, new_data, pred, plot, cov, bdd_danish_red)



### gamm ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group)

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", 
                "cholesterol_i", "marital_status_2cat_i", "education_i")

plot_adjusted_gamm <- map(POPs_group, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", round(p_value, 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels, covariates)

### gamm outliers ----
pollutant_labels <- set_names(
  c("Dioxin-like PCBs","Non-dioxin-like PCBs", "Most prevalent PCBs","HCB","ΣDDT","β-HCH","Σchlordane","ΣPBDE"), 
  POPs_group_outlier)

covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", 
                "cholesterol_i", "marital_status_2cat_i", "education_i")

plot_adjusted_gamm_outliers <- map(POPs_group_outlier, function(var) {
  
  formula <- as.formula(glue::glue("als ~ s({var}) + {paste(covariates, collapse = ' + ')}"))
  
  model <- gam(formula, family = binomial, method = "REML", data = bdd_danish)
  
  bdd_pred <- bdd_danish %>%                                                    # création bdd avec expo + covariables ramenées à leur moyenne
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))), 
                  .names = "adj_{.col}")) %>%
    select(all_of(var), starts_with("adj_")) %>%
    rename_with(~ gsub("adj_", "", .x)) 
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[1, "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "< 0.01", round(p_value, 2))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels[var] 
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect') + 
    theme_minimal()
  p
}) %>% 
  set_names(POPs_group)
rm(pollutant_labels, covariates)

## model 3 ----
### gamm ----
POPs_group_bis <- setdiff(POPs_group, "PCB_4")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_bis)
covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", 
                "cholesterol_i", "marital_status_2cat_i", "education_i")

model <- gam(als ~ s(PCB_DL) + s(PCB_NDL) + s(HCB) + s(ΣDDT) + 
               s(β_HCH) + s(Σchlordane) + s(ΣPBDE) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = "REML", 
             data = bdd_danish)


plot_copollutant_gamm <- map(POPs_group_bis, function(var) {
  
  bdd_pred <- bdd_danish %>%
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))))) %>%
    mutate(across(setdiff(POPs_group_bis, var), 
                  ~ mean(., na.rm = TRUE))) %>%                                 # Fixe les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates), all_of(setdiff(POPs_group_bis, var)))  
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "p-value < 0.01", sprintf("p-value = %.2f", p_value))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  
  # edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  # p_value <- model_summary$s.table[1, "p-value"]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
}) %>% set_names(POPs_group_bis)
rm(POPs_group_bis, pollutant_labels_bis)


### gamm outliers ----
POPs_group_outlier_bis <- setdiff(POPs_group_outlier, "PCB_4_outlier")
pollutant_labels_bis <- set_names(
  c("Dioxin-like PCBs", "Non-dioxin-like PCBs", 
    "HCB", "ΣDDT", "β-HCH", "Σchlordane", "ΣPBDE"), 
  POPs_group_outlier_bis)
covariates <- c("sex", "baseline_age", "smoking_2cat_i", "bmi", 
                "cholesterol_i", "marital_status_2cat_i", "education_i")

model <- gam(als ~ s(PCB_DL_outlier) + s(PCB_NDL_outlier) + s(HCB_outlier) + s(ΣDDT_outlier) + 
               s(β_HCH_outlier) + s(Σchlordane_outlier) + s(ΣPBDE_outlier) + 
               sex + baseline_age + 
               smoking_2cat_i + bmi + cholesterol_i + marital_status_2cat_i + education_i, 
             family = binomial, 
             method = "REML", 
             data = bdd_danish)


plot_copollutant_gamm_outlier <- map(POPs_group_outlier_bis, function(var) {
  
  bdd_pred <- bdd_danish %>%
    mutate(across(all_of(covariates), 
                  ~ if (is.numeric(.)) mean(., na.rm = TRUE) else names(which.max(table(.))))) %>%
    mutate(across(setdiff(POPs_group_outlier_bis, var), 
                  ~ mean(., na.rm = TRUE))) %>%                                 # Fixe les autres POPs à leur moyenne
    select(all_of(var), all_of(covariates), all_of(setdiff(POPs_group_outlier_bis, var)))  
  
  pred <- predict(model, newdata = bdd_pred, type = "link", se.fit = TRUE)
  
  bdd_pred <- bdd_pred %>%
    mutate(prob = plogis(pred$fit),                                             # plogit does exp(pred$fit) / (1 + exp(pred$fit))
           prob_lower = plogis(pred$fit - 1.96 * pred$se.fit),
           prob_upper = plogis(pred$fit + 1.96 * pred$se.fit))
  
  model_summary <- summary(model)
  edf <- format(model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "edf"], nsmall = 1, digits = 1) 
  p_value <- model_summary$s.table[rownames(model_summary$s.table) == paste0("s(", var, ")"), "p-value"]
  p_value_text <- ifelse(p_value < 0.01, "p-value < 0.01", sprintf("p-value = %.2f", p_value))
  x_max <- max(bdd_danish[[var]], na.rm = TRUE)
  x_label <- pollutant_labels_bis[var]
  
  
  # edf <- format(model_summary$s.table[1, "edf"], nsmall = 1, digits = 1) 
  # p_value <- model_summary$s.table[1, "p-value"]
  
  p1 <- ggplot(bdd_pred, aes(x = .data[[var]], y = prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "blue", alpha = 0.2) +
    labs(x = var, y = "Predicted probability of ALS") +
    annotate("text", x = x_max, y = Inf, label = paste("EDF: ", edf, "\np-value: ", p_value_text, sep = ""),
             hjust = 1, vjust = 1.2, size = 4, color = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(color = 'white'),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) 
  
  p2 <- ggplot(bdd_pred) +
    aes(x = "", y = .data[[var]]) +
    geom_boxplot(fill = "blue") +
    coord_flip() +
    ylab(x_label) + 
    xlab("") + 
    theme_minimal()
  
  p <- p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(10, 1),
                           guides = 'collect')
  p
}) %>% set_names(POPs_group_outlier_bis)
rm(POPs_group_outlier_bis, pollutant_labels_bis)



# Aline Davias
# 18/03/2025

# data loading - package loading ----
source("~/Documents/POP_ALS_2025_02_03/1_codes/2.3_analyses_POPs_ALS_survival.R")


# Table 1 - description of the subjects ----
# Description of the subject characteristics of the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263)
table_1 <- results_descriptive$comp$covar_comp_cases |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# Figure 1 - boxplot expo ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
figure_1 <- results_descriptive$comp$POPs_boxplot_comp_cases


# Table S1 - exposure distribution ----
# Distribution of pre-disease POP concentrations in ALS cases from the Danish EPIC, the FMC, the FMCF and the MFH Finnish cohorts (total sample size=263).
table_S1 <- results_descriptive$comp$POPs_table_comp_cases |>
  as_flex_table() |>
  flextable::font(fontname = "Calibri", part = "all") |> 
  fontsize(size = 10, part = "all") |>
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  merge_at(i = 1, j = 1, part = "header") |>  
  merge_at(i = 1, j = 2, part = "header")  

# footnotes:
## pg/ml.
## POPs were summed as follows: Dioxin-like PCBs corresponds to PCBs 118 and 156; non-dioxin-like PCBs corresponds to PCBs 28, 52, 74, 99, 101, 138, 153, 170, 180, 183, 187; most prevalent PCBs corresponds to PCBs 118, 138, 153, 180; ΣPBDE corresponds to PBDEs 47, 99, 153; ΣDDT corresponds to p,p’-DDT and p,p’-DDE and finally Σchlordane corresponds to trans-nonanchlor and oxychlordane.
## α-HCH had more than 95% of the observations below the limit of quantification (LOQ) and were therefore excluded from our analyses. 
## HCB, γ-HCH, β-HCH and PeCB were not included in any group and were therefore studied alone.'



# Figure S1 - heatmap of correlation between POP exposures ---- 
# Pearson correlations between pre-disease POP plasma concentrations in the Danish Diet, Cancer and Health study cohort (sample size: 498). 
plot.new()
tiff(filename = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_S1.tiff", units = "mm", width = 300, height = 270, res = 300)
corrplot(results_descriptive$comp$POPs_heatmap_cases, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         number.digits = 1,
         tl.cex = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()


# Export ----
table_1 <- read_docx() |> body_add_flextable(table_1) 
print(table_1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_1.docx")

table_S1 <- read_docx() |> body_add_flextable(table_S1)
print(table_S1, target = "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/table_S1.docx")

ggsave(
  "~/Documents/POP_ALS_2025_02_03/2_output/Article_POPs_ALS_survival/figure_1.tiff",
  figure_1,
  height = 4,
  width = 8,
  units = "in")




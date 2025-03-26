FattyAcid <- 
  read.csv2("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_FattyAcid.csv") |>
  rename(sample = PSEUDO_ID) |>
  mutate(sample = as.character(sample)) |>
  filter(!is.na(sample))


Lipid <- read.csv("/Volumes/shared/EOME/Weisskopf/POPs-ALS/Data/Lipids/Finnish Lipids/THLBB2020_24_Lipid.csv", sep=";") |>
  rename(sample = PSEUDO_ID) |>
  mutate(sample = as.character(sample)) |>
  filter(!is.na(sample)) |>
  pivot_wider(names_from = "ANALYSEX", 
              values_from = "MEANVAL")


test <- full_join(bdd_finnish, FattyAcid, by = 'sample')
test <- full_join(test, Lipid, by = c('sample', 'Barcode'))
test <- test |>
  filter(!sample == "9047039803")

tbl_merge(tbls = list(bdd_finnish |> select(sample) |> tbl_summary(), 
                      FattyAcid |> select(sample) |> tbl_summary()), 
          tab_spanner = c('bdd_finish', 'bdd_lipids'))

Lipid |>
  select(PSEUDO_ID) |>
  mutate(PSEUDO_ID = as.character(PSEUDO_ID)) |>
  tbl_summary()

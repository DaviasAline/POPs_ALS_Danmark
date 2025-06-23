
str(bdd_danish_proteomic)

bdd_danish_proteomic |>
  select(Index) |>
  unique() |>
  View()

bdd_danish_proteomic |>
  select(MissingFreq, Panel, Panel_Version, QC_Warning, MaxLOD, PlateLOD, NPX, Normalization) |>
  tbl_summary()

bdd_danish_proteomic |>
  filter(QC_Warning == "Warning") |>
  View()



descrip_num(data = bdd_danish, vars = proteomic)

bdd_danish_proteomic <- bdd_danish_proteomic |>
  mutate(NPX_MaxLOD = ifelse(NPX > MaxLOD, "> Max LOD", "< Max LOD"), 
         NPX_PlateLOD = ifelse(NPX > PlateLOD, "> Plate LOD", "< Plate LOD"))

bdd_danish_proteomic |> select(NPX_MaxLOD, NPX_PlateLOD) |> tbl_summary()

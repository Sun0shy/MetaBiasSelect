library(MetaBiasSelect)

res <- run_metabias_pipeline(
  object = T_cells_merged,
  cluster_col = "T_NK_celltype_refined",
  cancer_col = "orig.ident",
  sample_col = "sample_id",
  cancers = c("CRC", "PDAC", "STAD")
)

print(res$plot)

save_metabias_outputs(
  result = res,
  outdir = "MetaBiasSelect_output",
  prefix = "T_NK"
)

############################################
### plot top 5 robust immune TCs

setwd("Users/arkajyotibhattacharya/Google Drive/My Drive/Projects_27092022/Spatial transcriptomics/Results/Seurat_objects_single_cell_reference_dataset/")   ### only in RStudio

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

library(Seurat)
library(data.table)
OUTPUT_DIR = "Users/arkajyotibhattacharya/Google Drive/My Drive/Projects_27092022/Spatial transcriptomics/Results/Seurat_objects_single_cell_reference_dataset/"

############################################
### read in 


corrected_mix_mat = data.frame(fread("/Users/arkajyotibhattacharya/Google Drive/My Drive/Projects_27092022/Spatial transcriptomics/Results/TLS_TC_projection_on_single_cells/mixing_matrix_with_single_cell_annotation.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
visium = readRDS('/Users/arkajyotibhattacharya/Google Drive/My Drive/Projects_27092022/Spatial transcriptomics/Data/EC_annotated_reference_atlas.rds') 
visium@assays$RNA@counts = t(as.matrix(corrected_mix_mat))

saveRDS(visium, file = paste0("Single_cell_reference_data_TLS_TCs.RDS"))



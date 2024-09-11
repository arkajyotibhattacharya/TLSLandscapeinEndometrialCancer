
############################################
### plot top 5 robust immune TCs

setwd("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Seurat_objects_with_TCs/")   ### only in RStudio

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)
OUTPUT_DIR = "/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Seurat_objects_with_TCs/"

############################################
### read in 


corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_PAM_patient_using_TLS_TCs/Project_data_on_independent_components_{8638d230-f85f-4752-8f0d-6d89e8d8fdf6}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
visium = readRDS('/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/PAM_patient_Tcells_Arko.rds') 
visium@assays$RNA@counts = t(as.matrix(corrected_mix_mat))

saveRDS(visium, file = paste0("PAM_patient_Tcells__TLS_TCs.RDS"))



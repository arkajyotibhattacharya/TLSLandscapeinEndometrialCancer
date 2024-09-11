
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


corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_tdln_tonsil_using_tls_tc/TDLN_EC4/Project_data_on_independent_components_{c81e4842-0e71-4dec-800c-82f67b9d03b6}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
visium = readRDS('/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/TDLN_tonsil.rds') 


############################################
### plot functions

create_seurat_object_with_TLS_TCs = function(image_id, corrected_mix_mat = corrected_mix_mat, visium = visium) {
  mapping_file <- visium[[image_id]]$spatial_data@images[[image_id]]@coordinates
  
  if(image_id =="OC1") {
    mapping_file = subset(mapping_file, col<185)  ### outside image region
  }
  
  common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
  mapping_file = mapping_file[common_samples,]
  corrected_mix_mat_current = corrected_mix_mat[common_samples,]
  
  corrected_mix_mat_transposed = t(corrected_mix_mat_current)
  new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial" )
  new.seurat.object@images$image = new(
    Class = 'VisiumV1'
    ,assay = "spatial"
    ,key = "image_"
    ,coordinates = mapping_file
    ,image = visium[[image_id]]$spatial_data@images[[image_id]]@image
    ,scale.factors = visium[[image_id]]$spatial_data@images[[image_id]]@scale.factors
    ,spot.radius = visium[[image_id]]$spatial_data@images[[image_id]]@spot.radius
  )
  
  saveRDS(new.seurat.object, file = paste0(image_id, "_TLS_TCs.RDS"))
}


create_seurat_object_with_TLS_TCs(image_id = "TDLN_EC4", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_tdln_tonsil_using_tls_tc/TDLN_EC6/Project_data_on_independent_components_{b37ce3b2-613c-4350-b05c-84ff2d42ad5d}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "TDLN_EC6", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_tdln_tonsil_using_tls_tc/TDLN_EC7/Project_data_on_independent_components_{2b00046d-e56a-43bb-8158-b6d16a1d5772}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "TDLN_EC7", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_tdln_tonsil_using_tls_tc/TDLN_EC8/Project_data_on_independent_components_{d8481bd3-11be-4d38-b2af-c62f12f108ac}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "TDLN_EC8", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Results/Projection_tdln_tonsil_using_tls_tc/TONSIL/Project_data_on_independent_components_{15aa4615-c5fb-4550-9685-0dab703e1fad}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus_independent_component_", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "TONSIL", corrected_mix_mat = corrected_mix_mat, visium = visium)


TONSIL_TLS_TCs = readRDS(file = "TONSIL_TLS_TCs.RDS")
SpatialFeaturePlot(TONSIL_TLS_TCs
                   , features = "TC4390"
                   , stroke = NA
                   , slot = "counts"
                   , image.alpha = 0.2
                   # , pt.size.factor = 1.8
                   # , pt.size.factor = 1.4  ### for OC1 zoomed
                   # , pt.size.factor = pt_size
                   , alpha = 1
                   ,interactive = FALSE)&NoAxes()

TDLN_EC4_TLS_TCs = readRDS(file = "TDLN_EC4_TLS_TCs.RDS")
SpatialFeaturePlot(TDLN_EC4_TLS_TCs
                   , features = "TC1622"
                   , stroke = NA
                   , slot = "counts"
                   , image.alpha = 0.2
                   # , pt.size.factor = 1.8
                   # , pt.size.factor = 1.4  ### for OC1 zoomed
                   # , pt.size.factor = pt_size
                   , alpha = 1
                   ,interactive = FALSE)&NoAxes()





############################################
### plot top 5 robust immune TCs

setwd("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/Seurat_objects_with_TCs/")   ### only in RStudio

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)
OUTPUT_DIR = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/Seurat_objects_with_TCs/"

############################################
### read in 


corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC1/Project_data_on_independent_components_{a46a91f0-06bb-408c-9b46-ae7e1d986d66}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
visium = readRDS('/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/ALL_VISIUM_FILES_MERGED.rds') 


# corrected_mix_mat_v1 = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_ovarian_TC/OC2/Project_data_on_independent_components_{5863f5bf-8a01-4070-9625-4a9f56296b3c}/mixing_matrix.tsv"), row.names = 1)
# colnames(corrected_mix_mat_v1)[1] = "V1"
# 
# corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, corrected_mix_mat_v1))
# 
# corrected_mix_mat_v1 = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_ovarian_TC/OC3/Project_data_on_independent_components_{8daedad5-3184-41f5-ac02-1735f2ee4e05}/mixing_matrix.tsv"), row.names = 1)
# colnames(corrected_mix_mat_v1)[1] = "V1"
# 
# corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, corrected_mix_mat_v1))
# 
# visium = readRDS('/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/ALL_VISIUM_FILES_MERGED.rds') 
# visium = readRDS('/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/NEW_EC_samples_ALL_VISIUM_FILES_MERGED.rds')

TLS_TCs = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/TLS_multiple_testing_correction/Multiple_testing_results.txt"))
tcs_to_plot = paste0("TC", which(TLS_TCs$MVP_reject_H0.nPerm..10000.FDR..0.05.conf..0.8=="TRUE"))
                    

############################################
### plot functions

create_seurat_object_with_TLS_TCs = function(image_id, corrected_mix_mat = corrected_mix_mat, visium = visium) {
  mapping_file <- visium@images[[image_id]]@coordinates
  
  if(image_id =="OC1") {
    mapping_file = subset(mapping_file, col<185)  ### outside image region
  }
  
  common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
  mapping_file = mapping_file[common_samples,]
  corrected_mix_mat_current = corrected_mix_mat[common_samples,]
  
  corrected_mix_mat_transposed = t(corrected_mix_mat_current)
  corrected_mix_mat_transposed  = corrected_mix_mat_transposed[tcs_to_plot,]
  new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial" )
  new.seurat.object@images$image = new(
    Class = 'VisiumV1'
    ,assay = "spatial"
    ,key = "image_"
    ,coordinates = mapping_file
    ,image = visium@images[[image_id]]@image
    ,scale.factors = visium@images[[image_id]]@scale.factors
    ,spot.radius = visium@images[[image_id]]@spot.radius
  )
  
  saveRDS(new.seurat.object, file = paste0(image_id, "_TLS_TCs.RDS"))
}


create_seurat_object_with_TLS_TCs(image_id = "EC1", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC2/Project_data_on_independent_components_{835982b4-294c-4c41-b539-b507d9d0fd67}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC2", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC3/Project_data_on_independent_components_{d6724aee-9a95-4eaa-bb86-71639e7bbc4f}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC3", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC3/Project_data_on_independent_components_{d6724aee-9a95-4eaa-bb86-71639e7bbc4f}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC3", corrected_mix_mat = corrected_mix_mat, visium = visium)

visium = readRDS('/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/NEW_EC_samples_ALL_VISIUM_FILES_MERGED.rds') 

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC4/Project_data_on_independent_components_{2d2d339b-b7a7-4bbb-8322-8303f81e0c16}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC4", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC5/Project_data_on_independent_components_{9d87362a-e960-44f5-ae79-399e9c3f7095}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC5", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC6/Project_data_on_independent_components_{41828a1a-38e6-4b45-98ba-c4dbf063eb0e}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC6", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC7/Project_data_on_independent_components_{fdd44f44-0ff2-4e6a-ac2f-b1936df04ab1}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC7", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC8/Project_data_on_independent_components_{28409068-263c-4cf0-b6cb-14166f40e65e}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC8", corrected_mix_mat = corrected_mix_mat, visium = visium)

corrected_mix_mat = data.frame(fread("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC9/Project_data_on_independent_components_{6558338f-d315-4b7a-87c6-e08c8ee16e40}/mixing_matrix.tsv"), row.names = 1)
colnames(corrected_mix_mat)  = gsub("consensus.independent.component.", "TC", colnames(corrected_mix_mat))
create_seurat_object_with_TLS_TCs(image_id = "EC9", corrected_mix_mat = corrected_mix_mat, visium = visium)






plot_tc_image = function(image_id, tc_to_plot, corrected_mix_mat = corrected_mix_mat, visium = visium) {
  mapping_file <- visium@images[[image_id]]@coordinates
  
  if(image_id =="OC1") {
    mapping_file = subset(mapping_file, col<185)  ### outside image region
  }
  
  common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
  mapping_file = mapping_file[common_samples,]
  corrected_mix_mat_current = corrected_mix_mat[common_samples,]
  
  # corrected_mix_mat_current$ = range01(corrected_mix_mat_current$) ## cna per image
  corrected_mix_mat_transposed = t(corrected_mix_mat_current)
  
  if(tc_to_plot==""){
    scale_fill_custom = ggplot2::scale_fill_gradient(limits = c(0,1), low = "white", high = "red" )
  } else {
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-20,20), breaks=seq(-20, 20, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-10,10), breaks=seq(-10, 10, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-15,15), breaks=seq(-15, 15, by = 2), low = "#00006b", mid = "gray", high = "#a70000", midpoint = 0, oob = scales::squish)
    scale_fill_custom = ggplot2::scale_fill_gradientn(limits = c(-6,6), breaks=seq(-6, 6, by = 3), colours = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100), oob = scales::squish)
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-15,15), breaks=seq(-15, 15, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
  }
  
  new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial")
  new.seurat.object@images$image = new(
    Class = 'VisiumV1'
    ,assay = "spatial"
    ,key = "image_"
    ,coordinates = mapping_file
    ,image = visium@images[[image_id]]@image
    ,scale.factors = visium@images[[image_id]]@scale.factors
    ,spot.radius = visium@images[[image_id]]@spot.radius
  )
  
  # pt_size = round(max(max(mapping_file$row), max(mapping_file$col))/120, 2) ### does not work
  
  p = SpatialFeaturePlot(new.seurat.object
                         , features = tc_to_plot
                         , stroke = NA
                         , slot = "counts"
                         , image.alpha = 0.2
                         # , pt.size.factor = 1.8
                         # , pt.size.factor = 1.4  ### for OC1 zoomed
                         # , pt.size.factor = pt_size
                         , alpha = 1
                         ,interactive = FALSE)&scale_fill_custom&NoAxes()
  # ,interactive = FALSE)&NoAxes()
  p
}




############################################
### iterate through all 

images_to_plot = sort(grep("OC", names(visium@images), value = T, invert = F))  

plot_list_tcs = lapply(images_to_plot, function(curr_img){
  print(paste0("### start ", curr_img))
  plot_list_tcs_curr_tc = lapply(c("text", 'raw',tcs_to_plot), function(curr_tc){
    print(paste0("### ", curr_tc))
    p = NULL
    if(curr_tc =='text') {
      p = ggplot() + geom_text(aes(label = I(curr_img), y=0.5, x=0.2), hjust = 0) + xlim(0, 1) + ylim(0,1) + theme_void()
    } else {
      if(curr_tc=='raw'){
        p = plot_raw_image(curr_img, corrected_mix_mat = corrected_mix_mat, visium = visium)
        p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
      } else {
        p = plot_tc_image(curr_img, tc_to_plot=curr_tc,corrected_mix_mat = corrected_mix_mat, visium = visium)
        p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
      }
    }
    p
  })
  plot_list_tcs_curr_tc
})
plot_list_tcs = unlist(plot_list_tcs, recursive = F)

# png(paste0('overview_spatial_top5_robust_20_color15.png'), width = 1*length(tcs_to_plot), height = 1*length(images_to_plot), units = "in", res = 1800)
# plot_grid(plotlist = plot_list_tcs, nrow=length(images_to_plot), rel_widths = c(0.6,rep(1,length(tcs_to_plot)+2)), scale = 1.1)
# dev.off()


pdf(paste0('Survival_tree_TCs_on_OC_123.pdf'), width = 3*length(tcs_to_plot), height = 2*length(images_to_plot))
plot_grid(plotlist = plot_list_tcs, nrow=length(images_to_plot), rel_widths = c(0.6,rep(1,length(tcs_to_plot)+2)), scale = 1)
dev.off()



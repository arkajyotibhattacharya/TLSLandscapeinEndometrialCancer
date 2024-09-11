############## Plotting TCs on spatial samples ####################
################ Transformation to S4 objects #########################

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")


library(RColorBrewer)
library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)

OUTPUT_DIR = '/home/iro/Projection_results/Plots'


setwd(OUTPUT_DIR)

############################################

### read in
mix_mat_bc1 = data.frame(fread("/home/iro/Projection_results/BC1/Project_data_on_independent_components_{865ba4ef-97fe-4299-bd00-c0f4f7fa2571}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc2 = data.frame(fread("//home/iro/Projection_results/BC2/Project_data_on_independent_components_{ac65f7d4-6fbc-4de3-b046-50eab769705f}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc3 = data.frame(fread("/home/iro/Projection_results/BC3/Project_data_on_independent_components_{eade2a1a-9903-49fc-864b-bf0dceca5c17}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc4 = data.frame(fread("/home/iro/Projection_results/BC4/Project_data_on_independent_components_{63ab3898-a12d-41c6-a4ba-679f91bf6884}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc5 = data.frame(fread("/home/iro/Projection_results/BC5/Project_data_on_independent_components_{c83598c9-4fd3-4f99-980a-4cf6c34ba7fb}/mixing_matrix.tsv"), row.names = 1)


# Combine the mixing matrices into one dataframe
corrected_mix_mat = as.data.frame(rbind(mix_mat_bc1, mix_mat_bc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc3))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc4))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc5))

# Remove redundant objects from the workspace
rm(mix_mat_bc2)
rm(mix_mat_bc3)
rm(mix_mat_bc4)
rm(mix_mat_bc5)

#	Read the spatial data 
visium = readRDS("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/ALL_VISIUM_FILES_MERGED.rds")

colnames(corrected_mix_mat) = gsub(pattern = "consensus.independent.component.", replacement = "TC", x = colnames(corrected_mix_mat))


tcs_to_plot = colnames(corrected_mix_mat)[1:20]


tcs_to_plot_v1 = tcs_to_plot

############################################

### plot functions


plot_raw_image = function(image_id, corrected_mix_mat = corrected_mix_mat, visium = visium) {
  
  mapping_file <- visium@images[[image_id]]@coordinates
  
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
    
    ,assay = "Spatial"
    
    ,key = "image_"
    
    ,coordinates = mapping_file
    
    ,image = visium@images[[image_id]]@image
    
    ,scale.factors = visium@images[[image_id]]@scale.factors
    
    ,spot.radius = visium@images[[image_id]]@spot.radius
    
  )
  
  
  TCs_to_plot_index = colnames(corrected_mix_mat)[1]
  
  p = SpatialFeaturePlot(new.seurat.object
                         
                         , features =TCs_to_plot_index
                         , slot = "counts"
                         
                         , stroke = NA
                         
                         , image.alpha = 1
                         
                         , pt.size.factor = 1.8
                         
                         , alpha = 0
                         
                         # , interactive = FALSE)&NoLegend()&NoAxes()
                         
                         , interactive = FALSE)&NoAxes()
  
  p
  
}



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
    
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-6,6), breaks=seq(-6, 6, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    
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
                         ,slot = "counts"
                         , stroke = NA
                         
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

images_to_plot = paste("BC", c(1:5), sep = "")



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



pdf(paste0('TC_TCGA_100EV_BC1-5.pdf'), width = 3*length(tcs_to_plot), height = 2*length(images_to_plot))

plot_grid(plotlist = plot_list_tcs, nrow=length(images_to_plot), rel_widths = c(0.6,rep(1,length(tcs_to_plot)+2)), scale = 1)

dev.off()







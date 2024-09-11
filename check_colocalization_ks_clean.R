
############################################
### create colocalization corr between TCs

setwd(dirname(rstudioapi::getSourceEditorContext()$path))   ### only in RStudio

source("../shared_functions.R")
source("../config.R")
source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

library(Seurat)
library(BiocParallel)
library(ks)


# ### adapted from:
# https://github.com/sameelab/STANN/blob/main/downstream_analysis/run_kde_estimates.R
# https://github.com/sameelab/STANN/blob/main/downstream_analysis/run_pearson_correlation.py


############################################
### read in 

# DATASET_TYPE = 'tcga'
DATASET_TYPE = 'mrna'

corrected_mix_mat = load_dataframe(paste0('../../data/spatial_transcriptomics/all_spatial_lognorm_pc1removed_',DATASET_TYPE,'_mm.tsv'))
visium = readRDS('../../data/spatial_transcriptomics/ALL_VISIUM_FILES_MERGED.rds')

BBPARAM_multicore = MulticoreParam(50)

OUTPUT_DIR = file.path('../../plots/spatial_colocalization')


############################################
### read immune and CNA

immune_ic_list = colnames(load_dataframe(paste0("../../data/immune_ica_flipped/",DATASET_TYPE,"_ica_ic_immune_flipped.tsv")))
cna_ic_list = colnames(load_dataframe(paste0('../../data/tacna/',DATASET_TYPE,'_ica_ic_tacna.tsv')))

cna_df = load_dataframe(paste0('../../data/spatial_transcriptomics/all_spatial_',DATASET_TYPE,'_cna_burden.tsv'))
corrected_mix_mat$cnaburden_01 = cna_df[rownames(corrected_mix_mat), "cna_burden_01"]
corrected_mix_mat$cnaburden = cna_df[rownames(corrected_mix_mat), "cna_burden"]

CNABURDEN_CUTOFF = 0.5


############################################
### functions

create_mapping_mm = function(image_id, corrected_mix_mat = corrected_mix_mat_01, visium = visium){
  mapping_file <- visium@images[[image_id]]@coordinates
  
  common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
  mapping_file = mapping_file[common_samples,]
  corrected_mix_mat_current = corrected_mix_mat[common_samples,]
  mapping_mm_df = cbind(mapping_file, corrected_mix_mat_current)
  mapping_mm_df
}

matrix_quantile_cutoff = function(x, q=0.5) {
  quantile_cutoff = quantile(x, q)
  x>quantile_cutoff
}



############################################
### calculate kde estimate

images_to_plot = grep("EC", names(visium@images), value = T, invert = T)  ### EC is inhouse

tcs_to_plot = c(unique(c(immune_ic_list, cna_ic_list)),"cnaburden_01")


#job::job({
  dir.create(file.path(OUTPUT_DIR, "ks_plots_neg3"), showWarnings = F)
  # dir.create(file.path(OUTPUT_DIR, "ks_plots_pos3"), showWarnings = F)
  img_tc_kernels_positive = sapply(images_to_plot, function(curr_img) {
    print(curr_img)
    curr_mapping_df_mm = create_mapping_mm(curr_img, corrected_mix_mat = corrected_mix_mat, visium = visium)
    # curr_mapping_df_mm$cnaburden_per_image = range01(curr_mapping_df_mm$cnaburden_01)
    
    OUTPUT_DIR_IMG = file.path(OUTPUT_DIR, "ks_plots_neg3", curr_img)
    # OUTPUT_DIR_IMG = file.path(OUTPUT_DIR, "ks_plots_pos3", curr_img)
    dir.create(OUTPUT_DIR_IMG, showWarnings = F)
    
    # tc_kde_list = lapply(tcs_to_plot, function(curr_tc) {
    tc_kde_list = bplapply(tcs_to_plot, function(curr_tc) {
      if(curr_tc =="cnaburden_01") {
        x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] > CNABURDEN_CUTOFF ,]
      } else {
        # mm_cutoff = 3
        # x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] > mm_cutoff ,]  ### pos only
        mm_cutoff = -3
        x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] < mm_cutoff , ]  ### neg only
      }
      
      ### at least 50 spots must be there to identify pattern
      if(nrow(x)<50){ 
        x_kde = matrix(NA, nrow = max(curr_mapping_df_mm$row), ncol = max(curr_mapping_df_mm$col))
        return(x_kde)
      } else {
        
        ### bandwith optimisation start - avoids bandwith close to 0
        h_start <- matrix(0, nrow = 2, ncol = 2)
        diag(h_start) <- 9
        
        ### Hlscv.diag uses all threads - no idea how to limit
        x_hpi = Hlscv.diag(x=x[,c("row","col")], Hstart=h_start)
        diag(x_hpi) = pmax(diag(x_hpi),3)  ### rare edge cases: sometimes bandwith optimisation is close to 0
        
        x_kde = kde(x=x[,c("row","col")], w=abs(x[,curr_tc]), H=x_hpi, xmin = c(0, 0), xmax = c(max(curr_mapping_df_mm$row), max(curr_mapping_df_mm$col)))
        
        ### plot output
        png(file.path(OUTPUT_DIR_IMG, paste0(curr_tc,".png")), width=1050, height = 400)
        par(mfrow=c(1,3))
        plot(x$row, x$col, ann=FALSE,axes=F,frame.plot=T, ylim=c(0,max(curr_mapping_df_mm$col)), xlim=c(0, max(curr_mapping_df_mm$row)))
        image(x_kde$estimate, useRaster = T, axes = FALSE)
        
        x_kde_mask = matrix_quantile_cutoff(x_kde$estimate, 0.75)
        x_kde_plot = x_kde$estimate
        x_kde_plot[!x_kde_mask] = NA
        image(x_kde_plot)
        dev.off()
        
        return(x_kde$estimate)
      }
      
    },BPPARAM = BBPARAM_multicore)
    # })
    
    names(tc_kde_list) = tcs_to_plot
    tc_kde_list
    
  }, USE.NAMES = T, simplify=F)
  
  # saveRDS(img_tc_kernels_positive, file.path(OUTPUT_DIR, "spatial_kernels_tc_positive3_obj.rds"))
  saveRDS(img_tc_kernels_positive, file.path(OUTPUT_DIR, "spatial_kernels_tc_negative3_obj.rds"))
  
#})


############################################
### calculate colocalization correlation


### calc corr all pos|neg possibilties ~8 hours
#job::job({
  img_tc_kernels = list('pos'=readRDS(file.path(OUTPUT_DIR, "spatial_kernels_tc_positive3_obj.rds")),
                        'neg'=readRDS(file.path(OUTPUT_DIR, "spatial_kernels_tc_negative3_obj.rds")) )
  curr_img=names(img_tc_kernels$pos)[1] ### just need one example
  
  par_grid = expand.grid('image_id'=names(img_tc_kernels$pos),
                         'tc_1' = names(img_tc_kernels$pos[[curr_img]]),'tc_2' = names(img_tc_kernels$pos[[curr_img]]),
                         'tc1_direction' = c('pos','neg'), 'tc2_direction' = c('pos','neg'), stringsAsFactors=F)
  
  corr_list = bplapply(1:nrow(par_grid), function(curr_row) {
    # print(curr_row)
    curr_par_grid = par_grid[curr_row,]
    
    tc_1 = img_tc_kernels[[curr_par_grid$tc1_direction]][[curr_par_grid$image_id]][[curr_par_grid$tc_1]]
    tc_2 = img_tc_kernels[[curr_par_grid$tc2_direction]][[curr_par_grid$image_id]][[curr_par_grid$tc_2]]
    
    if(all(is.na(tc_1)) | all(is.na(tc_2))) {
      curr_par_grid$pearson_corr = NA
      curr_par_grid$pearson_corr_pvalue = NA
      return(curr_par_grid)
    }
    
    tc_1_mask = matrix_quantile_cutoff(tc_1, 0.75)  ### 0.75 not like 0.5 like in orginal
    tc_2_mask = matrix_quantile_cutoff(tc_2, 0.75)
    kde_mask_common = tc_1_mask | tc_2_mask  ### union not intersect like in original
    corr_obj = cor.test(tc_1[kde_mask_common], tc_2[kde_mask_common], method="pearson")
    
    curr_par_grid$pearson_corr = corr_obj$estimate
    curr_par_grid$pearson_corr_pvalue = corr_obj$p.value
    return(curr_par_grid)
  }, BPPARAM = BBPARAM_multicore)
  corr_df = do.call('rbind',corr_list)
  save_dataframe(corr_df, file.path(OUTPUT_DIR,paste0('spatial_kernels_correlation_',DATASET_TYPE,'.tsv')))
#})



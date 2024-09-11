library(Matrix)
library(Seurat)

sample_directories <- readLines("//zkh/appdata/research/OG/202200193_TRISS/Raw_data/Raw_and_cell_ranger_processed_data/sample_directories.txt")
sample_index <- read.table("//zkh/appdata/research/OG/202200193_TRISS/Raw_data/Raw_and_cell_ranger_processed_data/sample_index.txt")
sample_directories  = sample_directories[10:14]
sample_index  = as.data.frame(sample_index[10:14,])
colnames(sample_index) = c("directory", "sample_name")

process_multiple_samples <- function(sample_directories, sample_index) {
  processed_data <- list()
  total_samples_processed <- 0
  
  for (sample_dir in sample_directories) {
    sample_name <- sample_index[sample_index$directory == sample_dir, "sample_name"]
    
    if (length(sample_name) > 0) {
      sample_name <- as.character(sample_name)
      
      matrix_dir <- file.path(sample_dir, "filtered_feature_bc_matrix/")
      barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
      features.path <- file.path(matrix_dir, "features.tsv.gz")
      matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
      
      mat_filtered <- readMM(file = matrix.path)
      
      feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
      barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
      
      colnames(mat_filtered) <- barcode.names$V1
      rownames(mat_filtered) <- feature.names$V1
      
      tissue_positions <- read.csv(file.path(sample_dir, "spatial/tissue_positions.csv"), header = TRUE)
      
      spatial_data <- Load10X_Spatial(sample_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = sample_name, filter.matrix = TRUE, to.upper = FALSE)
      
      processed_data[[sample_name]] <- list(
        mat_filtered = mat_filtered,
        feature_names = feature.names,
        barcode_names = barcode.names,
        tissue_positions = tissue_positions,
        spatial_data = spatial_data
      )
      
      total_samples_processed <- total_samples_processed + 1
    } else {
      cat("Sample name not found for directory:", sample_dir, "\n")
    }
  }
  
  cat("Total samples processed:", total_samples_processed, "\n")
  
  return(processed_data)
}
processed_samples <- process_multiple_samples(sample_directories, sample_index)

saveRDS(processed_samples, file = "F://My Drive//Projects_27092022//Spatial transcriptomics//Data//visium_rds_files//TDLN_tonsil.rds")

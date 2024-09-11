
# PACKAGE MANAGER  
pacman::p_load(pacman, dplyr, Seurat, patchwork, ggplot2, RColorBrewer, cowplot, data.table)

# load RDS file 
rds_file <- readRDS("/home/feijasomefun/TLS_PP_R_Files/Spatial_part_3/spatial_TLS_n6_integrated.rds")
# load meta_data (with TLS info)
meta_deta <- rds_file@meta.data

# select cluster 24 = the cluster with TLS 
# this shows the TLS locations for all the ECs    
meta_data_cluster_24 = meta_deta[which(meta_deta$seurat_clusters==24),]
table(meta_data_cluster_24$orig.ident)


# LOADING EC (component) directory: 1,3,4,5,6,7 
#file_path <- "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/"
EC_file_path <- c(  
  EC1 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC1/Project_data_on_independent_components_{a46a91f0-06bb-408c-9b46-ae7e1d986d66}/mixing_matrix.tsv",
  EC3 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC3/Project_data_on_independent_components_{d6724aee-9a95-4eaa-bb86-71639e7bbc4f}/mixing_matrix.tsv",
  EC4 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC4/Project_data_on_independent_components_{2d2d339b-b7a7-4bbb-8322-8303f81e0c16}/mixing_matrix.tsv",
  EC5 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC5/Project_data_on_independent_components_{9d87362a-e960-44f5-ae79-399e9c3f7095}/mixing_matrix.tsv",
  EC6 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC6/Project_data_on_independent_components_{41828a1a-38e6-4b45-98ba-c4dbf063eb0e}/mixing_matrix.tsv",
  EC7 = "/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_of_TCGA_components_100EV_on_all_EC_samples/Results/EC7/Project_data_on_independent_components_{fdd44f44-0ff2-4e6a-ac2f-b1936df04ab1}/mixing_matrix.tsv")

EC_strings <- c("EC1","EC3","EC4","EC5","EC6","EC7")


# DATA PROCESSESING: 
# load EC + metadata files + process 
# function: remove dashes and extra 1's
# load metadata for relevant EC 
f_EC_processing <- function(f_EC_file_path, f_EC_string, f_meta_data){
  
  # row.names = 1 specifies first column used as row names 
  EC <- data.frame(fread(f_EC_file_path), row.names = 1) 
  # remove dashes with 1's at end of ACATCG... 
  
  rownames_EC_without_dashes = sapply(rownames(EC), function(x){
    y = strsplit(x, "-")[[1]][1]
  })
  rownames(EC) = rownames_EC_without_dashes
  
  # line selects rows from f_meta_data where column orig.ident is equal to "EC
  meta_data_EC = f_meta_data[which(f_meta_data$orig.ident==f_EC_string),]
  # remove dashes 
  rownames_meta_deta_EC_without_dashes = sapply(rownames(meta_data_EC), function(x){
    y = strsplit(x, "-")[[1]][1]
  })
  
  rownames(meta_data_EC) = rownames_meta_deta_EC_without_dashes
  length(which(rownames(EC)%in%rownames(meta_data_EC)))
  result_list <- list(EC = EC, metadata = meta_data_EC)
  return(result_list)}

# Apply the function to all ECs (mapply for multiple arguments)
l_EC_processing <- mapply(f_EC_processing, EC_file_path, EC_strings, MoreArgs = list(f_meta_data = meta_deta), SIMPLIFY = FALSE)

# remove metadata part for descriptive stats
l_EC_without_meta <- lapply(l_EC_processing, function(ec_list) ec_list$EC)

# DESCRIPTIVE stats ECs

# function to generate distribution of:
# mean, median, IQR, Max, Min, SD, NAs per row and column/component
run_descriptive_analysis <- function(f_file, f_starting_column, f_column_path,
                                     f_row_path) {
  # load packages for subseting and boxplots
  pacman::p_load(pacman, dplyr, ggplot2, RColorBrewer, cowplot, data.table)
  
  # load mixing matrix
  df_mixing_matrix <- f_file
  
  # function histograms 
  f_summary_histograms <- function(f_data_frame) {
    colnames(f_data_frame) <- c("mean activity score", "median activity score",
                                "minimum activity score", "maximum activity score", 
                                "activity scores standard deviation", 
                                "activity scores IQR", "missing values")
    
    # Define colors for each histogram
    hist_colors <- c("#d0a1f0", "#84d9ad", "#fca368", "#59c3d4", "#fa89b5", "#a1c7f0", "#f0f0f0")
    
    theme_histo <- theme_classic() + theme(
      title = element_text(size = 11),
      axis.line = element_line(color = "#413d45")
    )
    
    # List to store individual histogram plots
    histogram_list <- list()
    
    # Loop through each column in the data frame
    for (i in seq_along(colnames(f_data_frame))) {
      column <- colnames(f_data_frame)[i]
      
      # Create the histogram plot for the current column
      histogram_plot <- ggplot(f_data_frame, aes(x = !!sym(column))) +
        geom_histogram(color = "#130f17",
                       fill = hist_colors[i]) +
        labs(title = paste("Distribution of", column), 
             x = paste(column), y = "Number of components") +
        theme_histo
      
      # Store the plot in the list
      histogram_list[[column]] <- histogram_plot
    }
    
    # Arrange the histograms into a grid using cowplot
    combined_plots <- plot_grid(plotlist = histogram_list, ncol = 2, align = "hv")
    return(combined_plots)
  }
  
  # Assuming df_descriptive_data is a data.table
  df_descriptive_stats <- setDT(df_mixing_matrix)
  
  # Define a function to calculate summary statistics for a column or row = dimension
  calculate_summary <- function(dimension) {
    c(
      Mean = mean(dimension, na.rm = TRUE),
      Median = median(dimension, na.rm = TRUE),
      Min = min(dimension, na.rm = TRUE),
      Max = max(dimension, na.rm = TRUE),
      SD = sd(dimension, na.rm = TRUE),
      IQR = IQR(dimension, na.rm = TRUE),
      NA_Count = sum(is.na(dimension))
    )
  }
  
  # COLUMNS 
  # Apply the function to each column using sapply
  df_component_summary <- data.table(t(sapply(df_descriptive_stats[, f_starting_column:ncol(df_descriptive_stats), with = FALSE], calculate_summary)))
  output_file_path_comp <- f_column_path
  # Run distributions
  output_plots_comp <- f_summary_histograms(df_component_summary)
  # Save the combined plot to a PDF file
  ggsave(output_file_path_comp, plot = output_plots_comp, width = 12, height = 8)
  
  # ROWS 
  # Select the columns of interest
  selected_columns <- df_descriptive_stats[, f_starting_column:ncol(df_descriptive_stats), with = FALSE]
  # change y variable:
  
  # Apply the function to each row (1) using data.table's `apply`
  df_row_summary <- data.table(t(apply(selected_columns, 1, calculate_summary)))
  
  output_file_path_row <- f_row_path
  output_plots_row <- f_summary_histograms(df_row_summary)
  ggsave(output_file_path_row, plot = output_plots_row, width = 12, height = 8)
}

# No NA values 
sum(is.na(l_EC_without_meta))

# function to apply descriptive stats function to each EC 
f_apply_descriptive_stats <- function(EC_data, f_EC_strings) {
  
  generate_output_paths <- function(EC_names, base_path = "/home/feijasomefun/TLS_PP_R_Files/Spatial_part_3/Spatial_plots/") {
    lapply(EC_names, function(ec_name) {
      list(
        components = paste0(base_path, ec_name, "_components.pdf"),
        rows = paste0(base_path, ec_name, "_rows.pdf")
      )
    })
  }
  
  # Generate output paths
  output_paths <- generate_output_paths(EC_names = f_EC_strings)
  
  # lapply to iterate over each EC dataframe
  lapply(seq_along(EC_data), function(i) {
    EC_name <- names(EC_data)[i]
    
    # Extract output paths for the specific EC
    output_path_components <- output_paths[[i]]$components
    output_path_rows <- output_paths[[i]]$rows
    
    # Run descriptive analysis for the current EC
    run_descriptive_analysis(
      f_file = EC_data[[i]],
      f_starting_column = 1,
      f_column_path = output_path_components,
      f_row_path = output_path_rows
    )
  })
}
results_descriptive_stats <- f_apply_descriptive_stats(
  EC_data = l_EC_without_meta, f_EC_strings = EC_strings)

# RETRIEVING SIGNIFICANT COMPONENTS: 39 
# directory for significant components  
t_significiant_components <- fread("/home/feijasomefun/table_comparing_multiple_testing.csv")

#store component names
t_component_names <- t_significiant_components$c_component_number
# Replace spaces with dots and convert to character
component_names_select <- gsub(" ", ".", as.character(t_component_names))
print(component_names_select)

# SUBSET SIGNIFICANT COMPONENTS 
f_subset_components_EC <- function(f_data, f_columns) {
  dt_subset <- f_data[, f_columns, drop = FALSE]
  return(dt_subset)
}

# Apply the function to each element in the list
l_subset_components_EC <- lapply(l_EC_without_meta, f_subset_components_EC, f_columns = component_names_select)

# CHECK FOR TLS 

#EC <- data.frame(fread(f_EC_file_path), row.names = 1)
# load data EC log normalized PC 1 removed 

# Function to read data for a given EC
f_read_EC_data <- function(EC_number) {
  # sprintf %s is the placeholder for the EC numbers 
  file_path <- sprintf("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/spatial_counts_pc1removed/ncbi/EC%s_lognorm_pc1removed_ncbi.tsv", EC_number)
  ec_data <- data.frame(fread(file_path), row.names = 1)  # first column is rownames
  
  # Apply substring removal function to exclude the first row
  rownames_EC_without_dashes <- sapply(rownames(ec_data)[-1], function(x) {
    strsplit(x, "-")[[1]][1]
  })
  
  # Update row names, including the first row name
  rownames(ec_data) <- c(rownames(ec_data)[1], rownames_EC_without_dashes)
  return(ec_data)
}
# List of EC numbers
EC_numbers <- c("1", "3", "4", "5", "6", "7")
# Read data for each EC, apply substring removal, and store in a list
l_EC_data <- lapply(EC_numbers, f_read_EC_data)

# for each EC: 
# from metadata create a list with rownames where seuret cluster = 24 
# in the l_subset_components check the rownames 
# if the rownames is any of the rownames in the list
# create a new column and add True to column 

# remove components part --> only metadata 
l_EC_processed_meta <- lapply(l_EC_processing, function(ec_list) ec_list$metadata)

f_identify_TLS <- function(f_metadata, f_components, f_matching_name) {
  selected_rows <- f_metadata$seurat_cluster == 24
  # select the rows from the metadata 
  selected_row_names <- rownames(f_metadata[selected_rows, ])#, drop = FALSE]) # for as dataframe
  print(selected_row_names)
  
  # Check whether rownames of component are in selected row names 
  match_rows <- rownames(f_components) %in% selected_row_names
  # add match to new matching colummn in components dataframe 
  f_components[f_matching_name] <- match_rows
  
  return(f_components)
}

df_components_matched <- mapply(f_identify_TLS, 
                                l_EC_processed_meta, 
                                l_subset_components_EC, 
                                MoreArgs = list(f_matching_name = "TLS_presence"),
                                SIMPLIFY = FALSE)


check_components <- l_subset_components_EC[["EC1"]]
check <- f_identify_TLS(l_EC_processing[["EC1"]][["metadata"]],check_components,"TLS_presence")
print(check)



# check 
# for l_EC_data 
# use meta_data subset from processing files  
# make a list of the values that have seuret cluster 24 
# use arkos trick to have the first column be the row names 
# check if any of the values in the list match with the value in the first column



# EC1 
# meta_data_EC1
# TLS cluster 24, what is the cluster ?
# seperate files and load projected mixing matrix 
#  1) 6 mixing matrix weights; 6*7 graphs (combine; cowplot) 14 PDFs,  
# (function for columns and rows)
# 6 EC samples explore descriptive statistics 
# everything per sample 
# script; seperate metadata for EC1 
# from metadata TLS yes or NO (one column, yes when it is 24, no when not 24)
# 1 and 0
# load mixing matrix weights of all the components for EC1 (EC1 file)
# make a subset for 39 known components 
# EC1 39 columns add the yes or no TLS column 
# match the rows of the EC1 (TLS location) with component
# uncomment some things of the script 
#kernal smoothing
# EC 

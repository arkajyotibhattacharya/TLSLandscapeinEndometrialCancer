
# STEP 1: IMPORTING 
# import all the datasets

# load the data.table package for the fread function 
library(data.table)
library(parallel)
# load the tsv dataset (9709)

# original Table of mixing matrix file 
T_mixing_matrix <- fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Students/Feija\ Somefun/Data/ica_flipped_mixing_matrix_consensus.tsv")

# variable for the dataframe with TLS and non TLS identification of samples
T_TLS_identification <- fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Students/Feija\ Somefun/Data/TCGA_UCEC_metadata.txt")


common_samples = intersect(substr(T_mixing_matrix$V1,1,12), T_TLS_identification$`TCGA Participant Barcode`)

T_mixing_matrix = T_mixing_matrix[which(substr(T_mixing_matrix$V1,1,12)%in%common_samples),]
T_TLS_identification = T_TLS_identification[which(T_TLS_identification$`TCGA Participant Barcode`%in%common_samples),]



# load pacman and dplyr package for subseting
pacman:: p_load(pacman, dplyr)


# STEP 2: SUBSETING 
# subsetting the TLS and NONE TLS data  


# SIMPLIFIED TABLES 

# this is a simple version with only barcodes and TLS presence  
T_TLS_ANY <-select(filter(T_TLS_identification, TLS == "ANY"), c("TCGA Participant Barcode", TLS))  
# there are 30 TLS participant codes with ANY TLS 

T_TLS_NONE <-select(filter(T_TLS_identification, TLS == "NONE"), c("TCGA Participant Barcode", TLS))  
# there are 256 participant codes with NONE TLS


# STEP 3:
# match the TLS subsets with the mixing matrix


# Create a copy of the mixing matrix with a different memory address
c_mixing_matrix <- data.frame(T_mixing_matrix)

# TLS MATCHED dataframe 
# Create new data frame (df) with copy of mixing matrix and a new column TLS_match
df_TLS_bcode_check <- cbind(c_mixing_matrix, TLS_match = NA)

# Perform the matching 
df_TLS_bcode_check$TLS_match <- substr(df_TLS_bcode_check$V1, 1, 12) %in% T_TLS_ANY$"TCGA Participant Barcode"

# code breakdown:
# %in% operator to perform element-wise comparisons 
# element wise comparison (the TLS barcode basic) should check all the value in the df_bcode_check 
# first 12 characters of df_TLS_bcode_check$V1 and the barcodes in T_TLS_ANY$TCGA Participant Barcode
# The result is a logical vector that is used to complete the df_TLS_bcode_check$TLS_match column.

# Create the TLS_matched dataframe using filter function: 26 entries 
df_TLS_matched_2 <- filter(df_TLS_bcode_check, TLS_match)


#NONE tls MATCHED dataframe 
# Create new data frame (df) with copy of mixing matrix and a new column TLS_match
df_NONE_bcode_check_2 <- cbind(c_mixing_matrix, NONE_match = NA)

# Perform the matching without using explicit loops
df_NONE_bcode_check_2$NONE_match <- substr(df_NONE_bcode_check_2$V1, 1, 12) %in% T_TLS_NONE$"TCGA Participant Barcode"

# create the new NONE TLS matched dataframe: 192 entries 
df_NONE_matched_2 <- filter(df_NONE_bcode_check_2, NONE_match)


# STEP 4 Mann Whitney test 

# for each of the (9709) components compare the activity scores of TLS (26 entries) with NON TLS (192 entries) 




# DATA FRAME COLUMNS FOR MANN WHITNEY PROCESSING 

# create new data frame
# create column variable for storing component number as a vector 
# where the length of the vector is equal to the number of components = 9709 
# first column 1 = name last column 9711 = matching info

new_column <- vector("numeric", 9709)

column_list <- list(
  c_component_number = new_column,
  # create column for p value 
  c_p_value = new_column,
  # column negative - log10 of p value
  c_neg_log_p = new_column,
  # create column vector for median difference 
  c_median_difference = new_column, 
  # column for sign median difference 
  c_sign_median_difference = new_column,
  # column for processed p 
  c_adj_p = new_column,
  # column for correction of multiple testing values 
  c_correction = new_column,
  c_significant_p = new_column
)

# create a variable for corrected last column ID  
# column 9711 contains matching information 

last_col <- ncol(df_TLS_matched_2)
print(last_col)
update_last_col <- last_col - 1
print(update_last_col)

# dataframe to process the folowing calculation 
# - Log10(P-value ) * sign (median TLS ANY - median TLS None) 


# new dataframe for results of data processing   
df_results_MW <- data.frame(column_list) 


mwu_all_tcs = function(x, df_TLS_matched_2_local, df_NONE_matched_2_local){
  # x represents the column, x stars from column 2 (column 1 = barcodes) 
  # ends at column 9710 = independent component 9709, column 9711 contains matching info  
  # print(x)
  # [x] is column while [[x]] is to create a vector of the column x
  # wilxoc needs a vector as input 
  WR_TLS_col <- df_TLS_matched_2_local[[x]]
  WR_NONE_col <- df_NONE_matched_2_local[[x]]
  
  # Wilcox mann whitney test 
  wr_test_results <- wilcox.test(WR_TLS_col, WR_NONE_col)
  
  # Median of TLS and NONE column 
  median_TLS <- median(df_TLS_matched_2_local[[x]])
  median_NONE <- median(df_NONE_matched_2_local[[x]])
  
  # Median difference 
  median_dif <- median_TLS - median_NONE
  # Sign of median difference (know which side there is a difference)
  sign_median_dif <- sign(median_dif) 
  
  # print the column name of the dataset of index x 
  # print(paste("TLS component", colnames(df_TLS_matched_2_local)[x]))
  # print(paste("NONE component", colnames(df_NONE_matched_2_local)[x]))
  # print(wr_test_results)
  
  # add test results to dataframe Mann whitney results 
  df_results_MW$c_component_number[x-1] <- colnames(df_TLS_matched_2_local)[x]
  
  # TC as column 1 P value as column 2 
  df_results_MW$c_p_value[x-1] <- wr_test_results$p.value
  # -log10 of p value 
  neg_log_10_p <- -log10(wr_test_results$p.value)
  df_results_MW$c_neg_log_p[x-1] <- neg_log_10_p
  
  # add median difference to median difference column 
  df_results_MW$c_median_difference[x-1] <- median_dif
  
  # add sign of median difference to column
  df_results_MW$c_sign_median_difference[x-1] <- sign_median_dif
  
  # - Log10(P-value ) * sign (median TLS ANY - median TLS None) 
  full_calc <- neg_log_10_p * sign_median_dif
  return(full_calc)
}


start.time <- Sys.time()

# Your R code here
result <- sapply(c(2:update_last_col), function(x){mwu_all_tcs(x,df_TLS_matched_2_local = df_TLS_matched_2, df_NONE_matched_2_local = df_NONE_matched_2 )})

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

original_log_transformed_p_value_with_sign = result

multiple_testing_mwu = function(i){
  
  perm_data = T_TLS_identification
  set.seed(12345+i)										                                                    #makes sure the RNG in "sample" does exactly the same when rerunning the code (so, purely for reproducibility)
  perm_data$TLS = sample(perm_data$TLS)
  # this is a simple version with only barcodes and TLS presence  
  T_TLS_ANY <-select(filter(perm_data, TLS == "ANY"), c("TCGA Participant Barcode", TLS))  
  # there are 30 TLS participant codes with ANY TLS 
  
  T_TLS_NONE <-select(filter(perm_data, TLS == "NONE"), c("TCGA Participant Barcode", TLS))  
  # there are 256 participant codes with NONE TLS
  
  
  # STEP 3:
  # match the TLS subsets with the mixing matrix
  
  # TLS MATCHED dataframe 
  # Create new data frame (df) with copy of mixing matrix and a new column TLS_match
  df_TLS_bcode_check <- cbind(c_mixing_matrix, TLS_match = NA)
  
  # Perform the matching 
  df_TLS_bcode_check$TLS_match <- substr(df_TLS_bcode_check$V1, 1, 12) %in% T_TLS_ANY$"TCGA Participant Barcode"
  
  # code breakdown:
  # %in% operator to perform element-wise comparisons 
  # element wise comparison (the TLS barcode basic) should check all the value in the df_bcode_check 
  # first 12 characters of df_TLS_bcode_check$V1 and the barcodes in T_TLS_ANY$TCGA Participant Barcode
  # The result is a logical vector that is used to complete the df_TLS_bcode_check$TLS_match column.
  
  
  # Create the TLS_matched dataframe using filter function: 26 entries 
  df_TLS_matched_2 <- filter(df_TLS_bcode_check, TLS_match)
  
  
  #NONE tls MATCHED dataframe 
  # Create new data frame (df) with copy of mixing matrix and a new column TLS_match
  df_NONE_bcode_check_2 <- cbind(c_mixing_matrix, NONE_match = NA)
  
  # Perform the matching without using explicit loops
  df_NONE_bcode_check_2$NONE_match <- substr(df_NONE_bcode_check_2$V1, 1, 12) %in% T_TLS_NONE$"TCGA Participant Barcode"
  
  # create the new NONE TLS matched dataframe: 192 entries 
  df_NONE_matched_2 <- filter(df_NONE_bcode_check_2, NONE_match)
  
  
  # STEP 4 Mann Whitney test 
  
  # for each of the (9709) components compare the activity scores of TLS (26 entries) with NON TLS (192 entries) 
  
  
  # DATA FRAME COLUMNS FOR MANN WHITNEY PROCESSING 
  
  
  
  # create a variable for corrected last column ID  
  
  last_col <- ncol(df_TLS_matched_2)
  print(last_col)
  update_last_col <- last_col - 1
  print(update_last_col)
  
  # dataframe to process the folowing calculation 
  # - Log10(P-value ) * sign (median TLS ANY - median TLS None) 
  
  
  start.time <- Sys.time()
  
  # Your R code here
  result <- sapply(c(2:update_last_col), function(x){mwu_all_tcs(x,df_TLS_matched_2_local = df_TLS_matched_2, df_NONE_matched_2_local = df_NONE_matched_2 )})
  
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(time.taken)
  
  return(result)
  
}


permutation_log_transformed_p_value_with_sign_v1 = multiple_testing_mwu(1)

permutation_log_transformed_p_value_with_sign_v2 = multiple_testing_mwu(2)

plot(permutation_log_transformed_p_value_with_sign_v1, permutation_log_transformed_p_value_with_sign_v2)

nPerm = 10000
FDR = 0.05
conf_level = 0.8

time1 = proc.time()[3]
no_cores = 10
cl <- makeCluster(no_cores, type = "FORK")

outputs_combined = parLapply(cl, 1:nPerm,multiple_testing_mwu)
# outputs_combined = parLapply(cl, 1:20,multiple_testing_survival)
stopCluster(cl)
print((proc.time()[3] - time1)/60)

adjusted_Pval_combined = matrix(NA,length(outputs_combined[[1]]), nPerm)

for(i in 1:nPerm)
{
  adjusted_Pval_combined[,i] = outputs_combined[[i]]
  
}

# Calculate MVP threshold
time1=Sys.time()
original_mlogpval_ordered = sort(abs(original_log_transformed_p_value_with_sign),decreasing=TRUE )			#becomes a sorted list, no longer in the DF
adjusted_Pval_combined_mlogpval = abs(adjusted_Pval_combined)
adjusted_Pval_combined_mlogpval_v1 = apply(abs(adjusted_Pval_combined_mlogpval),2,function(x) sort(x,decreasing=TRUE))      #Should remain a DF in this step. Watch out, any problems in the earlier steps (eg. warnings) may introduce NA's, forcing this into layered list that cannot be processed in the subsequent script.

cutoffs = array(0,nPerm)									#creates a list including all cutoff results
for(j in 1:nPerm)
{
  for(row_number in 1:dim(adjusted_Pval_combined_mlogpval_v1)[1])
  {
    if(length(which(adjusted_Pval_combined_mlogpval_v1[,j]>original_mlogpval_ordered[row_number]))/row_number>=FDR)
    {
      cutoffs[j] = original_mlogpval_ordered[row_number]				#Corrected: now takes p-value from original p-values instead of permutation-p-values
      break
    }
  }
}												

#Calculate cutoff based on confidence level and use this to create a results matix with clear FDR TRUE/FALSE based on cutoff
list_of_cutoffs = quantile(cutoffs, conf_level)						#uses this list (in combination with the confidence level to mark a single cutoff value), works as expected
output = ifelse(abs(original_log_transformed_p_value_with_sign)<list_of_cutoffs,"False","True")			#checks every -log10pvalue in the original data for the cutoff, only values more extreme should become "TRUE"
results_final = as.data.frame(cbind(original_log_transformed_p_value_with_sign, output))
results_final$original_log_transformed_p_value_with_sign = as.numeric(results_final$original_log_transformed_p_value_with_sign)

colnames(results_final)[2]=c(paste("MVP_reject_H0", "nPerm=", nPerm, "FDR=", FDR, "conf=", conf_level))



# create new dataframe with only significant p values 
df_sig_p <- select(filter(df_results_MW, c_significant_p == "TRUE"), c_component_number, c_p_value, c_adj_p )
View(df_sig_p)


# final processed mann whitney data frame with the component number and the final calculation
# column final calculation  = - Log10(P-value) * sign (median TLS ANY - median TLS None) 
# Sign of (median TLS – NONE median) = +1 when TLS is more active than NONE
# sign of (median TLS - NONE median) = -1 NONE is more active than TLS 

df_mw_processed_results <- select(df_results_MW,"c_p_value","c_adj_p")
View(df_mw_processed_results)

# export calculation dataframe to csv
write.csv(df_results_MW,file="/home/feijasomefun/Downloads/df_results_MW.csv",fileEncoding = "UTF-8")

# export final dataframe to CSV 
write.csv(df_mw_processed_results,file="/home/feijasomefun/Downloads/df_raw_calculations.csv",fileEncoding = "UTF-8")





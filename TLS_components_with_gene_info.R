library(data.table)
library(readxl)
ica= data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Students/Feija\ Somefun/Data/ica_flipped_independent_components_consensus.tsv"), row.names = 1)

TLS_components = as.data.frame(read_excel("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Students/Feija\ Somefun/Results/Multiple_testing_results_TLS.xlsx"))
TLS_components_numbers = paste("consensus.independent.component.",gsub("TC_", "", TLS_components$Components), sep = "")
ica = ica[,TLS_components_numbers]

genomic_mapping = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Databases/Genomic\ mapping\ files/Entrezid_mapping_using_org_Hs_eg_db_17012024.txt"), row.names = 1)

common_genes = intersect(rownames(genomic_mapping), rownames(ica))

ica = ica[common_genes,]
genomic_mapping = genomic_mapping[common_genes,]
ica_combined = as.data.frame(cbind(genomic_mapping,ica))
ica_combined = as.data.frame(cbind(rownames(ica_combined),ica_combined))
colnames(ica_combined)[1] = "Entrezid"
colnames(ica_combined) = gsub("\\.", "_", colnames(ica_combined))


write.table(ica_combined, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Spatial\ transcriptomics/Results/ica_flipped_independent_components_consensus_100EV_TCGA_TLS_associated_with_gene_info.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

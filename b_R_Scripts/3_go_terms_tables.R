# ------------------------------------------- GO TERMS OILS VS PLACEBO ---------------------------------------------#

# GO analysis of the GC-JD-7904 project. 

# This small script prepares the tables for the GO analysis.

# ---> LIBRARY UPLOAD:

library(plyr)
library(dplyr)

# ---> DATA LOAD: 

# Results from DGE: 
# p = paired samples, b = batch correction, f = filtered, wl = without low quality samples, no = No outliers based
# on PCA and oPLSDA. 

edge_p_b_f_no <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/edge_filter_NoOutlier_paired.txt",
  header = TRUE,
  sep = "\t")

edge_p_b_f <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/table_data_filter_paired.txt",
  header = TRUE,
  sep = "\t")

limma_p_b_f <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/limma_paired_filter.txt",
  header = TRUE,
  sep = "\t")

edge_b_f_wl <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/table_data_filter_nopaired_wl.txt",
  header = TRUE,
  sep = "\t")

limma_b_f_wl <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/limma_filter_nopaired_wl.txt",
  header = TRUE,
  sep = "\t")

# Ensembled names:

ensembl <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/ensembl_id.txt",
  header = TRUE,
  sep = "\t")

# Human GO terms: 

human_go <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/3_go_terms_tables/human_go_terms.txt",
  header = TRUE,
  sep = "\t")

# ---> DATA MANIPULATION: 

# Get only the gene names and p_values: 

edge_p_b_f_no <- data.frame(gene = rownames(edge_p_b_f_no),
                         p_value = edge_p_b_f_no$FDR)

edge_p_b_f <- data.frame(gene = rownames(edge_p_b_f),
                         p_value = edge_p_b_f$PValue)

limma_p_b_f <- data.frame(gene = rownames(limma_p_b_f),
                         p_value = limma_p_b_f$P.Value)

edge_b_f_wl <- data.frame(gene = rownames(edge_b_f_wl),
                         p_value = edge_b_f_wl$PValue)

limma_b_f_wl <- data.frame(gene = rownames(limma_b_f_wl),
                         p_value = limma_b_f_wl$P.Value)

# Keep only one row for gene: 

human_go_final <- ddply(human_go, .(Gene.name), summarize,
                        GO.term.accesion = paste(unique(GO.term.accession), collapse = ","))
colnames(human_go_final) <- c("gene", "GO_term")
human_go_final <- human_go_final[!human_go_final$GO_term == "", ]

# Remove extra comas (,) and add spaces: 

human_go_final$GO_term <- gsub(pattern = "^,", replacement = "", human_go_final$GO_term)
human_go_final$GO_term <- gsub(pattern = ",$", replacement = "", human_go_final$GO_term)
human_go_final$GO_term <- gsub(pattern = ",", replacement = ", ", human_go_final$GO_term)

# Merge ensembl_id with table with genes: 

ensembl_table <- merge(edge_p_b_f_no, ensembl)
ensembl_table <- ensembl_table[ensembl_table$p_value < 0.05, ]
ensembl_table$p_value <- NULL 
ensembl_table$gene_id <- gsub(pattern = "\\.\\d+", "", ensembl_table$gene_id)


# ---> SAVE TABLES: 

write.table(edge_p_b_f_no, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/edge_paired_filter_NOutliers.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(edge_p_b_f, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/edge_paired_filter.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(limma_p_b_f, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/limma_paired_filter.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(edge_b_f_wl, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/edge_paired_filter_wl.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(limma_b_f_wl, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/limma_paired_filter_wl.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(human_go_final, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/human_go_final.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

write.table(ensembl_table,
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/3_go_terms_tables/edge_paired_filter_NOutliers_gene_list.txt",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)


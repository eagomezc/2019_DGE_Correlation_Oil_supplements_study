# ----------------------------------------- GENE ONTOLOGY ENRICHMENT ANALYSIS  -------------------------------------------#

# The gene ontology terms represents gene product properties. These terms covers three domains: Cellular Component, molecular,
# function and biological process. Each term is related with a general or very specific function in the system. 

# The GO ontology is structured as a directed acyclic graph, and each term has defined relationships to one or more other
# terms in the same domain, and sometimes to other domains. The GO vocabulary is designed to be species-neutral, and 
# includes terms applicable to prokaryotes and eukaryotes, single and multicellular organisms.

# With a DGE analysis we can find the genes that are down or upregulated between different conditions; the GO enrichment 
# analysis allows to understand a little more if there is a connection between all the differential expressed genes and 
# what are their functions in the organism. 

# ---> LIBRARIES: 

library(topGO)
library(devtools)
library(ggplot2)
library(ggpubr) 

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/hhy270/Documents/GitHub/2019_DGE_Correlation_Oil_supplements_study/a_Data/2_GO_analysis_(TopGO)/"
output <- "C:/Users/hhy270/Documents/GitHub/2019_DGE_Correlation_Oil_supplements_study/c_Expected_Outputs/2_GO_analysis_(TopGO)/"

# !!!! IMPORTANT: For this script to work the READ COUNTS FILE has to be called: 2_GO_DGE_p_values_list.txt
# !!!! IMPORTANT: For this script to work the SAMPLES INFO FILE has to be called: 2_GO_Human_go_terms.txt

#---> DATA LOAD: 

# Open the file with the genes and the Adjust p values (BH correction).

# The tab delimited table consists in two columns "gene" and "p_value", showing all the genes and their respective Adj p
# values obtained from the DGE. 
# See a_Data/2_GO_analysis_(TopGO)/2_GO_DGE_p_values_list.txt  

genes_p_value <- read.table(
  file = paste(input, "2_GO_DGE_p_values_list.txt", sep = ""),
  header = TRUE,
  sep = "\t")

# Open the file with the genes and their related GO terms for the organism of interest (e.g: Homo Sapiens).

# The tab delimited table consists in two columns "gene" and "GO_term", showing all the genes for a  genome and their 
# associated GO terms (separated by colons ","; e.g: GO:0005576, GO:0005615, GO:0070062,...).

## Step One
## GO annotations
gene_to_go_mapping_file <- "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/4_go_terms/human_go_final.txt"

## file of significant genes (2 column file: i.e. gene id and pvalue) file:
de_genes_file <- paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/4_go_terms/", pv_file, ".txt", sep = "")

output_directory <- paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/4_go_terms/", pv_file, "/", sep = "")

## Read in input file:
de_genes <- read.table(de_genes_file, header = TRUE)

colnames(de_genes) <- c("locus", "p_value")

## Read in GO annotations:  
gene_to_go_mapping <- readMappings(file = gene_to_go_mapping_file)

## Convert into topgo's genelist format:
topgo_genelist        <- de_genes$p_value
names(topgo_genelist) <- de_genes$locus

for (go_category in c("BP", "MF", "CC")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = function(x) { # fails to run without this
                      return(x <= 0.05)
                    },
                    allGenes    = topgo_genelist,
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 5) # Modify to reduce/increase stringency.
  
  # STEP THREE
  
  ## Calculate ks test using 'weight01' algorithm:
  result_weight_ks <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "ks")
  
  
  ## Combine results from statistical tests:
  result_weight_output <- GenTable(object       = my_go_data,
                                   weight_ks    = result_weight_ks,
                                   orderBy      = "weight_ks",
                                   numChar      = 42, # Get not full GO terms names
                                   topNodes     = length(score(result_weight_ks)))
  
  ## Subset calls with significance higher than expected:
  result_weight_output_sig <- subset(x      = result_weight_output,
                                     subset = (weight_ks <= 0.05))
  
  ## Correct for multiple testing:
  result_weight_output_sig$weight_ks_adjusted <- p.adjust(p      = result_weight_output_sig$weight_ks,
                                                          method = c("BH"),
                                                          n      = length(result_weight_output_sig$weight_ks))
  
  ## Get the significance genes id for each Go Term
  sig_genes <- de_genes[de_genes$p_value <= 0.05, ]
  genes <- lapply(result_weight_output_sig$GO.ID, function(x) as.character(unlist(genesInTerm(object = my_go_data, whichGO = x))))
  go_genes <- lapply(genes, function(x) intersect(x, sig_genes$locus))
  sig_go_genes <- lapply(go_genes, function(x) paste(x, collapse = ","))
  
  result_weight_output_sig$genes_id <- unlist(sig_go_genes)
  
  result_weight_output_sig$Term <- gsub(" [a-z]*\\.\\.\\.$", "", result_weight_output_sig$Term)
  result_weight_output_sig$Term <- gsub("\\.\\.\\.$", "", result_weight_output_sig$Term)
  
  result_weight_output_sig <- result_weight_output_sig[!(result_weight_output_sig$genes_id == ""), ]
  
  ## Print out all GO terms of interest:
  pdf(file = paste(output_directory,"/", go_category, "_distribution.pdf", sep = ""))
  for (go_term in 1:nrow(result_weight_output_sig)){
    print(showGroupDensity(object  = my_go_data,
                           whichGO = result_weight_output_sig[go_term, "GO.ID"],
                           ranks   = TRUE,
                           rm.one  = FALSE))
  }
  dev.off()
  ## Write to output:
  write.table(x         = result_weight_output,
              file      = file.path(output_directory,
                                    paste(go_category, "no_sig_ks_shortname.tsv", sep = "_")),
              row.names = FALSE,
              sep       = "\t")
  
  write.table(x         = result_weight_output_sig,
              file      = file.path(output_directory,
                                    paste(go_category, "sig_ks_shortname.tsv", sep = "_")),
              row.names = FALSE,
              sep       = "\t")
}

# }
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

# !!!! IMPORTANT: For this script to work the GENES P VALUE has to be called: 2_GO_DGE_p_values_list.txt
# !!!! IMPORTANT: For this script to work the GENOME GO TERMS has to be called: 2_GO_Human_go_terms.txt

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
# See a_Data/2_GO_analysis_(TopGO)/2_GO_Human_go_terms.txt

genome_GO_terms <- readMappings(file = paste(input, "2_GO_Human_go_terms.txt", sep = ""))

#---> DATA PREPARATION: 

# Transform genes and p values table in the genelist format requieres for TopGO:

topgo_genes <- genes_p_value$p_value
names(topgo_genes) <- genes_p_value$gene

#---> GENE ONTOLOGY ANALYSIS: 

# The gene ontology analysis runs for all the GO categories: Biological process, Molecullar Function and Cellular
# Component. 

for (go_category in c("BP", "MF", "CC")) {

  # Make the GOdata object:
  
  go_object  <- new("topGOdata",  # Name of the object 
                    description = paste("GOtest", go_category, sep = "_"), # Description 
                    ontology = go_category, # Which category is analysing
                    geneSel = function(x) {
                      return(x <= 0.05) # Fails to run without this even if it's not need it
                      }, 
                    allGenes = topgo_genes,
                    gene2GO = genome_GO_terms,
                    annot = annFUN.gene2GO,
                    nodeSize = 5) # Modify to reduce/increase stringency (Each GO term has to be associated with
                                  # at least 5 genes)
  
  # Kolmogorov-Smirnov test is used to identified the enriched GO terms: 
  
  ks_results <- runTest(object = go_object, algorithm = "weight01", statistic = "ks")
  
  # Produce the table with the results:
  
  ks_results_table <- GenTable(object = go_object, 
                             weight_ks = ks_results, 
                             orderBy = "weight_ks", 
                             numChar = 42, # Get not full GO terms names but the first _ characters. 
                             topNodes = length(score(ks_results)))
  
  # Separate the results statistical significance (ks test <= 0.05):
  
  ks_results_sig <- subset(ks_results_table, subset = (weight_ks <= 0.05))
  
  # Adjust p value (FDR) for multiple comparisons. In this case using the BH correction:
  
  ks_results_sig$weight_ks_adjusted <- p.adjust(p = ks_results_sig$weight_ks, method = c("BH"),
                                                n = length(ks_results_sig$weight_ks))

  # Get the significance genes id for each Go Term:
  
  sig_genes <- genes_p_value[genes_p_value$p_value <= 0.05, ]
  genes <- lapply(ks_results_sig$GO.ID, function(x) as.character(unlist(genesInTerm(object = go_object, whichGO = x))))
  go_genes <- lapply(genes, function(x) intersect(x, sig_genes$gene))
  sig_go_genes <- lapply(go_genes, function(x) paste(x, collapse = ","))
  
  ks_results_sig$genes_id <- unlist(sig_go_genes)
  
  ks_results_sig$Term <- gsub(" [a-z]*\\.\\.\\.$", "", ks_results_sig$Term)
  ks_results_sig$Term <- gsub("\\.\\.\\.$", "", ks_results_sig$Term)
  
  ks_results_sig <- ks_results_sig[!(ks_results_sig$genes_id == ""), ]
  
  # ---> OUTPUT:
  
  # Print all GO terms of interest:
  
  pdf(file = paste(output, "2_", go_category, "_density_plot.pdf", sep = ""))
  
  for (go_term in 1:nrow(ks_results_sig)){
    
    print(showGroupDensity(object  = go_object,
                           whichGO = ks_results_sig[go_term, "GO.ID"],
                           ranks   = TRUE,
                           rm.one  = FALSE))
  }
  
  dev.off()
  
  # Output tables: 
  
  write.table(ks_results_sig,
              file = paste(output, "2_", go_category, "_significant_enriched_GO_terms.tsv", sep = ""),
              row.names = FALSE,
              sep  = "\t")
}

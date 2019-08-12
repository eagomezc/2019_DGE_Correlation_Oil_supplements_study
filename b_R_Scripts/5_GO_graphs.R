# The following script takes in three tab-delimited text files output from the R package, topGO, combines the gene ontologies
# and plots as a barplot. 

## Load libraries:
libraries <- c("devtools", "ggplot2", "ggpubr")
for (lib in libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkgs = lib)
    library(lib, character.only = TRUE )
  }
}

## Density check table:
density_check <- read.table(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/density check.tsv",
                            header = TRUE,
                            sep = "\t")
density_check$BP <- gsub(pattern = " ", replacement = "", density_check$BP)
density_check$CC <- gsub(pattern = " ", replacement = "", density_check$CC)
density_check$MF <- gsub(pattern = " ", replacement = "", density_check$MF)

## Read table in for combined plot:
bp_go_terms_df <- read.csv(file="C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/BP_sig_ks_shortnames.tsv",
                           header=TRUE, sep="\t")

bp_go_terms_df <- bp_go_terms_df[bp_go_terms_df$GO.ID %in% density_check$BP, ]

## Add status column:
bp_go_terms_df$status <- "Biological Processes"

mf_go_terms_df <- read.csv(file="C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/MF_sig_ks_shortnames.tsv",
                           header=TRUE, sep ="\t")

mf_go_terms_df <- mf_go_terms_df[mf_go_terms_df$GO.ID %in% density_check$MF, ]

## Add status column:
mf_go_terms_df$status <- "Molecular Function"

cc_go_terms_df <- read.csv(file="C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/CC_sig_ks_shortnames.tsv",
                           header=TRUE, sep = "\t")

cc_go_terms_df <- cc_go_terms_df[cc_go_terms_df$GO.ID %in% density_check$CC, ]

## Add status column:
cc_go_terms_df$status <- "Cellular Component"

## Combine terms and filter for significance:
combined_go_terms_df<- rbind(bp_go_terms_df,
                             mf_go_terms_df,
                             cc_go_terms_df)

## Filter by significance:
combined_go_terms_df_filtered <- subset(combined_go_terms_df, 
                                        weight_ks_adjusted < 0.05)

## Log transform adjusted p values:
combined_go_terms_df_filtered$log10 <- -log10(combined_go_terms_df_filtered$weight_ks_adjusted)

## Create a new name term for plotting:
combined_go_terms_df_filtered$new_go_terms <- combined_go_terms_df_filtered$Term

## Generate barchart:

combined_go_terms_df_filtered$status <- factor(combined_go_terms_df_filtered$status, 
      levels = c("Molecular Function", "Cellular Component", "Biological Processes"))

combined_go_terms_df_filtered$new_go_terms <- gsub("(?<=^|; )([a-z])", "\\U\\1", 
                                                 combined_go_terms_df_filtered$new_go_terms, perl = T)

combined_go_terms_df_filtered$final_go_terms <- paste(combined_go_terms_df_filtered$new_go_terms,
                                                      " (", combined_go_terms_df_filtered$GO.ID, ")", sep = "")

bar_plot <- ggbarplot(combined_go_terms_df_filtered, x = "final_go_terms", y = "log10",
                      position = position_dodge(0.1),
                      fill = "status",# change fill color by mpg_level
                      color = NULL, #"white"            # Set bar border colors to white
                      palette = "jco",# jco journal color palett. see ?ggpar
                      sort.val = "desc",          # Sort the value in descending order
                      sort.by.groups = TRUE,     # Don't sort inside each group
                      ylab = "-log10(p)",
                      xlab = "Gene ontology (GO) term",
                      legend.title = "GO categories",
                      lab.col = "black",
                      lab.size = 4,
                      lab.vjust = 0.5,
                      lab.hjust = 1,
                      legend = "right",
                      rotate = TRUE,
                      ggtheme = theme_minimal()
)

## Make the axes labels bigger:
bar_plot <- bar_plot + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text=element_text(size=20),
        axis.title.x = element_text(size=20,face="bold"),
        axis.title.y = element_text(size=20,face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin=unit(c(0.5, 15, 0.5 , 0.5), "cm"))

## Update colours for plotting:
pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/5_GO_graphs/Go_terms_shortnames.pdf",
    width = 25, height = 12, onefile = TRUE)

bar_plot + 
  scale_fill_manual(values = c("dodgerblue2", 
                               "seagreen3", 
                               "brown2")) +
  geom_hline(yintercept = 1.301, linetype="dashed", colour="black") 

dev.off()

# ---------------------------------------> VIOLIN PLOTS:

# DATA UPLOAD:

# Normalized reads:

normalized_reads <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/normalized_reads_NOOut.tsv",
  header = TRUE,
  sep = "\t")

# Immune genes:

immune_genes <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/immune_genes.txt",
  header = TRUE,
  sep = "\t")

# DATA MANIPULATION:

immune_reads <- normalized_reads[rownames(normalized_reads) %in% immune_genes$Genes, ]
immune_reads_log <- log(immune_reads+1)

# VIOLIN PLOTS: 

immune_reads_log <- as.data.frame(t(immune_reads_log))
groups <- c(rep("Placebo", 18), rep("Oils", 18))

immune_reads <- as.data.frame(t(immune_reads))
immune_reads$groups <- groups


violin_table <- data.frame(groups = rep(groups, 15))
genes <- NULL
counts <- NULL

for (i in 1:ncol(immune_reads_log)) {
  genes[i] <- list(rep(colnames(immune_reads_log)[i], 36))
  counts[i] <- list(immune_reads_log[, i])
}

violin_table$genes <- unlist(genes)
violin_table$counts <- unlist(counts)
  
violin_table$genes <- factor(violin_table$genes, levels = colnames(immune_reads_log))

violin_table <- violin_table[violin_table$genes %in% c("PARK7", "S100A12", "GPX1", "CD74",
                                                       "IFIT3", "IFITM1", "IFITM3", "CD3D", "RPS19", "RAC2"), ]

# VIOLIN PLOTS:

ggplot(violin_table, aes(y = counts, x = groups)) +
  geom_violin(aes(fill = groups), trim = FALSE) +
  facet_wrap("genes", ncol = 5)

for (j in 1:(ncol(immune_reads)-1)) {
  pdf(file = paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/5_GO_graphs/", colnames(immune_reads)[j],
                   ".pdf", sep = ""), width = 25, height = 12)
  print(ggplot(immune_reads, aes(y = immune_reads[,j], x = groups)) +
  geom_violin(aes(fill = groups), trim = FALSE) +
  scale_y_continuous(name = colnames(immune_reads)[j]))
  dev.off()
}




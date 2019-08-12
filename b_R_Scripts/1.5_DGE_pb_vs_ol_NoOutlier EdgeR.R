# ---------------------------- RNA-Seq of Placebo vs Oil project with batch effects -----------------------------------#

# DGE analysis of the GC-JD-7904 project. 

# This was bulk RNA sequencing from 22 blood samples (paired), so is not single cell data. We used the QIASeq UPX 3' 
# Transcriptome kit. The libraries were sequenced on the NextSeq Mid-output kit with R1 = 151bp and R2 = 27bp. 
# The way the library is constructed means that we only get a single fastq for the whole library, rather than one 
# per sample and we have used the Qiagen analysis tool to demultiplex the data.

# ---> LIBRARIES: 

library(edgeR)
library(ggplot2)
library(gplots) 
library(ggrepel)

# ---> LOADING THE DATA:

# Open table with raw read counts:

read_counts <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/1.5_DGE_pb_vs_ol_NoOutlier/read_counts_Esteban.txt",
  header = TRUE,
  sep = "\t")

# Some genes are repeat (perhaps because a name problem during the mapping or some dupplication), but it's only
# a couple of them so it's sage to get rid off the duplicates. 

read_counts <- read_counts[!duplicated(read_counts$X), ] # Eliminate duplicates
row.names(read_counts) <- read_counts$X
read_counts$X <- NULL

# Open the txt file with the classification table (groups and batchs): 
# Since the RNA sequencing was made in two batches is important to have this consideration present at the moment
# of doing the DGE. 

classification <- read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/1.5_DGE_pb_vs_ol_NoOutlier/sample_info.txt",
  header = TRUE,
  sep = "\t")

# Based on the PCA and oPLS-DA analysis (see graphs in the same folder of this script) we identify a couple of 
# samples that, because of RNA sequencing quality or because they were natural outliers, werent located in the 
# right group (Placebo vs 4.5 supplement oil) and was affecting the results. The classiffication list has a 
# quality column with this information. 

classification <- classification[classification$quality == "ok", ]

# ---> DATA MANIPULATION:

# Keep only good samples: 

read_counts <- read_counts[ ,as.vector(classification$sample)]

# Create factors with groups, batchs and subjects:
# We need the group, batchs and subjects as factors for the design after in the DGE. 

groups <- factor(classification$group) 
names(groups) <- names(read_counts) # Keep the sample names in the group factor.
groups <- relevel(groups, ref = "Placebo") # Make sure that Placebo (control) becomes the reference factor.
                                           # This helps for not making mistakes at the moment of doing the 
                                           # comparison that you want to make.

batchs <- factor(classification$batch) # Batch as a factor.
subjects <- factor(classification$subject) # Subjects as a factor for the paired analysis. 

# Create the DDGEList object:
# DGEList takes raw counts. Not normalization step is requiered at this point. 

read_data <- DGEList(counts = read_counts, group = groups)

# Filtering and normalization: 

# Gene to be expressed at a reasonable level in a sample if it has two counts per each million mapped reads in that 
# sample. Gene should be expressed in at least 18 (half of the samples since is a paired analysis) to be in at least
# one of the conditions. 

keep <- rowSums(cpm(read_data)>2) >= 18 # Save only the genes that fullfil the condition. 
filter_read_data <- read_data[keep, , keep.lib.sizes = FALSE] # New DGEList object with only the filtered genes.

# Normalization: 

# Counts per million: 
# This normalization of counts per million has to be made before the other normalization, because it will be like
# normalizing over something that has been already normalized. This is visualization: scatter and violin plots. 

normalized_reads <- cpm(read_data)
normalized_filter_reads <- cpm(filter_read_data)

# TMN normalizacion: 
# Trimmed mean of Mvalues (TMM) between each pair of samples. We call the product of the original library size and 
# the scaling factor the effective library size. The effective library size replaces the original library size in all 
# downsteam analyses. TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes are 
# believed not differentially expressed between any pair of the samples. 

read_data <- calcNormFactors(read_data) 
filter_read_data <- calcNormFactors(filter_read_data)

# Save the CPM from the filtered data. 

write.table(normalized_filter_reads, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/5_GO_graphs/normalized_reads_NOOut.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# ---> DIFFERENTIAL GENE ANALYSIS (EDGER): 

# The design define how the response variable (DGE for each gene) is explained by the explained variables. In this
# case we are considering the group, the batch and the subject (because is paired sample experiment). 

design <- model.matrix(~subjects+batchs+groups) # Making the comparison by subjects. 
rownames(design) <- colnames(read_data)

# Estimate dispersion and have in consideration possible outliers. 

read_data_d <- estimateDisp(read_data, design, robust = TRUE) # Estimate the dispersion. 
filter_read_data_d <- estimateDisp(filter_read_data, design, robust = TRUE) # Robustified against potential outlier
                                                                            # genes.
# Here we run the statistcal analysys that check the DGE:
# We have to options: Likelihood ratio test and quasi-likelihood F test. The consensus said that the second one
# is best. The first one is recommended for SINGLE CELL RNA SEQ and NOT REPLICATES.

fit_data <- glmQLFit(read_data_d, design, robust = TRUE)
fit_filter_data <- glmQLFit(filter_read_data_d, design, robust = TRUE)

qlf_data <- glmQLFTest(fit_data)
qlf_data_filter <- glmQLFTest(fit_filter_data)

# Generates a table with the results of the DGE analysis. Then you can correct the p values using the desire
# methodology (BH my personal favorite). 

table_data <- as.data.frame(topTags(qlf_data, n = Inf, adjust.method = "BH", sort.by = "PValue"))
table_data_filter <- as.data.frame(topTags(qlf_data_filter, n = Inf, adjust.method = "BH", sort.by = "PValue"))

# Save the tables: 
# With the filtered data the results were more fructifery.

write.table(table_data, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/edge_NoOutlier_paired.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

write.table(table_data_filter, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/edge_filter_NoOutlier_paired.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)



edge_sig_all <- table_data_filter[(table_data_filter$PValue < 0.05), ] # A table with only the statistical significance genes. 
edge_sig_all$gene <- rownames(edge_sig_all)
edge_sig <- head(edge_sig_all, n = 10)
neg_edge_sig <- edge_sig_all[(edge_sig_all$logFC < 0) & (edge_sig_all$FDR < 0.05), ]
edge_sig_pv <- table_data_filter[(table_data_filter$FDR < 0.05), ]

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/volcano_edge_filter_NoOutlier.pdf",
    width = 25, height = 12, onefile = TRUE)

ggplot(table_data_filter, aes(logFC,-log10(PValue))) +  #LogFC as x and -log10(P.value) as y.  
  geom_point(colour = "gray16") +
  scale_y_continuous(name = "-Log10(p value)") +
  scale_x_continuous(name = "Log2(Fold Change)") +
  geom_text_repel(data = edge_sig, position = "identity", aes(label = gene), size = 8) +
  geom_text_repel(data = neg_edge_sig, position = "identity", aes(label = gene), size = 8) +
  {if (nrow(edge_sig_pv) > 0) geom_point(data = edge_sig_pv[edge_sig_pv$logFC > 0, ], 
                                         aes(logFC,-log10(PValue)), color = "red3") } +
  {if (nrow(edge_sig_pv) > 0) geom_point(data = edge_sig_pv[edge_sig_pv$logFC < 0, ], 
                                         aes(logFC,-log10(PValue)), color = "dodgerblue2") } + # Color only significant genes. 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 1) + 
  geom_text(aes(1.70, -log10(0.05), label = "p value = 0.05", vjust = 1.2), color = "gray0", size = 7.5) + 
  theme(legend.position = "none",
        axis.title = element_text(size = 25),
        axis.text.x  =  element_text(size = 20, hjust = 1, colour = "black"),
        axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1, 3, 1, 3), "cm"))  

dev.off()

# ---> DIFFERENTIAL GENE ANALYSIS (LIMMA): 

design_limma <- model.matrix(~subjects+batchs+groups)
# design_limma <- model.matrix(~batchs+groups)

limma_data <- voom(read_data, design_limma, normalize="quantile")
limma_data_filter <- voom(filter_read_data, design_limma, normalize="quantile")

# cor_data_filter <- duplicateCorrelation(limma_data_filter, design_limma) # Reduce correlation between different factors

fit_limma <- lmFit(limma_data, design_limma)  # Organized the values according to the design.
fit_limma_filter <- lmFit(limma_data_filter, design_limma)
# fit_limma_filter <- lmFit(limma_data_filter, design_limma, cor_data_filter$consensus.correlation) # No difference at all

fit2_limma <- eBayes(fit_limma, robust = TRUE)
fit2_limma_filter <- eBayes(fit_limma_filter, robust = TRUE)

limma_table <- topTable(fit2_limma, sort.by = "p",number = Inf, adjust.method="BH", coef = "groupsOils")
limma_table_filter <- topTable(fit2_limma_filter, sort.by = "p",number = Inf, adjust.method="BH", coef = "groupsOils")

write.table(limma_table, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/limma_paired_NoOutlier.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

write.table(limma_table_filter, 
            file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/limma_filter_paired_NoOutlier.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

edge_sig_l <- limma_table_filter[(limma_table_filter$P.Value < 0.05), ] # A table with only the statistical significance genes. 
edge_sig_l$gene <- rownames(edge_sig_l)
edge_sig_l <- head(edge_sig_l, n = 6)
edge_sig_pv <- limma_table_filter[(limma_table_filter$adj.P.Val < 0.05), ]

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/volcano_limma_filter_paired_NoOutlier.pdf",
    width = 25, height = 12, onefile = TRUE)

ggplot(limma_table_filter, aes(logFC,-log10(P.Value))) +  #LogFC as x and -log10(P.value) as y.  
  geom_point() +
  scale_y_continuous(name = "-Log10(P value)") +
  scale_x_continuous(name = "Log2(Fold Change)") +
  geom_text_repel(data = edge_sig_l, position = "identity", aes(label = gene), size = 5) +
  {if (nrow(edge_sig_pv) > 0) geom_point(data = edge_sig_pv, 
                                         aes(logFC,-log10(P.Value), 
                                             colour = "red")) } + # Color only significant genes. 
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed", size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 1),
        axis.text.y  = element_text(size = 20, hjust = 1))  

dev.off()

# ---> HEATMAP: 

# This analysis will be only concentrated in the genes that were differentiated expressed. I will be using 
# different approachs to see which one is the best: heatmap with LogFC, Normalized reads and raw reads. 

# Table preparation: 

# Vector with DE genes:
sign_genes <- edge_sig_all$gene
fdr_sign_genes <- rownames(edge_sig_pv)

# Creates matrix with raw and normalized reads of the GE genes (no corrected significant genes): 
matrix_reads <- read_counts[rownames(read_counts) %in% sign_genes, ]
matrix_reads <- as.matrix(matrix_reads)

matrix_norm_reads <- normalized_filter_reads[rownames(normalized_filter_reads) %in% sign_genes, ]

# Creates matrix with raw and normalized reads of the GE genes (corrected significant genes): 
matrix_fdr_reads <- read_counts[rownames(read_counts) %in% fdr_sign_genes, ]

matrix_fdr_norm_reads <- normalized_filter_reads[rownames(normalized_filter_reads) %in% fdr_sign_genes, ]

# Define a specific color for each of the populations.
placebo_index <- which((groups == "Placebo") == TRUE)
oils_index <-which((groups =="Oils") == TRUE)

colors <- NULL
colors[placebo_index] <- "blue"
colors[oils_index] <- "red"

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/norm_heatmap.pdf",
    width = 25, height = 12, onefile = TRUE)
par(xpd = TRUE)
heatmap.2(matrix_norm_reads,
          hclustfun = function(matrix_norm_reads) hclust(dist(t(matrix_norm_reads),method = "euclidean"),method = "ward.D2"),
          #breaks = breaks, 
          col = bluered(75),
          scale = "row",
          key.par = list(cex = 1.5),
          dendrogram = c("both"),
          colCol = colors, # Colors from columns
          density.info = "none",
          trace = "none",
          cexCol = 1.5,
          labRow = FALSE,
          margins = c(9, 9))
legend(x = 0.85, y =  0.99, inset=.05, title="Treatment:",  bty = "n",
       c("Placebo", "Oils"), fill=c("blue", "red"),
       text.font = 1, cex = 1.3)
rect(0.85, 0.89, 0.93, 0.99)
dev.off()

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/norm_heatmapsig.pdf",
    width = 25, height = 12, onefile = TRUE)
par(xpd = TRUE)
heatmap.2(matrix_fdr_norm_reads,
          hclustfun = function(matrix_fdr_norm_reads) hclust(dist(t(matrix_fdr_norm_reads),method = "euclidean"),method = "ward.D2"),
          scale = "row",
          #breaks = breaks, 
          col = bluered(75),
          key.par = list(cex = 1.5),
          dendrogram = c("both"),
          colCol = colors, # Colors from columns
          density.info = "none",
          trace = "none",
          cexCol = 1.5,
          cexRow = 0.7,
          margins = c(9, 9))
legend(x = 0.85, y =  0.99, inset=.05, title="Treatment:",  bty = "n",
       c("Placebo", "Oils"), fill=c("blue", "red"),
       text.font = 1, cex = 1.3)
rect(0.85, 0.89, 0.93, 0.99)
dev.off()

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/heatmap.pdf",
    width = 25, height = 12, onefile = TRUE)
par(xpd = TRUE)
heatmap.2(matrix_reads,
          hclustfun = function(matrix_reads) hclust(dist(t(matrix_reads),method = "euclidean"),method = "ward.D2"),
          #breaks = breaks, 
          col = bluered(75),
          scale = "row",
          key.par = list(cex = 1.5),
          dendrogram = c("both"),
          colCol = colors, # Colors from columns
          density.info = "histogram",
          trace = "none",
          cexCol = 1.5,
          labRow = FALSE,
          margins = c(9, 9))
legend(x = 0.85, y =  0.99, inset=.05, title="Treatment:",  bty = "n",
       c("Placebo", "Oils"), fill=c("blue", "red"),
       text.font = 1, cex = 1.3)
rect(0.85, 0.89, 0.93, 0.99)
dev.off()

pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/1.5_DGE_pb_vs_ol_NoOutlier/heatmapsig.pdf",
    width = 25, height = 12, onefile = TRUE)
par(xpd = TRUE)
heatmap.2(matrix_fdr_reads,
          hclustfun = function(matrix_fdr_reads) hclust(dist(t(matrix_fdr_reads),method = "euclidean"),method = "ward.D2"),
          scale = "row",
          #breaks = breaks, 
          col = bluered(75),
          key.par = list(cex = 1.5),
          dendrogram = c("both"),
          colCol = colors, # Colors from columns
          density.info =  "histogram",
          trace = "none",
          cexCol = 1.5,
          cexRow = 0.7,
          margins = c(9, 9))
legend(x = 0.85, y =  0.99, inset=.05, title="Treatment:",  bty = "n",
       c("Placebo", "Oils"), fill=c("blue", "red"),
       text.font = 1, cex = 1.3)
rect(0.85, 0.89, 0.93, 0.99)
dev.off()

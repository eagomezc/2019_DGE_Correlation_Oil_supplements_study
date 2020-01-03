#-------------------------------------- Spearman Correlation Analysis -------------------------------------------------#

# The Spearman correlation is the non-parametric version of the Pearson correlation. It measures the strengh and direction
# of association between two variables based on a ranking system. 

# This script is made to analysis the correlation between lipid mediators concentration againts changes in the expression 
# of peripheral blood cell activation markers 24h post fish oil supplementation.

# ---> LIBRARIES:

library(corrplot)
library(ggplot2)

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/hhy270/Documents/GitHub/2019_DGE_Correlation_Oil_supplements_study/a_Data/3_Correlation_analysis_(Corrplot)/"
output <- "C:/Users/hhy270/Documents/GitHub/2019_DGE_Correlation_Oil_supplements_study/c_Expected_Outputs/3_Correlation_analysis_(Corrplot)/"

# !!!! IMPORTANT: For this script to work the CORRELATION TABLE FILE has to be called: 3_Correlation_analysis_data.txt

#---> DATA LOAD: 

# Open the correlation table file.

# The table consists in a series of columns as next:
# First column "dataname" contains the name of the different samples (each row will contain the data of a specific sample)
# Second column "condition" contains inf about treatment that the patient (sample) received.
# From the third column we have the different features that are going to be compare (e.g: lipid mediators concentrations,
# cell counts, clinical scores, etc). 
# See a_Data/3_Correlation_analysis_(Corrplot)/3_Correlation_analysis_data.txt

data <-  read.table(
  file = paste(input, "3_Correlation_analysis_data.txt", sep = ""),
  header = TRUE,
  sep = "\t")

#---> DATA MANIPULATION:

# Separated the Lipid mediators concentration for the cell culture info:

lipids_mediators <- data[, c(1:27)]

# Get the cell culture info:

percentages <- data[, c(1, 2, 28: ncol(data))]

#--> CORRELATION ANALYSIS (SPEARMAN): 

for (cond in data[!duplicated(data$condition), ]$condition) {
  
  # Lipids mediators table:
  
  lm <- lipids_mediators[lipids_mediators$condition == cond, ] # Only takes one condition per time
  rownames(lm) <- lm$dataname # Define row names
  lm[, c(1,2)] <- NULL # Delete unnecesary colums
  
  # Cell culture table: 
  
  cells <- percentages[percentages$condition == cond, ]
  rownames(cells) <- cells$dataname
  cells[, c(1,2)] <- NULL
  
  # Correlation (Spearman):
  # It's going to take the columns from the first table and compare them with all the columns of the second table. 
  cor_final <- cor(lm, cells, method = "spearman")
  
  # P values table: 
  
  p_value_vector <- NULL # Creates empty vector
  p_value_data_frame <- data.frame(lm = colnames(lm)) # Creates p value data frame
  p_adj_value_data_frame <- data.frame(lm =  colnames(lm)) # Creates adj p value data frame
  
  for (i in 1:ncol(cells)) {
    x <- cells[, i]
    for (j in 1:ncol(lm)) {
      y <- lm[, j]
      correlation <- cor.test(x, y, method = "spearman", exact = TRUE)
      p_value_vector[j] <- correlation$p.value
    }
    p_adjust_vector <- p.adjust(p_value_vector, method = "BH")
    p_value_data_frame[, colnames(cells[i])] <- p_value_vector
    p_adj_value_data_frame[, colnames(cells[i])] <- p_adjust_vector
  }
  
  # Update table with the p values associated to every correlation. 
  
  rownames(p_value_data_frame) <- p_value_data_frame$lm
  p_value_data_frame$lm <- NULL
  p_value_data_frame <- as.matrix(p_value_data_frame)
  
  # Update table with the adj p values associated to every correlation. 
  
  rownames(p_adj_value_data_frame) <- p_adj_value_data_frame$lm
  p_adj_value_data_frame$lm <- NULL
  p_adj_value_data_frame <- as.matrix(p_adj_value_data_frame)
  
  #---> OUTPUT: 
  
  write.table(cor_final,
              file = paste(output, "3_", cond, "_correlation.tsv", sep = ""),
              row.names = TRUE,
              sep  = "\t")
  
  write.table(p_adj_value_data_frame,
              file = paste(output, "3_", cond, "_adj_p_values.tsv", sep = ""),
              row.names = TRUE,
              sep  = "\t")
  
  pdf(file = paste(output, "3_", cond, "_correlation_matrix.pdf", sep = ""), 
      width = 25, height = 12, onefile = TRUE)
  
  # Creates the correlation plot showing only significant results.
  
  adjust_p_value_analysis <- corrplot(cor_final, method = "circle", p.mat = p_adj_value_data_frame, sig.level = 0.1, 
                                      insig = "blank", tl.col = "black", tl.srt = 45, cl.pos = "r")
  
  print(adjust_p_value_analysis)
  
  dev.off()

}


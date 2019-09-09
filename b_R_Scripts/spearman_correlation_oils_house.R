#-------------------------------------- Spearman Correlation Analysis -------------------------------------------------#

# This script is to make Spearman correlation analysis between an specific column (lipid mediator) of the table againts 
# a set of columns with bacterial information.

# ---> LIBRARIES:

library("corrplot")
library("psych")
library("ggplot2")
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  tot <- sum(!is.na(var_name))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  message("Outliers identified: ", na2 - na1, " from ", tot, " observations")
  message("Proportion (%) of outliers: ", (na2 - na1) / tot*100)
  message("Mean of the outliers: ", mo)
  m2 <- mean(var_name, na.rm = T)
  message("Mean without removing outliers: ", m1)
  message("Mean if we remove outliers: ", m2)
  if((na2 - na1) > 0){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    message("Outliers successfully removed", "\n")
    return(invisible(dt))
  } else{
    message("Nothing changed", "\n")
    return(invisible(var_name))
  }
}

# ---> DATA UPLOAD:

datas <-  read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/7_correlation_analysis/correlation_data.txt",
  header = TRUE,
  sep = "\t")
datas <- datas[, -c(27:30, 37, 38)] # Delete some conditions columns because of Jesmond doing! 

#outlierKD(datas, RvD1)
#outlierKD(datas, AT.RvD1)

#datas <- datas[complete.cases(datas), ]

data_conditions <-  read.table(
  file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/input/7_correlation_analysis/condition_info.txt",
  header = TRUE,
  sep = "\t")

# ---> DATA MANIPULATION:

concentrations <- c(45:66)
cond <- "3g"

# Changes data names with info of the conditions:

data <- datas

data$dataname <- paste(data$dataname, data_conditions$condition, sep = "_")

# Table with lm:
lm <- data[, c(1:26)]
rownames(lm) <- lm$dataname
lm$dataname <- NULL
lm <- lm[concentrations, ]

# Table with difference percentages:
percentages <- data[, c(27: ncol(data))]
percentages <- percentages[concentrations, ]
rownames(percentages) <- rownames(lm)


# ---> CORRELATION ANALYSIS: 

# Correlation (Spearman):
# It's going to take the columns from the first table and compare them with the columns of the second table. 
cor_final <- cor(lm, percentages, method = "spearman")


p_value_vector <- NULL
p_value_data_frame <- data.frame(del = c(1:25))
p_adj_value_data_frame <- data.frame(del = c(1:25))

for (i in 1:ncol(percentages)) {
  x <- percentages[, i]
  for (j in 1:ncol(lm)) {
    y <- lm[, j]
    correlation <- cor.test(x, y, method = "spearman", exact = TRUE)
    p_value_vector[j] <- correlation$p.value
  }
  p_adjust_vector <- p.adjust(p_value_vector, method = "BH")
  p_value_data_frame[, colnames(percentages[i])] <- p_value_vector
  p_adj_value_data_frame[, colnames(percentages[i])] <- p_adjust_vector
}

rownames(p_value_data_frame) <- colnames(lm)
p_value_data_frame$del <- NULL
p_value_data_frame <- as.matrix(p_value_data_frame)

rownames(p_adj_value_data_frame) <- colnames(lm)
p_adj_value_data_frame$del <- NULL
p_adj_value_data_frame <- as.matrix(p_adj_value_data_frame)

# P values:

no_pvalue_analysis <- corrplot(cor_final, method = "circle", tl.col = "black", tl.srt = 45, title = "no p value")
p_value_analysis <- corrplot(cor_final, method = "circle", p.mat = p_value_data_frame, title = "me no adjust",
                             sig.level = 0.1, insig = "blank", tl.col = "black", tl.srt = 45)

# Correlation graph use!!!!!: 
#pdf(file = "C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/7_correlation_analysis/correlation_matrix.pdf",
#    width = 25, height = 12, onefile = TRUE)

#col3 <- colorRampPalette(c("red", "white", "blue")) 
adjust_p_value_analysis <- corrplot(cor_final, method = "circle", p.mat = p_adj_value_data_frame, sig.level = 0.1, 
                                    insig = "blank",col = col3, tl.col = "black", tl.srt = 45, tl.pos = "n",
                                    cl.pos = "b") # CHANGE THIS ONE

print(adjust_p_value_analysis)
#dev.off()

write.table(cor_final, 
            file = paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/7_correlation_analysis/By groups/",
                          cond, "_correlation_data_out_cd41&pa.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

write.table(p_adj_value_data_frame, 
            file = paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/7_correlation_analysis/By groups/",
                          cond, "_adjvalues_cd41&pa.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# ---> CORRELATION PLOTS: 

# Based on the correlation graph (adjust p value) call the interest columns from both tables (lm and percentages)
# to make the correlation plots, having in consideration the different concentrations: 

lm_significants <- c(3, 6, 8, 8, 18, 18, 18, 19)
percentajes_significants <- c(14, 14, 13, 14, 14, 15, 18, 13)

datas$condition <- factor(data_conditions$condition, levels = c("3g")) 

datas <- datas[c(45:66), ]

for (i in 1:length(lm_significants)) {

df_plot <- data.frame(lm = datas[, lm_significants[i]+1],
                      percentages = percentages[, percentajes_significants[i]],
                      condition = datas$condition)

outlierKD(df_plot, lm)
outlierKD(df_plot, percentages)

df_plot <- df_plot[complete.cases(df_plot), ]

model1=lm(df_plot$lm~df_plot$percentages)

p <- qplot(df_plot$percentages,
      df_plot$lm,
      geom = "point",
      col = df_plot$condition) + 
  geom_line(data = fortify(model1), aes(x = df_plot$percentages, y = .fitted), color = "blue")

pdf(file = paste("n_out/",
                 colnames(datas)[lm_significants[i]+1], "_vs_",
                 colnames(percentages)[percentajes_significants[i]], ".pdf", sep = ""),
    width = 25, height = 12, onefile = TRUE)
print(p)
dev.off()

}
  

for (i in 1:length(lm_significants)) {
  df_plot <- data.frame(lm = rank(datas[, lm_significants[i]+1], ties.method = "min"),
                        percentages = rank(percentages[, percentajes_significants[i]], ties.method = "min"),
                        condition = datas$condition)
  
  model1=lm(df_plot$lm~df_plot$percentages)
  
  p <- ggplot(data = df_plot) +  
    geom_point(aes(x=percentages, y = lm, color = condition)) +
    geom_line(data = fortify(model1), aes(x = df_plot$percentages, y = .fitted), color = "blue")
    #geom_smooth(aes(x=percentages, y = lm), method = "auto", se = FALSE)  
  
  pdf(file = paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/7_correlation_analysis/rank line/rank_line_",
                   colnames(datas)[lm_significants[i]+1], "_vs_",
                   colnames(percentages)[percentajes_significants[i]], ".pdf", sep = ""),
      width = 25, height = 12, onefile = TRUE)
  print(p)
  dev.off()
  
}

for (i in 1:length(lm_significants)) {
  df_plot <- data.frame(lm = log2(datas[, lm_significants[i]+1] + 94),
                        percentages = log2(percentages[, percentajes_significants[i]] + 94),
                        condition = datas$condition)
  
  outlierKD(df_plot, lm)
  outlierKD(df_plot, percentages)
  
  df_plot <- df_plot[complete.cases(df_plot), ]
  
  model1=lm(df_plot$lm~df_plot$percentages)
  
  p <- ggplot(data = df_plot) +  
    geom_point(aes(x=percentages, y = lm, color = condition)) +
    geom_line(data = fortify(model1), aes(x = df_plot$percentages, y = .fitted), color = "blue")
    #geom_smooth(aes(x=percentages, y = lm), method = "auto", se = FALSE)  
  
  pdf(file = paste("C:/Users/hhy270/Dropbox/GC-JD-7904 project/output/7_correlation_analysis/groups/log_NOout_3g/",
                   colnames(datas)[lm_significants[i]+1], "_vs_",
                   colnames(percentages)[percentajes_significants[i]], ".pdf", sep = ""),
      width = 25, height = 12, onefile = TRUE)
  print(p)
  dev.off()
  
}

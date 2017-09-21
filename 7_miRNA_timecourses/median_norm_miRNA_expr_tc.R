# The MIT License
# 
# Copyright (c) 2017 Piero Dalle Pezze
# 
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to 
# deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom 
# the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




# Plot time courses of all miRNA. 
# - X is Time
# - Y is median of standardised miRNA expression
# - colour is PAM clustering



library(reshape2)
library(ggplot2)

library(ggfortify) # autoplot
library(cluster) # pam

library(grid)
library(gridExtra)

source('../utilities/plots.R')



# Run PAM clustering (Partition around medoids)
run_pam <- function(df, k, name, scale., centre) {

  pamx <- pam(df, k, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
  # Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
  write.csv(pamx$clustering, file=paste0(name,k,"_labels.csv"), quote=FALSE)

  return (pamx)
}


run_pca <- function(df, name, scale., centre) {
  
  # COMPUTE PCA
  pca <- prcomp(t(df), scale.=scale., center=centre)
  # Write the calculated PCA `rotation` (PCA load). This is the table miRNAs (rows) vs PCAs (cols).
  write.csv(pca$rotation, file=paste0(name, '_pca_rotation.csv'), quote=FALSE)
  
  # Plot the variance of PCA components
  png(paste0(name,"_pca_variances",".png"), width=2000, height=1500, res=300)
  plot(pca, type = "l", main='Variance of PCA components')
  dev.off()
  
  return(pca)
}


filt_pca <- function(pca, component=1, threshold=0.05) {
  pca.rot <- pca$rotation
  # filter the table
  pca.rot.filt <- pca.rot[ pca.rot[, component] < threshold | pca.rot[, component] > threshold, ]
  
  # return the list of filtered ids
  return(rownames(pca.rot.filt))
}




################
# Load data sets
################

clusters <- 2

# select the files containing the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered_scaled"
filename.data.table <- "summarised_mirna_counts_after_mapping_filtered_data_table"
suffix <-".csv"

# load counts
counts.scaled <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=T,header=T,row.names=1)

# load counts data table (this contains meta data)
counts.data.table <- read.table(paste0(location,"/",filename.data.table,suffix), sep=",",fill=T,header=T)




#############################
# Prepare data frame to plot
#############################

# retrieve the experiment groups (without repeat number). We only need to drop the last two characters from 
# the name of each column. 
groups <- gsub('.{2}$', '', colnames(counts.scaled))

# remove duplicates
groups <- groups[!duplicated(groups)]


# calculate the median for each group of repeats.
counts.scaled.median <- matrix(0, ncol=length(groups), nrow = nrow(counts.scaled))

# extract the columns of replicates for each group and calculate the median by row. This new column will form a new data set with groups as column names.
for(i in seq(1:length(groups))) {
  # extract the columns of replicates
  counts.g <- counts.scaled[, grepl(groups[i], colnames(counts.scaled))]
  # calculate the median per row
  counts.g.median <- apply(counts.g, 1, median)
  # add this column to the median data frame
  counts.scaled.median[, i] <- counts.g.median
}

counts.scaled.median <- data.frame(counts.scaled.median)
rownames(counts.scaled.median) <- rownames(counts.scaled)
colnames(counts.scaled.median) <- groups




######################################################
# Separate counts.scaled.median into three data frames
######################################################

# Separate miWT, miA66, miA66_noEGF. We also add miRNA names as a new column (we need this for melting)
counts.scaled.miWT <- data.frame(miRNA=rownames(counts.scaled), counts.scaled.median[, grepl("miWT", colnames(counts.scaled.median)), drop=FALSE])
counts.scaled.miA66 <- data.frame(miRNA=rownames(counts.scaled), counts.scaled.median[, grepl("miA66", colnames(counts.scaled.median)), drop=FALSE])
counts.scaled.miA66.noEGF <- data.frame(miRNA=rownames(counts.scaled), counts.scaled.median[, grepl("noEGF", colnames(counts.scaled.median)), drop=FALSE])

# remove column `miA66_noEGF`` from miA66 data frame as this is treated separately.
counts.scaled.miA66 <- subset(counts.scaled.miA66, select = -miA66_noEGF)

# rename columns so that they only contain times (useful for melting later on)
colnames(counts.scaled.miWT) <- gsub('miWT_', '', colnames(counts.scaled.miWT))
colnames(counts.scaled.miA66) <- gsub('miA66_', '', colnames(counts.scaled.miA66))
# this is renamed manually
colnames(counts.scaled.miA66.noEGF) <- c('miRNA', '300')




###################################################
# Plot all miRNA time courses by WT, A66, and noEGF
###################################################

# plot
plot_expr_tc(counts.scaled.miWT, "miRNA_tc__miWT", "(WT)", line=TRUE)
plot_expr_tc(counts.scaled.miA66, "miRNA_tc__miA66", "(A66)", line=TRUE)
plot_expr_tc(counts.scaled.miA66.noEGF, "miRNA_tc__miA66noEGF", "(A66 noEGF)", line=FALSE)




#########################
# ADD PAM clustering info
#########################

# `scale` and `centre` must be FALSE because we already scaled and centred our data.
scale. <- FALSE
centre <- FALSE

# run PAM clustering for WT, A66, A66noEGF
pam.miWT <- run_pam(subset(counts.scaled.miWT, select = -c(miRNA), drop=FALSE), clusters, "median_counts_WT_samples__pam", scale., centre)
pam.miA66 <- run_pam(subset(counts.scaled.miA66, select = -c(miRNA), drop=FALSE), clusters, "median_counts_A66_samples__pam", scale., centre)
pam.miA66.noEGF <- run_pam(subset(counts.scaled.miA66.noEGF, select = -c(miRNA), drop=FALSE), clusters, "median_counts_A66noEGF_samples__pam", scale., centre)

# add PAM labels
counts.scaled.miWT.pam <- data.frame(counts.scaled.miWT, pam=factor(pam.miWT$clustering))
counts.scaled.miA66.pam <- data.frame(counts.scaled.miA66, pam=factor(pam.miA66$clustering))
counts.scaled.miA66.noEGF.pam <- data.frame(counts.scaled.miA66.noEGF, pam=factor(pam.miA66.noEGF$clustering))

# plot
plot_expr_tc_w_pam(counts.scaled.miWT.pam, "miRNA_tc__miWT_pam_clust", "(WT)", line=TRUE)
plot_expr_tc_w_pam(counts.scaled.miA66.pam, "miRNA_tc__miA66_pam_clust", "(A66)", line=TRUE)
plot_expr_tc_w_pam(counts.scaled.miA66.noEGF.pam, "miRNA_tc__miA66noEGF_pam_clust", "(A66 noEGF)", line=FALSE)




#####################################################################################################
# Median of the miRNAs belonging to the same cluster class. This for miWT, miA66, and miA66.noEGF
#####################################################################################################

# plot
plot_median_expr_tc_w_pam(subset(counts.scaled.miWT.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miWT_pam_clust_median', '(WT)', line=TRUE)
plot_median_expr_tc_w_pam(subset(counts.scaled.miA66.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miA66_pam_clust_median', '(A66)', line=TRUE)
plot_median_expr_tc_w_pam(subset(counts.scaled.miA66.noEGF.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miA66noEGF_pam_clust_median', '(A66 noEGF)', line=FALSE)



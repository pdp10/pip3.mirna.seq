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
# - Y is mean of standardised miRNA expression
# - colour is PAM clustering



library(reshape2)
library(ggplot2)

library(cluster) # pam

library(grid)
library(gridExtra)


source('../utilities/clustering.R')
source('../utilities/plots.R')





################
# Load data sets
################

suffix <-".csv"

# select the files containing the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered_rlog_scaled"
filename.mean <- "summarised_mirna_counts_after_mapping_filtered_rlog_mean_scaled"


# load counts
counts <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=T,header=T,row.names=1)

# load counts means
counts.mean <- read.table(paste0(location,"/",filename.mean,suffix), sep=",",fill=T,header=T,row.names=1)



# clustering
location.pam <- "../6_clustering"
filename.pam <- "pam_clustering_labels"
# load scaled counts
pam.labels <- read.table(paste0(location.pam,"/",filename.pam, suffix), sep=",",fill=T,header=T,row.names=1)

# extract the number of clusters
pam.classes <- unique(pam.labels[,1])


# combine the counts matrix with PAM clustering
counts.mean.wpam <- cbind(counts.mean, pam=pam.labels[,1])




###########
# Filtering
# THIS IS THE PART THAT IS DIFFERENT FROM 7_time_courses
###########

### Load significant miRNA calculated using DESeq2:strain (padj<=0.05, |lfc|>0.1) && PCA:PC1 (>0.05). These are 99 miRNA (see Venn diagram)
# select the files containing the data
location.venn.intersect <- "../6_clustering"
filename.venn.intersect <- "venn_diagram_intersect__filt_pca_PC1_VS_signif_deseq_strain"
# load the significant miRNA names 
miRNA.venn.intersect <- rownames(read.table(paste0(location.venn.intersect,"/",filename.venn.intersect,suffix), sep=",",fill=T,header=T, row.names=1))


# NOTE: WE DO NOT RUN PAM ON THIS NEW SUBSET, BUT RE-USE PAM LABELS AS COMPUTED PREVIOUSLY.
# filter the counts table with the Venn diagram intersection
counts.mean.wpam <- subset(counts.mean.wpam, rownames(counts.mean.wpam) %in% miRNA.venn.intersect)












###################################################
# Split the count matrix depending on cluster class
###################################################
# split_counts_matrix_by_pam.class() returns a list of lists
l <- split_counts_matrix_by_pam_class(df.wpam=counts.mean.wpam, pam.classes=pam.classes)
counts.clusters.wt.a66 <- l[[1]]
counts.clusters.a66.noEGF <- l[[2]]


##################
# Plot data frames
##################

for(pam.class in pam.classes) { 
  # melt the data.table
  df.wt.a66.melt <- melt(counts.clusters.wt.a66[[pam.class]], id=c('id', 'colour'))
  df.a66.noEGF.melt <- melt(counts.clusters.a66.noEGF[[pam.class]], id=c('id', 'colour'))
  
  plot_all_expr_tc_wcolour(df.line=df.wt.a66.melt, 
                           df.point=df.a66.noEGF.melt, 
                           filename=paste0('miRNA_tc__clust_', pam.class, '.png'),
                           title=paste0('miRNA expression - clust: ', pam.class))
}





####################################################
# Calculate the mean for each PAM class and strain
####################################################
# split_counts_matrix_by_pam.class() returns a list of lists
l <- mean_for_each_pam_class_and_strain(df.wt.a66=counts.clusters.wt.a66, df.a66.noEGF=counts.clusters.a66.noEGF, pam.classes=pam.classes)
counts.clusters.wt.a66.mean <- l[[1]]
counts.clusters.a66.noEGF.mean <- l[[2]]


#############################
# Plot the mean data frames
#############################

for(pam.class in pam.classes) {
  # melt the data.table
  df.wt.a66.melt <- melt(counts.clusters.wt.a66.mean[[pam.class]], id=c('colour'))
  df.a66.noEGF.melt <- melt(counts.clusters.a66.noEGF.mean[[pam.class]], id=c('colour'))
  
  plot_mean_expr_tc_wcolour(df.line=df.wt.a66.melt,
                            df.point=df.a66.noEGF.melt,
                            filename=paste0('miRNA_tc__mean_clust_', pam.class, '.png'),
                            title=paste0('miRNA expression - clust: ', pam.class))
}



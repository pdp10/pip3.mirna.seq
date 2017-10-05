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


source('../utilities/clustering.R')
source('../utilities/plots.R')





################
# Load data sets
################

suffix <-".csv"

# select the files containing the data
location <- "../data"
filename.scaled <- "summarised_mirna_counts_after_mapping_filtered_scaled"
filename.median.scaled <- "summarised_mirna_counts_after_mapping_filtered_median_scaled"
#filename.data.table <- "summarised_mirna_counts_after_mapping_filtered_data_table"

# load scaled counts
counts.scaled <- read.table(paste0(location,"/",filename.scaled,suffix), sep=",",fill=T,header=T,row.names=1)

# load median scaled counts
counts.median.scaled <- read.table(paste0(location,"/",filename.median.scaled,suffix), sep=",",fill=T,header=T,row.names=1)



# clustering
location.pam <- "../6_clustering"
filename.pam <- "pam_clustering_labels"
# load scaled counts
pam.labels <- read.table(paste0(location.pam,"/",filename.pam, suffix), sep=",",fill=T,header=T,row.names=1)

# extract the number of clusters
pam.classes <- unique(pam.labels[,1])


# combine the counts matrix with PAM clustering
counts.median.scaled.wpam <- cbind(counts.median.scaled, pam=pam.labels[,1])






###################################################
# Split the count matrix depending on cluster class
###################################################
# split_counts_matrix_by_pam.class() returns a list of lists
l <- split_counts_matrix_by_pam_class(df.wpam=counts.median.scaled.wpam, pam.classes=pam.classes)
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
# Calculate the median for each PAM class and strain
####################################################
# split_counts_matrix_by_pam.class() returns a list of lists
l <- median_for_each_pam_class_and_strain(df.wt.a66=counts.clusters.wt.a66, df.a66.noEGF=counts.clusters.a66.noEGF, pam.classes=pam.classes)
counts.clusters.wt.a66.median <- l[[1]]
counts.clusters.a66.noEGF.median <- l[[2]]


#############################
# Plot the median data frames
#############################

for(pam.class in pam.classes) {
  # melt the data.table
  df.wt.a66.melt <- melt(counts.clusters.wt.a66.median[[pam.class]], id=c('colour'))
  df.a66.noEGF.melt <- melt(counts.clusters.a66.noEGF.median[[pam.class]], id=c('colour'))

  plot_median_expr_tc_wcolour(df.line=df.wt.a66.melt,
                              df.point=df.a66.noEGF.melt,
                              filename=paste0('miRNA_tc__median_clust_', pam.class, '.png'),
                              title=paste0('miRNA expression - clust: ', pam.class))
}



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


# Calculate and plot the correlation matrix using pearson correlation coefficients


source('../utilities/plots.R')


################
# Load data sets
################

# select the files containing the data
location <- "../data"
filename.counts.log <- "summarised_mirna_counts_after_mapping_filtered_norm_log"
filename.counts.median.log <- "summarised_mirna_counts_after_mapping_filtered_median_log"
filename.counts.rlog <- "summarised_mirna_counts_after_mapping_filtered_norm_rlog"
suffix <-".csv"

# load counts (log)
counts.log <- read.table(paste0(location,"/", filename.counts.log, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)

# load log of median counts
counts.median.log <- read.table(paste0(location,"/", filename.counts.median.log, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)

# load counts (rlog)
counts.rlog <- read.table(paste0(location,"/", filename.counts.rlog, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


### Load significant miRNA calculated using DESeq2:strain (padj<=0.05, |lfc|>0.1) && PCA:PC2 (>0.05). These are 89 miRNA (see Venn diagram)
# select the files containing the data
location.venn.intersect <- "../6_clustering"
filename.venn.intersect <- "venn_diagram_intersect__filt_pca_pc2_VS_signif_deseq_strain"
# load the significant miRNA names 
miRNA.venn.intersect <- rownames(read.table(paste0(location.venn.intersect,"/",filename.venn.intersect,suffix), sep=",",fill=T,header=T, row.names=1))






###########
# Filtering
# THIS IS THE PART THAT IS DIFFERENT FROM 7_time_courses
###########

# filter the counts table with the Venn diagram intersection
counts.log <- subset(counts.log, rownames(counts.log) %in% miRNA.venn.intersect)
counts.median.log <- subset(counts.median.log, rownames(counts.median.log) %in% miRNA.venn.intersect)
counts.rlog <- subset(counts.rlog, rownames(counts.rlog) %in% miRNA.venn.intersect)







#########################################################
# Compute heatplot of the counts matrix and counts density
#########################################################
print('plot_counts_matrix_heatmap')

# heatplot for counts.log.scaled (heatplot scales for us)
# Compute the log first, then apply a z-trasform to N(0, 1) using scale().
# This works better than just using scale().
# Z transformation
plot_counts_matrix_heatplot(counts.log, filename="miRNA_counts_matrix_heatplot__log_scaled_by_row.png", scale='row')

# heatplot for counts.median (heatplot scales for us) 
plot_counts_matrix_heatplot(counts.median.log, filename="miRNA_counts_matrix_heatplot__median_log_scaled_by_row.png", scale='none')

# heatplot for rlog(counts) (heatplot scales for us)
plot_counts_matrix_heatplot(counts.rlog, filename="miRNA_counts_matrix_heatplot__rlog_scaled_by_row.png", scale='row')


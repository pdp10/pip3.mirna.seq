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
filename.counts <- "summarised_mirna_counts_after_mapping_filtered_rlog_scaled"
filename.counts.mean <- "summarised_mirna_counts_after_mapping_filtered_rlog_mean_scaled"
suffix <-".csv"

# load counts (log)
counts <- read.table(paste0(location,"/", filename.counts, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)

# load log of mean counts
counts.mean <- read.table(paste0(location,"/", filename.counts.mean, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)





###########
# Filtering
# THIS IS THE PART THAT IS DIFFERENT FROM 7_time_courses
###########

### Load significant miRNA calculated using DESeq2:strain (padj<=0.05, |lfc|>0.1) && PCA:PC1 (>0.05). These are 89 miRNA (see Venn diagram)
# select the files containing the data
location.venn.intersect <- "../6_clustering"
filename.venn.intersect <- "venn_diagram_intersect__filt_pca_PC1_VS_signif_deseq_strain"
# load the significant miRNA names 
miRNA.venn.intersect <- rownames(read.table(paste0(location.venn.intersect,"/",filename.venn.intersect,suffix), sep=",",fill=T,header=T, row.names=1))

# filter the counts table with the Venn diagram intersection
counts <- subset(counts, rownames(counts) %in% miRNA.venn.intersect)
counts.mean <- subset(counts.mean, rownames(counts.mean) %in% miRNA.venn.intersect)







#########################################################
# Compute heatplot of the counts matrix and counts density
#########################################################
print('plot_counts_matrix_heatmap')

plot_counts_matrix_heatplot(counts, filename="miRNA_counts_matrix_heatplot__rlog_scaled.png", scale='none')

# heatplot for samples (mean of repeats) (heatplot scales for us) 
plot_counts_matrix_heatplot(counts.mean, filename="miRNA_counts_matrix_heatplot__rlog_mean_scaled.png", scale='none')


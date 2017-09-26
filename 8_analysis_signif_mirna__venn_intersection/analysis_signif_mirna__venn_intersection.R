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



source('../utilities/plots.R')



################
# Load data sets
################

suffix <-".csv"

### Load data sets

# select the files containing the data
location.data <- "../data"
filename.scaled <- "summarised_mirna_counts_after_mapping_filtered_scaled"
filename.median.scaled <- "summarised_mirna_counts_after_mapping_filtered_median_scaled"

# load scaled counts
counts.scaled <- read.table(paste0(location.data,"/",filename.scaled,suffix), sep=",",fill=T,header=T,row.names=1)

# load median scaled counts
counts.median.scaled <- read.table(paste0(location.data,"/",filename.median.scaled,suffix), sep=",",fill=T,header=T,row.names=1)



### Load significant miRNA calculated using DESeq2:strain (padj<=0.05) && PCA:PC2 (>0.05). These are 89 miRNA (see Venn diagram)

# select the files containing the data
location.venn.intersect <- "../5_compare_deseq_signif_miRNA_vs_pca_loading"
filename.venn.intersect <- "summarised_mirna_counts_after_mapping_filtered_scaled__Venn_Diagram_intersect_miRNA"

# load the significant miRNA names 
miRNA.venn.intersect <- rownames(read.table(paste0(location.venn.intersect,"/",filename.venn.intersect,suffix), sep=",",fill=T,header=T, row.names=1))



###########################################
# Create counts matrix of significant miRNA
###########################################

# Filtered counts matrix of significant miRNA
counts.scaled.venn.intersect <- subset(counts.scaled, rownames(counts.scaled) %in% miRNA.venn.intersect)

# Filtered counts matrix of significant miRNA, with median of sample replicates 
counts.median.scaled.venn.intersect <- subset(counts.median.scaled, rownames(counts.median.scaled) %in% miRNA.venn.intersect)




#######################################
# Plot heatplot for these mini-matrices
#######################################

plot_counts_matrix_heatplot(counts.scaled.venn.intersect, filename="miRNA_counts_matrix_heatplot__venn_intersect.png", scale='none', method="ward.D2", labRow=FALSE)

plot_counts_matrix_heatplot(counts.median.scaled.venn.intersect, filename="miRNA_counts_matrix_heatplot__venn_intersect_median.png", scale='none', method="ward.D2", labRow=FALSE)





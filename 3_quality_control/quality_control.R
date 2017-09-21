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


library(reshape2)
library(ggplot2)

source('../utilities/plots.R')



################
# Load data sets
################

# select the files containing the data
location <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered"
filename.counts.data.table <- "summarised_mirna_counts_after_mapping_filtered_data_table"
suffix <-".csv"


# load counts
counts <- read.table(paste0(location,"/", filename.counts, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)

# load counts data table
counts.data.table <- read.table(paste0(location,"/", filename.counts.data.table, suffix), sep=",", fill=TRUE, header=TRUE)



################################################
# Compute and plot correlation matrix of samples
################################################
print('plot_samples_correlation_matrix')
p <- plot_corr_matrix(counts, "pearson")

# Add block separators (black lines)
separators <- data.frame(x=c(17.5,17.5,35.5,35.5), y=c(17.5,.5,.5,35.5), xend=c(37.5,17.5,35.5,37.5), yend=c(17.5,17.5,35.5,35.5))
p <- p + geom_segment(data=separators, aes(x, y, xend=xend, yend=yend), size=0.5, inherit.aes=F)

ggsave(paste0("miRNA_samples_pears_corr_matrix.png"), width=4, height=4, dpi=300)




###################################################
# Compute read density and cumulative distributions
###################################################
print('plot_read_density')
plot_read_density(counts.data.table)

print('plot_read_ecdf')
plot_read_ecdf(counts.data.table)



#############################################
# Extract miRNA with counts greater than 400k
#############################################
print('pca')
source('pca.R')
rm(list = ls())



###############################################
# Calculate coefficient of variance for Samples
###############################################
print('samples_coeff_var')
source('samples_coeff_var.R')
rm(list = ls())



#############################################
# Extract miRNA with counts greater than 400k
#############################################
print('extract_ids_gt_k')
source('extract_ids_gt_k.R')
extract_ids_gt_k <- function(counts, k=400000)
rm(list = ls())

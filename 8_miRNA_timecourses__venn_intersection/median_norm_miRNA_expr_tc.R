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

clusters <- 3

# select the files containing the data
location <- "../data"
filename.scaled <- "summarised_mirna_counts_after_mapping_filtered_scaled"
filename.median.scaled <- "summarised_mirna_counts_after_mapping_filtered_median_scaled"
#filename.data.table <- "summarised_mirna_counts_after_mapping_filtered_data_table"
suffix <-".csv"

# load scaled counts
counts.scaled <- read.table(paste0(location,"/",filename.scaled,suffix), sep=",",fill=T,header=T,row.names=1)

# load median scaled counts
counts.median.scaled <- read.table(paste0(location,"/",filename.median.scaled,suffix), sep=",",fill=T,header=T,row.names=1)


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
counts.scaled <- subset(counts.scaled, rownames(counts.scaled) %in% miRNA.venn.intersect)
counts.median.scaled <- subset(counts.median.scaled, rownames(counts.median.scaled) %in% miRNA.venn.intersect)







######################################################
# Separate counts.median.scaled into three data frames
######################################################

# Separate miWT, miA66, miA66_noEGF. We also add miRNA names as a new column (we need this for melting)
counts.scaled.miWT <- data.frame(miRNA=rownames(counts.scaled), counts.median.scaled[, grepl("miWT", colnames(counts.median.scaled)), drop=FALSE])
counts.scaled.miA66 <- data.frame(miRNA=rownames(counts.scaled), counts.median.scaled[, grepl("miA66", colnames(counts.median.scaled)), drop=FALSE])
counts.scaled.miA66.noEGF <- data.frame(miRNA=rownames(counts.scaled), counts.median.scaled[, grepl("noEGF", colnames(counts.median.scaled)), drop=FALSE])

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

# melt the data.table
df.melt <- melt(counts.scaled.miWT, id=c('miRNA'))
# plot
plot_expr_tc(df.melt, filename="miRNA_tc__miWT.png", line=TRUE, title='miRNA expression (WT)', xlab='time [m]', ylab='median standardised expression')

# melt the data.table
df.melt <- melt(counts.scaled.miA66, id=c('miRNA'))
# plot
plot_expr_tc(df.melt, filename="miRNA_tc__miA66.png", line=TRUE, title='miRNA expression (A66)', xlab='time [m]', ylab='median standardised expression')

# melt the data.table
df.melt <- melt(counts.scaled.miA66.noEGF, id=c('miRNA'))
# plot
plot_expr_tc(df.melt, filename="miRNA_tc__miA66noEGF.png", line=FALSE, title='miRNA expression (A66 noEGF)', xlab='time [m]', ylab='median standardised expression')




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
counts.scaled.miWT.pam <- data.frame(counts.scaled.miWT, colour=factor(pam.miWT$clustering), check.names=FALSE)
counts.scaled.miA66.pam <- data.frame(counts.scaled.miA66, colour=factor(pam.miA66$clustering), check.names=FALSE)
counts.scaled.miA66.noEGF.pam <- data.frame(counts.scaled.miA66.noEGF, colour=factor(pam.miA66.noEGF$clustering), check.names=FALSE)

# melt & plot
df.melt <- melt(counts.scaled.miWT.pam, id=c('miRNA', 'colour'))
plot_expr_tc_wcolour(df.melt, filename="miRNA_tc__miWT_pam_clust.png", line=TRUE, gradient=FALSE, title='miRNA expression (WT)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')
# melt & plot
df.melt <- melt(counts.scaled.miA66.pam, id=c('miRNA', 'colour'))
plot_expr_tc_wcolour(df.melt, filename="miRNA_tc__miA66_pam_clust.png", line=TRUE, gradient=FALSE, title='miRNA expression (A66)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')
# melt & plot
df.melt <- melt(counts.scaled.miA66.noEGF.pam, id=c('miRNA', 'colour'))
plot_expr_tc_wcolour(df.melt, filename="miRNA_tc__miA66noEGF_pam_clust.png", line=FALSE, gradient=FALSE, title='miRNA expression (A66 noEGF)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')




#####################################################################################################
# Median of the miRNAs belonging to the same cluster class. This for miWT, miA66, and miA66.noEGF
#####################################################################################################

# calculate the median
dt.median <- median_df_by_colour_id(subset(counts.scaled.miWT.pam, select = -c(miRNA), drop=FALSE), filename='miRNA_tc_miWT_pam_clust_median.csv')
# we restore this column to make the plot work. We group by miRNA. Here is pointless really, but this avoids the creation of a new plot function
dt.median$miRNA <- dt.median$colour
# melt the data.table
dt.median.melt <- melt(dt.median, id=c('miRNA', 'colour'))
# plot
plot_expr_tc_wcolour(dt.median.melt, filename='miRNA_tc_miWT_pam_clust_median.png', 
                     line=TRUE, gradient=FALSE, title='miRNA expression (WT)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')

dt.median <- median_df_by_colour_id(subset(counts.scaled.miA66.pam, select = -c(miRNA), drop=FALSE), filename='miRNA_tc_miA66_pam_clust_median.csv')
# we restore this column to make the plot work. We group by miRNA. Here is pointless really, but this avoids the creation of a new plot function
dt.median$miRNA <- dt.median$colour
# melt the data.table
dt.median.melt <- melt(dt.median, id=c('miRNA', 'colour'))
# plot
plot_expr_tc_wcolour(dt.median.melt, filename='miRNA_tc_miA66_pam_clust_median.png', 
                     line=TRUE, gradient=FALSE, title='miRNA expression (A66)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')

dt.median <- median_df_by_colour_id(subset(counts.scaled.miA66.noEGF.pam, select = -c(miRNA), drop=FALSE), filename='miRNA_tc_miA66noEGF_pam_clust_median.csv')
# we restore this column to make the plot work. We group by miRNA. Here is pointless really, but this avoids the creation of a new plot function
dt.median$miRNA <- dt.median$colour
# melt the data.table
dt.median.melt <- melt(dt.median, id=c('miRNA', 'colour'))
# plot
plot_expr_tc_wcolour(dt.median.melt, filename='miRNA_tc_miA66noEGF_pam_clust_median.png', 
                     line=FALSE, gradient=FALSE, title='miRNA expression (A66noEGF)', xlab='time [m]', ylab='median standardised expression', colorlab='pam')



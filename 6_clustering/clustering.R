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




### Calculate clustering


library(factoextra) # eigenvalues
library(ggplot2)
library(cluster) # pam

library(grid)
library(gridExtra)

source('../utilities/plots.R')




################
# Load data sets
################

clusters <- 3
scale. <- FALSE
suffix <-".csv"

# Counts Matrix
location.counts <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered_rlog_scaled"
counts <- read.table(paste0(location.counts, "/", filename.counts, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# PCA rotation matrix for counts
location.pca <- "../3_quality_control"
filename.pca <- "summarised_mirna_counts_after_mapping_filtered_rlog_scaled_samples_PCA_rotation"
df.pca <- read.table(paste0(location.pca, "/", filename.pca, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# DESeq contrast: strain (using WALD Test (~t-test))
# DESeq2 volcano plot: summarised_mirna_counts_after_mapping_filtered_VolcanoPlot_results_DESeq__strain__WT_A66.png
location.deseq.strain <- "../4_deseq__strain"
filename.deseq.strain <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66"
df.deseq.strain <- read.table(paste0(location.deseq.strain, "/", filename.deseq.strain, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# DESeq contrast: time (using Likelihood Ratio Test - all time points together (~ANOVA))
# DESeq2 volcano plot: summarised_mirna_counts_after_mapping_filtered_VolcanoPlot_results_DESeq__time.png
location.deseq.time <- "../4_deseq__time"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__time"
df.deseq.time <- read.table(paste0(location.deseq.time, "/", filename.deseq.time, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)





###########
# Filtering
###########

#pca
pca.thres.val <- 0.05
#deseq strain
deseq.strain.padj.thres.val <- 0.05
deseq.strain.lfc.thres.val <- 0.1
#deseq time
deseq.time.padj.thres.val <- 0.05
deseq.time.lfc.thres.val <- 0.1


# Create a pca data frame of significant loadings for PC1.
# remove rows containing NA
df.pca.filt <- df.pca[complete.cases(df.pca), ]
# filter
df.pca.filt <- df.pca[df.pca$PC1 < -pca.thres.val | df.pca$PC1 > pca.thres.val,]
# extract names
df.pca.filt.miRNA <- rownames(df.pca.filt)
# Write list of miRNAs
write.csv(df.pca.filt, file=paste0("miRNA__filtered_pca_PC1__loadings_gt_", gsub('\\.','', pca.thres.val), suffix), row.names=T, quote=FALSE)

# Create a DESeq2:strain results data frame of significant miRNA (padj, lfc)
df.deseq.strain.filt <- df.deseq.strain[complete.cases(df.deseq.strain), ]
df.deseq.strain.filt <- df.deseq.strain.filt[df.deseq.strain.filt$log2FoldChange < -deseq.strain.lfc.thres.val | df.deseq.strain.filt$log2FoldChange > deseq.strain.lfc.thres.val,]
df.deseq.strain.filt <- df.deseq.strain.filt[df.deseq.strain.filt$padj < deseq.strain.padj.thres.val,]
df.deseq.strain.filt.miRNA <- rownames(df.deseq.strain.filt)
write.csv(df.deseq.strain.filt, file=paste0("miRNA__filtered_deseq_strain__padj_", gsub('\\.','', deseq.strain.padj.thres.val), '_lfc_', gsub('\\.','', deseq.strain.lfc.thres.val), suffix), row.names=T, quote=FALSE)

# Create a DESeq2:time results data frame of significant miRNA (padj, lfc)
df.deseq.time.filt <- df.deseq.time[complete.cases(df.deseq.time), ]
df.deseq.time.filt <- df.deseq.time.filt[df.deseq.time.filt$log2FoldChange < -deseq.time.lfc.thres.val | df.deseq.time.filt$log2FoldChange > deseq.time.lfc.thres.val,]
df.deseq.time.filt <- df.deseq.time.filt[df.deseq.time.filt$padj < deseq.time.padj.thres.val,]
df.deseq.time.filt.miRNA <- rownames(df.deseq.time.filt)
write.csv(df.deseq.time.filt, file=paste0("miRNA__filtered_deseq_time__padj_", gsub('\\.','', deseq.time.padj.thres.val), '_lfc_', gsub('\\.','', deseq.time.lfc.thres.val), suffix), row.names=T, quote=FALSE)








#############
# CLUSTERING  
#############

# Run clustering
# Partition around medoids
pamx <- pam(counts, clusters, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
# Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
write.csv(pamx$clustering, file=paste0("pam_clustering_labels.csv"), quote=FALSE)


# Create vectors of labels when the miRNA is significative ('' otherwise). 

# vector of mirna with filtered PC1
mirna.pca.PC1.filt.labels <- ifelse(rownames(counts) %in% df.pca.filt.miRNA, rownames(counts), '')
# vector of mirna with significative DESeq2:Strain
mirna.deseq.strain.filt.labels <- ifelse(rownames(counts) %in% df.deseq.strain.filt.miRNA, rownames(counts), '')
# vector of mirna with significative DESeq2:Time
mirna.deseq.time.filt.labels <- ifelse(rownames(counts) %in% df.deseq.time.filt.miRNA, rownames(counts), '')

# vector of mirna with filtered PC1 and significative DESeq2:Strain
mirna.pca.PC1.deseq.strain.filt.labels <- ifelse(mirna.deseq.strain.filt.labels %in% mirna.pca.PC1.filt.labels, mirna.deseq.strain.filt.labels, '')
# vector of mirna with filtered PC1 and significative DESeq2:Time
mirna.pca.PC1.deseq.time.filt.labels <- ifelse(mirna.deseq.time.filt.labels %in% mirna.pca.PC1.filt.labels, mirna.deseq.time.filt.labels, '')



# create summary data frames to plot
df.loading.pam <- data.frame(df.pca, cluster=as.factor(pamx$clustering), check.names = FALSE)

# PLOT
# PC1 and DESeq:Strain
df.loading.pam$label <- mirna.pca.PC1.deseq.strain.filt.labels
plot_pca_clustering(df=df.loading.pam, filename=paste0('pam_clustering_w_deseq_strain_pca_PC1_plot_wlabels.png'), show.labels=TRUE)

# PC1 and DESeq:Time
df.loading.pam$label <- mirna.pca.PC1.deseq.time.filt.labels
plot_pca_clustering(df=df.loading.pam, filename=paste0('pam_clustering_w_deseq_time_pca_PC1_plot_wlabels.png'), show.labels=TRUE)






###################
# Plot Venn Diagram
###################

# PCA:PC1 vs DESeq2:strain
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_strain_VS_filt_pca_PC1",
                  category.names=c("PCA:PC1", "DESeq:strain"),
                  colours=c('green', 'magenta'))
# Save the intersection of these two sets
df.PC1.deseq.strain.miRNA <- intersect(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA)
write.csv(df.PC1.deseq.strain.miRNA, file=paste0("venn_diagram_intersect__filt_pca_PC1_VS_signif_deseq_strain", suffix), row.names=F, quote=FALSE)



# PCA:PC1 vs DESeq2:time
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.time.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_time_VS_filt_pca_PC1",
                  category.names=c("PCA:PC1", "DESeq:time"),
                  colours=c('green', 'magenta'))
# Save the intersection of these two sets
df.PC1.deseq.time.miRNA <- intersect(df.pca.filt.miRNA, df.deseq.time.filt.miRNA)
write.csv(df.PC1.deseq.time.miRNA, file=paste0("venn_diagram_intersect__filt_pca_PC1_VS_signif_deseq_time", suffix), row.names=F, quote=FALSE)



# PCA:PC1 vs DESeq2:strain vs DESeq2:time
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA, df.deseq.time.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_strain_VS_signif_deseq_time_VS_filt_pca_PC1",
                  category.names=c("PCA:PC1", "DESeq:strain", "DESeq:time"),
                  colours=c('green', 'magenta', 'yellow'))
# Save the intersection of these two sets
df.PC1.deseq.strain.time.miRNA <- intersect(df.PC1.deseq.strain.miRNA, df.PC1.deseq.time.miRNA)
write.csv(df.PC1.deseq.strain.time.miRNA, file=paste0("venn_diagram_intersect__filt_pca_PC1_VS_signif_deseq_strain_VS_signif_deseq_time", suffix), row.names=F, quote=FALSE)





#############################################################################################
# Generate PCAs for statistically significant miRNA from DESeq:strain (contrast: WT vs A66 in strain)
#############################################################################################

# The purpose is to see whether the data separation detected on PCA:PC1 using all miRNAs 
# is still present within the group of significant miRNAs.

# Create a counts table of significant miRNAs
counts.deseq.strain.signif <- counts[rownames(counts) %in% df.deseq.strain.filt.miRNA,]

pca.reads <- prcomp(counts.deseq.strain.signif, center=TRUE, scale.=TRUE)
# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca.reads)     
# Plot the variance of PCA components
png(paste0(filename.counts,"_reads_pca_comp_variances",".png"), width=2000, height=1500, res=300)
plot(pca.reads, type = "l", main='Variance of PCA components')
dev.off()
# Write the calculated PCA `rotation` (PCA load). This is the table samples (rows) vs PCAs (cols).
write.csv(pca.reads$rotation, file=paste0(filename.counts,"_reads_PCA_rotation", suffix), quote=FALSE)



# prepare the data set for plotting
## Extract STRAIN
# use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(counts.deseq.strain.signif)), 'WT', ifelse(grepl('noEGF', colnames(counts.deseq.strain.signif)), 'A66 no EGF', 'A66'))
## Extract TIME
time <- gsub(".*[_]([^.]+)[_].*", "\\1", colnames(counts.deseq.strain.signif))
# now replace 'noEGF' with '300', as noEGF sample was taken at 300m. 
time <- gsub('noEGF', '300', time)
df <- data.frame(pca.reads$rotation, 
                 strain=strain, 
                 time=as.numeric(time), 
                 check.names = FALSE)

# Plot PCA (colour is strain, shape is time)
plot_pca(df, eigen, paste0(filename.counts, '_deseq_strain_signif.png'))



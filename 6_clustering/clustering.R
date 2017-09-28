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
library(ggfortify) # autoplot
library(cluster) # pam

library(grid)
library(gridExtra)

source('../utilities/plots.R')




################
# Load data sets
################

suffix <-".csv"

# Counts Matrix
location.counts <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered_scaled"
counts.scaled <- read.table(paste0(location.counts, "/", filename.counts, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# PCA rotation matrix for counts.scaled
location.pca <- "../3_quality_control"
filename.pca <- "summarised_mirna_counts_after_mapping_filtered_scaled_PCA_rotation"
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


# Create a pca data frame of significant loadings for PC2.
# remove rows containing NA
df.pca.filt <- df.pca[complete.cases(df.pca), ]
# filter
df.pca.filt <- df.pca[df.pca$PC2 < -pca.thres.val | df.pca$PC2 > pca.thres.val,]
# extract names
df.pca.filt.miRNA <- rownames(df.pca.filt)

# Create a DESeq2:strain results data frame of significant miRNA (padj, lfc)
df.deseq.strain.filt <- df.deseq.strain[complete.cases(df.deseq.strain), ]
df.deseq.strain.filt <- df.deseq.strain.filt[df.deseq.strain.filt$log2FoldChange < -deseq.strain.lfc.thres.val | df.deseq.strain.filt$log2FoldChange > deseq.strain.lfc.thres.val,]
df.deseq.strain.filt <- df.deseq.strain.filt[df.deseq.strain.filt$padj < deseq.strain.padj.thres.val,]
df.deseq.strain.filt.miRNA <- rownames(df.deseq.strain.filt)

# Create a DESeq2:time results data frame of significant miRNA (padj, lfc)
df.deseq.time.filt <- df.deseq.time[complete.cases(df.deseq.time), ]
df.deseq.time.filt <- df.deseq.time.filt[df.deseq.time.filt$log2FoldChange < -deseq.time.lfc.thres.val | df.deseq.time.filt$log2FoldChange > deseq.time.lfc.thres.val,]
df.deseq.time.filt <- df.deseq.time.filt[df.deseq.time.filt$padj < deseq.time.padj.thres.val,]
df.deseq.time.filt.miRNA <- rownames(df.deseq.time.filt)






#############
# CLUSTERING    --- TO WORK ON THIS! 
#############

clusters <- 3
scale.=FALSE

# this is a % (the higher, the more we filter)
pca.thres <- rep(pca.thres.val, nrow(counts.scaled))
# this is padj (the lower, the more we filter)
deseq.strain.thres <- rep(deseq.strain.padj.thres.val, nrow(counts.scaled))
deseq.time.thres <- rep(deseq.time.padj.thres.val, nrow(counts.scaled))

# Unify the data sets (is there a better way to pass the thresholds without adding a new column?? aes, aes_string seem not to work..)
df.unified.strain <- data.frame(counts.scaled, miRNA=rownames(counts.scaled), df.pca, df.deseq.strain, pca_thres=pca.thres, deseq_thres=deseq.strain.thres)
df.unified.time <- data.frame(counts.scaled, miRNA=rownames(counts.scaled), df.pca, df.deseq.time, pca_thres=pca.thres, deseq_thres=deseq.time.thres)

# PLOT
# PC2 and DESeq:Strain
plot_clustering(counts.scaled, df.unified.strain, clusters, paste0('pam_clustering_w_deseq_strain_pca_pc2'), labels=FALSE, labels.col="miRNA", scale.)
plot_clustering(counts.scaled, df.unified.strain, clusters, paste0('pam_clustering_w_deseq_strain_pca_pc2'), labels=TRUE, labels.col="miRNA", scale.)

# PC2 and DESeq:Time
plot_clustering(counts.scaled, df.unified.time, clusters, paste0('pam_clustering_w_deseq_time_pca_pc2'), labels=FALSE, labels.col="miRNA", scale.)
plot_clustering(counts.scaled, df.unified.time, clusters, paste0('pam_clustering_w_deseq_time_pca_pc2'), labels=TRUE, labels.col="miRNA", scale.)







###################
# Plot Venn Diagram
###################

# PCA:PC2 vs DESeq2:strain
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_strain_VS_filt_pca_pc2",
                  category.names=c("PCA:PC2", "DESeq:strain"),
                  colours=c('green', 'magenta'))
# Save the intersection of these two sets
df.pc2.deseq.strain.miRNA <- intersect(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA)
write.csv(df.pc2.deseq.strain.miRNA, file=paste0("venn_diagram_intersect__filt_pca_pc2_VS_signif_deseq_strain", suffix), row.names=F, quote=FALSE)



# PCA:PC2 vs DESeq2:time
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.time.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_time_VS_filt_pca_pc2",
                  category.names=c("PCA:PC2", "DESeq:time"),
                  colours=c('green', 'magenta'))
# Save the intersection of these two sets
df.pc2.deseq.time.miRNA <- intersect(df.pca.filt.miRNA, df.deseq.time.filt.miRNA)
write.csv(df.pc2.deseq.time.miRNA, file=paste0("venn_diagram_intersect__filt_pca_pc2_VS_signif_deseq_time", suffix), row.names=F, quote=FALSE)



# PCA:PC2 vs DESeq2:strain vs DESeq2:time
plot_venn_diagram(data=list(df.pca.filt.miRNA, df.deseq.strain.filt.miRNA, df.deseq.time.filt.miRNA), 
                  filename="venn_diagram_intersect__signif_deseq_strain_VS_signif_deseq_time_VS_filt_pca_pc2",
                  category.names=c("PCA:PC2", "DESeq:strain", "DESeq:time"),
                  colours=c('green', 'magenta', 'yellow'))
# Save the intersection of these two sets
df.pc2.deseq.strain.time.miRNA <- intersect(df.pc2.deseq.strain.miRNA, df.pc2.deseq.time.miRNA)
write.csv(df.pc2.deseq.strain.time.miRNA, file=paste0("venn_diagram_intersect__filt_pca_pc2_VS_signif_deseq_strain_VS_signif_deseq_time", suffix), row.names=F, quote=FALSE)





######################
# Write list of miRNAs
######################

write.csv(df.pca.filt, file=paste0("miRNA__filtered_pca_pc2__loadings_gt_", gsub('\\.','', pca.thres.val), suffix), row.names=T, quote=FALSE)
write.csv(df.deseq.strain.filt, file=paste0("miRNA__filtered_deseq_strain__padj_", gsub('\\.','', deseq.strain.padj.thres.val), '_lfc_', gsub('\\.','', deseq.strain.lfc.thres.val), suffix), row.names=T, quote=FALSE)
write.csv(df.deseq.time.filt, file=paste0("miRNA__filtered_deseq_time__padj_", gsub('\\.','', deseq.time.padj.thres.val), '_lfc_', gsub('\\.','', deseq.time.lfc.thres.val), suffix), row.names=T, quote=FALSE)





#############################################################################################
# Generate PCAs for statistically significant miRNA from DESeq:strain (contrast: WT vs A66 in strain)
#############################################################################################

# The purpose is to see whether the data separation detected on PCA:PC2 using all miRNAs 
# is still present within the group of significant miRNAs.


# Create a counts table of significant miRNAs
counts.scaled.deseq.strain.signif <- counts.scaled[rownames(counts.scaled) %in% df.deseq.strain.filt.miRNA,]

# Now we prepare PCA using this data set (the data set is already scaled, so we only need to compute the transpose as we are interested 
# in the PCA for the samples)
counts.scaled.deseq.strain.signif.t <- t(counts.scaled.deseq.strain.signif)

# Now, create a data set from df.counts.signif.t with `strain` 
counts.scaled.deseq.strain.signif.t.metadata <- counts.scaled.deseq.strain.signif.t
## Extract STRAIN. Use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(counts.scaled.deseq.strain.signif)), 'WT', 'A66')

# now we add the columns
counts.scaled.deseq.strain.signif.t.metadata <- cbind(counts.scaled.deseq.strain.signif.t.metadata, strain)




# Run PCA
pca <- prcomp(counts.scaled.deseq.strain.signif.t, center=TRUE, scale.=scale.)
# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca)     

# Plot the variance of PCA components
png(paste0("pca__comp_variances__deseq_strain_signif_mirna",".png"), width=2000, height=1500, res=300)
plot(pca, type = "l", main='Variance of PCA components')
dev.off()

# Write the calculated PCA `rotation` (PCA loadings). This is the table miRNAs (rows) vs PCAs (cols).
write.csv(pca$rotation, file=paste0("pca_loadings_filt_deseq_strain_signif_mirna", suffix), quote=FALSE)




# Plot PCA (colour is strain)
# PC1 vs PC2
# plot without labels
c1c2.strain <- autoplot(pca, data=counts.scaled.deseq.strain.signif.t.metadata, colour='strain', x=1, y=2, scale=scale.) +
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()

# PC2 vs PC3
# plot without labels
c2c3.strain <- autoplot(pca, data=counts.scaled.deseq.strain.signif.t.metadata, colour='strain', x=2, y=3, scale=scale.) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1.5) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()

## COMBINED plots
c1c2c3.combined <- arrangeGrob(c1c2.strain, c2c3.strain, ncol=2)
ggsave(paste0("pca_c1c2c3__deseq_strain_signif_mirna.png"), plot=c1c2c3.combined, width=7, height=7, dpi=300)





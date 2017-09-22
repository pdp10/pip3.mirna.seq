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




# Run PAM clustering (colour is cluster)
plot_clustering <- function(df.pca, df.full, k, filename, labels=FALSE, labels.col="miRNA", scale.=FALSE) {
  
  # Partition around medoids
  pamx <- pam(df.pca, k, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
  
  # Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
  write.csv(pamx$clustering, file=paste0("pam_clustering_labels.csv"), quote=FALSE)
  
  # Plot Clustering 
  # PC1 vs PC2
  # plot without labels
  c1c2.pam <- autoplot(pamx, data=df.full, x=1, y=2, scale=scale.) +
    ggtitle('PAM clustering of miRNAs') +
    #coord_fixed(ratio=0.5) +
    theme_basic()

  # PC2 vs PC3
  # plot without labels
  c2c3.pam <- autoplot(pamx, data=df.full, x=2, y=3, scale=scale.) +
    ggtitle('PAM clustering of miRNAs') +
    #coord_fixed(ratio=0.5) +
    theme_basic()

  if(labels) {
    c1c2.pam <- c1c2.pam +
      geom_text(aes(label=ifelse( abs(PC2) > pca_thres & padj < deseq_thres, as.character(miRNA),'')),hjust=0,vjust=1, color="black", size=1.5)
    c2c3.pam <- c2c3.pam +
      geom_text(aes(label=ifelse( abs(PC2) > pca_thres & padj < deseq_thres, as.character(miRNA),'')),hjust=0,vjust=1, color="black", size=1.5)
    
    # write this list of genes:
    signif.mir <- ifelse( abs(df.full$PC2) > df.full$pca_thres & df.full$padj < df.full$deseq_thres, as.character(df.full$miRNA),NA)
    signif.mir <- signif.mir[!is.na(signif.mir)]
    write.csv(signif.mir, file=paste0(filename, "_miRNA_labels.csv"), quote=FALSE)

    ## COMBINE plots
    c1c2c3.pam.combined <- arrangeGrob(c1c2.pam, c2c3.pam, ncol=2)
    ggsave(paste0(filename, "_plot_w_labels.png"), plot=c1c2c3.pam.combined, width=8, height=3.5, dpi=300)
  } else {
    ## COMBINED plots
    c1c2c3.pam.combined <- arrangeGrob(c1c2.pam, c2c3.pam, ncol=2)
    ggsave(paste0(filename, "_plot.png"), plot=c1c2c3.pam.combined, width=8, height=3.5, dpi=300)
  }  
}





################
# Load data sets
################

suffix <-".csv"

# Counts Matrix
location.counts <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered_scaled"
counts.scaled <- read.table(paste0(location.counts, "/", filename.counts, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# PCA rotation matrix for counts.scaled
location.pca.rot <- "../3_quality_control"
filename.pca.rot <- "summarised_mirna_counts_after_mapping_filtered_scaled_PCA_rotation"
df.pca.rot <- read.table(paste0(location.pca.rot, "/", filename.pca.rot, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# DESeq contrast: strain
location.deseq.strain <- "../4_deseq__strain"
filename.deseq.strain <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66"
df.deseq.strain <- read.table(paste0(location.deseq.strain, "/", filename.deseq.strain, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)


# DESeq contrast: time
location.deseq.time <- "../4_deseq__time"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__time"
df.deseq.time <- read.table(paste0(location.deseq.time, "/", filename.deseq.time, suffix), sep=",", fill=TRUE, header=TRUE, row.names=1)




# this is a % (the higher, the more we filter)
pca.thres <- rep(0.05, nrow(counts.scaled))
# this is padj (the lower, the more we filter)
deseq.strain.thres <- rep(0.05, nrow(counts.scaled))
deseq.time.thres <- rep(0.05, nrow(counts.scaled))




# Unify the data sets (is there a better way to pass the thresholds without adding a new column?? aes, aes_string seem not to work..)
df.unified.strain <- data.frame(counts.scaled, miRNA=rownames(counts.scaled), df.pca.rot, df.deseq.strain, pca_thres=pca.thres, deseq_thres=deseq.strain.thres)
df.unified.time <- data.frame(counts.scaled, miRNA=rownames(counts.scaled), df.pca.rot, df.deseq.time, pca_thres=pca.thres, deseq_thres=deseq.time.thres)









#############
# CLUSTERING
#############

clusters <- 2
scale.=FALSE

# PC2 and DESeq:Strain
plot_clustering(counts.scaled, df.unified.strain, clusters, paste0('pam_clustering_w_deseq_strain_pca_pc2'), labels=FALSE, labels.col="miRNA", scale.)
plot_clustering(counts.scaled, df.unified.strain, clusters, paste0('pam_clustering_w_deseq_strain_pca_pc2'), labels=TRUE, labels.col="miRNA", scale.)

# PC2 and DESeq:Time
plot_clustering(counts.scaled, df.unified.time, clusters, paste0('pam_clustering_w_deseq_time_pca_pc2'), labels=FALSE, labels.col="miRNA", scale.)
plot_clustering(counts.scaled, df.unified.time, clusters, paste0('pam_clustering_w_deseq_time_pca_pc2'), labels=TRUE, labels.col="miRNA", scale.)





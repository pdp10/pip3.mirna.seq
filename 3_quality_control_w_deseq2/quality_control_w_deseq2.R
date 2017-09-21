
# This script runs some quality control for the counts table
# Most of the quality control is performed on normalised counts tables.
# These tables are generated only if they don't already exist.


# download DESeq2 from: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library("DESeq2")
library('vsn')
library("dplyr")
library("ggplot2")
library(grid)
library(gridExtra)
library("hexbin")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")

source('../utilities/plots.R')


# select the file countaining the data
location <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered"
filename.counts.metadata <- "summarised_mirna_counts_after_mapping_filtered_metadata"
suffix <-".csv"



#####################################
# Read Counts Table and Samples Table
#####################################

# load counts
counts <- read.table(paste0(location,"/",filename.counts,suffix), sep=",",fill=T,header=T,row.names=1)

# load counts metadata
counts.metadata <- read.table(paste0(location,"/",filename.counts.metadata,suffix), sep=",",fill=T,header=T,row.names=1)
# cast the column `time` from numeric to factor
counts.metadata$time <- as.factor(counts.metadata$time)
# cast the column `replicate` from numeric to factor
counts.metadata$replicate <- as.factor(counts.metadata$replicate)


#######################
# Prepare DESeq dataset
#######################
# create the dataset for DESeq. We use this as an annotated SummarizedExperiment object
# use assay(se), colData(se), design(se) to read the most important info
se <- DESeqDataSetFromMatrix(countData = counts, colData=counts.metadata, design = ~ 1)



######################
# Estimate dispersions
######################
# estimate the size factors from counts, counts.metadata
se <- estimateSizeFactors(se)
se <- estimateDispersions(se)
#str(fitInfo(se))
png(paste0(filename.counts,"_DESeq_dispersions_estimates",".png"), width=2000, height=2000, pointsize=14, res=300)
plotDispEsts(se)
dev.off()



############################
# Normalise the counts table
############################

# See workflow: http://www.bioconductor.org/help/workflows/rnaseqGene/ 

# 1) transformations of the counts in order to visually explore sample relationships. 
# 2) original raw counts for statistical testing. 
# => This is critical because the statistical testing methods rely on original count data (not scaled or transformed) 
# for calculating the precision of measurements.


# Many common statistical methods for exploratory analysis of multidimensional data, 
# for example clustering and principal components analysis (PCA), work best for data 
# that generally has the same range of variance at different ranges of the mean values. 
# When the expected amount of variance is approximately the same across different mean 
# values, the data is said to be homoskedastic. For RNA-seq counts, however, the expected 
# variance grows with the mean. For example, if one performs PCA directly on a matrix of 
# counts or normalized counts (e.g. correcting for differences in sequencing depth), the 
# resulting plot typically depends mostly on the genes with highest counts because they 
# show the largest absolute differences between samples. A simple and often used strategy 
# to avoid this is to take the logarithm of the normalized count values plus a pseudocount of 1; 
# however, depending on the choice of pseudocount, now the genes with the very lowest counts 
# will contribute a great deal of noise to the resulting plot, because taking the logarithm 
# of small counts actually inflates their variance.


# Plot SD vs mean
# We plot the standard deviation of each row (genes) against the mean
p<-meanSdPlot(assay(se), ranks=FALSE)
ggsave(paste0(filename.counts, "_meanSdPlot_counts.png"), width=4, height=4, dpi=300)

# The logarithm with a small pseudocount amplifies differences when the values 
# are close to 0. The low count genes with low signal-to-noise ratio will overly 
# contribute to sample-sample distances and PCA plots.
# apply a log(counts + 1) transformation
norm.log<-normTransform(se)
meanSdPlot(assay(norm.log), ranks = FALSE)
ggsave(paste0(filename.counts, "_meanSdPlot_log_counts.png"), width=4, height=4, dpi=300)


# The logarithm with a small pseudocount amplifies differences when the values are close to 0. 
# The low count genes with low signal-to-noise ratio will overly contribute to sample-sample distances and PCA plots.
# As a solution, DESeq2 offers two transformations for count data that stabilize the variance 
# across the mean: the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014), 
# and the variance stabilizing transformation (VST) for negative binomial data with a 
# dispersion-mean trend (Anders and Huber 2010), implemented in the vst function.

# For genes with high counts, the rlog and VST (see below) will give similar result to the ordinary log2 
# transformation of normalized counts. 
# For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples. 
# The rlog-transformed or VST data then becomes approximately homoskedastic, and can be used directly for 
# computing distances between samples, making PCA plots, or as input to downstream methods which perform 
# best with homoskedastic data.
# Which transformation to choose? 
# The rlog tends to work well on small datasets (n < 30), sometimes outperforming the VST when there is a 
# large range of sequencing depth across samples (an order of magnitude difference). The VST is much faster 
# to compute and is less sensitive to high count outliers than the rlog. We therefore recommend the VST for 
# large datasets (hundreds of samples). You can perform both transformations and compare the meanSdPlot or 
# PCA plots generated, as described below.
# Note that the two transformations offered by DESeq2 are provided for applications other than differential 
# testing. For differential testing we recommend the DESeq function applied to raw counts, as described later 
# in this workflow, which also takes into account the dependence of the variance of counts on the mean value 
# during the dispersion estimation step.

# As alternative, compute rlog as implemented in DESeq2. The rlog tends to work well on small datasets (n < 30)
# In the below function calls, we specified blind = FALSE, which means that differences between cell lines and 
# treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
# The experimental design is not used directly in the transformation, only in estimating the global amount of 
# variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
rld <- rlog(se, blind = FALSE)
vst <- varianceStabilizingTransformation(se, blind = FALSE)

df <- bind_rows(
  #as_data_frame(log2(counts(se, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(norm.log)[, 1:2]) %>% mutate(transformation = "log2(x+1)"),  
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vst)[, 1:2]) %>% mutate(transformation = "vst"))
colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)  
ggsave(paste0(filename.counts, "_norm_compare.png"), width=5, height=2, dpi=300)



##################################
# Write normalised tables to files
###################################

# convert rownames to a new column
norm.log.counts <- assay(norm.log)
norm.log.counts <- cbind(miRNA = rownames(norm.log.counts), norm.log.counts)
write.csv(norm.log.counts,file=paste0(location,"/",filename.counts,"_norm_log",suffix),quote=FALSE,row.names=FALSE)

# convert rownames to a new column
rld.counts <- assay(rld)
rld.counts <- cbind(miRNA = rownames(rld.counts), rld.counts)
write.csv(rld.counts,file=paste0(location,"/",filename.counts,"_norm_rlog",suffix),quote=FALSE,row.names=FALSE)

# convert rownames to a new column
vst.counts <- assay(vst)
vst.counts <- cbind(miRNA = rownames(vst.counts), vst.counts)
write.csv(vst.counts,file=paste0(location,"/",filename.counts,"norm_vst",suffix),quote=FALSE,row.names=FALSE)




##################
# Quality Control
##################

# Sample distances 
# A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: 
# Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?
# We use the R function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the rlog-transformed data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples 
# to be rows of its argument, and different dimensions (here, genes) to be columns.

# In order to plot the sample distance matrix with the rows/columns arranged by the distances 
# in our distance matrix, we manually provide sampleDists to the clustering_distance argument 
# of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains 
# the data values themselves, and would calculate distances between the rows/columns of the distance 
# matrix, which is not desired. We also manually specify a blue color palette using the colorRampPalette 
# function from the RColorBrewer package.

# norm.log
sampleDists <- dist(t(assay(norm.log)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( norm.log$strain, norm.log$time, sep=" - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(paste0(filename.counts,"_strain_time_sample_distances_norm_log",".png"), width=2000, height=2000, res=300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colors)
dev.off()

# rlog
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
#print(sampleDistMatrix)
rownames(sampleDistMatrix) <- paste( rld$strain, rld$time, sep=" - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(paste0(filename.counts,"_strain_time_sample_distances_norm_rlog",".png"), width=2000, height=2000, res=300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colors)
dev.off()

# vst
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vst$strain, vst$time, sep=" - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(paste0(filename.counts,"_strain_time_sample_distances_norm_vst",".png"), width=2000, height=2000, res=300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colors)
dev.off()




# Heatmap of sample-to-sample distances using the original count matrix values and PoissonDistance.
# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), 
# implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes 
# the inherent variance structure of counts into consideration when calculating the distances 
# between samples. The PoissonDistance function takes the original count matrix (not normalized) 
# with samples as rows instead of columns, so we need to transpose the counts in se.
poisd <- PoissonDistance(t(counts(se)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$strain, rld$time, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
png(paste0(filename.counts,"_strain_time_sample_distances_PoissonDist",".png"), width=2000, height=2000, res=300)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()




# PCA with log2(x+1), rlog or vst (VarianceStabilizingTransformation). vst is generally used for datasets with many samples
# log

# we redefine the plot so that we can change the theme. To do so, we need to return the data.
# This plot redesign is done by the functions plotPrettyPCA1Var() and plotPrettyPCA2Var(). 
# p1 <- plotPCA(norm.log, intgroup = c("time"))
p1 <- plotPrettyPCA1Var(norm.log, intgroup = c("time"), plotColor=FALSE)
#ggsave(paste0(filename.counts, "_pca_time_norm_log.png"), width=4, height=4, dpi=300)

p2 <- plotPrettyPCA1Var(norm.log, intgroup = c("strain"))
#ggsave(paste0(filename.counts, "_pca_strain_norm_log.png"), width=4, height=4, dpi=300)

p3 <- plotPrettyPCA2Var(norm.log, intgroup = c("strain", "time"))
#ggsave(paste0(filename.counts, "_pca_time_strain_norm_log.png"), width=4, height=4, dpi=300)

p.combined <- arrangeGrob(p1,p2,p3,ncol=2)
ggsave(paste0(filename.counts, "_pca_combined_norm_log.png"), plot=p.combined, width=7, height=7, dpi=300)



# rlog
p1 <- plotPrettyPCA1Var(rld, intgroup = c("time"), plotColor=FALSE)
#ggsave(paste0(filename.counts, "_pca_time_norm_rlog.png"), width=4, height=4, dpi=300)

p2 <- plotPrettyPCA1Var(rld, intgroup = c("strain"))
#ggsave(paste0(filename.counts, "_pca_strain_norm_rlog.png"), width=4, height=4, dpi=300)

p3 <- plotPrettyPCA2Var(rld, intgroup = c("strain", "time"))
#ggsave(paste0(filename.counts, "_pca_time_strain_norm_rlog.png"), width=4, height=4, dpi=300)

p.combined <- arrangeGrob(p1,p2,p3,ncol=2)
ggsave(paste0(filename.counts, "_pca_combined_norm_rlog.png"), plot=p.combined, width=7, height=7, dpi=300)



# vst
p1 <- plotPrettyPCA1Var(vst, intgroup = c("time"), plotColor=FALSE)
#ggsave(paste0(filename.counts, "_pca_time_norm_vst.png"), width=4, height=4, dpi=300)

p2 <- plotPrettyPCA1Var(vst, intgroup = c("strain"))
#ggsave(paste0(filename.counts, "_pca_strain_norm_vst.png"), width=4, height=4, dpi=300)

p3 <- plotPrettyPCA2Var(vst, intgroup = c("strain", "time"))
#ggsave(paste0(filename.counts, "_pca_time_strain_norm_vst.png"), width=4, height=4, dpi=300)

p.combined <- arrangeGrob(p1,p2,p3,ncol=2)
ggsave(paste0(filename.counts, "_pca_combined_norm_vst.png"), plot=p.combined, width=7, height=7, dpi=300)


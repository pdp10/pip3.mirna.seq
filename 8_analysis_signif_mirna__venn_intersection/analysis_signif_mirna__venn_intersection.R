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
source('../utilities/clustering.R')


################
# Load data sets
################

suffix <-".csv"
clusters <- 3

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
# FILTER counts matrix of significant miRNA
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







##################
# Plot time course
##################

######################################################
# Separate counts.median.scaled.venn.intersect into three data frames
######################################################

# Separate miWT, miA66, miA66_noEGF. We also add miRNA names as a new column (we need this for melting)
counts.scaled.venn.intersect.miWT <- data.frame(miRNA=rownames(counts.scaled.venn.intersect), counts.median.scaled.venn.intersect[, grepl("miWT", colnames(counts.median.scaled.venn.intersect)), drop=FALSE])
counts.scaled.venn.intersect.miA66 <- data.frame(miRNA=rownames(counts.scaled.venn.intersect), counts.median.scaled.venn.intersect[, grepl("miA66", colnames(counts.median.scaled.venn.intersect)), drop=FALSE])
counts.scaled.venn.intersect.miA66.noEGF <- data.frame(miRNA=rownames(counts.scaled.venn.intersect), counts.median.scaled.venn.intersect[, grepl("noEGF", colnames(counts.median.scaled.venn.intersect)), drop=FALSE])

# remove column `miA66_noEGF`` from miA66 data frame as this is treated separately.
counts.scaled.venn.intersect.miA66 <- subset(counts.scaled.venn.intersect.miA66, select = -miA66_noEGF)

# rename columns so that they only contain times (useful for melting later on)
colnames(counts.scaled.venn.intersect.miWT) <- gsub('miWT_', '', colnames(counts.scaled.venn.intersect.miWT))
colnames(counts.scaled.venn.intersect.miA66) <- gsub('miA66_', '', colnames(counts.scaled.venn.intersect.miA66))
# this is renamed manually
colnames(counts.scaled.venn.intersect.miA66.noEGF) <- c('miRNA', '300')


###################################################
# Plot all miRNA time courses by WT, A66, and noEGF
###################################################

# plot
plot_expr_tc(counts.scaled.venn.intersect.miWT, "miRNA_tc__miWT", "(WT)", line=TRUE)
plot_expr_tc(counts.scaled.venn.intersect.miA66, "miRNA_tc__miA66", "(A66)", line=TRUE)
plot_expr_tc(counts.scaled.venn.intersect.miA66.noEGF, "miRNA_tc__miA66noEGF", "(A66 noEGF)", line=FALSE)


#########################
# ADD PAM clustering info
#########################

# `scale` and `centre` must be FALSE because we already scaled and centred our data.
scale. <- FALSE
centre <- FALSE

# run PAM clustering for WT, A66, A66noEGF
pam.miWT <- run_pam(subset(counts.scaled.venn.intersect.miWT, select = -c(miRNA), drop=FALSE), clusters, "median_counts_WT_samples__pam", scale., centre)
pam.miA66 <- run_pam(subset(counts.scaled.venn.intersect.miA66, select = -c(miRNA), drop=FALSE), clusters, "median_counts_A66_samples__pam", scale., centre)
pam.miA66.noEGF <- run_pam(subset(counts.scaled.venn.intersect.miA66.noEGF, select = -c(miRNA), drop=FALSE), clusters, "median_counts_A66noEGF_samples__pam", scale., centre)

# add PAM labels
counts.scaled.venn.intersect.miWT.pam <- data.frame(counts.scaled.venn.intersect.miWT, pam=factor(pam.miWT$clustering), check.names=FALSE)
counts.scaled.venn.intersect.miA66.pam <- data.frame(counts.scaled.venn.intersect.miA66, pam=factor(pam.miA66$clustering), check.names=FALSE)
counts.scaled.venn.intersect.miA66.noEGF.pam <- data.frame(counts.scaled.venn.intersect.miA66.noEGF, pam=factor(pam.miA66.noEGF$clustering), check.names=FALSE)

# plot
plot_expr_tc_w_pam(counts.scaled.venn.intersect.miWT.pam, "miRNA_tc__miWT_pam_clust", "(WT)", line=TRUE)
plot_expr_tc_w_pam(counts.scaled.venn.intersect.miA66.pam, "miRNA_tc__miA66_pam_clust", "(A66)", line=TRUE)
plot_expr_tc_w_pam(counts.scaled.venn.intersect.miA66.noEGF.pam, "miRNA_tc__miA66noEGF_pam_clust", "(A66 noEGF)", line=FALSE)


#####################################################################################################
# Median of the miRNAs belonging to the same cluster class. This for miWT, miA66, and miA66.noEGF
#####################################################################################################

# plot
plot_median_expr_tc_w_pam(subset(counts.scaled.venn.intersect.miWT.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miWT_pam_clust_median', '(WT)', line=TRUE)
plot_median_expr_tc_w_pam(subset(counts.scaled.venn.intersect.miA66.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miA66_pam_clust_median', '(A66)', line=TRUE)
plot_median_expr_tc_w_pam(subset(counts.scaled.venn.intersect.miA66.noEGF.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc_miA66noEGF_pam_clust_median', '(A66 noEGF)', line=FALSE)




















############################### Add information from DESeq2  #################################




################
# Load data sets
################

# select the files containing the data
location.time <- "../4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66_at_t"
suffix <-".csv"
padj.deseq.time <- "__padj_005"
time <- list("0", "15", "40", "90", "180", "300")


location.strain <- "../4_deseq__strain"
filename.deseq.strain <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66"
suffix <-".csv"
padj.strain.thres <- 0.01




#############################
# Merge time course data sets
#############################

# Load one DESeq data set in order to get the miRNA names.
miRNAs <- row.names(read.table(paste0(location.time,"/",filename.deseq.time, time[1], suffix), sep=",",fill=T,header=T,row.names=1))

# Add log2FoldChange_time0
df.tc <- data.frame(t0=rep(0, length(miRNAs)))
# Add row names
rownames(df.tc) <- miRNAs

# Add log2FoldChange for the available time points.
for(i in seq(1, length(time))) {
  df.tmp <- read.table(paste0(location.time,"/",filename.deseq.time, time[i], suffix), sep=",",fill=T,header=T,row.names=1)[,2]
  # Add log2FoldChange_timeX
  df.tc <- cbind(df.tc, df.tmp)
}

# Remove the first fake column t0.
df.tc <- subset(df.tc, select = -c(t0) )

# Set column names (time points)
colnames(df.tc) <- c(paste0(time))





#############################################
# Extract DESeq:strain for colour information
#############################################

# Load strain information
df.strain <- read.table(paste0(location.strain,"/",filename.deseq.strain, suffix), sep=",",fill=T,header=T,row.names=1)

# Extract the names of the significant miRNA. This will be used as a filter later on.
df.strain.signif.mirna <- rownames(subset(df.strain, padj < padj.strain.thres))
#write.csv(df.strain.signif.mirna, file=paste0(filename.deseq.strain ,"_DESeq_strain__WT_A66__padj_", gsub("\\.","",padj.strain.thres), suffix))


# Extract the log2 fold change
df.strain <- df.strain[,2, drop=FALSE]





###############################
# FILTER with significant miRNA
###############################

# df.tc
df.tc.venn.intersect <- subset(df.tc, rownames(df.tc) %in% miRNA.venn.intersect)
# df.strain
df.strain.venn.intersect <- subset(df.strain, rownames(df.strain) %in% miRNA.venn.intersect)
df.strain.venn.intersect <- df.strain.venn.intersect[,1]



##################
# Melt data frames
##################

df <- df.tc.venn.intersect
# add miRNA as new row
df$miRNA <- rownames(df)
# Add strain (colour)
df$strain <- df.strain.venn.intersect
# Melt times so that we have 'miRNA', 'strain', 'variable' (time-cols), 'value' (time-vals)
df.melt <- melt(df, id=c('miRNA','strain'))


# Due to the high number of miRNA, we also filter those having significant change in DESeq contrast `Strain`.
df.signif <- subset(df, rownames(df) %in% df.strain.signif.mirna)
df.melt.signif <- melt(df.signif, id=c('miRNA','strain'))



#####################################
# Plot miRNA using Strain for colours
#####################################

# Plot all miRNAs
g <- ggplot(data=df.melt, aes(x=variable, y=value, group=miRNA, color=strain)) + 
  geom_line() + 
  scale_colour_gradient2(low="red", mid="lightgrey", high="navyblue") +
  ggtitle('miRNA expression with strain log2 FC') +
  xlab('Time [m]') + 
  ylab('log2 fold change') +
  theme_basic()
ggsave(paste0("miRNA_log_fc_tc__deseq_strain_time_factor__w_strain.png"), width=4, height=4, dpi=300)


# Plot miRNAs with significant change in Strain
g <- ggplot(data=df.melt.signif, aes(x=variable, y=value, group=miRNA, color=strain)) +
  geom_line() +
  scale_colour_gradient2(low="red", mid="lightgrey", high="navyblue") +
  ggtitle('miRNA expression with strain log2 FC') +
  xlab('Time [m]') +
  ylab('log2 fold change') +
  theme_basic()
ggsave(paste0("miRNA_log_fc_tc__deseq_strain_time_factor__w_strain__padj_", gsub('\\.','',as.character(padj.strain.thres)), ".png"), width=4, height=4, dpi=300)


























# 
# 
# # Basic time course
# 
# # add miRNA names as a column
# counts.median.scaled.venn.intersect.mirna <- data.frame(miRNA=rownames(counts.median.scaled.venn.intersect), 
#                                                         counts.median.scaled.venn.intersect)
# 
# # rename columns so that they only contain times (useful for melting later on)
# colnames(counts.median.scaled.venn.intersect.mirna) <- gsub('miWT_', 'wt', colnames(counts.median.scaled.venn.intersect.mirna))
# colnames(counts.median.scaled.venn.intersect.mirna) <- gsub('miA66_', 'a66', colnames(counts.median.scaled.venn.intersect.mirna))
# colnames(counts.median.scaled.venn.intersect.mirna) <- gsub('noEGF', 'noEGF', colnames(counts.median.scaled.venn.intersect.mirna))
# 
# # plot
# plot_expr_tc2(counts.median.scaled.venn.intersect.mirna, "miRNA_tc", "(signif)", line=TRUE)
# 
# 
# 
# 
# # ADD PAM clustering info
# 
# clusters = 3
# # `scale` and `centre` must be FALSE because we already scaled and centred our data.
# scale. <- FALSE
# centre <- FALSE
# 
# # run PAM clustering for WT, A66, A66noEGF
# pam.miRNA <- run_pam(counts.median.scaled.venn.intersect, clusters, "signif miRNA__pam", scale., centre)
# 
# # add PAM labels
# counts.median.scaled.venn.intersect.mirna.pam <- data.frame(counts.median.scaled.venn.intersect.mirna, pam=factor(pam.miRNA$clustering), check.names=FALSE)
# 
# # plot
# plot_expr_tc_w_pam(counts.median.scaled.venn.intersect.mirna.pam, "miRNA_tc__pam_clust", "(signif)", line=TRUE)
# 
# 
# 
# # Median of the miRNAs belonging to the same cluster class. 
# plot_median_expr_tc_w_pam(subset(counts.median.scaled.venn.intersect.mirna.pam, select = -c(miRNA), drop=FALSE), 'miRNA_tc__pam_clust_median', '(signif)', line=TRUE)
# 
# 

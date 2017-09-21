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




######################################################
# DESeq2-based comparison of strain (A66 vs WT) vs time
######################################################


# This version uses the Likelihood Ratio Test. 
# "For analyses using the likelihood ratio test (using nbinomLRT), 
# the p-values are determined solely by the difference in deviance 
# between the full and reduced model formula. A log2 fold change is 
# included, which can be controlled using the name argument."


# If you want to test comparisons of pairs of levels of disease, you
# should be using the Wald test (the standard DESeq() workflow), instead
# of the LRT. The LRT is for testing all levels of disease at once.





## Issue, padj is very high.. 



# download DESeq2 from: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library('DESeq2')
library('ggplot2')
library('pheatmap')

source('../utilities/plots.R')


# Questions: 
# (1) what are the genes differentially expressed at each time for each treatment versus control 
# (e.g. treated versus control at time 15, treated versus control at time 40 ...etc)? 
# (2) What are the genes differentially expressed at each time for one treatment versus the other?




# select the file countaining the data
location <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered"
filename.counts.metadata <- "summarised_mirna_counts_after_mapping_filtered_metadata"
suffix <-".csv"
padj.thres <- 0.05
lfc.thres <- 1

par(mar=c(5,5,5,5))


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

# we will override the design later.
dds <- DESeqDataSetFromMatrix(countData=counts, colData=counts.metadata, design = ~ strain + time + strain:time)

# combine the factors of interest into a single factor with all combinations of the original factors
#dds$id <- factor(paste0(dds$strain, dds$time))
# add an id column to dds
#dds$id <- c(paste(dds$strain, dds$time, dds$egf, sep='_'))


###########
# Run DESeq
###########

# Additional filters for dds


# run DESeq using the likelihood ratio test, where we remove the strain-specific differences over time. mi-RNA with small p values 
# from this test are those which at one or more time points after time 0 showed a strain-specific effect. Notice: this will not give 
# small p values to mi-RNA that moved up or down over time in the same way in both A66/WT.

dds <- DESeq(dds, test="LRT", reduced = ~ strain + time)

# VLADIMIR (https://github.com/wikiselev/rnaseq.mcf10a/blob/master/data-raw/diff-expr.R)  !!!! Different results too..
#dds <- DESeq(dds, test="LRT", reduced = ~ time)


dds.res <- results(dds)

# write results
write.csv(dds.res[order(dds.res$padj),], file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_per_tp", suffix))

# write significant results 
# Notice: I need not less than 0.5 to filter something here...
dds.res.signif <- subset(dds.res, padj < padj.thres) 
write.csv(dds.res.signif[order(dds.res.signif$padj),], file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_per_tp__padj_", gsub("\\.","",padj.thres), suffix))




##############################################
# Analysis: WT vs A66 in strain per time point
##############################################

times <- c('0','15','40','90','180','300')

# plot time course of plot counts
plotcounts.data <- plotCounts(dds, which.min(dds.res$padj), intgroup = c("time","strain"), returnData=TRUE)
ggplot(plotcounts.data, aes(x=as.numeric(time), y=count, color=strain, group=strain)) + 
  geom_point() + 
  geom_smooth(se=FALSE, method="loess") + 
  scale_y_log10(breaks = c(50,100,200,400)) +
  ggtitle('DE time course of strains: WT vs A66') +
  scale_x_discrete(name ="Time [min]", limits=times) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(paste0(filename.counts, "_timecourse_strain_vs_time.png"), width=5, height=2, dpi=300)



# We use Wald tests for the log2 fold changes at individual time points
#> resultsNames(dds)
#  [1] "Intercept"        "strain_WT_vs_A66" "time_15_vs_0"     "time_40_vs_0"     "time_90_vs_0"    
#  [6] "time_180_vs_0"    "time_300_vs_0"    "strainWT.time15"  "strainWT.time40"  "strainWT.time90" 
# [11] "strainWT.time180" "strainWT.time300"

# the interaction terms are the difference between the two groups at a given time 
# after accounting for the difference at time 0.So, we remove 0 (the first element)
dds.res.min.padj <- which.min(dds.res$padj)

for(t in times[-1]) {
  name <- paste0('strainWT.time', t)
  dds.res.time <- results(dds, name=name, test="Wald")
  dds.res.time[dds.res.min.padj,]
  # write results
  write.csv(dds.res.time[order(dds.res.time$padj),], file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_at_t", t, suffix))
}




# Clustering of significant mi-RNA by their profiles. 
# We extract a matrix of the shrunken log2 fold changes using the coef function:
betas <- coef(dds)
colnames(betas)
# We can now plot the log2 fold changes in a heatmap.
topGenes <- head(order(dds.res$padj),50)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
png(paste0(filename.counts,"_strain_time_miRNA_clustering",".png"), width=2000, height=2000, res=300)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)
dev.off()






########################
# Generate Volcano plots
########################

# This works with p-value. But shouldn't our results be considered by padj to adjust pvalue and reduce false positives?
# results Strain vs Time
png(paste0(filename.counts,"_results_DESeq__strain_time_volcano.png"), width=2000, height=2000, pointsize=18, res=300)
signif.genes <- volcano_plot(dds.res, title='Strain vs Time', thres.padj=padj.thres, thres.lfc=lfc.thres)
dev.off()
write.csv(signif.genes, file=paste0(filename.counts ,"_results_DESeq__strain_time__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",lfc.thres), suffix))



###################
# Generate MA-plots
###################
# An MA plot is an application of a Bland–Altman plot for visual representation of genomic data.
# The plot visualises the differences between measurements taken in two samples, by transforming the
# data onto M (log ratio) and A (mean average) scales, then plotting these values.
# Here we test for differential expression between conditions Strain over Time.
# The red colour marks genes (here miRNAs) detected as differentially expressed at 10% false discovery rate
# when Benjamini-Hochberg multiple testing adjustment is used. The symbols at the upper and lower plot border
# indicate genes (here miRNAs) with very large or infinite log fold change.

# Plot mean expression vs log fold change for strain over time
png(paste0(filename.counts,"_MAplot_results_DESeq__strain_time.png"), width=2000, height=2000, pointsize=18, res=300)
plotMA(dds.res, main="Strain vs Time")
dev.off()



##############################
# Generate p-values histograms
##############################
   
# histogram of p-values for the contrast WT vs A66 per time point
png(paste0(filename.counts,"_pvalues_results_DESeq__strain_time.png"), width=2000, height=2000, pointsize=18, res=300)
hist(dds.res$pvalue, breaks=100, col="skyblue", border="slateblue", main="p-values for Strain vs Time")
dev.off()


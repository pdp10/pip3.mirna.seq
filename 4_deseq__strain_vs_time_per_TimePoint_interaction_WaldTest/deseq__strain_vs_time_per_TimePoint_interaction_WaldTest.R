
######################################################
#DESeq2-based comparison of strain (A66 vs WT) vs time
######################################################

# Group-specific condition effects, individuals nested within groups

## For details, see:
## http://master.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#nested-indiv 
## https://support.bioconductor.org/p/92932/#93108 


# download DESeq2 from: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library('DESeq2')
library('ggplot2')

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
lfc.thres <- 0.5

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

# Analysis of group-specific condition effects, while controlling for differences in individual. 
# NOTICE: 
# replicate must "restart" for each strain (group). Our `replicate` corresponds to `ind.n`, NOT `ind`
##    grp ind cnd ind.n
## 1    X   1   A     1
## 2    X   1   B     1
## 3    X   2   A     2
## 4    X   2   B     2
## 5    X   3   A     3
## 6    X   3   B     3
## 7    Y   4   A     1
## 8    Y   4   B     1
## 9    Y   5   A     2
## 10   Y   5   B     2
## 11   Y   6   A     3
## 12   Y   6   B     3
dds <- DESeqDataSetFromMatrix(countData=counts, colData=counts.metadata, design = ~ strain + strain:replicate + strain:time)

# Our DESeqDataSet has a design like ~ grp + grp:ind.n + grp:cnd. 
# This new design will result in the following model matrix:
#
model.matrix(~ strain + strain:replicate + strain:time, colData(dds))



###########
# Run DESeq
###########

# Additional filters for dds


# run DESeq
dds.res <- DESeq(dds)



###########################################
# Analysis: WT vs A66 in strain per time point
###########################################


times <- c('15','40','90','180','300')

for(time in times) {
  
  # results WT vs A66 in strain per time point
  results.strain.WT_A66 <- results(dds.res, contrast=list(paste0("strainWT.time",time),paste0("strainA66.time", time)))
  results.signif.strain.WT_A66 <- subset(results.strain.WT_A66, padj < padj.thres)
  write.csv(results.strain.WT_A66, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_at_t", time, suffix))
  write.csv(results.signif.strain.WT_A66, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_at_t", time, "__padj_", gsub("\\.","",padj.thres), suffix))
  

  ########################
  # Generate Volcano plots 
  ########################
  
  # results WT vs A66 in strain per time point
  png(paste0(filename.counts,"_VolcanoPlot_results_DESeq__strain__WT_A66_at_t", time, ".png"), width=2000, height=2000, pointsize=18, res=300)
  signif.genes <- volcano_plot(results.strain.WT_A66, title=paste0('WT vs A66 at ', time, 'm'), thres.padj=padj.thres, thres.lfc=lfc.thres)
  dev.off()
  write.csv(signif.genes, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66_at_t", time, "__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",lfc.thres), suffix))
  

  ###################
  # Generate MA-plots
  ###################
  # An MA plot is an application of a Bland–Altman plot for visual representation of genomic data. 
  # The plot visualises the differences between measurements taken in two samples, by transforming the 
  # data onto M (log ratio) and A (mean average) scales, then plotting these values.
  # Here we test for differential expression between conditions (WT vs A66). 
  # The red colour marks genes (here miRNAs) detected as differentially expressed at 10% false discovery rate 
  # when Benjamini-Hochberg multiple testing adjustment is used. The symbols at the upper and lower plot border 
  # indicate genes (here miRNAs) with very large or infinite log fold change.
  
  # Plot mean expression vs log fold change for the contrast WT vs A66 per time point
  png(paste0(filename.counts,"_MAplot_results_DESeq__strain__WT_A66_at_t", time, ".png"), width=2000, height=2000, pointsize=18, res=300)
  plotMA(results.strain.WT_A66, main=paste("Shrunk LFC at ", time, "m", sep=""))
  dev.off()

  
  
  ##############################
  # Generate p-values histograms
  ##############################
  
  # histogram of p-values for the contrast WT vs A66 per time point
  png(paste0(filename.counts,"_pvalues_results_DESeq__strain__WT_A66_at_t", time,".png"), width=2000, height=2000, pointsize=18, res=300)
  hist(results.strain.WT_A66$pvalue, breaks=100, col="skyblue", border="slateblue", main=paste("p-values for WT vs A66 at ", time, "m", sep=""))
  dev.off()
  
}








#############################################
#############################################
# Now analyse WT_vs_A66
#############################################
#############################################


# results WT vs A66 in strain
results.strain.WT_A66 <- results(dds.res, name="strain_WT_vs_A66")
results.signif.strain.WT_A66 <- subset(results.strain.WT_A66, padj < padj.thres)
write.csv(results.strain.WT_A66, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66", suffix))
write.csv(results.signif.strain.WT_A66, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66__padj_", gsub("\\.","",padj.thres), suffix))


########################
# Generate Volcano plots 
########################

# results WT vs A66 in strain
png(paste0(filename.counts,"_VolcanoPlot_results_DESeq__strain__WT_A66.png"), width=2000, height=2000, pointsize=18, res=300)
signif.genes <- volcano_plot(results.strain.WT_A66, title=paste0('WT vs A66'), thres.padj=padj.thres, thres.lfc=lfc.thres)
dev.off()
write.csv(signif.genes, file=paste0(filename.counts ,"_results_DESeq__strain__WT_A66__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",lfc.thres), suffix))


###################
# Generate MA-plots
###################
# Plot mean expression vs log fold change for the contrast WT vs A66
png(paste0(filename.counts,"_MAplot_results_DESeq__strain__WT_A66.png"), width=2000, height=2000, pointsize=18, res=300)
plotMA(results.strain.WT_A66, main=paste("Shrunk LFC at ", time, "m", sep=""))
dev.off()



######################################################
#DESeq2-based comparison of strain (A66 vs WT) vs time
######################################################

# Notice: the design does not take into account the interaction between time and strain.



# download DESeq2 from: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library("DESeq2")

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
lfc.thres <- 0.2

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
dds <- DESeqDataSetFromMatrix(countData = counts, colData=counts.metadata, design = ~ 1)

# combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$time, dds$strain))

design(dds) <- ~ group


###########
# Run DESeq
###########
# Additional filters for dds

# run DESeq with and without betaPrior
resultsDESeq <- DESeq(dds)

###########################################
# Analysis: WT vs A66 in strain per time point
###########################################

times <- c('0','15','40','90','180','300')

for(time in times) {
  
  group1 <- paste(time, 'WT', sep='')
  group2 <- paste(time, 'A66', sep='')
  
  # results WT vs A66 in strain per time point
  results.strain.WT_A66 <- results(resultsDESeq, contrast=c("group",group1,group2))
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





# ###########
# # Run DESeq (SAME THING WITHOUT betaPrior)
# ###########
# # Additional filters for dds
# 
# # run DESeq with and without betaPrior
# resultsDESeq.noPrior <- DESeq(dds, betaPrior=FALSE)
# 
# ###########################################
# # Analysis: WT vs A66 in strain per time point
# ###########################################
# 
# for(time in times) {
#   
#   group1 <- paste(time, 'WT', sep='')
#   group2 <- paste(time, 'A66', sep='')
#   
#   # Same for the results of DESeq without betaPrior
#   results.strain.WT_A66.noPrior <- results(resultsDESeq.noPrior, contrast=c("group",group1,group2))
#   results.signif.strain.WT_A66.noPrior <- subset(results.strain.WT_A66.noPrior, padj < padj.thres)
#   write.csv(results.strain.WT_A66.noPrior, file=paste0(filename.counts ,"_noprior_results_DESeq__strain__WT_A66_at_t", time, suffix))
#   write.csv(results.signif.strain.WT_A66.noPrior, file=paste0(filename.counts ,"_noprior_results_DESeq__strain__WT_A66_at_t", time, "__padj_", gsub("\\.","",padj.thres), suffix))
#   
#   ########################
#   # Generate Volcano plots 
#   ########################
#   # Same for the results of DESeq without betaPrior
#   png(paste0(filename.counts,"_noprior_VolcanoPlot_results_DESeq__strain__WT_A66_at_t", time, ".png"), width=2000, height=2000, pointsize=18, res=300)
#   signif.genes <- volcano_plot(results.strain.WT_A66.noPrior, title=paste0('WT vs A66 at ', time, 'm (no BetaPrior)'), thres.padj=padj.thres, thres.lfc=0.5)
#   write.csv(signif.genes, file=paste0(filename.counts ,"_noprior_results_DESeq__strain__WT_A66_at_t", time, "__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",0.5), suffix))
#   dev.off()
# 
#   ###################
#   # Generate MA-plots
#   ###################
#   # Plot mean expression vs log fold change for the contrast WT vs A66 per time point
#   png(paste0(filename.counts,"_MAplot_noprior_results_DESeq__strain__WT_A66_at_t", time, ".png"), width=2000, height=2000, pointsize=18, res=300)
#   plotMA(results.strain.WT_A66.noPrior, main=paste("Unshrunk LFC at ", time, "m (no prior)", sep=""))
#   dev.off()
# 
# }
# 

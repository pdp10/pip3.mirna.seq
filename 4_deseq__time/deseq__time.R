
###############################################
#DESeq2-based comparison of time (t1 vs t2)
###############################################




# download DESeq2 from: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library("DESeq2")

source('../utilities/plots.R')



# select the file countaining the data
location <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered"
filename.counts.metadata <- "summarised_mirna_counts_after_mapping_filtered_metadata"
suffix <-".csv"
padj.thres <- 0.05
lfc.thres <- 0.1

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
dds <- DESeqDataSetFromMatrix(countData = counts, colData=counts.metadata, design = ~ time)







###########
# Run DESeq
###########

# Additional filters for dds


# run DESeq using Likelihood Ratio Test to determine DE miRNAs over all time points
dds.res <- DESeq(dds, test='LRT', reduced = ~ 1)


############################
# Save results
############################
results.time <- results(dds.res)
results.signif.time <- subset(results.time, padj < padj.thres)
write.csv(results.time, file=paste0(filename.counts ,"_results_DESeq__time", suffix))
write.csv(results.signif.time, file=paste0(filename.counts ,"_results_DESeq__time", "__padj_", gsub("\\.","",padj.thres), suffix))

########################
# Generate Volcano plots 
########################
png(paste0(filename.counts,"_VolcanoPlot_results_DESeq__time.png"), width=2000, height=2000, pointsize=18, res=300)
signif.genes <- volcano_plot(results.time, title=paste0("DESeq over time points (LRT)"), thres.padj=padj.thres, thres.lfc=lfc.thres)
dev.off()
write.csv(signif.genes, file=paste0(filename.counts ,"_results_DESeq__time", "__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",lfc.thres), suffix))

###################
# Generate MA-plots
###################
# Plot mean expression vs log fold change
png(paste0(filename.counts,"_MAplot_results_DESeq__time.png"), width=2000, height=2000, pointsize=18, res=300)
plotMA(results.time, main="Shrunk LFC")
dev.off()

##############################
# Generate p-values histograms
##############################
png(paste0(filename.counts,"_pvalues_results_DESeq__time.png"), width=2000, height=2000, pointsize=18, res=300)
hist(results.time$pvalue, breaks=100, col="skyblue", border="slateblue", main=paste0("p-values for time"))
dev.off()









# run DESeq using Wald Test to determine DE miRNA between pairs of time points
dds.res <- DESeq(dds)



time <- c('0','15','40','90','180','300')

for (i in 1:(length(time)-1)) {
  for(j in (i+1):length(time)) {

    ############################
    # Analysis: time[i] vs time[j] in time
    ############################
    print(paste('Contrast: time', time[i], time[j]))
    
    # results time[i] vs time[j] in time
    results.time <- results(dds.res, contrast=c("time",time[i],time[j]))
    results.signif.time <- subset(results.time, padj < padj.thres)
    write.csv(results.time, file=paste0(filename.counts ,"_results_DESeq__time__", time[i], "_", time[j], suffix))
    write.csv(results.signif.time, file=paste0(filename.counts ,"_results_DESeq__time__", time[i], "_", time[j], "__padj_", gsub("\\.","",padj.thres), suffix))
    
    
    ########################
    # Generate Volcano plots 
    ########################
    # results time[i] vs time[j] in time
    png(paste0(filename.counts,"_VolcanoPlot_results_DESeq__time__", time[i], "_", time[j], ".png"), width=2000, height=2000, pointsize=18, res=300)
    signif.genes <- volcano_plot(results.time, title=paste0(time[i], " vs ", time[j]), thres.padj=padj.thres, thres.lfc=lfc.thres)
    dev.off()
    write.csv(signif.genes, file=paste0(filename.counts ,"_results_DESeq__time__", time[i], "_", time[j], "__padj_", gsub("\\.","",padj.thres), '_lfc_', gsub("\\.","",lfc.thres), suffix))

    
    ###################
    # Generate MA-plots
    ###################
    # Plot mean expression vs log fold change for the contrast time[i] vs time[j] in time
    png(paste0(filename.counts,"_MAplot_results_DESeq__time__", time[i], "_", time[j], ".png"), width=2000, height=2000, pointsize=18, res=300)
    plotMA(results.time, main="Shrunk LFC")
    dev.off()

    
    ##############################
    # Generate p-values histograms
    ##############################
    # histogram of p-values for the contrast time[i] vs time[j] in time
    png(paste0(filename.counts,"_pvalues_results_DESeq__time__", time[i], "_", time[j], ".png"), width=2000, height=2000, pointsize=18, res=300)
    hist(results.time$pvalue, breaks=100, col="skyblue", border="slateblue", main=paste0("p-values for ", time[i], " vs ", time[j]))
    dev.off()

  }
}

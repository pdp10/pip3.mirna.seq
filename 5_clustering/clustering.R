
### Calculate clustering


library(factoextra) # eigenvalues
library(ggfortify) # autoplot
library(cluster) # pam

library(grid)
library(gridExtra)

source('../utilities/plots.R')


# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered_scaled"
suffix <-".csv"


# increase the left margin to avoid cropping the label
par(mar=c(5.1, 5.1, 4.1, 2.1))



#####################################
# Read Counts Table and Samples Table
#####################################

# load counts.scaled, the count matrix scaled using scale()
counts.scaled <- read.table(paste0(location, "/", filename,suffix), sep=",", fill=T, header=T, row.names=1)




#############
# CLUSTERING
#############

clusters <- seq(2, 10)
scale=TRUE

for(k in clusters) {

  ########################
  # Run PAM clustering (colour is cluster)
  ########################
  
  # Partition around medoids
  pamx <- pam(counts.scaled, k, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
  
  # Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
  write.csv(pamx$clustering, file=paste0(filename, "_pam",k,"_labels", suffix), quote=FALSE)
  
  # Plot Clustering 
  # PC1 vs PC2
  # plot without labels

  c1c2.pam <- autoplot(pamx, data=counts.scaled, x=1, y=2, scale=scale) +
    ggtitle('PAM clustering of miRNAs') +
    #coord_fixed(ratio=0.5) +
    theme_basic()
  #ggsave(paste0(filename, "_pam",k,"_pca_c1c2.png"), width=4, height=3.5, dpi=300)
  
  # PC2 vs PC3
  # plot without labels
  c2c3.pam <- autoplot(pamx, data=counts.scaled, x=2, y=3, scale=scale) + 
    ggtitle('PAM clustering of miRNAs') +
    #coord_fixed(ratio=0.5) +
    theme_basic()
  #ggsave(paste0(filename, "_pam",k,"_pca_c2c3.png"), width=4, height=3.5, dpi=300)
  
  ## COMBINED plots
  c1c2c3.pam.combined <- arrangeGrob(c1c2.pam, c2c3.pam, ncol=2)
  ggsave(paste0(filename, "_pam",k,"_pca_c1c2c3_combined.png"), plot=c1c2c3.pam.combined, width=8, height=3.5, dpi=300)

}



#### TODO 1: 
#Add miRNA labels for the 42 miRNA which are significant using 
#- DESeq (contrast: strain),                       (6_strain)  
#- DESeq (contrast: time),                         (6_time)
#- DESeq (contrast: strain) AND PCA:PC2 (42 miRNA). (6_strain)


#### TODO 2: 
#Plot time courses using clustering information 









# ### HEATMAP using colours from clustering
# 
# 
# 
# 
# # Read the list of cells in each cluster
# clusters <- read.csv(paste0(filename, "_pam3_labels",suffix),header=TRUE, row.names=1)
# 
# ordclusters <- clusters[order(row.names(clusters)),,drop=FALSE]
# 
# 
# ##### SKIP THIS FOR NOW
# # # read the DESeq results
# # finalRes_C1_C2<-read.csv("../DESeq-ByCluster/DESeq-allEPI-Clusters-C1-C2-7Apr2017.csv",row.names=1)
# # finalRes_C1_C3<-read.csv("../DESeq-ByCluster/DESeq-allEPI-Clusters-C1-C3-7Apr2017.csv",row.names=1)
# # finalRes_C2_C3<-read.csv("../DESeq-ByCluster/DESeq-allEPI-Clusters-C2-C3-7Apr2017.csv",row.names=1)
# # 
# # # filter the DE genes based on base expression, logfoldchange and padj
# # filt_C1_C2<-subset(finalRes_C1_C2,baseMean>100 & abs(log2FoldChange) >1 & padj<1e-5)
# # filt_C1_C3<-subset(finalRes_C1_C2,baseMean>100 & abs(log2FoldChange) >1 & padj<1e-5)
# # filt_C2_C3<-subset(finalRes_C1_C2,baseMean>100 & abs(log2FoldChange) >1 & padj<1e-5)
# # 
# # # Extract the subset of genes which are diff expressed between two clusters, FOR ALL CELLS. 
# # DE_C1_C2<-subset(counts,rownames(counts) %in% rownames(filt_C1_C2))
# # DE_C1_C3<-subset(counts,rownames(counts) %in% rownames(filt_C1_C2))
# # DE_C2_C3<-subset(counts,rownames(counts) %in% rownames(filt_C1_C2))
# # 
# # # extract the data only for cells in the relevant clusters
# # DE_C1_C2_subset <- DE_C1_C2[,(which(ordclusters$x %in% c(1,2)))]
# # 
# # toplot<-data.matrix(DE_C1_C2_subset)
# 
# 
# ordcounts <- counts[order(row.names(counts)), , drop=FALSE]
# toplot <- data.matrix(ordcounts)
# 
# 
# 
# # HEATMAP of differentially expressed genes
# # Uses heatplot instead of the default heatmap. Better for gene expression data
# #################################################
# 
# # build the color bar identifying the cells from different clusters
# totcolbar <- "black"
# totcolbar[grep("1",ordclusters$x)]<-"red"  # red
# totcolbar[grep("2",ordclusters$x)]<-"black"  # black
# totcolbar[grep("3",ordclusters$x)]<-"green" # green
# 
# # extract a vector corresponding the relevant cells
# colbar <- totcolbar[(which(ordclusters$x %in% c(1,2,3)))]
# 
# par(mar = c(5, 3, 0.5, 1))
# 
# 
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("made4")
# require(made4)
# heatplot(data=toplot, ColSideColors=colbar, key=TRUE, labRow="", lhei=c(1,4))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Something is wrong. I dont get the same cluster on different days. And the bar colour is wrong.
# 
# ## get the cluster members
# #hc.rows<- hclust(dist(t(toplot)), method="average")
# #plot(hc.rows)
# #ct<- cutree(hc.rows, k=3)
# #table(ct)
# #tableclust<- data.frame(t(toplot),ct)
# 
# 
# 
# # 
# # 
# # # volcanoPlot
# # ################################################################################################
# # 
# # source('plot_volcano.r')
# # 
# # # filter the DE genes based on base expression
# # filt4v_C1_C2<-subset(finalRes_C1_C2,baseMean>100)
# # 
# # data4plot<-data.frame(dataLabels=rownames(filt4v_C1_C2),logFoldChange=filt4v_C1_C2$log2FoldChange,Pvalue=filt4v_C1_C2$padj,baseMean=filt4v_C1_C2$baseMean)
# # 
# # hitsC1_left<-finalRes_C1_C2[c("WT1","NPNT","CDH3","BNC1","PODXL"),]
# # hitsC1_right<-finalRes_C1_C2[c("CDH1"),]
# # hitsC2_left<-finalRes_C1_C2[c("FN1","THBS1"),]
# # hitsC2_right<-finalRes_C1_C2[c("THY1","TCF21","ADORA2B","CDH7"),]
# # 
# # par(mfrow=c(1,1))
# # plot.volcano(data4plot,fold.change=100,p.level=1e-100,plot.legend = F,fade_style=1,lines=F,
# #              colour2 = c(180, 0, 0),colour1 = c(0, 120, 60))
# # 
# # with(hitsC1_left, points(log2FoldChange, -log10(padj), pch=21, bg=rgb(0, 120, 60,maxColorValue=255),col="black"))
# # with(hitsC1_right, points(log2FoldChange, -log10(padj), pch=21, bg=rgb(0, 120, 60,maxColorValue=255),col="black"))
# # with(hitsC2_left, points(log2FoldChange, -log10(padj), pch=21, bg=rgb(180, 0, 0,maxColorValue=255),col="black"))
# # with(hitsC2_right, points(log2FoldChange, -log10(padj), pch=21, bg=rgb(180, 0, 0,maxColorValue=255),col="black"))
# # text(hitsC1_left$log2FoldChange,-log10(hitsC1_left$padj), labels=rownames(hitsC1_left), pos=2, cex=0.5)
# # text(hitsC1_right$log2FoldChange,-log10(hitsC1_right$padj), labels=rownames(hitsC1_right), pos=4, cex=0.5)
# # text(hitsC2_left$log2FoldChange,-log10(hitsC2_left$padj), labels=rownames(hitsC2_left), pos=2, cex=0.5)
# # text(hitsC2_right$log2FoldChange,-log10(hitsC2_right$padj), labels=rownames(hitsC2_right), pos=4, cex=0.5)
# # #identify(finalRes_C1_C2$log2FoldChange, -log10(finalRes_C1_C2$padj), labels=rownames(finalRes_C1_C2),offset=0)
# # 
# # # # Scatter plot of all genes in gray and significantly different in red. I label a few transmembrane markers
# # # #################################################
# # 
# # # average counts for all genes in C1 and C2
# # meanValues<-cbind(rowMeans(Allepi[,which(ordclusters$cluster == 1)]),rowMeans(Allepi[,which(ordclusters$cluster == 2)]))
# # colnames(meanValues)<-c("C1","C2")
# # 
# # # average counts for diff expressed genes (NB: this is clumsy. one should extract them from meanValues)
# # DEmeanValues_C1_C2<-data.matrix(meanValues[rownames(Allepi) %in% rownames(finalRes_C1_C2),])
# # 
# # smoothScatter(log(meanValues[,"C1"]),log(meanValues[,"C2"]),pch=16,col="grey",xlab="log2 (in C1)",ylab="log2 (in C2)",cex=0.5)
# # points(log(DEmeanValues_C1_C2[,"C1"]),log(DEmeanValues_C1_C2[,"C2"]),pch=16,cex=1,col="red")
# # text(log(meanValues[c("TCF21","THY1","ADORA2B","WT1","NPR1","GPC5","LRP2","CDH3","BNC1"),1]),
# #      log(meanValues[c("TCF21","THY1","ADORA2B","WT1","NPR1","GPC5","LRP2","CDH3","BNC1"),2]),
# #      label=row.names(meanValues[c("TCF21","THY1","ADORA2B","WT1","NPR1","GPC5","LRP2","CDH3","BNC1"),]),col="darkgreen",pos=3)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 

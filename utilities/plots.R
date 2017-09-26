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


library(reshape2)
library(ggplot2)
library(data.table)


#source("https://bioconductor.org/biocLite.R")
#biocLite("made4")
library(made4)




################
# ggplot2 themes
################

# A simple theme without grid
theme_basic <- function(){
  theme_bw() %+replace%
  theme(
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank()
    )
}





##############
# time courses
##############

# Plot time courses
plot_expr_tc <- function(df, name, type, line=TRUE) {
  # melt the data.table
  df.melt <- melt(df, id=c('miRNA'))
  
  g <- ggplot(df.melt, aes(x=variable, y=value, group=miRNA)) + 
    ggtitle(paste0('miRNA expression ', type)) +
    xlab('Time [m]') + 
    ylab('median standardised expression') +
    theme_basic()
  if(line) {
    g <- g + geom_line()
  } else {
    g <- g + geom_point()
  }
  ggsave(paste0(name, ".png"), width=4, height=4, dpi=300)
}

# Plot time courses based on pam class
plot_expr_tc_w_pam <- function(df, name, type, line=TRUE) {
  # melt the data.table
  df.melt <- melt(df, id=c('miRNA', 'pam'))
  
  # Plot miRNA using PAM clustering for colours
  g <- ggplot(df.melt, aes(x=variable, y=value, group=miRNA, color=pam)) + 
    ggtitle(paste0('miRNA expression ', type)) +
    xlab('Time [m]') + 
    ylab('median standardised expression') +
    theme_basic()
  if(line) {
    g <- g + geom_line()
  } else {
    g <- g + geom_point()
  }
  ggsave(paste0(name, ".png"), width=4, height=4, dpi=300)
}

# Plot median time courses based on pam class
plot_median_expr_tc_w_pam <- function(df, name, type, line=TRUE) {
  # convert data.frame to data.table
  dt <- data.table(df)  
  # median every column by `pam` class (and for each time point)
  dt.median <- dt[, lapply(.SD, median), by=pam]
  # save the data.table
  write.csv(dt.median, file=paste0(name, '.csv'), row.names=F, quote=FALSE)
  # melt the data.table
  dt.median.melt <- melt(dt.median, id=c('pam'))
  
  # Plot miRNA using PAM clustering for colours
  g <- ggplot(dt.median.melt, aes(x=variable, y=value, group=pam, color=pam)) + 
    ggtitle(paste0('miRNA expression ', type)) +
    xlab('Time [m]') + 
    ylab('median standardised expression') +
    theme_basic()
  if(line) {
    g <- g + geom_line()
  } else {
    g <- g + geom_point()
  }
  ggsave(paste0(name, ".png"), width=4, height=4, dpi=300)
}





#####
# PCA
#####

# Plot PCA nicely
plotPrettyPCA1Var <- function(x, intgroup=c("condition"), ntop=500, plotColor=TRUE) {
  x.data <- plotPCA(x, intgroup=intgroup, ntop=ntop, returnData=TRUE)
  percentVar <- round(100 * attr(x.data, "percentVar"))
  # Depending on plotColor, we decide if plotting values using different colours or shapes.
  shape <- NULL
  color <- NULL
  if(plotColor) { 
    color <- intgroup[1] 
  } else { 
    shape <- intgroup[1] 
  }
  ggplot(x.data, aes_string(x="PC1", y="PC2", color=color, shape=shape)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
}

plotPrettyPCA2Var <- function(x, intgroup=c("condition1", "condition2"), ntop=500) {
  x.data <- plotPCA(x, intgroup=intgroup, ntop=ntop, returnData=TRUE)
  percentVar <- round(100 * attr(x.data, "percentVar"))
  ggplot(x.data, aes_string(x="PC1", y="PC2", color=intgroup[1], shape=intgroup[2])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
}




##########
# HEATPLOT
##########

# use these functions for clustering like heatplot.
dist.pear <- function(x) as.dist(1-cor(t(x)))
hclust.ave <- function(x) hclust(x, method="ward")

# plot the counts matrix using heatmap.2
plot_counts_matrix_heatmap <- function(df, filename="heatmap.png", palette=redgreen(50), dendrogram="both", scale="none", trace="none", labRow=FALSE, ylab='miRNA', xlab='Samples') {
  png(file=filename, width=8, height=8, units="in", bg="white", res=300)
  heatmap.2(as.matrix(df), col=palette,
            dendrogram=dendrogram, scale=scale, trace=trace,
            distfun=dist.pear, hclustfun=hclust.ave,
            margins = c(9, 3),
            labRow = labRow, ylab=ylab, xlab=xlab)
  dev.off()
}


# plot the counts matrix using heatplot
plot_counts_matrix_heatplot <- function(df, filename="heatplot.png", dendrogram="both", scale="row", method="ward", labRow=FALSE, ylab='miRNA', xlab='Samples') {
  png(file=filename, width=8, height=8, units="in", bg="white", res=300)  
  heatplot(df, 
           dend=dendrogram, scale=scale, 
           cols.default=FALSE, lowcol="blue", highcol="yellow", 
           method=method,
           margins = c(9, 3),
           keysize=1, #key.par = list(cex=0.5)
           labRow = labRow, ylab=ylab, xlab=xlab)
  dev.off()
}




####################
# CORRELATION MATRIX
####################

# Compute and plot correlation matrix of samples
plot_corr_matrix <- function(counts, method) {
  
  # calculate the correlation matrix using pearson correlation coefficients
  counts.corr.mat <- cor(counts, method='pearson')
  
  # cut off the lower triangle
  counts.corr.mat[lower.tri(counts.corr.mat)]<- NA
  
  # Melt the correlation matrix
  counts.corr.mat.melt <- melt(counts.corr.mat, na.rm = TRUE)
  
  # Heatmap
  ggheatmap <- ggplot(data = counts.corr.mat.melt, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.975, limit = c(0.95,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    scale_y_discrete(position="right") + 
    #theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.3, 0.4))+
    coord_fixed()
  return(ggheatmap)
}





############
# STATISTICS
############

plot_read_density <- function(d) {
  d$samples <- paste(d$strain, d$time, d$rep, sep = "_")
  p <- ggplot(as.data.frame(d), aes(value, color = samples)) +
    geom_density() +
    theme_basic() + 
    theme(legend.key.size = unit(0.4, "cm")) +
    scale_x_log10()
  ggsave("miRNA_read_density_per_sample.png", width=6, height=4, dpi=300)
  p <- p + facet_grid(strain ~ time)
  ggsave("miRNA_read_density_per_sample_facet.png", width=6, height=4, dpi=300)    
}

plot_read_ecdf <- function(d) {
  d$samples <- paste(d$strain, d$time, d$rep, sep = "_")
  p <- ggplot(as.data.frame(d), aes(value, color = samples)) +
    stat_ecdf() +
    theme_basic() + 
    theme(legend.key.size = unit(0.4, "cm")) +
    scale_x_log10() 
  ggsave("miRNA_read_ecdf_per_sample.png", width=6, height=4, dpi=300)
  p <- p + facet_grid(strain ~ time)
  ggsave("miRNA_read_ecdf_facet_per_sample.png", width=6, height=4, dpi=300)  
}

plot_counts_density <- function(df, filename) {
  png(file=filename, width=8, height=8, units="in", bg="white", res=300)
  plot(density(as.matrix(df)))
  dev.off()
}


###############
# Volcano plots
###############

# Generate a volcano plot and return the list of significant genes
volcano_plot <- function(deseq2.res, thres.padj = 0.05, thres.lfc = 0.5, title="padj versus fold change", show.x.annot=TRUE, show.y.annot=TRUE) {
  df <- data.frame(log2FoldChange=deseq2.res$log2FoldChange, negLog10Padj=-log10(deseq2.res$padj))
  rownames(df) <- rownames(deseq2.res)
  # plot
  par(mar = c(5,5,5,4))
  signif.genes <- ( df$negLog10Padj > -log10(thres.padj) & abs(df$log2FoldChange) > thres.lfc )
  plot(df, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~padj),
       pch=16, cex=0.5)
  
  points(df, pch=16, cex=0.5, col="grey")
  points(df[signif.genes & df$log2FoldChange>thres.lfc, ], pch=16, cex=0.5, col="red") 
  points(df[signif.genes & df$log2FoldChange<thres.lfc, ], pch=16, cex=0.5, col="blue") 
  
  abline(h = -log10(thres.padj), col="black", lty=2, lwd=1.5)
  abline(v = c(-thres.lfc, thres.lfc), col="black", lty=2, lwd=1.5)
  if(show.y.annot) {
    mtext(paste0("padj=",thres.padj), side=4, at=-log10(thres.padj), cex=0.8, line=0.5, las=1)
  }
  if(show.x.annot) {
    mtext(c(paste0("-", thres.lfc), paste0("+", thres.lfc)), side=3, at=c(-thres.lfc-0.05,thres.lfc+0.05), cex=0.8, line=0.2)
  }
  title(title)
  return(df[signif.genes,])
}





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





### Calculate PCA and T-SNE for the filtered counts data set


library(factoextra) # eigenvalues
library(rgl) # 3d plot
library(Rtsne)

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



########################
# Prepare data sets
########################


# Apply a log transformation for the counts, and transpose this for PCA (we cluster samples, not miRNA here.)
counts.scaled.t <- t(counts.scaled)

# Now, create a data set from counts.scaled.t with `strain`, `time`, `egf`, and `replicate` 
counts.scaled.t.metadata <- counts.scaled.t


## Extract STRAIN
# use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(counts.scaled)), 'WT', 'A66')

## Extract TIME
time <- gsub(".*[_]([^.]+)[_].*", "\\1", colnames(counts.scaled))
# now replace 'noEGF' with '300', as noEGF sample was taken at 300m. 
time <- gsub('noEGF', '300', time)

## Extract EGF
# use grepl to create a logical vector whether noEGF appears or not in n.
egf <- ifelse(grepl('noEGF', colnames(counts.scaled)), 'noEGF', 'EGF')

## Extract REPLICATE
replicate <- gsub(".*[_]([^.]+)", "\\1", colnames(counts.scaled))
replicate <- paste0(replicate)

# now we add the columns
counts.scaled.t.metadata <- cbind(counts.scaled.t.metadata, strain)
counts.scaled.t.metadata <- cbind(counts.scaled.t.metadata, time)
counts.scaled.t.metadata <- cbind(counts.scaled.t.metadata, egf)
counts.scaled.t.metadata <- cbind(counts.scaled.t.metadata, replicate)




##########
# Run PCA
##########

# Compute PCA

# `scale` and `centre` must be FALSE because we already scaled and centred our data.
scale <- FALSE
centre <- FALSE
pca <- prcomp(counts.scaled.t, center=centre, scale.=scale)

# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca)     

# Plot the variance of PCA components
png(paste0(filename,"_pca_comp_variances",".png"), width=2000, height=1500, res=300)
plot(pca, type = "l", main='Variance of PCA components')
dev.off()

# Show PCA summary
summary(pca)

# Write the calculated PCA `rotation` (PCA load). This is the table miRNAs (rows) vs PCAs (cols).
write.csv(pca$rotation, file=paste0(filename,"_PCA_rotation", suffix), quote=FALSE)


#############################
# Plot PCA (colour is strain)
#############################

# PC1 vs PC2
# plot without labels
c1c2.strain <- autoplot(pca, data=counts.scaled.t.metadata, colour='strain', x=1, y=2, scale=scale) +
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_strain_c1c2.png"), width=4, height=4, dpi=300)


# PC2 vs PC3
# plot without labels
c2c3.strain <- autoplot(pca, data=counts.scaled.t.metadata, colour='strain', x=2, y=3, scale=scale) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_strain_c2c3.png"), width=4, height=4, dpi=300)


#############################
# Plot PCA (shape is time)
#############################

# PC1 vs PC2
# plot without labels
c1c2.time <- autoplot(pca, data=counts.scaled.t.metadata, shape='time', x=1, y=2, scale=scale) +
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_time_c1c2.png"), width=4, height=4, dpi=300)


# PC2 vs PC3
# plot without labels
c2c3.time <- autoplot(pca, data=counts.scaled.t.metadata, shape='time', x=2, y=3, scale=scale) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_time_c2c3.png"), width=4, height=4, dpi=300)


#############################
# Plot PCA (colour is strain, shape is time)
#############################

# PC1 vs PC2
# plot without labels
c1c2.strain.time <- autoplot(pca, data=counts.scaled.t.metadata, colour='strain', shape='time', x=1, y=2, scale=scale) +
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_strain_time_c1c2.png"), width=4, height=4, dpi=300)


# PC2 vs PC3
# plot without labels
c2c3.strain.time <- autoplot(pca, data=counts.scaled.t.metadata, colour='strain', shape='time', x=2, y=3, scale=scale) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()
#ggsave(paste0(filename, "_pca_strain_time_c2c3.png"), width=4, height=4, dpi=300)



## COMBINED plots
c1c2.combined <- arrangeGrob(c1c2.strain, c1c2.time, c1c2.strain.time, ncol=2)
ggsave(paste0(filename, "_pca_c1c2_combined.png"), plot=c1c2.combined, width=7, height=7, dpi=300)

c2c3.combined <- arrangeGrob(c2c3.strain, c2c3.time, c2c3.strain.time, ncol=2)
ggsave(paste0(filename, "_pca_c2c3_combined.png"), plot=c2c3.combined, width=7, height=7, dpi=300)







#############################################################

## EMULATE ggplot2 colours 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colours = gg_color_hue(n=2)


###################################
# plot the third dimension (3D PCA)
###################################
pc3d <- cbind(pca$x[,1], pca$x[,2], pca$x[,3])
plot3d(pc3d, col=colours, xlab="PC1", ylab="PC2",zlab="PC3", main="PCA of samples: pc1 vs pc2 vs pc3", type="s",size=2,scale=0.4)
rgl.snapshot(filename=paste0(filename,"_3d_pca_c1c2c3.png"),"png")




####################
# t-SNE
####################

rtsnecb_out <- Rtsne(as.matrix(counts.scaled.t), verbose = TRUE, dims=2, perplexity=10, max_iter=1000)

plot(rtsnecb_out$Y,col=colours,pch=16, cex=1.5,cex.lab=1.5,cex.axis=1.5,main="t-SNE of samples: dim1 vs dim2",xlab="dim1",ylab="dim2")
dev.copy(png,paste0(filename,"_tSNE_dim1_dim2.png"))
dev.off()



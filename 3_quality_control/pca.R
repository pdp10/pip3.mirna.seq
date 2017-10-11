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

library(ggplot2)

library(grid)
library(gridExtra)

source('../utilities/plots.R')



#########
# MY DATA
#########

suffix <-".csv"
scale <- TRUE
centre <- TRUE

# select the file countaining the data
location <- "../data"
# the data must be log-transformed first. We also scale it.
filename <- "summarised_mirna_counts_after_mapping_filtered_rlog_scaled"
filename.mean <- "summarised_mirna_counts_after_mapping_filtered_rlog_mean_scaled"

# increase the left margin to avoid cropping the label
#par(mar=c(5.1, 5.1, 4.1, 2.1))



#####################################
# Read Counts Table and Samples Table
#####################################

# load counts
counts <- read.table(paste0(location, "/", filename,suffix), sep=",", fill=T, header=T, row.names=1)

counts.mean <- read.table(paste0(location, "/", filename.mean,suffix), sep=",", fill=T, header=T, row.names=1)




###################
# Prepare data sets (used for ggplot2 colours, shape)
###################

## Extract STRAIN
# use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(counts)), 'WT', ifelse(grepl('noEGF', colnames(counts)), 'A66 no EGF', 'A66'))

## Extract TIME
time <- gsub(".*[_]([^.]+)[_].*", "\\1", colnames(counts))
# now replace 'noEGF' with '300', as noEGF sample was taken at 300m. 
time <- gsub('noEGF', '300', time)

## Extract EGF
# use grepl to create a logical vector whether noEGF appears or not in n.
egf <- ifelse(grepl('noEGF', colnames(counts)), 'noEGF', 'EGF')

## Extract REPLICATE
replicate <- gsub(".*[_]([^.]+)", "\\1", colnames(counts))
replicate <- paste0(replicate)



####### SAME FOR counts.mean  ########## 

## Extract STRAIN
# use grepl to create a logical vector whether WT appears or not in n.
strain.mean <- ifelse(grepl('WT', colnames(counts.mean)), 'WT', ifelse(grepl('noEGF', colnames(counts.mean)), 'A66 no EGF', 'A66'))

## Extract TIME
time.mean <- gsub(".*[_]([^.]+).*", "\\1", colnames(counts.mean))
# now replace 'noEGF' with '300', as noEGF sample was taken at 300m. 
time.mean <- gsub('noEGF', '300', time.mean)

## Extract EGF
# use grepl to create a logical vector whether noEGF appears or not in n.
egf.mean <- ifelse(grepl('noEGF', colnames(counts.mean)), 'noEGF', 'EGF')




########################
# Run PCA of the samples
########################

# we need to traspose the counts
pca.samples <- prcomp(t(counts), center=centre, scale.=scale)
# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca.samples)     
# Plot the variance of PCA components
png(paste0(filename,"_samples_pca_comp_variances",".png"), width=2000, height=1500, res=300)
plot(pca.samples, type = "l", main='Variance of PCA components')
dev.off()
# Write the calculated PCA `rotation` (PCA load). This is the table miRNAs (rows) vs PCAs (cols).
write.csv(pca.samples$rotation, file=paste0(filename,"_samples_PCA_rotation", suffix), quote=FALSE)




###############################
# Run and plot PCA of the reads
###############################

pca.reads <- prcomp(counts, center=centre, scale.=scale)
# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca.reads)     
# Plot the variance of PCA components
png(paste0(filename,"_reads_pca_comp_variances",".png"), width=2000, height=1500, res=300)
plot(pca.reads, type = "l", main='Variance of PCA components')
dev.off()
# Write the calculated PCA `rotation` (PCA load). This is the table samples (rows) vs PCAs (cols).
write.csv(pca.reads$rotation, file=paste0(filename,"_reads_PCA_rotation", suffix), quote=FALSE)

# prepare the data set for plotting
df <- data.frame(pca.reads$rotation, 
                 strain=strain, time=as.numeric(time), egf=egf, replicate=as.numeric(replicate), 
                 check.names = FALSE)
# Plot PCA (colour is strain, shape is time)
plot_pca(df, eigen, filename=paste0(filename, "_pca_c1c2c3_combined.png"))




###############################
# Run and plot PCA of the reads (mean of the replicates)
###############################

pca.reads.mean <- prcomp(counts.mean, center=centre, scale.=scale)
eigen <- get_eig(pca.reads.mean)     
png(paste0(filename.mean,"_reads_mean_pca_comp_variances",".png"), width=2000, height=1500, res=300)
plot(pca.reads.mean, type = "l", main='Variance of PCA components')
dev.off()
write.csv(pca.reads.mean$rotation, file=paste0(filename.mean,"_reads_mean_PCA_rotation", suffix), quote=FALSE)

df.mean <- data.frame(pca.reads.mean$rotation, 
                 strain=strain.mean, time=as.numeric(time.mean), egf=egf.mean, 
                 check.names = FALSE)
# Plot PCA (colour is strain, shape is time)
plot_pca(df.mean, eigen, filename=paste0(filename.mean, "_pca_c1c2c3_combined.png"))














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
pc3d <- cbind(pca.samples$x[,1], pca.samples$x[,2], pca.samples$x[,3])
plot3d(pc3d, col=colours, xlab="PC1", ylab="PC2",zlab="PC3", main="PCA on samples", type="s",size=2,scale=0.4)
rgl.snapshot(filename=paste0(filename,"_3d_pca_c1c2c3.png"),"png")




####################
# t-SNE
####################

rtsnecb_out <- Rtsne(as.matrix(t(counts)), verbose = TRUE, dims=2, perplexity=10, max_iter=1000)

plot(rtsnecb_out$Y,col=colours,pch=16, cex=1.5,cex.lab=1.5,cex.axis=1.5,main="t-SNE of samples: dim1 vs dim2",xlab="dim1",ylab="dim2")
dev.copy(png,paste0(filename,"_tSNE_dim1_dim2.png"))
dev.off()



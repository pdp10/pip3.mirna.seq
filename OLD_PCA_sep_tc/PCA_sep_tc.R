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



library(factoextra)
library(rgl)
library(Rtsne)

# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping"
suffix <-"csv"

# build filename containing the counts
counts.filename <- paste(location,"/",filename,".",suffix,sep="")

# load counts
counts <- read.table(counts.filename, sep=",",fill=T,header=T,row.names=1)


########################
#### counts filtering
########################

# remove lowly expressed genes (with less than 10 reads in every samples)
maxExp <- do.call(pmax, counts)
Thresh <- counts[maxExp>1,]

# remove genes which variation across samples is less than two-fold
maxExpT <- do.call(pmax, Thresh)
minExpT <- do.call(pmin, Thresh)
Filt <- Thresh[(maxExpT+1)/(minExpT+1)>2,]

########################
#### PCA for all samples
########################

# transpose data for PCA
tdata<-t(log(Filt+0.1))

# # PCA
# PCA<-prcomp(tdata)
# eigen<-get_eig(PCA)     # This provide a list of components with the respective variance
# 
# write.csv(PCA$rotation, file=paste("PCA_",filename,"_rot.csv",sep=""), quote=FALSE)
# 
# # increase the left margin to avoid cropping the label
# par(mar=c(5.1, 5.1, 4.1, 2.1))
# 
# # build color vector
# colours<-rep(1,length(colnames(counts)))     # default colour black
# colours[grep("miWT",colnames(counts))]<-"blue"    
# colours[grep("miA66",colnames(counts))]<-"magenta"
# colours[grep("unmap",colnames(counts))]<-"orange"
# colours[grep("base",colnames(counts))]<-"green"
# # print(colours)
# 
# # build symbol vector
# symbols<-rep(3, length(colnames(counts)))       # default symbol crosses
# symbols[grep("miWT",colnames(counts))]<-15    # square
# symbols[grep("miA66",colnames(counts))]<-16     # circle
# symbols[grep("unmap",colnames(counts))]<-17     # triangle
# symbols[grep("base",colnames(counts))]<-18     # rhombus
# 
# # aggregate plots 
# #par(mfrow=c(2,2))
# 
# plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
#      xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
# text(PCA$x[,1], PCA$x[,2], labels=colnames(counts), pos=1, cex=0.7)
# dev.copy(png, paste("PCA_",filename,"_ALL_pc1_pc2_w_names",".png",sep=""))
# dev.off()
# plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
#      xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
# dev.copy(png,(paste("PCA_",filename,"_ALL_pc1_pc2",".png",sep="")))
# dev.off()
# 
# plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
#      xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
# text(PCA$x[,2], PCA$x[,3], labels=colnames(counts), pos=1, cex=0.7)
# dev.copy(png,paste("PCA_",filename,"_ALL_pc2_pc3_w_names",".png",sep=""))
# dev.off()
# plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
#      xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
# #text(PCA$x[,2], PCA$x[,3], labels=colnames(counts), pos=1, cex=0.7)
# dev.copy(png,paste("PCA_",filename,"_ALL_pc2_pc3",".png",sep=""))
# dev.off()
# 
# # plot the third dimension
# pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
# plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
# dev.copy(png,paste("PCA_",filename,"_ALL_pc1_pc2_pc3",".png",sep=""))
# dev.off()


########################
#### PCA without unmapped
########################

#tdata.No.unmapped<-tdata[rownames(tdata) != "unmap", ]
#counts.No.unmapped<-counts[rownames(counts) != "unmap", ]

tdata.No.unmapped<-tdata[-grep('unmap', rownames(tdata)), ]
counts.No.unmapped<-counts[, -grep('unmap', colnames(counts))]



# PCA
PCA<-prcomp(tdata.No.unmapped)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoUnmapped_sepTC_rot.csv",sep=""), quote=FALSE)

# build color vector
colours<-rep(1,length(colnames(counts.No.unmapped)))     # default colour black
colours[grep("miWT",colnames(counts.No.unmapped))]<-"blue"    
colours[grep("miA66",colnames(counts.No.unmapped))]<-"magenta"
colours[grep("unmap",colnames(counts.No.unmapped))]<-"orange"
colours[grep("base",colnames(counts.No.unmapped))]<-"green"
# print(colours)

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.unmapped)))       # default symbol crosses
symbols[grep("_0_",colnames(counts.No.unmapped))]<-0
symbols[grep("_15_",colnames(counts.No.unmapped))]<-1
symbols[grep("_40_",colnames(counts.No.unmapped))]<-2
symbols[grep("_90_",colnames(counts.No.unmapped))]<-3
symbols[grep("_180_",colnames(counts.No.unmapped))]<-4
symbols[grep("_300_",colnames(counts.No.unmapped))]<-5
symbols[grep("_noEGF_",colnames(counts.No.unmapped))]<-15
symbols[grep("unmap",colnames(counts.No.unmapped))]<-20
symbols[grep("base",colnames(counts.No.unmapped))]<-20

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples by tc w/o unmapped",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.unmapped), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoUnmapped_sepTC_pc1_pc2_w_text",".png",sep=""))
dev.off()
plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples by tc w/o unmapped",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
#text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.unmapped), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoUnmapped_sepTC_pc1_pc2",".png",sep=""))
dev.off()

plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples by tc w/o unmapped",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.unmapped), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoUnmapped_sepTC_pc2_pc3_w_text",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples by tc w/o unmapped",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
#text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.unmapped), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoUnmapped_sepTC_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoUnmapped_sepTC_pc1_pc2_pc3",".png",sep=""))
dev.off()


########################
#### t-SNE without unmap
########################

rtsnecb_out <- Rtsne(as.matrix(tdata.No.unmapped), verbose = TRUE, dims=2, perplexity=10, max_iter=1000)

plot(rtsnecb_out$Y,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="t-SNE by tc w/o unmapped",xlab="dim1",ylab="dim2")
dev.copy(png,paste("tSNE_",filename,"_NoUnmapped_sepTC_dim1_dim2",".png",sep=""))
dev.off()

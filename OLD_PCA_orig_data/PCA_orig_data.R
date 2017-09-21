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
filename <- "summarised_mirna_counts"
suffix <-"csv"
typecounts <- "batch corrected"      # needed for legend of figure

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

# PCA
PCA<-prcomp(tdata)
eigen<-get_eig(PCA)     # This provide a list of components with the respective variance

write.csv(PCA$rotation, file=paste("PCA_",filename,"_rot.csv",sep=""), quote=FALSE)

# increase the left margin to avoid cropping the label
par(mar=c(5.1, 5.1, 4.1, 2.1))

# build color vector
colours<-rep(1,length(colnames(counts)))    # default colour black
colours[grep("lane1",colnames(counts))]<-4  # blue
colours[grep("lane2",colnames(counts))]<-6  # magenta
colours[grep("lane3",colnames(counts))]<-3  # green
# print(colours)

# build symbol vector


plot(PCA$x,col=colours,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_pc1_pc2",".png",sep=""))
dev.off()

plot(PCA$x[,2],PCA$x[,3],pch=16,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA all samples",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_pc1_pc2_pc3",".png",sep=""))
dev.off()


########################
#### PCA without GTATGC
########################

#tdata.No.GTATGC<-tdata[rownames(tdata) != "GTATGC", ]
#counts.No.GTATGC<-counts[rownames(counts) != "GTATGC", ]

tdata.No.GTATGC<-tdata[-grep('GTATGC', rownames(tdata)), ]
counts.No.GTATGC<-counts[, -grep('GTATGC', colnames(counts))]



# PCA
PCA<-prcomp(tdata.No.GTATGC)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoGTATGC_rot.csv",sep=""), quote=FALSE)

# build color vector 
colours<-rep(1,length(colnames(counts.No.GTATGC)))    # default colour black
colours[grep("lane1",colnames(counts.No.GTATGC))]<-4  # blue
colours[grep("lane2",colnames(counts.No.GTATGC))]<-6  # magenta
colours[grep("lane3",colnames(counts.No.GTATGC))]<-3  # green

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.GTATGC)))       # default symbol crosses
symbols[grep("lane2",colnames(counts.No.GTATGC))]<-16    # disks
symbols[grep("lane3",colnames(counts.No.GTATGC))]<-6     # triangles

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.GTATGC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_pc1_pc2",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.GTATGC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_pc1_pc2_pc3",".png",sep=""))
dev.off()



########################
#### PCA without GTCGAG
########################

#tdata.No.GTCGAG<-tdata[rownames(tdata) != "GTCGAG", ]
#counts.No.GTCGAG<-counts[rownames(counts) != "GTCGAG", ]

tdata.No.GTCGAG<-tdata[-grep('GTCGAG', rownames(tdata)), ]
counts.No.GTCGAG<-counts[, -grep('GTCGAG', colnames(counts))]



# PCA
PCA<-prcomp(tdata.No.GTCGAG)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoGTCGAG_rot.csv",sep=""), quote=FALSE)

# build color vector 
colours<-rep(1,length(colnames(counts.No.GTCGAG)))    # default colour black
colours[grep("lane1",colnames(counts.No.GTCGAG))]<-4  # blue
colours[grep("lane2",colnames(counts.No.GTCGAG))]<-6  # magenta
colours[grep("lane3",colnames(counts.No.GTCGAG))]<-3  # green

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.GTCGAG)))       # default symbol crosses
symbols[grep("lane2",colnames(counts.No.GTCGAG))]<-16    # disks
symbols[grep("lane3",colnames(counts.No.GTCGAG))]<-6     # triangles

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTCGAG",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.GTCGAG), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTCGAG_pc1_pc2",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTCGAG",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.GTCGAG), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTCGAG_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoGTCGAG_pc1_pc2_pc3",".png",sep=""))
dev.off()



########################
#### PCA without GTGGCC
########################

#tdata.No.GTGGCC<-tdata[rownames(tdata) != "GTGGCC", ]
#counts.No.GTGGCC<-counts[rownames(counts) != "GTGGCC", ]

tdata.No.GTGGCC<-tdata[-grep('GTGGCC', rownames(tdata)), ]
counts.No.GTGGCC<-counts[, -grep('GTGGCC', colnames(counts))]



# PCA
PCA<-prcomp(tdata.No.GTGGCC)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoGTGGCC_rot.csv",sep=""), quote=FALSE)

# build color vector 
colours<-rep(1,length(colnames(counts.No.GTGGCC)))    # default colour black
colours[grep("lane1",colnames(counts.No.GTGGCC))]<-4  # blue
colours[grep("lane2",colnames(counts.No.GTGGCC))]<-6  # magenta
colours[grep("lane3",colnames(counts.No.GTGGCC))]<-3  # green

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.GTGGCC)))       # default symbol crosses
symbols[grep("lane2",colnames(counts.No.GTGGCC))]<-16    # disks
symbols[grep("lane3",colnames(counts.No.GTGGCC))]<-6     # triangles

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTGGCC",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.GTGGCC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTGGCC_pc1_pc2",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTGGCC",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.GTGGCC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTGGCC_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoGTGGCC_pc1_pc2_pc3",".png",sep=""))
dev.off()




#####################################
#### PCA without GTATGC,GTCGAG
#####################################

#print(colnames(counts))

tdata.No.GTATGC.GTCGAG<-tdata.No.GTATGC[-grep('GTCGAG', rownames(tdata.No.GTATGC)), ]
counts.No.GTATGC.GTCGAG<-counts.No.GTATGC[,-grep('GTCGAG', colnames(counts.No.GTATGC))]
#print(colnames(counts.No.GTATGC.GTCGAG))

# PCA
PCA<-prcomp(tdata.No.GTATGC.GTCGAG)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_rot.csv",sep=""), quote=FALSE)

# build color vector 
colours<-rep(1,length(colnames(counts.No.GTATGC.GTCGAG)))    # default colour black
colours[grep("lane1",colnames(counts.No.GTATGC.GTCGAG))]<-4  # blue
colours[grep("lane2",colnames(counts.No.GTATGC.GTCGAG))]<-6  # magenta
colours[grep("lane3",colnames(counts.No.GTATGC.GTCGAG))]<-3  # green

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.GTATGC.GTCGAG)))       # default symbol crosses
symbols[grep("lane2",colnames(counts.No.GTATGC.GTCGAG))]<-16    # disks
symbols[grep("lane3",colnames(counts.No.GTATGC.GTCGAG))]<-6     # triangles

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC,GTCGAG",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.GTATGC.GTCGAG), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_pc1_pc2",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC,GTCGAG",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.GTATGC.GTCGAG), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_pc1_pc2_pc3",".png",sep=""))
dev.off()

########################
#### t-SNE without GTATGC,GTCGAG
########################

rtsnecb_out <- Rtsne(as.matrix(tdata.No.GTATGC.GTCGAG), verbose = TRUE, dims=2, perplexity=10, max_iter=1000)

plot(rtsnecb_out$Y,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="t-SNE w/o GTATGC,GTCGAG",xlab="dim1",ylab="dim2")
dev.copy(png,paste("tSNE_",filename,"_NoGTATGC_NoGTCGAG_dim1_dim2",".png",sep=""))
dev.off()







#####################################
#### PCA without GTATGC,GTCGAG,GTGGCC
#####################################

#print(colnames(counts))
#print(colnames(counts.No.GTATGC))

tdata.No.GTATGC.GTCGAG.GTGGCC<-tdata.No.GTATGC[-grep('GTCGAG', rownames(tdata.No.GTATGC)), ]
counts.No.GTATGC.GTCGAG.GTGGCC<-counts.No.GTATGC[,-grep('GTCGAG', colnames(counts.No.GTATGC))]
#print(colnames(counts.No.GTATGC.GTCGAG.GTGGCC))

tdata.No.GTATGC.GTCGAG.GTGGCC<-tdata.No.GTATGC.GTCGAG.GTGGCC[-grep('GTGGCC', rownames(tdata.No.GTATGC.GTCGAG.GTGGCC)), ]
counts.No.GTATGC.GTCGAG.GTGGCC<-counts.No.GTATGC.GTCGAG.GTGGCC[,-grep('GTGGCC', colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]
#print(colnames(counts.No.GTATGC.GTCGAG.GTGGCC))




# PCA
PCA<-prcomp(tdata.No.GTATGC.GTCGAG.GTGGCC)
eigen<-get_eig(PCA)      # library(factoextra)

write.csv(PCA$rotation, file=paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_NoGTGGCC_rot.csv",sep=""), quote=FALSE)

# build color vector 
colours<-rep(1,length(colnames(counts.No.GTATGC.GTCGAG.GTGGCC)))    # default colour black
colours[grep("lane1",colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]<-4  # blue
colours[grep("lane2",colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]<-6  # magenta
colours[grep("lane3",colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]<-3  # green

# build symbol vector
symbols<-rep(3, length(colnames(counts.No.GTATGC.GTCGAG.GTGGCC)))       # default symbol crosses
symbols[grep("lane2",colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]<-16    # disks
symbols[grep("lane3",colnames(counts.No.GTATGC.GTCGAG.GTGGCC))]<-6     # triangles

plot(PCA$x,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC,GTCGAG,GTGGCC",
     xlab=sprintf("PC1 %.1f %%",eigen[1,2]),ylab=sprintf("PC2 %.1f %%",eigen[2,2]))
text(PCA$x[,1], PCA$x[,2], labels=colnames(counts.No.GTATGC.GTCGAG.GTGGCC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_NoGTGGCC_pc1_pc2",".png",sep=""))
dev.off()
plot(PCA$x[,2],PCA$x[,3],pch=symbols,col=colours,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="PCA samples w/o GTATGC,GTCGAG,GTGGCC",
     xlab=sprintf("PC2 %.1f %%",eigen[2,2]),ylab=sprintf("PC3 %.1f %%",eigen[3,2]))
text(PCA$x[,2], PCA$x[,3], labels=colnames(counts.No.GTATGC.GTCGAG.GTGGCC), pos=1, cex=0.7)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_NoGTGGCC_pc2_pc3",".png",sep=""))
dev.off()

# plot the third dimension
pc3d<-cbind(PCA$x[,1], PCA$x[,2], PCA$x[,3])
plot3d(pc3d, col=colours,xlab="PC1", ylab="PC2",zlab="PC3",type="s",size=1,scale=0.2)
dev.copy(png,paste("PCA_",filename,"_NoGTATGC_NoGTCGAG_NoGTGGCC_pc1_pc2_pc3",".png",sep=""))
dev.off()




########################
#### t-SNE without GTATGC,GTCGAG,GTGGCC
########################

rtsnecb_out <- Rtsne(as.matrix(tdata.No.GTATGC.GTCGAG.GTGGCC), verbose = TRUE, dims=2, perplexity=10, max_iter=1000)

plot(rtsnecb_out$Y,col=colours,pch=symbols,cex=1.5,cex.lab=1.5,cex.axis=1.5,main="t-SNE w/o GTATGC,GTCGAG,GTGGCC",xlab="dim1",ylab="dim2")
dev.copy(png,paste("tSNE_",filename,"_NoGTATGC_NoGTCGAG_NoGTGGCC_dim1_dim2",".png",sep=""))
dev.off()

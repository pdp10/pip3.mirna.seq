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



# load necessary libraries
library(Rtsne)

# select the file countaining the data
location <- "../data"
prefix <- "Allepi-CPM-clean-filt-log-batCor"
suffix <-"csv"

# genes of interest TODO: WRITE LOOPS!
GeneName1<-"TCF21"
GeneName2<-"BNC1"
GeneName3<-"WT1"

# build filename containing the counts
countfile <- paste(location,"/",prefix,".",suffix,sep="")

# load counts
counts <- read.table(countfile, sep=",",fill=T,header=T,row.names=1)

# transpose the count matrix for the tSNE
tcounts <- t(counts)

########################
# tSNE
#######################################################################

rtsne_out <- Rtsne(as.matrix(tcounts), verbose = TRUE, dims=3, perplexity=30, max_iter=2000)

# Get the coordinates from the t-SNE
coordinates <- rtsne_out$Y

# Attach the name of the cells to the coordinates
rownames(coordinates) <- rownames(tcounts)
write.csv(coordinates, file=paste("tSNE-coordinates-",prefix,".csv",sep=""),quote=FALSE)

# Get the maximum expression of interesting genes for normalisation
# NB: data is log. So expr(X)-expr(Y) is in fact log(x)-log(y) = log(x/y), 
# and therefore represent logfoldchange

maxGene1<-max(tcounts[,GeneName1])
maxGene2<-max(tcounts[,GeneName2])
maxGene3<-max(tcounts[,GeneName3])

# color with one gene

goi<-tcounts[,GeneName1]
rbPal <- colorRampPalette(c('blue',"white",'red'))
datcol <- rbPal(10)[as.numeric(cut(goi,breaks = 10))]
plot(rtsne_out$Y,col=datcol,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main=paste("tSNE", GeneName1," (red = high expression)"),
     xlab="dim1",ylab="dim2")

dev.copy(png,paste("tSNE-",prefix,"-",GeneName1,".png",sep=""))
dev.off()

goi<-tcounts[,GeneName2]
rbPal <- colorRampPalette(c('blue',"white",'red'))
datcol <- rbPal(10)[as.numeric(cut(goi,breaks = 10))]
plot(rtsne_out$Y,col=datcol,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main=paste("tSNE", GeneName2," (red = high expression)"),
     xlab="dim1",ylab="dim2")

dev.copy(png,paste("tSNE-",prefix,"-",GeneName2,".png",sep=""))
dev.off()

goi<-tcounts[,GeneName3]
rbPal <- colorRampPalette(c('blue',"white",'red'))
datcol <- rbPal(10)[as.numeric(cut(goi,breaks = 10))]
plot(rtsne_out$Y,col=datcol,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main=paste("tSNE", GeneName3," (red = high expression)"),
     xlab="dim1",ylab="dim2")

dev.copy(png,paste("tSNE-",prefix,"-",GeneName3,".png",sep=""))
dev.off()

# Colour with the expression of two genes

goi<-tcounts[,GeneName2]/maxGene2-tcounts[,GeneName1]/maxGene1

rbPal <- colorRampPalette(c('darkred',"yellow",'darkcyan'))
datcol <- rbPal(10)[as.numeric(cut(goi,breaks = 10))]
plot(rtsne_out$Y,col=datcol,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main=paste("tSNE Colours = ", GeneName2," (green) , ",GeneName1," (red)"),xlab="dim1",ylab="dim2")

dev.copy(png,paste("tSNE-",prefix,"-",GeneName2,"-",GeneName1,".png",sep=""))
dev.off()


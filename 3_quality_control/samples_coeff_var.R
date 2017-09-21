# Calculate the coefficient of variance vs mean for treatment / time point.


# For each miRNA/condition, get the means and SD between replicates. Then plot X=mean, Y=CV
# You should get plenty of points, that should follow a classical square decay.
# Most of the variation is at very low mean values because noise decreases as the sqrt(signal).
# Two plots:
# - scatterplot with colors representing Ctrl, EGF-treated and EGF+A66
# - scatterplot with colors representing the times.




# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered"
suffix <-".csv"


# load counts
counts <- read.table(paste(location,"/",filename,suffix,sep=""), sep=",",fill=TRUE,header=TRUE,row.names=1)


# increase the left margin to avoid cropping the label
par(mar=c(5.1, 5.1, 4.1, 2.1))


############################################
# For each miRNA/condition, get the means and SD between replicates. Then plot X=mean, Y=CV
# - scatterplot with colors representing Ctrl, EGF-treated and EGF+A66
############################################

# separate groups
# CTRL
counts.ctrl <- counts[, grep('noEGF', colnames(counts))]
counts.ctrl <- unique(counts.ctrl)
counts.ctrl.mean <- apply(counts.ctrl, 1, mean)
counts.ctrl.sd <- apply(counts.ctrl, 1, sd)
counts.ctrl.cv <- counts.ctrl.sd/counts.ctrl.mean
# EGF-trtm
counts.egf.trtm <- counts[, grep('miWT', colnames(counts))]
counts.egf.trtm.mean <- apply(counts.egf.trtm, 1, mean)
counts.egf.trtm.sd <- apply(counts.egf.trtm, 1, sd)
counts.egf.trtm.cv <- counts.egf.trtm.sd/counts.egf.trtm.mean
# EGF-a66
counts.egf.a66 <- counts[, grep("miA66_[[:digit:]]", colnames(counts))]
counts.egf.a66.mean <- apply(counts.egf.a66, 1, mean)
counts.egf.a66.sd <- apply(counts.egf.a66, 1, sd)
counts.egf.a66.cv <- counts.egf.a66.sd/counts.egf.a66.mean

symbols=20

# compute log as we plot log(mean)
counts.ctrl.mean.log <- log(counts.ctrl.mean)
counts.egf.trtm.mean.log <- log(counts.egf.trtm.mean)
counts.egf.a66.mean.log <- log(counts.egf.a66.mean)


plot(x=counts.egf.trtm.mean.log, y=counts.egf.trtm.cv, 
     xlim=c(0, max(counts.ctrl.mean.log, counts.egf.trtm.mean.log, counts.egf.a66.mean.log, na.rm=TRUE)),
     ylim=c(0, max(counts.ctrl.cv, counts.egf.trtm.cv, counts.egf.a66.cv, na.rm=TRUE)),     
     pch=symbols,col='green',
     cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main="CV samples by treatment type",xlab=sprintf("log(mean)"),ylab=sprintf("CV"))

points(x=counts.egf.a66.mean.log, y=counts.egf.a66.cv,
     pch=symbols,col='magenta',
     cex=1.5,cex.lab=1.5,cex.axis=1.5)

points(x=counts.ctrl.mean.log, y=counts.ctrl.cv,
     pch=symbols,col='blue',
     cex=1.5,cex.lab=1.5,cex.axis=1.5)

legend(x="topright", legend = c('ctrl','egf a66','egf trtm'), col=c("blue","magenta","green"), pch=symbols)

dev.copy(png, paste("CV_",filename,"_samples_by_trtm_type",".png",sep=""))
dev.off()









############################################
# For each miRNA/condition, get the means and SD between replicates. Then plot X=mean, Y=CV
# - scatterplot with colors representing the times.
############################################

# separate groups
# tc: 0
counts.0 <- counts[, grep('_0_', colnames(counts))]
counts.0.mean <- apply(counts.0, 1, mean)
counts.0.sd <- apply(counts.0, 1, sd)
counts.0.cv <- counts.0.sd/counts.0.mean
# tc: 15
counts.15 <- counts[, grep('_15_', colnames(counts))]
counts.15.mean <- apply(counts.15, 1, mean)
counts.15.sd <- apply(counts.15, 1, sd)
counts.15.cv <- counts.15.sd/counts.15.mean
# tc: 40
counts.40 <- counts[, grep('_40_', colnames(counts))]
counts.40.mean <- apply(counts.40, 1, mean)
counts.40.sd <- apply(counts.40, 1, sd)
counts.40.cv <- counts.40.sd/counts.40.mean
# tc: 90
counts.90 <- counts[, grep('_90_', colnames(counts))]
counts.90.mean <- apply(counts.90, 1, mean)
counts.90.sd <- apply(counts.90, 1, sd)
counts.90.cv <- counts.90.sd/counts.90.mean
# tc: 180
counts.180 <- counts[, grep('_180_', colnames(counts))]
counts.180.mean <- apply(counts.180, 1, mean)
counts.180.sd <- apply(counts.180, 1, sd)
counts.180.cv <- counts.180.sd/counts.180.mean
# tc: 300
counts.300 <- counts[, grep('_300_', colnames(counts))]
counts.300.mean <- apply(counts.300, 1, mean)
counts.300.sd <- apply(counts.300, 1, sd)
counts.300.cv <- counts.300.sd/counts.300.mean



# Extract miRNA with log(mean)>=5 and CV>=0.5
counts.tc <- counts[, grep('_[[:digit:]]{1,3}_', colnames(counts))]
counts.tc.mean <- apply(counts.tc, 1, mean)
counts.tc.sd <- apply(counts.tc, 1, sd)
counts.tc.cv <- counts.tc.sd/counts.tc.mean
# filter by CV>=0.5
miRNA.highMean.highCV <- counts[sapply(counts.tc.cv, function(x) { x >= 0.5}), ]
# filter by log(mean)>=5
miRNA.highMean.highCV <- miRNA.highMean.highCV[sapply(log(counts.tc.mean), function(x) {x >= 5}), ]
write.csv(miRNA.highMean.highCV, file=paste(filename,"_miRNA_highMean_highCV.csv", sep=""), quote=FALSE)




counts.0.mean.log <- log(counts.0.mean)
counts.15.mean.log <- log(counts.15.mean)
counts.40.mean.log <- log(counts.40.mean)
counts.90.mean.log <- log(counts.90.mean)
counts.180.mean.log <- log(counts.180.mean)
counts.300.mean.log <- log(counts.300.mean)

symbols=20

plot(x=counts.0.mean.log, y=counts.0.cv, 
     xlim=c(0, max(counts.0.mean.log, counts.15.mean.log, counts.40.mean.log, counts.90.mean.log, counts.180.mean.log, counts.300.mean.log, na.rm=TRUE)),
     ylim=c(0, max(counts.0.cv, counts.15.cv, counts.40.cv, counts.90.cv, counts.180.cv, counts.300.cv, na.rm=TRUE)),     
     pch=symbols,col='blue',
     cex=1.5,cex.lab=1.5,cex.axis=1.5,
     main="CV samples by time point",xlab=sprintf("log(mean)"),ylab=sprintf("CV"))

points(x=counts.15.mean.log, y=counts.15.cv,
       pch=symbols,col='magenta',
       cex=1.5,cex.lab=1.5,cex.axis=1.5)
points(x=counts.40.mean.log, y=counts.40.cv,
       pch=symbols,col='green',
       cex=1.5,cex.lab=1.5,cex.axis=1.5)
points(x=counts.90.mean.log, y=counts.90.cv,
       pch=symbols,col='orange',
       cex=1.5,cex.lab=1.5,cex.axis=1.5)
points(x=counts.180.mean.log, y=counts.180.cv,
       pch=symbols,col='cyan',
       cex=1.5,cex.lab=1.5,cex.axis=1.5)
points(x=counts.300.mean.log, y=counts.300.cv,
       pch=symbols,col='brown',
       cex=1.5,cex.lab=1.5,cex.axis=1.5)

legend(x="topright", legend = c('0m','15m','40m','90m','180m','300m'), col=c("blue","magenta","green", "orange", "cyan", "brown"), pch=symbols)

dev.copy(png, paste("CV_",filename,"_samples_by_time_course",".png",sep=""))
dev.off()


# Mark and label the points of the plot manually
#identify(x=c(log(counts.0.mean),log(counts.15.mean),log(counts.40.mean),log(counts.90.mean),log(counts.180.mean),log(counts.300.mean)), 
#         y=c(counts.0.cv,counts.15.cv,counts.40.cv,counts.90.cv,counts.180.cv,counts.300.cv))


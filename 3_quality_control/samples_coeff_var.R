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




# Calculate the coefficient of variance vs mean for treatment / time point.


# For each miRNA/condition, get the means and SD between replicates. Then plot X=mean, Y=CV
# You should get plenty of points, that should follow a classical square decay.
# Most of the variation is at very low mean values because noise decreases as the sqrt(signal).
# Two plots:
# - scatterplot with colors representing noEGF, EGF-treated and EGF+A66
# - scatterplot with colors representing the times.



source('../utilities/plots.R')
library(ggplot2)



# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered"
suffix <-".csv"


# load counts
counts <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=TRUE,header=TRUE,row.names=1)



############################################
# For each miRNA/condition, get the means and SD between replicates. Then plot X=mean, Y=CV
# - scatterplot with colors representing noEGF, EGF WT and EGF A66
############################################

# separate groups
# noEGF
counts.noEGF <- counts[, grep('noEGF', colnames(counts))]
counts.noEGF <- unique(counts.noEGF)
counts.noEGF.mean <- apply(counts.noEGF, 1, mean)
counts.noEGF.sd <- apply(counts.noEGF, 1, sd)
counts.noEGF.cv <- counts.noEGF.sd/counts.noEGF.mean
# EGF-wt
counts.egf.wt <- counts[, grep('miWT', colnames(counts))]
counts.egf.wt.mean <- apply(counts.egf.wt, 1, mean)
counts.egf.wt.sd <- apply(counts.egf.wt, 1, sd)
counts.egf.wt.cv <- counts.egf.wt.sd/counts.egf.wt.mean
# EGF-a66
counts.egf.a66 <- counts[, grep("miA66_[[:digit:]]", colnames(counts))]
counts.egf.a66.mean <- apply(counts.egf.a66, 1, mean)
counts.egf.a66.sd <- apply(counts.egf.a66, 1, sd)
counts.egf.a66.cv <- counts.egf.a66.sd/counts.egf.a66.mean

# compute log as we plot log(mean)
counts.noEGF.mean.log <- log1p(counts.noEGF.mean)
counts.egf.wt.mean.log <- log1p(counts.egf.wt.mean)
counts.egf.a66.mean.log <- log1p(counts.egf.a66.mean)


# Create a data frame containing what we need.
df.samples.by.trmt <- data.frame(mean=c(counts.noEGF.mean.log, 
                                        counts.egf.wt.mean.log, 
                                        counts.egf.a66.mean.log), 
                                 cv=c(counts.noEGF.cv, 
                                      counts.egf.wt.cv, 
                                      counts.egf.a66.cv),
                                 treatment=factor(c(rep('no EGF', length(counts.noEGF.mean.log)), 
                                                    rep('EGF WT', length(counts.egf.wt.mean.log)), 
                                                    rep('EGF A66', length(counts.egf.a66.mean.log)))))

# Plot
g <- ggplot(data=df.samples.by.trmt, aes(x=mean, y=cv, color=treatment)) +
     geom_point(size=0.8) +
     ggtitle('CV samples by treatment type') + xlab('log(mean)') + ylab('cv') +
     theme_basic() + 
     theme(axis.text=element_text(size=12),
           axis.title.x=element_text(size=14),
           axis.title.y=element_text(size=14))
ggsave(paste0("CV_",filename,"_samples_by_strain_type",".png"), width=5, height=4, dpi=300)







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


# compute log as we plot log(mean)
counts.0.mean.log <- log1p(counts.0.mean)
counts.15.mean.log <- log1p(counts.15.mean)
counts.40.mean.log <- log1p(counts.40.mean)
counts.90.mean.log <- log1p(counts.90.mean)
counts.180.mean.log <- log1p(counts.180.mean)
counts.300.mean.log <- log1p(counts.300.mean)


# Create a data frame containing what we need.
df.samples.by.tc <- data.frame(mean=c(counts.0.mean.log, 
                                      counts.15.mean.log,
                                      counts.40.mean.log, 
                                      counts.90.mean.log, 
                                      counts.180.mean.log, 
                                      counts.300.mean.log), 
                                cv=c(counts.0.cv, 
                                     counts.15.cv, 
                                     counts.40.cv, 
                                     counts.90.cv, 
                                     counts.180.cv, 
                                     counts.300.cv),
                                time=factor(c(rep('0m', length(counts.0.mean.log)), 
                                              rep('15m', length(counts.15.mean.log)), 
                                              rep('40m', length(counts.40.mean.log)), 
                                              rep('90m', length(counts.90.mean.log)), 
                                              rep('180m', length(counts.180.mean.log)), 
                                              rep('300m', length(counts.300.mean.log)))))

# Plot
g <- ggplot(data=df.samples.by.tc, aes(x=mean, y=cv, color=time)) +
    geom_point(size=0.8) +
    ggtitle('CV samples by time course') + xlab('log(mean)') + ylab('cv') +
    theme_basic() + 
    theme(axis.text=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))
ggsave(paste0("CV_",filename,"_samples_by_time_course",".png"), width=5, height=4, dpi=300)






##############################################
# Extract miRNA with log(mean)>=5 and CV>=0.5
##############################################

counts.tc <- counts[, grep('_[[:digit:]]{1,3}_', colnames(counts))]
counts.tc.mean <- apply(counts.tc, 1, mean)
counts.tc.sd <- apply(counts.tc, 1, sd)
counts.tc.cv <- counts.tc.sd/counts.tc.mean
# filter by CV>=0.5
miRNA.highMean.highCV <- counts[sapply(counts.tc.cv, function(x) { x >= 0.5}), ]
# filter by log(mean)>=5
miRNA.highMean.highCV <- miRNA.highMean.highCV[sapply(log(counts.tc.mean), function(x) {x >= 5}), ]
write.csv(miRNA.highMean.highCV, file=paste(filename,"_miRNA_highMean_highCV.csv", sep=""), quote=FALSE)





# Mark and label the points of the plot manually
#identify(x=c(log(counts.0.mean),log(counts.15.mean),log(counts.40.mean),log(counts.90.mean),log(counts.180.mean),log(counts.300.mean)), 
#         y=c(counts.0.cv,counts.15.cv,counts.40.cv,counts.90.cv,counts.180.cv,counts.300.cv))


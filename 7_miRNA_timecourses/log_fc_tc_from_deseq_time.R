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




# Plot time courses of all miRNA. 
# - X is Time
# - Y is log fold change as calculated by DESeq (contrasts: Time)
# - colour is (A66 vs WT) or (PAM clustering)



library(reshape2)
library(ggplot2)

source('../utilities/plots.R')


################
# Load data sets
################

# select the files containing the data
location.time <- "../4_deseq__time"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__time__0_"
suffix <-".csv"
padj.deseq.time <- "__padj_005"
time <- list("15", "40", "90", "180", "300")


location.strain <- "../4_deseq__strain"
filename.deseq.strain <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66"
suffix <-".csv"
padj.strain.thres <- 0.01


clusters <- 2
location.pam <- "../6_clustering"
filename.pam <- paste0("summarised_mirna_counts_after_mapping_filtered_scaled_pam", clusters, "_labels")
suffix <-".csv"



#############################
# Merge time course data sets
#############################

# Load one DESeq data set in order to get the miRNA names.
miRNAs <- row.names(read.table(paste0(location.time,"/",filename.deseq.time, time[1], suffix), sep=",",fill=T,header=T,row.names=1))

# Add log2FoldChange_time0
df.tc <- data.frame(t0=rep(0, length(miRNAs)))
# Add row names
rownames(df.tc) <- miRNAs

# Add log2FoldChange for the available time points.
for(i in seq(1, length(time))) {
  df.tmp <- read.table(paste0(location.time,"/",filename.deseq.time, time[i], suffix), sep=",",fill=T,header=T,row.names=1)[,2]
  # Add log2FoldChange_timeX
  df.tc <- cbind(df.tc, df.tmp)
}

# Set column names (time points)
colnames(df.tc) <- c('0', paste0(time))




#############################################
# Extract DESeq:strain for colour information
#############################################

# Load strain information
df.strain <- read.table(paste0(location.strain,"/",filename.deseq.strain, suffix), sep=",",fill=T,header=T,row.names=1)

# Extract the names of the significant miRNA. This will be used as a filter later on.
df.strain.signif.mirna <- rownames(subset(df.strain, padj < padj.strain.thres))
write.csv(df.strain.signif.mirna, file=paste0(filename.deseq.strain ,"_DESeq_strain__WT_A66__padj_", gsub("\\.","",padj.strain.thres), suffix))

# Extract the log2 fold change
df.strain <- df.strain[,2]


##################
# Melt data frames
##################

df <- df.tc
# add miRNA as new row
df$miRNA <- rownames(df)
# Add strain (colour)
df$strain <- df.strain
# Melt times so that we have 'miRNA', 'strain', 'variable' (time-cols), 'value' (time-vals)
df.melt <- melt(df, id=c('miRNA','strain'))


# Due to the high number of miRNA, we also filter those having significant change in DESeq contrast `Strain`.
df.signif <- df[df.strain.signif.mirna,]
df.melt.signif <- melt(df.signif, id=c('miRNA','strain'))



#####################################
# Plot miRNA using Strain for colours
#####################################

# Plot all miRNAs
g <- ggplot(data=df.melt, aes(x=variable, y=value, group=miRNA, color=strain)) + 
  geom_line() + 
  scale_colour_gradient2(low="red", mid="lightgrey", high="navyblue") +
  ggtitle('miRNA expression with strain log2 FC') +
  xlab('Time [m]') + 
  ylab('log2 fold change') +
  theme_basic()
ggsave(paste0("miRNA_log_fc_tc__deseq_time__w_strain.png"), width=4, height=4, dpi=300)


# Plot miRNAs with significant change in Strain
g <- ggplot(data=df.melt.signif, aes(x=variable, y=value, group=miRNA, color=strain)) + 
  geom_line() + 
  scale_colour_gradient2(low="red", mid="lightgrey", high="navyblue") +
  ggtitle('miRNA expression with strain log2 FC') +
  xlab('Time [m]') + 
  ylab('log2 fold change') +
  theme_basic()
ggsave(paste0("miRNA_log_fc_tc__deseq_time__w_strain__padj_", gsub('\\.','',as.character(padj.strain.thres)), ".png"), width=4, height=4, dpi=300)






##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



#############################################
# Extract PAM cluster labels for colour information
#############################################

# Load strain information
df.pam <- read.table(paste0(location.pam,"/",filename.pam, suffix), sep=",",fill=T,header=T,row.names=1)


##################
# Melt data frames
##################

df <- df.tc
# add miRNA as new row
df$miRNA <- rownames(df)
# Add PAM cluster (colour)
# pam is converted to factor so that it is not automatically interpreted as a continuous variable
df$pam <- factor(df.pam[,1])
# Melt times so that we have 'miRNA', 'pam', 'variable' (time-cols), 'value' (time-vals)
df.melt <- melt(df, id=c('miRNA','pam'))



#############################################
# Plot miRNA using PAM clustering for colours
#############################################

g <- ggplot(df.melt, aes(x=variable, y=value, group=miRNA, color=pam)) + 
  geom_line() + 
  ggtitle('miRNA expression with PAM clustering') +
  xlab('Time [m]') + 
  ylab('log2 fold change') +
  theme_basic()
ggsave(paste0("miRNA_log_fc_tc__deseq_time__w_pam_clust.png"), width=4, height=4, dpi=300)


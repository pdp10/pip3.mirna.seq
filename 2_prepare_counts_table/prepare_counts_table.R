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





# This script prepares the counts table. Here we filter the table.


library(data.table)
library(DESeq2)





# Create a data table from a counts matrix
counts_2_data_table <- function(m) {
  d <- as.data.frame(m)
  d$miRNA <- rownames(d)
  d <- melt(d)
  d <- as.data.table(d)
  d[, `:=`(strain, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 1))]
  d[, `:=`(time, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 2))]
  d[, `:=`(replicate, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 3))]
  d <- d[, list(miRNA, value, strain, time, replicate)]
  return(d)
}



log_counts <- function(counts, counts.metadata, filename) {
  
  # Prepare DESeq dataset
  #######################
  # create the dataset for DESeq. We use this as an annotated SummarizedExperiment object
  # use assay(se), colData(se), design(se) to read the most important info
  se <- DESeqDataSetFromMatrix(countData = counts, colData=counts.metadata, design = ~ 1)
  
  # Estimate dispersions
  ######################
  # estimate the size factors from counts, counts.metadata
  se <- estimateSizeFactors(se)
  se <- estimateDispersions(se)
  #str(fitInfo(se))
  png(filename, width=2000, height=2000, pointsize=14, res=300)
  plotDispEsts(se)
  dev.off()
  
  # Compute rlog as implemented in DESeq2. The rlog tends to work well on small datasets (n < 30)
  # In the below function calls, we specified blind = FALSE, which means that differences between cell lines and 
  # treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
  # The experimental design is not used directly in the transformation, only in estimating the global amount of 
  # variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
  rld <- rlog(se)
  # convert rownames to a new column
  rld.counts <- assay(rld)

  return(rld.counts)
}


# NOTE: counts must be log trasformed before running this function. 
# Alternatively, compute the median
mean_counts_by_reps <- function(counts) {
  
  # retrieve the experiment groups (without repeat number). We only need to drop the last two characters from 
  # the name of each column. 
  groups <- gsub('.{2}$', '', colnames(counts))
  
  # remove duplicates
  groups <- groups[!duplicated(groups)]
  
  # calculate the mean for each group of repeats.
  counts.mean <- matrix(0, ncol=length(groups), nrow = nrow(counts))
  
  # extract the columns of replicates for each group and calculate the mean by row. 
  # This new column will form a new data set with groups as column names.
  for(i in seq(1:length(groups))) {
    # extract the columns of replicates
    counts.g <- counts[, grepl(groups[i], colnames(counts))]
    # calculate the mean per row
    counts.g.mean <- apply(counts.g, 1, mean)
    # add this column to the mean data frame
    counts.mean[, i] <- counts.g.mean
  }
  
  counts.mean <- data.frame(counts.mean)
  rownames(counts.mean) <- rownames(counts)
  colnames(counts.mean) <- groups
  
  return(counts.mean)
}





#########
# My Data
#########

# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping"
filename.counts.metadata <- "summarised_mirna_counts_after_mapping_filtered_metadata"
suffix <-".csv"



#######################
# Load the counts table
#######################

# load counts
counts <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=TRUE, header=TRUE, row.names=1)
# load counts metadata
counts.metadata <- read.table(paste0(location,"/",filename.counts.metadata,suffix), sep=",",fill=T,header=T,row.names=1)
# cast the column `time` from numeric to factor
counts.metadata$time <- as.factor(counts.metadata$time)
# cast the column `replicate` from numeric to factor
counts.metadata$replicate <- as.factor(counts.metadata$replicate)



###########
# FILTERING
###########

# remove 'unmap' counts (unmapped)
counts<-counts[, -grep('unmap', colnames(counts))]

# remove 'base' counts (these seem to be replicates but have a wrong base)
counts <- counts[, -grep('base', colnames(counts))]

# Remove lowly expressed genes (with less than 1 reads in every samples)
# This requires some care. 
maxExp <- do.call(pmax, counts)
counts <- counts[maxExp>10,]

# remove all 0-expressed genes
#counts <- counts[rowSums(counts) > 0, ]



##############################################################
# Get the miRNA names after filtering and replace `-` with `_`
##############################################################
mirna <- gsub('-','_', rownames(counts))




########################
# Save the counts matrix
########################
print('filter counts matrix')
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.mirna <- cbind(miRNA = mirna, counts)
# write the filtered data to file
write.csv(counts.mirna, file=paste0(location, "/", filename, "_filtered", suffix), row.names=FALSE)




#############################
# Compute the log of the data
#############################
print('rlog counts matrix')
# We use DESeq2 rlog function.

# log the data using rlog
counts.rlog <- log_counts(counts, counts.metadata, filename=paste0(filename,"_DESeq_dispersions_estimates",".png"))
counts.rlog.mirna <- cbind(miRNA = mirna, counts.rlog)
write.csv(counts.rlog.mirna, file=paste0(location,"/",filename,"_filtered_rlog",suffix), row.names=FALSE)



###########################
# Scale the log of the data (clustering, pca, etc)
###########################
print('scale the rlog counts matrix')
# scale counts to N(0,1)
counts.rlog.scaled <- t(scale(t(counts.rlog), center=TRUE, scale=TRUE))
counts.rlog.scaled.mirna <- cbind(miRNA = mirna, counts.rlog.scaled)
# write the data to file
write.csv(counts.rlog.scaled.mirna, file=paste0(location, "/", filename, "_filtered_rlog_scaled", suffix), row.names=FALSE)





############### mean of the replicates #################



####################################
# Compute the mean of the log counts  (Mean of the replicates)
####################################
print('mean of the rlog counts matrix')
counts.rlog.mean <- mean_counts_by_reps(counts.rlog)
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.rlog.mean.mirna <- cbind(miRNA = mirna, counts.rlog.mean)
# write the data to file
write.csv(counts.rlog.mean.mirna, file=paste0(location, "/", filename, "_filtered_rlog_mean", suffix), row.names=FALSE)


###########################
# Scale the log mean of the data (clustering, pca, etc)
###########################
print('scale the mean of the rlog counts matrix')
# scale counts to N(0,1)
counts.rlog.mean.scaled <- t(scale(t(counts.rlog.mean), center=TRUE, scale=TRUE))
counts.rlog.mean.scaled.mirna <- cbind(miRNA = mirna, counts.rlog.mean.scaled)
# write the data to file
write.csv(counts.rlog.mean.scaled.mirna, file=paste0(location, "/", filename, "_filtered_rlog_mean_scaled", suffix), row.names=FALSE)






#########################
# Counts data table
#########################
print('create counts matrix table')
counts.data.table <- counts_2_data_table(counts)
# write the filtered data to file
write.csv(counts.data.table, file=paste0(location, "/", filename, "_filtered_data_table", suffix), row.names=FALSE)



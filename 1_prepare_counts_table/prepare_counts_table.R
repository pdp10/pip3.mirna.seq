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


# select the file countaining the data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping"
suffix <-".csv"



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





#######################
# Load the counts table
#######################

# load counts
counts <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=TRUE, header=TRUE, row.names=1)



###########
# FILTERING
###########

# remove 'unmap' counts (unmapped)
counts<-counts[, -grep('unmap', colnames(counts))]

# remove 'base' counts (these seem to be replicates but have a wrong base)
counts<-counts[, -grep('base', colnames(counts))]

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
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.mirna <- cbind(miRNA = mirna, counts)

# write the filtered data to file
write.csv(counts.mirna, file=paste0(location, "/", filename, "_filtered", suffix), row.names=FALSE)





#########################
# Log the counts matrix
#########################

counts.log <- log1p(counts)
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.log.mirna <- cbind(miRNA = mirna, counts.log)
# write the filtered data to file
write.csv(counts.log.mirna, file=paste0(location, "/", filename, "_filtered_log", suffix), row.names=FALSE)





#########################
# Scale the counts matrix (clustering, pca, etc)
#########################

# scale counts to N(0,1)
counts.scaled <- t(scale(t(counts), center=TRUE, scale=TRUE))
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.scaled.mirna <- cbind(miRNA = mirna, counts.scaled)
# write the filtered data to file
write.csv(counts.scaled.mirna, file=paste0(location, "/", filename, "_filtered_scaled", suffix), row.names=FALSE)







#########################
# Log and scaled the counts matrix
#########################

# scale counts to N(0,1)
counts.log.scaled <- t(scale(t(counts.log), center=TRUE, scale=TRUE))
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.log.scaled.mirna <- cbind(miRNA = mirna, counts.log.scaled)
# write the filtered data to file
write.csv(counts.log.scaled.mirna, file=paste0(location, "/", filename, "_filtered_log_scaled", suffix), row.names=FALSE)






#########################
# Counts data table
#########################

counts.data.table <- counts_2_data_table(counts)
# write the filtered data to file
write.csv(counts.data.table, file=paste0(location, "/", filename, "_filtered_data_table", suffix), row.names=FALSE)






################################
# Scale the median counts matrix
################################
# Because this is gene expression, one has to take the geometric mean, 
# which in the case of a lognormal law is equal to the median. 
# Then you can do all the transformations you want.


### 1) FIRST, WE CALCULATE THE MEDIAN

# retrieve the experiment groups (without repeat number). We only need to drop the last two characters from 
# the name of each column. 
groups <- gsub('.{2}$', '', colnames(counts))

# remove duplicates
groups <- groups[!duplicated(groups)]

# calculate the median for each group of repeats.
counts.median <- matrix(0, ncol=length(groups), nrow = nrow(counts))

# extract the columns of replicates for each group and calculate the median by row. 
# This new column will form a new data set with groups as column names.
for(i in seq(1:length(groups))) {
  # extract the columns of replicates
  counts.g <- counts[, grepl(groups[i], colnames(counts))]
  # calculate the median per row
  counts.g.median <- apply(counts.g, 1, median)
  # add this column to the median data frame
  counts.median[, i] <- counts.g.median
}

counts.median <- data.frame(counts.median)
rownames(counts.median) <- rownames(counts)
colnames(counts.median) <- groups

# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.median.mirna <- cbind(miRNA = mirna, counts.median)
# write the filtered data to file
write.csv(counts.median.mirna, file=paste0(location, "/", filename, "_filtered_median", suffix), row.names=FALSE)



### 2) SECOND, WE SCALE
# scale counts.median to N(0,1)
counts.median.scaled <- t(scale(t(counts.median), center=TRUE, scale=TRUE))


# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.median.scaled.mirna <- cbind(miRNA = mirna, counts.median.scaled)
# write the filtered data to file
write.csv(counts.median.scaled.mirna, file=paste0(location, "/", filename, "_filtered_median_scaled", suffix), row.names=FALSE)




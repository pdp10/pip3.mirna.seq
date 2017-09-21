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
  d$mir <- rownames(d)
  d <- melt(d)
  d <- as.data.table(d)
  d[, `:=`(strain, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 1))]
  d[, `:=`(time, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 2))]
  d[, `:=`(replicate, sapply(strsplit(as.character(d[, variable]), "_"), "[[", 3))]
  d <- d[, list(mir, value, strain, time, replicate)]
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



#########################
# Scale the counts matrix (clustering, pca, etc)
#########################

# scale counts to N(0,1)
counts.scaled <- t(scale(t(counts), center=TRUE, scale=TRUE))


#########################
# Log the counts matrix
#########################

counts.log <- log(counts+1)


#########################
# Log the counts matrix
#########################

counts.data.table <- counts_2_data_table(counts)



############
# Save files
############

mirna <- gsub('-','_', rownames(counts))

# COUNTS
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts <- cbind(miRNA = mirna, counts)

# write the filtered data to file
write.csv(counts, file=paste0(location, "/", filename, "_filtered", suffix), row.names=FALSE)



# COUNTS.SCALED
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.scaled <- cbind(miRNA = mirna, counts.scaled)

# write the filtered data to file
write.csv(counts.scaled, file=paste0(location, "/", filename, "_filtered_scaled", suffix), row.names=FALSE)



# LOG COUNTS
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.log <- cbind(miRNA = mirna, counts.log)

# write the filtered data to file
write.csv(counts.log, file=paste0(location, "/", filename, "_filtered_log", suffix), row.names=FALSE)





# COUNTS DATA TABLE
# convert rownames to a new column. NOTE: '-' is replaced with '_' in miRNA names
counts.log <- cbind(miRNA = mirna, counts.log)

# write the filtered data to file
write.csv(counts.data.table, file=paste0(location, "/", filename, "_filtered_data_table", suffix), row.names=FALSE)


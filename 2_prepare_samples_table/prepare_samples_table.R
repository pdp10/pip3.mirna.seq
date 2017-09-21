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





# This script generates the metadata table from the filtered counts table. 
# The metadata table is also called sample table in DESEq.


# select the file countaining the filtered counts data
location <- "../data"
filename <- "summarised_mirna_counts_after_mapping_filtered"
suffix <-".csv"


# load counts
counts <- read.table(paste0(location,"/",filename,suffix), sep=",",fill=T,header=T,row.names=1)


# Extract everything (time points, noEGF) which is between two '_' . Arbitrary strings before and after (.*)
# replace the whole string with the contents of the first capturing group: \\1 (we need to escape the backslash, hence the double backslash)
time <- gsub(".*[_]([^.]+)[_].*", "\\1", colnames(counts))
# now replace 'noEGF' with '300', as noEGF sample was taken at 300m. 
time <- gsub('noEGF', '300', time)

# use grepl to create a logical vector whether noEGF appears or not in n.
egf <- ifelse(grepl('noEGF', colnames(counts)), 'noEGF', 'EGF')

# use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(counts)), 'WT', 'A66')

replicate <- gsub(".*[_]([^.]+)", "\\1", colnames(counts))
replicate <- paste0(replicate)

# make the counts metadata table
counts.metadata <- data.frame(
    samples = colnames(counts),
    strain = strain,
    time = time,
    egf = egf, 
    replicate = replicate)


# we sort the metadata by treatment (descending) and time (ascending). 
# we use rank() otherwise the function order() complains as strain is a factor.
#counts.metadata <- counts.metadata[order(-rank(strain), time),]

#print(counts.metadata) 



# write the counts metadata to file
write.csv(counts.metadata, file=paste0(location,"/",filename, "_metadata", suffix), row.names=FALSE)


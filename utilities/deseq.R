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





# Create a data frame of Deseq2:log2FoldChange time courses. 
# `files` is a vector of Deseq2 results, each representing a time point. 
deseq_lfc_time_course_df <- function(files) {
  
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
  
  # Remove the first fake column t0.
  df.tc <- subset(df.tc, select = -c(t0) )
  
  # Set column names (time points)
  colnames(df.tc) <- c(paste0(time))
  
  return(df.tc)
}

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





# Extract and assemble the statistically significant miRNA computed from time course volcano plots 
# in one table. miRNA are sorted alphabetically. 
# The resulting summary table has columns: miRNA, log2foldchange, negLog10Padj, and time 



library("naturalsort")
source('../utilities/plots.R')


################
# Load data sets
################

# select the files containing the data
locations <- c(
  "../4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain",
  "../4_deseq__strain_vs_time_per_TimePoint_interaction_WaldTest",
  "../4_deseq__time"
)

filename.patterns <- c(
  "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66_at_t",
  "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66_at_t",
  "summarised_mirna_counts_after_mapping_filtered_results_DESeq__time__0_"
)

suffix <-".csv"

time1 <- list("0", "15", "40", "90", "180", "300")
time2 <- list("15", "40", "90", "180", "300")
time3 <- list("15", "40", "90", "180", "300")
times <- list(time1, time2, time3)

padj <- c("005", "005", "005")

lfc <- c("02", "05", "01")


#################################################################################
# Combine signif miRNA from DESeq WT vs A66 (strain contrast) for each time point
#################################################################################

for(i in 1:length(filename.patterns)) {
  
  # columns are miRNA, log2FoldChange, negLog10Padj, time
  df.complete <- data.frame(matrix(ncol=4, nrow = 0))
  
  # Read and stack together the signif miRNAs for each time point 
  for(j in 1:length(times[[i]])) {
    df <- read.table(paste0(locations[i], "/", filename.patterns[i], times[[i]][[j]], "__padj_", padj[i], "_lfc_", lfc[i], suffix), sep=",", fill=T, header=T)
    df <- na.omit(df)
    df$time <- rep(times[[i]][[j]], nrow(df))
    df.complete <- rbind(df.complete, df)
  }
  colnames(df.complete) <- c('miRNA','log2FoldChange', 'negLog10Padj', 'time')
  df.complete <- df.complete[naturalorder(df.complete$miRNA),]

  # rename rows
  row.names(df.complete) <- 1:nrow(df.complete)
  #print(df.complete)    
  write.csv(df.complete, file=paste0(substring(locations[i], 6), '__', filename.patterns[i],"_miRNA__padj_", padj[i], "_lfc_", lfc[i], suffix), quote=FALSE, row.names=FALSE)
  
}

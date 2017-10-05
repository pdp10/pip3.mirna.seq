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




# Run PAM clustering (Partition around medoids)
run_pam <- function(df, k, name, scale.=FALSE, centre=FALSE) {

  pamx <- pam(df, k, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
  # Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
  write.csv(pamx$clustering, file=paste0(name, "_labels.csv"), quote=FALSE)

  return (pamx)
}


# run pca
run_pca <- function(df, name, scale.=FALSE, centre=FALSE) {
  
  # COMPUTE PCA
  pca <- prcomp(t(df), scale.=scale., center=centre)
  # Write the calculated PCA `rotation` (PCA load). This is the table miRNAs (rows) vs PCAs (cols).
  write.csv(pca$rotation, file=paste0(name, '_pca_rotation.csv'), quote=FALSE)
  
  # Plot the variance of PCA components
  png(paste0(name,"_pca_variances",".png"), width=2000, height=1500, res=300)
  plot(pca, type = "l", main='Variance of PCA components')
  dev.off()
  
  return(pca)
}


# filter PCA rotation by threshold on a specified component
filt_pca <- function(pca, component=1, threshold=0.05) {
  pca.rot <- pca$rotation
  # filter the table
  pca.rot.filt <- pca.rot[ pca.rot[, component] < threshold | pca.rot[, component] > threshold, ]
  
  # return the list of filtered ids
  return(rownames(pca.rot.filt))
}




# Split the count matrix depending on cluster class
# Return a list of lists
split_counts_matrix_by_pam_class <- function(df.wpam, pam.classes) {
  
  df.cols <- colnames(df.wpam)
  df.wt.a66 <- list()
  df.a66.noEGF <- list()
  
  
  for(pam.class in pam.classes) {
    
    # Separate the miRNAs by pam.class
    ##################################
    # the inner subset filters the miRNAs with pam.class, whereas the outer subset drops the pam column as it is no longer necessary.
    df.wt.a66[[pam.class]] <- subset(subset(df.wpam, pam==pam.class), select= -c(pam))
    
    # Separate and stack WT, A66, and A66-noEGF
    ###########################################
    wt <- data.frame(df.wt.a66[[pam.class]][, grepl("miWT", df.cols), drop=FALSE])
    a66 <- data.frame(df.wt.a66[[pam.class]][, grepl("miA66", df.cols), drop=FALSE])
    a66.noEGF <- data.frame(df.wt.a66[[pam.class]][, grepl("noEGF", df.cols), drop=FALSE])
    
    # remove column `miA66_noEGF`` from miA66 data frame as this is treated separately.
    a66 <- subset(a66, select = -miA66_noEGF)
    
    # rename columns so that they only contain times (useful for melting later on)
    colnames(wt) <- gsub('miWT_', '', colnames(wt))
    colnames(a66) <- gsub('miA66_', '', colnames(a66))
    # this is renamed manually
    colnames(a66.noEGF) <- c('300')
    
    # add a row for the colour
    wt <- cbind (wt, colour=rep('WT', nrow(wt)))
    a66 <- cbind (a66, colour=rep('A66', nrow(a66)))
    a66.noEGF <- cbind (a66.noEGF, colour=rep('A66_noEGF', nrow(a66.noEGF)))
    
    # stack the two data frames (wt, a66) and (a66noEGF)
    #####################################
    df.wt.a66[[pam.class]] <- rbind(wt, a66)
    # save data frames
    write.csv(df.wt.a66[[pam.class]], file=paste0('miRNA_wt_a66_tc_cluster_', pam.class, '.csv'), row.names=T, quote=FALSE)
    # add id column
    df.wt.a66[[pam.class]]$id <- 1:nrow(df.wt.a66[[pam.class]])

    # we keep a66.noEGF separately because it has different columns
    df.a66.noEGF[[pam.class]] <- a66.noEGF
    # save data frame
    write.csv(df.a66.noEGF[[pam.class]], file=paste0('miRNA_a66noEGF_tc_cluster_', pam.class, '.csv'), row.names=T, quote=FALSE)
    # add id column
    df.a66.noEGF[[pam.class]]$id <- 1:nrow(df.a66.noEGF[[pam.class]])
  }
  return(list(df.wt.a66, df.a66.noEGF))
}



# Calculate the median for each PAM class and strain
median_for_each_pam_class_and_strain <- function(df.wt.a66, df.a66.noEGF, pam.classes) {
  
  df.wt.a66.median <- list()
  df.a66.noEGF.median <- list()
  
  for(pam.class in pam.classes) {
    dt.median <- median_df_by_colour_id(df=subset(df.wt.a66[[pam.class]], select=-c(id)), 
                                        filename=paste0('miRNA_wt_a66_tc_cluster_', pam.class, '_median.csv'))
    df.wt.a66.median[[pam.class]] <- data.frame(dt.median)
    colnames(df.wt.a66.median[[pam.class]]) <- gsub("^X", "",  colnames(df.wt.a66.median[[pam.class]]))

    dt.median <- median_df_by_colour_id(df=subset(df.a66.noEGF[[pam.class]], select=-c(id)), 
                                        filename=paste0('miRNA_a66noEGF_tc_cluster_', pam.class, '_median.csv'))
    df.a66.noEGF.median[[pam.class]] <- data.frame(dt.median)
    colnames(df.a66.noEGF.median[[pam.class]]) <- gsub("^X", "",  colnames(df.a66.noEGF.median[[pam.class]]))
  }
  
  return(list(df.wt.a66.median, df.a66.noEGF.median))
}



# create data table of the median of rows having the same the colour id (column) 
# return a data table
median_df_by_colour_id <- function(df, filename) {
  # convert data.frame to data.table
  dt <- data.table(df)  
  # median every column by `colour` class (and for each time point)
  dt.median <- dt[, lapply(.SD, median), by=colour]
  # save the data.table
  write.csv(dt.median, file=filename, row.names=F, quote=FALSE)
  return(dt.median)
}

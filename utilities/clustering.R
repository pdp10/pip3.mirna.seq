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
run_pam <- function(df, k, name, scale., centre) {

  pamx <- pam(df, k, diss=FALSE, metric="euclidean", stand=FALSE, do.swap=TRUE)
  # Write the calculated clustering. This is the table miRNAs (rows) vs ClusterLabels (cols).
  write.csv(pamx$clustering, file=paste0(name,k,"_labels.csv"), quote=FALSE)

  return (pamx)
}


# run pca
run_pca <- function(df, name, scale., centre) {
  
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


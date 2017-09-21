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




# Compare DESeq significant miRNA (contrast is strain: WT vs A66) against the selected miRNA with high 
# variance in PCA:PC2, the principal component where WT and A66 are separated.
# In the first part the counts table is filtered by DESeq significant miRNA. A PCA is run on this filtered counts table.
# This PCA also shows separation between WT and A66 along PC2.
# In the second part the miRNA with variance outside the interval [-pca.thres, pca.thres] in the PCA loadings table are selected. This list 
# is then compared with the significant miRNA from DESeq. Venn Diagram is generated.



library(factoextra) # eigenvalues
library(ggfortify) # autoplot
library(VennDiagram)

library("naturalsort")
source('../utilities/plots.R')




################
# Load data sets
################

# Counts file (we import the scaled data set)
location.counts <- "../data"
filename.counts <- "summarised_mirna_counts_after_mapping_filtered_scaled"

# DESeq file
location.deseq <- "../4_deseq__strain"
filename.deseq <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66__padj_005_lfc_01"

# PCA loadings file
location.pca <- "../3_quality_control"
filename.pca <- "summarised_mirna_counts_after_mapping_filtered_scaled_PCA_rotation"

suffix <-".csv"


# read tables
df.counts <- read.table(paste0(location.counts, "/", filename.counts, suffix), sep=",", fill=T, header=T)

df.deseq <- read.table(paste0(location.deseq, "/", filename.deseq, suffix), sep=",", fill=T, header=T)
df.deseq <- na.omit(df.deseq)
colnames(df.deseq)[1] <- 'miRNA'

df.pca <- read.table(paste0(location.pca, "/", filename.pca, suffix), sep=",", fill=T, header=T)
df.pca <- na.omit(df.pca)
colnames(df.pca)[1] <- 'miRNA'



#############################################################################################
# "=>" Generate PCAs from DESeq statistically significant miRNA (contrast: WT vs A66 in strain)
#############################################################################################

# The purpose is to see whether the data separation detected on PC2 using PCA for all miRNA 
# is still present within the selected significan miRNA.



# Create a counts table of significant miRNAs
df.counts.signif <- df.counts[df.counts$miRNA %in% df.deseq$miRNA,]
rownames(df.counts.signif) <- df.counts.signif$miRNA
df.counts.signif <- df.counts.signif[,-1]

# Write the count table filtered by DESeq significant miRNA
write.csv(df.counts.signif, file=paste0(filename.counts,"__counts_table_deseq_signif_mirna", suffix), quote=FALSE)
# Write the significant miRNA
df.counts.signif.miRNA <- rownames(df.counts.signif)
write.csv(df.counts.signif.miRNA, file=paste0(filename.counts,"__deseq_signif_mirna", suffix), row.names=F, quote=FALSE)



# Now we prepare PCA using this data set (the data set is already scaled, so we only need to compute the transpose as we are interested 
# in the PCA for the samples)
df.counts.signif.t <- t(df.counts.signif)

# Now, create a data set from df.counts.signif.t with `strain` 
df.counts.signif.t.metadata <- df.counts.signif.t
## Extract STRAIN. Use grepl to create a logical vector whether WT appears or not in n.
strain <- ifelse(grepl('WT', colnames(df.counts.signif)), 'WT', 'A66')

# now we add the columns
df.counts.signif.t.metadata <- cbind(df.counts.signif.t.metadata, strain)




# Run PCA
scale=FALSE
pca <- prcomp(df.counts.signif.t, center=TRUE, scale.=scale)

# This provide a list of components with the respective variance (used for plotting)
eigen <- get_eig(pca)     

# Plot the variance of PCA components
png(paste0(filename.counts,"_pca_comp_variances_signif_mirna",".png"), width=2000, height=1500, res=300)
plot(pca, type = "l", main='Variance of PCA components')
dev.off()

# Show PCA summary
summary(pca)

# Write the calculated PCA `rotation` (PCA load). This is the table miRNAs (rows) vs PCAs (cols).
write.csv(pca$rotation, file=paste0(filename.counts,"__PCA_rotation_filt_deseq_signif_mirna", suffix), quote=FALSE)




# Plot PCA (colour is strain)
# PC1 vs PC2
# plot without labels
c1c2.strain <- autoplot(pca, data=df.counts.signif.t.metadata, colour='strain', x=1, y=2, scale=scale) +
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()
#ggsave(paste0(filename.counts, "_pca_strain_c1c2_deseq_signif_mirna.png"), width=4, height=4, dpi=300)

# PC2 vs PC3
# plot without labels
c2c3.strain <- autoplot(pca, data=df.counts.signif.t.metadata, colour='strain', x=2, y=3, scale=scale) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1.5) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()
#ggsave(paste0(filename.counts, "_pca_strain_c2c3_deseq_signif_mirna.png"), width=4, height=4, dpi=300)

## COMBINED plots
c1c2c3.combined <- arrangeGrob(c1c2.strain, c2c3.strain, ncol=2)
ggsave(paste0(filename.counts, "_pca_strain_c1c2c3_deseq_signif_mirna.png"), plot=c1c2c3.combined, width=7, height=7, dpi=300)





# plot with labels
c1c2.strain <- autoplot(pca, data=df.counts.signif.t.metadata, colour='strain', x=1, y=2, scale=scale, label=TRUE, label.size=2) + 
  ggtitle('PCA of samples: pc1 vs pc2') +
  coord_fixed(ratio=2.5) + 
  xlab(sprintf("PC1 %.1f %%",eigen[1,2])) + 
  ylab(sprintf("PC2 %.1f %%",eigen[2,2])) +
  theme_basic()
#ggsave(paste0(filename.counts, "_pca_strain_c1c2_deseq_signif_mirna_w_labels.png"), width=4, height=4, dpi=300)

# plot with labels
c2c3.strain <- autoplot(pca, data=df.counts.signif.t.metadata, colour='strain', x=2, y=3, scale=scale, label=TRUE, label.size=2) + 
  ggtitle('PCA of samples: pc2 vs pc3') +
  coord_fixed(ratio=1.5) + 
  xlab(sprintf("PC2 %.1f %%",eigen[2,2])) + 
  ylab(sprintf("PC3 %.1f %%",eigen[3,2])) +
  theme_basic()
#ggsave(paste0(filename.counts, "_pca_strain_c2c3_deseq_signif_mirna_w_labels.png"), width=4, height=4, dpi=300)


## COMBINED plots
c1c2c3.combined <- arrangeGrob(c1c2.strain, c2c3.strain, ncol=2)
ggsave(paste0(filename.counts, "_pca_strain_c1c2c3_deseq_signif_mirna_w_labels.png"), plot=c1c2c3.combined, width=7, height=7, dpi=300)





#############################################################################################
# "<=" Extract the most variable miRNA from PCA-PC2 rotation and see whether these match with 
# the DESeq statistically significant miRNA (for contrast : WT vs A66) 
#############################################################################################

# A threshold for selecting miRNA in PCA-PC2. 
pca.thres <- 0.05


# Create a counts table of significant miRNAs. We select all miRNA outside the interval [-pca.thres, pca.thres].
df.pca.filt <- df.pca[df.pca$PC2 < -pca.thres | df.pca$PC2 > pca.thres,]
rownames(df.pca.filt) <- df.pca.filt$miRNA
df.pca.filt <- df.pca.filt[,-1]

# Write the significant count table
#write.csv(df.pca.filt, file=paste0(filename.counts,"__PCA_rotation_filt_pc2_thres_", pca.thres, suffix), quote=FALSE)

# Write the significant miRNA
df.pca.filt.miRNA <- rownames(df.pca.filt)
write.csv(df.pca.filt.miRNA, file=paste0(filename.counts,"__PCA_PC2_mirna_thres_", pca.thres, suffix), row.names=F, quote=FALSE)




###################
# Plot Venn Diagram
###################

# Save the intersection of these two sets of miRNA
df.signif.pca.filt.miRNA <- intersect(df.counts.signif.miRNA, df.pca.filt.miRNA)
write.csv(df.signif.pca.filt.miRNA, file=paste0(filename.counts,"__Venn_Diagram_intersect_miRNA", suffix), row.names=F, quote=FALSE)

# suppress the log file venn.diagram generates
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# Now we plot the Venn diagram for these two lists of miRNA. 
venn.diagram(
  x = list(df.counts.signif.miRNA, df.pca.filt.miRNA),
  category.names = c("signif miRNA by DESeq", paste0("miRNA filt by PC2 (t=", pca.thres, ")")),
  filename = paste0(filename.counts,"__Venn_Diagram__PC2_miRNA__DESeq_miRNA.png"),
  output = FALSE,
  imagetype="png",
  height = 1000, 
  width = 1000, 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = 'blank',
  fill = c('blue', 'magenta'),
  cex = 1,
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-15, 12),
  cat.dist = c(0.050, 0.055),
  cat.fontfamily = "sans"
)

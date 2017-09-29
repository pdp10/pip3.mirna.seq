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

source('../utilities/deseq.R')
source('../utilities/plots.R')


################
# Load data sets
################

suffix <-".csv"

location.strain <- "../4_deseq__strain"
filename.deseq.strain <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66"

location.strain.signif <- "../6_clustering"
filename.deseq.strain.signif <- "miRNA__filtered_deseq_strain__padj_005_lfc_01"

location.pam <- "../6_clustering"
filename.pam <- "pam_clustering_labels"


### Load significant miRNA calculated using DESeq2:strain (padj<=0.05, |lfc|>0.1) && PCA:PC2 (>0.05). These are 89 miRNA (see Venn diagram)
# select the files containing the data
location.venn.intersect <- "../6_clustering"
filename.venn.intersect <- "venn_diagram_intersect__filt_pca_pc2_VS_signif_deseq_strain"
# load the significant miRNA names 
miRNA.venn.intersect <- rownames(read.table(paste0(location.venn.intersect,"/",filename.venn.intersect,suffix), sep=",",fill=T,header=T, row.names=1))




####################################################
# DESeq2:time(log2FoldChange) wDESeq2:Strain or wPAM
####################################################

time <- list("15", "40", "90", "180", "300")
location.time <- "../4_deseq__time"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__time__0_"

plot_deseq_lfc_tc_wstrain(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                          deseq2.strain.file=paste0(location.strain,"/",filename.deseq.strain, suffix),
                          deseq2.strain.signif.file=paste0(location.strain.signif,"/",filename.deseq.strain.signif, suffix),
                          filter.vec=miRNA.venn.intersect,
                          filename.out="miRNA_log_fc_tc__deseq_time__w_strain__venn_diag_filt")

plot_deseq_lfc_tc_wpam(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                       pam.clust.file=paste0(location.pam,"/",filename.pam, suffix),
                       filter.vec=miRNA.venn.intersect,
                       filename.out="miRNA_log_fc_tc__deseq_time__w_pam_clust__venn_diag_filt")




##################################################################
# DESeq2:strain_time_factor(log2FoldChange) wDESeq2:Strain or wPAM
##################################################################

time <- list("0", "15", "40", "90", "180", "300")
location.time <- "../4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66_at_t"

plot_deseq_lfc_tc_wstrain(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                          deseq2.strain.file=paste0(location.strain,"/",filename.deseq.strain, suffix),
                          deseq2.strain.signif.file=paste0(location.strain.signif,"/",filename.deseq.strain.signif, suffix),
                          filter.vec=miRNA.venn.intersect,
                          filename.out="miRNA_log_fc_tc__deseq_strain_time_factor__w_strain__venn_diag_filt")

plot_deseq_lfc_tc_wpam(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                       pam.clust.file=paste0(location.pam,"/",filename.pam, suffix), 
                       filter.vec=miRNA.venn.intersect,
                       filename.out="miRNA_log_fc_tc__deseq_strain_time_factor__w_pam_clust__venn_diag_filt")




#######################################################################
# DESeq2:strain_time_interaction(log2FoldChange) wDESeq2:Strain or wPAM
#######################################################################

time <- list("15", "40", "90", "180", "300")
location.time <- "../4_deseq__strain_vs_time_per_TimePoint_interaction_WaldTest"
filename.deseq.time <- "summarised_mirna_counts_after_mapping_filtered_results_DESeq__strain__WT_A66_at_t"

plot_deseq_lfc_tc_wstrain(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                          deseq2.strain.file=paste0(location.strain,"/",filename.deseq.strain, suffix),
                          deseq2.strain.signif.file=paste0(location.strain.signif,"/",filename.deseq.strain.signif, suffix), 
                          filter.vec=miRNA.venn.intersect,
                          filename.out="miRNA_log_fc_tc__deseq_strain_time_interaction__w_strain__venn_diag_filt")

plot_deseq_lfc_tc_wpam(deseq2.tc.files=paste0(location.time,"/",filename.deseq.time, time, suffix),
                       pam.clust.file=paste0(location.pam,"/",filename.pam, suffix), 
                       filter.vec=miRNA.venn.intersect,
                       filename.out="miRNA_log_fc_tc__deseq_strain_time_interaction__w_pam_clust__venn_diag_filt")


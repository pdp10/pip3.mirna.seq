
#1
print('1_prepare_counts_table')
setwd('1_prepare_counts_table')
source('prepare_counts_table.R')
rm(list = ls())
setwd('../')

#2
print('2_prepare_samples_table')
setwd('2_prepare_samples_table')
source('prepare_samples_table.R')
rm(list = ls())
setwd('../')

#3
print('3_quality_control')
setwd('3_quality_control')
source('quality_control.R')
rm(list = ls())
setwd('../')

print('3_quality_control_w_deseq2')
setwd('3_quality_control_w_deseq2')
source('quality_control_w_deseq2.R')
rm(list = ls())
setwd('../')

#4
print('4_deseq__strain')
setwd('4_deseq__strain')
source('deseq__strain.R')
rm(list = ls())
setwd('../')

print('4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain')
setwd('4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain')
source('deseq__strain_vs_time_per_TimePoint__factor_TimeStrain.R')
rm(list = ls())
setwd('../')

print('4_deseq__strain_vs_time_per_TimePoint_interaction_WaldTest')
setwd('4_deseq__strain_vs_time_per_TimePoint_interaction_WaldTest')
source('deseq__strain_vs_time_per_TimePoint_interaction_WaldTest.R')
rm(list = ls())
setwd('../')

print('4_deseq__time')
setwd('4_deseq__time')
source('deseq__time.R')
rm(list = ls())
setwd('../')

#5
print('5_clustering')
setwd('5_clustering')
source('clustering.R')
rm(list = ls())
setwd('../')

#6
print('6_compare_deseq_signif_miRNA_over_time')
setwd('6_compare_deseq_signif_miRNA_over_time')
source('compare_deseq_signif_miRNA_over_time.R')
rm(list = ls())
setwd('../')

print('6_compare_deseq_signif_miRNA_vs_pca_loading')
setwd('6_compare_deseq_signif_miRNA_vs_pca_loading')
source('compare_deseq_signif_miRNA_vs_pca_loading.R')
rm(list = ls())
setwd('../')

#7
print('7_miRNA_timecourses')
setwd('7_miRNA_timecourses')
source('miRNA_timecourses.R')
rm(list = ls())
setwd('../')


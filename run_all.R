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
print('5_compare_deseq_signif_miRNA_over_time')
setwd('5_compare_deseq_signif_miRNA_over_time')
source('compare_deseq_signif_miRNA_over_time.R')
rm(list = ls())
setwd('../')

#6
print('6_clustering')
setwd('6_clustering')
source('clustering.R')
rm(list = ls())
setwd('../')

#7
print('7_miRNA_timecourses')
setwd('7_miRNA_timecourses')
source('miRNA_timecourses.R')
rm(list = ls())
setwd('../')

#8
print('8_analysis_signif_mirna__venn_intersection')
setwd('8_analysis_signif_mirna__venn_intersection')
source('analysis_signif_mirna__venn_intersection.R')
rm(list = ls())
setwd('../')

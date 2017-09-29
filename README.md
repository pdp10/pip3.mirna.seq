# Pipeline for processing miRNA sequencing data set

#### data
The data sets.

#### 1_prepare_counts_table: 
Prepare and filter data structures

#### 2_prepare_samples_table:
Create the samples table for DESeq

#### 3_quality_control:
Run quality control (pca, cv, read density, read ecdf, matrix corr).
It also calculate the heatplot of the counts matrix.

#### 3_quality_control_w_deseq2:
Run DESeq2 dispersions estimates, mean vs sd plot, PCA analyses per time point and strain

#### 4_deseq__strain:
Run DESeq2 contrasting strain (WT,A66)

#### 4_deseq__time:
Run DESeq2 contrasting time (0,15,40,90,180,300)

#### 4_deseq__strain_vs_time_per_TimePoint__factor_TimeStrain
Run DESeq2 contrasting the strain (WT,A66) for each time point. This uses a factor (0WT,15WT, 40WT, ..., 0A66, 15A66, 40A66, ..) as design formula. 

#### 4_deseq__strain_vs_time_per_TimePoint_interaction_WaldTest:
Run DESeq2 using strain + time + strain:time as design formula. The idea is to calculate the difference between A66 and WT along the time course.
Protocol: http://www.bioconductor.org/help/workflows/rnaseqGene/
Not sure our data can be used with this approach successfully. 

#### 4_deseq__strain_vs_egf:
Run DESeq2 contrasting the strain (WT,A66) for EGF condition (EGF, noEGF). This uses a factor (EGFWT,noEGFWT, EGFA66, noEGFA66) as design formula. 

#### 5_compare_deseq_signif_miRNA_over_time:
Compare DESeq2-based significant miRNA between time points

#### 6_clustering:
Run pam clustering with labels of significant miRNA from PCA:PC2 and DESeq2 (contrast: Strain or Time).
Compare DESeq2-based significant miRNA against PCA:PC2 loadings 
greater than a certain threshold. Calculate Venn diagram, and run 
PCA for the significant sets.

#### 7_miRNA_timecourses:
Generate time courses from the counts table. miRNA are plotted all together, all together with colour representing PAM clustering, or by 
median of each cluster time courses. Data are separated between WT, A66, and A66_noEGF. 
Additionally, time courses of DESeq2-based log2 fold change (lcf) are 
plotted using the three time-dependent DESeq2 experiments.

#### 8_miRNA_timecourses__venn_intersection: 
Like 7_miRNA_timecourses but using the Venn diagram intersection calculated in 6_clustering.
It also calculate the heatplot of the counts matrix for this selected miRNAs.



### Other folders: 

#### PCA_*: 
Preliminary PCA for the counts table.

#### external: 
External code.


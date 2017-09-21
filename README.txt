
- data: 
the data sets.

- 1_prepare_counts_table: 
Essentially data filtering.

- 2_prepare_samples_table:
Create the samples table for DESeq

- 3_quality_control:
Run DESeq dispersions estimates, mean vs sd plot, PCA analyses per time point and strain

- 4_deseq__strain:
Run DESeq contrasting strain (WT,A66)

- 4_deseq__time:
Run DESeq contrasting time (0,15,40,90,180,300)

- 4_deseq__strain_vs_time_per_TimePoint:
Run DESeq contrasting the strain (WT,A66) for each time point. This uses a factor (0WT,15WT, 40WT, ..., 0A66, 15A66, 40A66, ..) as design formula. 

- 4_deseq__strain_vs_time:
Run DESeq using strain + time + strain:time as design formula. The idea is to calculate the difference between A66 and WT along the time course.
Protocol: http://www.bioconductor.org/help/workflows/rnaseqGene/
Not sure our data can be used with this approach successfully. 

- 4_deseq__strain_vs_egf:
Run DESeq contrasting the strain (WT,A66) for EGF condition (EGF, noEGF). This uses a factor (EGFWT,noEGFWT, EGFA66, noEGFA66) as design formula. 





Other folders: 

- CV contains a coefficient of variance analysis

- PCA_* folders represent preliminary PCA for the counts table.

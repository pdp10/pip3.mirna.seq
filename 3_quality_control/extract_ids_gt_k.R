
### Extract miRNA with counts greater than k
extract_ids_gt_k <- function(counts, k=400000) {
  df <- counts[apply(counts, 1, function(x) any(x > k)), ]
  write.csv(df, file=paste0("counts_gt_", k,".csv"), quote=FALSE)
}

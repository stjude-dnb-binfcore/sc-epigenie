############################################################################################################
#' Function to export cluster peaks
#' @param seurat_obj
#' @param da_peaks
#' @param results_dir
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
#' 
# compute_fold_change_peaks.R

compute_fold_change_peaks <- function(seurat_obj, da_peaks, results_dir, assay = assay) {
  
  # 1. Get unique clusters
  clusters <- unique(da_peaks$cluster)
  
  # 2. Compute fold change for each cluster
  fc_list <- clusters |> 
    set_names() |> 
    map(~ FoldChange(seurat_obj, ident.1 = .x, ident.2 = NULL, assay = assay)) |> 
    map(~ arrange(.x, desc(avg_log2FC)))
  
  # 3. Combine into one data frame with cluster labels
  fc_df_combined <- imap_dfr(fc_list, ~ mutate(.x, cluster = .y))
  
  # 4. Add rownames as a column called "feature"
  fc_df_combined <- fc_df_combined |> 
    rownames_to_column("feature") |> 
    relocate(feature)
  
  # 5. Save to TSV
  output_file <- file.path(results_dir, "fc_df_combined.tsv")
  write_tsv(fc_df_combined, file = output_file)
  message("✔ Fold change results saved to: ", output_file)
  
  return(fc_df_combined)
}
############################################################################################################

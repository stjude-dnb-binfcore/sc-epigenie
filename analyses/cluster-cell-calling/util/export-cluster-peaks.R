############################################################################################################
#' Function to format peaks
#' @param peak_name A peak string in format "chr-start-end"
#'
#' @return A formatted peak string in "chr:start-end" format
#' @export
#'
#' @examples
#' 
format_peak <- function(peak_name) {
  peak_name <- as.character(peak_name)
  parts <- unlist(strsplit(peak_name, "-"))
  paste0(parts[1], ":", parts[2], "-", parts[3])
}



############################################################################################################
#' Function to export cluster peak plots and coverage browsers
#' @param seurat_obj Seurat object with ATAC assay
#' @param fc_df_combined Data frame with fold changes and features per cluster
#' @param clusters Vector of cluster labels
#' @param top_n Number of top peaks per cluster to plot
#' @param min.cutoff_value Numeric cutoff for FeaturePlot max.cutoff
#' @param plots_dir Directory to save static plots (CoveragePlot, VlnPlot+FeaturePlot)
#' @param clustering_column description
#'
#' @return NULL (saves files to disk)
#' @export
#'
#' @examples
#' 
export_cluster_peak_plots <- function(seurat_obj,
                                      fc_df_combined,
                                      clusters = NULL,
                                      top_n,
                                      min.cutoff_value,
                                      plots_dir,
                                      clustering_column) {
  
  # Use cluster levels from Seurat object if not provided
  if (is.null(clusters)) {
    clusters <- levels(factor(seurat_obj@meta.data[[clustering_column]]))
    message(glue::glue("Using clusters from column '{clustering_column}': {paste(clusters, collapse = ', ')}"))
  }
  
  for (clust in clusters) {
    message("Processing cluster: ", clust)
    
    top_peaks <- fc_df_combined %>%
      filter(cluster == clust) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = top_n) #%>%
     # mutate(clean_peak = sub("\\.\\.\\..*", "", feature)) %>%
     # filter(clean_peak %in% rownames(seurat_obj[["peaks"]]))
    
   
    for (i in seq_len(nrow(top_peaks))) {
      peak_name <- top_peaks$feature[i]
      region <- format_peak(peak_name)
      region <- as.character(peak_name)

      # Skip if peak not found
      if (!(peak_name %in% rownames(seurat_obj[["peaks"]]))) {
        warning("Peak not found in Seurat object: ",  peak_name)
        next
      }
      
      # pip install macs3
      
      # peaks <- CallPeaks(object = seurat_obj, group.by = clusters)
      
      # Set identity to desired clustering column
      Idents(seurat_obj) <- seurat_obj@meta.data[[clustering_column]]
      
      # === Static CoveragePlot ===
      p <- CoveragePlot(seurat_obj,
                        region = peak_name,
                        #ranges = peaks,
                        #sep = "-",
                        assay = "peaks",
                        extend.upstream = 40000,
                        extend.downstream = 20000,
                        group.by = clustering_column) +
        ggtitle(glue::glue("CoveragePlot - Cluster {clust}, Peak: {peak_name}"))
      print(p)
      plot_file <- file.path(plots_dir, paste0("CoveragePlot-cluster", clust, "_peak", i, ".png"))
      ggsave(plot_file, plot = p, width = 15, height = 8)
      
      # === Violin + FeaturePlot ===
      plot1 <- VlnPlot(seurat_obj, features = peak_name, group.by =  clustering_column) +
        ggtitle(glue::glue("Violin Plot - Cluster {clust}, Peak: {peak_name}"))
      plot2 <- FeaturePlot(seurat_obj, features = peak_name, max.cutoff = min.cutoff_value) +
        ggtitle(glue::glue("Feature Plot - Peak: {peak_name}"))
      combined_plot <- plot1 | plot2
      print(combined_plot) 

      vln_feature_file <- file.path(plots_dir, paste0("VlnPlot-FeaturePlot-cluster", clust, "-peak", i, ".png"))
      ggsave(vln_feature_file, plot = combined_plot, width = 15, height = 8)
      
      # === Interactive CoverageBrowser ===
      #cb <- suppressMessages(suppressWarnings(CoverageBrowser(seurat_obj, region = peak_name, assay = "peaks")))
      #browser_file <- file.path(plots_dir, paste0("coverage_browser_cluster", clust, "_peak", i, ".html"))
      #htmlwidgets::saveWidget(cb, file = browser_file, selfcontained = TRUE)
    }
  }
  message("Done Processing cluster: ", clust)
  
}
############################################################################################################

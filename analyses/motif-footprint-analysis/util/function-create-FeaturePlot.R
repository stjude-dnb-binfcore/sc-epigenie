#############################################################################################
#' This helper generates UMAP-based feature plots for chromVAR deviation Z-scores
#' using a custom ggplot implementation that mimics Seurat's FeaturePlot color behavior.
#' It loops through a set of motifs and saves plots to disk.
#'
#' Why?
#' Bug introduced due to package conflict between Seurat 4.4 vs SeuratObject version 5.0.0 (layer)
#  Error message:
#  FeaturePlot Error: ! The `slot` argument of `FetchData()` was deprecated in SeuratObject 5.0.0 and is now defunct. i Please use the `layer` argument instead. Run `rlang::last_trace()` to see where the error occurred.
#'
#' @param seurat_obj Seurat object containing chromVAR assay and UMAP reduction.
#' @param palette Data frame of color definitions with columns `color_names` and `hex_codes`.
#'        Must include entries for "gradient_3", "gradient_8", and "na_color".
#' @param min.cutoff Minimum cutoff for feature values. Can be numeric or quantile string (e.g., "q05").
#' @param max.cutoff Maximum cutoff for feature values. Can be numeric or quantile string (e.g., "q95").
#' @param ChromVAR_plots_dir Directory path where PNG plots will be saved.
#' @param split_variable Variable to split plots (if applicable).
#'
#' @details
#' - Requires UMAP embedding stored in `seurat_obj[["umap"]]`.
#' - Assumes chromVAR deviations are in `assay = "chromvar"` and layer = "data".
#' - Motifs to plot should exist in `rownames(seurat_obj[["chromvar"]])`.
#' - Compatible with Seurat v4/v5 (uses FetchData internally).
#'
#' @return Invisibly returns NULL. Side effect: saves PNG files for each motif.
#'
#' @examples
#' # Example usage:
#' create_featureplot_overall(
#'   seurat_obj = my_obj,
#'   palette = gradient_palette_df,
#'   min.cutoff = "q05",
#'   max.cutoff = "q95",
#'   ChromVAR_plots_dir = "plots/chromvar",
#'   split_variable = "ID"
#' )
#'
create_featureplot_variable <- function(seurat_obj, palette, min.cutoff, max.cutoff, ChromVAR_plots_dir, split_variable) {
  
  
  # Build UMAP once (aligned to cell barcodes)
  emb <- Embeddings(seurat_obj[["umap"]])
  umap_df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2])
  rownames(umap_df) <- rownames(emb)
  
  # Define colors
  low_color_df  <- dplyr::filter(palette, color_names == "gradient_3")
  high_color_df <- dplyr::filter(palette, color_names == "gradient_8")
  na_color_df   <- dplyr::filter(palette, color_names == "na_color")
  
  # Small ggplot function that mimics FeaturePlot and supports faceting
  featureplot_gg <- function(df_coords, values, title, pt.size = 0.1,
                             min.cutoff = NA, max.cutoff = NA,
                             palette = palette,
                             split.by = NULL,            # <— meta col name to facet by
                             meta = NULL) {              # <— pass seurat_obj@meta.data
    
    # clip by quantiles or numeric cutoffs, if provided
    clip_quant <- function(v, cmin, cmax) {
      if (!is.na(cmin)) {
        cminv <- if (is.character(cmin) && grepl("^q\\d{2}$", cmin)) {
          stats::quantile(v, as.numeric(sub("q", "", cmin)) / 100, na.rm = TRUE)
        } else cmin
        v[v < cminv] <- cminv
      }
      if (!is.na(cmax)) {
        cmaxv <- if (is.character(cmax) && grepl("^q\\d{2}$", cmax)) {
          stats::quantile(v, as.numeric(sub("q", "", cmax)) / 100, na.rm = TRUE)
        } else cmax
        v[v > cmaxv] <- cmaxv
      }
      v
    }
    
    values <- clip_quant(values, min.cutoff, max.cutoff)
    
    # Build df and align by cell barcodes
    df <- data.frame(
      UMAP_1        = df_coords$UMAP_1,
      UMAP_2        = df_coords$UMAP_2,
      feature_value = as.numeric(values)
    )
    rownames(df) <- rownames(df_coords)
    
    # Inject meta column as "ID" for faceting
    if (!is.null(split.by)) {
      stopifnot(is.character(split.by), length(split.by) == 1)
      stopifnot(split.by %in% colnames(meta))
      # align by rownames to ensure order matches embedding
      df$ID <- meta[rownames(df), split.by, drop = TRUE]
    }
    
    p <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2, color = feature_value)) +
      ggplot2::geom_point(size = 0.1) +
      ggplot2::scale_colour_gradient(
        low      = low_color_df$hex_codes[1],
        high     = high_color_df$hex_codes[1],
        na.value = na_color_df$hex_codes[1]
      ) +
      ggplot2::labs(title = title, color = title) +
      theme_Publication(base_size = 11) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
    
    # Facet if ID exists
    if (!is.null(split.by)) {
      p <- p + ggplot2::facet_wrap(~ ID)
    }
    
    p
  }
  
  
  # ----- Loop (v4/v5-safe) -----
  motifs_to_plot <- names(top_motifs)   # or a custom vector of motif IDs
  # pick assay and fetch values
  a <- "chromvar"
  print(a)
  
  for (m in motifs_to_plot) {
    if (!m %in% rownames(seurat_obj[["chromvar"]])) {
      warning("Motif ", m, " not found in chromvar; skipping.")
      next
    }
    
    # Fetch deviations from the chromVAR assay; use layer = "data" (v5)
    vals <- FetchData(seurat_obj, vars = m, assay = "chromvar", layer = "data")[, 1]
    
    p2 <- featureplot_gg(
      df_coords  = umap_df,
      values     = vals,
      title      = m,
      pt.size    = 0.1,
      min.cutoff = min.cutoff,     # e.g., "q05"
      max.cutoff = max.cutoff,     # e.g., "q95"
      split.by   = split_variable,     # <— facet by this meta column; change as needed
      meta       = seurat_obj@meta.data
    )
    
    p <- print(p1 + p2)
    
    plot_path <- file.path(ChromVAR_plots_dir, paste0("ChromVAR_FeaturePlot_", m, "_", split_variable, ".png"))
    ggsave(plot_path, plot = p, width = 15, height = 5)
    message("Plotting ", m,  " _", split_variable, ": Complete")
  }
  
}
############################################################################################################################################




############################################################################################################################################
# IGNORE this 
# plot-RunChromVAR-overall-FeaturePlot
# This code runs outside the container.
# It fails with the container as there is a bug introduced due to package conflict between Seurat 4.4 vs SeuratObject version 5.0.0 (layer)
# Error message:
# FeaturePlot Error: ! The `slot` argument of `FetchData()` was deprecated in SeuratObject 5.0.0 and is now defunct. i Please use the `layer` argument instead. Run `rlang::last_trace()` to see where the error occurred.
# Let's keep the code for future upgrades.

#motifs_to_plot <- names(top_motifs)   # or a custom vector of motif IDs

#for (m in motifs_to_plot) {
#  if (!m %in% rownames(seurat_obj[["chromvar"]])) {
#    warning("Motif ", m, " not found in chromvar", "; skipping.")
#    next
#    }
#  message("Plotting ", m)
#  p2 <- FeaturePlot(seurat_obj,
#                    features   = m,
#                    min.cutoff = min.cutoff_value,
#                    max.cutoff = max.cutoff_value,
#                    pt.size    = 0.1) + ggplot2::ggtitle(paste0(m))
#  p <- print(p1 + p2)

#  plot_path <- file.path(ChromVAR_plots_dir, paste0("ChromVAR_FeaturePlot_overall2_", m, ".png"))
#  ggsave(plot_path, plot = p, width = 15, height = 5)
#  message("Plotting ", m, ": Complete")
#  }
############################################################################################################################################



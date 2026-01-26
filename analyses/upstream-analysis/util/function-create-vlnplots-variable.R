# --------------------------------------------------------------------------
# Create QC violin plots (facet by feature), compatible with Seurat v4/v5
# --------------------------------------------------------------------------
#' Why?
#' Bug introduced due to package conflict between Seurat 4.4 vs SeuratObject version 5.0.0 (layer)
#  Error:
# ! The `slot` argument of `FetchData()` was deprecated in SeuratObject
#  5.0.0 and is now defunct.
# i Please use the `layer` argument instead.
# For more details, see: https://github.com/stjude-dnb-binfcore/sc-epigenie/issues/106
# --------------------------------------------------------------------------
#'
#' Generates a grid of violin plots for QC metrics (or other features), optionally grouped
#' by a metadata column, with a top-left title. Designed to work in mixed Seurat/SeuratObject
#' environments by preferring meta.data and using layer-aware FetchData only when needed.
#'
#' @param seurat_obj Seurat object.
#' @param features Character vector of feature names to plot. These are expected to be
#'   columns in `seurat_obj@meta.data` (e.g., `nCount_peaks`, `TSS.enrichment`).
#'   If a feature is not in meta.data, the function will attempt to fetch it from the
#'   active assay using `Seurat::FetchData()` (tries `layer=` first, then falls back).
#' @param out_path Character file path for saving the PNG. If `NULL`, no file is saved
#'   and the plot is only returned. Default: `NULL`.
#' @param split_variable Optional single character string naming a `meta.data` column
#'   used to group samples on the x-axis (e.g., `"ID"` or `"sample"`). If `NULL`,
#'   all cells are treated as one group. Default: `NULL`.
#' @param assay Character assay name to use if a feature is not found in `meta.data`.
#'   If `NULL`, `Seurat::DefaultAssay(seurat_obj)` is used. Default: `NULL`.
#' @param sample_title Optional title placed at the top-left via `patchwork::plot_annotation()`.
#'   Default: `NULL`.
#' @param ncol Integer number of columns in the plot grid. Default: `5`.
#' @param pt.size Numeric size for jitter points overlaid on violins. Default: `0.1`.
#' @param return_plot Logical; if `TRUE`, return the patchwork object; if `FALSE`,
#'   return `invisible(NULL)`. Default: `TRUE`.
#'
#' @details
#' - This function **does not** depend on Seurat's `VlnPlot()`, so it avoids
#'   the `slot` vs `layer` incompatibility between Seurat 4.x and SeuratObject 5.x.
#' - For QC features commonly stored in `meta.data`, reading directly is fastest and safest.
#' - If a feature isn’t in `meta.data`, we try `Seurat::FetchData(..., assay=, layer=)`
#'   and fall back to `Seurat::FetchData(..., assay=)` if the first call fails.
#' - Title is left-aligned using `plot_annotation()` + a theme that sets `hjust = 0`.
# --------------------------------------------------------------------------
create_vlnplots_variable <- function(seurat_obj, features, out_path,
                                     split_variable = NULL,  # optional meta column to group along x-axis (e.g., "ID")
                                     pt.size = 0.1, ncol = NULL, sample_title = NULL) {

  stopifnot(is.character(features), length(features) > 0)

  meta  <- seurat_obj@meta.data
  cells <- colnames(seurat_obj)

  # Pull values from meta.data (preferred for QC)
  get_vec <- function(f) {
    if (!(f %in% colnames(meta))) {
      warning("Feature '", f, "' not found in meta.data; skipping.")
      return(NULL)
    }
    v <- meta[cells, f, drop = TRUE]
    names(v) <- cells
    v
  }

  make_one_violin <- function(f) {
    vals <- get_vec(f)
    if (is.null(vals)) return(NULL)
    df <- data.frame(cell = names(vals), value = as.numeric(vals))

    if (!is.null(split_variable)) {
      stopifnot(split_variable %in% colnames(meta))
      df$group <- meta[df$cell, split_variable, drop = TRUE]
    } else {
      df$group <- "SeuratProject"
    }

    ggplot(df, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width", trim = TRUE, na.rm = TRUE) +
      geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85, na.rm = TRUE) +
      geom_jitter(size = pt.size, alpha = 0.25, width = 0.15, height = 0, na.rm = TRUE) +
      labs(title = f, x = NULL, y = NULL) +
      #theme_minimal(base_size = 10) +
      theme_Publication(base_size = 11) + 
      theme(
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        axis.text.x      = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
      )
  }

  plots <- Filter(Negate(is.null), lapply(features, make_one_violin))
  if (!length(plots)) stop("No valid features to plot.")

  pw <- wrap_plots(plots, ncol = ncol)

  # Title/subtitle/caption (left-aligned title)
  if (!is.null(sample_title) || !is.null(subtitle) || !is.null(caption)) {
    pw <- pw +
      patchwork::plot_annotation(
        title   = sample_title %||% NULL #,
        #subtitle = subtitle %||% NULL,
        #caption  = caption %||% NULL
      ) &
      theme(
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold")
      )
  }

  # Save
  ggsave(filename = out_path, plot = pw, dpi = 300)
  message("Saved: ", out_path)
  return(pw)  
}

#################################################################################
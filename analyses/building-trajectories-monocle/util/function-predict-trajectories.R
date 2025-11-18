#############################################################################################
#' Helper: Identify Dimensional Reductions Associated with an Assay
#'
#' This function mimics the Seurat v4 internal function for identifying
#' all dimensional reductions associated with a given assay in a Seurat object.
#'
#' @param object A Seurat object.
#' @param assay Name of the assay (default: active assay).
#'
#' @return A character vector of dimensional reduction names associated with the assay.
#' @examples
#' AssociatedDimReducs(seurat_obj)
#'
AssociatedDimReducs <- function(object, assay = NULL) {
  assay <- assay %||% DefaultAssay(object)
  reduc.names <- names(Filter(function(x) {
    x@assay.used == assay
  }, object@reductions))
  return(reduc.names)
}


#############################################################################################
#' Predict Trajectories and Add Pseudotime to Seurat Object
#'
#' @param seurat_obj A Seurat object
#' @param cell_type_name The name of the metadata column to use for lineage assignment
#' @param lineage_value The value in cell_type_name that defines the lineage root
#' @param plots_dir Directory to save plots
#'
#' @return A list with the updated Seurat object and the Monocle3 cell_data_set
#' @export
#'
#' @examples
#' result <- predict_trajectories(seurat_obj, "predicted.id", "Rods", plots_dir = "plots/")
#' seurat_obj <- result$seurat_obj
#' cds <- result$cds
predict_trajectories <- function(seurat_obj, cell_type_name, lineage_value, plots_dir) {
  
  set.seed(2025)
  
  cat("Programmatically rooting cells in trajectory analysis\n")
  
  # Convert Seurat to Monocle3 CDS
  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  
  # Helper: programmatically select root node for the specified lineage
  get_earliest_principal_node <- function(cds, cell_type_name, lineage_value) {
    cell_ids <- which(colData(cds)[, cell_type_name] == lineage_value)
    if (length(cell_ids) == 0) {
      stop("No cells found for lineage_value in cell_type_name column.")
    }
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
      as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
    ]
    root_pr_nodes
  }
  
  root_pr_nodes <- get_earliest_principal_node(cds, cell_type_name, lineage_value)
  cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
  
  # Plot trajectory colored by pseudotime
  traj_plot_path <- file.path(plots_dir, glue::glue("Plot-trajectory-programmatic-lineage.png"))
  p_traj <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
  print(p_traj)
  ggsave(filename = traj_plot_path, plot = p_traj, width = 6, height = 5, device = "png")
  
  # Add pseudotime to Seurat object metadata
  pseudotime_vec <- cds@principal_graph_aux@listData$UMAP$pseudotime
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = pseudotime_vec,
    col.name = glue::glue("{lineage_value}_pseudotime_programmatically")
  )
  
  # Plot pseudotime on Seurat UMAP
  pseudotime_plot_path <- file.path(plots_dir, glue::glue("Plot-pseudotime.png"))
  p_pseudotime <- FeaturePlot(seurat_obj, glue::glue("{lineage_value}_pseudotime_programmatically"), pt.size = 0.1) & scale_color_viridis_c()
  print(p_pseudotime)
  ggsave(filename = pseudotime_plot_path, plot = p_pseudotime, width = 6, height = 5, device = "png")
  
  
  # Return both objects for further use
  return(list(seurat_obj = seurat_obj, cds = cds))
}


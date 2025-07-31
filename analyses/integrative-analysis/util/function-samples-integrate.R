############################################################################################################
#' Function to integrate samples using Harmony
#' @param seurat_obj
#' @param variables_to_integrate
#' @param num_dim
#' @param assay
#' @param algorithm_value description
#'
#' @return
#' @export
#'
#' @examples
#' 
#'
harmony_integration <- function(seurat_obj, variables_to_integrate, num_dim, assay, algorithm_value){
  
  set.seed(1234) # Make code reproducible
  
  cat("Performing Harmony Integration\n")
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = variables_to_integrate, reduction.use = "lsi",
                           plot_convergence = TRUE, assay.use = assay, project.dim = FALSE)


  cat("Performing Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique\n")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 2:num_dim)
  
  cat("Generate and add metadata\n")
  reduction_names <- c(paste0("umap"), paste0("lsi")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  seurat_obj@meta.data <- merge_metadata(seurat_obj, metadata) 
  write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_{integration_method}.tsv"))) # Save metadata
  
  
  cat("Computing the k.param nearest neighbors for a given dataset and clusters\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 2:num_dim) %>%
    
    # https://stuartlab.org/signac/articles/mouse_brain_vignette#non-linear-dimension-reduction-and-clustering
    FindClusters(algorithm = algorithm_value) %>% 
    identity()
  
  cat("Returned Seurat object")
  return(seurat_obj)}
  
############################################################################################################
############################################################################################################
#' To visualize all the cells together, we can co-embed the scRNA-seq and scATAC-seq cells in the same low dimensional space
#'
#' @param reference_obj A Seurat object with known cell types in `reference_obj$celltype` (scRNA-seq)
#' @param seurat_obj.query A Seurat object (scATAC-seq) to label
#' @param transfer.anchors Anchors computed from `FindTransferAnchors()`
#' @param group_by_variable Variable to embed the cells
#'
#' @return An object
#' @export
#'
#' @examples
#' coembed <- co_embedding_cells(reference_obj = rna_reference, seurat_obj.query = atac_query, transfer.anchors = transfer.anchors, group_by_variable = group_by_variable)

co_embedding_cells <- function(reference_obj, seurat_obj.query, transfer.anchors, group_by_variable) {
  
  cat("Transferring imputed RNA expression into ATAC cells...\n")
  # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
  # full transcriptome if we wanted to
  genes.use <- VariableFeatures(reference_obj)
  refdata <- GetAssayData(reference_obj, assay = "RNA", layer = "data")[genes.use, ]
  
  # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
  # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
  imputation <- TransferData(anchorset = transfer.anchors, 
                             refdata = refdata, 
                             weight.reduction = seurat_obj.query[["lsi"]], 
                             dims = 2:30)
  
  # this line adds the imputed data matrix to the seurat_obj.query object
  seurat_obj.query[["RNA"]] <- imputation
  
  # Merge and co-embed
  coembed <- merge(x = reference_obj, y = seurat_obj.query)
  
  # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  #coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
  
  # Unify cell type labels
  if ("celltype" %in% colnames(coembed@meta.data)) {
    coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
  } else {
    coembed$celltype <- coembed$predicted.id
  }
  
  name1 <- paste0(plots_dir, "/", glue::glue("plot-coembed.png"))
  p1 <- DimPlot(coembed, group.by = group_by_variable)
  p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
  print(p1 + p2)
  ggsave(file = name1, width = 18, height = 6, device = "png")
  
  
  name2 <- paste0(plots_dir, "/", glue::glue("plot-coembed-predictions.png"))
  print(DimPlot(coembed, split.by = "seq_technology_assay", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend())
  ggsave(file = name2, width = 18, height = 6, device = "png")
  
  
  return(coembed)
}

############################################################################################################

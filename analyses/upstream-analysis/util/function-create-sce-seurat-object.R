########################################################################################################
#' -- internal utility: version check
#' --- Replace the old .is_SO_v5() with this robust version ---

#suppressPackageStartupMessages({
#  library(SingleCellExperiment)
#  library(SummarizedExperiment)
#  library(Matrix)
#})

.is_SO_v5 <- function() {
  v <- utils::packageVersion("SeuratObject")
  # Convert to character like "5.0.1", split on '.', take first token, coerce to integer
  maj <- as.integer(strsplit(as.character(v), "\\.")[[1L]][1L])
  isTRUE(maj >= 5L)
}
########################################################################################################


########################################################################################################
#' -- internal utility: read a matrix by layer/slot depending on SeuratObject major version
.get_mat <- function(so, assay, layer_or_slot, prefer_layer = TRUE) {
  if (.is_SO_v5()) {
    # SeuratObject v5+: must use 'layer='
    return(Seurat::GetAssayData(so, assay = assay, layer = layer_or_slot))
  } else {
    # SeuratObject v4.x: use 'slot='; translate common layer names to matching slots
    slot_name <- switch(
      layer_or_slot,
      "counts"      = "counts",
      "data"        = "data",
      "scale.data"  = "scale.data",
      layer_or_slot # pass through (for unusual names)
    )
    return(Seurat::GetAssayData(so, assay = assay, slot = slot_name))
  }
}
########################################################################################################


########################################################################################################
#' Convert a Seurat object to SingleCellExperiment without using as.SingleCellExperiment()
#'
#' @param so            A Seurat object.
#' @param assay         Character scalar; assay name to pull matrices from. Defaults to DefaultAssay(so).
#' @param counts_layer  Character scalar; which layer to use as counts assay in SCE (e.g. "counts").
#' @param expr_layer    Character scalar; optional expression/normalized layer to add as 'logcounts'
#'                       (e.g. "data"). Set NULL to skip.
#' @param add_reductions Logical; copy Seurat reductions into SCE reducedDims. Default TRUE.
#'
#' @return A SingleCellExperiment with assays(counts=...), optional logcounts, colData, rowData, reducedDims.
seurat_to_sce_safe <- function(
    so,
    assay         = NULL,
    counts_layer  = "counts",
    expr_layer    = "data",
    add_reductions = TRUE) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(so)
  
  # Validate that requested layers exist on the assay when on v5
  if (.is_SO_v5()) {
    lyr_avail <- tryCatch(SeuratObject::Layers(so[[assay]]), error = function(e) character(0))
    if (!(counts_layer %in% lyr_avail)) {
      stop(sprintf("Layer '%s' not found in assay '%s'. Available: %s",
                   counts_layer, assay, paste(lyr_avail, collapse = ", ")), call. = FALSE)
    }
    if (!is.null(expr_layer) && !(expr_layer %in% lyr_avail)) {
      warning(sprintf("Expr layer '%s' not found; will skip logcounts. Available: %s",
                      expr_layer, paste(lyr_avail, collapse = ", ")))
      expr_layer <- NULL
    }
  }
  
  # ---- counts matrix (required by scDblFinder)
  m_counts <- .get_mat(so, assay = assay, layer_or_slot = counts_layer)
  if (!inherits(m_counts, "dgCMatrix")) m_counts <- as(m_counts, "dgCMatrix")
  
  # ---- colData from meta.data (align to cell order)
  cells <- colnames(so)
  meta  <- so@meta.data
  cd <- DataFrame(meta[cells, , drop = FALSE])
  
  # ---- rowData (lightweight)
  rd <- DataFrame(row.names = rownames(m_counts))
  rd$feature <- rownames(m_counts)
  
  sce <- SingleCellExperiment(
    assays  = list(counts = m_counts),
    colData = cd,
    rowData = rd
  )
  
  # ---- optional expression layer -> logcounts
  if (!is.null(expr_layer)) {
    m_expr <- .get_mat(so, assay = assay, layer_or_slot = expr_layer)
    if (!inherits(m_expr, "dgCMatrix")) m_expr <- as(m_expr, "dgCMatrix")
    # keep shapes aligned to cells/rows
    m_expr <- m_expr[rownames(sce), colnames(sce), drop = FALSE]
    SummarizedExperiment::assay(sce, "logcounts") <- m_expr
  }
  
  # ---- reducedDims (optional)
  if (isTRUE(add_reductions) && length(so@reductions)) {
    rd_list <- list()
    for (red in names(so@reductions)) {
      emb <- try(Seurat::Embeddings(so[[red]]), silent = TRUE)
      if (!inherits(emb, "try-error")) {
        # Embeddings rows are cells; align to SCE columns
        common <- intersect(rownames(emb), colnames(sce))
        if (length(common)) {
          rd_list[[red]] <- as.matrix(emb[common, , drop = FALSE])[colnames(sce), , drop = FALSE]
        }
      }
    }
    if (length(rd_list)) {
      reducedDims(sce) <- rd_list
    }
  }
  
  sce
}
########################################################################################################



########################################################################################################
# --- dependencies
#suppressPackageStartupMessages({
#  library(Seurat)
#  library(SeuratObject)
#  library(SingleCellExperiment)
#  library(SummarizedExperiment)
#  library(Matrix)
#})

# --- set "data" (or any name) into an Assay using the right API for v4/v5
.set_data_matrix <- function(assay_obj, mat, name = "data") {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  if (.so_is_v5()) {
    # v5+: write into a layer
    assay_obj <- SeuratObject::SetAssayData(object = assay_obj, layer = name, new.data = mat)
  } else {
    # v4.x: write into a slot
    assay_obj <- SeuratObject::SetAssayData(object = assay_obj, slot  = name, new.data = mat)
  }
  assay_obj
}
########################################################################################################


########################################################################################################
#' Convert SingleCellExperiment to Seurat without using as.Seurat()
#'
#' @param sce         A SingleCellExperiment with at least a "counts" assay.
#' @param assay_name  Name for the Seurat assay to create (e.g., "peaks", "RNA").
#' @param data_from   SCE assay to copy into Seurat "data" (default "logcounts").
#'                    Use NULL to skip; if not found, falls back to log1p(counts).
#' @param add_reductions  If TRUE, copy reducedDims to Seurat reductions.
#' @return A Seurat object.
sce_to_seurat_safe <- function(sce, assay_name = "RNA", data_from = "logcounts", add_reductions = TRUE) {
  stopifnot(inherits(sce, "SingleCellExperiment"))
  
  ## 1) counts (required)
  if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
    stop("SCE must contain a 'counts' assay.", call. = FALSE)
  }
  counts_mat <- SummarizedExperiment::assay(sce, "counts")
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(counts_mat, "dgCMatrix")
  
  ## 2) create base Seurat object (this places counts in the assay already)
  so <- Seurat::CreateSeuratObject(
    counts   = counts_mat,
    assay    = assay_name,
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))
  )
  
  ## 3) optional "data" matrix
  if (!is.null(data_from)) {
    if (data_from %in% SummarizedExperiment::assayNames(sce)) {
      data_mat <- SummarizedExperiment::assay(sce, data_from)
    } else {
      # fallback: compute a simple log1p of counts if requested layer not present
      message("Assay '", data_from, "' not found in SCE; using log1p(counts) as 'data'.")
      data_mat <- as.matrix(Matrix::t(Matrix::t(counts_mat))) # keep sparse-friendly; we’ll log1p below
      data_mat <- as(data_mat, "dgCMatrix")
      data_mat@x <- log1p(data_mat@x)
    }
    # align to Seurat feature/cell order
    data_mat <- data_mat[rownames(so), colnames(so), drop = FALSE]
    # write to the assay using version-aware setter
    so[[assay_name]] <- .set_data_matrix(so[[assay_name]], data_mat, name = "data")
  }
  
  ## 4) reduced dimensions
  if (isTRUE(add_reductions) && length(SingleCellExperiment::reducedDims(sce)) > 0) {
    for (nm in names(SingleCellExperiment::reducedDims(sce))) {
      emb <- SingleCellExperiment::reducedDims(sce)[[nm]]
      # ensure rownames exist and align to cells
      if (is.null(rownames(emb))) {
        # try to infer from colnames(sce) length; if mismatch, skip
        if (nrow(emb) == ncol(so)) {
          rownames(emb) <- colnames(so)
        } else {
          warning("ReducedDim '", nm, "' has no rownames; skipping.")
          next
        }
      }
      common <- intersect(rownames(emb), colnames(so))
      if (!length(common)) next
      emb <- emb[colnames(so), , drop = FALSE]  # reorder to Seurat cells
      so[[nm]] <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(emb),
        key        = paste0(toupper(nm), "_"),
        assay      = assay_name
      )
    }
  }
  
  so
}
########################################################################################################


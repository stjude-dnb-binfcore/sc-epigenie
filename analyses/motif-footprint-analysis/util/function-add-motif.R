#########################################################################################################
#' Helper: Get a list of motif position frequency matrices from the JASPAR database
#' 
#' Get JASPAR motif position frequency matrices (PFMs)
#'
#' Retrieve a list of motif position frequency matrices (PFMs) from the
#' JASPAR database, typically used to build motif objects for downstream
#' analyses (e.g., chromatin accessibility or motif enrichment) in single-cell
#' workflows. Accepts a Seurat object in case you want to derive species or
#' genome context from metadata, though it is not strictly required.
#'
#' @param 
#'
#' @return A `TFBSTools::PFMatrixList` containing JASPAR motif PFMs.
#' @examples
#'
if (jaspar_library_version == "2018") {
  pfm <- getMatrixSet(x = JASPAR2018,
                      opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
  
  } else if (jaspar_library_version == "2020") {
    pfm <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
  
  } else if (jaspar_library_version == "2022") {
    pfm <- getMatrixSet(x = JASPAR2022,
                        opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
  
  } else if (jaspar_library_version == "2024") {
  # JASPAR2024 package does not seem to be supported by the TFBSTools package, which also causes an error when running getMatrixSet():
  # Error: unable to find an inherited method for function ‘getMatrixSet’ for signature ‘x = "function"’
  # The following resources indicate a good alternative solution to get the so desirable sweet peaks
  # https://github.com/stuart-lab/signac/discussions/1646
  # https://github.com/ge11232002/TFBSTools/issues/39
    
    jaspar <- JASPAR2024::JASPAR2024()
    pfm <- getMatrixSet(x = jaspar@db,
                        opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
}

#########################################################################################################


#########################################################################################################
#
#' Add motif analysis to a Seurat object (Signac)
#'
#' Harmonize peak ranges and counts for a specified ChromatinAssay, reconstruct
#' the assay to ensure strict alignment, and then add motif information using
#' a \pkg{BSgenome} reference and a list of JASPAR PFMs (via \pkg{TFBSTools}).
#' This helper standardizes seqlevels to the UCSC style, keeps only chromosomes
#' present in both the assay ranges and the genome, recreates the ChromatinAssay
#' with aligned ranges and counts, and finally calls \code{Signac::AddMotifs}.
#'
#' @param seurat_obj A \code{Seurat} object containing at least one
#'   \code{ChromatinAssay}. The specified \code{assay} must inherit
#'   \code{"ChromatinAssay"}.
#' @param assay A length-1 character string naming the ChromatinAssay to update.
#'   This function does not infer a default; you must specify the assay name.
#' @param genome A \code{BSgenome} object (e.g.,
#'   \code{BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}).
#'   A string genome name is not accepted.
#' @param pfm A \code{TFBSTools::PFMatrixList} of JASPAR motif position
#'   frequency matrices (PFMs), typically obtained via
#'   \code{TFBSTools::getMatrixSet}.
#'
#' @return The updated \code{Seurat} object with motifs added to the specified
#'   ChromatinAssay. Motifs can be accessed with \code{Signac::Motifs()} and
#'   used in downstream motif enrichment and accessibility analyses.
#'

add_motifs <- function(seurat_obj, assay, genome, pfm) {
  
  ############################################################
  # (1)  Identify the chromatin assay and get its ranges & counts
  # Find the first assay that is a ChromatinAssay (adjust if you have more than one)
  chrom_assays <- names(seurat_obj)[sapply(names(seurat_obj), function(a) inherits(seurat_obj[[a]], "ChromatinAssay"))]
  
  # Extract current ranges and counts
  old_atac <- seurat_obj[[assay]]
  gr       <- granges(old_atac)
  counts   <- GetAssayData(old_atac, slot = "counts")
  sep <- c("-", "-")
  
  ############################################################
  # (2) Harmonize seqlevels style and keep only chromosomes present in both
  # UCSC style on both objects
  seqlevelsStyle(gr)     <- "UCSC"
  seqlevelsStyle(genome) <- "UCSC"
  
  # Drop seqlevels not present in both; optionally keep only standard chromosomes
  common <- intersect(seqlevels(gr), seqlevels(genome))
  gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  
  # Create canonical peak name strings from GRanges using the detected separator style
  peak.names <- GRangesToString(gr, sep = sep)        # e.g., "chr1:100-200" or "chr1-100-200"
  # Ensure uniqueness (rare but can happen after pruning)
  if (anyDuplicated(peak.names) > 0) peak.names <- make.unique(peak.names)
  names(gr) <- peak.names
  
  ############################################################
  # (3) Subset & reorder the counts to exactly match gr
  # Peaks present in both objects, preserving counts' order
  shared <- rownames(counts)[ rownames(counts) %in% names(gr) ]
  
  # Build aligned objects
  counts_new <- counts[shared, , drop = FALSE]       # same order as 'shared'
  idx        <- match(shared, names(gr))             # numeric indices into gr
  gr_new     <- gr[idx]
  
  # Final alignment check
  stopifnot(identical(rownames(counts_new), names(gr_new)))
  
  ############################################################
  # (4) Recreate the ChromatinAssay (safer than slot editing)
  frag <- try(Fragments(old_atac), silent = TRUE)
  
  new_atac <- CreateChromatinAssay(counts = counts_new, ranges = gr_new, fragments = if (inherits(frag, "try-error")) NULL else frag)
  seurat_obj[[assay]] <- new_atac
  DefaultAssay(seurat_obj) <- assay
  
  ############################################################
  # (5) Run AddMotifs with the genome object)
  seurat_obj <- AddMotifs(object = seurat_obj,
                          genome = genome,  # BSgenome object, not a string
                          pfm    = pfm)
  
  ############################################################
  # (6) Sanity check
  # Retrieve the motif object
  motif_obj <- Motifs(seurat_obj)
  motif_obj
  
  # Inspect motif IDs and names
  head(rownames(motif_obj))        # motif IDs
  head(motif_obj@motif.names)      # motif names
  
  # Inspect motif-to-peak mapping
  head(motif_obj@data)  # sparse matrix: peaks × motifs
  
  # Save object
  saveRDS(motif_obj, file = paste0(results_dir, "/", "motif_obj.rds")) 
  
  #################################################################
  #return(seurat_obj)
  return(list(seurat_obj = seurat_obj, motif_obj = motif_obj))
}
#########################################################################################################


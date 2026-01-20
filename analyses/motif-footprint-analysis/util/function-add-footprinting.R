#############################################################################################
#' Get JASPAR motif position weight matrices (PWMs)
#'
#' Retrieve motif position **frequency** matrices (PFMs) from the JASPAR database
#' and convert them to **position weight matrices** (PWMs) for downstream scoring,
#' footprinting, and motif-centric analyses in single-cell workflows.
#'
#' @param species Character scalar specifying the species name (e.g., "Homo sapiens"),
#'   species common name, or numeric taxonomy id supported by JASPAR. You may also
#'   pass a vector (e.g., c("Homo sapiens", "Mus musculus")) to combine species.
#'   For broad sets, you can use \code{tax_group = "vertebrates"} via \code{opts}.
#'
#' @return A \code{TFBSTools::PWMatrixList} containing JASPAR motif PWMs.
#' @examples
#'
# extract position frequency matrices for the motifs
pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = species, all_versions = FALSE)) 



#############################################################################################
#

#' Add transcription factor footprint analysis to a Seurat object (Signac)
#'
#' Perform motif-centric footprinting by aggregating Tn5 insertion signal
#' around motif instances. This helper assumes motifs have already been added
#' to the target \code{ChromatinAssay} (e.g., via \code{Signac::AddMotifs})
#' and uses a \pkg{BSgenome} reference together with motif matrices
#' (PFM/PWM) to compute footprint profiles for selected motifs.
#'
#' @param seurat_obj A \code{Seurat} object containing at least one
#'   \code{ChromatinAssay}. The specified \code{assay} must inherit
#'   \code{"ChromatinAssay"} and have motif annotations present
#'   (see \code{Signac::AddMotifs}).
#' @param assay A length-1 character string naming the ChromatinAssay to update.
#'   This function does not infer a default; you must specify the assay name.
#' @param genome A \code{BSgenome} object (e.g.,
#'   \code{BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}).
#'   A string genome name is not accepted.
#' @param pwn A set of motif matrices used for footprinting. **Note:** This
#'   appears to be a typo for "PWM". Supply either a \code{TFBSTools::PFMatrixList}
#'   (PFMs) or \code{TFBSTools::PWMatrixList} (PWMs) compatible with Signac’s
#'   footprint functions. Prefer PWMs for scoring; PFMs can be converted.
#' @param motifs_to_fp A character vector of motif names/IDs to footprint
#'   (must match motif identifiers present in \code{Motifs(seurat_obj[[assay]])}).
#' @param cell_type_name A metadata column containing the cell type name.
#' @param top_n_value_footprinting Number of top motifs to footprint.

#'
#' @return The updated \code{Seurat} object with footprint data stored in the
#'   assay's motifs/feature-level metadata, suitable for downstream visualization
#'   (e.g., \code{Signac::PlotFootprint}) and comparative accessibility analyses.
#'
#'

add_footprinting <- function(seurat_obj, assay, genome, pfm, cell_type_name, top_n_value_footprinting) {
  
  ########################################################################################
  # (Step 1): Add motif information by using JASPAR motif position weight matrices (PWMs)
  # DefaultAssay(seurat_obj) <- assay

  # add motif information
  seurat_obj <- AddMotifs(seurat_obj, genome = genome, pfm = pwm)

  ############################################################
  # Sanity check
  # Retrieve the motif object
  motif_obj <- Motifs(seurat_obj)
  motif_obj
  
  # Inspect motif IDs and names
  head(rownames(motif_obj))        # motif IDs
  head(motif_obj@motif.names)      # motif names
  
  # Inspect motif-to-peak mapping
  head(motif_obj@data)  # sparse matrix: peaks × motifs
  
  ########################################################################################
  # (Step 2): Choose motif information 
  motif_names <- motif_obj@motif.names          # readable names like "CEBPB", "SPI1", etc.
  motif_ids   <- rownames(motif_obj)            # IDs like "MA0062.1"

  # Choose the motifs to footprint (examples)
  motifs_to_fp <- head(motif_names, top_n_value_footprinting)          # or use IDs: head(motif_ids, 6)

  #################################################################
  # OPTIONAL: if you want the top enriched motifs instead of arbitrary top 6
  #best6 <- enriched.motifs %>%
  #   arrange(desc(fold.enrichment), desc(percent.observed)) %>%
  #   slice_head(n = 6)
  # # map enriched IDs to motif_obj (tolerate version mismatches like ".1" vs ".2")
  #library(stringr)
  #requested_core <- str_replace(as.character(best6$motif), "\\.\\d+$", "")
  #avail_core     <- str_replace(motif_ids, "\\.\\d+$", "")
  #idx            <- match(requested_core, avail_core)
  #motifs_to_fp   <- motif_names[!is.na(idx)]  # names matched to IDs present
  #################################################################

  ########################################################################################
  # (Step 3): We need to attach the fragments files to the seurat object before any Motif footprinting steps.
  
  # Get sample name
  metadata_dir <- yaml$metadata_dir
  metadata_file <- yaml$metadata_file
  project_metadata_file <- file.path(metadata_dir, metadata_file) # metadata input file
  
  project_metadata <- read.csv(project_metadata_file, sep = "\t", header = TRUE)
  sample_name <- unique(as.character(project_metadata$ID))
  sample_name <- sort(sample_name, decreasing = FALSE)
  print(sample_name)
  
  # create-fragment-list
  # Attach appropriate fragments files to use for plotting later
  fragment_paths <- file.path(yaml$data_dir, sample_name, "outs", "fragments.tsv.gz")
  
  fragment_list <- mapply(function(sample, frag_path) {
    if (!file.exists(frag_path)) {
      warning("Fragment file not found for sample: ", sample)
      return(NULL)
    }
    matching_cells <- paste0(
      str_remove(grep(paste0("^", sample, ":"), colnames(seurat_obj), value = TRUE), pattern = paste0(sample, ":")),
      "-1"
    )
    names(matching_cells) <- grep(paste0("^", sample, ":"), colnames(seurat_obj), value = TRUE)
    
    tryCatch({
      CreateFragmentObject(path = frag_path, cells = matching_cells)
    }, error = function(e) {
      warning("Failed to create fragment object for sample: ", sample, " (", e$message, ")")
      return(NULL)
    })
  }, sample_name, fragment_paths, SIMPLIFY = FALSE)
  
  # Remove any NULL entries in case of failure
  fragment_list <- Filter(Negate(is.null), fragment_list)
  seurat_obj[[assay]]@fragments <- fragment_list
  # Fragments(seurat_obj)
  
  ########################################################################################
  # (Step 4): Compute footprints for the selected motifs ONCE ---
  seurat_obj <- Footprint(object     = seurat_obj,
                          motif.name = motifs_to_fp,     # *motif names or IDs*, NOT peak coordinates
                          genome     = genome)

  #################################################################
  return(list(seurat_obj = seurat_obj, motif_obj = motif_obj, motifs_to_fp = motifs_to_fp))
  }
#################################################################


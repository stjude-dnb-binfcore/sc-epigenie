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
if (jaspar_library_version == "2018") {
  pwm <- getMatrixSet(x = JASPAR2018,
                      opts = list(species = species, all_versions = FALSE))
  
  } else if (jaspar_library_version == "2020") {
    pwm <- getMatrixSet(x = JASPAR2020,
                        opts = list(species = species, all_versions = FALSE))
    
  } else if (jaspar_library_version == "2022") {
    pwm <- getMatrixSet(x = JASPAR2022,
                        opts = list(species = species, all_versions = FALSE))
    
  } else if (jaspar_library_version == "2024") {
    # JASPAR2024 package does not seem to be supported by the TFBSTools package, which also causes an error when running getMatrixSet():
    # Error: unable to find an inherited method for function ‘getMatrixSet’ for signature ‘x = "function"’
    # The following resources indicate a good alternative solution to get the so desirable sweet peaks
    # https://github.com/stuart-lab/signac/discussions/1646
    # https://github.com/ge11232002/TFBSTools/issues/39
    jaspar <- JASPAR2024::JASPAR2024()
    pwm <- getMatrixSet(x = jaspar@db,
                        opts = list(species = species, all_versions = FALSE))
}




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
#' @param motifs_to_fp_module_list

#'
#' @return The updated \code{Seurat} object with footprint data stored in the
#'   assay's motifs/feature-level metadata, suitable for downstream visualization
#'   (e.g., \code{Signac::PlotFootprint}) and comparative accessibility analyses.
#'
#'

add_footprinting <- function(seurat_obj, assay, genome, pfm, cell_type_name, top_n_value_footprinting, motifs_to_fp_module_list) {
  
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
  
  if (is.null(motifs_to_fp_module_list)) {
    message("Motifs for footprinting are being estimated based on the Z-scores from the chromVAR deviation matrix and we will use the top ", top_n_value_footprinting, " motifs", "\n")
 
    ######################################
    # Print best motifs for footprinting
    chromvar_mat <- seurat_obj[["chromvar"]]@data # Get the chromVAR deviation matrix
    avg_scores <- rowMeans(chromvar_mat) # Compute average deviation Z-score per motif across all cells
  
    # Sort and pick the top motifs
    # motifs_to_fp <- sort(avg_scores, decreasing = TRUE)[1:top_n_value_footprinting] 
    motifs_to_fp_scores <- sort(avg_scores, decreasing = TRUE)
  
    # Sanity check that these motifs exist in the motif_obj as well.
    motif_ids_in_obj <- sort(names(motif_names))            # "MA0004.1", "MA0006.1", ...
    motif_ids_fp     <- sort(names(motifs_to_fp_scores))           # "MA0060.3", "MA1644.1", ...
    overlap_ids <- intersect(motif_ids_fp, motif_ids_in_obj)
    #print(overlap_ids)

    # Choose the motifs to footprint (examples)
    #top_overlap_ids <- head(overlap_ids, top_n_value_footprinting)          # or use IDs: head(motif_ids, 6)
    #motifs_to_fp <- motif_names[top_overlap_ids]
    
    # Subset motifs scores to overlapping ids
    motifs_to_fp_scores_overlap <- motifs_to_fp_scores[overlap_ids]
    
    # Sort by average z-scores again
    motifs_to_fp_scores_overlap <- sort(motifs_to_fp_scores_overlap, decreasing = TRUE)
    print(head(motifs_to_fp_scores_overlap))
    # Get top n motifs to plot
    motifs_to_fp <- motif_names[head(names(motifs_to_fp_scores_overlap), top_n_value_footprinting)]

    } else if (length(motifs_to_fp_module_list) > 1) {
      message("Motifs for footprinting are being provided by the user.", "\n")
      motifs_to_fp_names <- as.character(motifs_to_fp_module_list)
      motifs_to_fp <- motif_names[motifs_to_fp_names]
    }
  
  print(motifs_to_fp)
  
  ######################################

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


############################################################################################################
#' Process and re-index scATAC-seq fragment files if needed
#'
#' This function renames cell barcodes in 10X Genomics fragment files to ensure compatibility
#' with Seurat and Signac functions (e.g., CallPeaks). It also re-indexes the files with Tabix
#' if required.
#'
#' @param yaml A list containing configuration values, including:
#'   \itemize{
#'     \item{\code{metadata_dir}}{ - Path to the directory containing metadata}
#'     \item{\code{metadata_file}}{ - File name of the metadata (TSV format)}
#'     \item{\code{data_dir}}{ - Path to the base directory with fragment files}
#'   }
#'
#' @return NULL (invisible). The function updates the fragment files in-place.
#' @export
#'
#' @examples
#' yaml <- list(
#'   metadata_dir = "data/metadata/",
#'   metadata_file = "project_metadata.tsv",
#'   data_dir = "data/processed/"
#' )
#' process_fragments(yaml)
#'
############################################################################################################

process_fragments <- function(yaml) {
  metadata_dir <- yaml$metadata_dir
  metadata_file <- yaml$metadata_file
  project_metadata_file <- file.path(metadata_dir, metadata_file)
  
  # Load metadata and extract sample names
  project_metadata <- read.csv(project_metadata_file, sep = "\t", header = TRUE)
  sample_name <- sort(unique(as.character(project_metadata$ID)), decreasing = FALSE)
  print(sample_name)
  
  for (sample in sample_name) {
    input_path <- file.path(yaml$data_dir, sample, "outs")
    fragments_gz <- file.path(input_path, "fragments.tsv.gz")
    index_file <- paste0(fragments_gz, ".tbi")
    backup_fragments_gz <- file.path(input_path, "fragments_original.tsv.gz")
    
    # Check if processed fragments and index exist
    if (file.exists(fragments_gz) && file.exists(index_file)) {
      message("Processed fragments and index exist for sample: ", sample, ". Skipping to next sample.")
      next
    }
    
    # If not, check for original backup and process if needed
    if (!file.exists(backup_fragments_gz)) {
      # Backup the original fragments file
      if (file.exists(fragments_gz)) {
        file.copy(fragments_gz, backup_fragments_gz)
        message("Backup created for sample: ", sample)
      } else {
        message("No fragments.tsv.gz found for sample: ", sample, ". Skipping.")
        next
      }
    } else {
      message("Backup already exists for sample: ", sample)
    }
    
    # Unzip, modify barcodes, save, and re-compress
    fragments_tmp <- file.path(input_path, "fragments.tsv")
    R.utils::gunzip(backup_fragments_gz, destname = fragments_tmp, remove = FALSE, overwrite = TRUE)
    dt <- data.table::fread(fragments_tmp, header = FALSE, sep = "\t")
    dt[[4]] <- paste0(sample, ":", gsub("-1$", "", dt[[4]]))  # Modify barcode column
    dt <- dt[, 1:5]  # Ensure only first 5 columns
    data.table::fwrite(dt, fragments_tmp, sep = "\t", col.names = FALSE)
    R.utils::gzip(fragments_tmp, destname = fragments_gz, overwrite = TRUE)
    message("Processed fragments for sample: ", sample)
    
    # Re-index the processed fragments file
    Rsamtools::bgzip(fragments_tmp, dest = fragments_gz, overwrite = TRUE)
    Rsamtools::indexTabix(fragments_gz, format = "bed")
    message("Re-indexed fragments for sample: ", sample)
  }
  
  message("Fragment processing and indexing complete.")
}

############################################################################################################


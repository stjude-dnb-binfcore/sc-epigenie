#################################################################################
# This will run the `01A-run-signac-qc.Rmd` script for multiple samples/libraries
# and save html report separetely for each one of them
# https://pkgs.rstudio.com/rmarkdown/reference/render.html
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(knitr)
  library(tidyverse)
  library(yaml)
  library(optparse)
})

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)

#################################################################################
# Set up directories and paths to file Inputs/Outputs
root_dir <- yaml$root_dir
metadata_dir <- yaml$metadata_dir
metadata_file = yaml$metadata_file
analysis_dir <- file.path(root_dir, "analyses", "upstream-analysis") 

# File path to plots directory for signac_qc
signac_qc_plots_dir <-
  file.path(plots_dir, "01_Signac_qc") 
if (!dir.exists(signac_qc_plots_dir)) {
  dir.create(signac_qc_plots_dir)
}

# Create signac_results_dir
signac_results_dir <- 
  file.path(module_results_dir, paste0("01_Signac_qc"))
if (!dir.exists(signac_results_dir)) {
  dir.create(signac_results_dir)
}

#######################################################
# Read metadata file and define `sample_name`
project_metadata_file <- file.path(metadata_dir, metadata_file) # metadata input file

# Read metadata file and define `sample_name`
project_metadata <- read.csv(project_metadata_file, sep = "\t", header = TRUE)
sample_name <- unique(as.character(project_metadata$ID))
sample_name <- sort(sample_name, decreasing = FALSE)
print(sample_name)

#####################################################################################
# Run markdown script per each library
for (i in seq_along(sample_name)){
  
  cat("Beginning to process sample:", sample_name[i], "\n")
  
  # Create directory to save html reports
  plots_dir <- file.path(signac_qc_plots_dir, sample_name[i]) 
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)}
    
  # Create results_dir per sample
  results_dir <- file.path(signac_results_dir, sample_name[i])
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)}
  
  # Render and save html
  rmarkdown::render("01A_run_signac_qc.Rmd", 
                    output_dir = file.path(plots_dir),
                    clean = TRUE, # using TRUE will clean intermediate files that are created during rendering
                    output_file = c(paste('Report-', 'signac-qc', '-', sample_name[i], '-', Sys.Date(), sep = '')),
                    output_format = 'all',
                    params = list(data_dir = yaml$data_dir,
                                  genome_name = yaml$genome_name_upstream,
                                  assay = yaml$assay_signac_qc,
                                  #annotation_Ensembl = "EnsDb.Mmusculus.v79",        
                                  annotation_Ensembl = yaml$annotation_Ensembl_upstream,
                                  species = yaml$species_upstream,
                                  object_species = yaml$object_species_upstream,
                                  min.cutoff_value = yaml$min.cutoff_value_upstream,
                                  #nCount_peaks_min = yaml$nCount_peaks_min_upstream,
                                  #nCount_peaks_max = yaml$nCount_peaks_max_upstream,
                                  #nFeature_peaks_min = yaml$nFeature_peaks_min_upstream,
                                  #nFeature_peaks_max = yaml$nFeature_peaks_max_upstream,
                                  #peak_region_fragments_min = yaml$peak_region_fragments_min_upstream,
                                  #peak_region_fragments_max = yaml$peak_region_fragments_max_upstream,
                                  pct_reads_in_peaks_value = yaml$pct_reads_in_peaks_value_upstream,
                                  blacklist_ratio_value = yaml$blacklist_ratio_value_upstream,
                                  #duplicate_value = yaml$duplicate_value_upstream,
                                  #mitochondrial_value = yaml$mitochondrial_value_upstream,
                                  TSS.enrichment_value = yaml$TSS.enrichment_value_upstream,
                                  nucleosome_signal_value = yaml$nucleosome_signal_value_upstream,
                                  use_threshold_filtering = yaml$use_threshold_filtering_upstream,
                                  condition_value1 = yaml$condition_value1,
                                  condition_value2 = yaml$condition_value2,
                                  condition_value3 = yaml$condition_value3,
                                  use_condition_split = yaml$use_condition_split_seurat_multiple_samples,
                                  print_pdf = yaml$print_pdf_seurat_multiple_samples,
                                  grouping = yaml$grouping,

                                  # the following parameters are the same across the module #
                                  PROJECT_NAME = yaml$PROJECT_NAME,
                                  PI_NAME = yaml$PI_NAME,
                                  TASK_ID = yaml$TASK_ID,
                                  PROJECT_LEAD_NAME = yaml$PROJECT_LEAD_NAME,
                                  DEPARTMENT = yaml$DEPARTMENT,
                                  LEAD_ANALYSTS = yaml$LEAD_ANALYSTS,
                                  GROUP_LEAD = yaml$GROUP_LEAD,
                                  CONTACT_EMAIL = yaml$CONTACT_EMAIL,
                                  PIPELINE = yaml$PIPELINE, 
                                  START_DATE = yaml$START_DATE,
                                  COMPLETION_DATE = yaml$COMPLETION_DATE))
}

#################################################################################

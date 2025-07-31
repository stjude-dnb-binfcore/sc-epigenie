#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
#.libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)})

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)

#################################################################################
# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "upstream-analysis") 
filename_filter_object = yaml$filename_filter_object_value
filename_summary_report = yaml$filename_summary_report_value

# File path to plots directory
plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Create module_results_dir
module_results_dir <- file.path(analysis_dir, "results")
if (!dir.exists(module_results_dir)) {
  dir.create(module_results_dir)
}


# Get the suffix or empty string if NULL
suffix <- if (!is.null(filename_filter_object)) {
  glue::glue("_{filename_filter_object}")
} else {
  ""
}

Signac_dir <- file.path(analysis_dir, "plots", "01_Signac_qc") 
scDblFinder_dir <- file.path(analysis_dir, "plots", "02_scDblFinder") 
Filter_object_dir <- file.path(analysis_dir, "plots", glue::glue('03_Filter_object{suffix}')) 
Final_summary_report_dir <- file.path(analysis_dir, "plots", glue::glue('04_Final_summary_report{suffix}'))

################################################################################################################
# Run Rmd scripts to process data per method
################################################################################################################
future_globals_value = 214748364800 #200 * 1024^3; # 150 * 1024^3; other options: 1000 * 1024^2 = 1048576000; 8000 * 1024^2 =8388608000
################################################################################################################

###############################################################################################################
# (1) Signac QC metrics
# Run the seurat_qc script for each sample/library and save html/pdf reports per each
source(paste0(analysis_dir, "/", "01B_run_signac_qc_multiple_samples.R"))

###############################################################################################################
# (2) Estimating and filtering out doublets
rmarkdown::render('02_run_scDblFinder.Rmd', 
                  clean = TRUE,
                  output_dir = file.path(scDblFinder_dir),
                  output_file = paste('Report-', 'scDblFinder', '-', Sys.Date(), sep = ''),
                  output_format = 'all',
                  params = list(assay = yaml$assay_signac_qc,
                    min.cutoff_value = yaml$min.cutoff_value_upstream,
                    root_dir = yaml$root_dir,
                    metadata_dir = yaml$metadata_dir,
                    metadata_file = yaml$metadata_file,
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


##############################################################################################################
# (3) Merging filtered data
rmarkdown::render('03_run_filter_object.Rmd', 
                  clean = TRUE,
                  output_dir = file.path(Filter_object_dir),
                  output_file = paste('Report-', 'Filter-object', '-', Sys.Date(), sep = ''),
                  output_format = 'all',
                  params = list(assay = yaml$assay_signac_qc,
                                min.cutoff_value = yaml$min.cutoff_value_upstream,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                use_condition_split = yaml$use_condition_split_filter_object,
                                print_pdf = yaml$print_pdf_filter_object,
                                grouping = yaml$grouping,
                                use_scDblFinder_filtering = yaml$use_scDblFinder_filtering_filter_object,
                                filename_filter_object = yaml$filename_filter_object_value,
                                root_dir = yaml$root_dir,
                                metadata_dir = yaml$metadata_dir,
                                metadata_file = yaml$metadata_file,
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
                                
################################################################################################################
# (4) Final QC summary report
rmarkdown::render('04_run_summary_report.Rmd', 
                  clean = TRUE,
                  output_dir = file.path(Final_summary_report_dir),
                  output_file = paste('Report-', 'Final-summary', '-', Sys.Date(), sep = ''),
                  output_format = 'all',
                  params = list(use_scDblFinder_filtering = yaml$use_scDblFinder_filtering_summary_report,
                    cellranger_parameters = yaml$cellranger_parameters,
                    filename_summary_report = yaml$filename_summary_report_value,
                    filename_filter_object = yaml$filename_filter_object_value,
                    root_dir = yaml$root_dir,
                    metadata_dir = yaml$metadata_dir,
                    metadata_file = yaml$metadata_file,
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

################################################################################################################   

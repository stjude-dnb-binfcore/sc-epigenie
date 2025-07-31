#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
#.libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)})

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
analysis_dir <- file.path(root_dir, "analyses", "cluster-cell-calling") 
report_dir <- file.path(analysis_dir, "plots") 
data_dir_module <- yaml$data_dir_module_name

################################################################################################################
# Run Rmd scripts 
################################################################################################################
future_globals_value = 214748364800 #200 * 1024^3; # 150 * 1024^3; other options: 1000 * 1024^2 = 1048576000; 8000 * 1024^2 =8388608000
resolution = yaml$resolution_find_markers
################################################################################################################

rmarkdown::render('03-find-peaks.Rmd', clean = TRUE,
                  output_dir = file.path(report_dir),
                  output_file = c(paste('Report-', 'find-peaks', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                assay = yaml$assay_clustering_module,
                                resolution_list = yaml$resolution_list_find_markers, 
                                n_value = yaml$n_value_find_markers,
                                genome_name = yaml$genome_name_cellranger,
                                ident.1_value = yaml$ident.1_value,
                                condition = yaml$condition_value1,
                                min.cutoff_value = yaml$min.cutoff_value_upstream,
                                OrgDb_value = yaml$OrgDb_value,
                                top_n_value_peaks = yaml$top_n_value_peaks,
                                root_dir = yaml$root_dir,
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


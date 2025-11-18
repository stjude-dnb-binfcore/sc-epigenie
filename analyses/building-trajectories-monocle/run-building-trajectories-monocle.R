#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
#.libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
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
# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
assay <- yaml$assay_annotation_module
analysis_dir <- file.path(root_dir, "analyses", "building-trajectories-monocle") 

module_plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)}

################################################################################################################
rmarkdown::render('01-building-trajectories-monocle.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'building-trajectories-monocle', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(assay = yaml$assay_annotation_module,
                                resolution_list = yaml$resolution_list_find_markers, 
                                
                                cell_type_name = yaml$cell_type_name_module,
                                lineage_number_value = yaml$lineage_number_value_module,
                                lineage_value = yaml$lineage_value_module,
                                lineage1_value = yaml$lineage1_value_module,
                                lineage1_value_vector = yaml$lineage1_value_vector_module,
                                lineage2_value = yaml$lineage2_value_module,
                                lineage2_value_vector = yaml$lineage2_value_vector_module,
                                
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

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
analysis_dir <- file.path(root_dir, "analyses", "integration-with-scrna-seq-data") 

module_plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)}

################################################################################################################
future_globals_value = 214748364800 # 200*1024^3; other options: 1000 * 1024^2 = 1048576000; 8000 * 1024^2 =8388608000
input_data = yaml$module_with_input_data
input_data_folder= yaml$input_data_folder_name

################################
# Set data_dir
# Caution! Sometimes this file will be located in the `cluster-cell-calling` module
# BUT if we had to remove contamination, then it will be located in the `cell-contamination-removal-analysis` module
data_dir_annotation_module <- file.path(root_dir, "analyses", input_data, "results", input_data_folder)
input_data_file <- file.path(data_dir_annotation_module, glue::glue("seurat_obj_gene_activity_matrix.rds"))

################################################################################################################
# Integrating with scRNA-seq
rmarkdown::render('01-integration-with-scrna-seq-data.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'integration-with-scrna-seq-data', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(reduction_value = yaml$reduction_value_annotation_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                data_file = input_data_file,
                                assay = yaml$assay_annotation_module,
                                nfeatures_value = yaml$nfeatures_value_module,
                                ct_palette_file = yaml$ct_palette_file_value,
                                group_by_variable = yaml$group_by_variable_module,

                                reference_dir = yaml$reference_dir_annotation_module,
                                reference_file_name = yaml$reference_file_name_annotation_module,
                                celltype_reference = yaml$celltype_reference_module,

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
future_globals_value = 214748364800 #200 * 1024^3; # 150 * 1024^3; other options: 1000 * 1024^2 = 1048576000; 8000 * 1024^2 =8388608000
resolution = yaml$resolution_find_markers
################################################################################################################

rmarkdown::render('02-peak-calling.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'peak-calling', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(assay = yaml$assay_clustering_module,
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


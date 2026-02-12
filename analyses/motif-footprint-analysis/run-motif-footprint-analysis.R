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
  library(BSgenome.Mmusculus.UCSC.mm39)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(BSgenome.Hsapiens.UCSC.hg38)
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
analysis_dir <- file.path(root_dir, "analyses", "motif-footprint-analysis") 

module_plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)}


################################################################################################################
# BSgenome package: 
# (1) For human: BSgenome.Hsapiens.UCSC.hg38 
# (2) For mouse: BSgenome.Mmusculus.UCSC.mm10 or BSgenome.Mmusculus.UCSC.mm39
genome = BSgenome.Mmusculus.UCSC.mm39 

# JASPAR database version to use for motif analysis. 
# It dependens on genome reference and versioning used in the project. 
# Please check the JASPAR database for more details: https://jaspar.genereg.net/downloads/.
jaspar_library_version = "2024" # Options: "2018", "2020", "2022" or "2024". 

#################################################################################
rmarkdown::render('01-motif-analysis.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'motif-analysis', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(data_dir_motif = yaml$data_dir_motif_module,
                                data_file_motif = yaml$data_file_motif_module,
                                assay = yaml$assay_annotation_module,
                                resolution_list = yaml$resolution_list_find_markers, 
                                min.cutoff_value = yaml$min.cutoff_value_motif_module,
                                max.cutoff_value = yaml$max.cutoff_value_motif_module,
                                top_n_value_motifs = yaml$top_n_value_motifs_module,
                                top_da_peak_names_motifs = yaml$top_da_peak_names_motifs_module,
                                cell_type_name = yaml$cell_type_name_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,

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


################################################################################################################
# JASPAR database doesn’t differentiate by genome build; it uses species identifiers.
# (1) For human: 9606 
# (2) For mouse (mm39, GRCm39, mm19, GRCm19): 10090
species = 10090

#################################################################################
rmarkdown::render('02-motif-footprinting-analysis.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'motif-footprinting-analysis', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(assay = yaml$assay_annotation_module,
                                cell_type_name = yaml$cell_type_name_module,
                                top_n_value_footprinting = yaml$top_n_value_footprinting_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                motifs_to_fp_module_list = yaml$motifs_to_fp_module_list_value,

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

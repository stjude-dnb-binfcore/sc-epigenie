#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(reticulate)
  library(gridExtra)
  library(glue)
  library(R.utils)
  library(shiny)
  library(shinyhelper)
  library(DT)
  library(magrittr)
  library(ggdendro)
  library(ShinyCell)
  library(RColorBrewer)
  library(Signac)
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
# Parameters
root_dir <- yaml$root_dir
PROJECT_NAME <- yaml$PROJECT_NAME
PI_NAME <- yaml$PI_NAME
assay <- yaml$assay_signac_qc
annotations_dir <- yaml$annotations_dir_rshiny_app
annotations_filename <- yaml$annotations_filename_rshiny_app
reduction_value <- yaml$reduction_value_annotation_module


# Set up directories and paths to root_dir and analysis_dir
analysis_dir <- file.path(root_dir, "analyses") 
module_dir <- file.path(analysis_dir, "rshiny-app") 

annotation_results_dir <- file.path(analysis_dir, "integration-with-scrna-seq-data", "results") 
annotations_all_results_dir <- file.path(annotation_results_dir, annotations_dir) 

# Input files
annotations_all_file <- file.path(annotations_all_results_dir, annotations_filename)

# Create results_dir
results_dir <- file.path(module_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}


################################################################################
### library(ShinyCell) ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
################################################################################
# It was downloaded from here: https://github.com/SGDDNB/ShinyCell/blob/master/R
# We modified the function to account for any type of assay
source(paste0(module_dir, "/util/makeShinyFiles_assay.R"))

################################################################################################################
### Generate R shiny app ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################################################################################################
cat("Beginning to process results from", "annotations_file", "\n")
seu1 <- readRDS(annotations_all_file)
DefaultAssay(seu1) <- assay

# Check if the specified dimensionality reduction exists in the Seurat object
if (reduction_value %in% names(seu1@reductions)) {
  
  # Retrieve the key prefix for the specified reduction (e.g., "PC_", "UMAP_")
  reduc_key <- seu1[[reduction_value]]@key
  message("Setting default reduction to ", reduction_value)
  
  # Construct names for dimension 1 and 2 (e.g., "PC_1", "PC_2")
  dim_1 <- paste0(reduc_key, "1")
  dim_2 <- paste0(reduc_key, "2")
  
  # Identify all metadata columns associated with the specified reduction
  cols_to_remove <- grep(reduc_key, colnames(seu1@meta.data), value = TRUE)
  
  # Remove the identified columns from the metadata
  seu1@meta.data <- seu1@meta.data[, !colnames(seu1@meta.data) %in% cols_to_remove]
  
} else {
  # If the reduction is not present, notify that the default will be set automatically
  message(reduction_value, " not found in Seurat object, default reduction will be determined by the function")
}

scConf1 <- createConfig(seu1, 
                        meta.to.include = NA, # Include all metadata (or specify if you want a subset)
                        maxLevels = 150)     # Use the number of unique levels

makeShinyFiles_assay(seu1, 
                     scConf1, 
                     shiny.prefix = "sc1", 
                     default.dimred = c(dim_1, dim_2),
                     shiny.dir = paste(results_dir, "shinyApp", sep = "/"))

cat("Make R shiny app for all files", "\n")
makeShinyCodesMulti(
  shiny.title = PROJECT_NAME, 
  shiny.footnotes = PI_NAME,
  shiny.prefix = c("sc1"),
  shiny.headers = c("Annotations"),
  shiny.dir = paste(results_dir, "shinyApp", sep = "/")) 
################################################################################################################   

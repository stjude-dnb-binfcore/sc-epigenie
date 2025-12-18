# Pipeline for Pipeline for Performing DNA sequence motif and motif footprinting analysis for sc-ATAC-Seq Analysis in 10X Genomics data

## Usage

`run-motif-footprint-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `run-motif-footprint-analysis.R`: User needs to specify the following:
   
   - `genome` to be used for the BSgenome package based on human or mouse experiment and genome build.
   - `species` to be used for the JASPAR database based on human or mouse experiment


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-motif-footprint-analysis.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains scripts tasked to performing:

- DNA sequence motif analysis in Signac and
- footprinting any motif that have positional information for.

For more information, see [Motif analysis with Signac](https://stuartlab.org/signac/articles/motif_vignette) and [Motif footprinting](https://stuartlab.org/signac/articles/footprint) .


## Folder structure 

The structure of this folder is as follows:

```
├── 01-motif-analysis.Rmd
├── 02-motif-footprinting-analysis.Rmd
├── lsf_script.txt
├── plots
|   ├── 01_motif_analysis
|   ├── 02_footprinting_analysis
|   ├── Report-motif-analysis-<Sys.Date()>.html
|   ├── Report-motif-analysis-<Sys.Date()>.pdf
|   ├── Report-motif-footprinting-analysis-<Sys.Date()>.html
|   └── Report-motif-footprinting-analysis-<Sys.Date()>.pdf
├── README.md
├── results
|   ├── 01_motif_analysis
|   └── 02_footprinting_analysis
├── run-motif-footprint-analysis.R
├── run-motif-footprint-analysis.sh
└── util
|   ├── function-add-footprinting.R
|___└── function-add-motif.R
```


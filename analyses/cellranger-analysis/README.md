# Pipeline for running and summarizing Cell Ranger count for single or multiple libraries for sc-/sn-ATAC-Seq Analysis in 10X Genomics data

## Usage

`submit-multiple-jobs.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

The `submit-multiple-jobs.sh` script is designed to run the following two steps: 
   - Step 1: To run the `j1.sh` script to align single or multiple libraries in parallel, i.e., `run-cellranger-analysis`. 
   - Step 2: To run `j2.sh` to summarize alignment results, i.e., `summarize-cellranger-analysis`. The latter script will be on hold and executed once all libraries are aligned and `j1.sh` is complete. This is been taken care of by `waiter.sh` script.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `../../project_parameters.Config.yaml`: define the `metadata_dir`, `genome_reference_path`, `cellranger_parameters`, and `genome_name_cellranger`. Please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues) to request the path to the reference genome of preference. Our team at the Bioinformatics core at DNB maintains the following genome references: `GRCh38`, `mm10`, and `GRCh38_mm10`, for human, mouse, and dual index genomes, respectively. Otherwise, specify the path to the reference genome of your preference.

- `waiter.sh`: Here, the user will need to replace `DST` with the Sample ID used for the samples of the project. This requires that all of the Sample IDs are the same, e.g., DST800, DST802, DST811 and so on.


### Run module on HPC

To run all of the scripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bsub < submit-multiple-jobs.sh
```

Please note that this will run the analysis module outside of the container while submitting lsf job on HPC. This is currently the only option of running the `cellranger-analysis` module. By default, we are using `cellranger-atac/2.1.0` as available on St Jude HPC.


## Folder content

This folder contains scripts tasked to run and summarize Cell Ranger count for single or multiple libraries for sc-/sn-ATAC-Seq Analysis in 10X Genomics data across the project. For more information and updates, please see [Cell Ranger support page](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).

This module uses CellRanger-atac v2.1.0 for the alignment.


## Folder structure 

The structure of this folder is as follows:

```
├── j1.sh
├── j2.sh
├── README.md
├── results
|   ├── 01_logs
|   ├── 02_cellranger_count
|   |   └── DefaultParameters
|   └── 03_cellranger_count_summary
|       └── DefaultParameters
├── submit-multiple-jobs.sh
└── waiter.sh
```

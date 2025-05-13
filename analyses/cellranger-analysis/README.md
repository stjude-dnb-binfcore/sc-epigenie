# Pipeline for running and summarizing Cell Ranger count for single or multiple libraries for single-cell ATAC (scATAC-Seq) data analysis in 10X Genomics data

## Usage

`submit-multiple-jobs.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

The `submit-multiple-jobs.sh` script is designed to run the following two steps: 
   - Step 1: To run the `j1.sh` script to align single or multiple libraries in parallel, i.e., `run-cellranger-analysis`. 
   - Step 2: To run `j2.sh` to summarize alignment results, i.e., `summarize-cellranger-analysis`. The latter script will be on hold and executed once all libraries are aligned and `j1.sh` is complete. This is been taken care of by `waiter.sh` script.


Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `../../project_parameters.Config.yaml`: define the `metadata_dir`, `genome_reference_path`, `cellranger_parameters`, and `genome_name_cellranger`. For a list of genome references maintained and supported by the Bioinformatics Core at DNB, please review the [Genome References wiki page](https://github.com/stjude-dnb-binfcore/sc-atac-seq/wiki/2.-Genome-References). Please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-atac-seq/issues) to request the path to the reference genome of preference. Otherwise, specify the path to the reference genome of your preference. 

User also need to define `sample_prefix` with the Sample ID used for the samples of the project. Sample IDs should follow a format like: PREFIX001 (e.g., DYE001, ABC-002, XYZ_003). You can specify multiple prefixes if your project uses more than one.

- `j1.sh`: 
  - `--force-cells`: User can add flags as necessary for their analyses and compare alignments with CellRanger, e.g., by using `--force-cells=8000` to constrain the number of cells to be used for alignment. We recommend running by default and, after careful assessment, editing parameters. We have found that the default parameters set up here work well for most cases.


### Handling Top-Ups and Technical Replicates in Cell Ranger Count

If a sample (`ID` column) has multiple technical replicates (i.e., multiple sequencing runs of the same library), list all corresponding FASTQ file paths in a single row of the metadata file (`FASTQ` column), separated by commas with no spaces. Use the same format for samples with multiple library names (`SAMPLE` column). Make sure that sample names and their FASTQ paths are listed in matching order (`SAMPLE` and `FASTQ` column, respectively). 


Example Metadata Format 1:

| ID | SAMPLE | FASTQ | 
:----------|:----------|:----------|
| DYE001 | seq_submission_code1_sample1 | /absolute_path/seq_submission_code1/replicate1,/absolute_path/seq_submission_code1/replicate2 | 

Example Metadata Format 2:

| ID | SAMPLE | FASTQ | 
:----------|:----------|:----------|
| DYE001 | seq_submission_code1_sample1,seq_submission_code2_sample2 | /absolute_path/seq_submission_code1,/absolute_path/seq_submission_code2 | 


Cell Ranger will automatically recognize these as multiple libraries for the same sample (`ID` column) and merge them into a single output for the sample during processing. There is no need to manually combine or rename the files—simply format them correctly in the metadata file, and the pipeline will handle the rest.


### Run module on HPC

To run all of the scripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bsub < submit-multiple-jobs.sh
```

Please note that this will run the analysis module outside of the container while submitting lsf job on HPC. This is currently the only option of running the `cellranger-analysis` module. By default, we are using `python/3.9.9` and `cellranger-atac/2.1.0` as available on St Jude HPC.


## Folder content

This folder contains scripts tasked to run and summarize Cell Ranger count for single or multiple libraries for single-cell ATAC (scATAC-Seq) data analysis in 10X Genomics data across the project. For more information and updates, please see [Cell Ranger support page](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).

This module uses CellRanger-atac v2.1.0 for the alignment.

For more information, on how to review  the web summary file in the output folder of the Cell Ranger ATAC analysis software, please see the [Interpreting Cell Ranger ATAC Web Summary Files for Single Cell ATAC Assay](https://www.10xgenomics.com/support/epi-atac/documentation/steps/sequencing/interpreting-cell-ranger-atac-web-summary-files-for-single-cell-atac-assay) file.


## Folder structure 

The structure of this folder is as follows:

```
├── j1.sh
├── j2.sh
├── README.md
├── results
|   ├── 01_logs
|   ├── 02_cellranger_count
|   |   └── ${cellranger_parameters}
|   └── 03_cellranger_count_summary
|       └── ${cellranger_parameters}
├── submit-multiple-jobs.sh
├── util
|   └── summarize_cellranger_results.py
└── waiter.sh
```

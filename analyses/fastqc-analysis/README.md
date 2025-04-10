# Pipeline for FastQC quality control tool for high throughput sequence data analysis

## Usage

`run-fastqc-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml`: define `metadata_dir`. FASTQ paths to the fastqc files with format: `path1/*R1-R3*.fastq.gz` are extracted from the `metadata_dir`.

If the module needs to be run more than one time, user will need to remove the `02-multiqc-reports` folder before rerunning the module or the code will give an error at that step. Files and folder related to the MultiQC step will be generated every time a new run is performed. Folder can be deleted manually or from the node as:

```
rm -r 02-multiqc-reports
```


### Run module on an interactive session on HPC within the container

To run the script on an interactive session on HPC, please run the following command from an interactive compute node (while within the container):

```
bash run-fastqc-analysis.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content

This folder contains a script tasked to run FastQC quality control tool for all libraries across the project. Each libary directory contains the following files:
- R1 (50bp forward read): This file contains the actual sequence of interest from the experiment. For ATAC-seq or similar applications, R1 typically corresponds to the sequencing of the fragment after the transposase action (for ATAC-seq, this is where you are looking for the inserted fragments). So R1 is crucial for the actual biological signal.

- R2 (16bp barcode sequence): This file contains the barcode sequences, which are used to identify and demultiplex samples. For 10x sequencing or similar platforms, barcodes are essential for assigning reads to the correct cell or sample. This sequence does not provide biological sequence information but is crucial for downstream analysis, such as mapping each read to its corresponding sample or cell.

- R3 (49bp reverse read): This file contains the reverse mate-pair read in a paired-end sequencing setup. This will provide complementary information to R1, allowing for a more complete understanding of the fragment and the sequencing of the other end of the DNA fragment. This file is important for creating a more complete picture of the sequencing fragment in paired-end setups.

We will run FastQC on:
- R1 (50bp forward read) to assess the quality of the actual sequencing data for your experiment;
- R3 (49bp reverse read) to check the reverse reads, especially if you're concerned about paired-end read quality.


For more information, please:
- Type from the command line: fastqc --help, or
- See [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Folder structure 

The structure of this folder is as follows:

```
├── lsf-script.txt
├── README.md
├── results
|   ├── 01-fastqc-reports
|   ├── 02-multiqc-reports
|   └──multiqc_report.html
└── run-fastqc-analysis.sh
```


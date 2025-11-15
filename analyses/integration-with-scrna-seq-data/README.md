# Pipeline for Peak calling per cell type following integration with scRNA-seq data for sc-ATAC-Seq Analysis in 10X Genomics data


## Usage

`run-integration-with-scrna-seq-data.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `future_globals_value` is hard-coded in each `run-integration-with-scrna-seq-data.R` script for each method. If necessary, user can increase/decrease resources.

We provide a fixed color palette for annotating cell types. Please ensure that the cell type names in your customized cell type marker reference list and/or reference dataset  (if you use any) exactly match those in the `.figures/palettes/cell_types_palette.tsv`. If any cell types are missing, please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues) with the list of missing types, and we will add them to the palette. Alternatively, users are welcome to use their own `cell_types_palette.tsv`, as long as the column names (`cell_type_names` and `hex_codes`), file name, and format match those in `./figures/palettes/cell_types_palette.tsv`. We reccomend that users use the same color palette files across pipelines and projects, i.e., same files as for the `sc-rna-seq-snap` pipeline.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-integration-with-scrna-seq-data.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains a script tasked to perform a data transfer method in the context of scATAC-seq to:

- Classify cells measured with scATAC-seq based on clustering results from scRNA-seq
- Co-embed scATAC-seq and scRNA-seq data

There are also scripts:
   (1) to perform peak calling and gene ontology analysis per cell type and 
   (2) to predict regions of the genome that are more likely to be in physical proximity in the nucleus, i.e., [Cicero’s co-accessibility analysis](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/).


## Folder structure 

The structure of this folder is as follows:

```
├── 01-integration-with-scrna-seq-data.Rmd
├── 02-peak-calling.Rmd
├── 03-finding-co-accessible-networks-cicero.Rmd
├── lsf_script.txt
├── plots
├── README.md
├── results
├── run-integration-with-scrna-seq-data.R
├── run-integration-with-scrna-seq-data.sh
└── util
|   ├── co-embedding-cells.R
|   ├── function_process_fragments.R
|___└── function-cell-type-fractions.R
```


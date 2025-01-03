# Pipeline for Data exploratory analysis


## Usage

`run-data-exploratory-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

User will need to install the R packages as specified in the Rmd script and for the R version specified in the `run-data-exploratory-analysis.sh` script.


### Run module on an interactive session on HPC within the container

To run this module on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-data-exploratory-analysis.sh
```


## Folder content
This folder contains a script tasked to analyze the distribution of samples across the project.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-run-data-exploratory-analysis.Rmd
├── input
├── plots
|   ├── Report-data-exploratory-analysis-<Sys.Date()>.html
|   └── Report-data-exploratory-analysis-<Sys.Date()>.pdf
├── README.md
├── run-data-exploratory-analysis.R
└── run-data-exploratory-analysis.sh
```
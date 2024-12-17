# Pipeline for Sample data analysis

## Usage

To run the script in this module from the command line sequentially, use:

```
bash run-sample-distribution-analysis.sh
```

`run-sample-distribution-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `run-sample-distribution-analysis.R`: define `root_dir`, and `analysis_dir`


## Folder content
This folder contains a script tasked to analyze the distribution of samples across the project.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-sample-distribution-analysis.Rmd
├── plots
    └── Report_<Sys.Date()>.html
├── README.md
└── run-sample-distribution-analysis.sh
```
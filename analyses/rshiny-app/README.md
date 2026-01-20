# Pipeline for generating an R shiny app for the project

## Usage

`run-rshiny-app.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.


### Run module on an interactive session on HPC within the container

To run the module on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-rshiny-app.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Launch the R Shiny app

To launch the R Shiny app, navigate to the script located at `./results/shinyApp/server.R`. Open it in RStudio on your local machine, and click the `Run App` button in the top-right corner to start the application.


## Folder content

This folder contains scripts designed to run and generate an R shiny app for the project.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-generate-rshiny-app.R
├── lsf-script.txt
├── README.md
├── run-rshiny-app.sh
└── util
|___└── makeShinyFiles_assay.R
```

# Pipeline for Building trajectories with Monocle 3 for sc-ATAC-Seq Analysis in 10X Genomics data

## Usage

`run-building-trajectories-monocle.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-building-trajectories-monocle.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains a script tasked to infer trajectories and pseudotime on single-cell ATAC-seq (scATAC-seq) data.

For more information, see [Building trajectories with Monocle 3](https://stuartlab.org/signac/articles/monocle) and [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/).


## Folder structure 

The structure of this folder is as follows:

```
├── 01-building-trajectories-monocle.Rmd
├── lsf_script.txt
├── plots
|   ├── 01_building_trajectories_monocle
|   ├── Report-building-trajectories-monocle-<Sys.Date()>.html
|   └── Report-building-trajectories-monocle-<Sys.Date()>.pdf
├── README.md
├── run-building-trajectories-monocle.R
├── run-building-trajectories-monocle.sh
└── util
|___└── monocle3_updated.R
```


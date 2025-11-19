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

We build trajectories with Monocle 3. By default, our pipeline roots cells programmatically—automatically selecting starting cells based on a specified cluster or cell type (for example, using the `cell_type_name` column and the `lineage_value`  parameter). This approach is useful for exploratory analyses or when the user does not wish to predefine lineages. For more details on this topic, see [Monocle3 issue](https://github.com/cole-trapnell-lab/monocle3/issues/328).

If the user prefers, they can specify one or two lineages in the parameters. The pipeline will then subset the data and build a separate trajectory for each specified lineage (up to two). 

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
|   ├── function-predict-trajectories.R
|___└── monocle3_updated.R
```


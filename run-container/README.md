# How to run the container for the Single cell/Single-nucleus ATAC Seq workflow (ScATACSeq)

We have generated Dockerfile and Definition file that contain all tools, packages, and dependencies necessary to run the code and analyses modules in the [sc-rna-seq-snap repository](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap). These are customized for `Rstudio/R v4.4.0` and `Seurat v4.4.0`. We are using the same Docker image for the `sc-atac-seq` repo. 

For more information, see [./sc-rna-seq-snap repository/run-container/README.md file](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/blob/main/run-container/README.md).


## To use the container in an R interactive session on HPC:

From an interactive node on HPC, the user can open their session. Please modify memory and resources as needed for the analysis module to run.
```
bsub -P hpcf_interactive -J hpcf_interactive -n 2 -q standard -R "rusage[mem=16G]" -Is "bash"
```

## Load specific version of Singularity

Please note that a version of Singularity is installed by default on all the cluster nodes at St Jude HPC. 
Otherwise the user needs to ensure and load Singularity module by running the following on HPC:
```
module load singularity/4.1.1
```


1. Pull the singularity container from the `sc-rna-seq-snap` root_dir
```
singularity pull docker://achronistjude/rstudio_4.4.0_seurat_4.4.0:latest
```


2. Start the singularity container

a. To run analysis modules via lsf

All analysis modules (except for `.analyses/cellranger-analysis`) are designed to be run while executing the container. User only needs to run the lsf script as described in the `README.md` files in each analysis module.


b. To run from the terminal

User can run analysis module while on interactive node after executing the container:

```
bash run-terminal.sh
```

Then user may navigate to their module of interest, `./sc-rna-seq-snap/analyses/<module_of_interest>`. For example:
```
cd ./sc-rna-seq-snap/analyses/upstream-analysis
bash run-upstream-analysis.sh
```


c. To run from Rstudio

User can also run analyses via Rstudio OnDemand after executing the container:

```
bash run-rstudio.sh
```

The `run-rstudio.sh` is running at `IP_ADDR:PORT`. When RStudio launches, please click "Session" -> "Restart R" (at the RStudio web session). For St Jude users, we advice to disconnect from CloudFlare WARP as this might lead to unstable behavior while on VPN.

Again, the user can navigate to their module of interest and explore/run their analyses.


If you encounter issues during this step related to RStudio Server and specifically to an invalid secure cookie error. This might be an issue with how the secure cookie is being handled during an HTTP request. In this case, please check if the following directories have been generated and if so, remove them:
```
rm -r .cache/
```
```
rm -r .config/
```
```
rm -r .local/
```
```
rm -r rstudio-container-tmp/
```

These folders cache history and user info. Then, kill the interactive session, start a new one, and hopefully, it works! 🎉



## Authors

Antonia Chroni, PhD ([@AntoniaChroni](https://github.com/AntoniaChroni)) and 
Walid Abu Al-Afia ([@walidabualafia](https://github.com/walidabualafia))


## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

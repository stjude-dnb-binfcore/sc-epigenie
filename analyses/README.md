# Using Analysis Modules in Single-cell ATAC-seq workflow (scEpiGenie)

This repository contains a collection of analysis modules designed to process and analyze single-cell ATAC (scATAC-Seq) data from 10X sequencing technology. 

**Analysis modules**

Each module is self-contained and can be executed independently or as part of a larger analysis pipeline. Below is a summary of each analysis module, including whether they are required or optional. Furthermore, the analysis modules should be run in the following recommended order:

1. `fastqc-analysis` module (description="Pipeline for FastQC quality control tool for high throughput sequence data analysis.", required=True)
2. `cellranger-analysis` module (description="Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.", required=True)
3. `upstream-analysis` module (description="Pipeline for estimating QC metrics and filtering low quality cells.", required=True)
4. `integrative-analysis` module (description="Pipeline for Integrative analysis.", required=False)
5. `cluster-cell-calling` module (description="Pipeline for cluster cell calling and differentially accessible peaks analysis.", required=True)
6. `integration-with-scrna-seq-data` module (description="Pipeline for Peak calling per cell type following integration with scRNA-seq data for sc-ATAC-Seq Analysis.", required=True)
7. `building-trajectories-monocle` module (description="Pipeline for Building trajectories with Monocle 3 for sc-ATAC-Seq Analysis.", required=False)
8. `rshiny-app` module (description="Pipeline for generating an R shiny app for the project.", required=False)
9. `project-updates` module (description="Pipeline for summarizing results from all modules and generating project reports.", required=False)


## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-epigenie/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

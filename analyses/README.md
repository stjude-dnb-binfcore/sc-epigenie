# Using Analysis Modules in scEpiGenie: A Workflow for Single-cell ATAC-seq (scATAC-seq)

This repository contains a collection of analysis modules designed to process and analyze single-cell ATAC (scATAC-Seq) data from 10X sequencing technology. 

**Analysis modules**

Each module is self-contained and can be executed independently or as part of a larger analysis pipeline. Below is a summary of each analysis module, including whether they are required or optional. Furthermore, the analysis modules should be run in the following recommended order:

1. `fastqc-analysis` module (description="Pipeline for FastQC quality control tool for high throughput sequence data analysis.", required=True)
2. `cellranger-analysis` module (description="Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.", required=True)
3. `upstream-analysis` module (description="Pipeline for estimating QC metrics and filtering low quality cells.", required=True)
4. `integrative-analysis` module (description="Pipeline for Integrative analysis.", required=False)
5. `cluster-cell-calling` module (description="Pipeline for cluster cell calling and gene marker analysis.", required=True)
________________________________________________________________________________________ 
🚧🚧🚧 Analysis Modules Under Development 🚧🚧🚧

The following analysis modules are currently under development. Please check back soon for updates and stay tuned for what's coming next!

6. `project-updates` module (description="Pipeline for summarizing results from all modules and generating project reports.", required=False)
7. `integration-with-scrna-seq-data` module (description="Pipeline for integrating scATAC-seq with scRNA-seq data.", required=True)

________________________________________________________________________________________


## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-epigenie/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

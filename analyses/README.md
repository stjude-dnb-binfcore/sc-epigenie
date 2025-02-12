# How to use analysis modules in the Single cell ATAC Seq workflow (ScATACSeq)

This repository contains a collection of analysis modules designed to process and analyze single cell ATAC (scATAC) data from 10X sequencing technology. 

**Analysis modules**

Each module is self-contained and can be executed independently or as part of a larger analysis pipeline. Below is a summary of each analysis module, including whether they are required or optional. Furthermore, the analysis modules should be run in the following recommended order:


🚧 ⚠️🚧 ⚠️🚧 ⚠️🚧 ⚠️🚧⚠️ **Currently under construction** 🚧 ⚠️ 🚧 ⚠️🚧 ⚠️🚧 ⚠️


1. `cellranger-analysis` module (description="Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.", required=True)

2. `upstream-analysis` module (description="Pipeline for estimating QC metrics and filtering low quality cells.", required=True)




## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-atac-seq/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

# Pipeline for estimating QC metrics for sc-ATAC-Seq Analysis in 10X Genomics data

## Usage

`run-upstream-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.
- `future.globals.maxSize` is hardwired coded in the `run-upstream-analysis.R`. If necessary, user can increase/decrease resources.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node (while within the container):

```
bash run-upstream-analysis.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content

This folder contains scripts tasked to:
(1) Infer QC metrics and associated plots to visually explore the quality of each library of the project.
(2) Evaluate QC metrics and set filters to remove low quality cells in 10x single-cell-ATAC-sequencing libraries.

## QC Steps and methods

The pipeline offers flexibility for users to include or exclude methods and adjust the workflow during QC based on factors such as sequence type, expected cell number, experiment type, and genome reference.

By default, the pipeline runs all methods from steps (1-2). Step (1) is mandatory for basic QC filtering, while integration of step (2) is optional. This can be configured in the `project_parameters.Config.yaml` file. However, we recommend reviewing all results, as this can provide valuable insights into the overall quality of each library.


### (1) Signac QC metrics

[Signac](https://stuartlab.org/signac/articles/pbmc_vignette) workflow is implemented to pre-process, filter and plot the ATAC-sequencing data. For more tutorials, see [Introduction to single cell ATAC data analysis in R](https://www.youtube.com/watch?v=e2396GKFMRY&ab_channel=Sanbomics) and [How to analyze single-cell ATAC-Seq data in R | Detailed Signac Workflow Tutorial](https://www.youtube.com/watch?v=yEKZJVjc5DY&ab_channel=Bioinformagician).

The CellRanger output from the `cellranger-analysis` module will be used for this step. User will have to define `params` as needed for their experiment. 
  - Calculate QC metrics: Nucleosome banding pattern, Transcriptional start site (TSS) enrichment score, Total number of fragments in peaks, Fraction of fragments in peaks, Ratio reads in genomic blacklist regions, etc.
  - Before and after filter: Plot "pct_reads_in_peaks", "peak_region_fragments", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal", "nCount_peaks", "nFeature_peaks", "pct_reads_in_peaks_promoters", "pct_reads_in_peaks_enhancers", "duplicate", "mitochondrial".
  - Data are normalized by using term frequency-inverse document frequency (TF-IDF) normalization. Then, we select the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells with the [FindTopFeatures()](https://stuartlab.org/signac/reference/findtopfeatures) function. We next run singular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above. This returns a reduced dimension representation of the object (for users who are more familiar with scRNA-seq, you can think of this as analogous to the output of PCA). The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI). After that the cells are embedded in a low-dimensional space we can run UMAP from the Seurat package.

Here, the user can select to implement the following strategies to remove low quality cells via the `run_QC` function:
- `approach 1` Filter cells based on various variables as defined in the `params`.
- `approach 2` Use percentile filtering. This will be defined in the `params` by setting `use_threshold_filtering_upstream: "NO"`.

Moreover, only libraries with more than 500 cells will be kept for merging and integration purposes. This value is a commonly used threshold for many single-cell ATAC-seq studies as a minimum for obtaining reliable and meaningful analysis. 
- Statistical Power: At least 500 cells are typically needed to ensure the analysis has enough statistical power to detect meaningful biological signals.
- Cell Diversity: With fewer cells, you may not capture sufficient cellular diversity, leading to incomplete or biased results.
- Clustering: Some clustering algorithms in single-cell RNA-seq require a minimum number of cells to create robust, meaningful clusters.

### Genome-to-Ensembl Mapping Options

We need tp assign gene annotations from EnsDb based on the genome reference used for the cohort. For more information on the various Ensembl releases, you can refer to the [Table of Assemblies](https://useast.ensembl.org/info/website/archives/assembly.html) and [Build Notes for Reference Packages](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps#ref-2020-a).

| Species | Genome Build | Ensembl Version(s)        | Annotation Package                              |
|---------|--------------|---------------------------|--------------------------------------------------|
| Human   | hg19 (GRCh37) | v75                       | `EnsDb.Hsapiens.v75`                             |
| Human   | hg38 (GRCh38) | v86 (older), v98 (newer)  | `EnsDb.Hsapiens.v86`, `EnsDb.Hsapiens.v98`       |
| Mouse   | mm9 (NCBI37)  | —                         | *No official* `EnsDb`; use `TxDb.Mmusculus.UCSC.mm9.knownGene` |
| Mouse   | mm10 (GRCm38) | v79                       | `EnsDb.Mmusculus.v79`                            |
| Mouse   | mm39 (GRCm39) | v104, v105                | `EnsDb.Mmusculus.v104`, `EnsDb.Mmusculus.v105`   |


### 🧪 Key QC Metrics and Suggested Thresholds (Per Cell/Nucleus)

| **Metric**                  | **Interpretation**                                                 | **Low-Quality Threshold**  | **High-Quality Range**    |
|-----------------------------|--------------------------------------------------------------------|----------------------------|---------------------------|
| **Fragments in peaks**      | Proportion of reads in peaks — reflects signal-to-noise                                | < 20–30%                   | > 30–50%                  |
| **TSS enrichment score**    | Enrichment of reads around transcription start sites                                   | < 4–6                      | > 6–10                    |
| **Nucleosome signal**       | Ratio of mono- to di-/tri-nucleosome fragments (chromatin state)                       | > 4                        | < 2–3                     |
| **Total fragments**         | Total number of fragments per cell. Avoid very low or very high.                       | < 1k or > 50k              | 3k–30k                    |
| **Blacklist ratio**         | Fraction of reads in ENCODE blacklist regions (noise -High = poor quality)             | > 0.01                     | < 0.01                    |
| **Mono-/multi-nucleosome**  | Fragment length distribution should show clear banding            | Poor/no banding    | Visible nucleosome bands   |
| **Doublet score**           | Estimate of multiplets (2+ cells) per barcode (Use to exclude multiplets)              | High                       | Low (filtered out)        |


### Quality control

[Baek and Lee, 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303019#s0015) emphasized the importance of filtering out barcodes corresponding to low-quality cells or doublets after processing sequencing read data. In general, single-cell sequencing QC relies on metrics like read counts (count depth) and feature counts per barcode ([Luecken and Theis, 2019](https://www.scopus.com/record/display.uri?eid=2-s2.0-85067863532&origin=inward&txGid=44862f35717b260b21d411030ddd82fe)). Barcodes with either abnormally low or high read/feature counts are typically flagged as low-quality cells or multiplets, respectively.

However, scATAC-seq data offers additional QC metrics that better reflect chromatin accessibility quality. Commonly used indicators include the fraction of reads in peaks (FRiP), promoter read ratios, blacklist region read ratios, and transcription start site (TSS) enrichment scores ([Buenrostro et al., 2015](https://www.nature.com/articles/nature14590), [Fang et al., 2021](https://www.nature.com/articles/s41467-021-21583-9), [Granja et al., 2021](https://www.nature.com/articles/s41588-021-00790-6)). Barcodes lacking characteristic nucleosome banding patterns—typical of high-quality ATAC-seq data—are also excluded ([Cusanovich et al., 2015](https://www.science.org/doi/epdf/10.1126/science.aab1601?src=getftr&utm_source=sciencedirect_contenthosting&getft_integrator=sciencedirect_contenthosting), [Cusanovich et al., 2018](https://www.cell.com/cell/fulltext/S0092-8674(18)30855-9)). Additionally, peaks located in blacklist regions or overlapping housekeeping genes may be removed during feature filtering ([Fang et al., 2021](https://www.nature.com/articles/s41467-021-21583-9)).

It’s important to note that no universal QC thresholds apply to all datasets. QC criteria should be adapted based on specific sample characteristics, such as data complexity, cellular heterogeneity, expected cell types, batch effects, or sequencing platform used ([Baek and Lee, 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303019#s0015)).

For more information on QC and other sc-ATAC-Seq related methods, see [Lei Xiong et al., 2019](https://www.nature.com/articles/s41467-019-12630-7#Sec10), [Li et al., 2021](https://www.nature.com/articles/s41467-021-26530-2), [Stuart et al., 2021](https://www.nature.com/articles/s41592-021-01282-5#Sec9), 
[Taavitsainen et al, 2021](https://www.nature.com/articles/s41467-021-25624-1#Sec10), [De Rop et al., 2023](https://www.nature.com/articles/s41587-023-01881-x), [Gamache et al., 2023](https://link.springer.com/article/10.1186/s13578-023-01120-5#Sec10)


#### Post alignment/cell quality filtering parameters

We recommend that the user use the following parameters for initial `scATAC` QC, and then adjust accordingly if necessary:

| Parameter                            | Suggested Value | Corrected Comment                                |
| ------------------------------------ | --------------- | ------------------------------------------------ |
| `peak_region_fragments_min_upstream` | `"100"`        | Cells with very few fragments in peaks likely represent background noise, empty droplets, or low-quality nuclei.           |
| `pct_reads_in_peaks_value_upstream`  | `"15"`          | % reads in peaks — **> 15%** (or > 20%)     |
| `blacklist_ratio_value_upstream`     | `"0.05"`        | % reads in ENCODE blacklist — **< 0.05**         |
| `mitochondrial_value_upstream`       | `"5"`           | % mitochondrial reads — **< 5%**                 |
| `TSS.enrichment_value_upstream`      | `"10"`          | TSS enrichment score — **> 2** |
| `nucleosome_signal_value_upstream`   | `"4"`           | Nucleosome signal — **< 4**                      |

           
### (2) Estimating and filtering out doublets

Popular approach of single-cell uses oil droplets or wells to isolate single cells along with barcoded beads. Depending on the cell density loaded, a proportion of reaction volumes (i.e. droplets or wells) will capture more than one cell, forming ‘doublets’ (or ‘multiplets’), i.e. two or more cells captured by a single reaction volume and thus sequenced as a single-cell artifact. 

The proportion of doublets is proportional to the number of cells captured. It is common in single-cell experiments to have 10-20% doublets, making accurate doublet detection critical.

Doublets are prevalent in single-cell sequencing data and can lead to artifactual findings. We will use a computational approach to calculate and remove doublets from the library. Here, we use [ScDblFinder](https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html) method for identifying doublets/multiplets in single-cell data.

The `seurat_obj_raw.rds` object from step (1) is used for this step.
 -  Summary table with doublet metrics and doublets prediction plot are generated.


### (3) Merging filtered data

Next, we merge count matrices from steps (1-2) after filtering out low quality cells and doublets (optional as defined in the `params`). Seurat object and metadata for the library along with UMAP embeddings are saved to be used for downstream analyses.

### (4) Final QC summary report

Lastly, we provide a final QC summary report containing graphs and summary tables across each QC step.

## Folder structure 

The structure of this folder is as follows:

```
├── 01A_run_signac_qc.Rmd
├── 01B_run_signac_qc_multiple_samples.R
├── 02_run_scDblFinder.Rmd
├── 03_run_filter_object.Rmd
├── 04_run_summary_report.Rmd
├── plots
├── lsf-script.txt
├── README.md
├── results
├── run-upstream-analysis.R
├── run-upstream-analysis.sh
└── util
|   ├── function-create-UMAP.R
|   ├── function-process-signac.R
|___└── function-run-QC.R
```

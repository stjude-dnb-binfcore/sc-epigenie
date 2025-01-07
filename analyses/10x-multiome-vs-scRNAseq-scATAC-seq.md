# 10x Genomics Multiome vs. 10x matched scRNA-seq and scATAC-seq

Multiome is a term for using different assays to probe different features in a biological space. The 10x Genomics Multiome kit directly combines scRNA-seq to look at gene expression and scATAC-seq to look at chromatin accessibility in the exact same cell.

The processing of 10x scATAC-seq with matched scRNA-seq versus the 10x multiomic pipeline (which combines both scRNA-seq and scATAC-seq in a single experiment) differs in how the data is generated, processed, and integrated. Both pipelines aim to analyze transcriptomic and epigenomic data at the single-cell level, but they have different workflows and considerations due to the way the data is handled.

# Key Differences

## 1. Experimental Setup

### 10x matched scRNA-seq and scATAC-seq

- In this setup, two separate experiments (or assays) are performed: one for scATAC-seq and one for scRNA-seq.
- These are matched data because some cells from the same sample are used for the scRNA-seq experiment and some for the scATAC-seq experiment.

### 10x Genomics Multiome

- The multiomic approach uses a single experiment where both scRNA-seq and scATAC-seq are captured from the same set of cells in a single run. It uses a combined protocol with unique barcodes for both RNA and ATAC-seq, enabling simultaneous profiling of both transcriptomic and epigenomic features.
- A more integrated experimental design, where both the RNA and ATAC data are collected at the same time, meaning that both types of data share the exact same cell-level barcodes. This allows linking the RNA data to the chromatin accessibility data for the same cell. Cells are processed in two distinct channels and the matching of barcodes (and cells) allows the integration of the data afterward.

## 2. Data Processing Workflow
The pipelines for scATAC-seq with matched scRNA-seq and 10x multiomic data follow different paths in terms of data processing, integration, and analysis.

### 10x matched scRNA-seq and scATAC-seq
This pipeline involves separate processing for each modality (RNA and ATAC), followed by the integration of the results. The steps are:

1. Preprocessing of scRNA-seq. As described in [sc-rna-seq-snap/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/tree/main/analyses).

2. Preprocessing of scATAC-seq. As described in [sc-atac-seq/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-atac-seq/tree/main/analyses).

3. Integration of both assays. As described in [sc-atac-seq/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-atac-seq/tree/main/analyses).

   > - Use integration methods such as Seurat’s canonical correlation analysis (CCA) or Harmony to link RNA and ATAC data and align the datasets in a common space.



### 10x Genomics Multiome
In the multiomic pipeline, both RNA and ATAC data are collected simultaneously, so the workflow is more integrated. The steps are:

1. Preprocessing

   > - Cells are processed in the same experiment, where both RNA and ATAC data are captured together with unique barcodes for both modalities.
   > - The demultiplexing step assigns both the RNA and ATAC data to the same cell barcode.
   > - RNA-seq: Align RNA reads to the transcriptome using Cell Ranger.
   > - ATAC-seq: Align ATAC reads to the genome using Cell Ranger ATAC.
   > - Perform QC on both the RNA and ATAC data. Filtering cells based on both gene counts (for RNA) and peak counts (for ATAC).

2. Feature Generation

   > - For RNA, generate the gene expression matrix based on UMIs.
   > - For ATAC, generate the chromatin accessibility matrix based on peaks or fragments.

3. Normalization and Correction

   > - Normalize both RNA and ATAC data.
   > - RNA normalization: Normalize the UMI count data, often using library size factors or other methods like log-normalization.
   > - ATAC normalization: Normalize the ATAC data for sequencing depth, peak distribution, and cell-specific biases.

4. Integration

   > - Multiomic integration is performed using methods like Seurat (multi-modal integration), which allows both RNA and ATAC data to be jointly analyzed and integrated based on shared cell barcodes.
   > - Dimensionality reduction (e.g., PCA, UMAP, t-SNE) can be done on the multi-modal data.
   > - The analysis also integrates information from both RNA and ATAC-seq to reveal coordinated changes in gene expression and chromatin accessibility.

5. Downstream Analysis

   > - Like the scATAC + scRNA-seq pipeline, clustering can be done based on both transcriptomic and epigenomic features.
   > - However, the major advantage is that the data are already matched at the cell level, so the analysis of joint gene expression and chromatin accessibility is more streamlined and integrated.


## Advantages & Disadvantages

|  | 10x matched scRNA-seq and scATAC-seq  | 10x Genomics Multiome |
|:-----------:|:----------:|:--------:|
| Advantages | Flexible and customizable: You can separately fine-tune RNA and ATAC processing pipelines. | More streamlined: Data are processed together in a single experiment, so the integration of RNA and ATAC data is inherently simpler. |
|  | Allows you to use existing pipelines for scRNA-seq and scATAC-seq separately, which may be advantageous if you need specialized tools for each modality. | Simultaneous profiling: You capture both the transcriptome and the epigenome in a single experiment, reducing potential mismatches or batch effects. |
| Disadvantages | Requires a post-processing integration step, which can be more challenging and may introduce potential mismatches or biases. | Less flexibility: You must use the tools and workflows that are compatible with the 10x multiomic platform (e.g., Cell Ranger and Seurat). |
|  | Separate experiments can sometimes lead to technical inconsistencies. | Requires a large amount of sequencing depth and high-quality data to capture both. |
|  |  | nuclei isolation is mandatory for 10x Multiome because it is a requisite for scATAC-seq’s tagmentation step. This contrasts with scRNA-seq, which can be performed on nuclei and whole cells. You can get an idea of how important the whole-cell transcriptome would be for your experiment in our informative blog on single-nucleus RNA sequencing. A workaround is to combine a standalone whole-cell scRNA-seq experiment with a standalone (single-nuclei) ATAC-seq experiment by dividing the sample for two separate analyses. |
|  |  | Compared to standalone scATAC-seq, 10x Genomics Multiome is currently outperformed in terms of sensitivity and library complexity. In a systematic benchmark study on peripheral blood mononuclear cells (De Rop et al., 2023), 10x Genomics Multiome produced half the unique fragment peaks as the most advanced 10x Single Cell ATAC protocol. |
|  |  | Additional costs while it is less sensitive and efficient in sequencing that standalone scATAC-seq. This has to be taken into account in designs for which scATAC-seq is the primary focus of a study. For these designs, 10x Genomics Single Cell ATAC may be the preferred option. |



# 10x Genomics Multiome: A Guide to research questions, pipeline design, and analyses modules

## Aim

### Deep characterization of cell populations
1. The combined data serves for cross-validation: To identify cell types in a tissue by grouping together nuclei with similar gene expression and chromatin accessibility profiles. 

2. To high-resolution mapping of cell fates in developmental biology and stem cell research. In addition, it is possible to identify transcription factors and regulatory mechanisms responsible for priming cells, as explained in more detail in another section.

3. Unearth novel cell types that are undistinguished by gene expression or chromatin accessibility alone, yet show a unique combination of gene expression and chromatin accessibility profile deviant from other cell types. 

4. How do chromatin accessibility and gene expression change over time? The technology is often applied to generating high-resolution roadmaps of cell fates in samples including: Developing tissues (e.g., Ma et al., 2020; Frazel et al. 2023); Stem cell-based model systems such as organoids (Lee et al., 2023); Immune cell lineages (Chopp et al., 2023); and Transition processes in cancer (Han et al., 2022).

### Map regulatory networks in single cells

5. What regulatory elements are active? To enhance a sample’s single-cell transcription profile with an additional layer of information about each cell’s regulatory elements.

6. Where are transcription factors binding? To identify known and novel sequence patterns occurring within peak regions. These can be binding motifs, short sequences of DNA to which transcription factors bind to regulate gene expression.

7. Which regulatory elements drive gene expression, how?

8. What are the unique gene networks in each cell type?

9. How do conditions influence regulatory elements and downstream pathways?


## Methods

[Signac package](https://stuartlab.org/signac/)

[Signac and Seurat vignette](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)

[Joint RNA and ATAC analysis: 10x multiomic](https://stuartlab.org/signac/articles/pbmc_multiomic)

[Analyzing adult mouse brain scATAC-seq](https://stuartlab.org/signac/articles/mouse_brain_vignette)

[Joint RNA and ATAC analysis: SNARE-seq](https://stuartlab.org/signac/articles/snareseq)


## Analysis modules

### 1. `cellranger-analysis` module (description="Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.", required=True)

#### 10x matched scRNA-seq and scATAC-seq

As described in [sc-rna-seq-snap/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/tree/main/analyses) and [sc-atac-seq/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-atac-seq/tree/main/analyses).

#### 10x Genomics Multiome
- CellRanger-arc count. This step will:
   a. Demultiplex the modalities
   b. Pass demultiplex QC?
   c. Obtain counts with CellRanger
   d. Pass CellRanger QC?

  > Sequencing read alignments of snATAC-seq and snMultiome-seq: To process sequenced snATAC-seq and snMutiome-seq data, we used the CellRanger-atac count (v.2.0, 10x Genomics) and CellRanger-arc count (v.2.0, 10x Genomics) pipelines, respectively. 

### 2. `upstream-analysis` module (description="Pipeline for estimating QC metrics and filtering low quality cells.", required=True)

#### 10x matched scRNA-seq and scATAC-seq

As described in [sc-rna-seq-snap/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/tree/main/analyses) and [sc-atac-seq/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-atac-seq/tree/main/analyses).

#### 10x Genomics Multiome
- Quality control, normalization, feature selection, dimensionality reduction and clustering of snMutiome-seq data
  > For snMultiome-seq data containing profiles of both snRNA- and snATAC-seq data, we first performed separate processing and filtering of cells using the same steps as were described for the processing of separate sc/snRNA-seq and snATAC-seq assays. To obtain the final list of barcodes, we retained the cells that passed the quality control filters in both the snRNA- and snATAC-seq assays. In the result, we obtained filtered gene- and peak-count matrices for the same set of cells. We then performed TF-IDF normalization of the peak-count matrix, followed by LSI dimensionality reduction using the RunTFIDF and RunSVD Signac functions. For normalization and dimensionality reduction of the gene-count matrix, we used the SCTransform and RunPCA functions of Seurat with the same parameters as used for regular sc/snRNA-seq data processing. We next computed the weighted nearest neighbour (WNN) graph with the FindMultiModalNeighbors function using both data modalities. We used 1:30 PCA components from snRNA-seq and 2:30 LSI components from snATAC-seq for this analysis. We performed nonlinear dimensionality reduction of the resulting WNN graph using the RunUMAP function of Seurat. Finally, we obtained clusters with the FindClusters function using the WNN graph, setting the argument algorithm = 3 (SLM).

- Identification of doublets in snATAC-seq and snMultiome-seq samples

- Merging of snATAC-seq data across samples (cohort-level objects)

- Merging of snATAC-seq data across cancers (pan-cancer-level objects)

### 3. `integrative-analysis` module (description="Pipeline for Integrative analysis.", required=True)

#### 10x matched scRNA-seq and scATAC-seq

As described in [sc-rna-seq-snap/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/tree/main/analyses) and [sc-atac-seq/analyses/README.md](https://github.com/stjude-dnb-binfcore/sc-atac-seq/tree/main/analyses).

#### 10x Genomics Multiome
- Inherently integrated since both RNA and ATAC data are generated from the same cells in a single experiment.
- Integration is more seamless because both data types share the same barcodes, allowing for more straightforward downstream analysis (such as joint clustering and multiomic visualization).
- The cell-level correspondence between RNA and ATAC data is ensured from the start.

  > Tools for Integration:
    > - Cell Ranger ATAC and Cell Ranger (for RNA) are typically used for preprocessing, but for integration: Seurat v4 (multi-modal integration) or MOFA are designed to handle multiomic data by leveraging the same barcodes.
    > - Signac: A tool that works seamlessly with Seurat for integrating scRNA-seq and scATAC-seq data in the multiomic pipeline.


### 4. `cell-types-annotation` module (description="Pipeline for annotating cell types.", required=True)

Cell type annotation of snATAC-seq and snMultiome-seq data.

### 5. `peak-calling` module (description="Pipeline for calling peaks and Motif Enrichment Analysis.", required=False)



## Summary of Analyses modules

| Analyses modules | 10x matched scRNA-seq and scATAC-seq  | 10x Genomics Multiome |
|:-----------:|:----------:|:--------:|
| cellranger-analysis | separate for each modality - different pipeline | CellRanger-arc count to align and demultiplex modalities |
| upstream-analysis | separate processing and filtering for each modality - different pipeline - identify doublets | separate processing and filtering for each modality - identify doublets |
| integrative-analysis | Seurat, Harmony, liger for scRNA-seq/Seurat, Harmony for scATAC-seq | Signac |
| cell-types-annotation | e.g., Integrate with scRNA-seq and label transfer | ... |
| peak-calling | yes | ... |
| differential-accessibility-analysis | yes | ... |
| plotting-genomic-regions-analysis | yes | ... |
| gene-ontology-enrichment-analysis | yes | ... |
| trajectory-analysis | yes | ... |


## References

- [10x Genomics Multiome vs. matched scRNA-seq and scATAC-seq](https://www.scdiscoveries.com/blog/10x-genomics-multiome-vs-scrna-seq-and-scatac-seq/).
- [Terekhanova et al., 2023](https://www.nature.com/articles/s41586-023-06682-5#Sec10) have sequenced, processed, and analyzed snATAC-seq and snMutiome-seq data. See also [PanCan_snATAC_publication GitHub repo](https://github.com/ding-lab/PanCan_snATAC_publication).
- Developing tissues (e.g., [Ma et al., 2020](https://www.sciencedirect.com/science/article/pii/S0092867420312538?via%3Dihub#sec4); Frazel et al. 2023)
- Stem cell-based model systems such as organoids by [Lee et al., 2023](https://www.nature.com/articles/s12276-023-01076-z).
- Immune cell lineages by [Chopp et al., 2023](https://www.science.org/doi/10.1126/sciimmunol.adi9066#sec-4).
- Transition processes in cancer by [Han et al., 2022](https://www.sciencedirect.com/science/article/pii/S1535610822005025?via%3Dihub). Code available at the [Single-cell-multi-omics repo](https://github.com/lifei176/Single-cell-multi-omics).
- [Bi et al., 2022](https://www.nature.com/articles/s41467-022-30924-1). Code available at the [Human-primed-to-naive-transition-analysis repo ](https://github.com/zftu/Human-primed-to-naive-transition-analysis/tree/main)
- [Analysis of Single-Cell Multiome ATAC + Gene Expression - Dr. Wayne Doyle](https://www.youtube.com/watch?v=y1mFdkDVc-c&ab_channel=ActiveMotif)



## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-atac-seq/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

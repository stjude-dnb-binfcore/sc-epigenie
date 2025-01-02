# 10x Genomics Multiome vs. matched scRNA-seq and scATAC-seq


The processing of 10x scATAC-seq with matched scRNA-seq versus the 10x multiomic pipeline (which combines both scRNA-seq and scATAC-seq in a single experiment) differs in how the data is generated, processed, and integrated. Both pipelines aim to analyze transcriptomic and epigenomic data at the single-cell level, but they have different workflows and considerations due to the way the data is handled.

# Key Differences Between 10x scATAC + scRNA-seq with Matching vs. 10x Multiomic Pipeline:

## 1. Experimental Setup:

### 10x scATAC + matched scRNA-seq:

- In this setup, two separate experiments (or assays) are performed: one for scATAC-seq and one for scRNA-seq.
- These are matched data because each cell in the scRNA-seq experiment is identified in the scATAC-seq experiment, usually by using the same cell barcode. This allows linking the RNA data to the chromatin accessibility data for the same cell.
- Cells are processed in two distinct channels and the matching of barcodes (and cells) allows the integration of the data afterward.

### 10x Multiomic Pipeline:

- The multiomic approach uses a single experiment where both scRNA-seq and scATAC-seq are captured from the same set of cells in a single run. It uses a combined protocol with unique barcodes for both RNA and ATAC-seq, enabling simultaneous profiling of both transcriptomic and epigenomic features.
- A more integrated experimental design, where both the RNA and ATAC data are collected at the same time, meaning that both types of data share the exact same cell-level barcodes.


## 2. Data Processing Workflow:
The pipelines for scATAC-seq with matched scRNA-seq and 10x multiomic data follow different paths in terms of data processing, integration, and analysis.

### scATAC-seq with Matched scRNA-seq Pipeline:
This pipeline involves separate processing for each modality (RNA and ATAC), followed by the integration of the results. The steps are:

1. Preprocessing of scRNA-seq:

- Demultiplex and assign RNA reads to cells based on the cell barcodes.
- Perform quality control (QC) on the RNA data (e.g., filtering low-quality cells based on the number of genes detected, number of UMIs, mitochondrial gene expression).
- Normalize RNA counts (e.g., TPM, CPM, or using tools like Seurat).

2. Preprocessing of scATAC-seq:

- Align ATAC-seq reads to the reference genome (e.g., Bowtie2 or BWA).
- Identify peaks (regions of open chromatin) using tools like MACS2 or HOMER.
- Count fragments in these peaks for each cell.
- Normalize ATAC-seq data for sequencing depth and other biases.

3. Integration:

- After preprocessing, integrate the scRNA-seq and scATAC-seq data.
- This is typically done using cell barcodes to match RNA and ATAC data at the single-cell level.
- Use integration methods such as Seurat’s canonical correlation analysis (CCA) or Harmony to link RNA and ATAC data and align the datasets in a common space.

4. Analysis:

- Dimensionality reduction (PCA, UMAP) is performed separately for the RNA and ATAC data, and then integrated.
- Identify cell-type clusters based on both transcriptomic and epigenomic features.
- Analyze the relationship between gene expression and chromatin accessibility. For example, look for cis-regulatory elements that correlate with gene expression changes.

### 10x Multiomic Pipeline:
In the multiomic pipeline, both RNA and ATAC data are collected simultaneously, so the workflow is more integrated. The steps are:

1. Preprocessing:

- Cells are processed in the same experiment, where both RNA and ATAC data are captured together with unique barcodes for both modalities.
- The demultiplexing step assigns both the RNA and ATAC data to the same cell barcode.
- Perform QC on both the RNA and ATAC data. Filtering cells based on both gene counts (for RNA) and peak counts (for ATAC).

2. Alignment:

- RNA-seq: Align RNA reads to the transcriptome using Cell Ranger.
- ATAC-seq: Align ATAC reads to the genome using Cell Ranger ATAC.

3.Feature Generation:

- For RNA, generate the gene expression matrix based on UMIs.
- For ATAC, generate the chromatin accessibility matrix based on peaks or fragments.

4. Normalization and Correction:

- Normalize both RNA and ATAC data.
- RNA normalization: Normalize the UMI count data, often using library size factors or other methods like log-normalization.
- ATAC normalization: Normalize the ATAC data for sequencing depth, peak distribution, and cell-specific biases.

5. Integration:

- Multiomic integration is performed using methods like Seurat v3/Seurat v4 (multi-modal integration), which allows both RNA and ATAC data to be jointly analyzed and integrated based on shared cell barcodes.
- Dimensionality reduction (e.g., PCA, UMAP, t-SNE) can be done on the multi-modal data.
- The analysis also integrates information from both RNA and ATAC-seq to reveal coordinated changes in gene expression and chromatin accessibility.

6. Analysis:

- Like the scATAC + scRNA-seq pipeline, clustering can be done based on both transcriptomic and epigenomic features.
- However, the major advantage is that the data are already matched at the cell level, so the analysis of joint gene expression and chromatin accessibility is more streamlined and integrated.


## Integration of Data:
### 10x scATAC-seq with matched scRNA-seq:
- Requires post-processing integration since the RNA and ATAC data come from separate experiments, even if matched at the cell level.
- The integration is usually performed after the individual data processing, and tools like Seurat’s CCA or Harmony can be used to align the RNA and ATAC datasets.
- Integration is not as straightforward as in the multiomic pipeline, since the two datasets are processed independently before integration.

### 10x Multiomic Pipeline:
- Inherently integrated since both RNA and ATAC data are generated from the same cells in a single experiment.
- Integration is more seamless because both data types share the same barcodes, allowing for more straightforward downstream analysis (such as joint clustering and multiomic visualization).
- The cell-level correspondence between RNA and ATAC data is ensured from the start.

## Tools for Integration:

### scATAC + matched scRNA-seq:

- Seurat: Seurat can integrate the two datasets via CCA (Canonical Correlation Analysis) or more recent multi-modal integration techniques (for version 4+).
- Harmony: Another integration tool that is commonly used for aligning datasets after normalization.

### 10x Multiomic Pipeline:

- Cell Ranger ATAC and Cell Ranger (for RNA) are typically used for preprocessing, but for integration:
Seurat v4 (multi-modal integration) or MOFA are designed to handle multiomic data by leveraging the same barcodes.
- Signac: A tool that works seamlessly with Seurat for integrating scRNA-seq and scATAC-seq data in the multiomic pipeline.


## Advantages & Disadvantages:

### scATAC-seq with Matched scRNA-seq:
- Advantages:
   - Flexible and customizable: You can separately fine-tune RNA and ATAC processing pipelines.
   - Allows you to use existing pipelines for scRNA-seq and scATAC-seq separately, which may be advantageous if you need specialized tools for each modality.
- Disadvantages:
   - Requires a post-processing integration step, which can be more challenging and may introduce potential mismatches or biases.
   - Separate experiments can sometimes lead to technical inconsistencies.

### 10x Multiomic Pipeline:
- Advantages:
   - More streamlined: Data are processed together in a single experiment, so the integration of RNA and ATAC data is inherently simpler.
   - Simultaneous profiling: You capture both the transcriptome and the epigenome in a single experiment, reducing potential mismatches or batch effects.
- Disadvantages:
    - Less flexibility: You must use the tools and workflows that are compatible with the 10x multiomic platform (e.g., Cell Ranger and Seurat).
    - Requires a large amount of sequencing depth and high-quality data to capture both.











---------------------------------------------------------------------------------------
**How 10x Multiome Works: A Guide**


**Aim**

For more information, please see [10x Genomics Multiome vs. matched scRNA-seq and scATAC-seq](https://www.scdiscoveries.com/blog/10x-genomics-multiome-vs-scrna-seq-and-scatac-seq/).

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



**Methods**



**References**

- [10x Genomics Multiome vs. matched scRNA-seq and scATAC-seq](https://www.scdiscoveries.com/blog/10x-genomics-multiome-vs-scrna-seq-and-scatac-seq/)


---------------------------------------------------------------------------------------

**Analysis modules**




## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-atac-seq/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

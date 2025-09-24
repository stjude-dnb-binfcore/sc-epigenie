###########################################################################
#' Function to run miQC
#' @param seurat_obj
#'
#' @return
#' @export
#'
#' @examples
create_qc_metrics <- function(seurat_obj) {
  
  # Compute nucleosome signal score per cell: calculates a ratio of the fragments that are between 147-294bp and 
  # those that are <147 because the average length to wrap around one histone is 147bp. So, this way we can calculate mono-nuclesomes and fragments that are nucleosome-free.
  seurat_obj <- NucleosomeSignal(object = seurat_obj)
  
  # Compute transcription start sites based on the gene location - TSS enrichment score per cell
  # The TSSPlot() function requires per-base enrichment profiles, which are only generated when `fast = FALSE`.
  # By default, fast = TRUE, which skips computing the enrichment matrix for speed (but also skips what TSSPlot() needs).
  seurat_obj <- TSSEnrichment(object = seurat_obj, fast = FALSE)
  
  # Blacklist region fragments: the number of fragments that map to Blacklist regions which are junk regions that have been annotated, 
  # they have high signal in different NGS experiments independent of cell type or experiment and so these are regions we know they are no good.
  # Let's calculate and add the `blacklist_ratio` column
  seurat_obj$blacklist_ratio <- seurat_obj$blacklist_region_fragments / seurat_obj$peak_region_fragments
  
  # Let's add fraction of reads in peaks
  seurat_obj$pct_reads_in_peaks <- seurat_obj$peak_region_fragments / seurat_obj$passed_filters * 100
  
  # https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/singlecell
  seurat_obj$pct_reads_in_peaks_promoters <- seurat_obj$promoter_region_fragments / seurat_obj$passed_filters * 100
  seurat_obj$pct_reads_in_peaks_enhancers <- seurat_obj$enhancer_region_fragments / seurat_obj$passed_filters * 100
  
  # Add more QC
  seurat_obj$high.tss <- ifelse(seurat_obj$TSS.enrichment > 2, 'High', 'Low')
  seurat_obj$nucleosome_group <- ifelse(seurat_obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  seurat_obj$nucleosome_group_category <- ifelse(seurat_obj$nucleosome_signal > 4, "High NS", "Low NS")
  
  return(seurat_obj)
}


###########################################################################
#' Function to run QC with default parameters
#' @param seurat_obj
#'
#' @return
#' @export
#'
#' @examples
run_QC <- function(seurat_obj) {
  if (use_threshold_filtering == "YES"){
    
    message("🔍 We will subset based on threshold values defined in the YAML.")
    
    seurat_obj <- subset(x = seurat_obj,
                         subset = 
                           #nCount_peaks > nCount_peaks_min &
                           #nCount_peaks < nCount_peaks_max &
                           #nFeature_peaks > nFeature_peaks_min &
                           #nFeature_peaks < nFeature_peaks_max &
                           #peak_region_fragments > peak_region_fragments_min &
                           #peak_region_fragments < peak_region_fragments_max &
                           pct_reads_in_peaks > pct_reads_in_peaks_value &
                           blacklist_ratio < blacklist_ratio_value &
                           #duplicate < duplicate_value &
                           #mitochondrial < mitochondrial_value &
                           nucleosome_signal < nucleosome_signal_value &
                           TSS.enrichment > TSS.enrichment_value)
    message("✅ Filtering based on fixed thresholds complete.")
    
  } else {
    
    message("🔍 We will subset based on percentile for filtering. That means that we will set thresholds via the quantile function that differentiates between bottom 2% and everything else.")
    
    # Calculate thresholds using quantiles (without print)
    # https://github.com/mousepixels/sanbomics_scripts/blob/main/scATAC_intro_R.Rmd
    low_prf <- print(quantile(seurat_obj[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02))
    hig_prf <- print(quantile(seurat_obj[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98))
    low_prp <- print(quantile(seurat_obj[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02))
    high_blr <- print(quantile(seurat_obj[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98))
    hig_ns <- print(quantile(seurat_obj[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98))
    low_ts <- print(quantile(seurat_obj[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02))
    
    seurat_obj <- subset(x = seurat_obj,
                         subset = peak_region_fragments > low_prf &
                           peak_region_fragments < hig_prf &
                           pct_reads_in_peaks > low_prp &
                           blacklist_ratio < high_blr &
                           nucleosome_signal < hig_ns &
                           TSS.enrichment > low_ts)
    message("✅ Filtering based on percentiles complete.")
    
  } 
  
  return(seurat_obj)
}



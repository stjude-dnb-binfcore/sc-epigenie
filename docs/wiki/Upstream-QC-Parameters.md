# Upstream QC Parameters

Guidance for **`./analyses/upstream-analysis`** QC settings in `project_parameters.Config.yaml`, based on parameter benchmarking of an internal mm10 scATAC-seq validation cohort (8 samples, unpublished data; 2026).

**Full validation report:** [upstream-qc-parameter-benchmark.md](https://github.com/stjude-dnb-binfcore/sc-epigenie/blob/main/docs/validation/upstream-qc-parameter-benchmark.md)

---

## Recommended defaults

Use these settings for a typical good-quality cohort. Do not change unless you have a specific reason (see decision tree below).

| Parameter | Default | Notes |
|-----------|---------|-------|
| `min.features_value_module` | `50` | Barcodes with fewer peaks are removed at object creation |
| `use_threshold_filtering_upstream` | `"YES"` | Fixed cutoffs (recommended over percentile filtering) |
| `pct_reads_in_peaks_value_upstream` | `20` | Lowering to 15 has negligible effect |
| `TSS.enrichment_value_upstream` | `4` | Do not lower below 4 without reviewing QC plots |
| `blacklist_ratio_value_upstream` | `0.05` | |
| `nucleosome_signal_value_upstream` | `3` | |
| `min.cutoff_value_upstream` | `"q0"` | See feature-selection note below |

---

## How `run_QC` works

QC is defined in `analyses/upstream-analysis/util/function-run-QC.R` and applied **per sample** before merging.

1. **`create_qc_metrics()`** — computes nucleosome signal, TSS enrichment, `blacklist_ratio`, `pct_reads_in_peaks`, etc.
2. **`run_QC()`** — filters cells using either:
   - **Threshold mode** (`use_threshold_filtering_upstream: "YES"`) — fixed YAML cutoffs on 4 metrics: `pct_reads_in_peaks`, `blacklist_ratio`, `nucleosome_signal`, `TSS.enrichment`
   - **Percentile mode** (`"NO"`) — per-sample 2nd / 98th quantile cutoffs on those 4 metrics **plus** `peak_region_fragments`

**Not in `run_QC`:** `min.features_value_module` (object creation) and `min.cutoff_value_upstream` (FindTopFeatures for DR — affects clustering/UMAP, not cell counts).

Percentile thresholds adapt to each sample's distribution, so poor-quality samples can pass with lower absolute QC than fixed defaults allow. Prefer threshold mode when samples should meet the same standards.

---

## Decision tree

```
Start with defaults above
        │
        ▼
Review Cell Ranger + upstream summary report
        │
        ├─ All samples look good ──────────────────► Keep defaults
        │
        ├─ One sample clearly low quality ───────────► Exclude sample OR set min.features = 200
        │                                              for that sample's impact; review UMAP
        │
        ├─ Need maximum cell recovery ───────────────► Set use_threshold_filtering = "NO"
        │                                              (percentile filter; retains ~30% more cells
        │                                               but lower median QC — use with caution)
        │
        ├─ Integration UMAP very noisy ──────────────► Try min.cutoff_value_upstream = "q75"
        │                                              (top 25% peaks for DR; does NOT change
        │                                               cell counts; affects clustering/UMAP)
        │
        └─ Unsure ───────────────────────────────────► Contact DNB Bioinformatics Core
```

---

## Quick reference: what matters

| Change | Cell count impact | Recommendation |
|--------|-------------------|----------------|
| Threshold → percentile filtering | **Large** (+~11k cells on 8-sample cohort) | Avoid as default; too lenient for low-quality samples |
| `min.features` 50 → 200 | **Large only for bad samples** | Use 200 when Cell Ranger/QC suggests poor sample quality |
| `TSS` 4 → 3 | Moderate (+~3.7k cells) | Not recommended; mostly rescues low-quality cells |
| `TSS` 4 → 2 | Small (+~1k cells) | **Do not use** |
| `% reads in peaks` 20 → 15 | Negligible (+5 cells) | Either is fine |
| `min.cutoff` q0 → q75/q95 | **None** (cell counts unchanged) | q75 may improve clustering; q95 too stringent |

---

## Important distinctions

**Cell filtering vs feature selection**

- Parameters like `TSS.enrichment`, `pct_reads_in_peaks`, and `use_threshold_filtering` control **which cells are kept**.
- `min.cutoff_value_upstream` controls **which peaks are used** for dimensionality reduction (FindTopFeatures). It does **not** change cell counts.

**Outlier samples**

Benchmark results were strongly influenced by one low-quality sample (`ATAC_7f2a`). For cohorts without outlier samples, `min.features = 50` vs `200` makes almost no difference. Always review per-sample Cell Ranger metrics before changing defaults.

> **St Jude internal:** Anonymous labels map to original sample IDs in `benchmarking-sc-atac-seq/comparison/SAMPLE_ID_MAPPING.md` on HPCF (do not publish).

---

## Pipeline notes

- **FindClusters (>46k cells):** If using percentile filtering on large cohorts, ensure `method = "igraph"` is set in `integrative-analysis/util/function-samples-integrate.R`.
- **Rerun scope:** After changing upstream parameters, rerun `upstream-analysis` and `integrative-analysis` at minimum.

---

Note: Benchmarking is currently ongoing. The recommendations and validated default parameters presented here are based on results from a single cohort and should be considered preliminary. Additional datasets and cohorts will be evaluated to further validate and refine these recommendations.

---

## References

- [Validation benchmark (repo doc)](https://github.com/stjude-dnb-binfcore/sc-epigenie/blob/main/docs/validation/upstream-qc-parameter-benchmark.md)
- [sc-epigenie README — Tutorial and Documentation](https://github.com/stjude-dnb-binfcore/sc-epigenie#tutorial-and-documentation)
- Benchmark data (St Jude internal): `common/pipeline_development_shared/.../benchmarking-sc-atac-seq/comparison/`

*Maintainers: DNB Bioinformatics Core. Last updated: 2026-08-28.*

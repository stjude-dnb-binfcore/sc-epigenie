# To replicate the clean environment behavior of R Markdown in an R script, consider clearing the workspace at the start of your script by calling
# rm(list = ls())

# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")

suppressPackageStartupMessages({
  
# Set up
  # 01A_run_seurat_qc.Rmd
  library(future)
  library(cowplot)
  library(devtools)
  library(forcats)
  library(GGally)
  library(stringr)
  library(ggpmisc)
  library(ggrepel)
  library(miQC) 
  library(flexmix) # to estimate mixtureModel for miQC
  library(scater) 
  library(Seurat) 
  library(SingleCellExperiment)
  library(irlba) # this solves the issue with RunUMAP code chunk
  library(scooter)
  library(tidyverse)
  library(fs) #file system
  library(RColorBrewer)
  
  # 02_run_SoupX.Rmd
  library(future)
  library(knitr)
  library(SoupX)
  library(Seurat)
  library(stringr)
  library(tidyverse)
  library(tinytex)
  library(hdf5r)
  

  # 03_run_scDblFinder.Rmd
  library(scDblFinder)
  library(Seurat)
  library(scater)
  library(future)
  library(tidyverse)
  library(grid)
  library(knitr)

  # 04_run_filter_object.Rmd
  library(devtools)
  library(future)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(ggthemes)
  library(scooter)
  library(RColorBrewer)
  library(knitr)

  # 05_run_summary_report.Rmd
  library(tidyverse)  
  library(knitr)
  library(patchwork)

  # 01-integrative-analysis.Rmd
  library(future)
  library(tidyverse)
  library(patchwork)
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(rliger)
  library(RcppPlanc)
  library(SeuratWrappers)
  library(scooter)
  library(reshape2)
  library(RColorBrewer)
  library(knitr)
  
  # run-cell-types-annotation.R
  library(yaml)
  library(tidyverse)
  library(celldex)
  
  # 01-cell-types-annotation-SingleR-broad.Rmd
  library(tidyverse)
  library(Seurat)
  library(SingleR)
  library(scooter)
  library(knitr)

  library(ShinyCell)
  library(shinyhelper)
  library(DT)
  library(ggdendro)

  
  library(scDblFinder)
  library(celldex)
  library(infercnv)
  library(numbat)
  
  ####################################################
  # sc-epigenie
  library(AnnotationHub)
  library(ensembldb)
  library(GenomeInfoDb) #translation between chromosome names
  
  library(GenomicRanges)
  library(rtracklayer)
  
  library(Matrix)
  library(htmlwidgets)
  
  library(enrichplot)
  
  library(org.Hs.eg.db)
  library(org.Mm.eg.db) # Mouse gene annotation (use org.Hs.eg.db for human)
  
  # Motif analysis with Signac
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Mmusculus.UCSC.mm10)
  
  # Transcription factor footprinting
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  library(monocle3)
  library(cicero)
  library(biovizBase)
  library(clusterProfiler)
  ####################################################
})
#############################################################################################
#' Function to process and cluster cell for an individual sample
#' @param seurat_obj
#' @param min.cutoff_value
#' @param assay
#' @param results_dir
#' @param plots_output
#' @param use_condition_split
#' @param condition1
#' @param condition2
#' @param condition3
#' @param print_pdf
#' @param grouping
#'
#' @return
#' @export
#'
#' @examples
#' 
process_signac <- function(seurat_obj, min.cutoff_value, assay, results_dir, 
                           plots_output, use_condition_split, condition1, condition2, 
                           condition3, print_pdf, grouping) {
  
  set.seed(1234) # Make code reproducible
  
  # Normalize data
  seurat_obj <- RunTFIDF(seurat_obj)
  seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = min.cutoff_value) 
  seurat_obj <- RunSVD(seurat_obj) # SVD applied to term-document matrix, i.e., peak-call count matrix afer TF-IDF transformation; analogous to PCA
  
  # Plot DepthCor
  #name <- paste0(plots_dir, "/", glue::glue("DepthCor.png"))
  #print(DepthCor(seurat_obj))
  #ggsave(file = name, width = 6, height = 5, device = "png")
  
  name <-  paste0(plots_dir, "/", glue::glue("DepthCor.pdf"))
  p <- DepthCor(seurat_obj)
  pdf(file = name, width = 6, height = 5)
  print(p)
  dev.off()
  
  cat("Run UMAP", "\n")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "lsi", dims = 2:30)
  #seurat_obj <- FindNeighbors(object = seurat_obj, reduction = 'lsi', dims = 2:30)
  #seurat_obj <- FindClusters(object = seurat_obj, verbose = FALSE, algorithm = 3)
  
  # Plot DimPlot
  #name <- paste0(plots_dir, "/", glue::glue("DimPlot.png"))
  #print(DimPlot(object = seurat_obj, label = TRUE) + NoLegend())
  #ggsave(file = name, width = 6, height = 5, device = "png")
  
  
  # Plot features
  #name <- paste0(plots_dir, "/", glue::glue("FeaturePlot.png"))
  #p1 <- FeaturePlot(seurat_obj, features = "nCount_peaks") + ggtitle("nCount_peaks")
  #p2 <- FeaturePlot(seurat_obj, features = "nFeature_peaks") + ggtitle("nFeature_peaks")
  #print(p1 + p2 + plot_layout(ncol=2))
  #ggsave(file = name, width = 8, height = 5, device = "png")
  
  # Save metadata
  reduction_names <- c(paste0("umap"), paste0("lsi")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  write_tsv(metadata, file = paste0(results_dir, "/", "metadata", ".tsv")) # Save metadata
  
  saveRDS(seurat_obj, file = paste0(results_dir, "/", "seurat_obj.rds")) # Save Seurat object
  
  cat("Plot UMAP", "\n")
  umap_values <- c(paste0("UMAP"))
  
  if (print_pdf == "YES") {
    for (umap_val in umap_values) {
      print(umap_val)
      set.seed(1234) # Make code reproducible for UMAPs
    
      # Save plots
      umap_out <- file.path(paste0(plots_output, "/", "UMAP", "/"))
      dir.create(umap_out)
      print(umap_out)
    
      # nFeature_peaks
      name <- paste0(umap_out, "01_", "nFeature_peaks.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$nFeature_peaks,
                                      palette = gradient_palette_df,
                                      title_name = "nFeature_peaks")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # nCount_peaks
      name <- paste0(umap_out, "02_", "nCount_peaks.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$nCount_peaks,
                                      palette = gradient_palette_df,
                                      title_name = "nCount_peaks")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # mitochondrial
      name <- paste0(umap_out, "03_", "mitochondrial.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$mitochondrial,
                                      palette = gradient_palette_df,
                                      title_name = "mitochondrial")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # orig.ident
      name <- paste0(umap_out, "04_",  "orig.ident.pdf")
      p <- create_UMAP_orig_ident(df = metadata,
                                  umap_val = umap_val,
                                  color_value = "orig.ident",
                                  title_name = "orig.ident")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # ID
      name <- paste0(umap_out, "05_",  "ID.pdf")
      p <- create_UMAP_ID(df = metadata,
                          umap_val = umap_val,
                          color_value = "ID",
                          title_name = "ID")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
     # condition1
     name <- paste0(umap_out, "06_", "condition1.pdf")
     p <- create_UMAP_condition(df = metadata,
                                umap_val = umap_val,
                                color_value = condition1)
     pdf(file = name, width = 6, height = 5)
     print(p)
     dev.off()
     
     # condition2
     if (!is.null(condition2)) {
       name <- paste0(umap_out, "07_", "condition2.pdf")
       p <- create_UMAP_condition(df = metadata,
                                  umap_val = umap_val,
                                  color_value = condition2)
       pdf(file = name, width = 6, height = 5)
       print(p)
       dev.off() 
     }
     
     # condition3
     if (!is.null(condition3)) {
       name <- paste0(umap_out, "08_", "condition3.pdf")
       p <- create_UMAP_condition(df = metadata,
                                  umap_val = umap_val,
                                  color_value = condition3)
       pdf(file = name, width = 6, height = 5)
       print(p)
       dev.off()
     }
     
     if (use_condition_split == "YES") {
       print("Print condition_split!")
       # condition_split1
       name <- paste0(umap_out, "09_", "condition_split1.pdf")
       p <- create_UMAP_condition_split(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata[[condition1]])
       pdf(file = name, width = 10, height = 6)
       print(p)
       dev.off()
       
       # condition_split2
       if (!is.null(condition2)) {
         name <- paste0(umap_out, "10_", "condition_split2.pdf")
         p <- create_UMAP_condition_split(df = metadata,
                                          umap_val = umap_val,
                                          color_value = metadata[[condition2]])
         pdf(file = name, width = 10, height = 6)
         print(p)
         dev.off()
       }
       
       # condition_split3
       if (!is.null(condition3)) {
         name <- paste0(umap_out, "11_", "condition_split3.pdf")
         p <- create_UMAP_condition_split(df = metadata,
                                          umap_val = umap_val,
                                          color_value = metadata[[condition3]])
         pdf(file = name, width = 10, height = 6)
         print(p)
         dev.off()
       }
       
       
       } else { 
         print("Do NOT Print condition_split!")
         next }
     }
    } else if (print_pdf == "NO") {
      for (umap_val in umap_values) {
        print(umap_val)
        set.seed(1234) # Make code reproducible for UMAPs
    
        # Save plots
        umap_out <- file.path(paste0(plots_output, "/", "UMAP", "/"))
        dir.create(umap_out)
        print(umap_out)
    
        # nFeature_peaks
        name <- paste0(umap_out, "01_", "nFeature_peaks.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$nFeature_peaks,
                                        palette = gradient_palette_df,
                                        title_name = "nFeature_peaks")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # nCount_peaks
        name <- paste0(umap_out, "02_", "nCount_peaks.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$nCount_peaks,
                                        palette = gradient_palette_df,
                                        title_name = "nCount_peaks")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # mitochondrial
        name <- paste0(umap_out, "03_", "mitochondrial.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$mitochondrial,
                                        palette = gradient_palette_df,
                                        title_name = "mitochondrial")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # orig.ident
        name <- paste0(umap_out, "04_", "orig.ident.png")
        p <- create_UMAP_orig_ident(df = metadata,
                                    umap_val = umap_val,
                                    color_value = "orig.ident",
                                    title_name = "orig.ident")
        ggsave(file = name, width = 6, height = 5, device = "png")
        
        
        # ID
        name <- paste0(umap_out, "05_", "ID.png")
        p <- create_UMAP_ID(df = metadata,
                            umap_val = umap_val,
                            color_value = "ID",
                            title_name = "ID")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # condition1
        name <- paste0(umap_out, "06_", "condition1.png")
        p <- create_UMAP_condition(df = metadata,
                                   umap_val = umap_val,
                                   color_value = condition1)
        ggsave(file = name, width = 6, height = 5, device = "png")

        
        # condition2
        if (!is.null(condition2)) {
          name <- paste0(umap_out, "07_", "condition2.png")
          p <- create_UMAP_condition(df = metadata,
                                     umap_val = umap_val,
                                     color_value = condition2)
          ggsave(file = name, width = 6, height = 5, device = "png")
          
        }
        
        # condition3
        if (!is.null(condition3)) {
          name <- paste0(umap_out, "08_", "condition3.png")
          p <- create_UMAP_condition(df = metadata,
                                     umap_val = umap_val,
                                     color_value = condition3)
          ggsave(file = name, width = 6, height = 5, device = "png")
          
        }

        if (use_condition_split == "YES"){
          print("Print condition_split!")
          # condition_split1
          name <- paste0(umap_out, "09_", "condition_split1.png")
          p <- create_UMAP_condition_split(df = metadata,
                                           umap_val = umap_val,
                                           color_value = condition1)
          ggsave(file = name, width = 10, height = 6, device = "png")
          
          # condition_split2
          if (!is.null(condition2)) {
            name <- paste0(umap_out, "10_", "condition_split2.png")
            p <- create_UMAP_condition_split(df = metadata,
                                             umap_val = umap_val,
                                             color_value = condition2)
            ggsave(file = name, width = 10, height = 6, device = "png")
          }
          
          # condition_split3
          if (!is.null(condition3)) {
            name <- paste0(umap_out, "11_", "condition_split3.png")
            p <- create_UMAP_condition_split(df = metadata,
                                             umap_val = umap_val,
                                             color_value = condition3)
            ggsave(file = name, width = 10, height = 6, device = "png")
          }
          
          
          } else { 
            print("Do NOT Print condition_split!")
            next }
    }
  }
  
  cat("Returned Seurat object", "\n")
  return(seurat_obj)
}



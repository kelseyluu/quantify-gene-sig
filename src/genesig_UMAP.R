genesig_UMAP <- function(seurat_obj, 
                         genesig_name,
                         group_colname,
                         img_format,
                         fontsize=16, 
                         percentile_thresh=0.98) {
  
  # UMAP by gene sig activity classification 
  DimPlot(seurat_obj, 
          reduction = "umap", 
          group.by = paste0(genesig_name, "_active"), 
          order = "active", 
          cols = c("gray", "red"), 
          pt.size = 0.4,
          split.by = group_colname) + 
     ggtitle("Gene signature activity UMAP") +
    theme_classic(base_size = fontsize) 
  

ggsave(paste0('./output/genesig_UMAP.', img_format), width = 12) %>% suppressMessages()

}
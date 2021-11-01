genesig_UMAP <- function(seurat_obj, 
                         genesig_name,
                         group_colname,
                         img_format) {
  
  # UMAP by gene sig activity classification 
  DimPlot(seurat_obj, 
          reduction = "umap", 
          group.by = paste0(genesig_name, "_active"), 
          order = "active", 
          cols = c("gray", "red"), 
          pt.size = 1,
          split.by = group_colname) + 
     ggtitle("Gene signature activity UMAP") +
    theme_classic(base_size = 16) 
  

ggsave(paste0('./output/genesig_UMAP.', img_format), width = 12) %>% suppressMessages()

}
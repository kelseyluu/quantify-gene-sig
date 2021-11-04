genesig_UMAP <- function(seurat_obj, 
                         genesig_name,
                         group_colname,
                         img_format,
                         out_dir) {
  
  # UMAP by gene sig activity classification 
  DimPlot(seurat_obj, 
          reduction = "umap", 
          group.by = paste0(genesig_name, "_active"), 
          order = "active", 
          cols = c("gray", "red"), 
          pt.size = 1,
          split.by = group_colname,
          ncol = 2) + 
     ggtitle("Gene signature activity UMAP") +
    theme_classic(base_size = 16) 
  

dir.create(here(out_dir), showWarnings = FALSE)    
ggsave(here(out_dir, paste0(genesig_name, '_UMAP.', img_format)), 
        width = 12, height = 6*ceiling(n_groups/2)) %>% suppressMessages()

}
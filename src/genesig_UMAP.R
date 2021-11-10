genesig_UMAP <- function(seurat_obj, 
                         genesig_name,
                         group_colname,
                         cluster_colname,
                         img_format,
                         out_dir) {

    # UMAP by gene sig activity classification 
    p <- DimPlot(seurat_obj, 
            reduction = "umap", 
            group.by = paste0(genesig_name, "_active"), 
            order = "active", 
            cols = c("gray", "red"), 
            pt.size = 1,
            split.by = group_colname,
            ncol = 2) +     
            ggtitle("Gene signature activity UMAP") +
            theme_classic(base_size = 16) 

    if (!is.null(cluster_colname)) {
        # get cluster center coordinates
        tmp <- DimPlot(seurat_obj, 
            reduction = "umap", 
            group.by = cluster_colname)

        umap_avg_df <- tmp$data %>% 
            group_by(cluster = get(cluster_colname)) %>% 
            summarise(UMAP1_avg=mean(UMAP_1), UMAP2_avg=mean(UMAP_2))
    
        p <- p + ggrepel::geom_text_repel(data = umap_avg_df,
                            aes(x=UMAP1_avg, 
                            y=UMAP2_avg, 
                            label=cluster),
                            size=3)
    }

    p

    dir.create(here(out_dir), showWarnings = FALSE)    
    ggsave(here(out_dir, paste0(genesig_name, '_UMAP.', img_format)), 
            width = 12, height = 6*ceiling(n_groups/2)) %>% suppressMessages()

}
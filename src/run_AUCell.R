run_AUCell <- function(seurat_obj, 
                       geneSet, 
                       cells_rankings,
                       save_data,
                       out_dir) {
                         
  genesig_name <- names(geneSet)
  
  # rank genes by expression in each cell:
  if (is.null(cells_rankings)) {
    cells_rankings <- AUCell_buildRankings(seurat_obj %>% GetAssayData(), plotStats = FALSE) 
  } 
  
  # calculate gene set scores
  cells_AUC <- AUCell_calcAUC(geneSet, cells_rankings, 
                              aucMaxRank=nrow(cells_rankings)*0.05) 
  
  # find best activity classification threshold
  set.seed(123)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                               nCores=1, assign=TRUE,
                                               plotHist=FALSE)
  thresh <- cells_assignment[[genesig_name]]$aucThr$selected %>% round(4)
  active_cells <- cells_assignment[[genesig_name]]$assignment
  cat(paste('Optimal gene signature threshold:', thresh, '\n'))
  cat(paste0('Cells passing threshold: ', length(active_cells), ' (', round(100 * length(active_cells) / length(Cells(seurat_obj)), 2), '%)\n'))

  # get scores per cell
  genesig_scores <- getAUC(cells_AUC)
  
  # add genesig scores to seurat metadata
  sig_name <- paste0(genesig_name, "_sig_score")
  seurat_obj <- AddMetaData(seurat_obj, 
                            t(genesig_scores),
                            col.name = sig_name)
  
  # add genesig binary active status to seurat metadata
  seurat_obj <- AddMetaData(seurat_obj, 
                            metadata = factor(ifelse(
                              Cells(seurat_obj) %in% active_cells,
                              "active", "inactive"), 
                              levels=c('active', 'inactive')),
                            col.name = paste0(genesig_name, "_active"))
  
  if (save_data) {
    dir.create(here(out_dir), showWarnings = FALSE)
    seurat_obj@meta.data %>% 
    select(all_of(sig_name)) %>% 
    write.csv(here(out_dir, paste0(sig_name, 's.csv')))
  }


  return(seurat_obj)
}
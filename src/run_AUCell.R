run_AUCell <- function(seurat_obj, geneSet, cells_rankings=NULL) {
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
  cat(paste('Optimal threshold:', cells_assignment[[genesig_name]]$aucThr$selected %>% round(4), '\n'))

  # get scores per cell
  genesig_scores <- getAUC(cells_AUC)
  
  # add genesig scores to seurat metadata
  seurat_obj <- AddMetaData(seurat_obj, 
                            t(genesig_scores),
                            col.name = paste0(genesig_name, "_sig_score"))
  
  # add genesig binary active status to seurat metadata
  seurat_obj <- AddMetaData(seurat_obj, 
                            metadata = ifelse(
                              Cells(seurat_obj) %in% cells_assignment[[genesig_name]]$assignment,
                              "active", "inactive") %>% as.factor,
                            col.name = paste0(genesig_name, "_active"))
  
  return(seurat_obj)
}
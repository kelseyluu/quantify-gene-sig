plot_ROC <- function(metadata, 
                             genesig_name,
                             id_var, 
                             response_var, 
                             response_label,
                             cell_level,
                             img_format,
                             out_dir,
                             save_data,
                             thresh
                             ) { 

    genesig_var <- paste0(genesig_name, "_sig_score")
    response_groups <- unique(metadata[[response_var]]) %>% as.character()
    nonresp_label <- response_groups[response_groups != response_label]

    # get patient level response table
    patient_summary <- metadata %>% group_by(patient_id=get(id_var)) %>%
    summarise(response = unique(get(response_var)),
            avg_genesig = mean(get(genesig_var))) %>%
    mutate(zscore = (avg_genesig - mean(avg_genesig)) / sd(avg_genesig),
        pred_response = ifelse(zscore > thresh, nonresp_label, response_label) %>% as.factor())   
    
    # print patient level performance stats 
    xtab <- table(patient_summary$pred_response, patient_summary$response)
    print(confusionMatrix(xtab, positive = response_label))

    # evaluate genesig performance with ROC
    patient_roc <- suppressMessages(roc(patient_summary$response, patient_summary$zscore))


    # visualizing ROC curve and reporting AUC:
    dir.create(here(out_dir), showWarnings = FALSE)    
    if (img_format == 'png') {
        png(here(out_dir, paste0(genesig_name, '_ROC.png')), height=600, width=600)
    } else {
        pdf(here(out_dir, paste0(genesig_name, '_ROC.pdf')), height=600, width=600)
    }

    plot(patient_roc, cex.lab = 2, cex.axis = 2)

    if (cell_level) {
        # get cell level response table
        cell_summary <- metadata %>%
            select(response=all_of(response_var), 
                genesig_score=all_of(genesig_var)) 
        
        # evaluate genesig performance with ROC
        cell_roc <- suppressMessages(roc(cell_summary$response, cell_summary$genesig_score))

        # add cell level ROC trace
        lines(cell_roc, col="blue")
        legend("bottomright", 
                cex = 2,
                legend=c(paste0("Patient level, AUC=", patient_roc$auc %>% round(4)), 
                        paste0("Cell level,  AUC=", cell_roc$auc %>% round(4))), 
                        col=c("black", "blue"), lty = c(1,1))
    } else {
        legend("bottomright", 
        cex = 2,
        legend=paste0("AUC=", 
        patient_roc$auc %>% round(4)), 
        col="black", ty = 1) 
    }

    out <- dev.off() 


    if (save_data) {
        dir.create(here(out_dir), showWarnings = FALSE)  
        write.csv(patient_summary, here(out_dir, 'patient_response_summary.csv'), row.names=F)
    }

}


predict_response <- function(metadata, 
                             genesig_name,
                             id_var, 
                             out_dir,
                             save_data, 
                             thresh) {

    genesig_var <- paste0(genesig_name, "_sig_score")

    patient_summary <- metadata %>% group_by(patient_id=get(id_var)) %>%
    summarise(avg_genesig = mean(get(genesig_var))) %>%
    mutate(zscore = (avg_genesig - mean(avg_genesig)) / sd(avg_genesig),
            pred_response = ifelse(zscore > thresh, 0, 1))  
    
    if (save_data) {
        dir.create(here(out_dir), showWarnings = FALSE)  
        write.csv(patient_summary, here(out_dir, 'patient_response_summary.csv'), row.names=F)
    }
}
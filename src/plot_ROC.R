plot_ROC <- function(metadata, 
                     genesig_name,
                     resp_label, 
                     response_colname,
                     img_format,
                     out_dir,
                     cell_level, 
                     id_colname=NULL
                     ) {

    genesig_colname = paste0(genesig_name, "_sig_score")

    # get patient/cell level response tables
    patient_response_table <- metadata %>% 
        select(all_of(id_colname), 
                response_group=all_of(response_colname)) %>% unique %>%
        inner_join(metadata %>% 
                    group_by_(id_colname) %>% 
                    summarise(genesig_score = mean(get(genesig_colname)))) %>% suppressMessages()

    if (cell_level) {
        cell_response_table <- metadata %>%
            select(response_group=all_of(response_colname), 
                   genesig_score=all_of(genesig_colname))
    }


    # evaluate genesig performance with ROC
    patient_roc_curve <- roc(response = ifelse(patient_response_table$response_group == resp_label, 1, 0), predictor = patient_response_table$genesig_score) %>% suppressMessages()

    if (cell_level) {
        cell_roc_curve <- roc(response = ifelse(cell_response_table$response_group == resp_label, 1, 0), predictor = cell_response_table$genesig_score) %>% suppressMessages()
    }

  
    # visualizing ROC curve and reporting AUC:
    dir.create(here(out_dir), showWarnings = FALSE)    
    if (img_format == 'png') {
        png(here(out_dir, paste0(genesig_name, '_ROC.png')))
    } else {
        pdf(here(out_dir, paste0(genesig_name, '_ROC.pdf')))
    }


    plot(patient_roc_curve, main = "Reponse group prediction ROC")
    if (cell_level) {
        lines(cell_roc_curve, col="blue")
        legend("bottomright", 
                legend=c(paste0("Patient level, AUC=", patient_roc_curve$auc %>% round(4)), 
                        paste0("Cell level,  AUC=", cell_roc_curve$auc %>% round(4))), 
                        col=c("blue", "blue"), lty = c(1,1))
    } else {
        legend("bottomright", 
                legend=paste0("AUC=", 
                patient_roc_curve$auc %>% round(4)), 
                col="black", ty = 1) 
    }

    out <- dev.off() 


}
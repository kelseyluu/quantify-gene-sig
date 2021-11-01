plot_ROC <- function(metadata, 
                     genesig_name,
                     resp_label, 
                     response_colname,
                     img_format,
                     patient_level=T, 
                     id_colname=NULL
                     ) {

    genesig_colname = paste0(genesig_name, "_sig_score")

    # get patient/cell level response tables
    if (patient_level) {
        response_table <- metadata %>% 
            select(all_of(id_colname), 
                   response_group=all_of(response_colname)) %>% unique %>%
            inner_join(metadata %>% 
                       group_by_(id_colname) %>% 
                       summarise(genesig_score = mean(get(genesig_colname)))) %>% suppressMessages()
    } else {
        response_table <- metadata %>%
            select(response_group=all_of(response_colname), 
                   genesig_score=all_of(genesig_colname))
    }


    # evaluate genesig performance with ROC
    roc_curve <- roc(response = ifelse(response_table$response_group == resp_label, 1, 0), predictor = response_table$genesig_score) %>% suppressMessages()

  
    # visualizing ROC curve and reporting AUC:
    if (img_format == 'png') {
        png(paste0('./output/ROC_', ifelse(patient_level, 'patient', 'cell'), '_level.png'))
    } else {
        pdf(paste0('./output/ROC_', ifelse(patient_level, 'patient', 'cell'), '_level.pdf'))
    }

    plot(roc_curve,
        main = paste(ifelse(patient_level, 'Patient', 'Cell'), "level ROC")) 
    legend("bottomright", 
            legend=paste0("AUC=", roc_curve$auc %>% round(4)),
                        col="black", 
                        lty = 1) 
    
    out <- dev.off() 
}
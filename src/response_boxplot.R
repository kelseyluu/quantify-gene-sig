response_boxplot <- function(metadata,
                            genesig_name,
                            group_colname,
                            id_colname,
                            comparisons,
                            cell_level,
                            img_format,
                            out_dir
                            ) {
  
genesig_colname = paste0(genesig_name, "_sig_score")

  # create df summarizing average genesig score per patient
  genesig_df <- metadata %>%
    group_by_(id_colname) %>%
    summarise(total_cells = n(),
              avg_geneSig_score = mean(get(genesig_colname))) 
  
  # merge genesig_df with response group 
  response_df <- suppressMessages(metadata %>%
    select(all_of(id_colname), all_of(group_colname)) %>% unique %>%
    inner_join(genesig_df))

  
  response_df %>% 
    ggplot(aes(x=get(group_colname), y=avg_geneSig_score, fill=get(group_colname))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(color="black", size=1.5) +
      stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                         comparisons = comparisons, 
                         method = "wilcox.test", 
                         method.args = list(alternative = 'two.sided'), 
                         size = 4.5) + 
      theme_classic(base_size = 18) +
      labs(title=paste(genesig_name, 'patient level boxplots'),
          x=NULL,
          y='average gene signature score') +
      NoLegend()


dir.create(here(out_dir), showWarnings = FALSE)    
ggsave(here(out_dir, paste0(genesig_name, '_patient_boxplot.', img_format))) %>% suppressMessages()



if (cell_level) {
  genesig_colname = paste0(genesig_name, "_sig_score")

  metadata %>% 
    ggplot(aes(x=get(group_colname), y=get(genesig_colname))) +
      geom_violin(aes(fill=get(group_colname)), 
                  scale='width', 
                  draw_quantiles=c(0.25, 0.5, 0.75)) +
      geom_boxplot(outlier.shape = NA, width=0.04) +
      stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                         comparisons = comparisons, 
                         method = "wilcox.test", 
                         method.args = list(alternative = 'two.sided'), 
                         size = 4.5) +
      theme_classic(base_size = 18) +
      NoLegend() + 
      labs(title=paste(genesig_name, 'cell level violin plot'), 
          x=NULL,
          y='gene signature score')


  dir.create(here(out_dir), showWarnings = FALSE)    
  ggsave(here(out_dir, paste0(genesig_name, '_cell_violinplot.', img_format))) %>% suppressMessages()
}




}

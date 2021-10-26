plot_patient_level_boxplot <- function(metadata,
                                       genesig_colname,
                                       group_colname,
                                       id_colname,
                                       comparisons,
                                       wilcox_direction = 'less') {
  
  # create df summarizing average genesig score per patient
  genesig_df <- metadata %>%
    group_by_(id_colname) %>%
    summarise(total_cells = n(),
              avg_geneSig_score = mean(get(genesig_colname)))
  
  # merge genesig_df with response group 
  response_df <- metadata %>%
    select(id_colname, group_colname) %>% unique %>%
    inner_join(genesig_df)

  
  # save plot 
  png(here('plots', 'patient_boxplot.png'))

  response_df %>% ggplot(aes(x=get(group_colname),
  y=avg_geneSig_score, fill=get(group_colname))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color="black", size=1.5) +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                       comparisons = comparisons, 
                       method = "wilcox.test", 
                       method.args = list(alternative = wilcox_direction), 
                       size = 4.5) + 
    theme_classic(base_size = 18) +
    labs(title='patient level boxplots',
         x=NULL,
         y='average gene signature score') +
    NoLegend()

  dev.off() 
}

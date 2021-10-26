cell_violinplot <- function(metadata,
                            genesig_name,
                            group_colname
                            #  comparisons,
                            #  wilcox_direction="less"
                             ) {
  
genesig_colname = paste0(genesig_name, "_sig_score")

metadata %>% 
  ggplot(aes(x=get(group_colname), y=get(genesig_colname))) +
    geom_violin(aes(fill=get(group_colname)), 
                scale='width', 
                draw_quantiles=c(0.25, 0.5, 0.75)) +
    geom_boxplot(outlier.shape = NA, width=0.04) +
    # stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
    #                    comparisons = comparisons, 
    #                    method = "wilcox.test", 
    #                    method.args = list(alternative = wilcox_direction), 
    #                    size = 4.5) +
    theme_classic(base_size = 18) +
    NoLegend() + 
    labs(title=paste(genesig_name, 'cell level violin plot'), 
        x=NULL,
        y='gene signature score')

ggsave('./output/cell_violinplot.png')

}
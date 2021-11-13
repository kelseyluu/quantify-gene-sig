#! /usr/bin/env Rscript

# define docopt arguments
options(warn=-1)
suppressMessages(library(docopt))

'Quantify gene signature in scRNA-seq data.
Usage:
    quantify_gene_sig.R [options] <seurat> <genesig> <id> <out_dir>

Arguments:
    seurat  Path to seurat object *.RMD file containing scRNA-seq data.
    genesig  Path to gene signature *.txt file, one gene name per line.
    id   Name of column in seurat metadata indicating unique patient id.
    out_dir   Directory to save output files. Will be created if directory does not exist.

Options:
    -r, --response=<response_var>   Name of column in seurat metadata indicating treatment response groups.
    -l, --label=<response_label>    Label denoting positive response patients (if response known and 2 response groups present). [Default: 1] 
    -c, --cluster=<cluster_var>  Name of column in seurat metadata indicating cell type clusters. Used to label UMAP clusters.
    -n, --name=<genesig_name>   Gene signature name to use for output file prefixes. [Default: gene_signature]
    -t, --thresh=<threshold>    Z-score threshold for patient response classification. [Default: -0.2468496]
    --rankings=<cell_rankings>  Precomputed AUCell cell rankings. 
    --cell  Output cell level plots in addition to patient level. 
    --wilcox_dir=<wilcox_dir>  Direction of Wilcoxon test comparisons between response groups. One of: "two.sided", "less", "greater". [Default: two.sided]
    -f, --format=<img_format>   Format to save plots. Either png or pdf. [Default: png]
    --save_data   Export gene signature scores.
    -h, --help  Show help screen.
    -v, --version  Show version.
' -> doc


arguments <- docopt(doc, version = 'gene_signature_quantification v1.0\n\n')


# load required packages and scripts
suppressMessages(library(here))
suppressMessages(library(tictoc))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(AUCell))
suppressMessages(library(pROC))
suppressMessages(library(caret))
suppressMessages(library(ggpubr))


source('./src/run_AUCell.R')
source('./src/response_boxplot.R')
source('./src/genesig_UMAP.R')
source('./src/predict_response.R')

div <- function() {
    cat('\n-----------------------\n\n')
}


# check arguments are valid
if (!arguments$format %in% c('png', 'pdf')) {
    stop('Invalid image format.')
}


cat('\nLoading data...\n')
tic('Done')
    # load seurat object and cell rankings
    seurat_obj <- readRDS(arguments$seurat)
    if (!is.null(arguments$rankings)) {
        cells_rankings <- readRDS(arguments$rankings)
    } else {
        cells_rankings <- NULL
    }
    n_cells <- length(Cells(seurat_obj))

    if (!is.null(arguments$response)) {
        response_groups <- unique(seurat_obj@meta.data[[arguments$response]])
        n_groups <- length(response_groups)
        combs_df <- t(combn(response_groups %>% as.character(), 2))
        comparisons <- split(combs_df, seq(nrow(combs_df)))

        cat(paste('Read Seurat object with', n_cells, 'cells and', n_groups, 'response groups:',  paste(response_groups, collapse = ', '), '\n'))
    } else {
        cat(paste('Read Seurat object with', n_cells, 'cells. No response groups specified. \n'))
    }



    # load in genesig
    genesig <- list(scan(arguments$genesig, 'character', quiet=TRUE)) 
    cat(paste('Read', length(genesig[[1]]), 'genes in signature.\n'))
    names(genesig) <- arguments$name
toc()
div()


cat('Quantifying signature...\n')
tic('Done')
    seurat_obj <- run_AUCell(seurat_obj, genesig, cells_rankings, arguments$save_data, arguments$out_dir) 
toc()
div()


if (!is.null(arguments$response)) {
    cat('Evaluating classifier performance...\n')
    tic('Done')
        response_boxplot(seurat_obj@meta.data,
                        arguments$name,
                        arguments$response,
                        arguments$id,
                        comparisons,
                        arguments$wilcox_dir,
                        arguments$cell,
                        arguments$format,
                        arguments$out_dir)


        if (!is.null(seurat_obj@reductions$umap)) {
            genesig_UMAP(seurat_obj, 
                        arguments$name,
                        arguments$response,
                        arguments$cluster,
                        arguments$format,
                        arguments$out_dir)
        }

        if (n_groups == 2) {
            plot_ROC(seurat_obj@meta.data,
                    arguments$name,
                    arguments$id,
                    arguments$response,
                    arguments$label,
                    arguments$cell,
                    arguments$format,
                    arguments$out_dir,
                    arguments$save_data,
                    as.numeric(arguments$thresh)
                    )
        }

    toc()

} else {
    cat('Generating response predictions...\n')
    tic('Done')
    predict_response(seurat_obj@meta.data,
            arguments$name,
            arguments$id,
            arguments$out_dir,
            arguments$save_data,
            arguments$thresh
            )
    toc()
}

cat('\n')


#! /usr/bin/env Rscript

# define docopt arguments
options(warn=-1)
suppressMessages(library(docopt))

'Quantify gene signature in scRNA-seq data.
Usage:
    quantify_gene_sig.R [options] <seurat> <genesig> <id_var> <group_var> <out_dir>

Arguments:
    seurat  path to seurat object *.RMD file containing scRNA-seq data
    genesig  path to gene signature *.txt file, one gene name per line
    id_var  name of column in seurat metadata indicating unique patient id
    group_var   name of column in seurat metadata indicating treatment response groups
    out_dir   directory to save output files 

Options:
    -h, --help  Show help screen.
    -v, --version  Show version.
    -n, --name=<genesig_name>   Gene signature name. [Default: gene_signature]
    -r, --rankings=<cell_rankings>  Precomputed AUCell cell rankings. 
    -l, --label=<response_label>    Label denoting positive response patients (if 2 response groups). [Default: 1]
    -c, --cell_level    Output cell level plots in addition to patient level. 
    -f, --format=<img_format>   Format to save plots. Either png or pdf. [Default: png]
' -> doc


INTERACTIVE = FALSE
if (INTERACTIVE) {
    arguments <- docopt(doc, version = 'gene_signature_quantification v1.0\n\n',
                        args = c('./data/subsampled_cd3min.RDS', './data/genesig.txt', 'abbr_group'))

    arguments
} else {
    arguments <- docopt(doc, version = 'gene_signature_quantification v1.0\n\n')
}


# load required packages and scripts
suppressMessages(library(here))
suppressMessages(library(tictoc))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(AUCell))
suppressMessages(library(pROC))
suppressMessages(library(ggpubr))

source('./src/run_AUCell.R')
source('./src/response_boxplot.R')
source('./src/genesig_UMAP.R')
source('./src/plot_ROC.R')

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
    response_groups <- unique(seurat_obj@meta.data[[arguments$group_var]])
    n_groups <- length(response_groups)

    cat(paste('Read Seurat object with', n_cells, 'cells and', n_groups, 'response groups:',  paste(response_groups, collapse = ', '), '\n'))


    # load in genesig
    genesig <- list(scan(arguments$genesig, 'character', quiet=TRUE)) 
    cat(paste('Read', length(genesig[[1]]), 'genes in signature.\n'))
    names(genesig) <- arguments$name
toc()
div()


cat('Quantifying signature...\n')
tic('Done')
    seurat_obj <- run_AUCell(seurat_obj, genesig, cells_rankings) 
toc()
div()


cat('Generating plots...\n')
tic('Done')
    response_boxplot(seurat_obj@meta.data,
                    arguments$name,
                    arguments$group_var,
                    arguments$id_var,
                    arguments$cell_level,
                    arguments$format,
                    arguments$out_dir)


    if (!is.null(seurat_obj@reductions$umap)) {
        genesig_UMAP(seurat_obj, 
                    arguments$name,
                    arguments$group_var,
                    arguments$format,
                    arguments$out_dir)
    }


    if (n_groups == 2) {
        plot_ROC(seurat_obj@meta.data,
                arguments$name,
                arguments$label,
                arguments$group_var,
                arguments$format,
                arguments$out_dir,
                arguments$cell_level,
                arguments$id_var)
    }

toc()

cat('\n')

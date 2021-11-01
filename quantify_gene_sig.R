#! /usr/bin/env Rscript

# define docopt arguments
options(warn=-1)
suppressMessages(library(docopt))

'Quantify gene signature in scRNA-seq data.
Usage:
    quantify_gene_sig.R [options] <seurat> <genesig> <id_var> <group_var> 
    
Options:
    -h, --help  Show help screen.
    -v, --version  Show version.
    -n, --name=<genesig_name>   Gene signature name. [Default: gene_signature]
    -r, --rankings=<cell_rankings>  Precomputed AUCell cell rankings. 
    -l, --label=<response_label>    Label denoting positive response patients. [Default: 1]
    -c, --cell_level    Output cell level ROC plot in addition to patient level. 
    -f, --format=<img_format>   Format to save plots. Either png or pdf. [Default: png]
    
Arguments:
    seurat  path to seurat object containing scRNA-seq data
    genesig  path to gene signature, one gene name per line
    group_var   name of column in seurat metadata indicating treatment response groups
    id_var  name of column in seurat metadata indicating unique patient id
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
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(AUCell))
suppressMessages(library(pROC))
suppressMessages(library(tictoc))

source('./src/run_AUCell.R')
source('./src/cell_violinplot.R')
source('./src/patient_boxplot.R')
source('./src/genesig_UMAP.R')
source('./src/plot_ROC.R')


# check arguments are valid
if (!arguments$format %in% c('png', 'pdf')) {
    stop('Invalid image format.')
}


tic('Loading data')
# load seurat object and cell rankings
seurat_obj <- readRDS(arguments$seurat)
if (!is.null(arguments$rankings)) {
    cells_rankings <- readRDS(arguments$rankings)
} else {
    cells_rankings <- NULL
}


# load in genesig
genesig <- list(scan(arguments$genesig, 'character', quiet=TRUE)) 
cat(paste('Read', length(genesig[[1]]), 'genes\n'))
names(genesig) <- arguments$name
toc()


tic('Quantifying signature')
seurat_obj <- run_AUCell(seurat_obj, genesig, cells_rankings) 
toc()


tic('Generating plots')
cell_violinplot(seurat_obj@meta.data, 
                arguments$name, 
                arguments$group_var, 
                arguments$format)

patient_boxplot(seurat_obj@meta.data,
                arguments$name,
                arguments$group_var,
                arguments$id_var,
                arguments$format)

genesig_UMAP(seurat_obj, 
             arguments$name,
             arguments$group_var,
             arguments$format)

plot_ROC(seurat_obj@meta.data,
         arguments$name,
         arguments$label,
         arguments$group_var,
         arguments$format,
         patient_level=T,
         arguments$id_var)

if (arguments$cell_level) {
   plot_ROC(seurat_obj@meta.data,
         arguments$name,
         arguments$label,
         arguments$group_var,
         arguments$format,
         patient_level=F) 
}

toc()


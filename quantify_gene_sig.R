#! /usr/bin/env Rscript

# define docopt arguments
options(warn=-1)
suppressMessages(library(docopt))

'Quantify gene signature in scRNA-seq data.
Usage:
    quantify_gene_sig.R [--name=<genesig_name>] <seurat> <genesig> <group_var> <id_var>
    
Options:
    -h --help  Show help screen.
    -v --version  Show version.
    -n --name=<genesig_name>   Gene signature name. [Default: gene_signature]
    

Arguments:
    seurat  path to seurat object containing scRNA-seq data
    genesig  path to gene signature, one gene name per line
    group_var   name of column in seurat metadata indicating groups to plot
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
suppressMessages(library(tictoc))

source('./src/run_AUCell.R')
source('./src/cell_violinplot.R')
source('./src/patient_boxplot.R')


tic('Loading data')
# load seurat object
seurat_obj <- readRDS(arguments$seurat)

# load in genesig
genesig <- list(scan(arguments$genesig, 'character', quiet=TRUE)) 
cat(paste('Read', length(genesig[[1]]), 'genes\n'))
names(genesig) <- arguments$name
toc()


tic('Quantifying signature')
seurat_obj <- run_AUCell(seurat_obj, genesig)
toc()


tic('Generating plots')
cell_violinplot(seurat_obj@meta.data, 
                arguments$name, 
                arguments$group_var)
patient_boxplot(seurat_obj@meta.data,
                arguments$name,
                arguments$group_var,
                arguments$id_var)
toc()


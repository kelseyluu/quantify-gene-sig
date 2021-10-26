# define docopt arguments
suppressMessages(library(docopt))

'Quantify gene signature in scRNA-seq data.
Usage:
    test.R [--name=<name>] <seurat> <genesig> <group_colname>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    -n --name=<genesig_name>   Gene signature name. [Default: gene_signature]
    

Arguments:
    seurat  path to seurat object containing scRNA-seq data
    genesig  path to gene signature, one gene name per line
    group_colname   name of column in seurat object indicating groups to plot
' -> doc


INTERACTIVE = FALSE
if (INTERACTIVE) {
    arguments <- docopt(doc, version = 'gene_signature_quantification v1.0\n\n',
                        args = c('./pbmc.RDS', '../results/rank_ordered_CRvsPD_C4D1_geneSet.txt'))

    arguments
} else {
    arguments <- docopt(doc, version = 'gene_signature_quantification v1.0\n\n')
}


# load required packages and scripts
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(AUCell))
suppressMessages(library(tictoc))

source('../capstone_scripts/run_AUCell.R')
source('../capstone_scripts/cell_violinplot.R')
source('../capstone_scripts/patient_boxplot.R')


tic('loading data')
# load seurat object
seurat_obj <- readRDS(arguments$seurat)

# load in genesig
genesig <- list(scan(arguments$genesig, 'character'))
names(genesig) <- arguments$name
toc()


tic('quantifying signature')
seurat_obj <- run_AUCell(seurat_obj, genesig)
toc()


genesig_colname = paste0(genesig_name, "_sig_score")
cell_violinplot(seurat_obj@metadata, genesig_colname, group_colname)
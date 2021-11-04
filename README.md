# QuantSig
Command line tool for quantifying and evaluating gene signatures in single-cell RNA-seq datasets for cancer immunotherapy response prediction. 

## Usage

``` txt
Usage: quantify_gene_sig.R [options] <seurat> <genesig> <id_var> <group_var> 

Arguments:
    seurat:     Path to seurat object RDS file containing scRNA-seq data.
    genesig:    Path to gene signature, one gene name per line.
    group_var:  Name of column in seurat metadata indicating treatment response groups.
    id_var:     Name of column in seurat metadata indicating unique patient id.
    
Options:
    -h, --help                        Show help screen.
    -v, --version                     Show version.
    -n, --name <genesig_name>         Gene signature name. [Default: gene_signature]
    -r, --rankings <cell_rankings>    Precomputed AUCell cell rankings. 
    -l, --label <response_label>      Label denoting positive response patients. [Default: 1]
    -c, --cell_level                  Output cell level ROC plot in addition to patient level. 
    -f, --format <img_format>         Format to save plots. Either png or pdf. [Default: png]
```

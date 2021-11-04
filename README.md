# QuantSig
Command line tool for quantifying and evaluating gene signatures in single-cell RNA-seq datasets for cancer immunotherapy response prediction. 

## Usage

``` txt
Usage: quantify_gene_sig.R [options] <seurat> <genesig> <id_var> <group_var> <out_dir>

Arguments:
    seurat:     Path to seurat object *.RDS file containing scRNA-seq data.
    genesig:    Path to gene signature *.txt file, one gene name per line.
    id_var:     Name of column in seurat metadata indicating unique patient id.
    group_var:  Name of column in seurat metadata indicating treatment response groups.
    out_dir:    Directory to save output files. Will be created if directory does not exist.
    
Options:
    -h, --help                        Show help screen.
    -v, --version                     Show version.
    -n, --name <genesig_name>         Gene signature name to use for output file prefixes. [Default: gene_signature]
    -r, --rankings <cell_rankings>    Precomputed AUCell cell rankings. 
    -l, --label <response_label>      Label denoting positive response patients (if 2 response groups). [Default: 1]
    -c, --cell_level                  Output cell level plots in addition to patient level.  
    -f, --format <img_format>         Format to save plots. Either png or pdf. [Default: png]
```

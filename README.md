# QuantSig
Command line tool for quantifying and evaluating gene signatures in single-cell RNA-seq datasets for cancer immunotherapy response prediction. 

## Usage

``` txt
Usage: quantify_gene_sig.R [options] <seurat> <genesig> <id> <out_dir>

Arguments:
    seurat:     Path to seurat object *.RDS file containing scRNA-seq data.
    genesig:    Path to gene signature *.txt file, one gene name per line.
    id:         Name of column in seurat metadata indicating unique patient id.
    out_dir:    Directory to save output files. Will be created if directory does not exist.
    
Options:
    -r, --response <response_var>   Name of column in seurat metadata indicating treatment response groups.
    -l, --label <response_label>    Label denoting positive response patients (if response known and 2 response groups present). [Default: 1] 
    -c, --cluster <cluster_var>     Name of column in seurat metadata indicating cell type clusters. Used to label UMAP clusters.
    -n, --name <genesig_name>       Gene signature name to use for output file prefixes. [Default: gene_signature]
    -t, --thresh <threshold>        Z-score threshold for patient response classification. [Default: -0.2468496]
    --rankings <cell_rankings>      Precomputed AUCell cell rankings. 
    --cell                          Output cell level plots in addition to patient level. 
    --wilcox_dir <wilcox_dir>       Direction of Wilcoxon test comparisons between response groups. One of: "two.sided", "less", "greater". [Default: two.sided]
    -f, --format <img_format>       Format to save plots. Either png or pdf. [Default: png]
    --save_data                     Export gene signature scores.
    -h, --help                      Show help screen.
    -v, --version                   Show version.
```

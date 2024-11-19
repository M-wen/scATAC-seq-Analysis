#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("-i", "--input", required=TRUE, help="pseudo multiomics rds file")
parser$add_argument("-o", "--outpath", required=TRUE, help="out put path")

args <- parser$parse_args()

library(tidyverse)
library(Seurat)
library(FigR)
library(optmatch)
library(data.table)
library(Pando)
library("BSgenome.Hsapiens.UCSC.hg38")
# muo_data <- read_rds('${project_path}/results/06_integrated/pbmc_multiomic_object.rds')
muo_data <- read_rds(args$input)
muo_data <- initiate_grn(
    muo_data,
    rna_assay = 'RNA',
    peak_assay = 'peaks',
    regions = phastConsElements20Mammals.UCSC.hg38,
    exclude_exons = TRUE
)
data('motifs')
data('motif2tf')
muo_data <- find_motifs(
    muo_data, 
    pfm = motifs, 
    motif_tfs = motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
)
genesUsed <- unique(rownames(muo_data@assays$RNA@data))
regions <- NetworkRegions(muo_data)
muo_data <- infer_grn(
    muo_data,
    peak_to_gene_method = 'Signac',
    genes = genesUsed,
    parallel = F
)
muo_data <- find_modules(
    muo_data, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# saveRDS(muo_data,"${project_path}/results/06_integrated/pbmc_GRN.rds")
saveRDS(muo_data,paste0(args$outpath,"/results/06_integrated/pbmc_GRN.rds"))
muo_data <- get_network_graph(muo_data, graph_name='umap_graph')
plot_network_graph(muo_data, graph='umap_graph')

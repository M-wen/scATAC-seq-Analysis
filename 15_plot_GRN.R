#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("-i", "--input", required=TRUE, help="input GRN rds file")
parser$add_argument("-o", "--outpath", required=TRUE, help="out put path")

args <- parser$parse_args()


library(Seurat)
library(Signac)
library(tidyverse)
library(AUCell)
library(SCENIC)
library(Pando)
# object <- readRDS("${project_path}/06_integrated/pbmc_GRN.rds")
setwd(args$outpath)
object <- readRDS(args$input)
grn_module <- object@grn@networks$glm_network@modules@meta
grn_module_2 <- grn_module[which(grn_module$n_genes > 1),]
tf <- unique(grn_module_2$tf)
tf_target_meta <- as.data.frame(grn_module_2)
tf_target <- lapply(unique(tf),function(x){
    targe_genes <- tf_target_meta[which(tf_target_meta$tf == x),"target"]
})

names(tf_target) <- tf
cellInfo <- data.frame(object@meta.data)
assaydata <- object@assays$RNA@data
assaydata <- assaydata[,which(colnames(assaydata) %in% rownames(cellInfo))]
cells_rankings <- AUCell_buildRankings(assaydata, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(tf_target, cells_rankings)
raw_rss <- calcRSS(AUC = cells_AUC, cellAnnotation = cellInfo[,"seurat_annotations"])
library(viridis)
bk = seq(0.1, 1,by = 0.01)
col_length = length(bk)
if_1 = inferno(col_length/2)
if _2 = rev(mako(col_length/2))
col= append(if_2,if_1)
library(pheatmap)
pheatmap(raw_rss, cluster_cols = T, cluster_rows = T,
         show_rownames = T,
         color = col,
         filename = "./results/06_integrated/celltype_target_gene_rss_heatmap.pdf",width=6,height=12
        )
p <- DotPlot(object,features=names(tf_target),group.by="seurat_annotations") +
coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave(paste0(args$outpath,"/results/06_integrated/celltype_tf_dotplot.pdf"),
      p,width = 6,height=12)


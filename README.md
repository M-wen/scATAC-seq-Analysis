# scATAC-seq-Analysis
Single-cell sequencing assay for transposase-accessible chromatin (scATAC-seq) is a powerful tool for chromatin accessibility research. The potential of this technology is attracting an increasing number of researchers to generate scATAC-seq data. Here, we develop a practical tutorial covering the key steps of a typical scATAC-seq analysis, from processing FASTQ files to downstream analysis, aimed at newcomers to the field. To illustrate the practical application of these steps, we have incorporated practice suggestions into a workflow and applied it to a publicly available dataset. First, we outline the steps for data pre-processing of scATAC-seq. Next, we detail the steps for downstream analysis. Finally, we describe the steps for multi-omics integration analysis.
## 01_cellranger
Fastq files were aligned to genome, and valid cell barcodes were called using the cellranger-atac count command with the default parameters. 
## 02_detection_doublets_Amulet
Remove doublets using Amulet.
## 03_clustering
QC, dimension reduction, batch correction, clustering, embedding, and cell type anntation.
## 04_peakcalling
For each cell cluster or cell type, we combined all fragments to generate a pseudobulk ATAC-seq dataset for individual biological replicates. Furthermore, we generated two pseudo-replicates which consisit of half of the fragments from pooled cell cluster or cell type. We called peak for each of the pseudobulk datasets and the pooled dataset of all replicates independently.
## 05_parsepeak
Two types of peaks were retained to generate a list of reproducible peaks:
1. The pooled peak set and overlapped by more than 50% of peak length with a peak in both individual replicates.
2. The pooled peak set and overlapped by 50% or more of peak length with a peak in both pseudo-replicates
## 06_interative_overlap_peak_merging
1Merge all cell type or cluster peak sets to a union peak set.
## 07_addPeakMatrix
Add the peak matrix to the ArchR project.
## 08_findDAR
Identifying differentially accessible region.
## 09_motifAnalysis
Add motif annotation.
## 10_Coaccessibility
calculate the co-accesssibility.
## 11_integrated_scRNA_scATAC
Integrated analysis with scRNA-seq data.
## 12_pairedCell
cells from scATAC-seq paired with scRNA-seq.
## 13_create_pseudo_multiomics_cells
create pseudo multiomics cells.
## 14_inferCNV
infer CNV with the pseudo multiomics using Pando.
## 15_plot_GRN
Plot TF and target gene expression of the eRegulon.





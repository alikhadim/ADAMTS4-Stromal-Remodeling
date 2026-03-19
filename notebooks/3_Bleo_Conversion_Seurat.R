library(reticulate)
use_python("/Users/aly/miniconda3/envs/spatial/bin/python", required = TRUE)
library(sceasy)
library(Seurat)

base_path <- "/Users/aly/Downloads/2024_Ali_scRNAseq/2024_Mahsa_scRNAseq/ScanpyObjects/241129_Bleo_Harmony_v2"
h5ad_file <- paste0(base_path, ".h5ad")
rds_file  <- paste0(base_path, ".rds")

sceasy::convertFormat(h5ad_file, from = "anndata", to = "seurat", outFile = rds_file)
seurat_object <- readRDS(rds_file)

adata <- import("anndata")$read_h5ad(h5ad_file)
seurat_object@misc$celltype_colors <- setNames(as.character(adata$uns$celltype_colors), levels(seurat_object$celltype))
DimPlot(seurat_object, group.by = "celltype", cols = seurat_object@misc$celltype_colors)
table(seurat_object$sample, seurat_object$celltype)
saveRDS(seurat_object, rds_file)
##

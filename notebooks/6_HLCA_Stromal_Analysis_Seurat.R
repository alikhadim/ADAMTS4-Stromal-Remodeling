################## Human HCLA Stromal Analysis ##################

suppressPackageStartupMessages({
  libs <- c("Seurat", "presto", "dplyr", "tidyr", "stringr", "tibble", "readxl", "ggplot2", 
            "viridis", "RColorBrewer", "scico", "ggridges", "gridExtra", "cowplot", 
            "patchwork", "ComplexHeatmap", "colorRamp2", "clusterProfiler", "DOSE", 
            "pathview", "enrichplot", "DESeq2", "limma", "fgsea", "openxlsx", 
            "gplots", "org.Hs.eg.db", "biomaRt", "writexl", "scales")
  lapply(libs, library, character.only = TRUE)
})
myorgdb <- org.Hs.eg.db

seurat_object <- readRDS("/Users/aly/Downloads/2024_Ali_scRNAseq/HumanLungCellAtlas/Objects/HCLA_full_scanpy_Stroma.rds")
out_dir <- "/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis/HLCA/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

# Filter Disease
exclude <- c("pulmonary sarcoidosis", "lung large cell carcinoma", "pleomorphic carcinoma", 
             "chronic rhinitis", "chronic obstructive pulmonary disease", "lung adenocarcinoma", 
             "squamous cell lung carcinoma", "pneumonia", "cystic fibrosis", 
             "lymphangioleiomyomatosis", "hypersensitivity pneumonitis", "COVID-19")

seurat_object <- subset(seurat_object, subset = !disease %in% exclude & !is.na(disease))
seurat_object$disease <- factor(seurat_object$disease, 
                                levels = c("normal", "non-specific interstitial pneumonia", "pulmonary fibrosis", "interstitial lung disease"))

# Rename and Reorder Clusters
cluster_map <- c("Adventitial fibroblasts" = "AF2.1", "Alveolar fibroblasts" = "AF1", 
                 "Mesothelium" = "Meso", "Myofibroblasts" = "MyoFBs", 
                 "Peribronchial fibroblasts" = "PeriFBs", "Pericytes" = "Pericytes", 
                 "SM activated stress response" = "StressAcSMCs", "Smooth muscle" = "SMCs", 
                 "Smooth muscle FAM83D+" = "FAM83D+SMCs", "Subpleural fibroblasts" = "SpFBs", "Unknown" = "AF2.2")

seurat_object$seurat_clusters <- factor(recode(seurat_object$ann_finest_level, !!!cluster_map), 
                                        levels = c("AF1", "AF2.1", "AF2.2", "Meso", "MyoFBs", "PeriFBs", "Pericytes", "StressAcSMCs", "SMCs", "FAM83D+SMCs", "SpFBs"))
Idents(seurat_object) <- seurat_object$seurat_clusters

# Plotting
colors <- c("AF2.1"="#4CAF50", "AF1"="#FF9800", "Meso"="#03A9F4", "MyoFBs"="gold", "PeriFBs"="#8BC34A", 
            "Pericytes"="#673AB7", "StressAcSMCs"="#F44336", "SMCs"="#795548", "FAM83D+SMCs"="#E91E63", 
            "SpFBs"="#FF5722", "AF2.2"="#1E9E9E")

p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, alpha = 0.3, raster = FALSE) + 
  scale_color_manual(values = colors) + theme_minimal() + 
  ggtitle("HLCA Stromal Cells", subtitle = paste("n =", ncol(seurat_object)))

p2 <- DimPlot(seurat_object, reduction = "umap", group.by = "disease", pt.size = 1, alpha = 0.1, raster = FALSE) + 
  theme_minimal() + ggtitle("By Disease")

combined <- p1 + p2
ggsave("1.UMAP_Seurat_combined_v3.png", combined, width = 15, height = 6, dpi = 600)

print(table(seurat_object$seurat_clusters, seurat_object$disease))


# Table Generation
cell_counts <- as.data.frame(table(Cluster = Idents(seurat_object), Disease = seurat_object$disease))
cell_counts_wide <- cell_counts %>%
  pivot_wider(names_from = Disease, values_from = Freq, values_fill = 0)
cell_counts_wide <- bind_rows(cell_counts_wide, summarise(cell_counts_wide, across(where(is.numeric), sum), Cluster = "Total"))
write_xlsx(cell_counts_wide, "Cell_Cluster_Disease_Counts_v2.xlsx")

# Stacked Barplot
cluster_props <- cell_counts %>%
  group_by(Disease) %>%
  mutate(Proportion = Freq / sum(Freq) * 100) %>%
  ungroup()

plot <- ggplot(cluster_props, aes(x = Disease, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(title = "Cell Proportions by Disease", x = "Disease", y = "Proportion (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Cluster_Proportions_by_Disease_Filtered_v2.pdf", plot, width = 4, height = 10)


############################
##### now plot and statistics for each disease group #####
library(Seurat)
library(dplyr)
library(readxl)
library(writexl)
library(pheatmap)
library(RColorBrewer)

# 1. Differential Expression Analysis
diseases <- setdiff(levels(seurat_object$disease), "normal")
gene_map <- seurat_object[["RNA"]]@meta.features$feature_name

for (d in diseases) {
  Idents(seurat_object) <- seurat_object$disease
  markers <- FindMarkers(seurat_object, ident.1 = d, ident.2 = "normal", logfc.threshold = 0, min.pct = 0) %>%
    mutate(feature_name = gene_map[match(rownames(.), rownames(seurat_object))], Ensembl_ID = rownames(.)) %>%
    filter(is.finite(avg_log2FC))
  
  write_xlsx(markers, paste0(d, "_vs_normal_DEGs_full.xlsx"))
  write_xlsx(markers %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC)), paste0(d, "_vs_normal_DEGs_filtered.xlsx"))
}

# 2. Signature Filtering and Matrix Preparation
sig_path <- "/Users/aly/Downloads/2024_Ali_scRNAseq/GSEA_signatures/fromGeorge_marker_genes.xlsx"
af1_sig <- toupper(read_excel(sig_path)[[1]])

deg_files <- list.files(pattern = "_vs_normal_DEGs_full.*\\.xlsx$")
res_list <- lapply(setNames(deg_files, gsub("_vs_normal_DEGs_full.*", "", deg_files)), function(f) {
  read_excel(f) %>% filter(toupper(feature_name) %in% af1_sig, avg_log2FC < -0.5, p_val_adj < 0.05)
})

common <- Reduce(intersect, lapply(res_list, function(df) toupper(df$feature_name))) %>% .[!grepl("^MT[-−]ND4$", .)]

h_mat <- sapply(res_list, function(df) df %>% filter(toupper(feature_name) %in% common) %>% arrange(match(toupper(feature_name), common)) %>% pull(avg_log2FC))
p_mat <- sapply(res_list, function(df) df %>% filter(toupper(feature_name) %in% common) %>% arrange(match(toupper(feature_name), common)) %>% pull(p_val_adj))
rownames(h_mat) <- rownames(p_mat) <- common

p_labels <- apply(p_mat, c(1,2), function(x) if(x < 0.001) "***" else if(x < 0.01) "**" else if(x < 0.05) "*" else "")

# 3. Optimized Heatmap
pheatmap(
  h_mat, 
  display_numbers = p_labels,
  color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize_row = 7, fontsize_number = 6,
  cellwidth = 50, cellheight = 10,
  main = "AF1 Downregulated Genes Across Lung Diseases",
  filename = "HCLA_AF1_Heatmap.pdf",
  width = 8, height = 14
)

################# Normal vs IPF etc enrichplots ######
print(table(seurat_object$disease))
### normal non-specific interstitial pneumonia pulmonary fibrosis interstitial lung disease

sig_path <- "/Users/aly/Downloads/2024_Ali_scRNAseq/GSEA_signatures/fromGeorge_marker_genes.xlsx"
af1_sig <- toupper(read_excel(sig_path)[[1]])
gene_map <- seurat_object[["RNA"]]@meta.features$feature_name
gene_info <- data.frame(Ensembl = rownames(seurat_object), Symbol = toupper(gene_map))
valid_ids <- gene_info %>% filter(Symbol %in% af1_sig) %>% pull(Ensembl)

id1 <- "pulmonary fibrosis"
id2 <- "normal"

Idents(seurat_object) <- "disease"
de_res <- FindMarkers(seurat_object, ident.1 = id1, ident.2 = id2, logfc.threshold = 0, min.pct = 0)
ranks <- de_res$avg_log2FC
names(ranks) <- rownames(de_res)
ranks <- ranks[is.finite(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

gsea_res <- fgsea(pathways = list(AF1 = valid_ids), stats = ranks, minSize = 10)
plotEnrichment(valid_ids, ranks) + 
  labs(title = paste("AF1 Signature:", id1, "vs", id2),
       subtitle = paste0("NES: ", round(gsea_res$NES, 2), " | p-adj: ", formatC(gsea_res$padj, format = "e", digits = 2)))

ggsave(paste0("GSEA_", gsub(" ", "_", id1), "_vs_", gsub(" ", "_", id2), ".pdf"), width = 6, height = 4)

####################################################################################################################################
#######################################################
library(Seurat)
library(dplyr)
library(pheatmap)

seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- "disease"

avg_expr <- AverageExpression(seurat_object, assays = "RNA", group.by = "disease")$RNA
gene_meta <- seurat_object[["RNA"]]@meta.features

disease_cols <- c("non-specific interstitial pneumonia", "pulmonary fibrosis", "interstitial lung disease")
normal_col <- "normal"

fc_normal <- log2((avg_expr[, normal_col] + 1) / (rowMeans(avg_expr[, disease_cols]) + 1))
fc_disease <- log2((rowMeans(avg_expr[, disease_cols]) + 1) / (avg_expr[, normal_col] + 1))

top_normal <- names(sort(fc_normal, decreasing = TRUE)[1:10])
top_disease <- names(sort(fc_disease, decreasing = TRUE)[1:10])
selected_ids <- c(top_normal, top_disease)

plot_matrix <- t(scale(t(avg_expr[selected_ids, ])))
colnames(plot_matrix) <- colnames(avg_expr)
rownames(plot_matrix) <- gene_meta[rownames(plot_matrix), "feature_name"]

pdf("Disease_Gene_Expression_Heatmap.pdf", width = 3.5, height = 7)
pheatmap(
  plot_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  main = "Top DEGs across diseases"
)
dev.off()

######################
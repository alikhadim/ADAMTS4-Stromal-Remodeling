################## Adamts4 manuscript code after conversion from h5ad ##################

suppressPackageStartupMessages({
  libs <- c("Seurat", "SeuratExtend", "presto", "tidyverse", "readxl", "openxlsx", "writexl", 
            "ggplot2", "viridis", "RColorBrewer", "scico", "ggridges", "patchwork", "cowplot", 
            "gridExtra", "ggpubr", "tidyplots", "ComplexHeatmap", "circlize", "tidyheatmaps", 
            "gplots", "pheatmap", "scales", "reshape2", "scales", "ggtext")
  lapply(libs, library, character.only = TRUE)
})

seurat_object <- readRDS("/Users/aly/Documents/2026_Adamts4_Manuscript/Objects/241129_Bleo_Harmony_v3.rds")
head(seurat_object@meta.data)
seurat_object

output_dir <- "/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis"
output_directory <- "/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  
}
setwd(output_directory)
color_vector <- c("AF1" = "#1f77b4", "AF2" = "#ff7f0e", "MyoFB" = "#279e68", "PeriFB" = "#d62728")
seurat_object$celltype <- seurat_object$leiden_0.6_annotation
seurat_object$seurat_clusters <- seurat_object$celltype
Idents(seurat_object) <- seurat_object$celltype

p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, alpha = 0.5, label = TRUE) + 
  scale_color_manual(values = color_vector) + theme_minimal() + ggtitle("UMAP by Cluster")

p2 <- DimPlot(seurat_object, reduction = "umap", group.by = "sample", pt.size = 1, alpha = 0.5) + 
  scale_color_brewer(palette = "Set2") + theme_minimal() + ggtitle("UMAP by Sample")

combined_plot <- p1 + p2
print(combined_plot)
table(seurat_object$sample, seurat_object$seurat_clusters)
ggsave("1.UMAP_Seurat_combined.png", combined_plot, width = 20, height = 8, dpi = 600)

s_names <- as.character(unique(seurat_object$sample))
u_list <- lapply(s_names, function(s) {
  DimPlot(seurat_object[, seurat_object$sample == s], reduction = "umap", group.by = "seurat_clusters", 
          label = T, pt.size = 2, alpha = 0.4, cols = color_vector) + 
    ggtitle(s) + theme_minimal() + NoLegend() + theme(panel.grid = element_blank())
})
ggsave("2.umap_2x2_nogrid.png", wrap_plots(u_list, ncol = 2, nrow = 2), width = 12, height = 10, dpi = 300)

df <- seurat_object@meta.data %>% count(sample, seurat_clusters) %>% group_by(sample) %>% mutate(p = n/sum(n)*100)
writexl::write_xlsx(df, "props.xlsx")
ggplot(df, aes(sample, p, fill = seurat_clusters)) + geom_col() + 
  scale_fill_manual(values = color_vector) + theme_minimal() + theme(panel.grid = element_blank())
ggsave("3.props_nogrid.pdf", width = 5, height = 6)

#########################

output_directory <- getwd()
outdir <- file.path(output_directory, "Seurat_Extended")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
if(length(VariableFeatures(seurat_object)) == 0) seurat_object <- FindVariableFeatures(seurat_object)
genes <- VariableFeatures(seurat_object)

Idents(seurat_object) <- "seurat_clusters"
names(color_vector) <- levels(seurat_object)

pdf(file.path(outdir, "4.DimPlot2.pdf"), width = 6, height = 4)
DimPlot(seurat_object, cols = color_vector) + theme_umap_arrows()
dev.off()

pdf(file.path(outdir, "5.ClusterDistrBar.pdf"), width = 8, height = 2)
ClusterDistrBar(seurat_object$sample, seurat_object$seurat_clusters, cols = color_vector)
dev.off()

genes.zscore <- CalcStats(seurat_object, features = genes, group.by = "seurat_clusters", order = "p", n = 30, assay = "RNA")
pdf(file.path(outdir, "6.Heatmap_genes_zscore.pdf"), width = 2.5, height = 15)
Heatmap(genes.zscore, lab_fill = "zscore")
dev.off()


for(s in unique(seurat_object$sample)){
  seu_sub <- subset(seurat_object, subset = sample == s)
  g <- genes[genes %in% rownames(seu_sub)]
  if(length(g)==0) next
  z <- CalcStats(seu_sub, features=g, group.by="seurat_clusters", order="p", n=10, assay="RNA")
  if(is.null(z) || nrow(z)==0) next
  pdf(file.path(outdir, paste0("7.Heatmap_genes_zscore_", s, ".pdf")), width=2.5, height=8)
  print(Heatmap(z, lab_fill="zscore"))
  dev.off()
  message("Saved Heatmap for sample: ", s)
}

toplot2 <- CalcStats(seurat_object, features = genes[1:500], method = "zscore", order = "p")
Heatmap(toplot2, lab_fill = "zscore", feature_text_subset = genes[1:40], expand_limits_x = c(-0.5, 11))

grouped_features <- list(
  AF1 = c("Wnt2", "Scube2", "Fgf10", "Inmt", "Gdf10"),
  AF2 = c("Pi16", "Adh7", "Scara5"),
  SMC = c("Actc1", "Myh11", "Notch3", "Acta2"),
  PeriFB = c("Fgf18", "P2ry14", "Aspn"),
  Misc = c("Gli1", "Mki67", "Adamts4"), 
  MyoFB = c("Cthrc1", "Nrep", "Piezo2", "Grem1", "Fst")
)
pdf(file.path(outdir, "8.DotPlot2_grouped_features.pdf"), width = 3, height = 6)
DotPlot2(seurat_object, features = grouped_features)
dev.off()


pdf(file.path(outdir, "9.DotPlot2_grouped_features_by_sample.pdf"), width=8, height=8)
print(DotPlot2(seurat_object, features=grouped_features, group.by="seurat_clusters", split.by="sample", show_grid=FALSE))
dev.off()

samples <- unique(seurat_object$sample)
for(s in samples){
  seu_sub <- subset(seurat_object, subset = sample == s)
  pdf(file.path(outdir, paste0("10.DotPlot2_grouped_features_", s, ".pdf")), width = 3, height = 6)
  p <- DotPlot2(seu_sub, features = grouped_features)
  print(p)
  dev.off()
  message("Saved DotPlot2 for sample: ", s)
}

genes_vln <- c("Fgf10", "Adamts4", "Gdf10")
cells <- WhichCells(seurat_object)
fname <- paste0("11.VlnPlot2_", paste(genes_vln, collapse = "_"), ".pdf")

pdf(file.path(outdir, fname), width = 8, height = 3)
VlnPlot2(
  seurat_object,
  features = genes_vln,
  group.by = "seurat_clusters",
  pt = FALSE,
  comparisons = list(c(1,2), c(1,3), c(3,4)),
  cells = cells,
  cols = color_vector,
  hide.ns = TRUE,
  stat.method = "wilcox.test"
)
dev.off()

samples <- unique(seurat_object$sample)
for (s in samples) {
  seurat_sub <- subset(seurat_object, subset = sample == s)
  fname <- paste0("12.VlnPlot2_", s, "_", paste(genes_vln, collapse = "_"), ".pdf")
  
  pdf(file.path(outdir, fname), width = 8, height = 3)
  print(
    VlnPlot2(
      seurat_sub,
      features = genes_vln,
      group.by = "seurat_clusters",
      pt = FALSE,
      comparisons = list(c(1,2), c(1,3), c(3,4)),
      cols = color_vector,
      hide.ns = TRUE,
      stat.method = "wilcox.test"
    )
  )
  dev.off()
}


genes <- c("Scube2", "Cthrc1", "Lrp1")  # Replace with any 3 genes
p <- FeaturePlot3(seurat_object, feature.1 = genes[1], feature.2 = genes[2], feature.3 = genes[3], pt.size = 2)
fname <- paste0("FeaturePlot3_", paste(genes, collapse = "_"), ".pdf")
ggsave(file.path(outdir, fname), p, width = 10, height = 8)
samples <- unique(seurat_object$sample)
for (s in samples) {
  seurat_sub <- subset(seurat_object, subset = sample == s)
  p <- FeaturePlot3(seurat_sub, 
                    feature.1 = genes[1], 
                    feature.2 = genes[2], 
                    feature.3 = genes[3], 
                    pt.size = 2)
  fname <- paste0("13.FeaturePlot3_", s, "_", paste(genes, collapse = "_"), ".pdf")
  ggsave(file.path(outdir, fname), p, width = 10, height = 8)
}



genes <- c("Scube2","Adamts4","Cthrc1","Pi16")

# Full object
p1 <- FeaturePlot3(seurat_object, feature.1=genes[1], feature.2=genes[2], feature.3=genes[3], pt.size=2)
p2 <- FeaturePlot3(seurat_object, feature.1=genes[1], feature.2=genes[2], feature.3=genes[4], pt.size=2)
fname <- paste0("14.FeaturePlot4_",paste(genes,collapse="_"),".pdf")
ggsave(file.path(outdir,fname), cowplot::plot_grid(p1,p2,ncol=2), width=20, height=8)

# Per sample
samples <- unique(seurat_object$sample)
for(s in samples){
  seurat_sub <- subset(seurat_object, subset=sample==s)
  p1 <- FeaturePlot3(seurat_sub, feature.1=genes[1], feature.2=genes[2], feature.3=genes[3], pt.size=2)
  p2 <- FeaturePlot3(seurat_sub, feature.1=genes[1], feature.2=genes[2], feature.3=genes[4], pt.size=2)
  fname <- paste0("15.FeaturePlot4_",s,"_",paste(genes,collapse="_"),".pdf")
  ggsave(file.path(outdir,fname), cowplot::plot_grid(p1,p2,ncol=2), width=20, height=8)
}

seurat_object <- GeneSetAnalysis(seurat_object, genesets = hall50$mouse)
matr <- seurat_object@misc$AUCell$genesets

idents <- c("AF1", "MyoFB")
pdf(file.path(outdir, paste0("16.Waterfall_", idents[1], "_vs_", idents[2], ".pdf")), width=8, height=10)
WaterfallPlot(matr, f = seurat_object$seurat_clusters, ident.1 = idents[1], ident.2 = idents[2])
dev.off()

genes <- VariableFeatures(seurat_object)[1:50]
idents <- c("MyoFB", "AF1")
pdf(file.path(outdir, paste0("17.Waterfall_sigs_", idents[1], "_vs_", idents[2], ".pdf")), width=18, height=3)
WaterfallPlot(
  seurat_object,
  group.by = "seurat_clusters",
  features = genes,
  ident.1 = idents[1],
  ident.2 = idents[2],
  length = "logFC"
)
dev.off()


########################


genes <- c("Adamts4","Scube2","Pi16","Cthrc1","Fgf18")
expr_mat <- GetAssayData(seurat_object, layer="data")[genes,]
samples <- seurat_object$sample
adamts4_cells <- expr_mat["Adamts4",] > 0

coexp_counts <- sapply(genes[-1], function(g) tapply(expr_mat[g,] > 0 & adamts4_cells, samples, sum))
coexp_frac <- sapply(genes[-1], function(g) tapply(expr_mat[g,] > 0 & adamts4_cells, samples, sum)/tapply(adamts4_cells, samples, sum) * 100)

df <- melt(coexp_counts, varnames=c("Sample","Gene"), value.name="Count")
df$Fraction <- melt(coexp_frac, varnames=c("Sample","Gene"), value.name="Fraction")$Fraction

p <- ggplot(df, aes(Gene, Sample, size=Count, fill=Fraction)) +
  geom_point(shape=21, color="gray20", stroke=0.7, alpha=0.9) +
  scale_size_continuous(range=c(5,15), name="Number of cells") +
  scale_fill_gradientn(colors=rev(brewer.pal(9,"RdBu")), name="Percentage of Adamts4+ cells") +
  theme_minimal(base_size=16) +
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1, face="bold", size=13, color="#333333"),
    axis.text.y=element_text(face="bold", size=13, color="#333333"),
    axis.title=element_blank(),
    panel.grid=element_blank(),
    legend.position="right",
    legend.title=element_text(face="bold", size=13),
    legend.text=element_text(size=12),
    panel.border=element_rect(color="black", linewidth=0.8, fill=NA), 
    plot.background=element_rect(fill="#f8f8f8", color=NA)
  ) +
  guides(fill=guide_colorbar(barwidth=1.5, barheight=6), size=guide_legend(order=1))

ggsave("18.Adamts4_coexpression_plot.pdf", p, width=6, height=8)


#######################################################################################
#subset_seurat <- seurat_object
# ---------------------------- #
# LIF / MYO scores
# ---------------------------- #
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
sig_data <- readxl::read_xlsx("/Users/aly/Downloads/2024_Ali_scRNAseq/2024_Mahsa_scRNAseq/genes.xlsx")

plot_sig <- function(gene_set, name, obj) {
  g <- intersect(gene_set, rownames(obj))
  if(length(g) == 0) return(NULL)
  obj[[name]] <- colMeans(GetAssayData(obj, layer = "data")[g, , drop = FALSE])
  
  # Global
  p_glob <- FeaturePlot(obj, features = name, pt.size = 2) + 
    scale_color_viridis_c(option = "inferno") + theme_minimal() + NoGrid()
  ggsave(paste0(name, "_UMAP.pdf"), p_glob, width = 10, height = 8)
  ggsave(paste0(name, "_UMAP.png"), p_glob, width = 10, height = 8, dpi = 300, bg = "white")
  
  # Sample-wise Combined
  s_names <- as.character(unique(obj$sample))
  p_list <- lapply(s_names, function(s) {
    FeaturePlot(obj[, obj$sample == s], features = name, pt.size = 1.5) + 
      scale_color_viridis_c(option = "inferno") + theme_minimal() + NoGrid() + ggtitle(s)
  })
  
  combined <- wrap_plots(p_list, ncol = length(p_list))
  ggsave(paste0(name, "_combined.pdf"), combined, width = 15, height = 3)
  ggsave(paste0(name, "_combined.png"), combined, width = 15, height = 3, dpi = 300, bg = "white")
}

# Run for both signatures
plot_sig(sig_data$LIF, "LIF_Score", seurat_object)
plot_sig(sig_data$MYO, "MYO_Score", seurat_object)
##############################################

#####################################
# ---------------------------- #
# DEGs excel file and later heatmaps
# ---------------------------- #

Idents(seurat_object) <- "seurat_clusters"
m <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

wb <- createWorkbook()
for (n in c(10, 20, 100, 300, 500, 1000)) {
  top <- m %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
  addWorksheet(wb, paste0("Top", n))
  writeData(wb, paste0("Top", n), top[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")])
}
saveWorkbook(wb, file.path(output_directory, "21.Top_Markers_DEGs.xlsx"), overwrite = TRUE)

#####################################
# ---------------------------- #
# Volcano plot now
# ---------------------------- #

print(table(seurat_object$seurat_clusters))
outdir <- file.path(getwd(), "6.Volcanos"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
clusters <- levels(Idents(seurat_object))

for (i in 1:(length(clusters)-1)) for (j in (i+1):length(clusters)) {
  c1 <- clusters[i]; c2 <- clusters[j]
  if (sum(Idents(seurat_object)==c1)<30 | sum(Idents(seurat_object)==c2)<30) next
  res <- tryCatch(FindMarkers(seurat_object, ident.1=c1, ident.2=c2, group.by="seurat_clusters", logfc.threshold=0.25), error=function(e) NULL)
  if (is.null(res) | nrow(res)==0) next
  res$gene <- rownames(res)
  res <- res[!grepl("^(B2m|H2-|Gm|Rik|G5|B5|E5|Hist)", res$gene), ]
  res <- res |> mutate(neg_log10_padj=-log10(p_val_adj), direction=if_else(avg_log2FC>1,"up","down"),
                       candidate=abs(avg_log2FC)>=1 & p_val_adj<0.05,
                       log2FC_adjusted=pmin(avg_log2FC,60)) |> filter(neg_log10_padj<=250)
  top <- bind_rows(res |> filter(direction=="up") |> arrange(p_val_adj) |> head(5),
                   res |> filter(direction=="down") |> arrange(p_val_adj) |> head(5))
  xmax <- min(ceiling(max(abs(res$log2FC_adjusted), na.rm=TRUE)),40)
  p <- res |> tidyplot(x=log2FC_adjusted, y=neg_log10_padj) |>
    add_data_points(data=filter_rows(!candidate), color="lightgrey", rasterize=TRUE) |>
    add_data_points(data=filter_rows(candidate,direction=="up"), color="#FF7777", alpha=0.6) |>
    add_data_points(data=filter_rows(candidate,direction=="down"), color="#7DA8E6", alpha=0.6) |>
    add_reference_lines(x=c(-0.5,0.5), y=-log10(0.05)) |>
    add_data_labels_repel(data=top, label=gene, color="#000000", background=TRUE,
                          box.padding=0.2, point.padding=0.2, min.segment.length=0.05, max.overlaps=50) |>
    adjust_x_axis_title("$Log[2]~fold~change$") |>
    adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") |>
    adjust_legend_position("none") +
    coord_cartesian(xlim=c(-xmax,xmax)) +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=10))
  ggsave(file.path(output_dir,paste0("volcano_plot_",c1,"_vs_",c2,".pdf")), p, width=4, height=4, bg="white", dpi=300)
  message("Saved: ", file.path(output_dir,paste0("22.volcano_plot_",c1,"_vs_",c2,".pdf")))
}


############################### Now the subsetting of data ###############################
#table(seurat_object$sample, seurat_object$seurat_clusters)

############ Now subsetting analysis ###########
# ---------------------------- #
# Step 1: Subset the Seurat Object for Cthrc1+ MyoFBs
# ---------------------------- #

table(seurat_object$sample, seurat_object$seurat_clusters)

subset_seurat <- subset(seurat_object, subset = seurat_clusters == "MyoFB") 
output_directory <- "/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis/MyoFB"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  
}
setwd(output_directory)
# ---------------------------- #
# Step 2: Run PCA on the subset
# ---------------------------- #
subset_seurat <- NormalizeData(subset_seurat)
subset_seurat <- FindVariableFeatures(subset_seurat, selection.method = "vst", nfeatures = 2000)
subset_seurat <- ScaleData(subset_seurat, features = VariableFeatures(subset_seurat))
subset_seurat <- RunPCA(subset_seurat, features = VariableFeatures(subset_seurat))
subset_seurat <- RunUMAP(subset_seurat, dims = 1:20)
subset_seurat <- FindNeighbors(subset_seurat, dims = 1:20)
subset_seurat <- FindClusters(subset_seurat, resolution = 0.4)
#subset_seurat <- subset(subset_seurat, subset = sample != "Bleo_d30")
#table(subset_seurat$sample, subset_seurat$seurat_clusters)
cols <- brewer.pal(10, "Set2")

# Subset UMAP
p1 <- DimPlot(subset_seurat, reduction = "umap", group.by = "seurat_clusters", pt.size = 2, label = T, label.size = 6, cols = cols) + 
  theme_minimal() + theme(panel.grid = element_blank()) + ggtitle("UMAP - MyoFB Subset")
ggsave("UMAP_MyoFB_Set2.png", p1, width = 10, height = 8, dpi = 600)

# Sample-wise UMAPs
u_list <- lapply(as.character(unique(subset_seurat$sample)), function(s) {
  DimPlot(subset_seurat[, subset_seurat$sample == s], reduction = "umap", group.by = "seurat_clusters", 
          label = T, label.size = 5, pt.size = 1, cols = cols) + 
    ggtitle(s) + theme_minimal() + NoLegend() + theme(panel.grid = element_blank())
})
ggsave("umap_by_sample.png", wrap_plots(u_list, ncol = 2), width = 18, height = 15, dpi = 300)

# Proportions
df <- subset_seurat@meta.data %>% count(sample, seurat_clusters) %>% group_by(sample) %>% mutate(p = n/sum(n)*100)
writexl::write_xlsx(df, "cluster_proportions.xlsx")

ggplot(df, aes(sample, p, fill = seurat_clusters)) + geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = cols) + 
  theme_minimal() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("cluster_proportions.pdf", width = 5, height = 6)

##########################################
# ---------------------------- #
# DEGs excel file and later heatmaps
# ---------------------------- #

# DEGs and Excel Export
Idents(subset_seurat) <- "seurat_clusters"
m <- FindAllMarkers(subset_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
m <- m[!grepl("^Neat1|^B2m|^H2-", m$gene), ]

wb <- createWorkbook()
for (n in c(10, 20, 100, 300, 500, 1000)) {
  top_n <- m %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
  sheet <- paste0("Top", n)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, top_n[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")])
}
saveWorkbook(wb, file.path(output_directory, "Top_Markers_DEGs.xlsx"), overwrite = TRUE)

# Seurat Heatmap
top10 <- m %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC) %>% pull(gene) %>% unique()
subset_seurat <- ScaleData(subset_seurat, features = top10)
hp <- DoHeatmap(subset_seurat, features = top10, group.colors = set2_colors) + scale_fill_viridis_c(option = "E")
ggsave(file.path(output_directory, "Subset_heatmap.png"), hp, width = 10, height = 10)

# Sample-Averaged Heatmap (Fixed Dash/Underscore issue)
top6 <- m %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC) %>% pull(gene) %>% unique()
avg_exp <- AggregateExpression(subset_seurat, features = top6, group.by = "sample", return.seurat = F)$RNA
colnames(avg_exp) <- gsub("-", "_", colnames(avg_exp))

# Order samples and scale
s_cols <- c("Saline" = "darkseagreen", "Bleo_d14" = "deepskyblue", "Bleo_d30" = "royalblue", "Bleo_d60" = "orchid")
avg_exp <- avg_exp[, intersect(names(s_cols), colnames(avg_exp))]
scaled_exp <- t(scale(t(avg_exp)))

ann <- data.frame(sample = colnames(scaled_exp), row.names = colnames(scaled_exp))
pdf("Subset_heatmap_sample_expression.pdf", width = 8, height = 12)
pheatmap(scaled_exp, annotation_col = ann, annotation_colors = list(sample = s_cols), 
         cluster_cols = F, scale = "none", border_color = "gray", cellwidth = 20, cellheight = 15)
dev.off()

###########################################


subset_seurat<-subset(seurat_object,subset=seurat_clusters=="AF1")
output_directory<-"/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis/AF1";if(!dir.exists(output_directory))dir.create(output_directory,T);setwd(output_directory)

subset_seurat <- NormalizeData(subset_seurat) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=2000) %>%
  ScaleData(features=VariableFeatures(.)) %>%
  RunPCA(features=VariableFeatures(.))

subset_seurat <- subset_seurat %>%
  RunUMAP(dims=1:15) %>%
  FindNeighbors(dims=1:20) %>%
  FindClusters(resolution=0.4)

umap_coords <- Embeddings(subset_seurat, "umap")
umap_coords[, 2] <- -umap_coords[, 2]
subset_seurat[["umap"]]@cell.embeddings <- umap_coords

set2_cols <- brewer.pal(10, "Paired")

p1 <- DimPlot(subset_seurat, reduction = "umap", group.by = "seurat_clusters", pt.size = 2, label = T, label.size = 6, cols = set2_cols) + 
  theme_minimal() + theme(panel.grid = element_blank()) + ggtitle("UMAP - AF1 Subset")
ggsave("UMAP_AF1_Set2.png", p1, width = 8, height = 6.5, dpi = 600)

# Sample-wise UMAPs
u_list <- lapply(as.character(unique(subset_seurat$sample)), function(s) {
  DimPlot(subset_seurat[, subset_seurat$sample == s], reduction = "umap", group.by = "seurat_clusters", 
          label = T, label.size = 5, pt.size = 1, cols = set2_cols) + 
    ggtitle(s) + theme_minimal() + NoLegend() + theme(panel.grid = element_blank())
})
ggsave("umap_by_sample.png", wrap_plots(u_list, ncol = 2), width = 18, height = 15, dpi = 300)

# Proportions
df <- subset_seurat@meta.data %>% count(sample, seurat_clusters) %>% group_by(sample) %>% mutate(p = n/sum(n)*100)
writexl::write_xlsx(df, "cluster_proportions.xlsx")

ggplot(df, aes(sample, p, fill = seurat_clusters)) + geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = set2_cols) + 
  theme_minimal() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("cluster_proportions.pdf", width = 5, height = 6)

# DEGs and Excel Export
Idents(subset_seurat) <- "seurat_clusters"
m <- FindAllMarkers(subset_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
m <- m[!grepl("^Neat1|^B2m|^H2-", m$gene), ]

wb <- createWorkbook()
for (n in c(10, 20, 100, 300, 500, 1000)) {
  top_n <- m %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
  sheet <- paste0("Top", n)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, top_n[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")])
}
saveWorkbook(wb, "Top_Markers_DEGs.xlsx", overwrite = TRUE)

# Heatmaps
top10 <- m %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC) %>% pull(gene) %>% unique()
subset_seurat <- ScaleData(subset_seurat, features = top10)
hp <- DoHeatmap(subset_seurat, features = top10, group.colors = set2_cols) + scale_fill_viridis_c(option = "E")
ggsave("Subset_heatmap.png", hp, width = 10, height = 10)

top6 <- m %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC) %>% pull(gene) %>% unique()
avg_exp <- AggregateExpression(subset_seurat, features = top6, group.by = "sample", return.seurat = F)$RNA
colnames(avg_exp) <- gsub("-", "_", colnames(avg_exp))

s_cols <- c("Saline" = "darkseagreen", "Bleo_d14" = "deepskyblue", "Bleo_d30" = "royalblue", "Bleo_d60" = "orchid")
avg_exp <- avg_exp[, intersect(names(s_cols), colnames(avg_exp))]
scaled_exp <- t(scale(t(avg_exp)))

pdf("Subset_heatmap_sample_expression.pdf", width = 8, height = 12)
pheatmap(scaled_exp, annotation_col = data.frame(sample = colnames(scaled_exp), row.names = colnames(scaled_exp)), 
         annotation_colors = list(sample = s_cols), cluster_cols = F, scale = "none", border_color = "gray", cellwidth = 20, cellheight = 15)
dev.off()

############################################################################

################## For Conversion back to scanpy h5ad for any downstream velocity measurements 

############################################################################

table(subset_seurat$seurat_clusters)
subset_seurat$seurat_clusters <- subset_seurat$seurat_clusters
Idents(subset_seurat) <- subset_seurat$seurat_clusters

table(Idents(subset_seurat))
table(subset_seurat$seurat_clusters)


save_dir <- "/Users/aly/Documents/2026_Adamts4_Manuscript/Analysis/MyoFB/Conversion"
dir.create(save_dir, recursive = T, showWarnings = F)

fwrite(data.frame(gene = rownames(subset_seurat)), file.path(save_dir, "genes.csv"))
fwrite(data.frame(barcode = colnames(subset_seurat)), file.path(save_dir, "barcodes.csv"))
fwrite(subset_seurat@meta.data, file.path(save_dir, "cell_metadata.csv"), row.names = T)

if ("umap" %in% names(subset_seurat@reductions)) {
  fwrite(Embeddings(subset_seurat, "umap"), file.path(save_dir, "umap.csv"), row.names = T)
}

cols <- c("0"="#A6CEE3","1"="#1F78B4","2"="#B2DF8A","3"="#33A02C","4"="#FB9A99",
          "5"="#E31A1C","6"="#FDBF6F","7"="#FF7F00","8"="#CAB2D6","9"="#6A3D9A")
fwrite(data.frame(cluster=names(cols), color=cols), file.path(save_dir, "cluster_colors.csv"))

for (lyr in c("spliced", "unspliced", "ambiguous")) {
  if (lyr %in% names(subset_seurat@assays)) {
    mat <- GetAssayData(subset_seurat, assay = lyr, layer = "counts")
    
    missing_g <- setdiff(rownames(subset_seurat), rownames(mat))
    if (length(missing_g) > 0) {
      mat <- rbind(mat, Matrix(0, length(missing_g), ncol(mat), dimnames=list(missing_g, colnames(mat)), sparse=T))
    }
    
    missing_c <- setdiff(colnames(subset_seurat), colnames(mat))
    if (length(missing_c) > 0) {
      mat <- cbind(mat, Matrix(0, nrow(mat), length(missing_c), dimnames=list(rownames(mat), missing_c), sparse=T))
    }
    
    mat <- mat[rownames(subset_seurat), colnames(subset_seurat)]
    writeMM(mat, file.path(save_dir, paste0(lyr, ".mtx")))
  }
}
############################################################################
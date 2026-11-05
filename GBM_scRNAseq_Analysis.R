# Youssef (Stay strong for yourself)
# Glioblastoma Single-Cell RNA-seq Analysis
# GSE162931 | GBM131 and GBM022 | Naive vs RT 2 Days
# Run sections sequentially. Only four samples are included.

# 1. Load Libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(scDblFinder)
library(SingleR)
library(celldex)
library(harmony)
library(miloR)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)

# 2. Load Four GBM Samples
# Download and extract the four original GEO sample matrices into data/extracted/.
# Each directory must contain matrix.mtx, barcodes.tsv, and features.tsv/genes.tsv
# (possibly .gz), directly inside the directory.

# GBM131 Naive
gbm131_naive_counts <- Read10X(data.dir = "data/extracted/GBM131_Naive")
gbm131_naive <- CreateSeuratObject(counts = gbm131_naive_counts, project = "GBM131_Naive")

# GBM131 RT 2 Days
gbm131_rt_counts <- Read10X(data.dir = "data/extracted/GBM131_RT2d")
gbm131_rt <- CreateSeuratObject(counts = gbm131_rt_counts, project = "GBM131_RT2d")

# GBM022 Naive
gbm022_naive_counts <- Read10X(data.dir = "data/extracted/GBM022_Naive")
gbm022_naive <- CreateSeuratObject(counts = gbm022_naive_counts, project = "GBM022_Naive")

# GBM022 RT 2 Days
# If the files are not extracted yet, uncomment the next two lines.
# dir.create("data/extracted/GBM022_RT2d", recursive = TRUE, showWarnings = FALSE)
# untar("data/raw/GSM4967244.tar.gz", exdir = "data/extracted/GBM022_RT2d")
# Check file layout before Read10X:
list.files("data/extracted/GBM022_RT2d", recursive = TRUE)
gbm022_rt_counts <- Read10X(data.dir = "data/extracted/GBM022_RT2d")
gbm022_rt <- CreateSeuratObject(counts = gbm022_rt_counts, project = "GBM022_RT2d")

# Verify All Four Samples
gbm131_naive
gbm131_rt
gbm022_naive
gbm022_rt
# ========================================================================
# 3. Quality Control

# Mitochondrial Percentage - gbm131_naive
gbm131_naive[["percent.mt"]] <- PercentageFeatureSet(gbm131_naive, pattern = "^MT-")

# Mitochondrial Percentage - gbm131_rt
gbm131_rt[["percent.mt"]] <- PercentageFeatureSet(gbm131_rt, pattern = "^MT-")

# Mitochondrial Percentage - gbm022_naive
gbm022_naive[["percent.mt"]] <- PercentageFeatureSet(gbm022_naive, pattern = "^MT-")

# Mitochondrial Percentage - gbm022_rt
gbm022_rt[["percent.mt"]] <- PercentageFeatureSet(gbm022_rt, pattern = "^MT-")
# ############################################
# 3.1. QC Visualization before Filtering
# QC Violin Plot - gbm131_naive
VlnPlot(gbm131_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC Violin Plot - gbm131_rt
VlnPlot(gbm131_rt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC Violin Plot - gbm022_naive
VlnPlot(gbm022_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC Violin Plot - gbm022_rt
VlnPlot(gbm022_rt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ===========================================================================
# 4. Cell Filtering

# Filter gbm131_naive
gbm131_naive <- subset(gbm131_naive, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Filter gbm131_rt
gbm131_rt <- subset(gbm131_rt, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Filter gbm022_naive
gbm022_naive <- subset(gbm022_naive, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Filter gbm022_rt
gbm022_rt <- subset(gbm022_rt, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#           =====================================================================
# 4.1.QC Visualization After Filtering
# GBM131 Naive
VlnPlot(gbm131_naive,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# GBM131 RT
VlnPlot(gbm131_rt,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# GBM022 Naive
VlnPlot(gbm022_naive,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# GBM022 RT
VlnPlot(gbm022_rt,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# =================================================================================
# 5. Doublet Detection
library(scDblFinder) # Doublet detection
library(SingleCellExperiment) # Single-cell data structures

# A. Convert gbm131_naive to SingleCellExperiment
sce_gbm131_naive <- as.SingleCellExperiment(gbm131_naive)

# A. Convert gbm131_rt to SingleCellExperiment
sce_gbm131_rt <- as.SingleCellExperiment(gbm131_rt)

# A. Convert gbm022_naive to SingleCellExperiment
sce_gbm022_naive <- as.SingleCellExperiment(gbm022_naive)

# A. Convert gbm022_rt to SingleCellExperiment
sce_gbm022_rt <- as.SingleCellExperiment(gbm022_rt)
# #################################################
# B. Run scDblFinder - gbm131_naive
set.seed(100)
sce_gbm131_naive <- scDblFinder(sce_gbm131_naive)

# B. Run scDblFinder - gbm131_rt
set.seed(100)
sce_gbm131_rt <- scDblFinder(sce_gbm131_rt)

# B. Run scDblFinder - gbm022_naive
set.seed(100)
sce_gbm022_naive <- scDblFinder(sce_gbm022_naive)

# B. Run scDblFinder - gbm022_rt
set.seed(100)
sce_gbm022_rt <- scDblFinder(sce_gbm022_rt)

set.seed(100)


# C. Add Doublet Results - gbm131_naive
gbm131_naive$doublet_score <- sce_gbm131_naive$scDblFinder.score
gbm131_naive$doublet_class <- sce_gbm131_naive$scDblFinder.class
table(gbm131_naive$doublet_class)

# C. Add Doublet Results - gbm131_rt
gbm131_rt$doublet_score <- sce_gbm131_rt$scDblFinder.score
gbm131_rt$doublet_class <- sce_gbm131_rt$scDblFinder.class
table(gbm131_rt$doublet_class)

# C. Add Doublet Results - gbm022_naive
gbm022_naive$doublet_score <- sce_gbm022_naive$scDblFinder.score
gbm022_naive$doublet_class <- sce_gbm022_naive$scDblFinder.class
table(gbm022_naive$doublet_class)

# C. Add Doublet Results - gbm022_rt
gbm022_rt$doublet_score <- sce_gbm022_rt$scDblFinder.score
gbm022_rt$doublet_class <- sce_gbm022_rt$scDblFinder.class
table(gbm022_rt$doublet_class)
# ========================================================================

# Doublets detection and visualization 

library(ggplot2)

# 1. GBM131 Naive

table(gbm131_naive$doublet_class)

ggplot(gbm131_naive@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM131 Naive",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

VlnPlot(gbm131_naive,
        features = "doublet_score",
        group.by = "doublet_class")


# 2. GBM131 RT 2 Days

table(gbm131_rt$doublet_class)

ggplot(gbm131_rt@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM131 RT 2 Days",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

VlnPlot(gbm131_rt,
        features = "doublet_score",
        group.by = "doublet_class")


# 3. GBM022 Naive

table(gbm022_naive$doublet_class)

ggplot(gbm022_naive@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM022 Naive",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

VlnPlot(gbm022_naive,
        features = "doublet_score",
        group.by = "doublet_class")


# 4. GBM022 RT 2 Days

table(gbm022_rt$doublet_class)

ggplot(gbm022_rt@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM022 RT 2 Days",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

VlnPlot(gbm022_rt,
        features = "doublet_score",
        group.by = "doublet_class")

# #################################################
# D. Remove Doublets - gbm131_naive
gbm131_naive <- subset(gbm131_naive, subset = doublet_class == "singlet")

# D. Remove Doublets - gbm131_rt
gbm131_rt <- subset(gbm131_rt, subset = doublet_class == "singlet")

# D. Remove Doublets - gbm022_naive
gbm022_naive <- subset(gbm022_naive, subset = doublet_class == "singlet")

# D. Remove Doublets - gbm022_rt
gbm022_rt <- subset(gbm022_rt, subset = doublet_class == "singlet")

# ########################################################

library(Seurat)
library(ggplot2)

# 1. Doublet Score Scatter Plots After Removal

# GBM131 Naive
ggplot(gbm131_naive@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM131 Naive - After Doublet Removal",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

# GBM131 RT
ggplot(gbm131_rt@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM131 RT - After Doublet Removal",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

# GBM022 Naive
ggplot(gbm022_naive@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM022 Naive - After Doublet Removal",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

# GBM022 RT
ggplot(gbm022_rt@meta.data,
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "GBM022 RT - After Doublet Removal",
       x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()

# ########################################
# 2. Violin Plots After Doublet Removal

# GBM131 Naive
VlnPlot(gbm131_naive,
        features = "doublet_score",
        group.by = "doublet_class") +
  ggtitle("GBM131 Naive - After Doublet Removal")

# GBM131 RT
VlnPlot(gbm131_rt,
        features = "doublet_score",
        group.by = "doublet_class") +
  ggtitle("GBM131 RT - After Doublet Removal")

# GBM022 Naive
VlnPlot(gbm022_naive,
        features = "doublet_score",
        group.by = "doublet_class") +
  ggtitle("GBM022 Naive - After Doublet Removal")

# GBM022 RT
VlnPlot(gbm022_rt,
        features = "doublet_score",
        group.by = "doublet_class") +
  ggtitle("GBM022 RT - After Doublet Removal")

# #################################
# 3. Verify Doublet Removal

table(gbm131_naive$doublet_class)
table(gbm131_rt$doublet_class)
table(gbm022_naive$doublet_class)
table(gbm022_rt$doublet_class)

# #################################
# 4. Check Remaining Cell Numbers

cell_counts <- data.frame(
  Sample = c("GBM131_Naive", "GBM131_RT",
             "GBM022_Naive", "GBM022_RT"),
  Remaining_Cells = c(
    ncol(gbm131_naive),
    ncol(gbm131_rt),
    ncol(gbm022_naive),
    ncol(gbm022_rt)
  )
)

print(cell_counts)

ggplot(cell_counts,
       aes(x = Sample, y = Remaining_Cells, fill = Sample)) +
  geom_col() +
  labs(title = "Remaining Cells After Doublet Removal",
       x = "Sample", y = "Number of Cells") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ================================================================================

# 6. Normalization

# Normalize gbm131_naive
gbm131_naive <- NormalizeData(gbm131_naive)

# Normalize gbm131_rt
gbm131_rt <- NormalizeData(gbm131_rt)

# Normalize gbm022_naive
gbm022_naive <- NormalizeData(gbm022_naive)

# Normalize gbm022_rt
gbm022_rt <- NormalizeData(gbm022_rt)

# 7. Find Variable Features

# Variable Genes - gbm131_naive
gbm131_naive <- FindVariableFeatures(gbm131_naive, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(gbm131_naive)

# Variable Genes - gbm131_rt
gbm131_rt <- FindVariableFeatures(gbm131_rt, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(gbm131_rt)

# Variable Genes - gbm022_naive
gbm022_naive <- FindVariableFeatures(gbm022_naive, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(gbm022_naive)

# Variable Genes - gbm022_rt
gbm022_rt <- FindVariableFeatures(gbm022_rt, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(gbm022_rt)
#              =========================================
library(Seurat)
library(ggplot2)

# 1. GBM131 Naive

gbm131_naive <- FindVariableFeatures(gbm131_naive, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(gbm131_naive)

head(VariableFeatures(gbm131_naive), 10)

LabelPoints(
  plot = VariableFeaturePlot(gbm131_naive),
  points = head(VariableFeatures(gbm131_naive), 10),
  repel = TRUE
)

# 2. GBM131 RT 2 Days

gbm131_rt <- FindVariableFeatures(gbm131_rt, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(gbm131_rt)

head(VariableFeatures(gbm131_rt), 10)

LabelPoints(
  plot = VariableFeaturePlot(gbm131_rt),
  points = head(VariableFeatures(gbm131_rt), 10),
  repel = TRUE
)

# 3. GBM022 Naive

gbm022_naive <- FindVariableFeatures(gbm022_naive, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(gbm022_naive)

head(VariableFeatures(gbm022_naive), 10)

LabelPoints(
  plot = VariableFeaturePlot(gbm022_naive),
  points = head(VariableFeatures(gbm022_naive), 10),
  repel = TRUE
)

# 4. GBM022 RT 2 Days

gbm022_rt <- FindVariableFeatures(gbm022_rt, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(gbm022_rt)

head(VariableFeatures(gbm022_rt), 10)

LabelPoints(
  plot = VariableFeaturePlot(gbm022_rt),
  points = head(VariableFeatures(gbm022_rt), 10),
  repel = TRUE
)


# 8. Scaling

# Scale gbm131_naive
gbm131_naive <- ScaleData(gbm131_naive, features = VariableFeatures(gbm131_naive))

# Scale gbm131_rt
gbm131_rt <- ScaleData(gbm131_rt, features = VariableFeatures(gbm131_rt))

# Scale gbm022_naive
gbm022_naive <- ScaleData(gbm022_naive, features = VariableFeatures(gbm022_naive))

# Scale gbm022_rt
gbm022_rt <- ScaleData(gbm022_rt, features = VariableFeatures(gbm022_rt))

# 9. PCA

# PCA - gbm131_naive
gbm131_naive <- RunPCA(gbm131_naive, features = VariableFeatures(gbm131_naive))
ElbowPlot(gbm131_naive)

# PCA - gbm131_rt
gbm131_rt <- RunPCA(gbm131_rt, features = VariableFeatures(gbm131_rt))
ElbowPlot(gbm131_rt)

# PCA - gbm022_naive
gbm022_naive <- RunPCA(gbm022_naive, features = VariableFeatures(gbm022_naive))
ElbowPlot(gbm022_naive)

# PCA - gbm022_rt
gbm022_rt <- RunPCA(gbm022_rt, features = VariableFeatures(gbm022_rt))
ElbowPlot(gbm022_rt)

# 10. Clustering and UMAP

# Clustering - gbm131_naive
gbm131_naive <- FindNeighbors(gbm131_naive, dims = 1:20)
gbm131_naive <- FindClusters(gbm131_naive, resolution = 0.5)
gbm131_naive <- RunUMAP(gbm131_naive, dims = 1:20)
DimPlot(gbm131_naive, reduction = "umap", label = TRUE)

# Clustering - gbm131_rt
gbm131_rt <- FindNeighbors(gbm131_rt, dims = 1:20)
gbm131_rt <- FindClusters(gbm131_rt, resolution = 0.5)
gbm131_rt <- RunUMAP(gbm131_rt, dims = 1:20)
DimPlot(gbm131_rt, reduction = "umap", label = TRUE)

# Clustering - gbm022_naive
gbm022_naive <- FindNeighbors(gbm022_naive, dims = 1:20)
gbm022_naive <- FindClusters(gbm022_naive, resolution = 0.5)
gbm022_naive <- RunUMAP(gbm022_naive, dims = 1:20)
DimPlot(gbm022_naive, reduction = "umap", label = TRUE)

# Clustering - gbm022_rt
gbm022_rt <- FindNeighbors(gbm022_rt, dims = 1:20)
gbm022_rt <- FindClusters(gbm022_rt, resolution = 0.5)
gbm022_rt <- RunUMAP(gbm022_rt, dims = 1:20)
DimPlot(gbm022_rt, reduction = "umap", label = TRUE)
# #=======================================================

# 1. GBM131 Naive

DimHeatmap(gbm131_naive, dims = 1:10, cells = 500, balanced = TRUE)

gbm131_naive <- FindNeighbors(gbm131_naive, dims = 1:10)

gbm131_naive <- FindClusters(gbm131_naive, resolution = 0.5)

gbm131_naive <- RunUMAP(gbm131_naive, dims = 1:10)


# 2. GBM131 RT 2 Days

DimHeatmap(gbm131_rt, dims = 1:10, cells = 500, balanced = TRUE)

gbm131_rt <- FindNeighbors(gbm131_rt, dims = 1:10)

gbm131_rt <- FindClusters(gbm131_rt, resolution = 0.5)

# 3. GBM022 Naive

DimHeatmap(gbm022_naive, dims = 1:10, cells = 500, balanced = TRUE)

gbm022_naive <- FindNeighbors(gbm022_naive, dims = 1:10)

gbm022_naive <- FindClusters(gbm022_naive, resolution = 0.5)


# 4. GBM022 RT 2 Days

DimHeatmap(gbm022_rt, dims = 1:10, cells = 500, balanced = TRUE)

gbm022_rt <- FindNeighbors(gbm022_rt, dims = 1:10)

gbm022_rt <- FindClusters(gbm022_rt, resolution = 0.5)

# ########################################################################
# 11. SingleR Annotation Before Integration
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

# SingleR - gbm131_naive
pred_gbm131_naive <- SingleR(test = as.SingleCellExperiment(gbm131_naive), ref = hpca_ref, labels = hpca_ref$label.main)
gbm131_naive$SingleR_annotation <- pred_gbm131_naive$labels
DimPlot(gbm131_naive, group.by = "SingleR_annotation", label = TRUE, repel = TRUE)

# SingleR - gbm131_rt
pred_gbm131_rt <- SingleR(test = as.SingleCellExperiment(gbm131_rt), ref = hpca_ref, labels = hpca_ref$label.main)
gbm131_rt$SingleR_annotation <- pred_gbm131_rt$labels
DimPlot(gbm131_rt, group.by = "SingleR_annotation", label = TRUE, repel = TRUE)

# SingleR - gbm022_naive
pred_gbm022_naive <- SingleR(test = as.SingleCellExperiment(gbm022_naive), ref = hpca_ref, labels = hpca_ref$label.main)
gbm022_naive$SingleR_annotation <- pred_gbm022_naive$labels
DimPlot(gbm022_naive, group.by = "SingleR_annotation", label = TRUE, repel = TRUE)

# SingleR - gbm022_rt
pred_gbm022_rt <- SingleR(test = as.SingleCellExperiment(gbm022_rt), ref = hpca_ref, labels = hpca_ref$label.main)
gbm022_rt$SingleR_annotation <- pred_gbm022_rt$labels
DimPlot(gbm022_rt, group.by = "SingleR_annotation", label = TRUE, repel = TRUE)

# 12. Add Sample Information

# Metadata - gbm131_naive
gbm131_naive$sample_id <- "GBM131_Naive"
gbm131_naive$patient_id <- "GBM131"
gbm131_naive$condition <- "Naive"

# Metadata - gbm131_rt
gbm131_rt$sample_id <- "GBM131_RT2d"
gbm131_rt$patient_id <- "GBM131"
gbm131_rt$condition <- "RT"

# Metadata - gbm022_naive
gbm022_naive$sample_id <- "GBM022_Naive"
gbm022_naive$patient_id <- "GBM022"
gbm022_naive$condition <- "Naive"

# Metadata - gbm022_rt
gbm022_rt$sample_id <- "GBM022_RT2d"
gbm022_rt$patient_id <- "GBM022"
gbm022_rt$condition <- "RT"

# 13. Merge GBM Samples
merged_gbm <- merge(gbm131_naive,
                    y = list(gbm131_rt, gbm022_naive, gbm022_rt),
                    add.cell.ids = c("131N", "131RT", "022N", "022RT"))

# 14. Normalize Merged Data
merged_gbm <- JoinLayers(merged_gbm)
merged_gbm <- NormalizeData(merged_gbm)
merged_gbm <- FindVariableFeatures(merged_gbm, nfeatures = 2000)
merged_gbm <- ScaleData(merged_gbm)
merged_gbm <- RunPCA(merged_gbm)

# 15. Harmony Integration
# Correct patient-associated variation; inspect whether biology is preserved.
merged_gbm <- RunHarmony(merged_gbm, group.by.vars = "patient_id")

# 16. Clustering After Integration
merged_gbm <- FindNeighbors(merged_gbm, reduction = "harmony", dims = 1:20)
merged_gbm <- FindClusters(merged_gbm, resolution = 0.5)
merged_gbm <- RunUMAP(merged_gbm, reduction = "harmony", dims = 1:20)

# 17. Integrated UMAP
DimPlot(merged_gbm, group.by = "sample_id")
DimPlot(merged_gbm, group.by = "condition")
DimPlot(merged_gbm, label = TRUE)

# 18. SingleR Annotation After Integration
# Reference labels are exploratory for malignant cultured GBM cells.
merged_gbm <- JoinLayers(merged_gbm)
ref_pred <- SingleR(test = as.SingleCellExperiment(merged_gbm), ref = hpca_ref, labels = hpca_ref$label.main)
merged_gbm$SingleR_annotation <- ref_pred$labels
DimPlot(merged_gbm, group.by = "SingleR_annotation", label = TRUE, repel = TRUE)

# 19. Cluster Marker Genes
Idents(merged_gbm) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(merged_gbm, only.pos = TRUE)
head(cluster_markers)
write.csv(cluster_markers, "GBM_cluster_markers.csv", row.names = FALSE)
# #=======================================================================
# 20. GBM State Markers - DotPlot, FeaturePlot, Violin Plot, Heatmap
DefaultAssay(merged_gbm) <- 'RNA'
merged_gbm <- JoinLayers(merged_gbm, assay = 'RNA')
gbm_genes <- c('SOX2','OLIG2','PDGFRA','EGFR','CD44','CHI3L1','VIM','NES')
gbm_genes <- intersect(gbm_genes, rownames(merged_gbm))
if (length(gbm_genes) == 0) stop('None of the GBM marker genes were found. Check gene identifiers.')
p20a <- DotPlot(merged_gbm, features = gbm_genes, group.by = 'seurat_clusters') + RotatedAxis() + ggtitle('GBM state markers by cluster')
print(p20a)
ggsave('GBM_plots/20_GBM_marker_DotPlot.png', p20a, width=11, height=6, dpi=300)
p20b <- FeaturePlot(merged_gbm, features = gbm_genes, reduction = 'umap', ncol=4)
print(p20b)
ggsave('GBM_plots/20_GBM_marker_FeaturePlot.png', p20b, width=15, height=9, dpi=300)
p20c <- VlnPlot(merged_gbm, features = gbm_genes[1:min(4,length(gbm_genes))], group.by = 'condition', ncol=2, pt.size=0)
print(p20c)
ggsave('GBM_plots/20_GBM_marker_Violin.png', p20c, width=12, height=8, dpi=300)

# 20A. SingleR Annotation UMAP - reference-style labeled plots
# IMPORTANT: HPCA annotations are NOT verified cell identities in cultured GBM.
# This section uses existing SingleR metadata; it does not invent labels.
annotation_columns <- intersect(c('SingleR_annotation','SingleR_annotation_integrated','old_annotation','new_annotation'),
                                colnames(merged_gbm@meta.data))
if (length(annotation_columns) > 0) {
  for (annotation_col in annotation_columns) {
    p_ann <- DimPlot(merged_gbm, reduction='umap', group.by=annotation_col,
                     label=TRUE, repel=TRUE, label.size=3) +
      ggtitle(paste('GBM integrated UMAP -', annotation_col))
    print(p_ann)
    ggsave(paste0('GBM_plots/20A_SingleR_',annotation_col,'.png'),p_ann,
           width=15,height=9,dpi=300)
  }
  if (length(annotation_columns) >= 2) {
    p_old <- DimPlot(merged_gbm,reduction='umap',group.by=annotation_columns[1],
                     label=TRUE,repel=TRUE,label.size=2.5) + ggtitle(annotation_columns[1])
    p_new <- DimPlot(merged_gbm,reduction='umap',group.by=annotation_columns[2],
                     label=TRUE,repel=TRUE,label.size=2.5) + ggtitle(annotation_columns[2])
    p_compare <- p_old + p_new + patchwork::plot_layout(ncol=2)
    print(p_compare)
    ggsave('GBM_plots/20A_SingleR_Before_After_Comparison.png',p_compare,
           width=19,height=8,dpi=300)
  }
} else {
  message('20A: No SingleR annotation metadata found. Run SingleR annotation first.')
}

# 21. Differential Abundance - MiloR and neighborhood visualizations
# Check the reduction name in the converted SingleCellExperiment object.
gbm_sce <- as.SingleCellExperiment(merged_gbm)
print(reducedDimNames(gbm_sce))
# Explicitly transfer Harmony coordinates with a predictable name.
reducedDim(gbm_sce, 'HARMONY') <- Embeddings(merged_gbm, 'harmony')
gbm_milo <- Milo(gbm_sce)
gbm_milo <- buildGraph(gbm_milo, k=30, d=20, reduced.dim='HARMONY')
gbm_milo <- makeNhoods(gbm_milo, prop=0.1, k=30, d=20, refined=TRUE, reduced_dims='HARMONY')
gbm_milo <- countCells(gbm_milo, meta.data=as.data.frame(colData(gbm_milo)), sample='sample_id')
gbm_milo <- calcNhoodDistance(gbm_milo, d=20, reduced.dim='HARMONY')
# Display neighborhoods and their size distributions.
p21 <- plotNhoodSizeHist(gbm_milo)
print(p21)
ggsave('GBM_plots/21_Milo_neighborhood_sizes.png', p21, width=8, height=5, dpi=300)

# 22. MiloR Design and cell-count composition
milo_design <- unique(merged_gbm@meta.data[, c('sample_id','patient_id','condition')])
rownames(milo_design) <- milo_design$sample_id
milo_design$patient_id <- factor(milo_design$patient_id)
milo_design$condition <- factor(milo_design$condition, levels=c('Naive','RT'))
print(milo_design)
cell_composition <- as.data.frame(table(merged_gbm$sample_id, merged_gbm$seurat_clusters))
colnames(cell_composition) <- c('sample_id','cluster','cells')
p22 <- ggplot(cell_composition, aes(x=sample_id, y=cells, fill=cluster)) + geom_col(position='fill') +
  scale_y_continuous(labels=scales::percent) + labs(x='Sample',y='Fraction of cells',fill='Cluster',title='Cluster composition per sample') +
  theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1))
print(p22)
ggsave('GBM_plots/22_Cluster_composition.png', p22, width=9, height=6, dpi=300)

# 23. Differential Abundance test and MiloR visualizations
# With two patient pairs, treat results as exploratory; significance is unstable.
milo_results <- testNhoods(gbm_milo, design=~patient_id+condition, design.df=milo_design)
print(head(milo_results))
write.csv(milo_results, 'GBM_MiloR_results.csv', row.names=FALSE)
p23a <- ggplot(as.data.frame(milo_results), aes(x=logFC, y=-log10(pmax(SpatialFDR, 1e-300)))) +
  geom_point(alpha=0.6, size=1) + theme_classic() +
  labs(title='Milo differential abundance (exploratory)',x='Neighborhood log2 fold change',y='-log10 Spatial FDR')
print(p23a)
ggsave('GBM_plots/23_Milo_DA_scatter.png', p23a, width=8, height=6, dpi=300)
# plotNhoodGraphDA requires UMAP coordinates in the Milo object.
reducedDim(gbm_milo, 'UMAP') <- Embeddings(merged_gbm, 'umap')
p23b <- plotNhoodGraphDA(gbm_milo, milo_results, layout='UMAP', alpha=0.1)
print(p23b)
ggsave('GBM_plots/23_Milo_DA_neighborhood_graph.png', p23b, width=9, height=7, dpi=300)

# 24. Pseudobulk Counts - sample library size visualization
Idents(merged_gbm) <- 'sample_id'
pb_counts <- AggregateExpression(merged_gbm, assays='RNA', group.by='sample_id', return.seurat=FALSE)$RNA
print(colnames(pb_counts))
pb_lib <- data.frame(sample_id=colnames(pb_counts), total_counts=colSums(pb_counts))
p24 <- ggplot(pb_lib, aes(x=sample_id,y=total_counts)) + geom_col(fill='steelblue') +
  theme_classic() + labs(title='Pseudobulk library sizes',x='Sample',y='Total UMI counts') +
  theme(axis.text.x=element_text(angle=45,hjust=1))
print(p24)
ggsave('GBM_plots/24_Pseudobulk_library_sizes.png', p24, width=8, height=5, dpi=300)

# 25. Paired Differential Expression - edgeR (exploratory)
pb_design <- unique(merged_gbm@meta.data[, c('sample_id','patient_id','condition')])
pb_design <- pb_design[match(colnames(pb_counts), pb_design$sample_id), , drop=FALSE]
stopifnot(!anyNA(pb_design$sample_id), nrow(pb_design)==ncol(pb_counts))
pb_design$patient_id <- factor(pb_design$patient_id)
pb_design$condition <- factor(pb_design$condition, levels=c('Naive','RT'))
dge <- DGEList(counts=pb_counts)
keep <- filterByExpr(dge, group=pb_design$condition)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~patient_id+condition, data=pb_design)
# Very limited residual df with only two patient pairs; interpret cautiously.
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef='conditionRT')
gbm_de <- topTags(qlf, n=Inf)$table
gbm_de$gene <- rownames(gbm_de)
write.csv(gbm_de, 'GBM_RT_vs_Naive_DEG.csv', row.names=FALSE)
p25 <- plotMDS(dge, labels=pb_design$sample_id, col=as.integer(pb_design$condition), plot=FALSE)
png('GBM_plots/25_Pseudobulk_MDS.png', width=1600, height=1200, res=200)
plot(p25$x, p25$y, pch=19, col=as.integer(pb_design$condition), xlab='MDS 1', ylab='MDS 2', main='Pseudobulk MDS by sample')
text(p25$x, p25$y, labels=pb_design$sample_id, pos=3, cex=0.75)
legend('topright', legend=levels(pb_design$condition), col=seq_along(levels(pb_design$condition)), pch=19)
dev.off()

# 26. Differential Expression - Volcano and MA plots
gbm_de$significance <- 'Not Significant'
gbm_de$significance[!is.na(gbm_de$FDR) & gbm_de$FDR<0.05 & gbm_de$logFC>0.25] <- 'Upregulated'
gbm_de$significance[!is.na(gbm_de$FDR) & gbm_de$FDR<0.05 & gbm_de$logFC< -0.25] <- 'Downregulated'
p26a <- ggplot(gbm_de, aes(x=logFC, y=-log10(pmax(FDR,1e-300)), color=significance)) +
  geom_point(size=1,alpha=0.6) + geom_vline(xintercept=c(-0.25,0.25),linetype='dashed') +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') + theme_classic() +
  labs(title='RT vs Naive volcano (exploratory)',x='log2 fold change',y='-log10 FDR')
print(p26a)
ggsave('GBM_plots/26_Volcano.png',p26a,width=9,height=6,dpi=300)
p26b <- ggplot(gbm_de, aes(x=logCPM,y=logFC,color=significance)) + geom_point(size=1,alpha=0.6) +
  geom_hline(yintercept=0,linetype='dashed') + theme_classic() +
  labs(title='MA plot: RT vs Naive',x='Average log2 CPM',y='log2 fold change')
print(p26b)
ggsave('GBM_plots/26_MA_plot.png',p26b,width=9,height=6,dpi=300)

# 26A. Labeled Volcano Plot - reference-style gene labels
# Use FDR when available. With only two patient pairs, treat labels as exploratory.
if (requireNamespace('ggrepel',quietly=TRUE)) {
  label_candidates <- gbm_de[is.finite(gbm_de$PValue) & is.finite(gbm_de$logFC),,drop=FALSE]
  label_candidates <- label_candidates[order(label_candidates$PValue),,drop=FALSE]
  label_candidates <- head(label_candidates,20)
  p26_labeled <- p26a +
    ggrepel::geom_text_repel(data=label_candidates,
                             aes(x=logFC,y=-log10(pmax(FDR,1e-300)),label=gene),
                             size=3,max.overlaps=Inf,show.legend=FALSE) +
    ggtitle('GBM RT vs Naive - Labeled Volcano (exploratory)')
  print(p26_labeled)
  ggsave('GBM_plots/26A_Labeled_Volcano.png',p26_labeled,width=12,height=8,dpi=300)
} else message('26A: Install ggrepel for gene labels: install.packages("ggrepel")')

# 27. Upregulated and Downregulated Genes - counts and top-gene heatmap
up_genes <- gbm_de$gene[!is.na(gbm_de$FDR) & gbm_de$FDR<0.05 & gbm_de$logFC>0.25]
down_genes <- gbm_de$gene[!is.na(gbm_de$FDR) & gbm_de$FDR<0.05 & gbm_de$logFC< -0.25]
write.csv(data.frame(gene=up_genes),'GBM_upregulated.csv',row.names=FALSE)
write.csv(data.frame(gene=down_genes),'GBM_downregulated.csv',row.names=FALSE)
p27a <- ggplot(data.frame(direction=c('Upregulated','Downregulated'), genes=c(length(up_genes),length(down_genes))),
               aes(x=direction,y=genes,fill=direction)) + geom_col() + theme_classic() +
  labs(title='Number of DE genes (exploratory)',x='',y='Gene count')
print(p27a)
ggsave('GBM_plots/27_DE_gene_counts.png',p27a,width=7,height=5,dpi=300)
# Heatmap uses normalized pseudobulk logCPM, not per-cell expression.
top_genes <- head(gbm_de$gene[order(gbm_de$PValue)], 20)
top_genes <- intersect(top_genes, rownames(dge))
if (length(top_genes)>1) {
  logcpm <- cpm(dge, log=TRUE, prior.count=2)
  z <- t(scale(t(logcpm[top_genes,,drop=FALSE])))
  z[!is.finite(z)] <- 0
  png('GBM_plots/27_Top_DE_gene_heatmap.png',width=1600,height=1400,res=200)
  heatmap(z, Rowv=NA, Colv=NA, scale='none', margins=c(10,9), main='Top genes (exploratory pseudobulk)')
  dev.off()
}
# 28.0. Load Required Libraries

library(CellChat)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Create Output Folder
dir.create("GBM_plots", showWarnings = FALSE)

# 28.1. Initialization & Database Setup

# Rename Clusters
merged_gbm$cellchat_group <- paste0(
  "Cluster_",
  as.character(merged_gbm$seurat_clusters)
)

# Add Sample Information
merged_gbm$samples <- as.character(merged_gbm$sample_id)

# Check Groups
table(merged_gbm$cellchat_group)

# Create CellChat Object
cellchat <- createCellChat(
  object = merged_gbm,
  group.by = "cellchat_group",
  assay = "RNA"
)

# Use Human CellChat Database
cellchat@DB <- CellChatDB.human

# Subset Expression Data
cellchat <- subsetData(cellchat)

# 28.2. Identify Overexpressed Genes & Interactions

# Identify Signaling Genes
cellchat <- identifyOverExpressedGenes(
  cellchat,
  do.fast = FALSE
)

# Identify Ligand-Receptor Interactions
cellchat <- identifyOverExpressedInteractions(cellchat)

# 28.3. Infer Cell-Cell Communication

# Compute Communication Probabilities
cellchat <- computeCommunProb(cellchat)

# Filter Small Cell Groups
cellchat <- filterCommunication(
  cellchat,
  min.cells = 10)

# 28.4. Build the Communication Network

# Compute Signaling Pathway Probabilities
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate Communication Network
cellchat <- aggregateNet(cellchat)


# 28.5.1. Circle Plot - Overall Communication Network

# Cell Group Sizes
groupSize <- as.numeric(table(cellchat@idents))

# Display Circle Plot
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE
)


# Save Circle Plot

png(
  "GBM_plots/28_CellChat_Communication_Network.png",
  width = 1600,
  height = 1300,
  res = 200
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE
)

dev.off()


# 28.5.2. Extract Ligand-Receptor Results

lr_results <- subsetCommunication(cellchat)

# 28.5.3. View Results Table

head(lr_results)

View(lr_results)

# Save Results
write.csv(
  lr_results,
  "GBM_CellChat_Ligand_Receptor_Results.csv",
  row.names = FALSE
)


# 28.6.1. First Detected Signaling Pathway

pathways <- cellchat@netP$pathways

if (length(pathways) >= 1) {
  
  pathway1 <- pathways[1]
  
  print(pathway1)
  
  netVisual_aggregate(
    cellchat,
    signaling = pathway1,
    layout = "circle"
  )
  
}

# 28.6.2. Second Detected Signaling Pathway

if (length(pathways) >= 2) {
  
  pathway2 <- pathways[2]
  
  print(pathway2)
  
  netVisual_aggregate(
    cellchat,
    signaling = pathway2,
    layout = "circle"
  )
  
}
# 28.7. Save CellChat Object

saveRDS(
  cellchat,
  "GBM_CellChat_Results.rds"
)
cellchat <- readRDS("GBM_CellChat_Results.rds")

# 29-34. GO and KEGG enrichment - dotplots and barplots
# Use all tested genes as the background universe for GO analysis.
all_ids <- bitr(gbm_de$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db)
background_ids <- unique(all_ids$ENTREZID)

# 29. GO Biological Process - Upregulated
if (length(up_genes)>0) {
  up_ids <- bitr(up_genes,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db)
  if (nrow(up_ids)>0) {
    go_up_bp <- enrichGO(gene=unique(up_ids$ENTREZID), universe=background_ids, OrgDb=org.Hs.eg.db,ont='BP',readable=TRUE)
    write.csv(as.data.frame(go_up_bp),'GBM_GO_BP_up.csv',row.names=FALSE)
    if(nrow(as.data.frame(go_up_bp))>0) {
      p29a <- dotplot(go_up_bp,showCategory=15)+ggtitle('GO BP: Upregulated')
      p29b <- barplot(go_up_bp,showCategory=15)+ggtitle('GO BP: Upregulated')
      print(p29a);print(p29b)
      ggsave('GBM_plots/29_GO_BP_up_dotplot.png',p29a,width=11,height=7,dpi=300)
      ggsave('GBM_plots/29_GO_BP_up_barplot.png',p29b,width=11,height=7,dpi=300)
    }
  }
}

# 30. GO Biological Process - Downregulated
if (length(down_genes)>0) {
  down_ids <- bitr(down_genes,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db)
  if (nrow(down_ids)>0) {
    go_down_bp <- enrichGO(gene=unique(down_ids$ENTREZID),universe=background_ids,OrgDb=org.Hs.eg.db,ont='BP',readable=TRUE)
    write.csv(as.data.frame(go_down_bp),'GBM_GO_BP_down.csv',row.names=FALSE)
    if(nrow(as.data.frame(go_down_bp))>0) {
      p30 <- dotplot(go_down_bp,showCategory=15)+ggtitle('GO BP: Downregulated')
      print(p30)
      ggsave('GBM_plots/30_GO_BP_down_dotplot.png',p30,width=11,height=7,dpi=300)
    }
  }
}

# 31. GO Molecular Function - Upregulated
if (exists('up_ids') && nrow(up_ids)>0) {
  go_up_mf <- enrichGO(gene=unique(up_ids$ENTREZID),universe=background_ids,OrgDb=org.Hs.eg.db,ont='MF',readable=TRUE)
  write.csv(as.data.frame(go_up_mf),'GBM_GO_MF_up.csv',row.names=FALSE)
  if(nrow(as.data.frame(go_up_mf))>0) {
    p31 <- dotplot(go_up_mf,showCategory=15)+ggtitle('GO MF: Upregulated')
    print(p31)
    ggsave('GBM_plots/31_GO_MF_up_dotplot.png',p31,width=11,height=7,dpi=300)
  }
}

# 32. GO Cellular Component - Upregulated
if (exists('up_ids') && nrow(up_ids)>0) {
  go_up_cc <- enrichGO(gene=unique(up_ids$ENTREZID),universe=background_ids,OrgDb=org.Hs.eg.db,ont='CC',readable=TRUE)
  write.csv(as.data.frame(go_up_cc),'GBM_GO_CC_up.csv',row.names=FALSE)
  if(nrow(as.data.frame(go_up_cc))>0) {
    p32 <- dotplot(go_up_cc,showCategory=15)+ggtitle('GO CC: Upregulated')
    print(p32)
    ggsave('GBM_plots/32_GO_CC_up_dotplot.png',p32,width=11,height=7,dpi=300)
  }
}

# 33. KEGG - Upregulated
if (exists('up_ids') && nrow(up_ids)>0) {
  kegg_up <- enrichKEGG(gene=unique(up_ids$ENTREZID),universe=background_ids,organism='hsa')
  write.csv(as.data.frame(kegg_up),'GBM_KEGG_up.csv',row.names=FALSE)
  if(nrow(as.data.frame(kegg_up))>0) {
    p33 <- dotplot(kegg_up,showCategory=15)+ggtitle('KEGG: Upregulated')
    print(p33)
    ggsave('GBM_plots/33_KEGG_up_dotplot.png',p33,width=11,height=7,dpi=300)
  }
}

# 34. KEGG - Downregulated
if (exists('down_ids') && nrow(down_ids)>0) {
  kegg_down <- enrichKEGG(gene=unique(down_ids$ENTREZID),universe=background_ids,organism='hsa')
  write.csv(as.data.frame(kegg_down),'GBM_KEGG_down.csv',row.names=FALSE)
  if(nrow(as.data.frame(kegg_down))>0) {
    p34 <- dotplot(kegg_down,showCategory=15)+ggtitle('KEGG: Downregulated')
    print(p34)
    ggsave('GBM_plots/34_KEGG_down_dotplot.png',p34,width=11,height=7,dpi=300)
  }
}
# Create Cell Annotations

gbm_annotations <- data.frame(
  cell = colnames(merged_gbm),
  group = paste(
    merged_gbm$patient_id,
    merged_gbm$condition,
    sep = "_"
  )
)


# 35. GBM Cell State Scoring - exploratory gene programs
# Curated small marker sets are proxies, not validated Neftel signatures.
gbm_state_sets <- list(
  OPC_like = c('PDGFRA','OLIG1','OLIG2','SOX10','CSPG4'),
  NPC_like = c('SOX2','NES','DLL3','ASCL1','DCX'),
  AC_like  = c('GFAP','SLC1A3','ALDOC','AQP4','S100B'),
  MES_like = c('CD44','CHI3L1','VIM','SERPINE1','ANXA1')
)

for (state in names(gbm_state_sets)) {
  
  genes <- intersect(
    gbm_state_sets[[state]],
    rownames(merged_gbm)
  )
  
  if (length(genes) >= 3) {


    merged_gbm <- AddModuleScore(
      merged_gbm,
      features = list(genes),
      name = paste0(state, "_score")
    )
    
    col <- paste0(state, "_score1")


    p <- VlnPlot(
      merged_gbm,
      features = col,
      group.by = "condition",
      pt.size = 0
    ) +
      ggtitle(paste(state, "by condition (descriptive)"))
    
    print(p)
    
    ggsave(
      paste0("GBM_plots/35_", state, "_Violin.png"),
      p,
      width = 7,
      height = 5,
      dpi = 300
    )


    # Compare sample-level averages, not cells as independent biological replicates.
    sample_scores <- aggregate(merged_gbm@meta.data[[col]],
                               by=list(sample_id=merged_gbm$sample_id,patient_id=merged_gbm$patient_id,condition=merged_gbm$condition),FUN=mean)
    names(sample_scores)[4] <- 'mean_score'
    p <- ggplot(sample_scores,aes(x=condition,y=mean_score,group=patient_id,color=patient_id)) +
      geom_line() + geom_point(size=3) + theme_classic() + labs(title=paste(state,'paired patient means'),y='Mean module score')
    print(p)
    ggsave(paste0('GBM_plots/35_',state,'_Paired.png'),p,width=7,height=5,dpi=300)
    write.csv(sample_scores,paste0('GBM_',state,'_sample_scores.csv'),row.names=FALSE)
  } else message('Skipping ',state,': fewer than 3 marker genes found')
}

# 36. Radiation Response Signatures - expression patterns
radiation_genes <- intersect(c('CDKN1A','GADD45A','DDB2','MDM2','BAX','RAD51','H2AX','H2AFX','ATM','ATR','CHEK1','CHEK2','BBC3'),rownames(merged_gbm))
if (length(radiation_genes)>0) {
  p <- DotPlot(merged_gbm,features=radiation_genes,group.by='sample_id') + RotatedAxis() + ggtitle('Radiation response genes by sample')
  print(p)
  ggsave('GBM_plots/36_Radiation_DotPlot.png',p,width=12,height=6,dpi=300)
  scaled_rad <- intersect(radiation_genes,rownames(LayerData(merged_gbm,assay='RNA',layer='scale.data')))
  if(length(scaled_rad)>0) {
    p <- DoHeatmap(merged_gbm,features=scaled_rad,group.by='sample_id') + ggtitle('Radiation genes (scaled expression)')
    print(p)
    ggsave('GBM_plots/36_Radiation_Heatmap.png',p,width=12,height=7,dpi=300)
  }
  # Raw pseudobulk count summary is descriptive; edgeR results above are inferential.
  rad_pb <- as.matrix(pb_counts[intersect(radiation_genes,rownames(pb_counts)),,drop=FALSE])
  write.csv(rad_pb,'GBM_radiation_response_pseudobulk_counts.csv')
}

# 37. Cell Cycle Scoring - G1/S/G2M descriptive analysis
cc_s <- intersect(Seurat::cc.genes.updated.2019$s.genes,rownames(merged_gbm))
cc_g2m <- intersect(Seurat::cc.genes.updated.2019$g2m.genes,rownames(merged_gbm))
if(length(cc_s)>=5 && length(cc_g2m)>=5) {
  merged_gbm <- CellCycleScoring(merged_gbm,s.features=cc_s,g2m.features=cc_g2m,set.ident=FALSE)
  p <- DimPlot(merged_gbm,group.by='Phase',reduction='umap') + ggtitle('Cell cycle phase')
  print(p)
  ggsave('GBM_plots/37_CellCycle_UMAP.png',p,width=8,height=6,dpi=300)
  phase_counts <- as.data.frame(table(sample_id=merged_gbm$sample_id,Phase=merged_gbm$Phase))
  p <- ggplot(phase_counts,aes(x=sample_id,y=Freq,fill=Phase)) + geom_col(position='fill') +
    scale_y_continuous(labels=scales::percent) + theme_classic() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) + labs(y='Fraction of cells',title='Cell cycle phase by sample')
  print(p)
  ggsave('GBM_plots/37_CellCycle_Composition.png',p,width=9,height=6,dpi=300)
  write.csv(phase_counts,'GBM_cell_cycle_counts.csv',row.names=FALSE)
} else message('Cell cycle gene lists not available or insufficient matches; check Seurat version and gene symbols')

# 38. GSEA - ranked exploratory pseudobulk results
# Requires clusterProfiler, org.Hs.eg.db; ranking uses edgeR logFC from the paired analysis.
ranked <- gbm_de[is.finite(gbm_de$logFC) & !is.na(gbm_de$gene),c('gene','logFC')]
ranked <- ranked[!duplicated(ranked$gene),]
rank_ids <- suppressMessages(bitr(ranked$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db))
ranked <- merge(ranked,rank_ids,by.x='gene',by.y='SYMBOL')
ranked <- ranked[!duplicated(ranked$ENTREZID),]
gene_rank <- ranked$logFC
names(gene_rank) <- ranked$ENTREZID
gene_rank <- sort(gene_rank,decreasing=TRUE)
if(length(gene_rank)>=100) {
  gsea_go <- tryCatch(gseGO(geneList=gene_rank,OrgDb=org.Hs.eg.db,ont='BP',minGSSize=10,maxGSSize=500,
                            pvalueCutoff=1,verbose=FALSE),error=function(e) {message('GSEA skipped: ',conditionMessage(e));NULL})
  if(!is.null(gsea_go)) {
    write.csv(as.data.frame(gsea_go),'GBM_GSEA_GO_BP.csv',row.names=FALSE)
    if(nrow(as.data.frame(gsea_go))>0) {
      p <- dotplot(gsea_go,showCategory=15) + ggtitle('GSEA GO Biological Process (exploratory)')
      print(p)
      ggsave('GBM_plots/38_GSEA_GO_DotPlot.png',p,width=12,height=8,dpi=300)
      p <- enrichplot::gseaplot2(gsea_go,geneSetID=1,title=gsea_go@result$Description[1])
      print(p)
      ggsave('GBM_plots/38_GSEA_Top_Pathway.png',p,width=10,height=7,dpi=300)
    }
  }
} else message('GSEA skipped: insufficient ranked genes')

# 39. Transcription Factor Expression (not direct activity inference)
# These plots show TF RNA abundance; activity inference requires a validated regulon model.
tf_genes <- intersect(c('STAT3','RELA','NFKB1','MYC','E2F1','TP53','SOX2','OLIG2'),rownames(merged_gbm))
if(length(tf_genes)>0) {
  p <- DotPlot(merged_gbm,features=tf_genes,group.by='sample_id') + RotatedAxis() +
    ggtitle('Transcription factor RNA expression (not activity)')
  print(p)
  ggsave('GBM_plots/39_TF_Expression_DotPlot.png',p,width=11,height=6,dpi=300)
}

# 40. Pseudotime - OPTIONAL exploratory Slingshot analysis
# Pseudotime does not establish that irradiation caused cell-state transitions.
# Requires: BiocManager::install('slingshot')
if (requireNamespace('slingshot',quietly=TRUE)) {
  sce_trajectory <- as.SingleCellExperiment(merged_gbm)
  # Only run if UMAP and multiple clusters are available.
  if ('UMAP' %in% reducedDimNames(sce_trajectory) && length(unique(merged_gbm$seurat_clusters))>=2) {
    sce_trajectory <- tryCatch(slingshot::slingshot(sce_trajectory,clusterLabels='seurat_clusters',reducedDim='UMAP'),
                               error=function(e) {message('Slingshot skipped: ',conditionMessage(e));NULL})
    if(!is.null(sce_trajectory)) {
      pseudo <- slingshot::slingPseudotime(sce_trajectory)
      if(ncol(pseudo)>0) {
        merged_gbm$pseudotime_1 <- pseudo[,1]
        p <- FeaturePlot(merged_gbm,features='pseudotime_1',reduction='umap') + ggtitle('Exploratory pseudotime')
        print(p)
        ggsave('GBM_plots/40_Pseudotime_UMAP.png',p,width=8,height=6,dpi=300)
      }
    }
  }
} else message('Optional Slingshot analysis skipped: package not installed')

# 41. Save Final Seurat Object
saveRDS(merged_gbm,'GBM_integrated_final.rds')
cat('Completed. Plots are in:',normalizePath('GBM_plots'),'\n')

# Load packages
library(Seurat)
library(SingleR)
library(SeuratWrappers)
library(monocle3)
library(reticulate)
library(Matrix)
library(ggplot2)
library(cowplot)
library(magrittr)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 300000 * 1024^2)

# Load and process expression data
Control.data <- Read10X(data.dir = "Control\\filtered_feature_bc_matrix/")
MI_1D.data <- Read10X(data.dir = "MI_1D\\filtered_feature_bc_matrix/")
MI_3D.data <- Read10X(data.dir = "MI_3D\\filtered_feature_bc_matrix/")
MI_5D.data <- Read10X(data.dir = "MI_5D\\filtered_feature_bc_matrix/")
MI_7D.data <- Read10X(data.dir = "MI_7D\\filtered_feature_bc_matrix/")

MI.control <- CreateSeuratObject(Control.data, project = "MI_CTR", min.cells = 10)
MI.1D <- CreateSeuratObject(MI_1D.data, project = "MI_1D", min.cells = 10)
MI.3D <- CreateSeuratObject(MI_3D.data, project = "MI_3D", min.cells = 10)
MI.5D <- CreateSeuratObject(MI_5D.data, project = "MI_5D", min.cells = 10)
MI.7D <- CreateSeuratObject(MI_7D.data, project = "MI_7D", min.cells = 10)

MI.control[['percent.mito']] <- PercentageFeatureSet(MI.control, pattern = "^mt-")
MI.1D[['percent.mito']] <- PercentageFeatureSet(MI.1D, pattern = "^mt-")
MI.3D[['percent.mito']] <- PercentageFeatureSet(MI.3D, pattern = "^mt-")
MI.5D[['percent.mito']] <- PercentageFeatureSet(MI.5D, pattern = "^mt-")
MI.7D[['percent.mito']] <- PercentageFeatureSet(MI.7D, pattern = "^mt-")

FeatureScatter(MI.control, "nCount_RNA", "nFeature_RNA")
FeatureScatter(MI.control, "nCount_RNA", "percent.mito")
FeatureScatter(MI.1D, "nCount_RNA", "nFeature_RNA")
FeatureScatter(MI.1D, "nCount_RNA", "percent.mito")
FeatureScatter(MI.3D, "nCount_RNA", "nFeature_RNA")
FeatureScatter(MI.3D, "nCount_RNA", "percent.mito")
FeatureScatter(MI.5D, "nCount_RNA", "nFeature_RNA")
FeatureScatter(MI.5D, "nCount_RNA", "percent.mito")
FeatureScatter(MI.7D, "nCount_RNA", "nFeature_RNA")
FeatureScatter(MI.7D, "nCount_RNA", "percent.mito")

VlnPlot(MI.control, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(MI.1D, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(MI.3D, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(MI.5D, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(MI.7D, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# Apply filter criteria
quantile(MI.control$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.1D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.3D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.5D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.7D$nFeature_RNA, probs = c(0.02,0.98))

MI.control <- subset(x = MI.control, subset = nFeature_RNA > 371 & nFeature_RNA < 5061 & percent.mito < 10)
MI.1D <- subset(x = MI.1D, subset = nFeature_RNA > 350 & nFeature_RNA < 5199 & percent.mito < 10)
MI.3D <- subset(x = MI.3D, subset = nFeature_RNA > 407 & nFeature_RNA < 7825 & percent.mito < 10)
MI.5D <- subset(x = MI.5D, subset = nFeature_RNA > 493 & nFeature_RNA < 6741 & percent.mito < 10)
MI.7D <- subset(x = MI.7D, subset = nFeature_RNA > 410 & nFeature_RNA < 6132 & percent.mito < 10)

# MergeData
experiment.aggregate <- merge(x = MI.control, y = c(MI.1D,MI.3D,MI.5D,MI.7D))

samplename = experiment.aggregate@meta.data$orig.ident
table(samplename)

batchorder = rep("MI",length(samplename))
batchorder[samplename %in% c("MI_CTR")] = "1"
batchorder[samplename %in% c("MI_1D")] = "2"
batchorder[samplename %in% c("MI_3D")] = "3"
batchorder[samplename %in% c("MI_5D")] = "4"
batchorder[samplename %in% c("MI_7D")] = "5"
names(batchorder) = colnames(experiment.aggregate)

batchid = rep("MI",length(samplename))
batchid[samplename %in% c("MI_CTR")] = "MI_CTR"
batchid[samplename %in% c("MI_1D")] = "MI_1D"
batchid[samplename %in% c("MI_3D")] = "MI_3D"
batchid[samplename %in% c("MI_5D")] = "MI_5D"
batchid[samplename %in% c("MI_7D")] = "MI_7D"
names(batchid) = colnames(experiment.aggregate)

experiment.aggregate <- AddMetaData(object = experiment.aggregate, metadata = batchorder, col.name = "batchorder")
table(experiment.aggregate@meta.data$batchorder)

experiment.aggregate <- AddMetaData(object = experiment.aggregate, metadata = batchid, col.name = "batchid")
table(experiment.aggregate@meta.data$batchid)

# Perform integration
MI.list <- SplitObject(experiment.aggregate, split.by = "batchid")

for (i in 1:length(MI.list)) {
  MI.list[[i]] <- NormalizeData(MI.list[[i]], verbose = FALSE)
  MI.list[[i]] <- FindVariableFeatures(MI.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

MI.anchors <- FindIntegrationAnchors(object.list = MI.list, dims = 1:30)
MI.integrated <- IntegrateData(anchorset = MI.anchors, dims = 1:30)
all.genes <- rownames(MI.integrated)
MI.integrated <- ScaleData(MI.integrated, features = all.genes, verbose = FALSE)
MI.integrated <- RunPCA(MI.integrated, verbose = FALSE)
ElbowPlot(object = MI.integrated, ndims = 50)
use.pcs = 1:40
MI.integrated <- RunUMAP(MI.integrated, reduction = "pca", dims = use.pcs)
MI.integrated <- FindNeighbors(object = MI.integrated, reduction = "pca", dims = use.pcs)
MI.integrated <- FindClusters(object = MI.integrated, resolution = seq(0.5,2,0.1))
sapply(grep("^integrated_snn_res",colnames(MI.integrated@meta.data),value = TRUE), function(x) length(unique(MI.integrated@meta.data[,x])))
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "integrated_snn_res.1.4", pt.size=0.5, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "batchid", pt.size=0.5)
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "orig.ident", split.by = "batchid",  pt.size=0.5, ncol = 3) + NoLegend()
MI.integrated <- SetIdent(MI.integrated, value = 'integrated_snn_res.1.4')

DefaultAssay(MI.integrated) <- "RNA"
MI.integrated <- NormalizeData(MI.integrated, verbose = FALSE)
all.genes <- rownames(MI.integrated)
MI.integrated <- ScaleData(MI.integrated, features = all.genes, verbose = FALSE)
markers_all_RNA <- FindAllMarkers(object = MI.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_preprocessing_cluster_allmarker_Genes_RNA_RES14.csv")
table(MI.integrated@meta.data$integrated_snn_res.1.4,MI.integrated@meta.data$batchid)

# SingleR
mrsd.se <- ImmGenData()
MI.integrated.se <- as.SingleCellExperiment(MI.integrated, assay = "RNA")
mrsd.common <- intersect(rownames(MI.integrated.se), rownames(mrsd.se))
mrsd.se <- mrsd.se[mrsd.common,]
MI.integrated.se <- MI.integrated.se[mrsd.common,]
MI.integrated.se <- logNormCounts(MI.integrated.se)
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.main)
MI.integrated[["mrsd.main"]] <- DC.mrsd.pred$labels
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.fine)
MI.integrated[["mrsd.fine"]] <- DC.mrsd.pred$labels

# Filter Endothelial/Fibroblast, double cell type clusteres and other technical noise[PMID:30471926] (CL 5, 10, 12, 23, 24, 27, 29, 30) 
MI.integrated.1st.filtered = subset(MI.integrated, idents = c(0,1,2,3,4,6,7,8,9,11,13,14,15,16,17,18,19,20,21,22,25,26,28,31,32,33))
DefaultAssay(MI.integrated.1st.filtered) <- "integrated"
MI.integrated.1st.filtered <- RunPCA(MI.integrated.1st.filtered, verbose = FALSE)
ElbowPlot(object = MI.integrated.1st.filtered, ndims = 50)
use.pcs = 1:40
MI.integrated.1st.filtered <- RunUMAP(MI.integrated.1st.filtered, reduction = "pca", dims = use.pcs)
MI.integrated.1st.filtered <- FindNeighbors(object = MI.integrated.1st.filtered, reduction = "pca", dims = use.pcs)
MI.integrated.1st.filtered <- FindClusters(object = MI.integrated.1st.filtered, resolution = seq(0.5,2,0.1))
sapply(grep("^integrated_snn_res",colnames(MI.integrated.1st.filtered@meta.data),value = TRUE), function(x) length(unique(MI.integrated.1st.filtered@meta.data[,x])))
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "integrated_snn_res.1.5", pt.size=1, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "batchid", pt.size=0.5)
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "orig.ident", split.by = "batchid",  pt.size=0.5, ncol = 3) + NoLegend()
MI.integrated.1st.filtered <- SetIdent(MI.integrated.1st.filtered, value = 'integrated_snn_res.1.5')

# Cluster rename
new.cluster.ids <- c("MAC", "MAC", "MAC", "MAC", "Neutrophil", "Neutrophil", "MAC", "MAC", "B", "MAC", "MAC", "Monocytes", "MAC", "NK", "MAC", "DC", "Neutrophil", "MAC", "T", "DC", "Monocytes", "Xcr1_DC", "Mature_DC", "DC", "ILC2", "MAC", "Plasma", "Mast")
names(x = new.cluster.ids) <- levels(x = MI.integrated.1st.filtered)
MI.integrated.1st.filtered <- RenameIdents(object = MI.integrated.1st.filtered, new.cluster.ids)
MI.integrated.1st.filtered <- StashIdent(MI.integrated.1st.filtered, save.name = "CellType")
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "CellType", pt.size=0.5)
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "CellType", split.by = "batchid",  pt.size=0.5, ncol = 2) + NoLegend()
markers_all_RNA <- FindAllMarkers(object = MI.integrated.1st.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_1st_filtered_cluster_allmarker_Genes_RNA_RES15.csv")
table(MI.integrated.1st.filtered@meta.data$integrated_snn_res.1.5,MI.integrated.1st.filtered@meta.data$batchid)

# Main_Figure1
DimPlot(object = MI.integrated.1st.filtered, reduction = 'umap', group.by = "CellType", pt.size=1) 
MI.integrated.1st.filtered <- SetIdent(MI.integrated.1st.filtered, value = 'batchid')
fig1bSS = subset(MI.integrated.1st.filtered, idents = "MI_7D")
DimPlot(object = fig1bSS, reduction = 'umap', group.by = "CellType", pt.size=1) + NoLegend()+ xlim(c(-13.5,13.5)) + ylim(c(-13,15))
fig1bD1 = subset(MI.integrated.1st.filtered, idents = "MI_7D")
DimPlot(object = fig1bD1, reduction = 'umap', group.by = "CellType", pt.size=1) + NoLegend()+ xlim(c(-13.5,13.5)) + ylim(c(-13,15))
fig1bD3 = subset(MI.integrated.1st.filtered, idents = "MI_7D")
DimPlot(object = fig1bD3, reduction = 'umap', group.by = "CellType", pt.size=1) + NoLegend()+ xlim(c(-13.5,13.5)) + ylim(c(-13,15))
fig1bD5 = subset(MI.integrated.1st.filtered, idents = "MI_7D")
DimPlot(object = fig1bD5, reduction = 'umap', group.by = "CellType", pt.size=1) + NoLegend()+ xlim(c(-13.5,13.5)) + ylim(c(-13,15))
fig1bD7 = subset(MI.integrated.1st.filtered, idents = "MI_7D")
DimPlot(object = fig1bD7, reduction = 'umap', group.by = "CellType", pt.size=1) + NoLegend()+ xlim(c(-13.5,13.5)) + ylim(c(-13,15))

Macrophage_genes <- c("Adgre1", "Cd68", "Csf1r")
Neutrophil_genes <- c("Csf3r", "S100a9", "S100a8")
B_genes <- c("Ms4a1", "Cd79a", "Ly6d")
Monocyte_genes <- c("Ly6c2", "Chil3", "F10")
NK_genes <- c("Nkg7", "Klrb1c", "Gzma")
CD209DC_genes <- c("Cd209a", "Klrd1", "Flt3")
T_genes <- c("Cd3e", "Cd3d", "Lef1")
XcrDC_genes <- c("Xcr1", "Ifi205", "Itgae")
MigratoryDC_genes <- c("Ccr7", "Fscn1", "Cacnb3")
ILC2_genes <- c("Rora", "Cxcr6", "Gata3")
PC_genes <- c("Jchain", "Iglv1", "Mzb1")
Mast_genes <- c("Cma1", "Kit", "Rab27b")

MI.integrated.1st.filtered <- SetIdent(MI.integrated.1st.filtered, value = 'CellType')
MI.integrated.1st.filtered@active.ident <- factor(MI.integrated.1st.filtered@active.ident, levels = c("Mast","Plasma","ILC2","Mature_DC","Xcr1_DC","T","DC","NK","Monocytes","B","Neutrophil","MAC"))
DotPlot(object = MI.integrated.1st.filtered, features = c(Macrophage_genes,Neutrophil_genes,B_genes,Monocyte_genes,NK_genes,CD209DC_genes,T_genes,XcrDC_genes,MigratoryDC_genes,ILC2_genes,PC_genes,Mast_genes), cols = c("White", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
FeaturePlot(MI.integrated.1st.filtered, features = "Adgre1",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "S100a9",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Cd79a",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Ly6c2",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Nkg7",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Cd209a",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Cd3e",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Xcr1",cols = c("gray", "red"), pt.size=1)
FeaturePlot(MI.integrated.1st.filtered, features = "Fscn1",cols = c("gray", "red"), pt.size=1)

LR_markers_MI_filtered_Final <- FindAllMarkers(object = MI.integrated.1st.filtered, test.use = "LR", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(LR_markers_MI_filtered_Final, file = "Cluster_allmarker_Genes_Integrated_Filtered_Final_12_cell_Types.csv")
table(MI.integrated.1st.filtered@meta.data$integrated_snn_res.0.9,MI.integrated.1st.filtered@meta.data$orig.ident)


# Macrophage subset analysis
MI.macrophage = subset(MI.integrated.1st.filtered, idents = c(0,1,2,3,6,7,9,10,11,12,14,17,20,25))
DefaultAssay(MI.macrophage) <- "integrated"
MI.macrophage.recluster <- RunPCA(object = MI.macrophage, verbose = FALSE)
ElbowPlot(object = MI.macrophage.recluster, ndims = 50)
use.pcs = 1:30
MI.macrophage.recluster <- FindNeighbors(object = MI.macrophage.recluster, reduction = "pca", dims = use.pcs)
MI.macrophage.recluster <- FindClusters(object = MI.macrophage.recluster, reduction.type = "pca", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
MI.macrophage.recluster <- RunUMAP(object = MI.macrophage.recluster, reduction = "pca", dims = use.pcs)
sapply(grep("^integrated_snn_res",colnames(MI.macrophage.recluster@meta.data),value = TRUE), function(x) length(unique(MI.macrophage.recluster@meta.data[,x])))
DimPlot(object = MI.macrophage.recluster, reduction = 'umap', group.by = "integrated_snn_res.1.5", pt.size=0.1, label = TRUE, repel = TRUE) + NoLegend()
MI.macrophage.recluster <- SetIdent(MI.macrophage.recluster, value = 'integrated_snn_res.1.5')

# 1st round filter: filtering clusters other than macrophages
MI.macrophage.recluster = subset(MI.macrophage.recluster, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18))
DefaultAssay(MI.macrophage.recluster) <- "integrated"
MI.macrophage.recluster <- RunPCA(object = MI.macrophage.recluster, verbose = FALSE)
ElbowPlot(object = MI.macrophage.recluster, ndims = 50)
use.pcs = 1:40
MI.macrophage.recluster <- FindNeighbors(object = MI.macrophage.recluster, reduction = "pca", dims = use.pcs)
MI.macrophage.recluster <- FindClusters(object = MI.macrophage.recluster, reduction.type = "pca", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
MI.macrophage.recluster <- RunUMAP(object = MI.macrophage.recluster, reduction = "pca", dims = use.pcs)
sapply(grep("^integrated_snn_res",colnames(MI.macrophage.recluster@meta.data),value = TRUE), function(x) length(unique(MI.macrophage.recluster@meta.data[,x])))
DimPlot(object = MI.macrophage.recluster, reduction = 'umap', group.by = "integrated_snn_res.1.5", pt.size=0.1, label = TRUE, repel = TRUE) + NoLegend()
MI.macrophage.recluster <- SetIdent(MI.macrophage.recluster, value = 'integrated_snn_res.1.5')

# 2nd round filter and Figure 3: neutrophil subcluster filtering
MI.macrophage.recluster = subset(MI.macrophage.recluster, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16))
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19)
new.cluster.ids <- c(1,0,3,4,2,5,6,7,0,9,2,10,11,12,8,13,8,14,15)
MI.macrophage.recluster@active.ident <- plyr::mapvalues(x = MI.macrophage.recluster@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MI.macrophage.recluster <- StashIdent(MI.macrophage.recluster, save.name = "FigureCluster")
DimPlot(object = MI.macrophage.recluster, reduction = 'umap', group.by = "FigureCluster", pt.size=1, label = TRUE, repel = TRUE)
MI.macrophage.recluster <- SetIdent(MI.macrophage.recluster, value = 'FigureCluster')
MI.macrophage.recluster@active.ident <- factor(MI.macrophage.recluster@active.ident, ordered = T, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

DefaultAssay(MI.macrophage.recluster) <- "RNA"
markers_all_RNA <- FindAllMarkers(object = MI.macrophage.recluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_MF_subset_Figure_Genes.csv")
table(MI.macrophage.recluster@meta.data$FigureCluster,MI.macrophage.recluster@meta.data$batchid)

# Figure 3a right
top10 <- markers_all_RNA %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MI.macrophage.recluster, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

# Supplementary Figure 2
VlnPlot(MI.macrophage.recluster, features = c("Lyve1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Cd163"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Timd4"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("H2-Eb1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("H2-Ab1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Cd74"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Ccr2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Ly6c2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Plac8"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Cd68"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Fcgr1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Itgam"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Ace"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Ear2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Itgal"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Fcrls"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Rgs10"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Adgre1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Trem2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Gpnmb"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.macrophage.recluster, features = c("Spp1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")

#cell cycle analysis and Supplementary Figure 4
s.genes <- c("Usp1", "Casp8ap2", "Brip1", "Prim1", "Ung", "Atad2", "Gmnn", "Ubr7", "Rpa2", "Mcm2", "Clspn", "Pola1", "Fen1", "Dscc1", "Uhrf1", "Chaf1b", "Dtl", "Blm", "Pold3", "Rrm2", "Rad51ap1", "Tipin", "Ccne2", "Pcna", "Wdr76", "Msh2", "Nasp", "Rad51", "Hells", "Mcm6", "Mcm4", "Mcm5", "Rfc2", "Slbp", "Tyms", "Cdc45", "E2f8", "Cenpu", "Rrm1", "Cdca7", "Cdc6", "Exo1", "Gins2")
g2m.genes <- c("Ube2c","Lbr","Ctcf","Cdc20","Cbx5","Kif11","Anp32e","Birc5","Cdk1","Anln", "Aurka", "Aurkb", "Bub1", "Ccnb2", "Cdc25c", "Cdca2", "Cdca3", "Cdca8", "Cenpa", "Cenpe", "Cenpf", "Ckap2", "Ckap2l", "Ckap5", "Cks1b", "Cks2", "Dlgap5", "Ect2", "Gas2l3", "Gtse1", "Hjurp", "Hmgb2", "Hmmr", "Kif20b", "Kif23", "Kif2c", "Mki67", "Ncapd2", "Ndc80", "Nek2", "Nuf2", "Nusap1", "Psrc1", "Rangap1", "Smc4", "Tacc3", "Tmpo", "Top2a", "Tpx2", "Ttk", "Tubb4b")
MI.macrophage.recluster <- CellCycleScoring(MI.macrophage.recluster, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(MI.macrophage.recluster[[]])
FeaturePlot(MI.macrophage.recluster, features = "S.Score",cols = c("yellow", "red"), reduction = 'umap', pt.size=1) + DarkTheme()
FeaturePlot(MI.macrophage.recluster, features = "G2M.Score",cols = c("yellow", "red"), reduction = 'umap', pt.size=1) + DarkTheme()

# Pseudotime analysis and Figure 3c-d
cds <- as.cell_data_set(MI.macrophage.recluster, assay = "RNA")
cds <- estimate_size_factors(cds)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(MI.macrophage.recluster[["RNA"]])
plot_cells(cds, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size=7)
plot_cells(cds, color_cells_by = "DEGCluster", label_cell_groups=FALSE, label_leaves=TRUE, label_branch_points=FALSE, graph_label_size=10)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=5, show_trajectory_graph = TRUE, cell_size = 1.5, trajectory_graph_color = "grey", trajectory_graph_segment_size = 2.5)


MF_genes <- c("Hmox1", "Ltc4s", "Saa3", "Clec4e", "Chil3", "Ccr2")

MF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% MF_genes,
                      colData(cds)$DEGCluster %in% c("LateMF","EarlyMF","IFN","IntMF","cMono")]
plot_genes_in_pseudotime(MF_lineage_cds,
                         color_cells_by="DEGCluster",
                         min_expr=0.5, cell_size = 2, panel_order = c("Ccr2","Chil3","Clec4e","Hmox1","Ltc4s","Saa3")) + scale_color_manual(breaks = c("cMono", "EarlyMF", "IntMF", "IFN", "LateMF"), values=c("#ed68ed", "#00bfc4", "#00be68", "#0db802", "#f8766d"))

MF_genes <- c("Mgl2", "Lyve1", "Folr2", "Fcrls", "Rgs10", "Trem2")
MF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% MF_genes,
                       colData(cds)$DEGCluster %in% c("LateMF","EarlyMF","IFN","IntMF","cMono")]
plot_genes_in_pseudotime(MF_lineage_cds,
                         color_cells_by="DEGCluster",
                         min_expr=0.5, cell_size = 2, panel_order = c("Trem2","Rgs10","Fcrls","Folr2","Lyve1","Mgl2")) + scale_color_manual(breaks = c("cMono", "EarlyMF", "IntMF", "IFN", "LateMF"), values=c("#ed68ed", "#00bfc4", "#00be68", "#0db802", "#f8766d"))


# Neutrophil subset analysis (Figure 4a)
MI.neutrophil = subset(MI.integrated.1st.filtered, idents = c(4,5,16))
DefaultAssay(MI.neutrophil) <- "integrated"
MI.neutrophil.recluster <- RunPCA(object = MI.neutrophil, verbose = FALSE)
ElbowPlot(object = MI.neutrophil.recluster, ndims = 50)
use.pcs = 1:15
MI.neutrophil.recluster <- FindNeighbors(object = MI.neutrophil.recluster, reduction = "pca", dims = use.pcs)
MI.neutrophil.recluster <- FindClusters(object = MI.neutrophil.recluster, reduction.type = "pca", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
MI.neutrophil.recluster <- RunUMAP(object = MI.neutrophil.recluster, reduction = "pca", dims = use.pcs)
sapply(grep("^integrated_snn_res",colnames(MI.neutrophil.recluster@meta.data),value = TRUE), function(x) length(unique(MI.neutrophil.recluster@meta.data[,x])))
DimPlot(object = MI.neutrophil.recluster, reduction = 'umap', group.by = "integrated_snn_res.0.5", pt.size=2, label = FALSE, repel = TRUE, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend()
MI.neutrophil.recluster <- SetIdent(MI.neutrophil.recluster, value = 'integrated_snn_res.0.5')

DefaultAssay(MI.neutrophil.recluster) <- "RNA"
markers_all_RNA <- FindAllMarkers(object = MI.neutrophil.recluster, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_Neutrophil_subset_Figure_Genes.csv")
table(MI.neutrophil.recluster@meta.data$integrated_snn_res.0.5,MI.neutrophil.recluster@meta.data$batchid)

# Figure 4a right
all.genes <- rownames(MI.neutrophil.recluster)
MI.neutrophil.recluster <- ScaleData(MI.neutrophil.recluster, features = all.genes, verbose = FALSE)
top10 <- markers_all_RNA %>% group_by(cluster) %>% top_n(n = 10, wt = p_val_adj)
DoHeatmap(MI.neutrophil.recluster, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))
+ NoLegend() 

# Figure 4c
VlnPlot(MI.neutrophil.recluster, features = "Fpr1", pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.neutrophil.recluster, features = "Tnf", pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.neutrophil.recluster, features = "Ifitm2", pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.neutrophil.recluster, features = "Isg15", pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.neutrophil.recluster, features = "Manf", pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300")) + NoLegend() + geom_boxplot(width=0.1, fill="white")


# DC subset analysis (Figure 4d)
MI.DC = subset(MI.integrated.1st.filtered, idents = c(15,19,21,22,23))
DefaultAssay(MI.DC) <- "integrated"
MI.DC.recluster <- RunPCA(object = MI.DC, verbose = FALSE)
ElbowPlot(object = MI.DC.recluster, ndims = 50)
use.pcs = 1:30
MI.DC.recluster <- FindNeighbors(object = MI.DC.recluster, reduction = "pca", dims = use.pcs)
MI.DC.recluster <- FindClusters(object = MI.DC.recluster, reduction.type = "pca", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
MI.DC.recluster <- RunUMAP(object = MI.DC.recluster, reduction = "pca", dims = use.pcs)
sapply(grep("^integrated_snn_res",colnames(MI.DC.recluster@meta.data),value = TRUE), function(x) length(unique(MI.DC.recluster@meta.data[,x])))
DimPlot(object = MI.DC.recluster, reduction = 'umap', group.by = "integrated_snn_res.0.2", pt.size=2, label = FALSE, repel = TRUE, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend()
MI.DC.recluster <- SetIdent(MI.DC.recluster, value = 'integrated_snn_res.0.2')

DefaultAssay(MI.DC.recluster) <- "RNA"
markers_all_RNA <- FindAllMarkers(object = MI.DC.recluster, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_DC_subset_Figure_Genes.csv")
table(MI.DC.recluster@meta.data$integrated_snn_res.0.2,MI.DC.recluster@meta.data$batchid)

# Figure 4d right
all.genes <- rownames(MI.DC.recluster)
MI.DC.recluster <- ScaleData(MI.DC.recluster, features = all.genes, verbose = FALSE)
top10 <- markers_all_RNA %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MI.DC.recluster, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

# Figure 4f
VlnPlot(MI.DC.recluster, features = "Ms4a7", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.DC.recluster, features = "Ccr2", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.DC.recluster, features = "Cd209a", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.DC.recluster, features = "Xcr1", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.DC.recluster, features = "Fscn1", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(MI.DC.recluster, features = "Siglech", pt.size = 0, assay = "RNA", y.max = 4, cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702")) + NoLegend() + geom_boxplot(width=0.1, fill="white")

# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(SPOTlight)

# Load and process expression data
VISIUM_MI_1D <- Load10X_Spatial(data.dir = "Visium\\Raw_Data\\1D/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", filter.matrix = TRUE, to.upper = FALSE)
VISIUM_MI_3D <- Load10X_Spatial(data.dir = "Visium\\Raw_Data\\3D/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", filter.matrix = TRUE, to.upper = FALSE)
VISIUM_MI_5D <- Load10X_Spatial(data.dir = "Visium\\Raw_Data\\5D/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", filter.matrix = TRUE, to.upper = FALSE)
VISIUM_MI_7D <- Load10X_Spatial(data.dir = "Visium\\Raw_Data\\7D/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", filter.matrix = TRUE, to.upper = FALSE)

VISIUM_MI_1D <- SCTransform(VISIUM_MI_1D, assay = "Spatial", verbose = FALSE)
VISIUM_MI_3D <- SCTransform(VISIUM_MI_3D, assay = "Spatial", verbose = FALSE)
VISIUM_MI_5D <- SCTransform(VISIUM_MI_5D, assay = "Spatial", verbose = FALSE)
VISIUM_MI_7D <- SCTransform(VISIUM_MI_7D, assay = "Spatial", verbose = FALSE)

# Data QC
plot1 <- VlnPlot(VISIUM_MI_1D, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_1D, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_1D, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_1D, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_3D, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_3D, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_3D, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_3D, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_5D, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_5D, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_5D, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_5D, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_7D, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_7D, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(VISIUM_MI_7D, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(VISIUM_MI_7D, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Dimensionality reduction, clustering and visualization
VISIUM_MI_1D <- RunPCA(VISIUM_MI_1D, assay = "SCT", verbose = FALSE)
ElbowPlot(object = VISIUM_MI_1D, ndims = 50)
VISIUM_MI_1D <- FindNeighbors(VISIUM_MI_1D, reduction = "pca", dims = 1:15)
VISIUM_MI_1D <- FindClusters(VISIUM_MI_1D, verbose = FALSE, resolution = seq(0.1,1,0.1))
VISIUM_MI_1D <- RunUMAP(VISIUM_MI_1D, reduction = "pca", dims = 1:15)

VISIUM_MI_3D <- RunPCA(VISIUM_MI_3D, assay = "SCT", verbose = FALSE)
ElbowPlot(object = VISIUM_MI_3D, ndims = 50)
VISIUM_MI_3D <- FindNeighbors(VISIUM_MI_3D, reduction = "pca", dims = 1:15)
VISIUM_MI_3D <- FindClusters(VISIUM_MI_3D, verbose = FALSE, resolution = seq(0.1,1,0.1))
VISIUM_MI_3D <- RunUMAP(VISIUM_MI_3D, reduction = "pca", dims = 1:15)

VISIUM_MI_5D <- RunPCA(VISIUM_MI_5D, assay = "SCT", verbose = FALSE)
ElbowPlot(object = VISIUM_MI_5D, ndims = 50)
VISIUM_MI_5D <- FindNeighbors(VISIUM_MI_5D, reduction = "pca", dims = 1:15)
VISIUM_MI_5D <- FindClusters(VISIUM_MI_5D, verbose = FALSE, resolution = seq(0.1,1,0.1))
VISIUM_MI_5D <- RunUMAP(VISIUM_MI_5D, reduction = "pca", dims = 1:15)

VISIUM_MI_7D <- RunPCA(VISIUM_MI_7D, assay = "SCT", verbose = FALSE)
ElbowPlot(object = VISIUM_MI_7D, ndims = 50)
VISIUM_MI_7D <- FindNeighbors(VISIUM_MI_7D, reduction = "pca", dims = 1:15)
VISIUM_MI_7D <- FindClusters(VISIUM_MI_7D, verbose = FALSE, resolution = seq(0.1,1,0.1))
VISIUM_MI_7D <- RunUMAP(VISIUM_MI_7D, reduction = "pca", dims = 1:15)

sapply(grep("^SCT_snn_res",colnames(VISIUM_MI_1D@meta.data),value = TRUE), function(x) length(unique(VISIUM_MI_1D@meta.data[,x])))
sapply(grep("^SCT_snn_res",colnames(VISIUM_MI_3D@meta.data),value = TRUE), function(x) length(unique(VISIUM_MI_3D@meta.data[,x])))
sapply(grep("^SCT_snn_res",colnames(VISIUM_MI_5D@meta.data),value = TRUE), function(x) length(unique(VISIUM_MI_5D@meta.data[,x])))
sapply(grep("^SCT_snn_res",colnames(VISIUM_MI_7D@meta.data),value = TRUE), function(x) length(unique(VISIUM_MI_7D@meta.data[,x])))

SpatialDimPlot(VISIUM_MI_1D, label = TRUE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "SCT_snn_res.0.8")
SpatialDimPlot(VISIUM_MI_3D, label = TRUE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "SCT_snn_res.0.8")
SpatialDimPlot(VISIUM_MI_5D, label = TRUE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "SCT_snn_res.0.5")
SpatialDimPlot(VISIUM_MI_7D, label = TRUE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "SCT_snn_res.0.6")

VISIUM_MI_1D <- SetIdent(VISIUM_MI_1D, value = 'SCT_snn_res.0.8')
VISIUM_MI_3D <- SetIdent(VISIUM_MI_3D, value = 'SCT_snn_res.0.8')
VISIUM_MI_5D <- SetIdent(VISIUM_MI_5D, value = 'SCT_snn_res.0.5')
VISIUM_MI_7D <- SetIdent(VISIUM_MI_7D, value = 'SCT_snn_res.0.6')

# Signature scores (Supplementary Figure 1)
Cardiomyocyte_genes <- list(c("Myh7","Myh6","Actn2","Nkx2-5","Tnni3","Tnnt2"))
Endothelial_genes <- list(c("Cdh5", "Ly6c1", "Kdr"))
Fibro_genes <- list(c("Col1a1", "Pdgfra", "Lamc1"))

VISIUM_MI_1D <- AddModuleScore(object = VISIUM_MI_1D, features = Cardiomyocyte_genes, name = 'Cardiomyocyte')
VISIUM_MI_3D <- AddModuleScore(object = VISIUM_MI_3D, features = Cardiomyocyte_genes, name = 'Cardiomyocyte')
VISIUM_MI_5D <- AddModuleScore(object = VISIUM_MI_5D, features = Cardiomyocyte_genes, name = 'Cardiomyocyte')
VISIUM_MI_7D <- AddModuleScore(object = VISIUM_MI_7D, features = Cardiomyocyte_genes, name = 'Cardiomyocyte')
SpatialFeaturePlot(VISIUM_MI_1D, features = c("Cardiomyocyte1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = c("Cardiomyocyte1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = c("Cardiomyocyte1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = c("Cardiomyocyte1"), alpha = c(0.1,1.5))

VISIUM_MI_1D <- AddModuleScore(object = VISIUM_MI_1D, features = Endothelial_genes, name = 'Endothelial')
VISIUM_MI_3D <- AddModuleScore(object = VISIUM_MI_3D, features = Endothelial_genes, name = 'Endothelial')
VISIUM_MI_5D <- AddModuleScore(object = VISIUM_MI_5D, features = Endothelial_genes, name = 'Endothelial')
VISIUM_MI_7D <- AddModuleScore(object = VISIUM_MI_7D, features = Endothelial_genes, name = 'Endothelial')
SpatialFeaturePlot(VISIUM_MI_1D, features = c("Endothelial1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = c("Endothelial1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = c("Endothelial1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = c("Endothelial1"), alpha = c(0.1,1.5))

VISIUM_MI_1D <- AddModuleScore(object = VISIUM_MI_1D, features = Fibro_genes, name = 'Fibro')
VISIUM_MI_3D <- AddModuleScore(object = VISIUM_MI_3D, features = Fibro_genes, name = 'Fibro')
VISIUM_MI_5D <- AddModuleScore(object = VISIUM_MI_5D, features = Fibro_genes, name = 'Fibro')
VISIUM_MI_7D <- AddModuleScore(object = VISIUM_MI_7D, features = Fibro_genes, name = 'Fibro')
SpatialFeaturePlot(VISIUM_MI_1D, features = c("Fibro1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = c("Fibro1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = c("Fibro1"), alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = c("Fibro1"), alpha = c(0.1,1.5))

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
new.cluster.ids <- c(0,0,0,1,0,0,1,0,0)
VISIUM_MI_1D@active.ident <- plyr::mapvalues(x = VISIUM_MI_1D@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VISIUM_MI_1D <- StashIdent(VISIUM_MI_1D, save.name = "FigureCluster")
SpatialDimPlot(VISIUM_MI_1D, label = FALSE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "FigureCluster", cols = c("0" = "#F8736C", "1" = "#619CFF"))

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
new.cluster.ids <- c(0,0,0,0,2,2,1,0,0)
VISIUM_MI_3D@active.ident <- plyr::mapvalues(x = VISIUM_MI_3D@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VISIUM_MI_3D <- StashIdent(VISIUM_MI_3D, save.name = "FigureCluster")
SpatialDimPlot(VISIUM_MI_3D, label = FALSE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "FigureCluster", cols = c("0" = "#F8736C", "1" = "#619CFF", "2" = "#00BA38"))

current.cluster.ids <- c(0,1,2,3,4,5)
new.cluster.ids <- c(0,0,0,0,0,2)
VISIUM_MI_5D@active.ident <- plyr::mapvalues(x = VISIUM_MI_5D@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VISIUM_MI_5D <- StashIdent(VISIUM_MI_5D, save.name = "FigureCluster")
SpatialDimPlot(VISIUM_MI_5D, label = FALSE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "FigureCluster", cols = c("0" = "#F8736C", "2" = "#00BA38"))

current.cluster.ids <- c(0,1,2,3,4,5,6,7)
new.cluster.ids <- c(0,0,2,2,0,0,0,0)
VISIUM_MI_7D@active.ident <- plyr::mapvalues(x = VISIUM_MI_7D@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VISIUM_MI_7D <- StashIdent(VISIUM_MI_7D, save.name = "FigureCluster")
SpatialDimPlot(VISIUM_MI_7D, label = FALSE, label.size = 10, pt.size.factor = 1.5, alpha = c(0.1, 1), group.by = "FigureCluster", cols = c("0" = "#F8736C", "2" = "#00BA38"))



# Figure 3e and Supplementary Figure 10
SpatialFeaturePlot(VISIUM_MI_1D, features = "Clec4e", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Clec4e", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Clec4e", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Clec4e", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "Hmox1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Hmox1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Hmox1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Hmox1", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "Trem2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Trem2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Trem2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Trem2", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "Lyve1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Lyve1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Lyve1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Lyve1", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "F13a1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "F13a1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "F13a1", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "F13a1", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "Cbr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Cbr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Cbr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Cbr2", alpha = c(0.1,1.5))

SpatialFeaturePlot(VISIUM_MI_1D, features = "Folr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_3D, features = "Folr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_5D, features = "Folr2", alpha = c(0.1,1.5))
SpatialFeaturePlot(VISIUM_MI_7D, features = "Folr2", alpha = c(0.1,1.5))



# SPOTlight analyses
# Day1 Post-MI
set.seed(1243)
sample_id <- "MI1D"
clust_vr <- "CellType"
trn <- "MI"
cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}

# Define cell types of interest
ct_interest <- c("Monocytes", "Neutrophil", "MAC")

# ST-seq Data
se_obj <- "{an_MI}/{robj_dir}/Processed_ST_{sample_id}.RDS" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(file = .)


# Create a name/color dataframe
col_pal_v <- c("#f8766d", "#de8c00", "#b79f00", "#7cae00", "#00ba38", "#00c08b", "#00bfc4", "#00b4f0", "#619cff", "#c77cff", "#f564e3", "#ff64b0")
MI_order <- levels(MI.integrated.1st.filtered)
col_df <- data.frame(ct_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                    x = MI_order,
                                    perl = TRUE),
                     plt_name = MI_order,
                     ct_col = col_pal_v[1:length(MI_order)],
                     stringsAsFactors = FALSE)

col_df <- arrange(col_df, ct_name) 

# Deonvolution
decon_mtrx_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = MI.integrated.1st.filtered,
  counts_spatial = se_obj@assays$Spatial@counts,
  clust_vr = "specific_cell_type_mod",
  cluster_markers = CellType_Markers1,
  cl_n = cl_n,
  hvg = hvg,
  ntop = ntop,
  transf = transf,
  method = method,
  min_cont = min_cont,
  assay = "RNA",
  slot = "counts")

"{an_MI}/{robj_dir}/decon_mtrx_atlas_{sample_id}_{spotlight_id}_{clust_vr}.rds" %>%
  glue::glue() %>%
  here::here() %>%
  saveRDS(object = decon_mtrx_ls, file = .)

decon_mtrx <- decon_mtrx_ls[[2]]
decon_mtrx <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]

# Set as 0 cell types predicted to be under 3 % of the spot
decon_mtrx[decon_mtrx < 0.03] <- 0

# Change names to original ones
new_colnames <- data.frame(ct_name = colnames(decon_mtrx), stringsAsFactors = FALSE) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)
colnames(decon_mtrx) <- new_colnames

# Add deconvolution matrix to Seurat object metadata
se_obj@meta.data <- cbind(se_obj@meta.data, round(decon_mtrx, 3))

# Visualization
h <- NMF::coef(decon_mtrx_ls[[1]][[1]])

train_labs <- data.frame(ct_name = decon_mtrx_ls[[1]][[2]]) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)

rownames(h) <- paste("Topic", 1:nrow(h), sep = " ")

profile_plt <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = train_labs)

plt_2 <- profile_plt[[2]] +
  ggplot2::scale_x_discrete(limits = unique(train_labs)) +
  ggplot2::scale_y_discrete(labels = glue::glue("Topic {1:nrow(h)}")) +
  ggplot2::labs(x = "", y = "") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(hjust = 1)
  )

profile_plt[[1]]

"{an_MI}/{plt_dir}/all_ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = profile_plt[[1]],
    base_height = 16,
    base_width = 16)

"{an_MI}/{plt_dir}/ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = plt_2 +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 18)),
    base_height = 16,
    base_width = 16)

basis_spotlight <- data.frame(NMF::basis(decon_mtrx_ls[[1]][[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(decon_mtrx_ls[[1]][[2]], width = 30))

basis_spotlight %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")


scaleFUN <- function(x) sprintf("%.2f", x)

ct_all <- colnames(decon_mtrx)

ct_plt_ls <- lapply(ct_all, function(ct) {
  tmp_plt <- Seurat::SpatialFeaturePlot(object = se_obj,
                                        features = ct,
                                        alpha = c(0, 1))
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0, 1))
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = ct) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::heat.colors(10, rev = TRUE),
      # Same number of breaks for all plots
       breaks = seq(min(se_obj@meta.data[, ct]),
                    max(se_obj@meta.data[, ct]),
                    length.out = 3),
      # 2 decimals in the legend
      labels = scaleFUN
      # limits = c(0, 1)
    )
  return(tmp_plt)
})

ct_grid <- cowplot::plot_grid(
  plotlist = ct_plt_ls,
  axis = "trbl",
  align = "hv",
  nrow = 3,
  ncol = 4)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = ct_grid,
    base_height = 12,
    base_width = 16)

ct_plt_ls <- lapply(ct_interest, function(ct) {
  print(ct)
  tmp_plt <- Seurat::SpatialPlot(
    object = se_obj,
    features = ct,
    alpha = c(0, 1)) +
    theme(
      legend.title = ggplot2::element_blank())
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::heat.colors(10, rev = TRUE),
        breaks = seq(min(se_obj@meta.data[, ct]),
                     max(se_obj@meta.data[, ct]),
                     length.out = 3),
        labels = scaleFUN
      )
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::labs(title = ct) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10))
  return(tmp_plt)
})

interest_grid <- cowplot::plot_grid(plotlist = ct_plt_ls,
                                    axis = "trbl",
                                    align = "hv",
                                    nrow = 1,
                                    ncol = 3)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = interest_grid,
    base_height = 10,
    base_width = 16)

# Spatial scatterpie
spsct_plt1 <- SPOTlight::spatial_scatterpie(
  se_obj = se_obj,
  cell_types_all = ct_all,
  img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
  pie_scale = 0.4,
  slice = sample_id) + 
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))  


"{an_MI}/{plt_dir}/{sample_id}_Spatial_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = spsct_plt1,
    base_width = 12,
    base_height = 9)


sct_plt1 <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                       cell_types_all = ct_all,
                                       pie_scale = 0.6,
                                       slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt1,
    base_width = 12,
    base_height = 9)

metadata_subset <- se_obj@meta.data[, ct_all]

keep_0.75 <- colSums(se_obj@meta.data[, ct_all] > 0) < 0.75 * ncol(se_obj)
keep_g0 <- colSums(se_obj@meta.data[, ct_all] > 0) > 0

ct_var <- c(colnames(decon_mtrx)[keep_0.75 & keep_g0])

sct_plt_int <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                          cell_types_all = ct_var,
                                          pie_scale = 0.4) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::labs(fill = "") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int,
    base_width = 12,
    base_height = 9)


sct_plt_int1 <- SPOTlight::spatial_scatterpie(se_obj = se_obj,
                                             cell_types_all = ct_var,
                                             img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
                                             pie_scale = 0.4,
                                             slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 12))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest1_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int1,
    base_width = 12,
    base_height = 9)

if (sample_id == "MI1D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI3D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI5D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI5D1") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.5"] #5
} else if (sample_id == "MI5D2") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI7D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
}

# Format the data into long format
metadata_long <- se_obj@meta.data %>% 
  tidyr::pivot_longer(cols = all_of(ct_all),
                      names_to = "immune_key",
                      values_to = "immune_val") %>%
  dplyr::left_join(col_df,
                   by = c("immune_key" = "ct_name")) %>%
  dplyr::mutate(immune_val = dplyr::if_else(immune_val > 0.001, immune_val, 0))

"{an_MI}/{robj_dir}/metadata_long_{sample_id}.csv" %>%
  glue::glue() %>%
  here::here() %>%
  write.csv(metadata_long, file = .)

# Visualization
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1 + 2 * s)
    round(seq(min(x) + d, max(x) - d, length = n), 2)
  }
}

keep_ct <- metadata_long %>%
  dplyr::group_by(immune_key) %>%
  dplyr::summarise(prop_sum = sum(immune_val)) %>% 
  dplyr::filter(prop_sum > 0) %>%
  dplyr::pull(immune_key)

bplt_tmp <- metadata_long %>%
  dplyr::filter(immune_key %in% keep_ct) %>%
  dplyr::mutate(
    immune_key = stringr::str_wrap(string = immune_key,
                                   width = 8)) %>%
  ggpubr::ggboxplot(data = .,
                    x = "Stratification",
                    y = "immune_val",
                    facet.by = "immune_key",
                    color = "Stratification",
                    fill = "Stratification",
                    add.params = list(size=0), #????
                    add = "jitter",
                    scales = "free",
                    repel = TRUE,
                    outlier.shape = NA,
                    alpha = 0.6,
                    palette = c("#f8766d", "#de8c00", "#7cae00", "#00bfc4", "#c77cff") # ("0" = "#f8766d", "1" = "#de8c00", "2" = "#7cae00", "3" = "#00bfc4", "4" = "#c77cff", "5" = "#b79f00", "6" = "#00ba38", "7" = "#00c08b", "8" = "#00b4f0", "9" = "#619cff", "10" = "#f564e3", "11" = "#ff64b0")
  ) +
  ggplot2::scale_y_continuous(
    breaks = equal_breaks(n = 3, s = 0.05),
    expand = c(0.05, 0),
    labels = function(x) sapply(x, FUN = function(i) format(x = i, nsmall = 0))) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 12),
    strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = "Stratification",
                y = "Proportion",
                color = "Regions",
                fill = "Regions")

# Add P values
#bplt_tmp <- bplt_tmp +
#  ggpubr::stat_compare_means(method = "anova", size = 6)

"{an_MI}/{plt_dir}/strat_bplot_oro_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = bplt_tmp,
    base_height = 5,
    base_width = 8)



# Day3 Post-MI
set.seed(1243)
sample_id <- "MI3D"
clust_vr <- "CellType"
trn <- "MI"
cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}

# Define cell types of interest
ct_interest <- c("Monocytes", "Neutrophil", "MAC")

# ST-seq Data
se_obj <- "{an_MI}/{robj_dir}/Processed_ST_{sample_id}.RDS" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(file = .)


# Create a name/color dataframe
col_pal_v <- c("#f8766d", "#de8c00", "#b79f00", "#7cae00", "#00ba38", "#00c08b", "#00bfc4", "#00b4f0", "#619cff", "#c77cff", "#f564e3", "#ff64b0")
MI_order <- levels(MI.integrated.1st.filtered)
col_df <- data.frame(ct_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                    x = MI_order,
                                    perl = TRUE),
                     plt_name = MI_order,
                     ct_col = col_pal_v[1:length(MI_order)],
                     stringsAsFactors = FALSE)

col_df <- arrange(col_df, ct_name) 

# Deonvolution
decon_mtrx_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = MI.integrated.1st.filtered,
  counts_spatial = se_obj@assays$Spatial@counts,
  clust_vr = "specific_cell_type_mod",
  cluster_markers = CellType_Markers1,
  cl_n = cl_n,
  hvg = hvg,
  ntop = ntop,
  transf = transf,
  method = method,
  min_cont = min_cont,
  assay = "RNA",
  slot = "counts")

"{an_MI}/{robj_dir}/decon_mtrx_atlas_{sample_id}_{spotlight_id}_{clust_vr}.rds" %>%
  glue::glue() %>%
  here::here() %>%
  saveRDS(object = decon_mtrx_ls, file = .)

decon_mtrx <- decon_mtrx_ls[[2]]
decon_mtrx <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]

# Set as 0 cell types predicted to be under 3 % of the spot
decon_mtrx[decon_mtrx < 0.03] <- 0

# Change names to original ones
new_colnames <- data.frame(ct_name = colnames(decon_mtrx), stringsAsFactors = FALSE) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)
colnames(decon_mtrx) <- new_colnames

# Add deconvolution matrix to Seurat object metadata
se_obj@meta.data <- cbind(se_obj@meta.data, round(decon_mtrx, 3))

# Visualization
h <- NMF::coef(decon_mtrx_ls[[1]][[1]])

train_labs <- data.frame(ct_name = decon_mtrx_ls[[1]][[2]]) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)

rownames(h) <- paste("Topic", 1:nrow(h), sep = " ")

profile_plt <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = train_labs)

plt_2 <- profile_plt[[2]] +
  ggplot2::scale_x_discrete(limits = unique(train_labs)) +
  ggplot2::scale_y_discrete(labels = glue::glue("Topic {1:nrow(h)}")) +
  ggplot2::labs(x = "", y = "") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(hjust = 1)
  )

profile_plt[[1]]

"{an_MI}/{plt_dir}/all_ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = profile_plt[[1]],
    base_height = 16,
    base_width = 16)

"{an_MI}/{plt_dir}/ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = plt_2 +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 18)),
    base_height = 16,
    base_width = 16)

basis_spotlight <- data.frame(NMF::basis(decon_mtrx_ls[[1]][[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(decon_mtrx_ls[[1]][[2]], width = 30))

basis_spotlight %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")


scaleFUN <- function(x) sprintf("%.2f", x)

ct_all <- colnames(decon_mtrx)

ct_plt_ls <- lapply(ct_all, function(ct) {
  tmp_plt <- Seurat::SpatialFeaturePlot(object = se_obj,
                                        features = ct,
                                        alpha = c(0, 1))
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0, 1))
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = ct) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::heat.colors(10, rev = TRUE),
      # Same number of breaks for all plots
      breaks = seq(min(se_obj@meta.data[, ct]),
                   max(se_obj@meta.data[, ct]),
                   length.out = 3),
      # 2 decimals in the legend
      labels = scaleFUN
      # limits = c(0, 1)
    )
  return(tmp_plt)
})

ct_grid <- cowplot::plot_grid(
  plotlist = ct_plt_ls,
  axis = "trbl",
  align = "hv",
  nrow = 3,
  ncol = 4)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = ct_grid,
    base_height = 12,
    base_width = 16)

ct_plt_ls <- lapply(ct_interest, function(ct) {
  print(ct)
  tmp_plt <- Seurat::SpatialPlot(
    object = se_obj,
    features = ct,
    alpha = c(0, 1)) +
    theme(
      legend.title = ggplot2::element_blank())
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::heat.colors(10, rev = TRUE),
        breaks = seq(min(se_obj@meta.data[, ct]),
                     max(se_obj@meta.data[, ct]),
                     length.out = 3),
        labels = scaleFUN
      )
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::labs(title = ct) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10))
  return(tmp_plt)
})

interest_grid <- cowplot::plot_grid(plotlist = ct_plt_ls,
                                    axis = "trbl",
                                    align = "hv",
                                    nrow = 1,
                                    ncol = 3)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = interest_grid,
    base_height = 10,
    base_width = 16)

# Spatial scatterpie
spsct_plt1 <- SPOTlight::spatial_scatterpie(
  se_obj = se_obj,
  cell_types_all = ct_all,
  img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
  pie_scale = 0.4,
  slice = sample_id) + 
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))  


"{an_MI}/{plt_dir}/{sample_id}_Spatial_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = spsct_plt1,
    base_width = 12,
    base_height = 9)


sct_plt1 <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                       cell_types_all = ct_all,
                                       pie_scale = 0.6,
                                       slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt1,
    base_width = 12,
    base_height = 9)

metadata_subset <- se_obj@meta.data[, ct_all]

keep_0.75 <- colSums(se_obj@meta.data[, ct_all] > 0) < 0.75 * ncol(se_obj)
keep_g0 <- colSums(se_obj@meta.data[, ct_all] > 0) > 0

ct_var <- c(colnames(decon_mtrx)[keep_0.75 & keep_g0])

sct_plt_int <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                          cell_types_all = ct_var,
                                          pie_scale = 0.4) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::labs(fill = "") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int,
    base_width = 12,
    base_height = 9)


sct_plt_int1 <- SPOTlight::spatial_scatterpie(se_obj = se_obj,
                                              cell_types_all = ct_var,
                                              img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
                                              pie_scale = 0.4,
                                              slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 12))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest1_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int1,
    base_width = 12,
    base_height = 9)

if (sample_id == "MI1D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI3D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI5D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI5D1") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.5"] #5
} else if (sample_id == "MI5D2") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI7D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
}

# Format the data into long format
metadata_long <- se_obj@meta.data %>% 
  tidyr::pivot_longer(cols = all_of(ct_all),
                      names_to = "immune_key",
                      values_to = "immune_val") %>%
  dplyr::left_join(col_df,
                   by = c("immune_key" = "ct_name")) %>%
  dplyr::mutate(immune_val = dplyr::if_else(immune_val > 0.001, immune_val, 0))

"{an_MI}/{robj_dir}/metadata_long_{sample_id}.csv" %>%
  glue::glue() %>%
  here::here() %>%
  write.csv(metadata_long, file = .)

# Visualization
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1 + 2 * s)
    round(seq(min(x) + d, max(x) - d, length = n), 2)
  }
}

keep_ct <- metadata_long %>%
  dplyr::group_by(immune_key) %>%
  dplyr::summarise(prop_sum = sum(immune_val)) %>% 
  dplyr::filter(prop_sum > 0) %>%
  dplyr::pull(immune_key)

bplt_tmp <- metadata_long %>%
  dplyr::filter(immune_key %in% keep_ct) %>%
  dplyr::mutate(
    immune_key = stringr::str_wrap(string = immune_key,
                                   width = 8)) %>%
  ggpubr::ggboxplot(data = .,
                    x = "Stratification",
                    y = "immune_val",
                    facet.by = "immune_key",
                    color = "Stratification",
                    fill = "Stratification",
                    add.params = list(size=0), #????
                    add = "jitter",
                    scales = "free",
                    repel = TRUE,
                    outlier.shape = NA,
                    alpha = 0.6,
                    palette = c("#f8766d", "#de8c00", "#7cae00", "#00bfc4", "#c77cff") # ("0" = "#f8766d", "1" = "#de8c00", "2" = "#7cae00", "3" = "#00bfc4", "4" = "#c77cff", "5" = "#b79f00", "6" = "#00ba38", "7" = "#00c08b", "8" = "#00b4f0", "9" = "#619cff", "10" = "#f564e3", "11" = "#ff64b0")
  ) +
  ggplot2::scale_y_continuous(
    breaks = equal_breaks(n = 3, s = 0.05),
    expand = c(0.05, 0),
    labels = function(x) sapply(x, FUN = function(i) format(x = i, nsmall = 0))) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 12),
    strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = "Stratification",
                y = "Proportion",
                color = "Regions",
                fill = "Regions")

# Add P values
#bplt_tmp <- bplt_tmp +
#  ggpubr::stat_compare_means(method = "anova", size = 6)

"{an_MI}/{plt_dir}/strat_bplot_oro_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = bplt_tmp,
    base_height = 5,
    base_width = 8)



# Day5 Post-MI
set.seed(1243)
sample_id <- "MI5D"
clust_vr <- "CellType"
trn <- "MI"
cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}

# Define cell types of interest
ct_interest <- c("Monocytes", "Neutrophil", "MAC")

# ST-seq Data
se_obj <- "{an_MI}/{robj_dir}/Processed_ST_{sample_id}.RDS" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(file = .)


# Create a name/color dataframe
col_pal_v <- c("#f8766d", "#de8c00", "#b79f00", "#7cae00", "#00ba38", "#00c08b", "#00bfc4", "#00b4f0", "#619cff", "#c77cff", "#f564e3", "#ff64b0")
MI_order <- levels(MI.integrated.1st.filtered)
col_df <- data.frame(ct_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                    x = MI_order,
                                    perl = TRUE),
                     plt_name = MI_order,
                     ct_col = col_pal_v[1:length(MI_order)],
                     stringsAsFactors = FALSE)

col_df <- arrange(col_df, ct_name) 

# Deonvolution
decon_mtrx_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = MI.integrated.1st.filtered,
  counts_spatial = se_obj@assays$Spatial@counts,
  clust_vr = "specific_cell_type_mod",
  cluster_markers = CellType_Markers1,
  cl_n = cl_n,
  hvg = hvg,
  ntop = ntop,
  transf = transf,
  method = method,
  min_cont = min_cont,
  assay = "RNA",
  slot = "counts")

"{an_MI}/{robj_dir}/decon_mtrx_atlas_{sample_id}_{spotlight_id}_{clust_vr}.rds" %>%
  glue::glue() %>%
  here::here() %>%
  saveRDS(object = decon_mtrx_ls, file = .)

decon_mtrx <- decon_mtrx_ls[[2]]
decon_mtrx <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]

# Set as 0 cell types predicted to be under 3 % of the spot
decon_mtrx[decon_mtrx < 0.03] <- 0

# Change names to original ones
new_colnames <- data.frame(ct_name = colnames(decon_mtrx), stringsAsFactors = FALSE) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)
colnames(decon_mtrx) <- new_colnames

# Add deconvolution matrix to Seurat object metadata
se_obj@meta.data <- cbind(se_obj@meta.data, round(decon_mtrx, 3))

# Visualization
h <- NMF::coef(decon_mtrx_ls[[1]][[1]])

train_labs <- data.frame(ct_name = decon_mtrx_ls[[1]][[2]]) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)

rownames(h) <- paste("Topic", 1:nrow(h), sep = " ")

profile_plt <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = train_labs)

plt_2 <- profile_plt[[2]] +
  ggplot2::scale_x_discrete(limits = unique(train_labs)) +
  ggplot2::scale_y_discrete(labels = glue::glue("Topic {1:nrow(h)}")) +
  ggplot2::labs(x = "", y = "") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(hjust = 1)
  )

profile_plt[[1]]

"{an_MI}/{plt_dir}/all_ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = profile_plt[[1]],
    base_height = 16,
    base_width = 16)

"{an_MI}/{plt_dir}/ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = plt_2 +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 18)),
    base_height = 16,
    base_width = 16)

basis_spotlight <- data.frame(NMF::basis(decon_mtrx_ls[[1]][[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(decon_mtrx_ls[[1]][[2]], width = 30))

basis_spotlight %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")


scaleFUN <- function(x) sprintf("%.2f", x)

ct_all <- colnames(decon_mtrx)

ct_plt_ls <- lapply(ct_all, function(ct) {
  tmp_plt <- Seurat::SpatialFeaturePlot(object = se_obj,
                                        features = ct,
                                        alpha = c(0, 1))
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0, 1))
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = ct) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::heat.colors(10, rev = TRUE),
      # Same number of breaks for all plots
      breaks = seq(min(se_obj@meta.data[, ct]),
                   max(se_obj@meta.data[, ct]),
                   length.out = 3),
      # 2 decimals in the legend
      labels = scaleFUN
      # limits = c(0, 1)
    )
  return(tmp_plt)
})

ct_grid <- cowplot::plot_grid(
  plotlist = ct_plt_ls,
  axis = "trbl",
  align = "hv",
  nrow = 3,
  ncol = 4)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = ct_grid,
    base_height = 12,
    base_width = 16)

ct_plt_ls <- lapply(ct_interest, function(ct) {
  print(ct)
  tmp_plt <- Seurat::SpatialPlot(
    object = se_obj,
    features = ct,
    alpha = c(0, 1)) +
    theme(
      legend.title = ggplot2::element_blank())
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::heat.colors(10, rev = TRUE),
        breaks = seq(min(se_obj@meta.data[, ct]),
                     max(se_obj@meta.data[, ct]),
                     length.out = 3),
        labels = scaleFUN
      )
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::labs(title = ct) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10))
  return(tmp_plt)
})

interest_grid <- cowplot::plot_grid(plotlist = ct_plt_ls,
                                    axis = "trbl",
                                    align = "hv",
                                    nrow = 1,
                                    ncol = 3)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = interest_grid,
    base_height = 10,
    base_width = 16)

# Spatial scatterpie
spsct_plt1 <- SPOTlight::spatial_scatterpie(
  se_obj = se_obj,
  cell_types_all = ct_all,
  img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
  pie_scale = 0.4,
  slice = sample_id) + 
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))  


"{an_MI}/{plt_dir}/{sample_id}_Spatial_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = spsct_plt1,
    base_width = 12,
    base_height = 9)


sct_plt1 <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                       cell_types_all = ct_all,
                                       pie_scale = 0.6,
                                       slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt1,
    base_width = 12,
    base_height = 9)

metadata_subset <- se_obj@meta.data[, ct_all]

keep_0.75 <- colSums(se_obj@meta.data[, ct_all] > 0) < 0.75 * ncol(se_obj)
keep_g0 <- colSums(se_obj@meta.data[, ct_all] > 0) > 0

ct_var <- c(colnames(decon_mtrx)[keep_0.75 & keep_g0])

sct_plt_int <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                          cell_types_all = ct_var,
                                          pie_scale = 0.4) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::labs(fill = "") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int,
    base_width = 12,
    base_height = 9)


sct_plt_int1 <- SPOTlight::spatial_scatterpie(se_obj = se_obj,
                                              cell_types_all = ct_var,
                                              img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
                                              pie_scale = 0.4,
                                              slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 12))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest1_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int1,
    base_width = 12,
    base_height = 9)

if (sample_id == "MI1D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI3D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI5D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI5D1") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.5"] #5
} else if (sample_id == "MI5D2") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI7D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
}

# Format the data into long format
metadata_long <- se_obj@meta.data %>% 
  tidyr::pivot_longer(cols = all_of(ct_all),
                      names_to = "immune_key",
                      values_to = "immune_val") %>%
  dplyr::left_join(col_df,
                   by = c("immune_key" = "ct_name")) %>%
  dplyr::mutate(immune_val = dplyr::if_else(immune_val > 0.001, immune_val, 0))

"{an_MI}/{robj_dir}/metadata_long_{sample_id}.csv" %>%
  glue::glue() %>%
  here::here() %>%
  write.csv(metadata_long, file = .)

# Visualization
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1 + 2 * s)
    round(seq(min(x) + d, max(x) - d, length = n), 2)
  }
}

keep_ct <- metadata_long %>%
  dplyr::group_by(immune_key) %>%
  dplyr::summarise(prop_sum = sum(immune_val)) %>% 
  dplyr::filter(prop_sum > 0) %>%
  dplyr::pull(immune_key)

bplt_tmp <- metadata_long %>%
  dplyr::filter(immune_key %in% keep_ct) %>%
  dplyr::mutate(
    immune_key = stringr::str_wrap(string = immune_key,
                                   width = 8)) %>%
  ggpubr::ggboxplot(data = .,
                    x = "Stratification",
                    y = "immune_val",
                    facet.by = "immune_key",
                    color = "Stratification",
                    fill = "Stratification",
                    add.params = list(size=0), #????
                    add = "jitter",
                    scales = "free",
                    repel = TRUE,
                    outlier.shape = NA,
                    alpha = 0.6,
                    palette = c("#f8766d", "#de8c00", "#7cae00", "#00bfc4", "#c77cff") # ("0" = "#f8766d", "1" = "#de8c00", "2" = "#7cae00", "3" = "#00bfc4", "4" = "#c77cff", "5" = "#b79f00", "6" = "#00ba38", "7" = "#00c08b", "8" = "#00b4f0", "9" = "#619cff", "10" = "#f564e3", "11" = "#ff64b0")
  ) +
  ggplot2::scale_y_continuous(
    breaks = equal_breaks(n = 3, s = 0.05),
    expand = c(0.05, 0),
    labels = function(x) sapply(x, FUN = function(i) format(x = i, nsmall = 0))) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 12),
    strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = "Stratification",
                y = "Proportion",
                color = "Regions",
                fill = "Regions")

# Add P values
#bplt_tmp <- bplt_tmp +
#  ggpubr::stat_compare_means(method = "anova", size = 6)

"{an_MI}/{plt_dir}/strat_bplot_oro_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = bplt_tmp,
    base_height = 5,
    base_width = 8)




# Day7 Post-MI
set.seed(1243)
sample_id <- "MI7D"
clust_vr <- "CellType"
trn <- "MI"
cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}

# Define cell types of interest
ct_interest <- c("Monocytes", "Neutrophil", "MAC")

# ST-seq Data
se_obj <- "{an_MI}/{robj_dir}/Processed_ST_{sample_id}.RDS" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(file = .)


# Create a name/color dataframe
col_pal_v <- c("#f8766d", "#de8c00", "#b79f00", "#7cae00", "#00ba38", "#00c08b", "#00bfc4", "#00b4f0", "#619cff", "#c77cff", "#f564e3", "#ff64b0")
MI_order <- levels(MI.integrated.1st.filtered)
col_df <- data.frame(ct_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                    x = MI_order,
                                    perl = TRUE),
                     plt_name = MI_order,
                     ct_col = col_pal_v[1:length(MI_order)],
                     stringsAsFactors = FALSE)

col_df <- arrange(col_df, ct_name) 

# Deonvolution
decon_mtrx_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = MI.integrated.1st.filtered,
  counts_spatial = se_obj@assays$Spatial@counts,
  clust_vr = "specific_cell_type_mod",
  cluster_markers = CellType_Markers1,
  cl_n = cl_n,
  hvg = hvg,
  ntop = ntop,
  transf = transf,
  method = method,
  min_cont = min_cont,
  assay = "RNA",
  slot = "counts")

"{an_MI}/{robj_dir}/decon_mtrx_atlas_{sample_id}_{spotlight_id}_{clust_vr}.rds" %>%
  glue::glue() %>%
  here::here() %>%
  saveRDS(object = decon_mtrx_ls, file = .)

decon_mtrx <- decon_mtrx_ls[[2]]
decon_mtrx <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]

# Set as 0 cell types predicted to be under 3 % of the spot
decon_mtrx[decon_mtrx < 0.03] <- 0

# Change names to original ones
new_colnames <- data.frame(ct_name = colnames(decon_mtrx), stringsAsFactors = FALSE) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)
colnames(decon_mtrx) <- new_colnames

# Add deconvolution matrix to Seurat object metadata
se_obj@meta.data <- cbind(se_obj@meta.data, round(decon_mtrx, 3))

# Visualization
h <- NMF::coef(decon_mtrx_ls[[1]][[1]])

train_labs <- data.frame(ct_name = decon_mtrx_ls[[1]][[2]]) %>%
  dplyr::left_join(col_df, by = "ct_name") %>%
  dplyr::pull(plt_name)

rownames(h) <- paste("Topic", 1:nrow(h), sep = " ")

profile_plt <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = train_labs)

plt_2 <- profile_plt[[2]] +
  ggplot2::scale_x_discrete(limits = unique(train_labs)) +
  ggplot2::scale_y_discrete(labels = glue::glue("Topic {1:nrow(h)}")) +
  ggplot2::labs(x = "", y = "") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(hjust = 1)
  )

profile_plt[[1]]

"{an_MI}/{plt_dir}/all_ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = profile_plt[[1]],
    base_height = 16,
    base_width = 16)

"{an_MI}/{plt_dir}/ct_profiles_{spotlight_id}_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = plt_2 +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 18)),
    base_height = 16,
    base_width = 16)

basis_spotlight <- data.frame(NMF::basis(decon_mtrx_ls[[1]][[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(decon_mtrx_ls[[1]][[2]], width = 30))

basis_spotlight %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")


scaleFUN <- function(x) sprintf("%.2f", x)

ct_all <- colnames(decon_mtrx)

ct_plt_ls <- lapply(ct_all, function(ct) {
  tmp_plt <- Seurat::SpatialFeaturePlot(object = se_obj,
                                        features = ct,
                                        alpha = c(0, 1))
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0, 1))
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = ct) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::heat.colors(10, rev = TRUE),
      # Same number of breaks for all plots
      breaks = seq(min(se_obj@meta.data[, ct]),
                   max(se_obj@meta.data[, ct]),
                   length.out = 3),
      # 2 decimals in the legend
      labels = scaleFUN
      # limits = c(0, 1)
    )
  return(tmp_plt)
})

ct_grid <- cowplot::plot_grid(
  plotlist = ct_plt_ls,
  axis = "trbl",
  align = "hv",
  nrow = 3,
  ncol = 4)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = ct_grid,
    base_height = 12,
    base_width = 16)

ct_plt_ls <- lapply(ct_interest, function(ct) {
  print(ct)
  tmp_plt <- Seurat::SpatialPlot(
    object = se_obj,
    features = ct,
    alpha = c(0, 1)) +
    theme(
      legend.title = ggplot2::element_blank())
  
  if (sum(se_obj@meta.data[, ct]) == 0) {
    tmp_plt <- tmp_plt + ggplot2::scale_alpha(range = c(0,0))
  } else {
    tmp_plt <- tmp_plt +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::heat.colors(10, rev = TRUE),
        breaks = seq(min(se_obj@meta.data[, ct]),
                     max(se_obj@meta.data[, ct]),
                     length.out = 3),
        labels = scaleFUN
      )
  }
  
  tmp_plt <- tmp_plt +
    ggplot2::labs(title = ct) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10))
  return(tmp_plt)
})

interest_grid <- cowplot::plot_grid(plotlist = ct_plt_ls,
                                    axis = "trbl",
                                    align = "hv",
                                    nrow = 1,
                                    ncol = 3)

"{an_MI}/{plt_dir}/{sample_id}_scRNA_reference_arrangement_interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = interest_grid,
    base_height = 10,
    base_width = 16)

# Spatial scatterpie
spsct_plt1 <- SPOTlight::spatial_scatterpie(
  se_obj = se_obj,
  cell_types_all = ct_all,
  img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
  pie_scale = 0.4,
  slice = sample_id) + 
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))  


"{an_MI}/{plt_dir}/{sample_id}_Spatial_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = spsct_plt1,
    base_width = 12,
    base_height = 9)


sct_plt1 <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                       cell_types_all = ct_all,
                                       pie_scale = 0.6,
                                       slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_all, "ct_col"],
    breaks = ct_all) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt1,
    base_width = 12,
    base_height = 9)

metadata_subset <- se_obj@meta.data[, ct_all]

keep_0.75 <- colSums(se_obj@meta.data[, ct_all] > 0) < 0.75 * ncol(se_obj)
keep_g0 <- colSums(se_obj@meta.data[, ct_all] > 0) > 0

ct_var <- c(colnames(decon_mtrx)[keep_0.75 & keep_g0])

sct_plt_int <- SPOTlight::scatterpie_plot(se_obj = se_obj,
                                          cell_types_all = ct_var,
                                          pie_scale = 0.4) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::labs(fill = "") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int,
    base_width = 12,
    base_height = 9)


sct_plt_int1 <- SPOTlight::spatial_scatterpie(se_obj = se_obj,
                                              cell_types_all = ct_var,
                                              img_path = here::here(sprintf("F:/Single_Cell_R/MI/Visium/Raw_Ref2020/%s/spatial/tissue_lowres_image.png", sample_id)),
                                              pie_scale = 0.4,
                                              slice = sample_id) +
  ggplot2::scale_fill_manual(
    values = col_df[col_df$plt_name %in% ct_var, "ct_col"],
    breaks = col_df[col_df$plt_name %in% ct_var, "plt_name"]) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 12))

"{an_MI}/{plt_dir}/{sample_id}_Scatterpie_Interest1_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = sct_plt_int1,
    base_width = 12,
    base_height = 9)

if (sample_id == "MI1D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI3D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
} else if (sample_id == "MI5D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI5D1") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.5"] #5
} else if (sample_id == "MI5D2") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #4
} else if (sample_id == "MI7D") {
  se_obj@meta.data[, "Stratification"] <- se_obj@meta.data[, "SCT_snn_res.0.25"] #5
}

# Format the data into long format
metadata_long <- se_obj@meta.data %>% 
  tidyr::pivot_longer(cols = all_of(ct_all),
                      names_to = "immune_key",
                      values_to = "immune_val") %>%
  dplyr::left_join(col_df,
                   by = c("immune_key" = "ct_name")) %>%
  dplyr::mutate(immune_val = dplyr::if_else(immune_val > 0.001, immune_val, 0))

"{an_MI}/{robj_dir}/metadata_long_{sample_id}.csv" %>%
  glue::glue() %>%
  here::here() %>%
  write.csv(metadata_long, file = .)

# Visualization
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1 + 2 * s)
    round(seq(min(x) + d, max(x) - d, length = n), 2)
  }
}

keep_ct <- metadata_long %>%
  dplyr::group_by(immune_key) %>%
  dplyr::summarise(prop_sum = sum(immune_val)) %>% 
  dplyr::filter(prop_sum > 0) %>%
  dplyr::pull(immune_key)

bplt_tmp <- metadata_long %>%
  dplyr::filter(immune_key %in% keep_ct) %>%
  dplyr::mutate(
    immune_key = stringr::str_wrap(string = immune_key,
                                   width = 8)) %>%
  ggpubr::ggboxplot(data = .,
                    x = "Stratification",
                    y = "immune_val",
                    facet.by = "immune_key",
                    color = "Stratification",
                    fill = "Stratification",
                    add.params = list(size=0), #????
                    add = "jitter",
                    scales = "free",
                    repel = TRUE,
                    outlier.shape = NA,
                    alpha = 0.6,
                    palette = c("#f8766d", "#de8c00", "#7cae00", "#00bfc4", "#c77cff") # ("0" = "#f8766d", "1" = "#de8c00", "2" = "#7cae00", "3" = "#00bfc4", "4" = "#c77cff", "5" = "#b79f00", "6" = "#00ba38", "7" = "#00c08b", "8" = "#00b4f0", "9" = "#619cff", "10" = "#f564e3", "11" = "#ff64b0")
  ) +
  ggplot2::scale_y_continuous(
    breaks = equal_breaks(n = 3, s = 0.05),
    expand = c(0.05, 0),
    labels = function(x) sapply(x, FUN = function(i) format(x = i, nsmall = 0))) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 12),
    strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = "Stratification",
                y = "Proportion",
                color = "Regions",
                fill = "Regions")

# Add P values
#bplt_tmp <- bplt_tmp +
#  ggpubr::stat_compare_means(method = "anova", size = 6)

"{an_MI}/{plt_dir}/strat_bplot_oro_{sample_id}_{clust_vr}.pdf" %>%
  glue::glue() %>%
  here::here() %>%
  cowplot::save_plot(
    filename = .,
    plot = bplt_tmp,
    base_height = 5,
    base_width = 8)
# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file=file.path(project_root, "results", "Q1_log.txt"), append=TRUE)
}

log("加载T/NK细胞...")
tnk <- readRDS(file.path(project_root, "results", "objects", "tnk_cells.rds"))

log("重新归一化和特征选择...")
tnk <- NormalizeData(tnk)
tnk <- FindVariableFeatures(tnk, nfeatures = 2000)
tnk <- ScaleData(tnk)

log("PCA...")
tnk <- RunPCA(tnk, npcs = 30)

log("Harmony批次校正...")
tnk <- RunHarmony(tnk, group.by.vars = "patient")

log("UMAP...")
tnk <- RunUMAP(tnk, reduction = "harmony", dims = 1:20)

log("聚类...")
tnk <- FindNeighbors(tnk, reduction = "harmony", dims = 1:20)
tnk <- FindClusters(tnk, resolution = 0.4)
log(paste("T/NK簇数:", length(unique(tnk$seurat_clusters))))

log("计算Tpex stemness score...")
stemness_genes <- c("TCF7", "SLAMF6", "CCR7", "IL7R", "SELL", "LEF1")
exhaustion_genes <- c("PDCD1", "TOX", "HAVCR2", "TIGIT", "LAG3", "CTLA4")
cytotox_genes <- c("GZMB", "PRF1", "IFNG", "GNLY", "GZMA", "GZMH")

tnk <- AddModuleScore(tnk, features = list(stemness_genes), name = "stemness_score")
tnk <- AddModuleScore(tnk, features = list(exhaustion_genes), name = "exhaustion_score")
tnk <- AddModuleScore(tnk, features = list(cytotox_genes), name = "cytotox_score")

log("检查关键marker基因表达...")
markers_check <- c("TCF7", "PDCD1", "CXCR5", "TOX", "GZMB", "CD8A", "CD4", "FOXP3", "NCAM1")
expr_check <- AverageExpression(tnk, features = markers_check, group.by = "seurat_clusters")
write.csv(expr_check$RNA, file.path(project_root, "results", "Q1_marker_expression.csv"))

log("保存T/NK重聚类结果...")
saveRDS(tnk, file.path(project_root, "results", "objects", "tnk_clustered.rds"))

log("Q1完成!")

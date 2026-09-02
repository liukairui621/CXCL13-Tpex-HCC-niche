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
      file=file.path(project_root, "results", "Q0_full_log.txt"), append=TRUE)
}

log("加载Seurat对象...")
seurat_obj <- readRDS(file.path(project_root, "results", "objects", "seurat_raw.rds"))

log("质控过滤...")
seurat_obj <- subset(seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           nCount_RNA < 80000 &
           percent.mt < 15
)
log(paste("过滤后细胞数:", ncol(seurat_obj)))

log("归一化...")
seurat_obj <- NormalizeData(seurat_obj)

log("找高变基因...")
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)

log("缩放...")
seurat_obj <- ScaleData(seurat_obj)

log("PCA...")
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

log("Harmony批次校正...")
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "patient")

log("UMAP降维...")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)

log("构建近邻图...")
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)

log("聚类...")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
log(paste("簇数:", length(unique(seurat_obj$seurat_clusters))))

log("保存结果...")
saveRDS(seurat_obj, file.path(project_root, "results", "objects", "seurat_clustered.rds"))

log("全部完成!")

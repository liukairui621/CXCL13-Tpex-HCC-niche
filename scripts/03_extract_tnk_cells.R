# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)

log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file=file.path(project_root, "results", "Q0_annotation_log.txt"), append=TRUE)
}

log("加载聚类结果...")
seurat_obj <- readRDS(file.path(project_root, "results", "objects", "seurat_clustered.rds"))

log("查看各簇细胞类型分布...")
cluster_celltype <- table(seurat_obj$seurat_clusters, seurat_obj$celltype)
write.csv(cluster_celltype, file.path(project_root, "results", "cluster_celltype_table.csv"))

log("查看各簇组织来源分布...")
cluster_site <- table(seurat_obj$seurat_clusters, seurat_obj$site)
write.csv(cluster_site, file.path(project_root, "results", "cluster_site_table.csv"))

log("提取T/NK细胞...")
tnk_cells <- subset(seurat_obj, subset = celltype == "T/NK")
log(paste("T/NK细胞数:", ncol(tnk_cells)))

log("查看T/NK细胞组织来源分布...")
print(table(tnk_cells$site))

log("保存T/NK细胞对象...")
saveRDS(tnk_cells, file.path(project_root, "results", "objects", "tnk_cells.rds"))

log("完成!")

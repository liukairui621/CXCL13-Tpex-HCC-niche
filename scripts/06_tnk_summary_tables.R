# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(dplyr)

log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file=file.path(project_root, "results", "Q1_detail_log.txt"), append=TRUE)
}

log("加载数据...")
tnk <- readRDS(file.path(project_root, "results", "objects", "tnk_clustered.rds"))

log("检查关键簇的详细marker...")
key_genes <- c("TCF7", "PDCD1", "CXCR5", "TOX", "SLAMF6", "GZMB", 
               "PRF1", "HAVCR2", "TIGIT", "CD8A", "CD8B", "SELL",
               "CCR7", "IL7R", "CXCL13", "MKI67")

expr <- AverageExpression(tnk, features = key_genes, group.by = "seurat_clusters")
write.csv(expr$RNA, file.path(project_root, "results", "Q1_detail_markers.csv"))

log("检查各簇在患者间的分布...")
patient_dist <- table(tnk$seurat_clusters, tnk$patient)
write.csv(patient_dist, file.path(project_root, "results", "Q1_patient_distribution.csv"))

log("检查关键簇在不同组织来源的分布...")
site_dist <- table(tnk$seurat_clusters, tnk$site)
write.csv(site_dist, file.path(project_root, "results", "Q1_site_distribution.csv"))

log("完成!")

# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(dplyr)

log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file=file.path(project_root, "results", "Q1_final_log.txt"), append=TRUE)
}

log("加载数据...")
tnk <- readRDS(file.path(project_root, "results", "objects", "tnk_clustered.rds"))

log("标注T/NK细胞亚群...")
tnk$cellstate <- dplyr::recode(as.character(tnk$seurat_clusters),
  "0" = "Tn_Tcm",
  "1" = "Tcm",
  "2" = "Terminal_Tex",
  "3" = "Early_Tex",
  "4" = "Tpex_like",
  "5" = "Treg",
  "6" = "Cycling_T",
  "7" = "Tumor_T1",
  "8" = "Tumor_T2",
  "9" = "Tumor_T3",
  "10" = "Tumor_T4",
  "11" = "Tumor_T5",
  "12" = "Tumor_T6",
  "13" = "Tumor_T7",
  "14" = "NK",
  "15" = "NK_dim"
)

log("Tpex-like细胞统计...")
tpex_cells <- tnk$cellstate == "Tpex_like"
log(paste("Tpex-like总细胞数:", sum(tpex_cells)))
log(paste("肿瘤内Tpex-like:", sum(tpex_cells & tnk$site == "Tumor")))

log("各患者Tpex-like分布...")
tpex_meta <- tnk@meta.data[tpex_cells, ]
patient_tpex <- table(tpex_meta$patient, tpex_meta$site)
write.csv(patient_tpex, file.path(project_root, "results", "Q1_tpex_patient_dist.csv"))

log("保存标注后的T/NK对象...")
saveRDS(tnk, file.path(project_root, "results", "objects", "tnk_annotated.rds"))

log("Q1全部完成!")

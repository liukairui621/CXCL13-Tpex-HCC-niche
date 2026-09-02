# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(data.table)
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(dplyr)

log_file <- file.path(project_root, "results", "Q0_log.txt")
cat("开始时间:", as.character(Sys.time()), "\n", file=log_file)

cat("步骤1: 读取count矩阵...\n", file=log_file, append=TRUE)
counts <- fread(
  file.path(project_root, "data", "GSE149614", "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"),
  sep = "\t"
)
cat("count矩阵读取完成\n", file=log_file, append=TRUE)

cat("步骤2: 转换稀疏矩阵...\n", file=log_file, append=TRUE)
genes <- counts[[1]]
counts[, V1 := NULL]
counts_sparse <- Matrix(as.matrix(counts), sparse = TRUE)
rownames(counts_sparse) <- genes
rm(counts)
gc()
cat("稀疏矩阵转换完成\n", file=log_file, append=TRUE)

cat("步骤3: 读取元数据...\n", file=log_file, append=TRUE)
metadata <- read.table(
  gzfile(file.path(project_root, "data", "GSE149614", "GSE149614_HCC.metadata.updated.txt.gz")),
  header = TRUE, row.names = 1, sep = "\t"
)
cat("元数据读取完成\n", file=log_file, append=TRUE)

cat("步骤4: 创建Seurat对象...\n", file=log_file, append=TRUE)
seurat_obj <- CreateSeuratObject(
  counts = counts_sparse,
  meta.data = metadata,
  min.cells = 3,
  min.features = 200
)
cat("Seurat对象创建完成\n", file=log_file, append=TRUE)
cat("细胞数:", ncol(seurat_obj), "\n", file=log_file, append=TRUE)
cat("基因数:", nrow(seurat_obj), "\n", file=log_file, append=TRUE)

cat("步骤5: 计算线粒体基因比例...\n", file=log_file, append=TRUE)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

cat("步骤6: 保存初始对象...\n", file=log_file, append=TRUE)
saveRDS(seurat_obj, file.path(project_root, "results", "objects", "seurat_raw.rds"))

cat("Q0分析完成! 结束时间:", as.character(Sys.time()), "\n", file=log_file, append=TRUE)

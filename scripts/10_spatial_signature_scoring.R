# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(singscore)

log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file=file.path(project_root, "results", "Q3_log.txt"), append=TRUE)
}

log("定义样本信息...")
samples <- list(
  HCC1R  = list(path=file.path(project_root, "data", "GSE238264", "HCC1R"),  response="Responder"),
  HCC2R  = list(path=file.path(project_root, "data", "GSE238264", "HCC2R"),  response="Responder"),
  HCC3R  = list(path=file.path(project_root, "data", "GSE238264", "HCC3R"),  response="Responder"),
  HCC4R  = list(path=file.path(project_root, "data", "GSE238264", "HCC4R"),  response="Responder"),
  HCC5NR = list(path=file.path(project_root, "data", "GSE238264", "HCC5NR"), response="NonResponder"),
  HCC6NR = list(path=file.path(project_root, "data", "GSE238264", "HCC6NR"), response="NonResponder"),
  HCC7NR = list(path=file.path(project_root, "data", "GSE238264", "HCC7NR"), response="NonResponder")
)

log("定义signature基因...")
tpex_up    <- c("TCF7","SLAMF6","IL7R","CCR7","SELL","CD8A","CD8B")
tpex_down  <- c("GZMB","PRF1","HAVCR2","TIGIT")
cxcl13_up  <- c("CXCL13","CXCR3","LAG3","PDCD1","TOX","TIGIT","HAVCR2")
tls_up     <- c("CCL19","CCL21","CXCL13","LTB","TNFSF13B","CR2","SELL","CD19","MS4A1")
spp1_up    <- c("SPP1","CD163","TREM2","LGALS3","FN1","VEGFA")

log("开始逐样本处理...")
all_scores <- list()

for(sample_id in names(samples)) {
  log(paste("处理样本:", sample_id, samples[[sample_id]]$response))
  
  tryCatch({
    # 读取Visium数据
    obj <- Load10X_Spatial(
      data.dir = samples[[sample_id]]$path,
      filename = "filtered_feature_bc_matrix.h5"
    )
    
    # 归一化
    obj <- NormalizeData(obj)
    
    # 提取表达矩阵
    expr_mat <- as.matrix(GetAssayData(obj, layer="data"))
    
    # 计算signature scores
    rankData <- rankGenes(expr_mat)
    
    tpex_sc <- simpleScore(rankData,
                            upSet   = intersect(tpex_up,   rownames(expr_mat)),
                            downSet = intersect(tpex_down, rownames(expr_mat)))
    cxcl13_sc <- simpleScore(rankData,
                              upSet = intersect(cxcl13_up, rownames(expr_mat)))
    tls_sc <- simpleScore(rankData,
                           upSet = intersect(tls_up, rownames(expr_mat)))
    spp1_sc <- simpleScore(rankData,
                            upSet = intersect(spp1_up, rownames(expr_mat)))
    
    # 汇总
    score_df <- data.frame(
      spot        = rownames(tpex_sc),
      sample      = sample_id,
      response    = samples[[sample_id]]$response,
      tpex_score  = tpex_sc$TotalScore,
      cxcl13_score= cxcl13_sc$TotalScore,
      tls_score   = tls_sc$TotalScore,
      spp1_score  = spp1_sc$TotalScore
    )
    
    all_scores[[sample_id]] <- score_df
    log(paste("  完成，spot数:", nrow(score_df)))
    
  }, error = function(e) {
    log(paste("  错误:", e$message))
  })
}

log("合并所有样本...")
combined <- do.call(rbind, all_scores)
write.csv(combined, file.path(project_root, "results", "Q3_all_scores.csv"), row.names=FALSE)
log(paste("总spot数:", nrow(combined)))

log("计算各样本平均score...")
summary_df <- combined %>%
  group_by(sample, response) %>%
  summarise(
    n_spots      = n(),
    mean_tpex    = mean(tpex_score),
    mean_cxcl13  = mean(cxcl13_score),
    mean_tls     = mean(tls_score),
    mean_spp1    = mean(spp1_score),
    .groups = "drop"
  )
write.csv(summary_df, file.path(project_root, "results", "Q3_sample_summary.csv"), row.names=FALSE)
print(summary_df)

log("计算Responder vs NonResponder差异...")
for(sc in c("tpex_score","cxcl13_score","tls_score","spp1_score")) {
  r  <- combined[combined$response=="Responder",    sc]
  nr <- combined[combined$response=="NonResponder", sc]
  tt <- wilcox.test(r, nr)
  log(paste(sc, "R vs NR p:", round(tt$p.value, 4),
            "R mean:", round(mean(r),4),
            "NR mean:", round(mean(nr),4)))
}

log("计算spot级别相关...")
cor_tpex_cxcl13 <- cor(combined$tpex_score, combined$cxcl13_score, method="spearman")
cor_tpex_tls    <- cor(combined$tpex_score, combined$tls_score,    method="spearman")
cor_tpex_spp1   <- cor(combined$tpex_score, combined$spp1_score,   method="spearman")
log(paste("Tpex vs CXCL13 Spearman rho:", round(cor_tpex_cxcl13, 3)))
log(paste("Tpex vs TLS Spearman rho:",    round(cor_tpex_tls,    3)))
log(paste("Tpex vs SPP1 Spearman rho:",   round(cor_tpex_spp1,   3)))

log("Q3完成!")

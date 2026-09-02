# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

# Task 1: Rerun GSE202069 analysis with proper label audit
# Outputs are written under results/ in HCC_PROJECT_ROOT.

library(data.table)
library(singscore)

logf <- file(file.path(project_root, "results", "logs", "rerun_q4b_gse202069.log"), open="wt")
lg <- function(msg) { cat(paste0("[", Sys.time(), "] ", msg, "\n"), file=logf); cat(msg, "\n") }

lg("=== Task 1: GSE202069 Rerun ===")

# --- 1. Load expression data ---
lg("Loading expression data...")
expr <- fread(file.path(project_root, "data", "GSE202069", "GSE202069_tpm.txt.gz"), sep="\t")
gene_names <- expr[[1]]
expr_mat <- as.matrix(expr[, -1, with=FALSE])
rownames(expr_mat) <- gene_names
lg(paste("Genes:", nrow(expr_mat), "Samples:", ncol(expr_mat)))

# --- 2. Reproduce FAD inference (audit) ---
lg("Reproducing FAD inference for audit...")
fad_genes <- c("ACSL6","ACSL5","CPT1C","ADH7","ACADVL","ADH5","ACOX3",
               "ACSL3","ALDH1B1","ACSL1","ACAA2","ACAT1","ACADM","ACADL",
               "ACAT2","CPT1B","CPT1A","ACADS","GCDH","ACAA1","ALDH9A1",
               "HADHA","HADHB","CPT2","ALDH2","ALDH3A2","ACOX1","ECI1",
               "ECHS1","HADH","ACADSB","ECI2","ALDH7A1","ACSL4","CYP4A22",
               "CYP4A11","ADH1A","ADH1B","ADH6","EHHADH","ADH1C","ADH4")
fad_available <- intersect(fad_genes, rownames(expr_mat))
lg(paste("FAD genes available:", length(fad_available), "/", length(fad_genes)))

fad_score <- colMeans(expr_mat[fad_available, ], na.rm=TRUE)
pt_samples <- names(fad_score)[grepl("^PT", names(fad_score))]
pt_scores <- sort(fad_score[pt_samples])

tertile <- quantile(pt_scores, c(1/3, 2/3))
pt_subtype <- ifelse(pt_scores <= tertile[1], "F1",
               ifelse(pt_scores <= tertile[2], "F2", "F3"))
pt_response <- ifelse(pt_subtype == "F1", "FAD_inferred_Responder",
                ifelse(pt_subtype == "F3", "FAD_inferred_NonResponder", "Uncertain"))

label_audit <- data.frame(
  sample = names(pt_scores),
  FAD_score = round(pt_scores, 4),
  FAD_subtype = pt_subtype,
  inferred_group = pt_response,
  label_source = "FAD_tertile_inference",
  is_clinical_response = FALSE,
  notes = ifelse(pt_subtype == "F2",
    "Excluded: middle tertile, uncertain classification",
    "FAD-inferred subgroup, NOT clinical RECIST response")
)
write.csv(label_audit, file.path(project_root, "results", "tables", "GSE202069_label_audit.csv"), row.names=FALSE)
lg("Label audit saved")

# --- 3. Score calculation ---
lg("Computing scores...")
responders <- label_audit$sample[label_audit$FAD_subtype == "F1"]
nonresponders <- label_audit$sample[label_audit$FAD_subtype == "F3"]

tpex_up   <- c("TCF7","SLAMF6","IL7R","CCR7","SELL","CD8A","CD8B")
tpex_down <- c("GZMB","PRF1","HAVCR2","TIGIT")
cxcl13_up <- c("CXCL13","CXCR3","LAG3","PDCD1","TOX","TIGIT","HAVCR2")

tpex_up_avail   <- intersect(tpex_up, rownames(expr_mat))
tpex_down_avail <- intersect(tpex_down, rownames(expr_mat))
cxcl13_up_avail <- intersect(cxcl13_up, rownames(expr_mat))
lg(paste("Tpex up:", paste(tpex_up_avail, collapse=","), "| down:", paste(tpex_down_avail, collapse=",")))
lg(paste("CXCL13 up:", paste(cxcl13_up_avail, collapse=",")))

rankData <- rankGenes(expr_mat)
tpex_sc   <- simpleScore(rankData, upSet=tpex_up_avail, downSet=tpex_down_avail)
cxcl13_sc <- simpleScore(rankData, upSet=cxcl13_up_avail)

score_df <- data.frame(
  sample       = rownames(tpex_sc),
  tpex_score   = tpex_sc$TotalScore,
  cxcl13_score = cxcl13_sc$TotalScore
)
score_df$combined_score <- score_df$tpex_score + score_df$cxcl13_score
score_df$group <- ifelse(score_df$sample %in% responders, "FAD_Responder",
                  ifelse(score_df$sample %in% nonresponders, "FAD_NonResponder", NA))

ici_df <- score_df[!is.na(score_df$group), ]
lg(paste("FAD_Responder:", sum(ici_df$group=="FAD_Responder"),
         "FAD_NonResponder:", sum(ici_df$group=="FAD_NonResponder")))

# --- 4. Statistical tests with effect sizes ---
lg("Statistical tests...")
library(effsize)
library(coin)

run_test <- function(df, var, grp_var="group", grp1="FAD_Responder", grp2="FAD_NonResponder") {
  x1 <- df[[var]][df[[grp_var]] == grp1]
  x2 <- df[[var]][df[[grp_var]] == grp2]
  wt <- wilcox.test(x1, x2, conf.int=TRUE)
  cd <- cohen.d(x1, x2)
  n1 <- length(x1); n2 <- length(x2)
  r_effect <- abs(qnorm(wt$p.value/2)) / sqrt(n1 + n2)  # rank-biserial approx
  data.frame(
    variable = var,
    mean_responder = mean(x1),
    mean_nonresponder = mean(x2),
    median_responder = median(x1),
    median_nonresponder = median(x2),
    mann_whitney_p = wt$p.value,
    mann_whitney_CI_lower = wt$conf.int[1],
    mann_whitney_CI_upper = wt$conf.int[2],
    cohens_d = cd$estimate,
    cohens_d_CI_lower = cd$conf.int[1],
    cohens_d_CI_upper = cd$conf.int[2],
    n_responder = n1,
    n_nonresponder = n2,
    direction = ifelse(mean(x1) > mean(x2), "Responder > NonResponder", "NonResponder > Responder")
  )
}

stats_202 <- rbind(
  run_test(ici_df, "tpex_score"),
  run_test(ici_df, "cxcl13_score"),
  run_test(ici_df, "combined_score")
)
stats_202$cohort <- "GSE202069"
stats_202$label_type <- "FAD_inferred"

write.csv(ici_df, file.path(project_root, "results", "tables", "GSE202069_rerun_scores.csv"), row.names=FALSE)

# --- 5. Load GSE140901 for cross-cohort ---
lg("Loading GSE140901 scores...")
ici140 <- read.csv(file.path(project_root, "results", "tables", "GSE140901_scores_original_method.csv"))
ici140_filt <- ici140[ici140$response %in% c("Responder","NonResponder"), ]
ici140_filt$combined_score <- ici140_filt$tpex_score + ici140_filt$cxcl13_score
ici140_filt$group <- ici140_filt$response

stats_140 <- rbind(
  run_test(ici140_filt, "tpex_score", "group", "Responder", "NonResponder"),
  run_test(ici140_filt, "cxcl13_score", "group", "Responder", "NonResponder"),
  run_test(ici140_filt, "combined_score", "group", "Responder", "NonResponder")
)
stats_140$cohort <- "GSE140901"
stats_140$label_type <- "clinical_RECIST"

all_stats <- rbind(stats_202, stats_140)
write.csv(all_stats, file.path(project_root, "results", "stats", "Figure6_stats_summary.csv"), row.names=FALSE)

lg("=== KEY RESULTS ===")
for(i in 1:nrow(all_stats)) {
  s <- all_stats[i,]
  lg(paste(s$cohort, s$variable, ": d=", round(s$cohens_d, 3),
           "p=", round(s$mann_whitney_p, 4),
           "direction:", s$direction))
}

close(logf)
lg("Task 1 data complete. Proceeding to plot script.")

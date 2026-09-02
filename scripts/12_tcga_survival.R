# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

# Task 3: TCGA Cox modeling rerun
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(ggplot2)
library(singscore)

logf <- file(file.path(project_root, "results", "logs", "rerun_tcga_cox.log"), open="wt")
lg <- function(msg) { cat(paste0("[", Sys.time(), "] ", msg, "\n"), file=logf); cat(msg, "\n") }

lg("=== Task 3: TCGA Cox Modeling Rerun ===")

# --- 1. Load data (reproduce from scratch) ---
lg("Loading expression data...")
expr <- fread(file.path(project_root, "data", "TCGA_LIHC", "TCGA_LIHC_expression.gz"), sep="\t")
gene_names <- expr[[1]]
expr_mat <- as.matrix(expr[, -1, with=FALSE])
rownames(expr_mat) <- gene_names

lg("Loading survival data...")
surv <- read.table(file.path(project_root, "data", "TCGA_LIHC", "TCGA_LIHC_survival.txt"),
                   header=TRUE, sep="\t", stringsAsFactors=FALSE)

lg("Computing scores...")
tpex_up   <- c("TCF7","SLAMF6","IL7R","CCR7","SELL","CD8A","CD8B")
tpex_down <- c("GZMB","PRF1","HAVCR2","TIGIT")
cxcl13_up <- c("CXCL13","CXCR3","LAG3","PDCD1","TOX","TIGIT","HAVCR2")

rankData <- rankGenes(expr_mat)
tpex_sc   <- simpleScore(rankData, upSet=intersect(tpex_up, rownames(expr_mat)),
                          downSet=intersect(tpex_down, rownames(expr_mat)))
cxcl13_sc <- simpleScore(rankData, upSet=intersect(cxcl13_up, rownames(expr_mat)))

score_df <- data.frame(
  sample = colnames(expr_mat),
  tpex_score = tpex_sc$TotalScore,
  cxcl13_score = cxcl13_sc$TotalScore,
  combined_score = tpex_sc$TotalScore + cxcl13_sc$TotalScore
)

merged <- merge(score_df, surv, by.x="sample", by.y="sample")
lg(paste("Matched samples:", nrow(merged)))
merged <- merged[!is.na(merged$OS.time) & !is.na(merged$OS), ]
lg(paste("With OS data:", nrow(merged)))

# --- 2. Continuous univariate Cox (PRIMARY) ---
lg("=== PRIMARY: Continuous Univariate Cox ===")

cox_cont <- list()
cont_results <- data.frame()

for(var in c("tpex_score","cxcl13_score","combined_score")) {
  fml <- as.formula(paste("Surv(OS.time, OS) ~", var))
  cx <- coxph(fml, data=merged)
  s <- summary(cx)
  ci <- s$conf.int
  cox_cont[[var]] <- cx

  res <- data.frame(
    variable = var,
    analysis = "continuous_univariate",
    HR = ci[1,1],
    HR_lower = ci[1,3],
    HR_upper = ci[1,4],
    coef = s$coefficients[1,1],
    se = s$coefficients[1,3],
    z = s$coefficients[1,4],
    p = s$coefficients[1,5],
    n = cx$n,
    n_events = cx$nevent,
    concordance = s$concordance[1],
    ph_test_p = cox.zph(cx)$table[1,"p"]
  )
  cont_results <- rbind(cont_results, res)
  lg(paste(var, ": HR=", round(ci[1,1],3),
           " (", round(ci[1,3],3), "-", round(ci[1,4],3), ")",
           " p=", format.pval(s$coefficients[1,5], digits=3),
           " PH test p=", round(cox.zph(cx)$table[1,"p"], 3)))
}

write.csv(cont_results, file.path(project_root, "results", "tables", "TCGA_cox_results_continuous.csv"), row.names=FALSE)

# --- 3. Dichotomized Cox (SENSITIVITY) ---
lg("=== SENSITIVITY: Dichotomized Cox ===")

merged$tpex_grp <- ifelse(merged$tpex_score > median(merged$tpex_score), "High", "Low")
merged$cxcl13_grp <- ifelse(merged$cxcl13_score > median(merged$cxcl13_score), "High", "Low")
merged$combined_grp <- ifelse(merged$combined_score > median(merged$combined_score), "High", "Low")

dich_results <- data.frame()

for(grp_var in c("tpex_grp","cxcl13_grp","combined_grp")) {
  fml <- as.formula(paste("Surv(OS.time, OS) ~", grp_var))
  cx <- coxph(fml, data=merged)
  s <- summary(cx)
  ci <- s$conf.int

  res <- data.frame(
    variable = grp_var,
    analysis = "dichotomized_univariate_median_split",
    HR = ci[1,1],
    HR_lower = ci[1,3],
    HR_upper = ci[1,4],
    coef = s$coefficients[1,1],
    se = s$coefficients[1,3],
    z = s$coefficients[1,4],
    p = s$coefficients[1,5],
    n = cx$n,
    n_events = cx$nevent,
    concordance = s$concordance[1],
    ph_test_p = cox.zph(cx)$table[1,"p"]
  )
  dich_results <- rbind(dich_results, res)
  lg(paste(grp_var, ": HR=", round(ci[1,1],3),
           " (", round(ci[1,3],3), "-", round(ci[1,4],3), ")",
           " p=", format.pval(s$coefficients[1,5], digits=3)))
}

write.csv(dich_results, file.path(project_root, "results", "tables", "TCGA_cox_results_dichotomized.csv"), row.names=FALSE)

# --- 4. Multivariate sensitivity ---
lg("=== SENSITIVITY: Multivariate Cox ===")

clin <- read.delim(file.path(project_root, "data", "TCGA_LIHC", "LIHC_clinicalMatrix.txt"),
                   stringsAsFactors=FALSE)
clin$sample_id <- substr(clin$sampleID, 1, 15)
merged$sample_id <- substr(merged$sample, 1, 15)

clin_clean <- clin %>%
  select(sample_id, age_at_initial_pathologic_diagnosis, gender,
         pathologic_stage, neoplasm_histologic_grade) %>%
  mutate(
    age = as.numeric(age_at_initial_pathologic_diagnosis),
    sex = ifelse(gender == "MALE", 1, 0),
    stage = case_when(
      pathologic_stage %in% c("Stage I") ~ 1,
      pathologic_stage %in% c("Stage II") ~ 2,
      pathologic_stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC") ~ 3,
      pathologic_stage %in% c("Stage IV","Stage IVA","Stage IVB") ~ 4,
      TRUE ~ NA_real_
    ),
    grade = case_when(
      neoplasm_histologic_grade == "G1" ~ 1,
      neoplasm_histologic_grade == "G2" ~ 2,
      neoplasm_histologic_grade == "G3" ~ 3,
      neoplasm_histologic_grade == "G4" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  select(sample_id, age, sex, stage, grade)

df_mv <- merged %>% left_join(clin_clean, by="sample_id")
df_cc <- df_mv[complete.cases(df_mv[,c("OS","OS.time","tpex_score","age","sex","stage")]),]
lg(paste("Complete cases (age+sex+stage):", nrow(df_cc)))

df_cc_full <- df_mv[complete.cases(df_mv[,c("OS","OS.time","tpex_score","age","sex","stage","grade")]),]
lg(paste("Complete cases (age+sex+stage+grade):", nrow(df_cc_full)))

mv_results <- data.frame()

# Model 1: score + age + sex + stage
for(var in c("tpex_score","cxcl13_score","combined_score")) {
  fml <- as.formula(paste("Surv(OS.time, OS) ~", var, "+ age + sex + stage"))
  cx <- coxph(fml, data=df_cc)
  s <- summary(cx)
  ci <- s$conf.int[var,]
  pv <- s$coefficients[var, "Pr(>|z|)"]
  ph <- cox.zph(cx)$table[var, "p"]

  res <- data.frame(
    variable = var,
    model = "score + age + sex + stage",
    HR = ci[1], HR_lower = ci[3], HR_upper = ci[4],
    p = pv, n = cx$n, n_events = cx$nevent,
    ph_test_p = ph,
    concordance = s$concordance[1]
  )
  mv_results <- rbind(mv_results, res)
  lg(paste("MV1", var, ": HR=", round(ci[1],3),
           " (", round(ci[3],3), "-", round(ci[4],3), ")",
           " p=", format.pval(pv, digits=3),
           " n=", cx$n))
}

# Model 2: + grade
for(var in c("tpex_score","cxcl13_score","combined_score")) {
  fml <- as.formula(paste("Surv(OS.time, OS) ~", var, "+ age + sex + stage + grade"))
  cx <- coxph(fml, data=df_cc_full)
  s <- summary(cx)
  ci <- s$conf.int[var,]
  pv <- s$coefficients[var, "Pr(>|z|)"]
  ph <- cox.zph(cx)$table[var, "p"]

  res <- data.frame(
    variable = var,
    model = "score + age + sex + stage + grade",
    HR = ci[1], HR_lower = ci[3], HR_upper = ci[4],
    p = pv, n = cx$n, n_events = cx$nevent,
    ph_test_p = ph,
    concordance = s$concordance[1]
  )
  mv_results <- rbind(mv_results, res)
  lg(paste("MV2", var, ": HR=", round(ci[1],3),
           " (", round(ci[3],3), "-", round(ci[4],3), ")",
           " p=", format.pval(pv, digits=3),
           " n=", cx$n))
}

write.csv(mv_results, file.path(project_root, "results", "tables", "TCGA_multivariate_sensitivity.csv"), row.names=FALSE)

# --- 5. Plots ---
lg("Generating plots...")

theme_pub <- function() {
  theme_classic() +
  theme(
    text         = element_text(family="sans", color="black"),
    axis.text    = element_text(size=10, color="black"),
    axis.title   = element_text(size=11, color="black"),
    legend.text  = element_text(size=9),
    plot.title   = element_text(size=11, hjust=0.5),
    plot.subtitle= element_text(size=8, color="grey50", hjust=0.5, face="italic"),
    axis.line    = element_line(color="black", linewidth=0.5),
    axis.ticks   = element_line(color="black", linewidth=0.5)
  )
}

# Fig5D: Continuous univariate Cox forest
forest_data <- cont_results
forest_data$label <- c("Tpex-like Score", "CXCL13-associated exhaustion score", "Combined Score")
forest_data$label <- factor(forest_data$label, levels=rev(forest_data$label))
forest_data$sig <- ifelse(forest_data$p < 0.05, "p < 0.05", "ns")

p5D <- ggplot(forest_data, aes(x=HR, y=label)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=0.6) +
  geom_errorbar(aes(xmin=HR_lower, xmax=HR_upper), width=0.2, linewidth=0.8, color="#333333") +
  geom_point(aes(color=sig), size=4) +
  scale_color_manual(values=c("p < 0.05"="#E41A1C", "ns"="#999999")) +
  theme_pub() +
  labs(x="Hazard Ratio (95% CI)", y="",
       title="Univariate Cox Regression (TCGA-LIHC)",
       subtitle=paste0("Continuous scores, N=", cont_results$n[1], ", events=", cont_results$n_events[1])) +
  theme(legend.position="none") +
  geom_text(aes(label=paste0("HR=", round(HR,2), " p=", format.pval(p, digits=2))),
            hjust=-0.1, size=3, nudge_x=0.02)

tiff(file.path(project_root, "results", "panels", "Fig5D_TCGA_continuous_univariate_cox.tiff"),
     width=2600, height=1600, res=300, compression="lzw")
print(p5D); dev.off()

# FigS: Multivariate sensitivity forest
mv_plot <- mv_results
mv_plot$label <- paste0(gsub("_score","", mv_plot$variable), " | ", mv_plot$model)
mv_plot$label <- factor(mv_plot$label, levels=rev(mv_plot$label))
mv_plot$sig <- ifelse(mv_plot$p < 0.05, "p < 0.05", "ns")

pMV <- ggplot(mv_plot, aes(x=HR, y=label)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=0.6) +
  geom_errorbar(aes(xmin=HR_lower, xmax=HR_upper), width=0.2, linewidth=0.7, color="#333333") +
  geom_point(aes(color=sig), size=3.5) +
  scale_color_manual(values=c("p < 0.05"="#E41A1C", "ns"="#999999")) +
  theme_pub() +
  labs(x="Hazard Ratio (95% CI)", y="",
       title="Multivariate Cox Sensitivity (TCGA-LIHC)",
       subtitle="Adjusting for clinical covariates") +
  theme(legend.position="none", axis.text.y=element_text(size=8)) +
  geom_text(aes(label=paste0("p=", format.pval(p, digits=2))),
            hjust=-0.1, size=2.5, nudge_x=0.02)

tiff(file.path(project_root, "results", "panels", "FigS_TCGA_multivariate_sensitivity.tiff"),
     width=3200, height=2400, res=300, compression="lzw")
print(pMV); dev.off()

close(logf)
cat("Task 3 complete.\n")

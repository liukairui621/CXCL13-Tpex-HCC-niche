# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

# Task 4: Spatial validation at patient level
library(data.table)
library(ggplot2)
library(ggpubr)
library(effsize)
library(lme4)

logf <- file(file.path(project_root, "results", "logs", "rerun_spatial.log"), open="wt")
lg <- function(msg) { cat(paste0("[", Sys.time(), "] ", msg, "\n"), file=logf); cat(msg, "\n") }

lg("=== Task 4: Spatial Patient-Level Rerun ===")

# --- 1. Load spot-level data ---
lg("Loading spot-level data...")
spots <- fread(file.path(project_root, "results", "Q3_all_scores.csv"))
spots$combined_score <- spots$tpex_score + spots$cxcl13_score + spots$tls_score
lg(paste("Total spots:", nrow(spots), "Samples:", length(unique(spots$sample))))

# --- 2. Patient/sample-level summary ---
lg("Computing patient-level summaries...")
patient_summary <- spots[, .(
  n_spots = .N,
  mean_tpex = mean(tpex_score),
  mean_cxcl13 = mean(cxcl13_score),
  mean_tls = mean(tls_score),
  mean_spp1 = mean(spp1_score),
  combined_score = mean(tpex_score + cxcl13_score + tls_score),
  sd_tpex = sd(tpex_score),
  sd_cxcl13 = sd(cxcl13_score),
  sd_tls = sd(tls_score)
), by=.(sample, response)]

write.csv(patient_summary, file.path(project_root, "results", "tables", "spatial_sample_summary.csv"), row.names=FALSE)
lg("Patient summary:")
print(patient_summary)

# --- 3. Patient-level statistics ---
lg("Patient-level statistical tests...")

run_patient_test <- function(df, var) {
  r <- df[[var]][df$response == "Responder"]
  nr <- df[[var]][df$response == "NonResponder"]
  wt <- wilcox.test(r, nr, conf.int=TRUE, exact=FALSE)
  cd <- cohen.d(r, nr)
  data.frame(
    variable = var,
    level = "patient",
    mean_R = mean(r), mean_NR = mean(nr),
    median_R = median(r), median_NR = median(nr),
    diff_means = mean(r) - mean(nr),
    wilcox_p = wt$p.value,
    wilcox_CI_lower = wt$conf.int[1],
    wilcox_CI_upper = wt$conf.int[2],
    cohens_d = cd$estimate,
    d_CI_lower = cd$conf.int[1],
    d_CI_upper = cd$conf.int[2],
    n_R = length(r), n_NR = length(nr),
    direction = ifelse(mean(r) > mean(nr), "R > NR", "NR > R")
  )
}

patient_stats <- rbind(
  run_patient_test(patient_summary, "mean_tpex"),
  run_patient_test(patient_summary, "mean_cxcl13"),
  run_patient_test(patient_summary, "mean_tls"),
  run_patient_test(patient_summary, "mean_spp1"),
  run_patient_test(patient_summary, "combined_score")
)

write.csv(patient_stats, file.path(project_root, "results", "stats", "spatial_patient_level_stats.csv"), row.names=FALSE)

lg("=== PATIENT-LEVEL RESULTS ===")
for(i in 1:nrow(patient_stats)) {
  s <- patient_stats[i,]
  lg(paste(s$variable, ": d=", round(s$cohens_d, 3),
           " p=", round(s$wilcox_p, 4),
           " direction:", s$direction,
           " N_R=", s$n_R, " N_NR=", s$n_NR))
}

# --- 4. Spot-level descriptive stats (for comparison) ---
lg("Spot-level descriptive statistics...")

run_spot_test <- function(df, var) {
  r <- df[[var]][df$response == "Responder"]
  nr <- df[[var]][df$response == "NonResponder"]
  wt <- wilcox.test(r, nr, conf.int=TRUE)
  cd <- cohen.d(r, nr)
  data.frame(
    variable = var,
    level = "spot",
    mean_R = mean(r), mean_NR = mean(nr),
    median_R = median(r), median_NR = median(nr),
    diff_means = mean(r) - mean(nr),
    wilcox_p = wt$p.value,
    wilcox_CI_lower = wt$conf.int[1],
    wilcox_CI_upper = wt$conf.int[2],
    cohens_d = cd$estimate,
    d_CI_lower = cd$conf.int[1],
    d_CI_upper = cd$conf.int[2],
    n_R = length(r), n_NR = length(nr),
    direction = ifelse(mean(r) > mean(nr), "R > NR", "NR > R")
  )
}

spot_stats <- rbind(
  run_spot_test(spots, "tpex_score"),
  run_spot_test(spots, "cxcl13_score"),
  run_spot_test(spots, "tls_score"),
  run_spot_test(spots, "spp1_score"),
  run_spot_test(spots, "combined_score")
)

write.csv(spot_stats, file.path(project_root, "results", "stats", "spatial_spot_level_descriptive_stats.csv"), row.names=FALSE)

lg("=== SPOT-LEVEL RESULTS (descriptive only, pseudo-replication warning) ===")
for(i in 1:nrow(spot_stats)) {
  s <- spot_stats[i,]
  lg(paste(s$variable, ": d=", round(s$cohens_d, 3),
           " p=", format.pval(s$wilcox_p, digits=3),
           " direction:", s$direction,
           " N_R=", s$n_R, " N_NR=", s$n_NR,
           " *** PSEUDO-REPLICATION: N inflated by spot count ***"))
}

# --- 5. Mixed effects sensitivity (if lme4 available) ---
lg("Attempting mixed-effects sensitivity analysis...")
tryCatch({
  for(var in c("tpex_score","cxcl13_score","tls_score","combined_score")) {
    spots$response_bin <- ifelse(spots$response == "Responder", 1, 0)
    fml <- as.formula(paste(var, "~ response_bin + (1|sample)"))
    m <- lmer(fml, data=spots)
    s <- summary(m)
    coef_row <- s$coefficients["response_bin",]
    lg(paste("Mixed-effects", var, ": beta=", round(coef_row[1], 4),
             " SE=", round(coef_row[2], 4),
             " t=", round(coef_row[3], 3)))
  }
  lg("Mixed-effects analysis completed successfully.")
}, error = function(e) {
  lg(paste("Mixed-effects failed:", e$message))
  lg("Skipping mixed-effects. Patient-level Wilcoxon remains the appropriate test.")
})

# --- 6. Plots ---
lg("Generating patient-level plots...")

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

resp_colors <- c("Responder"="#F4A460","NonResponder"="#20B2AA")
patient_summary$response <- factor(patient_summary$response, levels=c("Responder","NonResponder"))

make_patient_box <- function(data, y_var, y_lab, title_str) {
  ggplot(data, aes_string(x="response", y=y_var, fill="response")) +
    geom_boxplot(width=0.45, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
    geom_jitter(width=0.1, size=3, alpha=0.9) +
    scale_fill_manual(values=resp_colors) +
    theme_pub() +
    labs(x="", y=y_lab, title=title_str,
         subtitle="Patient-level (n=4R vs 3NR) | GSE238264 Visium") +
    theme(legend.position="none", axis.text.x=element_text(size=11)) +
    stat_compare_means(method="wilcox.test", label="p.format",
                       label.x=1.5, label.y.npc=0.96, size=3.5,
                       exact=FALSE)
}

panels <- list(
  list(var="mean_tpex",     lab="Mean Tpex-like Score",    title="Tpex-like Signature",    fn="Fig4A_patient_level_tpex_boxplot.tiff"),
  list(var="mean_cxcl13",   lab="Mean CXCL13-associated exhaustion score", title="CXCL13-associated exhaustion signature", fn="Fig4B_patient_level_cxcl13_boxplot.tiff"),
  list(var="mean_tls",      lab="Mean TLS Score",          title="TLS Signature",          fn="Fig4C_patient_level_tls_boxplot.tiff"),
  list(var="combined_score", lab="Mean Combined Score",     title="Combined Score",         fn="Fig4D_patient_level_combined_boxplot.tiff")
)

for(pnl in panels) {
  p <- make_patient_box(patient_summary, pnl$var, pnl$lab, pnl$title)
  tiff(file.path(file.path(project_root, "results", "panels"), pnl$fn),
       width=1600, height=2000, res=300, compression="lzw")
  print(p); dev.off()
  cat(paste("Saved:", pnl$fn, "\n"))
}

close(logf)
cat("Task 4 complete.\n")

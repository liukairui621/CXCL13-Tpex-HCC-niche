suppressPackageStartupMessages({
  library(data.table)
  library(singscore)
  library(effsize)
  library(ggplot2)
})

project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork = FALSE)
base_dir <- file.path(project_root, "results", "external_validation")
for (subdir in c("tables", "stats", "panels", "logs")) {
  dir.create(file.path(base_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

log_file <- file(file.path(base_dir, "logs", "external_validation.log"), open = "wt")
log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message), file = log_file, append = TRUE)
  message(message)
}

tpex_up <- c("TCF7", "SLAMF6", "IL7R", "CCR7", "SELL", "CD8A", "CD8B")
tpex_down <- c("GZMB", "PRF1", "HAVCR2", "TIGIT")
cxcl13_genes <- c("CXCL13", "CXCR3", "LAG3", "PDCD1", "TOX", "TIGIT", "HAVCR2")

run_test <- function(scores, variable, dataset, treatment) {
  responder <- scores[[variable]][scores$response == "Responder"]
  nonresponder <- scores[[variable]][scores$response == "Nonresponder"]
  wilcox <- wilcox.test(responder, nonresponder, conf.int = TRUE, exact = FALSE)
  effect <- cohen.d(responder, nonresponder)
  data.frame(
    dataset = dataset,
    variable = variable,
    treatment = treatment,
    label_source = "GEO_clinical_metadata",
    mean_R = mean(responder),
    mean_NR = mean(nonresponder),
    median_R = median(responder),
    median_NR = median(nonresponder),
    wilcox_p = wilcox$p.value,
    wilcox_CI_lower = wilcox$conf.int[1],
    wilcox_CI_upper = wilcox$conf.int[2],
    cohens_d = effect$estimate,
    d_CI_lower = effect$conf.int[1],
    d_CI_upper = effect$conf.int[2],
    n_R = length(responder),
    n_NR = length(nonresponder),
    direction = ifelse(mean(responder) > mean(nonresponder), "R > NR", "NR > R")
  )
}

analyze_cohort <- function(expression_file, metadata_file, dataset, treatment) {
  expression <- fread(expression_file)
  genes <- expression[[1]]
  matrix <- as.matrix(expression[, -1, with = FALSE])
  rownames(matrix) <- genes
  metadata <- fread(metadata_file)

  ranked <- rankGenes(matrix)
  tpex <- simpleScore(
    ranked,
    upSet = intersect(tpex_up, rownames(matrix)),
    downSet = intersect(tpex_down, rownames(matrix))
  )
  cxcl13 <- simpleScore(ranked, upSet = intersect(cxcl13_genes, rownames(matrix)))
  scores <- data.frame(
    sample = rownames(tpex),
    tpex_score = tpex$TotalScore,
    cxcl13_associated_exhaustion_score = cxcl13$TotalScore,
    stringsAsFactors = FALSE
  )
  scores$response <- metadata$response[match(scores$sample, metadata$sample_id)]
  scores <- scores[!is.na(scores$response), ]
  scores$combined_score <- scores$tpex_score + scores$cxcl13_associated_exhaustion_score

  write.csv(
    scores,
    file.path(base_dir, "tables", paste0(dataset, "_signature_scores.csv")),
    row.names = FALSE
  )
  statistics <- rbind(
    run_test(scores, "tpex_score", dataset, treatment),
    run_test(scores, "cxcl13_associated_exhaustion_score", dataset, treatment),
    run_test(scores, "combined_score", dataset, treatment)
  )
  write.csv(
    statistics,
    file.path(base_dir, "stats", paste0(dataset, "_stats.csv")),
    row.names = FALSE
  )
  log_message(sprintf("%s: %d responders, %d nonresponders", dataset,
                      sum(scores$response == "Responder"),
                      sum(scores$response == "Nonresponder")))
  statistics
}

stats_215011 <- analyze_cohort(
  file.path(base_dir, "processed", "GSE215011_expression_matrix.tsv"),
  file.path(base_dir, "tables", "GSE215011_sample_metadata.csv"),
  "GSE215011",
  "Nivolumab (anti-PD-1)"
)
stats_279750 <- analyze_cohort(
  file.path(base_dir, "processed", "GSE279750_expression_matrix.tsv"),
  file.path(base_dir, "tables", "GSE279750_sample_metadata.csv"),
  "GSE279750",
  "anti-PD-L1"
)

stats_140901 <- read.csv(file.path(project_root, "results", "stats", "Figure6_stats_summary.csv"))
stats_140901 <- stats_140901[stats_140901$cohort == "GSE140901", ]
stats_140901$variable <- ifelse(
  stats_140901$variable == "cxcl13_score",
  "cxcl13_associated_exhaustion_score",
  stats_140901$variable
)
formatted_140901 <- data.frame(
  dataset = "GSE140901",
  variable = stats_140901$variable,
  treatment = "Nivolumab (anti-PD-1)",
  label_source = "GEO_clinical_RECIST",
  mean_R = stats_140901$mean_responder,
  mean_NR = stats_140901$mean_nonresponder,
  median_R = stats_140901$median_responder,
  median_NR = stats_140901$median_nonresponder,
  wilcox_p = stats_140901$mann_whitney_p,
  wilcox_CI_lower = NA_real_,
  wilcox_CI_upper = NA_real_,
  cohens_d = stats_140901$cohens_d,
  d_CI_lower = stats_140901$cohens_d_CI_lower,
  d_CI_upper = stats_140901$cohens_d_CI_upper,
  n_R = stats_140901$n_responder,
  n_NR = stats_140901$n_nonresponder,
  direction = stats_140901$direction
)

cross_cohort <- rbind(stats_215011, stats_279750, formatted_140901)
write.csv(cross_cohort, file.path(base_dir, "stats", "cross_cohort_effect_sizes.csv"), row.names = FALSE)

plot_forest <- function(data, variable, title, color, filename) {
  subset <- data[data$variable == variable, ]
  subset$label <- factor(subset$dataset, levels = rev(c("GSE140901", "GSE279750", "GSE215011")))
  plot <- ggplot(subset, aes(x = cohens_d, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = d_CI_lower, xmax = d_CI_upper), width = 0.2) +
    geom_point(size = 3.5, color = color) +
    theme_classic() +
    labs(x = "Cohen's d (responder vs nonresponder)", y = NULL, title = title)
  ggsave(file.path(base_dir, "panels", filename), plot, width = 7, height = 4.5, dpi = 600)
}

plot_forest(
  cross_cohort,
  "cxcl13_associated_exhaustion_score",
  "CXCL13-associated exhaustion score across clinical cohorts",
  "#E41A1C",
  "cross_cohort_CXCL13_associated_exhaustion.tiff"
)
plot_forest(
  cross_cohort,
  "tpex_score",
  "Tpex-like score across clinical cohorts",
  "#377EB8",
  "cross_cohort_Tpex.tiff"
)

close(log_file)

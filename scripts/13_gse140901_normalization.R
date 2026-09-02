# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

# Task 2: GSE140901 NanoString normalization rerun
library(GEOquery)
library(singscore)
library(effsize)
library(ggplot2)
library(ggpubr)

logf <- file(file.path(project_root, "results", "logs", "rerun_gse140901.log"), open="wt")
lg <- function(msg) { cat(paste0("[", Sys.time(), "] ", msg, "\n"), file=logf); cat(msg, "\n") }

lg("=== Task 2: GSE140901 NanoString Rerun ===")

# --- 1. Load metadata ---
lg("Loading GEO metadata...")
gse <- getGEO(filename=file.path(project_root, "data", "GSE140901", "GSE140901_series_matrix.txt.gz"))
pd <- pData(gse)
pd$response <- ifelse(pd[,"best_response:ch1"] == "PR", "Responder",
               ifelse(pd[,"best_response:ch1"] == "PD", "NonResponder", "SD"))

# --- 2. Parse ALL RCC files with full sections ---
lg("Parsing RCC files with controls...")
rcc_dir <- file.path(project_root, "data", "GSE140901", "rcc_files")
rcc_files <- list.files(rcc_dir, pattern="\\.RCC\\.gz$", full.names=TRUE)

parse_rcc_full <- function(filepath) {
  lines <- readLines(gzcon(file(filepath, "rb")))
  start <- grep("<Code_Summary>", lines) + 1
  end   <- grep("</Code_Summary>", lines) - 1
  data_lines <- lines[start:end]
  data_mat <- do.call(rbind, strsplit(data_lines, ","))
  df <- data.frame(
    CodeClass = data_mat[,1],
    Name  = data_mat[,2],
    Count = as.numeric(data_mat[,4]),
    stringsAsFactors = FALSE
  )
  return(df)
}

all_rcc <- lapply(rcc_files, parse_rcc_full)
sample_names <- sapply(strsplit(basename(rcc_files), "_"), `[`, 1)

# --- 3. Build matrices by class ---
build_matrix <- function(rcc_list, class_name, snames) {
  gene_list <- lapply(rcc_list, function(x) {
    d <- x[x$CodeClass == class_name, ]
    setNames(d$Count, d$Name)
  })
  all_genes <- unique(unlist(lapply(gene_list, names)))
  mat <- do.call(cbind, lapply(gene_list, function(x) {
    v <- rep(NA, length(all_genes)); names(v) <- all_genes; v[names(x)] <- x; v
  }))
  colnames(mat) <- snames
  return(mat)
}

endo_mat <- build_matrix(all_rcc, "Endogenous", sample_names)
hk_mat   <- build_matrix(all_rcc, "Housekeeping", sample_names)
pos_mat  <- build_matrix(all_rcc, "Positive", sample_names)
neg_mat  <- build_matrix(all_rcc, "Negative", sample_names)

lg(paste("Endogenous genes:", nrow(endo_mat), "Housekeeping:", nrow(hk_mat)))
lg(paste("Positive controls:", nrow(pos_mat), "Negative controls:", nrow(neg_mat)))

# --- 4A. Original method: log2(count+1) ---
lg("Method A: log2(count+1)...")
expr_orig <- log2(endo_mat + 1)

# --- 4B. NanoString normalization ---
lg("Method B: NanoString normalization...")

# Step 1: Positive control normalization
pos_geo_mean <- apply(pos_mat, 2, function(x) exp(mean(log(x[x > 0]))))
pos_norm_factor <- mean(pos_geo_mean) / pos_geo_mean
lg(paste("Positive control norm factors range:",
         round(min(pos_norm_factor), 3), "-", round(max(pos_norm_factor), 3)))

expr_pos_norm <- sweep(endo_mat, 2, pos_norm_factor, "*")

# Step 2: Negative control background subtraction
neg_means <- colMeans(neg_mat)
neg_sds   <- apply(neg_mat, 2, sd)
bg_threshold <- neg_means + 2 * neg_sds
lg(paste("Background thresholds range:",
         round(min(bg_threshold), 1), "-", round(max(bg_threshold), 1)))

expr_bg_sub <- sweep(expr_pos_norm, 2, bg_threshold, "-")
expr_bg_sub[expr_bg_sub < 1] <- 1

# Step 3: Housekeeping gene normalization
hk_pos_norm <- sweep(hk_mat, 2, pos_norm_factor, "*")
hk_geo_mean <- apply(hk_pos_norm, 2, function(x) exp(mean(log(x[x > 0]))))
hk_norm_factor <- mean(hk_geo_mean) / hk_geo_mean
lg(paste("HK norm factors range:",
         round(min(hk_norm_factor), 3), "-", round(max(hk_norm_factor), 3)))

expr_normalized <- sweep(expr_bg_sub, 2, hk_norm_factor, "*")
expr_normalized <- log2(expr_normalized + 1)

# --- 5. Score both methods ---
lg("Computing scores for both methods...")

tpex_up   <- c("TCF7","SLAMF6","IL7R","CCR7","SELL","CD8A","CD8B")
tpex_down <- c("GZMB","PRF1","HAVCR2","TIGIT")
cxcl13_up <- c("CXCL13","CXCR3","LAG3","PDCD1","TOX","TIGIT","HAVCR2")

compute_scores <- function(emat, method_name) {
  tu <- intersect(tpex_up, rownames(emat))
  td <- intersect(tpex_down, rownames(emat))
  cu <- intersect(cxcl13_up, rownames(emat))
  lg(paste(method_name, "- Tpex up:", paste(tu, collapse=","), "down:", paste(td, collapse=",")))
  lg(paste(method_name, "- CXCL13 up:", paste(cu, collapse=",")))

  rd <- rankGenes(emat)
  ts <- simpleScore(rd, upSet=tu, downSet=td)
  cs <- simpleScore(rd, upSet=cu)

  df <- data.frame(
    sample = rownames(ts),
    tpex_score = ts$TotalScore,
    cxcl13_score = cs$TotalScore,
    method = method_name
  )

  pd_m <- pd[match(df$sample, rownames(pd)), ]
  df$response <- pd_m$response
  df$best_response <- pd_m[, "best_response:ch1"]
  df$combined_score <- df$tpex_score + df$cxcl13_score
  return(df)
}

scores_orig <- compute_scores(expr_orig, "log2_count_plus1")
scores_norm <- compute_scores(expr_normalized, "NanoString_normalized")

write.csv(scores_orig, file.path(project_root, "results", "tables", "GSE140901_scores_original_method.csv"), row.names=FALSE)
write.csv(scores_norm, file.path(project_root, "results", "tables", "GSE140901_scores_normalized_method.csv"), row.names=FALSE)

# --- 6. Compare methods ---
lg("Comparing methods...")

run_test_140 <- function(df, var) {
  filt <- df[df$response %in% c("Responder","NonResponder"),]
  r <- filt[[var]][filt$response=="Responder"]
  nr <- filt[[var]][filt$response=="NonResponder"]
  wt <- wilcox.test(r, nr, conf.int=TRUE)
  cd <- cohen.d(r, nr)
  data.frame(
    variable=var, method=df$method[1],
    mean_R=mean(r), mean_NR=mean(nr),
    median_R=median(r), median_NR=median(nr),
    wilcox_p=wt$p.value,
    cohens_d=cd$estimate,
    d_CI_lower=cd$conf.int[1], d_CI_upper=cd$conf.int[2],
    n_R=length(r), n_NR=length(nr),
    direction=ifelse(mean(r) > mean(nr), "R>NR", "NR>R")
  )
}

comp <- rbind(
  run_test_140(scores_orig, "tpex_score"),
  run_test_140(scores_orig, "cxcl13_score"),
  run_test_140(scores_orig, "combined_score"),
  run_test_140(scores_norm, "tpex_score"),
  run_test_140(scores_norm, "cxcl13_score"),
  run_test_140(scores_norm, "combined_score")
)
write.csv(comp, file.path(project_root, "results", "stats", "GSE140901_method_comparison.csv"), row.names=FALSE)

lg("=== METHOD COMPARISON ===")
for(i in 1:nrow(comp)) {
  s <- comp[i,]
  lg(paste(s$method, s$variable, "d=", round(s$cohens_d,3), "p=", round(s$wilcox_p,4), s$direction))
}

# --- 7. Plots ---
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

resp_colors <- c("Responder"="#F4A460","NonResponder"="#20B2AA")

# FigS6C: GSE140901 CXCL13 boxplot (normalized method)
filt_norm <- scores_norm[scores_norm$response %in% c("Responder","NonResponder"),]
filt_norm$response <- factor(filt_norm$response, levels=c("Responder","NonResponder"))

pS6C <- ggplot(filt_norm, aes(x=response, y=cxcl13_score, fill=response)) +
  geom_boxplot(width=0.45, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.12, size=2.5, alpha=0.8) +
  scale_fill_manual(values=resp_colors) +
  theme_pub() +
  labs(x="", y="CXCL13-associated exhaustion score",
       title="GSE140901: CXCL13-associated exhaustion (NanoString Normalized)",
       subtitle="PR=6, PD=8 | Clinical RECIST labels") +
  theme(legend.position="none", axis.text.x=element_text(size=11)) +
  stat_compare_means(method="wilcox.test", label="p.format",
                     label.x=1.5, label.y.npc=0.96, size=3.5)

tiff(file.path(project_root, "results", "panels", "FigS6C_GSE140901_cxcl13_boxplot.tiff"),
     width=1800, height=2200, res=300, compression="lzw")
print(pS6C); dev.off()

# FigS6E: Sensitivity comparison
both_filt <- rbind(
  scores_orig[scores_orig$response %in% c("Responder","NonResponder"),],
  scores_norm[scores_norm$response %in% c("Responder","NonResponder"),]
)
both_filt$response <- factor(both_filt$response, levels=c("Responder","NonResponder"))
both_filt$method <- factor(both_filt$method,
  levels=c("log2_count_plus1","NanoString_normalized"),
  labels=c("log2(count+1)","NanoString Norm"))

pS6E <- ggplot(both_filt, aes(x=response, y=cxcl13_score, fill=response)) +
  geom_boxplot(width=0.45, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.12, size=2, alpha=0.7) +
  facet_wrap(~method) +
  scale_fill_manual(values=resp_colors) +
  theme_pub() +
  labs(x="", y="CXCL13-associated exhaustion score",
       title="Normalization Sensitivity: GSE140901",
       subtitle="Comparing log2(count+1) vs standard NanoString normalization") +
  theme(legend.position="none", strip.text=element_text(size=10, face="bold")) +
  stat_compare_means(method="wilcox.test", label="p.format",
                     label.x=1.5, label.y.npc=0.96, size=3)

tiff(file.path(project_root, "results", "panels", "FigS6E_GSE140901_normalization_sensitivity.tiff"),
     width=3200, height=2200, res=300, compression="lzw")
print(pS6E); dev.off()

close(logf)
cat("Task 2 complete.\n")

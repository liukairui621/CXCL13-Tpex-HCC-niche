# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(monocle3)
library(SeuratWrappers)

theme_pub <- function() {
  theme_classic() +
  theme(
    text         = element_text(family="sans", color="black"),
    axis.text    = element_text(size=10, color="black"),
    axis.title   = element_text(size=11, color="black"),
    legend.text  = element_text(size=9),
    legend.title = element_text(size=10),
    plot.title   = element_text(size=11, hjust=0.5),
    axis.line    = element_line(color="black", linewidth=0.5),
    axis.ticks   = element_line(color="black", linewidth=0.5)
  )
}

cellstate_colors <- c(
  "Tn_Tcm"="steelblue","Tcm"="#63C5C8","Tpex_like"="#E41A1C",
  "Early_Tex"="#F4A460","Terminal_Tex"="#8B4513","Cycling_T"="#4DAF4A",
  "Treg"="#984EA3","Tumor_T1"="#BDBDBD","Tumor_T2"="#969696",
  "Tumor_T3"="#737373","Tumor_T4"="#525252","Tumor_T5"="#252525",
  "Tumor_T6"="#D9D9D9","Tumor_T7"="#F0F0F0","NK"="#E8601C","NK_dim"="#FDAE61"
)

cat("Loading data...\n")
tnk <- readRDS(file.path(project_root, "results", "objects", "tnk_annotated.rds"))
cd8_states <- c("Tn_Tcm","Tcm","Tpex_like","Early_Tex","Terminal_Tex","Cycling_T")
cd8 <- subset(tnk, subset=cellstate %in% cd8_states)

cat("Building pseudotime...\n")
cds <- as.cell_data_set(cd8)
cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds)
root_cells <- colnames(cd8)[cd8$cellstate=="Tn_Tcm"]
cds <- order_cells(cds, root_cells=root_cells[1:min(10,length(root_cells))])

cat("Drawing Fig2A...\n")
p2A <- plot_cells(cds, color_cells_by="cellstate",
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE, label_branch_points=FALSE,
                  cell_size=0.6, group_label_size=3.5) +
  scale_color_manual(values=cellstate_colors) +
  theme_pub() + ggtitle("CD8\u207a T Cell Trajectory")
tiff(file.path(project_root, "results", "figures", "Fig2A.tiff"), width=2800, height=2400,
     res=300, compression="lzw")
print(p2A)
dev.off()
cat("Fig2A done\n")

cat("Drawing Fig2B...\n")
p2B <- plot_cells(cds, color_cells_by="pseudotime",
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE, label_branch_points=FALSE,
                  cell_size=0.6) +
  scale_color_gradientn(colors=c("#2166AC","#92C5DE","#F7F7F7",
                                  "#F4A582","#D6604D","#B2182B")) +
  theme_pub() + ggtitle("Pseudotime")
tiff(file.path(project_root, "results", "figures", "Fig2B.tiff"), width=2800, height=2400,
     res=300, compression="lzw")
print(p2B)
dev.off()
cat("Fig2B done\n")

cat("Drawing Fig2C...\n")
genes2C <- intersect(c("TCF7","SLAMF6","IL7R","PDCD1","TOX","GZMB","PRF1"),
                     rownames(tnk))
pt_vals <- pseudotime(cds)
pt_vals <- pt_vals[is.finite(pt_vals)]
expr_mat <- GetAssayData(cd8, layer="data")[genes2C, names(pt_vals)]
pt_df <- data.frame(pseudotime=pt_vals)
for(g in genes2C) pt_df[[g]] <- as.numeric(expr_mat[g,])
pt_long <- tidyr::pivot_longer(pt_df, cols=all_of(genes2C),
                                names_to="gene", values_to="expression")
gene_colors <- c("TCF7"="#2196A6","SLAMF6"="#63C5C8","IL7R"="steelblue",
                  "PDCD1"="#E41A1C","TOX"="#F4A460","GZMB"="#8B4513","PRF1"="#984EA3")
p2C <- ggplot(pt_long, aes(x=pseudotime, y=expression, color=gene)) +
  geom_smooth(method="loess", se=FALSE, span=0.5, linewidth=1) +
  scale_color_manual(values=gene_colors) +
  theme_pub() +
  labs(x="Pseudotime", y="Expression", color="Gene")
tiff(file.path(project_root, "results", "figures", "Fig2C.tiff"), width=2800, height=2000,
     res=300, compression="lzw")
print(p2C)
dev.off()
cat("Fig2C done\n")

cat("Drawing Fig2D...\n")
meta <- tnk@meta.data
meta_cd8 <- meta[meta$cellstate %in% cd8_states,]
meta_cd8$cellstate <- factor(meta_cd8$cellstate, levels=cd8_states)

p2D_stem <- ggplot(meta_cd8, aes(x=cellstate, y=stemness_score1, fill=cellstate)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3, linewidth=0.4) +
  scale_fill_manual(values=cellstate_colors) +
  theme_pub() + labs(x="", y="Stemness Score", title="Stemness") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1,size=9))

p2D_cyto <- ggplot(meta_cd8, aes(x=cellstate, y=cytotox_score1, fill=cellstate)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3, linewidth=0.4) +
  scale_fill_manual(values=cellstate_colors) +
  theme_pub() + labs(x="", y="Cytotoxicity Score", title="Cytotoxicity") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1,size=9))

p2D <- plot_grid(p2D_stem, p2D_cyto, ncol=2)
tiff(file.path(project_root, "results", "figures", "Fig2D.tiff"), width=3600, height=2000,
     res=300, compression="lzw")
print(p2D)
dev.off()
cat("Fig2D done\nFigure 2 complete\n")

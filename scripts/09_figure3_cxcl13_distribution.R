# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

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
expr_data <- GetAssayData(tnk, layer="data")
cxcl13_expr <- as.numeric(expr_data["CXCL13",])
meta <- tnk@meta.data
meta$CXCL13_expr <- cxcl13_expr
meta$CXCL13_pos  <- cxcl13_expr > 0

cat("Drawing Fig3A - CXCL13 dot plot...\n")
seurat_obj <- readRDS(file.path(project_root, "results", "objects", "seurat_clustered.rds"))
p3A <- DotPlot(seurat_obj, features="CXCL13", group.by="celltype",
               dot.scale=8, cols=c("#E8F4F8","#E41A1C")) +
  theme_pub() +
  labs(x="", y="Cell Type") +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))
tiff(file.path(project_root, "results", "figures", "Fig3A.tiff"), width=2000, height=2000,
     res=300, compression="lzw")
print(p3A)
dev.off()
cat("Fig3A done\n")

cat("Drawing Fig3B - CXCL13 violin...\n")
cellstate_order <- c("Treg","Cycling_T","Tpex_like","Early_Tex","Terminal_Tex",
                     "Tn_Tcm","Tcm","Tumor_T1","Tumor_T2","Tumor_T3",
                     "Tumor_T4","Tumor_T5","Tumor_T6","Tumor_T7","NK","NK_dim")
meta$cellstate <- factor(meta$cellstate, levels=cellstate_order)

p3B <- ggplot(meta, aes(x=cellstate, y=CXCL13_expr, fill=cellstate)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.08, fill="white", outlier.size=0.2, linewidth=0.4) +
  scale_fill_manual(values=cellstate_colors) +
  theme_pub() +
  labs(x="", y="CXCL13 Expression") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1,size=8))
tiff(file.path(project_root, "results", "figures", "Fig3B.tiff"), width=3600, height=2000,
     res=300, compression="lzw")
print(p3B)
dev.off()
cat("Fig3B done\n")

cat("Drawing Fig3C - CXCL13+ proportion...\n")
prop_df <- meta %>%
  group_by(cellstate) %>%
  summarise(pct_pos=mean(CXCL13_pos)*100, .groups="drop") %>%
  arrange(desc(pct_pos))
prop_df$cellstate <- factor(prop_df$cellstate, levels=prop_df$cellstate)

p3C <- ggplot(prop_df, aes(x=cellstate, y=pct_pos, fill=cellstate)) +
  geom_bar(stat="identity", color="black", linewidth=0.3, width=0.75) +
  scale_fill_manual(values=cellstate_colors) +
  theme_pub() +
  labs(x="", y="CXCL13\u207a cells (%)") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1,size=8))
tiff(file.path(project_root, "results", "figures", "Fig3C.tiff"), width=3000, height=2000,
     res=300, compression="lzw")
print(p3C)
dev.off()
cat("Fig3C done\n")

cat("Drawing Fig3D - descriptive Treg/Cycling_T comparison...\n")
meta$group <- ifelse(meta$cellstate %in% c("Treg","Cycling_T"),
                     as.character(meta$cellstate), "Other T/NK")
meta$group <- factor(meta$group, levels=c("Treg","Cycling_T","Other T/NK"))

p3D <- ggplot(meta, aes(x=group, y=CXCL13_expr, fill=group)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.2, linewidth=0.4) +
  scale_fill_manual(values=c("Treg"="#984EA3","Cycling_T"="#4DAF4A",
                              "Other T/NK"="#BDBDBD")) +
  theme_pub() +
  labs(x="", y="CXCL13 Expression") +
  theme(legend.position="none",
        axis.text.x=element_text(size=11))
tiff(file.path(project_root, "results", "figures", "Fig3D.tiff"), width=2000, height=2400,
     res=300, compression="lzw")
print(p3D)
dev.off()
ggsave(file.path(project_root, "results", "figures", "Fig3D.pdf"), p3D,
       width=5, height=5, device=cairo_pdf, bg="white")
cat("Fig3D done\nFigure 3 complete\n")

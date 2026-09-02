# Public reproduction configuration. Set HCC_PROJECT_ROOT to the repository root.
project_root <- normalizePath(Sys.getenv("HCC_PROJECT_ROOT", "."), mustWork=FALSE)
for (subdir in c("", "objects", "figures", "logs", "tables", "stats", "panels")) {
  dir.create(file.path(project_root, "results", subdir), recursive=TRUE, showWarnings=FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)

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

cat("Drawing Fig1A...\n")
p1A <- DimPlot(tnk, group.by="cellstate", label=TRUE, label.size=3,
               repel=TRUE, cols=cellstate_colors, pt.size=0.3) +
  ggtitle("T/NK Cell Subpopulations") + theme_pub()
tiff(file.path(project_root, "results", "figures", "Fig1A.tiff"), width=3000, height=2400,
     res=300, compression="lzw")
print(p1A)
dev.off()
cat("Fig1A done\n")

cat("Drawing Fig1B...\n")
markers <- c("TCF7","SLAMF6","IL7R","SELL","CCR7",
             "PDCD1","TOX","HAVCR2","LAG3","TIGIT",
             "GZMB","PRF1","IFNG","FOXP3","IL2RA","NKG7","GNLY")
p1B <- DotPlot(tnk, features=markers, group.by="cellstate",
               dot.scale=5, cols=c("#E8F4F8","#E41A1C")) +
  coord_flip() + theme_pub() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        axis.title=element_blank())
tiff(file.path(project_root, "results", "figures", "Fig1B.tiff"), width=3600, height=2400,
     res=300, compression="lzw")
print(p1B)
dev.off()
cat("Fig1B done\n")

cat("Drawing Fig1C...\n")
genes1C <- c("TCF7","PDCD1","SLAMF6","TOX","IL7R")
plist <- FeaturePlot(tnk, features=genes1C, ncol=3,
                     order=TRUE, combine=FALSE,
                     cols=c("#F0F0F0","#E41A1C"))
plist <- lapply(plist, function(p) {
  p + theme_pub() +
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          title=element_text(size=10, face="italic"))
})
p1C <- plot_grid(plotlist=plist, ncol=3)
tiff(file.path(project_root, "results", "figures", "Fig1C.tiff"), width=3600, height=2400,
     res=300, compression="lzw")
print(p1C)
dev.off()
cat("Fig1C done\n")

cat("Drawing Fig1D...\n")
meta <- tnk@meta.data
tpex_pt <- meta %>%
  group_by(patient, cellstate) %>% summarise(n=n(), .groups="drop") %>%
  group_by(patient) %>% mutate(pct=n/sum(n)*100) %>%
  filter(cellstate=="Tpex_like")
pt_colors <- brewer.pal(10,"Set3")
p1D <- ggplot(tpex_pt, aes(x=reorder(patient,-pct), y=pct, fill=patient)) +
  geom_bar(stat="identity", color="black", linewidth=0.3, width=0.7) +
  scale_fill_manual(values=pt_colors) +
  geom_hline(yintercept=mean(tpex_pt$pct), linetype="dashed",
             color="#E41A1C", linewidth=0.7) +
  theme_pub() +
  labs(x="Patient", y="Tpex-like cells (% of T/NK)") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1))
tiff(file.path(project_root, "results", "figures", "Fig1D.tiff"), width=2400, height=2000,
     res=300, compression="lzw")
print(p1D)
dev.off()
cat("Fig1D done\nFigure 1 complete\n")

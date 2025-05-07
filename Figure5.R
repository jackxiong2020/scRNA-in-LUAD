## library packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(CellChat)
library(survival)
library(AUCell)
library(RColorBrewer)
## load data
Tumor_cell = readRDS("~/data/Tumor_cell.rds")
############################## > Fig5a  ##################################
umap_plot <- UMAPPlot(Tumor_cell, group.by = "subtype", label = TRUE, label.size = 6) +
  NoLegend() +
  ggtitle("") +
  scale_color_manual(values = pal_d3("category20")(20)[c(18, 2, 18, 18)]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans")
  )

ggsave("~/figure5/Fig5a.pdf", umap_plot, width = 3.2, 
       height = 3.2, device = "pdf", bg = "transparent")

############################## >Fig5b  #################################
# Generate LAMP3 expression violin plot
p <- VlnPlot(
  Tumor_cell,
  features = "KLF4",
  group.by = "subtype",
  cols = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")  # D3 palette [4,2,3,1]
) + 
  labs(title = "KLF4", x = "", y = "Expression level") +
  theme(
    plot.title = element_text(face = "italic", size = 14),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid = element_blank()
  )

# Save plot
ggsave("~/figure5/Fig5_KLF4.pdf", p, width = 3.76, height = 3.45, bg = "transparent")

############################## >Fig5c   ##################################
# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "C2", group.by = "seurat_clusters") +
  labs(x = "", y = "Expression level") +
  NoLegend() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "italic"),
    panel.grid = element_blank()
  )
ggsave("Fig5c_C2.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "KLF4", group.by = "seurat_clusters") +
  labs(x = "", y = "Expression level") +
  NoLegend() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "italic"),
    panel.grid = element_blank()
  )

ggsave("Fig5c_KLF4.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

############################## > Fig5e  ##################################
# Calculate AUC scores and prepare data
score <- AUCell_run(Tumor_cell@assays$RNA@counts, genelist)@assays@data$AUC
data <- data.frame(t(score)) %>% 
  mutate(cluster = Tumor_cell$subtype) %>%
  select(cluster, AT1.cell, AT2.cell) %>%
  gather(celltype, AUCell, -cluster) %>%
  mutate(celltype = gsub("\\.", " ", celltype),
         cluster = factor(cluster, levels = c("C2", "C3", "C4", "C1")))

# Filter AT1 cells and set colors
data_AT1 <- data %>% 
  filter(celltype == "AT1 cell", AUCell < 0.55)
colors <- pal_d3("category20")(20)[c(2, 18, 18, 18)]

# Create and save boxplot
p <- ggpubr::ggboxplot(
  data_AT1, x = "cluster", y = "AUCell",
  fill = "cluster", color = "black",
  palette = colors, width = 0.4,
  add = "jitter", add.params = list(size = 0.4, alpha = 0.5)
) +
  ggtitle("AT1 score") +
  ylab("AUCell score") + xlab("") +
  ylim(0, 0.75) +
  scale_fill_manual(values = colors) +
  stat_compare_means(
    comparisons = list(c("C4", "C1"), c("C3", "C4"), c("C3", "C2")),
    method = "wilcox", label = "p.signif", size = 5
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14),
    legend.position = "right"
  )

ggsave("~/figure5/Fig5e.pdf", p, width = 2.5, height = 3.45, bg = "transparent")

############################## > Fig5f  ##################################
# Set cell identities and find markers
Tumor_cell$subtype <- factor(Tumor_cell$subtype, levels = c("C1","C2","C3","C4"))
Idents(Tumor_cell) <- "subtype"
Deg <- FindAllMarkers(Tumor_cell, logfc.threshold = 0.25) %>% 
  filter(avg_log2FC > 0.5, cluster == "C2")

# Perform GO enrichment
Go_bp <- enrichGO(
  gene = Deg$gene,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",
  ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05
)

# Prepare top 10 results
result <- Go_bp@result %>% 
  mutate(rank = -log10(p.adjust)) %>% 
  arrange(p.adjust) %>% 
  head(10)

# Create and save barplot
p <- ggplot(result, aes(x = reorder(Description, Count), y = Count, fill = rank)) +
  geom_col(alpha = 0.7, width = 0.8) +
  scale_fill_gradient(
    low = "orange", high = "#FF0000",breaks = c(6.5, 7.5, 8.5),name = "-log10(p adjust)"
  ) +
  coord_flip() +
  labs(
    x = "",y = "Gene number",title = "Enrichment GO Terms for C2"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.6, fill = NA),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    plot.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave("fig5f.pdf", p, width = 7.78, height = 4, bg = "transparent")

############################## >Fig5g  #################################
# Load and prepare gene sets
geneset <- read.gmt("~/data/h.all.v7.5.1.symbols.gmt") %>%
  mutate(term = gsub("HALLMARK_", "", term)) %>%
  split(.$term) %>%
  lapply("[[", 2)

# Calculate AUC scores
AUC.df <- AUCell_run(Tumor_cell@assays$RNA@counts, geneset)@assays@data$AUC %>%
  t() %>%
  data.frame() %>%
  mutate(cluster = Tumor_cell$subtype)

# Prepare P53 pathway data
data1 <- AUC.df %>%
  gather(celltype, AUCell, -cluster) %>%
  mutate(celltype = gsub("\\.", " ", celltype)) %>%
  filter(celltype == "P53 PATHWAY") %>%
  mutate(cluster = factor(cluster, levels = c("C2", "C3", "C4", "C1")))

# Create and save boxplot
p <- ggpubr::ggboxplot(
  data1, x = "cluster", y = "AUCell",
  fill = "cluster", color = "black",
  palette = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), # D3 colors [2,3,1,4]
  width = 0.4, add = "jitter",
  add.params = list(size = 0.4, alpha = 0.5)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c("C2","C1"), c("C2","C4"), c("C2","C3")),
    method = "wilcox", label = "p.signif", size = 5, vjust = 0.3
  ) +
  labs(
    title = "P53 PATHWAY",
    y = "AUCell score", 
    x = ""
  ) +
  ylim(0, 0.3) +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14)
  )

ggsave("~/figure5/Fig5g.pdf", p, width = 2.5, height = 3.45, bg = "transparent")

############################## > Fig5h  ##################################
# Perform cell cycle scoring
Tumor_cell <- CellCycleScoring(
  Tumor_cell,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

# Prepare data for plotting
Prop_data <- Tumor_cell@meta.data %>%
  dplyr::select(subtype, Phase)

# Create and save bar plot
p <- ggstatsplot::ggbarstats(
  data = Prop_data,x = Phase,y = subtype,
  results.subtitle = FALSE,bf.message = FALSE,
  proportion.test = FALSE,label.args = list(size = 3, fill = "white", alpha = 0.85),
  perc.k = 2,title = "",xlab = "",ylab = "Percentage",legend.title = "Phase"
) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C")) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

# Remove text layers if needed (uncomment if necessary)
 p <- delete_layers(p, match_type = "GeomText")

ggsave("~/figure5/Fig5h.pdf", p, width = 3.5, height = 3.6, bg = "transparent")

############################## > Fig5i  ##################################
# Load and prepare tumor suppressor genes
TSG <- read.table("~/data/LUAD_regulated_TSgenes.txt", header = FALSE)
colnames(TSG) <- c("Entrez", "Gene")

# Find marker genes
Idents(tumor) <- "subtype"
Deg <- FindAllMarkers(tumor, logfc.threshold = 0.25)

# Filter and combine significant markers
Deg1 <- Deg %>% 
  filter(avg_log2FC > 0.5) %>% 
  filter(cluster %in% c("C1", "C2", "C3", "C4"))

# Identify TSGs in each cluster
TSG_merge <- lapply(paste0("C", 1:4), function(x) {
  intersect(TSG$Gene, filter(Deg, cluster == x)$gene)
}) %>% 
  unlist() %>% 
  unique() %>% 
  sort()

# Create heatmap
heatmap_data <- AverageExpression(
  tumor,
  features = TSG_merge,
  slot = "data",
  group.by = "subtype"
)$SCT[, c("C2", "C1", "C3", "C4")]

# Define color scale
bk <- c(seq(-1, -0.1, 0.01), seq(0, 2, 0.01))
color_palette <- c(
  colorRampPalette(c("#2166ac", "#f7fbff"))(length(bk)/2),
  colorRampPalette(c("#f7fbff", "#b2182b"))(length(bk)/2)
)

# Generate and save heatmap
pdf("~/figure5/Fig5i.pdf", width = 3.38, height = 5.2)
ComplexHeatmap::pheatmap(
  heatmap_data,scale = "row",cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames = TRUE,show_rownames = TRUE,border = TRUE,
  color = color_palette,breaks = bk,name = "Expression",
  heatmap_legend_param = list(
    at = seq(-1, 2, 1),
    labels = seq(-1, 2, 1)
  ),
  row_names_gp = gpar(fontsize = 7, fontface = "italic"),
  column_names_gp = gpar(fontsize = 12, fontface = "plain")
)
dev.off()



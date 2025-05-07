## library packages
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(AUCell)
library(GSVA)
library(irGSEA)
library(ggpubr)
############################## > Fig2a  ##################################
# Load and prepare merged data
merge_data <- readRDS("~/data/merge_MIA_IA_Advanced_Metastatic_tumor_AT2.rds")

# Filter and relabel stages
merge_data$stage <- factor(
  merge_data$stage,
  levels = c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic"),
  labels = c(
    "AT2 cell (5,498)", "MIA (1,190)", "IA (8,463)", 
    "Advanced (5,298)", "Metastatic (12,141)"
  )
)

# Define color palette
stage_colors <- c("#479D88", "#6CB8D2", "#F08080", "#89558D", "#D1352B")

# Create UMAP plot
stage_umap <- UMAPPlot(
  merge_data,
  group.by = "stage",
  label = FALSE,
  cols = stage_colors
) +
  ggtitle("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(size = 14, family = "sans",face = "plain",color = "black"),
    axis.title = element_text(size = 16, family = "sans", face = "plain"),
    plot.title = element_text(size = 16, family = "sans", face = "plain"),
    legend.text = element_text(size = 14, family = "sans", face = "plain"),
    legend.title = element_text(size = 14, family = "sans", face = "plain")
  )

# Save plot
output_dir <- "～/Figure2/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig2a.pdf"),
  plot = stage_umap,
  width = 5.86,
  height = 3.73
)

############################## > Fig2b  ##################################
# Define Hallmark gene set categories
hallmark_categories <- list(
  "Death and stress" = sort(c("APOPTOSIS", "P53 PATHWAY", "HYPOXIA")),
  "Differentiation" = sort(c("EPITHELIAL MESENCHYMAL TRANSITION", "NOTCH SIGNALING", 
                             "WNT BETA CATENIN SIGNALING", "TGF BETA SIGNALING")),
  "Immune" = sort(c("IL6 JAK STAT3 SIGNALING", "INFLAMMATORY RESPONSE", "TNFA SIGNALING VIA NFKB")),
  "Metabolism" = sort(c("OXIDATIVE PHOSPHORYLATION", "GLYCOLYSIS", "FATTY ACID METABOLISM")),
  "Proliferation" = sort(c("MYC TARGETS V1", "MYC TARGETS V2", "E2F TARGETS", "MITOTIC SPINDLE"))
)

hallmark_use <- unlist(hallmark_categories, use.names = FALSE)

# Calculate pathway scores using AUCell
gsea_result <- irGSEA.score(
  object = merge_data,assay = "SCT",slot = "data",
  seeds = 123,ncores = 36,min.cells = 3,
  min.feature = 0,msigdb = TRUE,
  species = "Homo sapiens",category = "H",
  geneid = "symbol",method = "AUCell"
)

# Perform differential analysis
diff_result <- irGSEA.integrate(
  object = gsea_result,
  group.by = "stage",
  method = "AUCell"
)

# Process p-values
pval_matrix <- diff_result$AUCell %>%
  dplyr::select(Name, p_val_adj, cluster) %>%
  tidyr::spread(key = "cluster", value = "p_val_adj") %>%
  mutate(Name = gsub("HALLMARK |-", " ", Name)) %>%
  filter(Name %in% hallmark_use) %>%
  arrange(match(Name, hallmark_use)) %>%
  column_to_rownames("Name") %>%
  dplyr::select("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")

# Process AUC scores
auc_matrix <- diff_result$AUCell %>%
  dplyr::select(Name, avg_diff, cluster) %>%
  tidyr::spread(key = "cluster", value = "avg_diff") %>%
  mutate(Name = gsub("HALLMARK |-", " ", Name)) %>%
  filter(Name %in% hallmark_use) %>%
  arrange(match(Name, hallmark_use)) %>%
  column_to_rownames("Name") %>%
  dplyr::select("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")

# Create row annotations
row_anno <- data.frame(
  Genesets = rep(
    names(hallmark_categories),
    times = sapply(hallmark_categories, length)
  ),
  row.names = hallmark_use
)

# Define colors
anno_colors <- list(
  Genesets = c(
    "Death and stress" = "#1f77b4",
    "Differentiation" = "#ff7f0e",
    "Immune" = "#2ca02c",
    "Metabolism" = "#9467bd",
    "Proliferation" = "#d62728"
  )
)

# Custom function for significance stars
sig_stars <- function(j, i, x, y, width, height, fill) {
  if (auc_matrix[i, j] > 0.3 & pval_matrix[i, j] < 0.001) {
    grid.text("**", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
  }
}

# Create heatmap
heatmap_plot <- ComplexHeatmap::pheatmap(
  auc_matrix,name = "AUCell score",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames = TRUE,show_rownames = TRUE,
  border = TRUE,scale = "row",
  color = colorRamp2(
    seq(-2, 2, length.out = 100),
    colors = c("#2166ac", "#f7fbff", "#b2182b")
  ),
  annotation_row = row_anno,annotation_colors = anno_colors,
  cell_fun = sig_stars,row_names_gp = gpar(fontsize = 8)
)

# Save plot
output_dir <- "～/Figure2/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf(
  file.path(output_dir, "Fig2b.pdf"),
  width = 5.5,
  height = 4.10
)
draw(heatmap_plot)
dev.off()

############################## > Fig2c  ##################################
# Load and prepare gene sets
geneset <- read.gmt("~/data/h.all.v7.5.1.symbols.gmt") %>%
  mutate(term = str_remove_all(term, "HALLMARK_") %>%
           mutate(term = str_replace_all(term, "_", " ")) %>%
           split(.$term) %>%
           lapply(function(x) x[[2]]))
         
# Calculate AUC scores
score <- AUCell_run(merge_data@assays$RNA@counts, geneset)@assays@data$AUC
AUC.df <- data.frame(t(score)) %>%
  mutate(cluster = merge_data$stage) %>%
  gather(Pathway, AUCell, -cluster) %>%
  mutate(Pathway = str_replace_all(Pathway, "_", " "))

### prepare P53 pathway data
p53_data <- AUC.df %>%
  filter(Pathway == "P53 PATHWAY",
         cluster %in% c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")) %>%
  mutate(cluster = factor(cluster, levels = c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")))

# Define color palette
stage_colors <- c("AT2 cell" = "#BBDD78", 
                  "MIA" = "#6CB8D2",
                  "IA" = "#D55640", 
                  "Advanced" = "#89558D", 
                  "Metastatic" = "#D1352B")

# Create and save boxplot
comparisons <- list(
  c("MIA", "AT2 cell"), 
  c("IA", "Metastatic"),
  c("IA", "Advanced"),
  c("MIA", "IA")
)

p53_plot <- ggboxplot(
  tnfa_data, x = "cluster", y = "AUCell",
  color = "black",fill = "cluster",width = 0.4,
  palette = "npg",xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)) +
  stat_compare_means(
    comparisons = comparisons,method = "wilcox.test",label = "p.format",
  size = 4,bracket.size = 0.5,
  tip.length = 0,label.y = c(0.3, 0.28, 0.35, 0.40))+
  scale_fill_manual(values = stage_colors) +
  labs(y = "P53 PATHWAY score", x = "") +
  coord_cartesian(ylim = c(0, 0.4)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1, 
      size = 12, family = "sans", color = "black"
    ),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save plot
output_dir <- "~/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "Fig2c_p53.pdf"), p53_plot, width = 3.16, height = 4.55)

### Prepare TNF-α signaling data
tnfa_data <- data %>%
  filter(celltype == "TNFA SIGNALING VIA NFKB",
         cluster %in% c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")) %>%
  mutate(cluster = factor(cluster, levels = c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")))

# Define  comparisons
stage_colors <- c("AT2 cell" = "#BBDD78", 
                  "MIA" = "#6CB8D2",
                  "IA" = "#D55640", 
                  "Advanced" = "#89558D", 
                  "Metastatic" = "#D1352B")

comparisons <- list(
  c("MIA", "AT2 cell"),
  c("IA", "Metastatic"),
  c("IA", "Advanced"), 
  c("MIA", "IA")
)

# Create the boxplot
tnfa_plot <- ggboxplot(
  tnfa_data, x = "cluster", y = "AUCell",
  color = "black",fill = "cluster",width = 0.4,
  palette = "npg",xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)
) +
  labs(y = "TNFA SIGNALING VIA NFKB score", x = "") +
  stat_compare_means(
    comparisons = comparisons,label = "p.format",method = "wilcox.test",
    size = 4,vjust = -0.2,bracket.size = 0.5,
    tip.length = 0,label.y = c(0.37, 0.35, 0.42, 0.47)
  ) +
  scale_fill_manual(values = stage_colors) +
  coord_cartesian(ylim = c(0, 0.45)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1, 
      size = 12, family = "sans", color = "black"
    ),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save the plot
output_dir <- "~/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "Fig2c_tnfa.pdf"), p53_plot, width = 3.16, height = 4.55)

### Prepare WNT signaling data
wnt_data <- data %>%
  filter(celltype == "WNT.BETA.CATENIN.SIGNALING",
         cluster %in% c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")) %>%
  mutate(cluster = factor(cluster, levels = c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")))

# Define  comparisons
comparisons <- list(
  c("MIA", "AT2 cell"),
  c("MIA", "IA"),
  c("IA", "Advanced"),
  c("IA", "Metastatic")
)

# Create the boxplot
wnt_plot <- ggboxplot(
  wnt_data,x = "cluster",y = "AUCell",
  color = "black",fill = "cluster",width = 0.4,
  xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)
) +
  labs(y = "WNT BETA CATENIN SIGNALING score", x = "") +
  stat_compare_means(
    comparisons = comparisons,method = "wilcox.test",
    label = "p.format",size = 4,vjust = -0.2,
    bracket.size = 0.5,tip.length = 0,
    label.y = c(0.15, 0.17, 0.19, 0.22)
  ) +
  scale_fill_manual(values = stage_colors) +
  coord_cartesian(ylim = c(0, 0.2)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1,
      size = 12, family = "sans", color = "black"
    ),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save the plot
output_dir <- "~/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(
  file.path(output_dir, "Fig2c_wnt.pdf"),wnt_plot,width = 3.16,height = 4.55
)

# Prepare MYC TARGETS V1 data
myc_data <- data %>%
  filter(celltype == "MYC TARGETS V1",
         cluster %in% c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")) %>%
  mutate(cluster = factor(cluster, levels = c("AT2 cell", "MIA", "IA", "Advanced", "Metastatic")))

# statistical comparisons
comparisons <- list(
  c("MIA", "AT2 cell"),
  c("MIA", "IA"),
  c("IA", "Advanced"),
  c("IA", "Metastatic")
)

# Create the boxplot visualization
myc_plot <- ggboxplot(
  myc_data,x = "cluster",y = "AUCell",
  color = "black",fill = "cluster",width = 0.4,
  xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)
) +
  labs(y = "MYC TARGETS V1 score", x = "") +
  stat_compare_means(
    comparisons = comparisons,method = "wilcox.test",label = "p.format",
    size = 4,vjust = -0.2,bracket.size = 0.5,tip.length = 0,
    label.y = c(0.5, 0.60, 0.67, 0.74)
  ) +
  scale_fill_manual(values = stage_colors) +
  coord_cartesian(ylim = c(0, 0.8)) +  # Adjusted y-axis limit to better show all comparisons
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1,
                               size = 12, family = "sans", color = "black"),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save the plot
output_dir <- "~/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(
  file.path(output_dir, "Fig2c_myc.pdf"),  # Fixed filename with underscore
  myc_plot,
  width = 3.16,
  height = 4.55
)

############################## > Fig2d  ##################################
# Load and prepare scHCL data
hcl_result <- readRDS("~/data/scHCL_merge_data.rds")
merge_data$scHCL <- hcl_result$scHCL
data <- data.frame(scHCL = merge_data$scHCL, Subtype = merge_data$stage)

# Define cell type classifications
cell_type_mapping <- list(
  "AT1 cell" = c("AT1.cell..Adult.Lung1.", "AT1.cell.Adult.Lung2.", "AT1.cell.Adult.Lung3."),
  "AT2 cell" = c("AT2.cell.Adult.Lung1.", "AT2.cell.Adult.Lung2.", "AT2.cell.Adult.Lung3."),
  "Club cell" = c("Club.cell_KLK11.high.Adult.Lung2.", "Club.cell.Adult.Lung1.", "Club.cell.Adult.Lung3."),
  "Basal cell" = "Basal.cell_S100A2.high.Adult.Trachea2.",
  "Ciliated cell" = "Ciliated.cell.Adult.Lung3.",
  "Epithelium placenta like cell" = c("Epi.Placenta_VentoTormo.", "VCT2.Placenta_VentoTormo.", "VCT3.Placenta_VentoTormo."),
  "Goblet cell" = c("Goblet.cell.Adult.Esophagus2.", "Goblet.cell.Adult.Trachea2."),
  "Epithelium luminal like cell" = c("Luminal.cell_AGR2.high.Breast.Epithelium_Nguyen.", 
                                     "Luminal.cell_CD74.high.Breast.Epithelium_Nguyen.",
                                     "Luminal.cell_KRT23.high.Breast.Epithelium_Nguyen.",
                                     "Luminal.cell_SAA2.high.Breast.Epithelium_Nguyen."),
  "Alveolar bipotent intermediate cell" = "Alveolar.bipotent.intermediate.cell.Adult.Lung1."
)

# Apply cell type classification
data$scHCL1 <- "Other"
for (cell_type in names(cell_type_mapping)) {
  data$scHCL1[data$scHCL %in% cell_type_mapping[[cell_type]]] <- cell_type
}

# Set factor levels
data$Subtype <- factor(data$Subtype, 
                       levels = c("AT1 cell", "AT2 cell", "MIA", "IA", "Advanced", "Metastatic"))
data$scHCL1 <- factor(data$scHCL1,
                      levels = rev(c("AT1 cell", "AT2 cell", "Club cell", "Basal cell", 
                                     "Ciliated cell", "Epithelium luminal like cell",
                                     "Epithelium placenta like cell", 
                                     "Alveolar bipotent intermediate cell",
                                     "Goblet cell", "Other")))

# Calculate proportions and counts
Prop_data <- data %>% 
  group_by(Subtype, scHCL1) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(Freq = n/sum(n)) %>% 
  ungroup()

Num_data <- data %>% 
  group_by(Subtype, scHCL1) %>% 
  summarise(Number = n(), .groups = "drop") %>% 
  mutate(scHCL1 = factor(scHCL1, levels = c("AT1 cell", "AT2 cell", "Club cell", "Basal cell",
                                            "Ciliated cell", "Epithelium luminal like cell",
                                            "Epithelium placenta like cell",
                                            "Alveolar bipotent intermediate cell",
                                            "Goblet cell", "Other")))

# Define colors
col_cell <- pal_d3("category20")(20)[c(6,3,10,11,12,14,15,16,17,18)]

# Create stacked bar plot (top panel)
p_top <- ggplot(Num_data, aes(x = Subtype, y = Number, fill = scHCL1)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_cell, name = "scHCL inferred") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    plot.margin = unit(c(1,0,0,1), "lines"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Create dot plot (center panel)
p_center <- ggplot(Prop_data, aes(x = Subtype, y = scHCL1)) +
  geom_point(aes(size = Freq, color = scHCL1)) +
  scale_color_manual(values = rev(col_cell)) +
  labs(size = "Proportion", color = "Cell type") +
  guides(color = "none") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, size = 12, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_blank(),
    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
    plot.margin = unit(c(0,0,1,1), "lines"),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Combine and save plots
output_dir <- "~/figure2/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf(file.path(output_dir, "Fig2d.pdf"), width = 8.16, height = 9)
ggarrange(
  p_top,
  p_center,
  nrow = 2,
  ncol = 1,
  align = "v",
  heights = c(2,5),
  common.legend = FALSE
)
dev.off()

############################## > Fig2e  ##################################
# Load and prepare epithelial markers
Epi_markers <- read.csv("~/data/Epi_markers_from_normal_lung.csv")
Epi_markers$celltype <- factor(Epi_markers$celltype,
                               levels = c("AT2 cell", "AT1 cell", "Club cell", "Basal cell", "Ciliated cell"))

# Prepare gene lists for AUCell
genelist <- Epi_markers %>%
  dplyr::select(gene, celltype) %>%
  split(.$celltype) %>%
  lapply(function(x) x[[1]])

# Calculate AUC scores
score <- AUCell_run(merge_data@assays$RNA@counts, genelist)
AUC.df <- data.frame(t(score@assays@data$AUC))
AUC.df$stage <- merge_data@meta.data[,"stage"]  # Changed from 'object' to 'merge_data' for consistency

# Prepare plotting data
plot_data <- AUC.df %>%
  dplyr::select(stage, AT1.cell, AT2.cell) %>%
  tidyr::gather(celltype, AUCell, -stage) %>%
  mutate(
    celltype = gsub("\\.", " ", celltype),
    stage = factor(stage, levels = c("AT1 cell", "AT2 cell", "MIA", "IA", "Advanced", "Metastatic"))
  )

### AT2 
data_AT2 <- plot_data %>% filter(celltype == "AT2 cell")

# Define plot parameters
at2_colors <- pal_d3("category20")(20)[c(6,3,10,11,12,14,15,16,17,18)]
comparisons <- list(
  c("IA", "Metastatic"),
  c("IA", "Advanced"), 
  c("MIA", "IA"),
  c("AT2 cell", "MIA")
)

# Create AT2 plot
at2_plot <- ggboxplot(
  data_AT2,x = "stage", y = "AUCell",
  color = "black",fill = "stage",width = 0.4,
  xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)
) +
  stat_compare_means(
    comparisons = comparisons,method = "wilcox.test",  label = "p.format",     
    size = 4,vjust = -0.2,bracket.size = 0.5,
    tip.length = 0,label.y = c(0.80, 0.9, 1.0, 1.1)
  ) +
  scale_fill_manual(values = at2_colors) +
  labs(y = "AT2 score", x = "") +
  coord_cartesian(ylim = c(0, 1.1)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1,
      size = 14, family = "sans", color = "black"
    ),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save plot
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(
  file.path(output_dir, "Fig2e_AT2.pdf"),  # Added .pdf extension
  at2_plot,
  width = 3.46,
  height = 4.55
)

### AT1
at1_colors <- pal_d3("category20")(20)[c(3,6,10,11,12,14,15,16,17,18)]
at1_plot <- ggboxplot(
  data_AT1,x = "stage", y = "AUCell",
  color = "black",fill = "stage",width = 0.4,
  xlab = FALSE,alpha = 0.6,bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.2,add.params = list(color = "black", size = 0.4, alpha = 0.5)
) +
  stat_compare_means(
    comparisons = comparisons,method = "wilcox.test",  label = "p.format",     
    size = 4,vjust = -0.2,bracket.size = 0.5,
    tip.length = 0,label.y = label.y = c(0.65,0.75,0.85,0.95)
  ) +
  scale_fill_manual(values = at2_colors) +
  labs(y = "AT1 score", x = "") +
  coord_cartesian(ylim = c(0, 1.1)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1,
      size = 14, family = "sans", color = "black"
    ),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

ggsave(
  file.path(output_dir, "Fig2e_AT2.pdf"),  
  at2_plot,
  width = 3.46,
  height = 4.55
)

############################## > Fig2f  ##################################
# Prepare data for Monocle analysis
count_matrix <- as(as.matrix(merge_data@assays$RNA@counts), 'sparseMatrix')

# Create metadata objects
gene_metadata <- data.frame(
  gene_short_name = rownames(merge_data@assays$RNA),
  row.names = rownames(merge_data@assays$RNA)
)
cell_metadata <- merge_data@meta.data

# Create CellDataSet object
cds <- newCellDataSet(
  cellData = count_matrix,
  phenoData = new('AnnotatedDataFrame', data = cell_metadata),
  featureData = new('AnnotatedDataFrame', data = gene_metadata),
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit = 0.5
)

# Quality control and normalization
cds <- detectGenes(cds, min_expr = 0.5)
cds <- cds[Biobase::fData(cds)$num_cells_expressed > 10, ]
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, cores = 24, relative_expr = TRUE)

# Select highly variable genes for ordering
disp_table <- dispersionTable(cds)
ordering_genes <- subset(
  disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit
)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

# Dimensionality reduction and trajectory construction
cds <- reduceDimension(
  cds,method = "DDRTree",max_components = 2,norm_method = "log",pseudo_expr = 1
)
cds <- orderCells(cds)

# Define plotting parameters
pseudotime_colors <- scale_color_gradientn(values = seq(0, 1, 0.2),colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))

# Create pseudotime trajectory plot
pseudotime_plot <- plot_cell_trajectory(
  cds,color_by = "Pseudotime",cell_size = 0.5,show_branch_points = FALSE
) +
  pseudotime_colors +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_blank()
  )

# Save plot
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig2f_pseudotime.pdf"),
  pseudotime_plot,
  width = 4.03,
  height = 3.18
)
# 分段

# Define stage colors 
stage_colors <- pal_d3("category20")(20)[c(6, 3, 1, 4, 2)]
names(stage_colors) <- levels(cds$stage)  

# Create the trajectory plot by stage
stage_plot <- plot_cell_trajectory(
  cds,color_by = "stage",cell_size = 0.5,
  show_branch_points = FALSE
) +
  scale_color_manual(
    values = stage_colors,name = "Stage"  
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 3),ncol = 1  
  ))

# Save plot with consistent formatting
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  file.path(output_dir, "Fig2f_stage.pdf"), 
  stage_plot,width = 4.03, 
  height = 3.18,device = "pdf"
)

# Add AT1 scores to the CellDataSet object
cds$AT1 <- data_AT1[,"AUCell"]

# Define the color gradient for AT1 scores
at1_color_scale <- scale_color_gradientn(
  values = seq(0, 1, 0.2),
  colors = c("#330066", "#336699", "#66CC66", "#FFCC33"),
  name = "AT1 Score"  # Clear legend title
)

# Create the trajectory plot colored by AT1 score
at1_plot <- plot_cell_trajectory(
  cds,color_by = "AT1",cell_size = 0.5,
  show_branch_points = FALSE,show_tree = TRUE  
) +
  at1_color_scale +
  theme(
    legend.position = "right",legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),legend.key.height = unit(0.5, "cm")  
  ) +
  labs(title = NULL)  

# Ensure output directory exists
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the plot with high quality settings
ggsave(
  filename = file.path(output_dir, "Fig2f_AT1.pdf"),
  plot = at1_plot,
  width = 4.03,
  height = 3.18,
  device = cairo_pdf,  # For better font handling
  dpi = 300
)

# Add AT2 scores to the CellDataSet object
cds$AT2 <- data_AT2[,"AUCell"]

# Define the color gradient for AT2 scores (using same scheme as AT1 for consistency)
at2_color_scale <- scale_color_gradientn(
  values = seq(0, 1, 0.2),
  colors = c("#330066", "#336699", "#66CC66", "#FFCC33"),
  name = "AT2 Score"  # Clear legend title
)

# Create the trajectory plot colored by AT2 score
at2_plot <- plot_cell_trajectory(
  cds,color_by = "AT2",cell_size = 0.5,show_branch_points = FALSE,
  show_tree = TRUE  # Ensure trajectory tree is visible
) +
  at2_color_scale +
  theme(
    legend.position = c(0.85, 0.2),  # Adjusted position to avoid overlap
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", color = "grey80")
  ) +
  labs(title = NULL)  # Remove any default title

# Ensure output directory exists
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the plot with high quality settings
ggsave(
  filename = file.path(output_dir, "Fig2f_AT2.pdf"),
  plot = at2_plot,
  width = 4.03,
  height = 3.18,
  device = cairo_pdf,  # For better font handling
  dpi = 300
)

############################## > Fig2h  ##################################
## Expression and matrix transfer into S3 data
# Load TCGA data and prepare gene sets
TCGA <- readRDS("~/data/GDC_LUAD_early_stage.rds")
geneset <- list(
  AT2 = Epi_markers$gene[Epi_markers$celltype == "AT2 cell"],
  AT1 = Epi_markers$gene[Epi_markers$celltype == "AT1 cell"]
)

# Perform GSVA analysis
Exp <- as.matrix(GetAssayData(TCGA))
ssgsea <- gsva(Exp, gset.idx.list = geneset, kcdf = "Poisson", parallel.sz = 26)

# Add scores to object metadata
TCGA$AT1 <- as.numeric(ssgsea["AT1", ])
TCGA$AT2 <- as.numeric(ssgsea["AT2", ])

# Classify samples into High/Low groups
TCGA$AT_level <- case_when(
  TCGA$AT1 >= quantile(TCGA$AT1, 2/3) & TCGA$AT2 >= quantile(TCGA$AT2, 2/3) ~ "High",
  TCGA$AT1 <= quantile(TCGA$AT1, 1/3) & TCGA$AT2 <= quantile(TCGA$AT2, 1/3) ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
TCGA <- TCGA[, TCGA$AT_level %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ AT_level, data = TCGA@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = TCGA@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "AT1 & AT2 score",legend.labs = c("High", "Low"),
  legend = c(0.8, 0.8),surv.scale = "percent",xlab = "Months",
  palette = c("red", "black"),size = 0.5,
  risk.table.y.text = FALSE,pval.size = 6,ggtheme = theme_classic(),
  font = list(
    x = c(16, "plain"),
    y = c(16, "plain"),
    tickslab = c(14, "plain"),
    legend = c(16),
    legend.title = c(14, "italic")
  )
)

# Adjust risk table appearance
surv_plot$table <- surv_plot$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 12)
  )

# Save plot
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig2h_TCGA.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)

## GSE72094
GSE72094=readRDS("~/data/GSE72094_LUAD_early_stage.rds")
geneset <- list(
  AT2 = Epi_markers$gene[Epi_markers$celltype == "AT2 cell"],
  AT1 = Epi_markers$gene[Epi_markers$celltype == "AT1 cell"]
)

# Perform GSVA analysis
Exp <- as.matrix(GetAssayData(GSE72094))
ssgsea <- gsva(Exp, gset.idx.list = geneset, kcdf = "Poisson", parallel.sz = 26)

# Add scores to object metadata
GSE72094$AT1 <- as.numeric(ssgsea["AT1", ])
GSE72094$AT2 <- as.numeric(ssgsea["AT2", ])

# Classify samples into High/Low groups
GSE72094$AT_level <- case_when(
  GSE72094$AT1 >= quantile(GSE72094$AT1, 2/3) & GSE72094$AT2 >= quantile(GSE72094$AT2, 2/3) ~ "High",
  GSE72094$AT1 <= quantile(GSE72094$AT1, 1/3) & GSE72094$AT2 <= quantile(GSE72094$AT2, 1/3) ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
GSE72094<- GSE72094[, GSE72094$AT_level %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ AT_level, data = GSE72094@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = GSE72094@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "AT1 & AT2 score",legend.labs = c("High", "Low"),
  legend = c(0.8, 0.8),surv.scale = "percent",xlab = "Months",
  palette = c("red", "black"),size = 0.5,
  risk.table.y.text = FALSE,pval.size = 6,ggtheme = theme_classic(),
  font = list(
    x = c(16, "plain"),
    y = c(16, "plain"),
    tickslab = c(14, "plain"),
    legend = c(16),
    legend.title = c(14, "italic")
  )
)

# Adjust risk table appearance
surv_plot$table <- surv_plot$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 12)
  )

# Save plot
output_dir <- "~/figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig2h_GSE72094.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)



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

############################## > Fig4a  ##################################
# Generate and save UMAP plot
umap_plot <- UMAPPlot(Tumor_cell, group.by = "subtype", label = TRUE, label.size = 6) +
  NoLegend() +
  ggtitle("") +
  scale_color_manual(values = pal_d3("category20")(20)[c(4, 18, 3, 18)]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans")
  )

ggsave("~/figure4/Fig4a.pdf", umap_plot, width = 3.2, 
       height = 3.2, device = "pdf", bg = "transparent")

############################## > Fig4c  ##################################
# Calculate cell type proportions by patient
prop_data <- Tumor_cell@meta.data %>%
  dplyr::select(patient, stage, subtype) %>%
  dplyr::rename(group = subtype) %>%
  dplyr::count(patient, group, name = "n") %>%
  dplyr::group_by(patient) %>%
  dplyr::mutate(Freq = n/sum(n)) %>%
  dplyr::mutate(stage = ifelse(patient %in% MIA, "MIA", "IA")) 

prop_data$stage=factor(prop_data$stage,levels=c("MIA","IA"))
                                
# Create and save boxplot
p <- ggboxplot(prop_data, 
                       x = "group", y = "Freq", color = "stage",
                       palette = c("IA" = "#1F77B4", "MIA" = "#FF7F0E"),  # D3 category20 colors
                       add = "jitter", width = 0.4, alpha = 0.5, size = 0.4,
                       xlab = FALSE, legend = "right") +
  stat_compare_means(
    aes(group = stage), 
    method = "wilcox", 
    label = "p.signif",
    label.y = 0.95, size = 5) +
  ylim(0, 1) +
  labs(y = "Proportion", color = "Stage") +
 theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

ggsave("~/figure4/Fig4c.pdf", p, width = 3.85, height = 3.25, bg = "transparent")

############################## > Fig4d  ##################################
# Load and integrate datasets
eLUAD_GSE189357 <- readRDS("~/data/GSE189357.rds")
AdLUAD_MetLUAD <- merge_data[, merge_data$stage %in% c("Advanced","Metastatic")]

anchors <- FindIntegrationAnchors(list(eLUAD_GSE189357, AdLUAD_MetLUAD))
integrated_data <- IntegrateData(anchors)
integrated_data$stage <- factor(integrated_data$stage, levels=c("Early","Advanced","Metastatic"))

# Find marker genes
Tumor_cell$subtype <- factor(Tumor_cell$subtype, levels=c("C1","C2","C3","C4"))
Idents(Tumor_cell) <- "subtype"
Deg <- FindAllMarkers(Tumor_cell, logfc.threshold=0.25)
Deg <- dplyr::filter(Deg, avg_log2FC>0.5, cluster %in% paste0("C",1:4))
geneset = list(C1= Deg[,Deg$cluster %in% "C1"]$gene,
               C2= Deg[,Deg$cluster %in% "C2"]$gene,
               C3= Deg[,Deg$cluster %in% "C3"]$gene,
               C4= Deg[,Deg$cluster %in% "C4"]$gene)

# Calculate AUC scores
score <- AUCell_run(integrated_data@assays$RNA@counts, geneset)@assays@data$AUC
AUC.df <- data.frame(t(score)) %>% 
  dplyr::mutate(cluster=integrated_data$stage) %>%
  dplyr::select(cluster, everything())

# Prepare plot data
plot_data <- AUC.df %>% 
  tidyr::gather(celltype, AUCell, -cluster) %>%
  dplyr::mutate(celltype=gsub("\\."," ", celltype),
                cluster=factor(cluster, levels=c("Early","Advanced","Metastatic")))

# Create and save plot
p <- ggpubr::ggboxplot(plot_data, 
                       x="celltype", y="AUCell", fill="cluster",
                       palette=c("#6CB8D2","#89558D","#D1352B"),  # Blue, purple, red
                       width=0.4, xlab=FALSE, legend="right",
                       add="jitter", add.params=list(size=0.4, alpha=0.5)) +
  ggpubr::stat_compare_means(aes(group=cluster), method="t.test", label="p.signif") +
  ggplot2::ylim(0,1) +
  ggplot2::labs(y="AUCell score", x="", fill="Stage") +
  ggplot2::theme(
    panel.grid=element_blank(),
    panel.background=element_rect(fill="transparent"),
    axis.line=element_line(colour="black"),
    axis.text.x=element_text(angle=50, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14)
  )

ggsave("~/figure4/Fig4d.pdf", p, width=5, height=3.3, bg="transparent")

############################## > Fig4e  ##################################
# Create violin plot for CRABP2 expression
p <- Seurat::VlnPlot(
  Tumor_cell,
  features = "CRABP2",
  group.by = "subtype",
  cols = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")  # D3 category20 colors [4,2,3,1]
) +
  ggplot2::ggtitle("CRABP2") +
  ggplot2::xlab("") +
  ggplot2::ylab("Expression level") +
  ggplot2::theme(
    plot.title = element_text(face = "italic", size = 14, family = "sans"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans")
  )

# Save plot
ggsave(
  filename = "~/figure4/Fig4_CRABP2.pdf",
  plot = p,
  width = 3.76,
  height = 3.45,
  bg = "transparent"
)

# Generate LAMP3 expression violin plot
p <- VlnPlot(
  Tumor_cell,
  features = "LAMP3",
  group.by = "subtype",
  cols = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")  # D3 palette [4,2,3,1]
) + 
  labs(title = "LAMP3", x = "", y = "Expression level") +
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
ggsave("~/figure4/Fig4_LAMP3.pdf", p, width = 3.76, height = 3.45, bg = "transparent")

############################## > Fig4f  ##################################
# Process data
eLUAD_GSE189357 <- eLUAD_GSE189357 %>%
  SCTransform(return.only.var.genes = TRUE) %>%
  RunPCA(assay = "SCT") %>%
  RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.2)

# Calculate AUC scores
auc_scores <- AUCell_run(eLUAD_GSE189357@assays$RNA@counts, geneset)@assays@data$AUC
for(cluster in c("C1","C2","C3","C4")) {
  eLUAD_GSE189357[[cluster]] <- auc_scores[cluster,]
}

# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "CRABP2", group.by = "seurat_clusters") +
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

ggsave("Fig4f_CRABP2.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "LAMP3", group.by = "seurat_clusters") +
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

ggsave("Fig4f_LAMP3.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "C1", group.by = "seurat_clusters") +
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

ggsave("Fig4f_C1.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

# Create and save plot
p <- VlnPlot(eLUAD_GSE189357, features = "C3", group.by = "seurat_clusters") +
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

ggsave("Fig4f_C3.pdf", p, width = 3.14, height = 3.39, bg = "transparent")

############################## > Fig4h  ##################################
# Load data and set up cell groups
cellchat <- readRDS("~/data/LUAD_cellchat.rds")
tumor_cells <- "C1"
other_cells <- setdiff(unique(cellchat@idents), c("C2","C3","C4"))

# Extract tumor-specific interactions
net_matrix <- cellchat@net$count
tumor_interactions <- matrix(0, nrow(net_matrix), ncol(net_matrix), 
                             dimnames = dimnames(net_matrix))
tumor_interactions[tumor_cells, other_cells] <- net_matrix[tumor_cells, other_cells]
tumor_interactions <- tumor_interactions[-c(2:4), -c(2:4)]  # Remove C2-C4

# Create and save plot
p <- netVisual_circle(
  tumor_interactions,
  vertex.weight = rowSums(tumor_interactions),
  weight.scale = TRUE,
  label.edge = TRUE,
  vertex.label.cex = 1,
  arrow.size = 0,
  edge.label.cex = 1,
  title.name = "",
  color.use = c("#1F77B4", "#8C564B", "#9467BD", "#E377C2", 
                "#FF7F0E", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78")  # D3 palette [4,13,6,10,2,5,12,1,14]
) +
  ggtitle("Number of interactions") +
  theme(plot.title = element_text(size = 14, face = "plain"))

ggsave("fig4_Cir_C1_main.pdf", p, width = 5.4, height = 4.7)

tumor_cells <- "C3"
other_cells <- setdiff(unique(cellchat@idents), c("C1","C2","C4"))

# Extract tumor-specific interactions
net_matrix <- cellchat@net$count
tumor_interactions <- matrix(0, nrow(net_matrix), ncol(net_matrix), 
                             dimnames = dimnames(net_matrix))
tumor_interactions[tumor_cells, other_cells] <- net_matrix[tumor_cells, other_cells]
tumor_interactions <- tumor_interactions[-c(1,2,4), -c(1,2,4)]  # Remove C2-C4

# Create and save plot
p <- netVisual_circle(
  tumor_interactions,
  vertex.weight = rowSums(tumor_interactions),
  weight.scale = TRUE,
  label.edge = TRUE,
  vertex.label.cex = 1,
  arrow.size = 0,
  edge.label.cex = 1,
  title.name = "",
  color.use = c("#1F77B4", "#8C564B", "#9467BD", "#E377C2", 
                "#FF7F0E", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78")  # D3 palette [4,13,6,10,2,5,12,1,14]
) +
  ggtitle("Number of interactions") +
  theme(plot.title = element_text(size = 14, face = "plain"))

ggsave("fig4_Cir_C3_main.pdf", p, width = 5.4, height = 4.7)

############################## > Fig4i  #################################
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
geneset = list(C1= Deg[,Deg$cluster %in% "C1"]$gene,
               C2= Deg[,Deg$cluster %in% "C2"]$gene,
               C3= Deg[,Deg$cluster %in% "C3"]$gene,
               C4= Deg[,Deg$cluster %in% "C4"]$gene)
# Perform GSVA analysis
Exp <- as.matrix(GetAssayData(GSE72094))
ssgsea <- gsva(Exp, gset.idx.list = geneset, kcdf = "Poisson", parallel.sz = 26)

# Add scores to object metadata
GSE72094$C1 <- as.numeric(ssgsea["C1", ])
GSE72094$C3 <- as.numeric(ssgsea["C3", ])

# Classify samples into High/Low groups
GSE72094$C1_level <- case_when(
  GSE72094$C1 >= quantile(GSE72094$C1, 2/3) ~ "High",
  GSE72094$C1 <= quantile(GSE72094$C1, 1/3)  ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
GSE72094_C1<- GSE72094[, GSE72094$C1 %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ C1_level, data = GSE72094_C1@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = GSE72094_C1@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "C1 score",legend.labs = c("High", "Low"),
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
ggsave(
  file.path(output_dir, "Fig4i_GSE72094_C1.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)

# Classify samples into High/Low groups
GSE72094$C3_level <- case_when(
  GSE72094$C3 >= quantile(GSE72094$C3, 2/3) ~ "High",
  GSE72094$C3 <= quantile(GSE72094$C3, 1/3)  ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
GSE72094_C3<- GSE72094[, GSE72094$C3 %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ C3_level, data = GSE72094_C3@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = GSE72094_C3@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "C3 score",legend.labs = c("High", "Low"),
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

ggsave(
  file.path(output_dir, "Fig4i_GSE72094_C3.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)

## TCGA
TCGA=readRDS("~/data/GDC_LUAD_early_stage.rds")
geneset = list(C1= Deg[,Deg$cluster %in% "C1"]$gene,
               C2= Deg[,Deg$cluster %in% "C2"]$gene,
               C3= Deg[,Deg$cluster %in% "C3"]$gene,
               C4= Deg[,Deg$cluster %in% "C4"]$gene)
# Perform GSVA analysis
Exp <- as.matrix(GetAssayData(TCGA))
ssgsea <- gsva(Exp, gset.idx.list = geneset, kcdf = "Poisson", parallel.sz = 26)

# Add scores to object metadata
TCGA$C1 <- as.numeric(ssgsea["C1", ])
TCGA$C3 <- as.numeric(ssgsea["C3", ])

# Classify samples into High/Low groups
TCGA$C1_level <- case_when(
  TCGA$C1 >= quantile(TCGA$C1, 2/3) ~ "High",
  TCGA$C1 <= quantile(TCGA$C1, 1/3)  ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
TCGA_C1<- TCGA[, TCGA$C1 %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ C1_level, data = TCGA_C1@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = TCGA_C1@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "C1 score",legend.labs = c("High", "Low"),
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
ggsave(
  file.path(output_dir, "Fig4i_TCGA_C1.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)

# Classify samples into High/Low groups
TCGA$C3_level <- case_when(
  TCGA$C3 >= quantile(TCGA$C3, 2/3) ~ "High",
  TCGA$C3 <= quantile(TCGA$C3, 1/3)  ~ "Low",
  TRUE ~ "Mid"
)

# Filter for High/Low groups only
TCGA_C3<- TCGA[, TCGA$C3 %in% c("High", "Low")]

# Survival analysis
fit <- surv_fit(Surv(time, status) ~ C3_level, data = TCGA_C3@meta.data)

# Create survival plot
surv_plot <- ggsurvplot(
  fit,data = TCGA_C3@meta.data,risk.table = TRUE,
  pval = TRUE,legend.title = "C3 score",legend.labs = c("High", "Low"),
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

ggsave(
  file.path(output_dir, "Fig4i_TCGA_C3.pdf"),
  surv_plot,
  width = 4.9,
  height = 5.25,
  device = cairo_pdf
)

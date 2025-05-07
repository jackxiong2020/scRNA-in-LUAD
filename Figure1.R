############################## > scRNA process and cell annotation  ##################################
# Load required libraries
library(Seurat)
library(harmony)
library(dplyr)
library(infercnv)
library(ggplot2)
library(pheatmap)

# Define sample names
sample_names <- c("P01", "P02", "03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12")

# Initialize list to store Seurat objects
sc_data <- list()

# Function to load and preprocess single-cell data
load_and_process_sc_data <- function(sample_name) {
  # Load 10X data and create Seurat object
  data_dir <- paste0("~/cellranger/", sample_name, "/outs/filtered_feature_bc_matrix")
  seurat_obj <- Read10X(data.dir = data_dir) %>%
    CreateSeuratObject(project = sample_name, min.cells = 3, min.features = 500)
  
  return(seurat_obj)
}

# Function to remove doublets using DoubletFinder
remove_doublets <- function(seurat_obj) {
  # Preprocess for DoubletFinder
  seurat_obj <- seurat_obj %>%
    SCTransform(return.only.var.genes = TRUE) %>%
    RunPCA(assay = "SCT") %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1.0)
  
  # Find optimal pK parameter
  sweep_res <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  optimal_pk <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Calculate expected doublet rate
  doublet_rate <- ncol(seurat_obj) * 8 * 1e-6
  homotypic_prop <- modelHomotypic(seurat_obj$seurat_clusters)
  n_exp_poi <- round(doublet_rate * ncol(seurat_obj))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(
    seurat_obj, 
    PCs = 1:20, 
    pN = 0.25, 
    pK = optimal_pk, 
    nExp = n_exp_poi_adj, 
    reuse.pANN = FALSE, 
    sct = TRUE
  )
  
  # Update metadata
  df_col <- colnames(seurat_obj@meta.data)[grepl("DF.classification", colnames(seurat_obj@meta.data))]
  seurat_obj@meta.data[, df_col] <- ifelse(seurat_obj@meta.data[, df_col] == "Doublet", "doublet", "singlet")
  seurat_obj$DoubletFinder <- seurat_obj@meta.data[, df_col]
  
  # Filter out doublets
  seurat_obj <- seurat_obj[, seurat_obj$DoubletFinder == "singlet"]
  
  return(seurat_obj)
}

# Function for quality control
apply_quality_control <- function(seurat_obj) {
  # Filter cells with high mitochondrial gene percentage
  seurat_obj <- PercentageFeatureSet(
    seurat_obj, 
    pattern = "^MT-", 
    col.name = "percent_mito"
  )
  seurat_obj <- seurat_obj[, seurat_obj$percent_mito < 20]
  
  # Filter cells with extremely high gene counts (top 2%)
  seurat_obj <- seurat_obj[, seurat_obj$nFeature_RNA < quantile(seurat_obj$nFeature_RNA, 0.98)]
  
  return(seurat_obj)
}

# Function for cell cycle correction
correct_cell_cycle <- function(seurat_obj) {
  s_genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features = s_genes,
    g2m.features = g2m_genes
  )
  
  return(seurat_obj)
}

# Main processing loop
for (sample in sample_names) {
  cat("Processing sample:", sample, "\n")
  
  # Step 1: Load data
  sc_data[[sample]] <- load_and_process_sc_data(sample)
  
  # Step 2: Remove doublets
  sc_data[[sample]] <- remove_doublets(sc_data[[sample]])
  
  # Step 3: Quality control
  sc_data[[sample]] <- apply_quality_control(sc_data[[sample]])
  
  # Step 4: Cell cycle correction
  sc_data[[sample]] <- correct_cell_cycle(sc_data[[sample]])
}

# Merge all Seurat objects into one dataset
scLUAD <- merge(
  x = sc_data[[1]],
  y = do.call(c, sc_data[2:12]),
  add.cell.ids = sample_names
)

# Add Patient information
scLUAD$Patient <- scLUAD$orig.ident

# Define clinical stages (MIA vs IA)
mia_patients <- c("P01", "P02", "P03", "P04", "P05", "P06")
ia_patients <- c("P07", "P08", "P09", "P10", "P11", "P12")  

# Add Stage information
scLUAD$Stage <- ifelse(
  scLUAD$Patient %in% mia_patients, 
  "MIA", 
  ifelse(scLUAD$Patient %in% ia_patients, "IA", NA)
)

# Process merged data: normalization, clustering, and dimensionality reduction
scLUAD <- scLUAD %>%
  SCTransform(return.only.var.genes = TRUE) %>%
  RunPCA(assay = "SCT") %>%
  RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.2)

# Define cell type markers
cell_markers <- list(
  Epithelial = c("EPCAM", "KRT19", "CDH1", "KRT18"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRAC"),
  B_cell = c("CD79A", "IGHM", "IGHG3", "IGHA2"),
  Myeloid = c("CD68", "MARCO", "FCGR3A", "LYZ"),
  NK_cell = c("NCAM1", "NKG7", "GNLY", "KLRD1"),
  Mast_cell = c("KIT", "MS4A2", "GATA2"),
  Fibroblast = c("DCN", "COL1A1", "COL1A2", "THY1"),
  Endothelial = c("PECAM1", "CLDN5", "FLT1", "RAMP2")
)

# Annotate cell types based on clusters
cluster_annotations <- list(
  "T cell" = c(0, 6, 13),
  "NK cell" = c(3),
  "Epithelial cell" = c(1, 11, 12),
  "B cell" = c(8, 9),
  "Myeloid cell" = c(2),
  "Fibroblast" = c(5, 10),
  "Mast cell" = c(7),
  "Endothelial cell" = c(4)
)

# Apply annotations
scLUAD$cell_type <- "Unknown"
for (cell_type in names(cluster_annotations)) {
  scLUAD$cell_type[scLUAD$seurat_clusters %in% cluster_annotations[[cell_type]]] <- cell_type
}

# Set factor levels for proper ordering in plots
scLUAD$cell_type <- factor(
  scLUAD$cell_type,
  levels = c(
    "Epithelial cell", "Endothelial cell", "Fibroblast",
    "T cell", "NK cell", "B cell", "Myeloid cell", "Mast cell"
  )
)
############################## > Fig1 b  ##################################
# Define groups and process metadata
mia_patients <- c("P01","P02","P03","P04","P05","P06")
ia_patients <- c("P07","P08","P09","P10","P11","P12")

ScLUAD$Patient <- ScLUAD$orig.ident
ScLUAD$Stage <- factor(ifelse(ScLUAD$Patient %in% mia_patients, "MIA", "IA"),
                       levels = c("MIA","IA"),
                       labels = sprintf("%s (n = %s)", c("MIA","IA"), 
                                        format(table(ScLUAD$Stage), big.mark = ",")))

# Create and save UMAP plot
umap_plot <- UMAPPlot(ScLUAD, group.by = "cell_type", split.by = "Stage",
                      label = F, pt.size = 0, raster = T,
                      cols = pal_d3("category20")(20)[c(4,6,10,2,5,12,1,14)]) + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, family = "sans"),
        strip.text = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 12, family = "sans"))

if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = T)
ggsave("~/figure1/Fig1b.pdf", umap_plot, width = 8, height = 3.33)

############################## > Fig1 c  ##################################
# Create and save dot plot
dot_plot <- DotPlot(ScLUAD, features = markers, group.by = "cell_type") +
  coord_flip() + 
  theme_bw() + 
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33")) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 12, family = "sans", face = "plain"),
    axis.text.y = element_text(size = 10, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    plot.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Save plot
if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = T)
ggsave("~/figure1/Fig1c.pdf", dot_plot, width = 4.85, height = 5.2)
############################## > Fig1 d  ##################################
# Create proportion plot
prop_data <-  as.data.frame(table(ScLUAD$Patient, ScLUAD$cell_type)) %>% 
  setNames(c("Samples", "Cell type", "Cell Number"))

prop_plot <- ggplot(
  data = prop_data,
  aes(x = Samples, weight = `Cell Number`, fill = `Cell type`)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = pal_d3("category20")(20)[c(4,6,10,2,5,12,1,14)]) +
  labs(x = "Patient", y = "Proportion") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 12, family = "sans", face = "plain"),
    axis.text.y = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    plot.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Save plot
if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = TRUE)
ggsave("~/figure1/Fig1d.pdf", prop_plot, width = 6, height = 3.3)
############################## > Fig1 e  ##################################
# Calculate cell type proportions by patient
prop_data <- ScLUAD@meta.data %>%
  dplyr::select(Patient, Stage, cell_type) %>%
  dplyr::group_by(Patient, cell_type) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(Freq = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Stage = ifelse(Patient %in% mia_patients, "MIA", "IA"))

# Create proportion boxplot
prop_plot <- ggpubr::ggboxplot(
  prop_data,x = "cell_type",y = "Freq",
  color = "Stage",add = "jitter",width = 0.4,
  palette = "npg",xlab = FALSE,legend = "right",
  alpha = 0.5,bxp.errorbar = TRUE,size = 0.4,
  bxp.errorbar.width = 0.5,
  add.params = list(size = 1, alpha = 1)
) +
  # Add statistical comparisons
  stat_compare_means(
    aes(group = Stage),
    label = "p.signif",
    method = "t.test",
    label.y = 0.9,
    size = 4
  ) +
  # Adjust plot aesthetics
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Proportion") +
  scale_color_manual(values = cols_cluster_stage[1:2]) +
  # Custom theme settings
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1,size = 12, family = "sans", face = "plain"),
    axis.text.y = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    plot.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = TRUE)
ggsave("~/figure1/Fig1e.pdf", prop_plot, width = 6, height = 3.83)
############################## > inferCNV  ##################################
setwd("~/infercnv")
# Rename cell types and prepare annotation
ScLUAD$cell_type[ScLUAD$cell_type == "Epithelial cell"] <- "Maligant"
normal_cells <- c("T cell","B cell","NK cell","Fibroblast","Mast cell","Endothelial cell","Myeloid cell")

# Run inferCNV
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = GetAssayData(ScLUAD, slot="counts", assay="RNA"),
  annotations_file = data.frame(ScLUAD$cell_type),
  delim = " ",
  gene_order_file = "gencode_v19_gene_pos.txt",
  ref_group_names = normal_cells
)

infercnv::run(infercnv_obj, cutoff=0.1, out_dir="infercnv", 
              cluster_by_groups=T, denoise=T, HMM=F, 
              num_threads=36, useRaster=FALSE)

# Post-processing
infercnv_t <- readRDS("~/infercnv/run.final.infercnv_obj")
expr <- infercnv_t@expr.data
gene_pos <- read.table("gencode_v19_gene_pos.txt", header=F, sep="\t")
rownames(gene_pos) <- gene_pos$V1

# Process common genes
common_genes <- intersect(rownames(expr), gene_pos$V1)
expr_filtered <- expr[common_genes,]
gene_pos_filtered <- gene_pos[common_genes,]

# Calculate CNV scores
normal_cnv <- colMeans((expr_filtered[,unlist(infercnv_t@reference_grouped_cell_indices)] - 1)^2)
tumor_cells <- infercnv_t@observation_grouped_cell_indices$Maligant

# Clustering analysis
set.seed(123)
kmeans_res <- kmeans(t(expr_filtered), centers=6)
cluster_df <- data.frame(
  CB = names(kmeans_res$cluster),
  kmeans_class = factor(kmeans_res$cluster),
  class = ifelse(names(kmeans_res$cluster) %in% tumor_cells, "Tumor", "Normal")
) %>% arrange(kmeans_class)

## Tumor cell indentify
tumor_df <- cluster_df %>% 
  filter(kmeans_class %in% c(1, 3, 4, 5, 6))

# Annotate tumor cells in Seurat object
ScLUAD$cell_type[colnames(ScLUAD) %in% tumor_df$CB] <- "Tumor cell"
############################## > Fig1 f  ##################################
# Prepare annotations
top_anno <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = "NA", col = "NA"), 
    labels = 1:22,
    labels_gp = gpar(cex = 1.0)
  ))
  
# Set color scheme
num_clusters <- length(unique(cluster_df$kmeans_class))
color_v <- RColorBrewer::brewer.pal(num_clusters, "Dark2")[1:num_clusters]
names(color_v) <- as.character(1:num_clusters)
  
left_anno <- rowAnnotation(
  df = cluster_df,
  col = list(
    class = c("Tumor" = "red"),
    kmeans_class = color_v
  )
)
  
# Create and save heatmap
ht <- Heatmap(
    t(expr)[rownames(cluster_df),],col = colorRamp2(c(0.8, 1, 1.2), c("#377EB8", "#F0F0F0", "#E41A1C")),
    cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,
    show_row_names = FALSE,column_split = factor(sub_geneFile$V2, levels = paste0("chr", 1:22)),
    column_gap = unit(2, "mm"),
    heatmap_legend_param = list(
      title = "Modified expression",
      direction = "vertical",
      title_position = "leftcenter-rot",
      at = c(0.8, 1, 1.2),
      legend_height = unit(3, "cm")
    ),
    top_annotation = top_anno,left_annotation = left_anno,
    row_title = NULL,column_title = NULL
  )
  
pdf("~/figure1/Fig1f.pdf", width = 20, height = 15)
  draw(ht, heatmap_legend_side = "right")
dev.off()
############################## > Fig1 g  ##################################
# Calculate CNV scores
CNV_score <- data.frame(
  CB = rownames(expr),
  score = colMeans((expr - 1)^2) * 10000
) %>%
  left_join(kmeans_df_s, by = "CB") %>%
  mutate(
    Class = ifelse(is.na(kmeans_class), "normal", as.character(kmeans_class)),
    Class = plyr::mapvalues(Class, 
                            from = c(1:6, "normal"), 
                            to = c("K1","K2","K3","K4","K5","K6","Normal")),
    Class = factor(Class, levels = c("K5","K1","K4","K3","K6","K2","Normal"))
  )

# Plot boxplot
color_pal <- c(pal_d3("category20")(20)[1:6], "red")
names(color_pal) <- c("K1","K2","K3","K4","K5","K6","Normal")

p <- ggpubr::ggboxplot(CNV_score, x = "Class", y = "score", fill = "Class") +
  stat_compare_means(ref.group = "Normal", label = "p.signif") +
  scale_fill_manual(values = color_pal) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(y = "InferCNV score", x = "Kmean class") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, 
                               size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    plot.title = element_text(size = 14, family = "sans")
  )

if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = TRUE)
ggsave("~/figure1/Fig1g.pdf", p, width = 3.2, height = 3.7)
############################## > Fig1 h  ##################################
## 肿瘤细胞
Tumor= ScLUAD[,which(ScLUAD$cell_type %in% "Tumor cell")]

# Data processing: Calculate tumor cell proportions per sample
tumor_prop_data <- ScLUAD[, ScLUAD$cell_type == "Tumor cell"]@meta.data %>%
  dplyr::count(Patient, Stage, cell_type, name = "n") %>%
  dplyr::group_by(Patient) %>% 
  dplyr::mutate(Freq = n / sum(n)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Stage = ifelse(Patient %in% mia_patients, "MIA", "IA"))

# Create boxplot 
tumor_prop_plot <- ggpubr::ggboxplot(
  tumor_prop_data, 
  x = "Stage", y = "Freq", color = "Stage",
  palette = "npg", add = "jitter", width = 0.4, 
  xlab = FALSE, legend = "right", alpha = 0.5,
  size = 0.4, bxp.errorbar = TRUE, 
  bxp.errorbar.width = 0.5,
  add.params = list(size = 1, alpha = 1)
) +
  # Add statistical test
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", 
    label.y = 0.9, 
    size = 4
  ) +
  # Set axes and labels
  coord_cartesian(ylim = c(0, 1)) + 
  labs(y = "Proportion") + 
  scale_color_manual(values = cols_cluster_stage[1:2]) +
  # Set theme
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 14, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

if(!dir.exists("~/figure1")) dir.create("~/figure1", recursive = TRUE)
ggsave("~/figure1/Fig1h.pdf", tumor_prop_plot, width = 3.05, height = 2.86)
## library packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(scHCL)
library(NMF)
library(SPOTlight)
library(SingleCellExperiment)
library(scater)
library(scran)
library(presto)
library(ggcorrplot)
library(Rcpp)
library(doParallel)
library(spacexr)
library(SpatialInferCNV)
library(ggpubr)

############################## > Fig3a  ##################################
# Subset and process Tumor_cell cells
Tumor_cell <- ScLUAD[, ScLUAD$cell_type == "Tumor_cell cell"] %>%
  SCTransform(assay = "RNA", return.only.var.genes = TRUE) %>%
  RunPCA(assay = "SCT") %>%
  RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.2)

# Annotate subtypes
cluster_mapping <- c("0" = "C1", "1" = "C2", "2" = "C3", "3" = "C4")
Tumor_cell$subtype <- cluster_mapping[as.character(Tumor_cell$seurat_clusters)]

# Prepare stage labels with counts
stage_counts <- table(Tumor_cell$stage)
Tumor_cell$stage <- factor(
  Tumor_cell$stage,
  levels = rev(c("MIA", "IA")),
  labels = rev(paste0(c("MIA", "IA"), " (n = ", format(stage_counts[c("MIA", "IA")], big.mark = ","), ")"))
)

# Define colors and create UMAP plot
cols_cluster <- pal_d3("category20")(20)[c(4, 2, 3, 1)]
umap_plot <- UMAPPlot(
  Tumor_cell,group.by = "subtype",split.by = "stage",
  label = TRUE,label.size = 5,
  pt.size = 0.1,cols = cols_cluster
) +
  ggtitle("") +
  theme(
    strip.text = element_text(size = 14, family = "sans", face = "plain"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.position = "none"
  )

# Save plot
output_dir <- "~/figure3/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig3a.pdf"),
  umap_plot,
  width = 4,
  height = 3.35,
  device = cairo_pdf
)

############################## > Fig3b  ##################################
# Prepare data for plotting
sample_table <- as.data.frame(table(
  Tumor_cell@meta.data$patient,
  Tumor_cell@meta.data$subtype
))
names(sample_table) <- c("Patient", "Subtype", "CellNumber")

# Define color palette
subtype_colors <- pal_d3("category20")(20)[c(4, 2, 3, 1, 5)]

# Create the proportion plot
prop_plot <- ggplot(sample_table, aes(x = Patient, y = CellNumber, fill = Subtype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_fill_manual(values = subtype_colors, name = "Subtype") +
  labs(x = "Patient", y = "Proportion") +
  scale_y_continuous(expand = c(0, 0)) + # Remove padding at bottom
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(
      size = 12, angle = 50, hjust = 1, vjust = 1,
      family = "sans", margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain"),
    plot.margin = unit(c(5, 5, 5, 5), "mm") # Balanced margins
  )

# Save plot
output_dir <- "~/figure3"

ggsave(
  file.path(output_dir, "Fig3b.pdf"),
  prop_plot,
  width = 4.47,
  height = 3.3,
  device = cairo_pdf
)

############################## > Fig3c  ##################################
# Set identities and factor levels
Idents(Tumor_cell) <- factor(Tumor_cell$subtype, levels = c("C1", "C2", "C3", "C4"))

# Find differentially expressed genes
Deg <- FindAllMarkers(
  Tumor_cell,logfc.threshold = 0.5,only.pos = TRUE,assay = "SCT"
)

# Filter significant markers
significant_markers <- Deg %>%
  dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create heatmap
heatmap_plot <- DoHeatmap(
  object = Tumor_cell,
  features = significant_markers$gene,
  group.by = "subtype",
  group.colors = pal_d3("category20")(20)[c(4, 2, 3, 1)],
  assay = "SCT",
  label = FALSE,
  size = 3.5  # Adjust size of color key
) + scale_fill_gradientn(colors = c("blue", "white", "red"))+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9.5,family = "sans",face = "italic"),
    legend.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 10, family = "sans"),
    legend.key.size = unit(0.4, "cm")
  ) 

# Save plot
output_dir <- "~/figure3"

ggsave(
  file.path(output_dir, "Fig3c.pdf"),
  heatmap_plot,
  width = 4.78,
  height = 4.96,
  device = cairo_pdf
)

############################## > Single cell spatial transcriptome process  ##################################
# Load and process spatial data
process_spatial <- function(object) {
  object %>%
    SCTransform(assay = "Spatial", return.only.var.genes = TRUE) %>%
    RunPCA(assay = "SCT") %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(reduction = "pca") %>%
    FindClusters(resolution = 0.6)
}

ST_P13 <- readRDS("~/data/P13_bin50.rds") %>% process_spatial()
ST_P14 <- readRDS("~/data/P14_bin50.rds") %>% process_spatial()

# RCTD analysis function with improved error handling
Run_RCTD <- function(ST_object, Ref_Sc_object, cluster = "cell_type") {

  # Prepare reference data
  counts_sc <- Ref_Sc_object@assays$RNA@counts
  clusters <- as.factor(Ref_Sc_object@meta.data[[cluster]])
  clusters <- droplevels(clusters)
  names(clusters) <- colnames(Ref_Sc_object) 
  nUMI <- Ref_Sc_object$nCount_RNA
  names(nUMI) <- colnames(Ref_Sc_object)
  
  # Build reference
  reference <- Reference(counts_sc, clusters, nUMI)
  
  # Prepare spatial data
  counts_st <- ST_object@assays$Spatial@counts
  coords <- GetTissueCoordinates(ST_object)
  colnames(coords) <- c("x", "y")
  coords <- coords[complete.cases(coords), ]  # Remove NA coordinates
  
  # Create SpatialRNA object
  query <- SpatialRNA(
    coords = coords,
    counts = counts_st,
    nUMI = colSums(counts_st))
    
    # Run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 48)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    # Add results to Seurat object
    ST_object <- AddMetaData(
      ST_object,
      metadata = RCTD@results$results_df
    )
    
    return(ST_object)
    }

# Run RCTD with error handling
ST_P13 <- Run_RCTD(ST_P13, LUAD, cluster = "cell_type")
ST_P14 <- Run_RCTD(ST_P14, LUAD, cluster = "cell_type")

############################## > Saptial infercnv  ##################################
# Set up analysis directory
setwd("~/ST_P13_infercnv/")

# Prepare cell annotations
ST_P13$cell_type[ST_P13$cell_type == "Epithelial cell"] <- "Malignant"  
cell_anno <- data.frame(row.names = colnames(ST_P13),cell_type = ST_P13$cell_type)

# Define normal reference cell types
normal_cells <- c("T cell", "B cell", "Fibroblast", "Mast cell", "Endothelial cell", "Myeloid cell")

# Get count matrix
gene_matrix <- GetAssayData(
  ST_P13,slot = "counts",assay = "Spatial"
)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = gene_matrix,annotations_file = cell_anno,delim = "\t", 
  gene_order_file = "gencode_v19_gene_pos.txt",ref_group_names = normal_cells
)

# Run inferCNV analysis
infercnv_results <- infercnv::run(
    infercnv_obj,cutoff = 0.1,cluster_by_groups = TRUE,
    denoise = TRUE,HMM = FALSE,num_threads = 36,
    useRaster = FALSE,plot_steps = FALSE, no_plot = FALSE,output_format = "pdf"
)

# Set up analysis directory
setwd("~/ST_P14_infercnv/")

# Prepare cell annotations
ST_P14$cell_type[ST_P14$cell_type == "Epithelial cell"] <- "Malignant"  
cell_anno <- data.frame(row.names = colnames(ST_P14),cell_type = ST_P14$cell_type)

# Define normal reference cell types
normal_cells <- c("T cell", "B cell", "Fibroblast", "Mast cell", "Endothelial cell", "Myeloid cell")

# Get count matrix
gene_matrix <- GetAssayData(
  ST_P14,slot = "counts",assay = "Spatial"
)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = gene_matrix,annotations_file = cell_anno,delim = "\t", 
  gene_order_file = "gencode_v19_gene_pos.txt",ref_group_names = normal_cells
)

# Run inferCNV analysis
infercnv_results <- infercnv::run(
  infercnv_obj,cutoff = 0.1,cluster_by_groups = TRUE,
  denoise = TRUE,HMM = FALSE,num_threads = 36,
  useRaster = FALSE,plot_steps = FALSE, no_plot = FALSE,output_format = "pdf"
)

############################## > Fig3d  ##################################
# Subset and annotate epithelial cells  P13
Epi_P13 <- ST_P13[, ST_P13$cell_type == "Epithelial cell"]
Epi_P13_RCTD <- Run_RCTD(Epi_P13, Tumor_cell, cluster = "subtype")

# Process cell type annotations
Epi_P13$subtype <- paste0(Epi_P13_RCTD$first_type, " cancer cell")

# Clean and consolidate cell types
Epi_P13 <- Epi_P13[, !(Epi_P13$subtype %in% "NA cancer cell")]  # Remove NA cells

Epi_P13$subtype <- case_when(
  Epi_P13$subtype %in% c("B cell", "Myeloid cell", "T cell") ~ "Immune cell",
  Epi_P13$subtype %in% c("Endothelial cell", "Fibroblast") ~ "Stromal cell",
  TRUE ~ Epi_P13$subtype
)

# Set factor levels with custom colors
subtype_levels <- c("C1 cancer cell", "C2 cancer cell", "C3 cancer cell", "C4 cancer cell", 
                    "Immune cell", "Stromal cell")
subtype_colors <- c("#D62728FF", "#FF7F0EFF", "#77dd77", "#1F77B4FF", 
                    "#C7C7C7FF", pal_d3("category20")(20)[6])

Epi_P13$subtype <- factor(Epi_P13$subtype, levels = subtype_levels)

# Create spatial plot
spatial_plot <- SpatialDimPlot(
  Epi_P13,group.by = "subtype",pt.size.factor = 1.2,  
  image.alpha = 0.8,     stroke = 0.1           
) +
  scale_fill_manual(
    values = subtype_colors,name = "Cell type",drop = FALSE         
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank()
  )

# Save plot
output_dir <- "~/figure3"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "Fig3d_P13.pdf"),spatial_plot,width = 7.78,height = 4.72,
  device = cairo_pdf
)

# Subset and annotate epithelial cells P14
Epi_P14 <- ST_P14[, ST_P14$cell_type == "Epithelial cell"]
Epi_P13_RCTD <- Run_RCTD(Epi_P14, Tumor_cell, cluster = "subtype")

# Process cell type annotations
Epi_P14$subtype <- paste0(Epi_P13_RCTD$first_type, " cancer cell")

# Clean and consolidate cell types
Epi_P14 <- Epi_P14[, !(Epi_P14$subtype %in% "NA cancer cell")]  # Remove NA cells

Epi_P14$subtype <- case_when(
  Epi_P14$subtype %in% c("B cell", "Myeloid cell", "T cell") ~ "Immune cell",
  Epi_P14$subtype %in% c("Endothelial cell", "Fibroblast") ~ "Stromal cell",
  TRUE ~ Epi_P14$subtype
)

# Set factor levels with custom colors
subtype_levels <- c("C1 cancer cell", "C2 cancer cell", "C3 cancer cell", "C4 cancer cell", 
                    "Immune cell", "Stromal cell")
subtype_colors <- c("#D62728FF", "#FF7F0EFF", "#77dd77", "#1F77B4FF", 
                    "#C7C7C7FF", pal_d3("category20")(20)[6])

Epi_P14$subtype <- factor(Epi_P14$subtype, levels = subtype_levels)

# Create spatial plot
spatial_plot <- SpatialDimPlot(
  Epi_P14,group.by = "subtype",pt.size.factor = 1.2,  
  image.alpha = 0.8,     stroke = 0.1           
) +
  scale_fill_manual(
    values = subtype_colors,name = "Cell type",drop = FALSE         
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank()
  )

# Save plot
output_dir <- "~/figure3"

ggsave(
  file.path(output_dir, "Fig3d_P14.pdf"),spatial_plot,width = 7.78,height = 4.72,
  device = cairo_pdf
)

############################## > Fig3e  ##################################
## 肿瘤细胞占比
# Merge spatial datasets
ST_LUAD <- merge(ST_P13, ST_P14)

# Calculate cell type proportions
sample_table <- as.data.frame(table(
  ST_LUAD@meta.data$orig.ident,
  ST_LUAD@meta.data$subtype
))
names(sample_table) <- c("Sample", "CellType", "CellCount")

# Define color scheme matching previous plots
subtype_colors <- c(
  pal_d3("category20")(20)[c(4, 2, 3, 1)],  # Cancer subtypes
  "#C7C7C7FF",                              # Immune cells
  pal_d3("category20")(20)[6]               # Stromal cells
)

# Create proportion plot
prop_plot <- ggplot(sample_table, aes(x = Sample, y = CellCount, fill = CellType)) +
  geom_bar(position = "fill", stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = subtype_colors,name = "Subtype",
    labels = c("C1", "C2", "C3", "C4", "Immune", "Stromal")  # Cleaner legend labels
  ) +
  labs(
    x = "Patient",y = "Proportion",title = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +  # Tight y-axis
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(
      size = 12,angle = 50,hjust = 1,
      vjust = 1,margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    plot.margin = unit(c(5, 5, 5, 5), "mm")
  )

ggsave(
  file.path(output_dir, "Fig3e.pdf"),
  prop_plot,
  width = 3.52,
  height = 3.3,
  device = cairo_pdf
)

############################## > Fig3f   ##################################
# Prepare data
merge_data$subtype <- factor(merge_data$subtype,
                             levels = c("AT1 cell", "AT2 cell", "C1", "C2", "C3", "C4", "Advanced", "Metastatic"))

data <- data.frame(
  scHCL = merge_data$scHCL,
  Subtype = merge_data$subtype
)

# Define cell type mappings using a more maintainable approach
cell_type_mapping <- list(
  "AT1 cell" = c("AT1.cell..Adult.Lung1.", "AT1.cell.Adult.Lung2.", "AT1.cell.Adult.Lung3."),
  "AT2 cell" = c("AT2.cell.Adult.Lung1.", "AT2.cell.Adult.Lung2.", "AT2.cell.Adult.Lung3."),
  "Club cell" = c("Club.cell_KLK11.high.Adult.Lung2.", "Club.cell.Adult.Lung1.", "Club.cell.Adult.Lung3."),
  "Basal cell" = c("Basal.cell_S100A2.high.Adult.Trachea2.", "Basal.cell_KRT6A.high.Adult.Trachea2."),
  "Ciliated cell" = "Ciliated.cell.Adult.Lung3."
)

# Apply cell type classification
data$scHCL1 <- "Other"
for (cell_type in names(cell_type_mapping)) {
  data$scHCL1[data$scHCL %in% cell_type_mapping[[cell_type]]] <- cell_type
}

# Calculate proportions and counts
Prop_data <- data %>% 
  group_by(Subtype, scHCL1) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(Freq = n/sum(n)) %>% 
  ungroup()

Num_data <- data %>% 
  group_by(Subtype, scHCL1) %>% 
  summarise(Number = n(), .groups = "drop")

# Set factor levels
cell_levels <- c("AT2 cell", "AT1 cell", "Club cell", "Basal cell", "Ciliated cell", "Other")
Prop_data$scHCL1 <- factor(Prop_data$scHCL1, levels = rev(cell_levels))
Num_data$scHCL1 <- factor(Num_data$scHCL1, levels = cell_levels)

# Define colors
col_cell <- c("#77dd77", pal_d3("category20")(20)[c(6,10,11,12,18)])

# Create top panel (stacked bar plot)
p_top <- ggplot(Num_data, aes(x = Subtype, y = Number, fill = scHCL1)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = col_cell, name = "scHCL inferred") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 12, family = "sans"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = unit(c(1,0,0,1), "lines")
  )

# Create center panel (dot plot)
p_center <- ggplot(Prop_data, aes(x = Subtype, y = scHCL1)) +
  geom_point(aes(size = Freq, color = scHCL1)) +
  scale_color_manual(values = rev(col_cell)) +
  scale_size_continuous(range = c(1, 5)) +
  labs(x = "", y = "", size = "Proportion") +
  guides(color = "none") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.7, color = "black"),
    axis.text.x = element_text(angle = 45, size = 12, hjust = 1, vjust = 1, 
                               family = "sans", color = "black"),
    axis.text.y = element_text(size = 12, family = "sans", color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.margin = unit(c(0,0,1,1), "lines"),
    legend.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

output_dir <- "~/figure3/"

ggsave(
  file.path(output_dir, "fig3f.pdf"),
  umap_plot,
  width = 4.85,
  height = 5.15,
  device = cairo_pdf
)

############################## > Fig3g  ##################################
# Load and prepare epithelial markers
epi_markers <- read.csv("~/data/Epi_markers_from_normal_lung.csv")

# Define cell type factor levels
celltype_levels <- c("AT2 cell", "AT1 cell", "Club cell", "Basal cell", "Ciliated cell")
epi_markers$celltype <- factor(epi_markers$celltype, levels = celltype_levels)

# Prepare gene lists for AUCell
gene_lists <- epi_markers %>%
  dplyr::select(gene, celltype) %>%
  split(.$celltype) %>%
  lapply(function(x) x[[1]])

# Calculate AUC scores 
auc_scores <- AUCell_run(merge_data@assays$RNA@counts, gene_lists)
auc_matrix <- auc_scores@assays@data$AUC

# Prepare data frame for plotting 
auc_df <- data.frame(t(auc_matrix)) %>%
  dplyr::mutate(cluster = merge_data$subtype) %>%
  dplyr::select(cluster, everything())

# Define colors 
# Base colors from pal_d3
col_cluster=c(
  "AT1"="#F5D2A8",
  "AT2"="#BBDD78",
  "C1"="#D55640",
  "C2"="#E69F84",
  "C3"="#479D88",
  "C4"="#6CB8D2",
  "Advanced"="#89558D",
  "Metastatic"="#D1352B"
  )

# Set cluster factor levels 
cluster_levels <- rev(c("AT1 cell", "AT2 cell", "C3", "C2", "C4", "C1", 
                        "Advanced", "Metastatic"))
auc_df$cluster <- factor(auc_df$cluster, levels = cluster_levels)

# Create and save ridge plot 
ridge_plot <- ggplot(auc_df, aes(x = AT2.cell, y = cluster, fill = cluster)) +
  ggridges::geom_density_ridges(scale = 2) +
  scale_fill_manual(values = col_cluster) +
  labs(x = "AT2 score", y = "", color = "Class") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Save plot
output_dir <- "~/figure3/"
ggsave(
  file.path(output_dir, "fig3g_AT2_ridge.pdf"),
  umap_plot,
  width = 3.6,
  height = 2.58,
  device = cairo_pdf
)

ridge_plot <- ggplot(auc_df, aes(x = AT2.cell, y = cluster, fill = cluster)) +
  ggridges::geom_density_ridges(scale = 2) +
  scale_fill_manual(values = col_cluster) +
  labs(x = "AT1 score", y = "", color = "Class") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain")
  )

# Save plot
output_dir <- "~/figure3/"
ggsave(
  file.path(output_dir, "fig3g_AT1_ridge.pdf"),
  umap_plot,
  width = 3.6,
  height = 2.58,
  device = cairo_pdf
)
############################## > Fig3h,3j,3l  ##################################
# Prepare data for Monocle analysis
Tumor_AT1_AT2_cell = merge_data[,which(merer_data$stage %in% c("Tumor cell","AT1 cell","AT2 cell"))]
count_matrix <- as(as.matrix(Tumor_AT1_AT2_cell@assays$RNA@counts), 'sparseMatrix')

# Create metadata objects
gene_metadata <- data.frame(
  gene_short_name = rownames(Tumor_AT1_AT2_cell@assays$RNA),
  row.names = rownames(Tumor_AT1_AT2_cell@assays$RNA)
)
cell_metadata <- Tumor_AT1_AT2_cell@meta.data

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

subtype_plot <- plot_cell_trajectory(
  cds,
  color_by = "subtype",
  cell_size = 1,
  show_branch_points = FALSE
) +
  facet_wrap(~subtype1, scales = "free") +
  scale_color_manual(values = c(col_AT1, col_AT2, col_tumor[1:4])) +
  labs(color = "Cell type") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain"),
    legend.position = "right"
  )

# save plot
ggsave(
  filename = "~/figure3/Fig3h_pseudotime.pdf",
  plot = pseudotime_plot,
  width = 5.6,
  height = 3.78,
  device = "pdf"
)

# Generate AT1 trajectory plot 
at1_plot <- plot_cell_trajectory(
  cds,
  color_by = "AT1",          
  cell_size = 0.6,           
  show_branch_points = FALSE  
) +
  # Define color gradient for AT1 scores
  scale_color_gradientn(
    values = seq(0, 1, 0.2),  
    colors = c("#330066", "#66CC66", "#FFCC33", "#FFCC33")  
  ) +
  labs(color = "AT1 score") + 
  # Custom theme settings
  theme(
    legend.position = "right",  
    panel.grid = element_blank(), 
    axis.text = element_text(size = 10, family = "sans"), 
    axis.title = element_text(size = 12, family = "sans")  
  )

# Save plot
ggsave(
  filename = "~/figure3/Fig3j_AT1.pdf", 
  plot = at1_plot,       
  width = 3.83,         
  height = 3.78,         
  device = "pdf",        
  units = "in",         
  useDingbats = FALSE   
)

# Generate plot object showing AT2 expression along trajectory
at2_plot <- plot_cell_trajectory(
  cds,  
  color_by = "AT2",  
  cell_size = 0.6,   
  show_branch_points = FALSE  
) +
  # Apply custom color gradient for AT2 scores
  scale_color_gradientn(
    values = seq(0, 1, 0.2),  
    colors = c("#330066", "#336699", "#66CC66", "#FFCC33") 
  ) +
  labs(color = "AT2 score") +  
  # Customize plot appearance
  theme(
    panel.grid = element_blank(),  
    axis.text = element_text(size = 10, family = "sans"),  
    axis.title = element_text(size = 12, family = "sans"),  
    legend.position = "right"  
  )

# Save plot
ggsave(
  filename = "~/figure3/Fig3j_AT2.pdf", 
  plot = at2_plot,       
  width = 3.83,         
  height = 3.78,        
  device = "pdf",      
  units = "in",         
  useDingbats = FALSE    
)

subtype_plot <- plot_cell_trajectory(
  cds,
  color_by = "subtype",
  cell_size = 1,
  show_branch_points = FALSE
) +
  facet_wrap(~subtype1, scales = "free") +
  scale_color_manual(values = c(col_AT1, col_AT2, col_tumor[1:4])) +
  labs(color = "Cell type") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 12, family = "sans", face = "plain"),
    axis.title = element_text(size = 14, family = "sans", face = "plain"),
    legend.text = element_text(size = 12, family = "sans", face = "plain"),
    legend.title = element_text(size = 12, family = "sans", face = "plain"),
    legend.position = "right"
  )

# save plot
ggsave(
  filename = "~/figure3/Fig3l_subtype.pdf",
  plot = subtype_plot,
  width = 5.6,       
  height = 3.78,    
  device = "pdf",    
  units = "in",      
  dpi = 300         
)

############################## > Fig3k  ##################################
LUAD_mouse= readRDS("~/data/LUAD_mouse_Marj.rds")

calculateDiversityScores <- function(seurat_object, 
                                     group_column, 
                                     n_cells = 50, 
                                     n_repeats = 100,
                                     removeOutlier = TRUE,
                                     nsd = 3,
                                     topPCs = 1:30) {
  # Initialize results dataframe
  results <- data.frame(
    Iteration = integer(),
    Group = character(),
    DiversityScore = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get PCA space (first 30 PCs by default)
  pca_embeddings <- seurat_object@reductions$pca@cell.embeddings[, topPCs]
  groups <- unique(seurat_object@meta.data[[group_column]])
  
  # Perform repeated sampling
  for (i in seq_len(n_repeats)) {
    message("Processing iteration ", i, " of ", n_repeats)
    
    for (group in groups) {
      # Get cells for current group
      group_cells <- rownames(seurat_object@meta.data)[
        seurat_object@meta.data[[group_column]] == group
      ]
      
      # Skip if not enough cells
      if (length(group_cells) < n_cells) {
        message("Skipping group ", group, " - insufficient cells")
        next
      }
      
      # Random sampling and diversity calculation
      sampled_cells <- sample(group_cells, n_cells)
      sample_pca <- pca_embeddings[sampled_cells, , drop = FALSE]
      sample_groups <- seurat_object@meta.data[sampled_cells, group_column]
      
      # Calculate diversity score
      div_score <- mean(
        calDiversityScore(sample_pca, sample_groups, removeOutlier, nsd),
        na.rm = TRUE
      )
      
      # Store results
      results <- rbind(results, data.frame(
        Iteration = i,
        Group = group,
        DiversityScore = div_score
      ))
    }
  }
  
  return(results)
}

# Calculate diversity scores
diversity_scores <- calculateDiversityScores(
  seurat_object = LUAD_mouse,
  group_column = "mouseID",  
  n_cells = 30,
  n_repeats = 100
)

# Prepare data for analysis
Mouse_data <- data.frame(
  AT2 = LUAD_mm$AT2_score,AT1 = LUAD_mm$AT1_score,Mouse = LUAD_mm$mouseID,Stage = LUAD_mm$Stage)%>% 
  group_by(Mouse) %>%
  summarise(AT2_mean = mean(AT2),AT1_mean = mean(AT1),Stage = first(Stage)) %>% 
  left_join(diversity_scores %>% rename(Mouse = Patient), by = "Mouse") %>%
  filter(Stage %in% c("KPT 12w", "KPT 2w", "KPT 20w", "KPT 30w"))

# Create AT1 diversity plot
p_AT1 <- ggplot(Mouse_data, aes(x = AT1_mean, y = diversity_scores)) +
  geom_point(size = 2, alpha = 0.3, color = "#6baed6") +
  geom_smooth(method = 'lm', formula = y ~ x, 
              color = "#756bb1", fill = "#cbc9e2") +
  labs(x = "AT1 score", y = "Transcriptional diversity") +
  stat_cor(size = 5, label.y = 0.6) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.position = 'none'
  )

# Save plot
ggsave("~/figure3/Fig3k_AT1.pdf", p_AT1, 
       width = 3.45, height = 3.45, bg = "transparent")

# Create AT2 diversity plot
p_AT2 <- ggplot(Mouse_data, aes(x = AT2_mean, y = diversity_scores)) +
  geom_point(size = 2, alpha = 0.3, color = "#6baed6") +
  geom_smooth(method = 'lm', formula = y ~ x, 
              color = "#756bb1", fill = "#cbc9e2") +
  labs(x = "AT2 score", y = "Transcriptional diversity") +
  stat_cor(size = 5, label.y = 0.6) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, family = "sans"),
    legend.position = 'none'
  )

# Save plot
ggsave("~/figure3/Fig3k_AT2.pdf", p_AT2, 
       width = 3.45, height = 3.45, bg = "transparent")

############################## > Fig3m  ##################################
# Create dot plot
p <- DotPlot(tumor, 
             features = c("CDKN1A", "CLDN4", "KRT8", "NDRG1", "SFN"),
             group.by = "subtype") +
  coord_flip() +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colors = c("#330066", "#336699", "#66CC66", "#FFCC33")
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    color = guide_colorbar(title = "Average expression"),
    size = guide_legend(title = "Percent expressed")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(
      angle = 50, hjust = 1, 
      size = 14, family = "sans", color = "black"
    ),
    axis.text.y = element_text(
      size = 14, family = "sans", 
      face = "italic", color = "black"
    ),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans")
  )

# Save plot
ggsave("~/Fig3/figure3m.pdf", p, width = 4.11, height = 3.48, bg = "transparent")



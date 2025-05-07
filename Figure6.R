library(Deseq2)
library(ggplot2)
library(Seurat)
library(GseaVis)
library(ggrepel)
############################## > Fig6i  ##################################
# Load and prepare data
RNA_result <- read.csv("~/data/A549p2_PCDH PK A549p2_P_KLF4.csv") %>%
  arrange(log2FoldChange)

Tumor_cell <- readRDS("~/data/Tumor_cell.rds")
Idents(Tumor_cell) <- "subtype"
Deg <- FindAllMarkers(Tumor_cell, logfc.threshold = 0.25) %>%
  filter(avg_log2FC > 1, cluster == "C2")

# Prepare gene list for coloring
diff_gene <- RNA_result %>%
  mutate(
    Diff = case_when(
      padj < 0.5 & abs(log2FoldChange) >= 0.5 & log2FoldChange >= 0.5 ~ "Up",
      padj < 0.5 & abs(log2FoldChange) >= 0.5 ~ "Down",
      TRUE ~ "Not"
    ),
    p.value = pmin(-log10(padj), 300)  # Cap p-values at 300
  )

# Select significant genes
select_genes <- diff_gene %>%
  filter(Gene %in% Deg$gene, log2FoldChange > 0.5)

# Create volcano plot
p <- ggplot(diff_gene, aes(log2FoldChange, p.value, color = Diff)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_point(
    data = select_genes, shape = 1, size = 1, stroke = 1, color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = select_genes,aes(label = Gene),color = "black",
    size = 4,box.padding = 0.35,point.padding = 0.3,show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("red", "blue", "grey"),
    limits = c("Up", "Down", "Not"),
    guide = "none"
  ) +
  labs(
    x = bquote(~Log[2]~"(fold change)"),
    y = bquote(~-Log[10]~italic("P-value"))
  ) +
  coord_cartesian(ylim = c(0, 320), xlim = c(-4, 4)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "longdash", color = "grey40") +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey40") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Save plot
ggsave("~/figure6/Fig6i.pdf", p, width = 3.7, height = 5.8)

############################## > Fig6j,k  ##################################
# Filter DEGs and prepare geneset
geneset <- Deg %>% 
  filter(avg_log2FC > 1, cluster == "C2") %>% 
  transmute(gs_name = "Tumor_clusters2_score", gene_symbol = gene)

# Run GSEA
egmt <- GSEA(
  geneList = geneList,TERM2GENE = geneset,
  minGSSize = 1,pvalueCutoff = 0.99,verbose = FALSE
)

# Create and save GSEA plot
pdf("~/figure6/Fig6j.pdf", width = 3.8, height = 3.6)
gseaNb(
  object = egmt,geneSetID = "Tumor_clusters2_score",
  subPlot = 2,addPval = TRUE,pvalX = 0.6,
  pvalY = 0.65,pCol = 'black',pHjust = 0,
  newHtCol = c("blue", "white", "red")
)
dev.off()

### AT1
epi_markers <- read.csv("~/data/Epi_markers_from_normal_lung.csv")
gene_lists <- epi_markers %>%
  dplyr::select(gene, celltype) %>%
  split(.$celltype) %>%
  lapply(function(x) x[[1]])

AT1_markers=data.frame(data.frame(gene_symbol=gene_lists[["AT1 cell"]],gs_name=rep("AT1 score",length(geneset[[2]]))))
AT1_markers =AT1_markers %>% dplyr::select("gs_name", everything())

egmt <- clusterProfiler::GSEA(geneList, TERM2GENE=AT1_markers, 
                              minGSSize = 1,
                              pvalueCutoff = 0.99,
                              verbose=FALSE)
display_geneset="AT1 score"
pdf(file="~/figure6/Fig6k_AT1.pdf",width = 3.8,height = 3.6)
gseaNb(object =egmt,
       geneSetID =display_geneset,
       subPlot = c(2),
       addPval = T,
       pvalX = 0.5,pvalY = 0.65,
       pCol = 'black',
       pHjust = 0,
       newHtCol = c("blue","white", "red"))
dev.off()

AT2_markers=data.frame(data.frame(gene_symbol=gene_lists[["AT2 cell"]],gs_name=rep("AT2 score",length(geneset[[2]]))))
AT2_markers =AT2_markers %>% dplyr::select("gs_name", everything())

egmt <- clusterProfiler::GSEA(geneList, TERM2GENE=AT2_markers, 
                              minGSSize = 1,
                              pvalueCutoff = 0.99,
                              verbose=FALSE)
display_geneset="AT2 score"
pdf(file="~/figure6/Fig6k_AT2.pdf",width = 3.8,height = 3.6)
gseaNb(object =egmt,
       geneSetID =display_geneset,
       subPlot = c(2),
       addPval = T,
       pvalX = 0.5,pvalY = 0.65,
       pCol = 'black',
       pHjust = 0,
       newHtCol = c("blue","white", "red"))
dev.off()


#!/usr/bin/env Rscript

# =============================================================================
# Author: Oliver Abinader
# Description:
#   Additional downstream analyses:
#   - Cell composition by sample and cluster
#   - Per-cluster SingleR cell type composition
#   - Dominant cell type assignment per cluster
#   - Optional exploratory gene-gene correlation analysis
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(writexl)
  library(grid)
})

# -----------------------------
# Paths
# -----------------------------
input_rds <- "/path/to/analysis_folder/de_and_annotation/integratedSeuratObject_celltypeIdentification.rds"
output_dir <- "/path/to/analysis_folder/exploratory"
composition_dir <- file.path(output_dir, "cluster_composition")
correlation_dir <- file.path(output_dir, "gene_correlation")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(composition_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
merged_seurat_filtered <- readRDS(input_rds)
meta <- merged_seurat_filtered@meta.data

# -----------------------------
# Cell proportions per sample within each cluster
# -----------------------------
cell_counts <- meta %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(Cell_Count = n(), .groups = "drop")

total_cells <- meta %>%
  group_by(orig.ident) %>%
  summarise(Total_Cells = n(), .groups = "drop")

cell_percentages <- cell_counts %>%
  left_join(total_cells, by = "orig.ident") %>%
  mutate(Percent = (Cell_Count / Total_Cells) * 100)

write.xlsx(
  x = cell_percentages,
  file = file.path(composition_dir, "cell_percentages_by_sample_and_cluster.xlsx"),
  rowNames = FALSE
)

# Plot one barplot per cluster
all_clusters <- sort(unique(cell_percentages$seurat_clusters))

for (cluster_of_interest in all_clusters) {
  cluster_data <- cell_percentages %>%
    filter(seurat_clusters == cluster_of_interest)

  g <- ggplot(cluster_data, aes(x = orig.ident, y = Percent, fill = orig.ident)) +
    geom_bar(stat = "identity", color = "black") +
    labs(
      title = paste("Cluster", cluster_of_interest, "- % of Cells per Sample"),
      x = "Sample",
      y = "Percentage of Cells (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white"),
      legend.position = "none"
    )

  ggsave(
    filename = file.path(
      composition_dir,
      paste0("cluster_", cluster_of_interest, "_percentage_plot.png")
    ),
    plot = g,
    device = "tiff",
    width = 8,
    height = 8,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )
}

# -----------------------------
# Per-cluster cell type percentages (SingleR fine labels)
# -----------------------------
if ("singleR_label.fine" %in% colnames(meta)) {
  cluster_celltype_percentage <- meta %>%
    group_by(seurat_clusters, singleR_label.fine) %>%
    summarise(cell_count = n(), .groups = "drop_last") %>%
    mutate(cluster_total = sum(cell_count)) %>%
    ungroup() %>%
    mutate(percentage = (cell_count / cluster_total) * 100)

  write.xlsx(
    x = as.data.frame(cluster_celltype_percentage),
    file = file.path(composition_dir, "cluster_celltype_percentage.xlsx"),
    rowNames = FALSE
  )

  # Dominant cell type per cluster
  top_celltype_per_cluster <- cluster_celltype_percentage %>%
    group_by(seurat_clusters) %>%
    slice_max(order_by = percentage, n = 1, with_ties = FALSE) %>%
    ungroup()

  dominant_celltype_per_cluster <- top_celltype_per_cluster[, c("seurat_clusters", "singleR_label.fine")]

  write.xlsx(
    x = dominant_celltype_per_cluster,
    file = file.path(composition_dir, "dominant_celltype_per_cluster.xlsx"),
    rowNames = FALSE
  )

  # Add dominant cell type label to each cell, should be in the metadata (optional)
  meta$dominant_celltype_per_cluster <- dominant_celltype_per_cluster$singleR_label.fine[
    match(meta$seurat_clusters, dominant_celltype_per_cluster$seurat_clusters)
  ]

  # Save UMAP by dominant cluster label
  p_dominant <- DimPlot(
    merged_seurat_filtered,
    reduction = "umap.harmony",
    group.by = "dominant_celltype_per_cluster",
    label = TRUE,
    repel = TRUE,
    label.size = 4
  ) + ggtitle("UMAP by Dominant Cell Type per Cluster")

  ggsave(
    filename = file.path(composition_dir, "UMAP_dominant_celltype_per_cluster.png"),
    plot = p_dominant,
    device = "tiff",
    width = 8,
    height = 8,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )

  # Stacked barplot
  g <- ggplot(
    cluster_celltype_percentage,
    aes(x = factor(seurat_clusters), y = percentage, fill = singleR_label.fine)
  ) +
    geom_bar(stat = "identity") +
    scale_x_discrete(name = "Seurat Cluster") +
    scale_y_continuous(name = "Percentage") +
    labs(fill = NULL) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )

  ggsave(
    filename = file.path(composition_dir, "cluster_celltype_percentage.png"),
    plot = g,
    device = "tiff",
    width = 8,
    height = 8,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )
}

# -----------------------------
# OPTIONAL: exploratory gene-gene correlation analysis
# -----------------------------
# This section is intentionally left as a template because gene sets and cluster groupings are project-specific.

run_gene_correlation <- FALSE

if (run_gene_correlation) {
  # Example setup:
  expr_data <- as.matrix(GetAssayData(merged_seurat_filtered, assay = "SCT", layer = "data"))
  # Note: if this series of code run out of error, make sure the expression data is in a dataframe format.

  # Define a "cluster group of interest", as an example: consider clusters 1,3,4,6,7 as one cluster group. 
  # This is from a typical single-cell clustering project that had like 12 clusters in total.
  set1_clusters <- c(1, 3, 4, 6, 7) # In this case, we picked these clusters as one "cluster group of interest". 
  # Here, we are not doing this cluster group vs. the remaining cells in all toher clusters.

  selected_cells <- rownames(meta)[meta$seurat_clusters %in% set1_clusters]
  set1_data <- expr_data[, colnames(expr_data) %in% selected_cells, drop = FALSE]

  # Replace these with project-specific genes
  genes_1 <- c("GENE_A", "GENE_B", "GENE_C")
  genes_2 <- c("GENE_X", "GENE_Y", "GENE_Z")

  result <- data.frame(
    Gene1 = character(),
    Gene2 = character(),
    Correlation = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )

  for (gene1 in genes_1) {
    for (gene2 in genes_2) {
      if (gene1 %in% rownames(set1_data) && gene2 %in% rownames(set1_data)) {
        gene1_data <- as.numeric(set1_data[gene1, ])
        gene2_data <- as.numeric(set1_data[gene2, ])

        res <- cor.test(gene1_data, gene2_data, method = "spearman", exact = FALSE)

        result <- rbind(
          result,
          data.frame(
            Gene1 = gene1,
            Gene2 = gene2,
            Correlation = unname(res$estimate),
            P_Value = res$p.value
          )
        )
      }
    }
  }

  result$Adjusted_P_Value <- p.adjust(result$P_Value, method = "BH", n = nrow(result))
  result <- result[order(-result$Correlation), ]

  write_xlsx(
    x = result,
    path = file.path(correlation_dir, "gene_correlation_results.xlsx")
  )
}

# Save updated object if dominant cluster labels were added to the metadata
saveRDS(
  merged_seurat_filtered,
  file = file.path(output_dir, "integratedSeuratObject_exploratory.rds")
)

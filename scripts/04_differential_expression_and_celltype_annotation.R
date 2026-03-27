#!/usr/bin/env Rscript

# =============================================================================
# Author: Oliver Abinader
# Description:
#   Downstream analysis on integrated Seurat object:
#   - Differential expression FindMarkers (MS vs HC by cluster)
#   - Global differential expression FindAllMarkers (MS vs HC)
#   - Cluster marker detection (FindAllMarkers) ??
#   - Gene annotation (full gene names + summaries)
#   - Cell type annotation using SingleR + celldex DICE reference as an example
#   - UMAP export with SingleR labels
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(mygene)
  library(SingleR)
  library(celldex)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
input_rds <- "/path/to/analysis_folder/integration/integratedSeuratObject_SCTharmony.rds"
output_dir <- "/path/to/analysis_folder/de_and_annotation"
de_dir <- file.path(output_dir, "de_results")
annot_dir <- file.path(output_dir, "celltype_annotation")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
merged_seurat_filtered <- readRDS(input_rds)

# Ensure harmony clusters are active for cluster-based analyses
if (!"harmony_clusters" %in% colnames(merged_seurat_filtered@meta.data)) {
  stop("harmony_clusters column not found in metadata.")
}
# But to actually check if it is active, do that: Idents(merged_seurat_filtered)

# -----------------------------
# Prepare for DE
# -----------------------------
merged_seurat_filtered <- PrepSCTFindMarkers(merged_seurat_filtered)

# -----------------------------
# DE: MS vs HC within each (Harmony) cluster
# -----------------------------
merged_seurat_filtered$Disease_Cluster <- paste0(
  merged_seurat_filtered$Disease_status,
  "_cluster",
  merged_seurat_filtered$harmony_clusters
)

Idents(merged_seurat_filtered) <- "Disease_Cluster"

clusters <- sort(unique(merged_seurat_filtered$seurat_clusters))
# Check that seurat_clusters is actually same as harmony_clusters

for (clust in clusters) {
  ident_1 <- paste0("MS_cluster", clust)
  ident_2 <- paste0("HC_cluster", clust)

  if (ident_1 %in% levels(Idents(merged_seurat_filtered)) &&
      ident_2 %in% levels(Idents(merged_seurat_filtered))) {

    de_result <- FindMarkers(
      merged_seurat_filtered,
      ident.1 = ident_1,
      ident.2 = ident_2,
      verbose = FALSE
    )

    de_result$gene <- rownames(de_result)

    write.xlsx(
      de_result,
      file = file.path(de_dir, paste0("cluster", clust, "_MSvsHC.xlsx")),
      rowNames = FALSE
    )
  }
}

# -----------------------------
# Global DE: MS vs HC
# -----------------------------
Idents(merged_seurat_filtered) <- "Disease_status"

de_global <- FindMarkers(
  merged_seurat_filtered,
  ident.1 = "MS",
  ident.2 = "HC",
  verbose = FALSE
)

de_global$gene <- rownames(de_global)

write.xlsx(
  de_global,
  file = file.path(de_dir, "AllClusters_MSvsHC.xlsx"),
  rowNames = FALSE
)

# -----------------------------
# Find all cluster markers
# -----------------------------
Idents(merged_seurat_filtered) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  merged_seurat_filtered,
  verbose = FALSE
)

# Add full gene names
all_markers$Full_Gene_Name <- mapIds(
  org.Hs.eg.db,
  keys = all_markers$gene,
  column = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Add gene summaries from MyGene.info
genes <- unique(all_markers$gene)

gene_info <- queryMany(
  genes,
  scopes = "symbol",
  fields = c("summary"),
  species = "human"
)

gene_info <- as.data.frame(gene_info)

if (all(c("query", "summary") %in% colnames(gene_info))) {
  gene_info <- gene_info[, c("query", "summary")]
  colnames(gene_info) <- c("gene", "RefSeq_Gene_Summary")

  all_markers <- left_join(all_markers, gene_info, by = "gene")
}

write.xlsx(
  x = all_markers,
  file = file.path(de_dir, "DE_results_allmarkers_annotated.xlsx"),
  rowNames = FALSE
)

# -----------------------------
# SingleR cell type annotation
# -----------------------------
ref <- celldex::DatabaseImmuneCellExpressionData() # Example: DICE reference dataset

# Use RNA counts for SingleR
test_matrix <- GetAssayData(
  object = merged_seurat_filtered,
  assay = "SCT", # SCT assay
  layer = "counts"
)

data.pred.main <- SingleR(
  test = test_matrix,
  ref = ref,
  assay.type.test = 1, # Use the count data and not the normalized data
  labels = ref$label.main # Gives broad/ higher-level cell types
)

data.pred.fine <- SingleR(
  test = test_matrix,
  ref = ref,
  assay.type.test = 1,
  labels = ref$label.fine # Gives more granular/ subtype-level annotations
)

# Confirm row order
stopifnot(all(rownames(data.pred.main) == colnames(merged_seurat_filtered)))
stopifnot(all(rownames(data.pred.fine) == colnames(merged_seurat_filtered)))

# Add labels to metadata
merged_seurat_filtered$singleR_label.main <- data.pred.main$labels
merged_seurat_filtered$singleR_label.fine <- data.pred.fine$labels

# Save UMAP with fine labels # For more granular annotations
p_singleR <- DimPlot(
  object = merged_seurat_filtered,
  reduction = "umap.harmony",
  group.by = "singleR_label.fine",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5
) + ggtitle("Integrated UMAP with SingleR Fine Labels")

ggsave(
  filename = file.path(annot_dir, "IntegratedUMAP_singleR_label_fine.png"),
  plot = p_singleR,
  device = "tiff",
  width = 8,
  height = 8,
  units = "in",
  dpi = 350,
  compression = "lzw"
)

# Save annotated object
saveRDS(
  merged_seurat_filtered,
  file = file.path(output_dir, "integratedSeuratObject_celltypeIdentification.rds")
)

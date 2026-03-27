#!/usr/bin/env Rscript

# =============================================================================
# Author: Oliver Abinader
#
# Description:
#   Reproducible scRNA-seq downstream preprocessing workflow in Seurat:
#   - Load Cell Ranger feature-barcode matrices from multiple samples
#   - Create and merge Seurat objects
#   - Add sample-level metadata
#   - Compute QC metrics
#   - Filter low-quality cells
#   - Perform SCTransform normalization
#   - Run PCA, clustering, and UMAP (unintegrated)
#   - Perform Harmony-based integration
#   - Run clustering and UMAP (integrated)
#   - Export QC plots and save integrated Seurat object
#
# Input:
#   Parent directory containing multiple sample folders, each with:
#     - matrix.mtx.gz
#     - features.tsv.gz
#     - barcodes.tsv.gz
#
# Output:
#   - QC violin plots
#   - Integrated UMAP plot
#   - Integrated Seurat object (.rds)
#   - Summary tables for cells per sample and cells per cluster
#
# Notes:
#   - Designed for Seurat v5-style workflows using IntegrateLayers + Harmony.
#   - Update sample names in the merge section to match your dataset.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(openxlsx)
})

# -----------------------------
# User-defined paths
# -----------------------------
input_dir <- "."
# NOTE: "." represents the current working directory.
# Example: it should be set to the parent directory that contains all sample folders, where each of the latter includes files like matrix, barcode, feature.

output_dir <- "/path/to/analysis_folder/integration"
qc_dir <- file.path(output_dir, "qc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Identify all sample directories containing Cell Ranger matrices
# -----------------------------
all_dirs <- list.dirs(path = input_dir, recursive = TRUE, full.names = TRUE)

sample_dirs <- all_dirs[
  sapply(all_dirs, function(x) file.exists(file.path(x, "matrix.mtx.gz")))
]

if (length(sample_dirs) == 0) {
  stop("No sample directories containing matrix.mtx.gz were found.")
}

message("Found ", length(sample_dirs), " sample directories.")

# -----------------------------
# Read all samples into Seurat objects
# -----------------------------
seurat_list <- list()

for (x in sample_dirs) {
  sample_name <- basename(x)
  message("Reading sample: ", sample_name)

  counts <- ReadMtx(
    mtx = file.path(x, "matrix.mtx.gz"),
    features = file.path(x, "features.tsv.gz"),
    cells = file.path(x, "barcodes.tsv.gz")
  )

  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3
  )

  seurat_list[[sample_name]] <- seurat_obj
}

# -----------------------------
# Merge all Seurat objects
# -----------------------------
# Generic merge approach for any number of samples
sample_names <- names(seurat_list)

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names
)

# -----------------------------
# Add metadata: You may add new columns describing more about the cells/ samples.
# -----------------------------
# Example 1: we are adding Disease_status column, which can be inferred from sample names starting with "MS"
merged_seurat$Disease_status <- ifelse(
  str_detect(merged_seurat$orig.ident, "^MS"),
  "MS",
  "HC"
)

# Example 2: Age: extract numeric value preceding "yo" if present
merged_seurat$Age <- as.numeric(
  str_extract(merged_seurat$orig.ident, "(\\d+)(?=yo)")
)

# Example 3: Extract barcode portion after the first underscore in merged cell names
merged_seurat$Barcode <- sapply(
  strsplit(colnames(merged_seurat), "_"),
  function(x) paste(x[-1], collapse = "_")
)

# Example 4: If we are intersted in adding raw counts of a gene called "FOXP3" to metadata by first checking if the gene is present
if ("FOXP3" %in% rownames(merged_seurat)) {
  merged_seurat$FOXP3_counts <- FetchData(
    merged_seurat,
    vars = "FOXP3",
    assay = "RNA",
    layer = "counts"
  )$FOXP3
} else {
  warning("FOXP3 not found in dataset. FOXP3_counts metadata column will not be added.")
  merged_seurat$FOXP3_counts <- NA
}

# -----------------------------
# QC metrics
# -----------------------------
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^MT-"
)

# QC violin plots (all samples)
p0 <- VlnPlot(
  merged_seurat,
  features = "nFeature_RNA",
  assay = "RNA",
  layer = "counts",
  raster = FALSE
) +
  NoLegend() +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  filename = file.path(qc_dir, "VlnPlot_nFeature_RNA.png"),
  plot = p0,
  device = "tiff",
  width = 8,
  height = 8, 
  units = "in",
  dpi = 350,
  compression = "lzw"
) 

p1 <- VlnPlot(
  merged_seurat,
  features = "nCount_RNA",
  assay = "RNA",
  layer = "counts",
  raster = FALSE
) +
  NoLegend() +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  filename = file.path(qc_dir, "VlnPlot_nCount_RNA.png"),
  plot = p1,
  device = "tiff",
  width = 8,
  height = 8, 
  units = "in",
  dpi = 350,
  compression = "lzw"
) 
         
p2 <- VlnPlot(
  merged_seurat,
  features = "percent.mt",
  assay = "RNA",
  layer = "counts",
  raster = FALSE
) +
  NoLegend() +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  filename = file.path(qc_dir, "VlnPlot_percent_mt.png"),
  plot = p2,
  device = "tiff",
  width = 8,
  height = 8, 
  units = "in",
  dpi = 350,
  compression = "lzw"
) 

# -----------------------------
# Filtering low-quality cells
# -----------------------------
# Default filtering:
#   - nFeature_RNA > 200
#   - nFeature_RNA < 6000
#   - nCount_RNA < 40000
#   - percent.mt < 15
#   - FOXP3_counts > 0 (only if FOXP3 exists). 
#
# Note:
# FOXP3 filtering is biologically specific and may not apply to all datasets.
# Keep or remove depending on study design.
# Also, note that the filtering numbers may differ from a study to another, based on the QC plots and the input from the experimentalist. 

if ("FOXP3" %in% rownames(merged_seurat)) {
  merged_seurat_filtered <- subset(
    merged_seurat,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 6000 &
      nCount_RNA < 40000 &
      percent.mt < 15 &
      FOXP3_counts > 0
  )
} else {
  merged_seurat_filtered <- subset(
    merged_seurat,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 6000 &
      nCount_RNA < 40000 &
      percent.mt < 15
  )
}

message("Cells retained after filtering: ", ncol(merged_seurat_filtered))

# -----------------------------
# Unintegrated workflow
# -----------------------------
merged_seurat_filtered <- SCTransform(
  merged_seurat_filtered,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

merged_seurat_filtered <- RunPCA(
  merged_seurat_filtered,
  verbose = FALSE
) # Reduction name is PCA by default

#Check the Elbow plot to find the "elbow"—the point where variance tapers off-

merged_seurat_filtered <- FindNeighbors(
  merged_seurat_filtered,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

merged_seurat_filtered <- FindClusters(
  merged_seurat_filtered,
  resolution = 0.5,
  cluster.name = "unintegrated_clusters",
  verbose = FALSE
)

merged_seurat_filtered <- RunUMAP(
  merged_seurat_filtered,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap.unintegrated",
  verbose = FALSE
)

# -----------------------------
# Harmony integration (Seurat v5)
# -----------------------------
merged_seurat_filtered <- IntegrateLayers(
  object = merged_seurat_filtered,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  verbose = FALSE
)

merged_seurat_filtered <- FindNeighbors(
  merged_seurat_filtered,
  dims = 1:30,
  reduction = "harmony",
  verbose = FALSE
)

merged_seurat_filtered <- FindClusters(
  merged_seurat_filtered,
  resolution = 0.5,
  cluster.name = "harmony_clusters",
  verbose = FALSE
)

merged_seurat_filtered <- RunUMAP(
  merged_seurat_filtered,
  dims = 1:30,
  reduction = "harmony",
  reduction.name = "umap.harmony",
  verbose = FALSE
)

message("Number of Seurat clusters: ",
        length(unique(merged_seurat_filtered$seurat_clusters)))
# Note that the clusters in harmony_clusters will be refelcted in the seurat_clusters column.
         
# -----------------------------
# Save integrated UMAP
# -----------------------------
p_umap_sample <- DimPlot(
  object = merged_seurat_filtered,
  reduction = "umap.harmony",
  group.by = "orig.ident",
  label = FALSE
) + ggtitle("Integrated UMAP by Sample")
# You can also generate other UMAP(s) by selecting appropriate variable(s) column from the metadata.

ggsave(
  filename = file.path(output_dir, "IntegratedUMAP_bySample.png"),
  plot = p_umap_sample,
  device = "tiff",
  width = 8,
  height = 8, 
  units = "in",
  dpi = 350,
  compression = "lzw"
) 

# -----------------------------
# Summary tables
# -----------------------------
cells_per_sample <- merged_seurat_filtered@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  arrange(orig.ident)

write.xlsx(
  x = cells_per_sample,
  file = file.path(output_dir, "cells_per_sample.xlsx"),
  rowNames = FALSE
)

cells_per_cluster <- merged_seurat_filtered@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  arrange(seurat_clusters)

write.xlsx(
  x = cells_per_cluster,
  file = file.path(output_dir, "cells_per_seurat_cluster.xlsx"),
  rowNames = FALSE
)

# -----------------------------
# Save Seurat object
# -----------------------------
saveRDS(
  merged_seurat_filtered,
  file = file.path(output_dir, "integratedSeuratObject_SCTharmony.rds")
)

message("Pipeline complete.")

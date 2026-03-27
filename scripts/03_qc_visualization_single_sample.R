#!/usr/bin/env Rscript

# =============================================================================
# Author: Oliver Abinader
# Description:
#   Create publication-style QC violin plots for a single Seurat object.
#   This script is intended for manual / sample-specific visualization and is
#   not yet automated across multiple samples.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
})

# -----------------------------
# User inputs
# -----------------------------
input_rds <- "/path/to/analysis_folder/integration/integratedSeuratObject_SCTharmony.rds"
output_dir <- "/path/to/analysis_folder/qc_single_sample"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load object
seurat_obj <- readRDS(input_rds)

# OPTIONAL:
# Replace this with a subset if you want one specific sample only:
# seurat_obj <- subset(seurat_obj, subset = orig.ident == "YourSampleName")

rna_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
y_labels <- c(
  "No. of Unique Genes Detected per Cell",
  "No. of Reads Detected per Cell",
  "Percent of Mitochondrial Reads (%)"
)

p_list <- VlnPlot(
  seurat_obj,
  features = rna_features,
  assay = "RNA",
  ncol = 3,
  layer = "counts"
)

p_list <- lapply(seq_along(rna_features), function(i) {
  p_list[[i]] +
    ylab(y_labels[i]) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
})

p_combined <- wrap_plots(p_list, ncol = 3)

ggsave(
  filename = file.path(output_dir, "VlnPlot_RNA_single_sample.png"),
  plot = p_combined,
  device = "tiff",
  width = 8,
  height = 8,
  units = "in",
  dpi = 350,
  compression = "lzw"
)

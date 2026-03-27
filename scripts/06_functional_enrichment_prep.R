# Author: Oliver Abinader

# Prepare significantly upregulated and downregulated gene lists from differential expression result files for downstream functional enrichment analysis (e.g., WebGestalt GO BP).

library(readxl)
library(openxlsx)
library(dplyr)

# -----------------------------
# User-defined paths
# -----------------------------

# Folder containing DE result Excel files
# Example: cluster0_MSvsHC.xlsx, cluster1_MSvsHC.xlsx, etc.
base_dir <- "/path/to/analysis_folder/de_and_annotation/de_results"

# Output folder for GO enrichment input files
output_dir <- "/path/to/analysis_folder/de_and_annotation/functional_enrichment_prep"

# Create output folder if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# Customzied Parameters
# -----------------------------

# Thresholds for defining significantly regulated genes
fdr_cutoff <- 0.05
log2fc_up_cutoff <- 0.1
log2fc_down_cutoff <- -0.1

# -----------------------------
# Find DE result files
# -----------------------------

files <- list.files(base_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .xlsx files found in base_dir. Check the DE results folder path.")
}

# -----------------------------
# Process each DE file
# -----------------------------

for (file in files) {

  # Extract file name without extension
  file_name <- tools::file_path_sans_ext(basename(file))

  # Read DE results
  df <- read_excel(file)

  # Check required columns
  required_cols <- c("p_val_adj", "avg_log2FC") # You have to make sure that these columns are numeric.
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "Skipping file: ", file_name,
        " | Missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
    next
  }

  # Filter significantly upregulated genes
  df_up <- df %>%
    filter(p_val_adj < fdr_cutoff & avg_log2FC >= log2fc_up_cutoff)

  # Filter significantly downregulated genes
  df_down <- df %>%
    filter(p_val_adj < fdr_cutoff & avg_log2FC <= log2fc_down_cutoff)

  # Write full filtered tables (keeps all DE columns)
  write.xlsx(
    df_up,
    file.path(output_dir, paste0(file_name, "_up.xlsx")),
    rowNames = FALSE
  )

  write.xlsx(
    df_down,
    file.path(output_dir, paste0(file_name, "_down.xlsx")),
    rowNames = FALSE
  )

  # Optional: write gene-symbol-only text files for easy copy/paste into WebGestalt
  # Assumes gene symbols are stored either in a "gene" column or rownames were exported as "Symbol"
  gene_col <- NULL
  if ("gene" %in% colnames(df)) {
    gene_col <- "gene"
  } else if ("Symbol" %in% colnames(df)) {
    gene_col <- "Symbol"
  }

  if (!is.null(gene_col)) {
    write.table(
      df_up[[gene_col]],
      file = file.path(output_dir, paste0(file_name, "_up_genes.txt")),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    write.table(
      df_down[[gene_col]],
      file = file.path(output_dir, paste0(file_name, "_down_genes.txt")),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  } else {
    warning(
      paste0(
        "No 'gene' or 'Symbol' column found in ", file_name,
        ". Wrote filtered .xlsx files only."
      )
    )
  }

  message(
    paste0(
      "Processed: ", file_name,
      " | Up: ", nrow(df_up),
      " | Down: ", nrow(df_down)
    )
  )
}

# -----------------------------
# Notes for downstream enrichment
# -----------------------------
# Example WebGestalt settings used in this project:
# - Analysis type: Over-Representation Analysis (ORA)
# - Organism: Homo sapiens
# - Functional database: Gene Ontology - Biological Process (GO BP)
# - Input: Gene symbols
# - Reference set: genome protein-coding
# - FDR cutoff: 5%

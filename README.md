# scRNAseq-analysis

## Author
Oliver Abinader

---

# Overview

This repository contains a **reproducible single-cell RNA-seq (scRNA-seq) workflow** based on **10x Genomics data and Seurat**.

The pipeline spans:

- Raw sequencing processing (BCL/CBCL → FASTQ)
- Cell Ranger processing to generate feature-barcode matrices
- Quality control and preprocessing in Seurat
- Downstream analysis including clustering, differential expression, annotation, and enrichment

This workflow is designed for **cancer and tumor microenvironment studies**, but is modular and adaptable to other scRNA-seq applications.


# Repository Structure

scripts/
├── qc/
│ └── check_filtered_matrix_md5.sh
├── 01_bcl_to_fastq_to_feature_matrix.sh
├── 02_seurat_preprocessing_integration.R
├── 03_qc_visualization_single_sample.R
├── 04_differential_expression_and_celltype_annotation.R
├── 05_cluster_composition_and_optional_correlation.R
└── 06_functional_enrichment_prep.R


# Workflow Overview

## Step 1 — BCL/CBCL → FASTQ

Raw 10x Genomics sequencing data is converted into FASTQ files using:

- `cellranger mkfastq`

## Step 2 — FASTQ organization

FASTQ files are separated into **sample-specific directories** to ensure independent processing.

## Step 3 — FASTQ → Cell Ranger matrix

Each sample is processed using:

- `cellranger count`

This generates: filtered_feature_bc_matrix/ which serves as the input for downstream Seurat analysis.


# QC Step — MD5 checksum generation

This pipeline includes a lightweight QC step that generates **MD5 checksums** for key Cell Ranger output files.

**Script:**
scripts/qc/check_filtered_matrix_md5.sh

### Purpose
Generates MD5 checksums for key Cell Ranger output files:
- barcodes.tsv.gz  
- features.tsv.gz  
- matrix.mtx.gz  

### Output
md5sums_filtered_matrix.txt

### Important note
This step:
- ✔ Generates integrity snapshots
- ❌ Does NOT compare against reference values


# Additional Scripts

## 01 — BCL/CBCL processing
**01_bcl_to_fastq_to_feature_matrix.sh**

- Converts raw 10x Genomics BCL/CBCL data into FASTQ files
- Runs Cell Ranger preprocessing
- Generates filtered feature-barcode matrices


## 02 — Seurat preprocessing & integration
**02_seurat_preprocessing_integration.R**

- Loads Cell Ranger outputs into Seurat
- Performs QC filtering
- SCTransform normalization
- Harmony-based integration (`IntegrateLayers`)
- PCA, clustering, and UMAP visualization


## 03 — QC visualization
**03_qc_visualization_single_sample.R**

- Generates QC violin plots per sample
- Evaluates:
  - nCount_RNA
  - nFeature_RNA
  - mitochondrial content


## 04 — Differential expression & annotation
**04_differential_expression_and_celltype_annotation.R**

- Differential expression analysis (cluster-level and global)
- Marker gene identification
- Cell type annotation using SingleR (DICE reference)


## 05 — Cluster composition & correlation (optional)
**05_cluster_composition_and_optional_correlation.R**

- Cluster composition analysis across samples
- Cell type proportion estimation
- Dominant cluster identification
- Optional gene-gene correlation analysis


## 06 — Functional enrichment preparation
**06_functional_enrichment_prep.R**

- Splits DE results into:
  - upregulated genes
  - downregulated genes
- Prepares input for Gene Ontology analysis (WebGestalt)


# Output Structure

results/
├── integration/
├── qc/
├── de_and_annotation/
├── de_results/
├── celltype_annotation/
├── exploratory/
├── cluster_composition/
└── gene_correlation/


# Requirements

- Cell Ranger (10x Genomics)
- Reference transcriptome (e.g., GRCh38)
- R (Seurat v5 compatible)
- Harmony
- Linux/HPC environment with sufficient compute resources


# Project-Specific Notes

## Seurat version
- Built for Seurat v5
- Uses `IntegrateLayers()` for Harmony integration

  
## FOXP3 filtering

FOXP3_counts > 0
This filtering is **specific to Treg-focused analysis** and should not be generalized.


## Cell Ranger settings
Check whether `--include-introns` is enabled depending on version and experiment design.


## Metadata consistency
Ensure consistent naming across:

- FASTQ directories
- Cell Ranger outputs
- Seurat metadata


# Summary

This pipeline provides a **modular and reproducible framework for scRNA-seq analysis**, from raw sequencing data through biological interpretation and functional enrichment.

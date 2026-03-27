## Author
Oliver Abinader

# scRNAseq-analysis
10x Genomics reproducible scRNA-seq preprocessing and Seurat workflow (spanning raw sequencing conversion through downstream analysis).

## Overview
This repository documents a single-cell RNA-seq (scRNA-seq) preprocessing and downstream analysis workflow developed for cancer research projects. 
The pipeline is designed to support reproducible analysis from raw 10x Genomics sequencing output through generation of sample-level filtered feature-barcode matrices and downstream analysis.




## Repository Structure
```bash
scripts/
  01_bcl_to_fastq_to_feature_matrix.sh
```





## Preprocessing Workflow
### Step 1: BCL/CBCL to FASTQ
If the experiment begins with raw 10x Genomics sequencing output (`BCL/CBCL`), FASTQ files are generated using `cellranger mkfastq`.

### Step 2: Organize FASTQs by Sample
After FASTQ generation, FASTQ files for each biological sample should be separated into individual sample-specific directories.

### Step 3: FASTQ to Filtered Feature-Barcode Matrix
Each sample is processed independently using `cellranger count`.
This step generates the `filtered_feature_bc_matrix/` directory for each sample, which will be used for downstream Seurat-based single-cell analysis.

### Scripts

- **01_bcl_to_fastq_to_feature_matrix.sh**  
  Converts raw 10x Genomics BCL output to FASTQ files and generates filtered feature-barcode matrices using Cell Ranger.

- **02_seurat_preprocessing_integration.R**  
  Loads Cell Ranger output into Seurat, performs QC, metadata curation, cell filtering, SCTransform normalization, Harmony integration, clustering, and UMAP visualization.

- **03_qc_visualization_single_sample.R**  
  Generates publication-style QC violin plots for a single Seurat object or subsetted sample.

- **04_differential_expression_and_celltype_annotation.R**  
  Performs differential expression analyses (global and per cluster), cluster marker detection, marker annotation, and cell type assignment using SingleR with the DICE immune reference.

- **05_cluster_composition_and_optional_correlation.R**  
  Computes cluster/sample composition, per-cluster cell type proportions, dominant cluster annotations, and includes an optional exploratory gene-gene correlation framework.

## Output Structure

Recommended output directories:

```bash
results/
├── integration/
│   ├── qc/
├── de_and_annotation/
│   ├── de_results/
│   └── celltype_annotation/
└── exploratory/
    ├── cluster_composition/
    └── gene_correlation/


Notes
The workflow is written for Seurat v5 and uses Harmony via IntegrateLayers.
Some filtering steps (e.g., FOXP3-positive cell retention) are project-specific and may need to be adjusted depending on the biological question.
The exploratory correlation module is intentionally optional and should be customized for each study.

# Important note on FOXP3
Since your current project is **Treg-focused**, the line:
FOXP3_counts > 0
is biologically meaningful for this dataset, but for GitHub/public sharing it should be documented as:
dataset-specific
optional depending on study design




## Requirements
- Cell Ranger
- 10x Genomics-compatible reference transcriptome (e.g., `refdata-gex-GRCh38-2020-A`)
- HPC or Linux environment with sufficient CPU and memory resources

## Notes
- For future projects, confirm whether `--include-introns` is enabled or defaulted in the Cell Ranger version being used.

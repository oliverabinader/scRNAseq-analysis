## Author
Oliver Abinader

# scRNAseq-analysis
10x Genomics reproducible scRNA-seq preprocessing and Seurat workflow

## Overview
This repository documents a single-cell RNA-seq (scRNA-seq) preprocessing and downstream analysis workflow developed for cancer research projects. 
The pipeline is designed to support reproducible analysis from raw 10x Genomics sequencing output through generation of sample-level filtered feature-barcode matrices and downstream analysis.


## Preprocessing Workflow
### Step 1: BCL/CBCL to FASTQ
If the experiment begins with raw 10x Genomics sequencing output (`BCL/CBCL`), FASTQ files are generated using `cellranger mkfastq`.

### Step 2: Organize FASTQs by Sample
After FASTQ generation, FASTQ files for each biological sample should be separated into individual sample-specific directories.

### Step 3: FASTQ to Filtered Feature-Barcode Matrix
Each sample is processed independently using `cellranger count`.
This step generates the `filtered_feature_bc_matrix/` directory for each sample, which will be used for downstream Seurat-based single-cell analysis.

## Repository Structure
```bash
scripts/
  01_bcl_to_fastq_to_feature_matrix.sh
```

## Requirements
- Cell Ranger
- 10x Genomics-compatible reference transcriptome (e.g., `refdata-gex-GRCh38-2020-A`)
- HPC or Linux environment with sufficient CPU and memory resources

## Notes
- For future projects, confirm whether `--include-introns` is enabled or defaulted in the Cell Ranger version being used.

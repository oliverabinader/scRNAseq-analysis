## Author
Oliver Abinader

# scRNAseq-analysis
- A reproducible 10x Genomics and Seurat-based single-cell RNA-seq (scRNA-seq) workflow spanning raw sequencing conversion through downstream analysis.


## Overview
- This repository documents a reproducible single-cell RNA-seq preprocessing and downstream analysis workflow developed for cancer-focused research projects.
- The pipeline supports analysis from raw 10x Genomics sequencing output through generation of sample-level filtered feature-barcode matrices, followed by downstream analysis in Seurat.
- The workflow is modular and designed to be adaptable for future single-cell studies, including immune-focused and tumor microenvironment applications.


## Repository Structure
scripts/
- qc/check_filtered_matrix_md5.sh
- 01_bcl_to_fastq_to_feature_matrix.sh
- 02_seurat_preprocessing_integration.R
- 03_qc_visualization_single_sample.R
- 04_differential_expression_and_celltype_annotation.R
- 05_cluster_composition_and_optional_correlation.R
- 06_functional_enrichment_prep.R


**Preprocessing Workflow**
**Step 1:** BCL/CBCL to FASTQ
- If the experiment begins with raw 10x Genomics sequencing output (BCL/CBCL), FASTQ files are generated using cellranger mkfastq.

**Step 2:** Organize FASTQs by Sample
- After FASTQ generation, FASTQ files for each biological sample should be separated into individual sample-specific directories.

**Step 3:** FASTQ to Filtered Feature-Barcode Matrix
- Each sample is processed independently using cellranger count.
- This step generates the filtered_feature_bc_matrix/ directory for each sample, which serves as the input for downstream Seurat-based single-cell analysis.
- This pipeline includes a quality control step that generates MD5 checksums for key Cell Ranger output files.


**Scripts**
- 01_bcl_to_fastq_to_feature_matrix.sh: Converts raw 10x Genomics BCL/CBCL output to FASTQ files and generates filtered feature-barcode matrices using Cell Ranger.

- 02_seurat_preprocessing_integration.R: Loads Cell Ranger output into Seurat, performs quality control (QC), cell filtering, SCTransform normalization, Harmony-based integration, dimensionality reduction, clustering, and UMAP visualization.

- 03_qc_visualization_single_sample.R: Generates publication-style QC violin plots for a single Seurat object or single sample to support its quality assessment.

- 04_differential_expression_and_celltype_annotation.R: Performs differential expression analysis (global and cluster-level), identifies marker genes, supports marker-based cluster interpretation, and assigns cell type labels using SingleR with the DICE immune reference.

- 05_cluster_composition_and_optional_correlation.R: Computes cluster/sample composition, per-cluster cell type proportions, dominant cluster annotations, and includes an optional exploratory framework for gene-gene correlation analysis.

- 06_functional_enrichment_prep.R: Filters differential expression result files into significantly upregulated and downregulated gene sets for downstream Gene Ontology enrichment analysis performed in WebGestalt.


**Output Structure**
results/
├──integration/
    ├──qc/
├──de_and_annotation/
    ├──de_results/
    ├──celltype_annotation/
├──exploratory/
    ├──cluster_composition/
    ├──gene_correlation/


**Requirements**
- Cell Ranger
- 10x Genomics-compatible reference transcriptome (e.g., refdata-gex-GRCh38-2020-A)
- R (compatible with Seurat v5 workflow)
- HPC or Linux environment with sufficient CPU and memory resources for Cell Ranger processing


**Project-Specific Notes**
- This workflow is written for Seurat v5 and uses Harmony integration via IntegrateLayers().
- Some filtering or feature-retention steps may be dataset-specific and should be reviewed before applying the workflow to other projects.
- The exploratory gene-gene correlation module in Script 05 is optional and intended as a customizable downstream extension rather than a required core scRNA-seq step.


**Important Note on FOXP3 Filtering**
- In the current project, the use of: FOXP3_counts > 0 is biologically motivated because the analysis is Treg-focused. However, this should be treated as a project-specific filtering criterion rather than a universal preprocessing rule.


**Additional Notes**
- For future projects, confirm whether --include-introns is enabled or defaulted in the Cell Ranger version being used.
- Sample naming conventions should remain consistent across FASTQ directories, Cell Ranger outputs, and Seurat metadata to simplify downstream merging and annotation.

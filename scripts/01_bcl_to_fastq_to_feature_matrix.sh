#!/bin/bash

set -euo pipefail

###############################################################################
# PART 1: BCL/CBCL -> FASTQ using cellranger mkfastq
#
# Use this section only if starting from raw sequencing output (BCL/CBCL files).
# If FASTQ files are already available, skip to PART 2.
###############################################################################

# Example command for converting BCL/CBCL files to FASTQ files:
# cellranger mkfastq \
#   --run=/path/to/sequencing_run_folder \
#   --csv=/path/to/sample_sheet.csv \
#   --output-dir=/path/to/output_fastq_dir \
#   --delete-undetermined \
#   --localcores=32 \
#   --localmem=256

# Real example used in project:
# cellranger mkfastq \
#   --run=/path/to/projects/scRNAseq_07.03.24/220209_A01277_0046_BHYNJWDRXY \
#   --csv=/path/to/projects/scRNAseq_07.03.24/220209_A01277_0046_BHYNJWDRXY/cellranger-tiny-bcl_2-9-22.csv \
#   --output-dir=/path/to/projects/scRNAseq_07.03.24/220209_A01277_0046_BHYNJWDRXY/fastq \
#   --delete-undetermined \
#   --localcores=32 \
#   --localmem=256

###############################################################################
# PART 2: Organize FASTQ files by sample
#
# IMPORTANT:
# After running mkfastq, move FASTQ files into one folder per sample.
# Example structure:
#
# /project_name/raw_data/
#   sample_1/
#     *.fastq.gz
#   sample_2/
#     *.fastq.gz
#   sample_3/
#     *.fastq.gz
#   sample_4/
#     *.fastq.gz
#
# This sample-specific organization is required before running cellranger count.
###############################################################################

###############################################################################
# PART 3: FASTQ -> filtered_feature_bc_matrix using cellranger count
#
# This section loops through sample directories and runs cellranger count independently for each sample.
#
# Usage:
#   bash 01_bcl_to_fastq_to_feature_matrix.sh <transcriptome_path> <base_fastq_dir> <output_parent_dir>
#
# Example:
#   bash 01_bcl_to_fastq_to_feature_matrix.sh \
#     /path/to/refdb_from10x/refdata-gex-GRCh38-2020-A \
#     /path/to/projects/PROJECT_NAME/raw_data \
#     /path/to/projects/PROJECT_NAME/analysis/cellranger_count
###############################################################################

# Check for the required arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <transcriptome_path> <base_fastq_dir> <output_parent_dir>"
    exit 1
fi

# User-provided inputs
TRANSCRIPTOME="$1"
BASE_FASTQ_DIR="$2"
OUTPUT_PARENT_DIR="$3"

# Create output parent directory if it does not exist
mkdir -p "${OUTPUT_PARENT_DIR}"

# Loop through each sample directory in the base FASTQ directory
for FASTQ_DIR in "${BASE_FASTQ_DIR}"/*/; do
    SAMPLE_ID=$(basename "${FASTQ_DIR}")
    SAMPLE_OUT_DIR="${OUTPUT_PARENT_DIR}/${SAMPLE_ID}"

    echo "Running cellranger count for sample: ${SAMPLE_ID}"
    echo "FASTQ directory: ${FASTQ_DIR}"
    echo "Output directory: ${SAMPLE_OUT_DIR}"

    # Run Cell Ranger count
    # Note:
    # --create-bam can be set to true if BAM output is needed for downstream use.
    # For future projects, verify whether --include-introns is already enabled by default in the Cell Ranger version in use.
    cellranger count \
        --id="${SAMPLE_ID}" \
        --fastqs="${FASTQ_DIR}" \
        --transcriptome="${TRANSCRIPTOME}" \
        --create-bam=false \
        --output-dir="${SAMPLE_OUT_DIR}" \
        > "${SAMPLE_OUT_DIR}_cellranger_count.log" 2>&1
done

echo "All cellranger count jobs completed successfully."

###############################################################################
# OPTIONAL: Example single-sample command
#
# Example:
# cellranger count \
#   --id=sample \
#   --fastqs=/path/to/projects/SC_032625-449242614/raw_data/sample/ \
#   --transcriptome=/path/to/refdb_from10x/refdata-gex-GRCh38-2020-A \
#   --create-bam=true \
#   --output-dir=/path/to/projects/SC_032625-449242614/analysis/out_1/ \
#   > ./cellranger_count.log 2>&1
###############################################################################

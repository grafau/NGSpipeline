#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p short
#SBATCH -J samtools_coverage
#SBATCH -t 0-5:00:00
#SBATCH -o logs_metrics/SamCov_%A_%a.out
#SBATCH -e logs_metrics/SamCov_%A_%a.err

# === CRITICAL: Update this number to match the TOTAL count of your BAM files ===
# Run `ls MERGED_DIR/*.RG.mapped.sorted.bam | wc -l` to get the count.
#SBATCH --array=1-50  

# Usage: bash QC_stats.sh

# Load environment
eval "$(conda shell.bash hook)"
conda activate ngs

# Directories
MERGED_DIR=~/projects/rbgk/projects/greenrice/read_processing/contemporary/leaf/mapping
DEDUP_DIR=~/projects/rbgk/projects/greenrice/read_processing/contemporary/leaf/deduplicated
QC_DIR=~/projects/rbgk/projects/greenrice/read_processing/contemporary/leaf/QC

# 3) Select BAM files matching BOTH patterns (merged and single)
# The pattern *.RG.mapped.sorted.bam captures both:
#   - HRG0790aE1.RG.mapped.sorted.bam
#   - HRG0791aE1_merged.RG.mapped.sorted.bam
mapfile -t MERGED_BAMS < <(ls "${MERGED_DIR}"/*.RG.mapped.sorted.bam)

NUM_SAMPLES=${#MERGED_BAMS[@]}

# Safety check for array ID
if [ "$SLURM_ARRAY_TASK_ID" -gt "$NUM_SAMPLES" ]; then
  echo "ERROR: Task ID $SLURM_ARRAY_TASK_ID > number of samples ($NUM_SAMPLES)"
  exit 1
fi

# Get specific BAM for this task
BAM="${MERGED_BAMS[$SLURM_ARRAY_TASK_ID-1]}"

# Extract ID dynamically
# This works for both cases:
#   Input: HRG0790aE1.RG.mapped.sorted.bam        -> ID: HRG0790aE1
#   Input: HRG0791aE1_merged.RG.mapped.sorted.bam -> ID: HRG0791aE1_merged
ID=$(basename "$BAM" .RG.mapped.sorted.bam)

echo "[$SLURM_ARRAY_TASK_ID/$NUM_SAMPLES] Processing sample: $ID"

# 4) General coverage (MAPPED reads)
samtools coverage "$BAM" \
  > "$QC_DIR/${ID}.GENcoverage.log"
echo "  - General coverage done: ${ID}.GENcoverage.log"

# --- Add Flagstat Step Here (Required for % Endogenous calculation) ---
samtools flagstat "$BAM" > "$QC_DIR/${ID}.flagstats.log"
echo "  - Flagstat done: ${ID}.flagstats.log"
# ----------------------------------------------------------------------

# 5) Usable coverage (DEDUP reads)
# Matches the input naming convention automatically
DEDUP="${DEDUP_DIR}/${ID}.RG.mapped.sorted.rmdup.bam"

if [ -f "$DEDUP" ]; then
  samtools coverage "$DEDUP" \
    > "$QC_DIR/${ID}.USEFULcoverage.log"
  echo "  - Usable coverage done: ${ID}.USEFULcoverage.log"
else
  echo "  ! Deduped BAM not found, skipping: $DEDUP"
fi

echo "Finished sample: $ID"

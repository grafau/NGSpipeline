#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p himem
#SBATCH -J a_read_processing_part2
#SBATCH -t 0-330:00:00
#SBATCH -o log/a_read_processing_part2_%A_%a.log
#SBATCH -e err/a_read_processing_part2_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE

# Module load & libraries
eval "$(conda shell.bash hook)"
conda activate ngs

# Variables
DAT=~/projects/rbgk/projects/greenrice/raw_fastq/embryo/unifIDlinks
REF=~/projects/rbgk/users_area/ykyungle/GREENrice/ref/IRGSP-1.0_genome.fasta
HOM=~/projects/rbgk/projects/greenrice/read_processing/historical/embryo
BIN=~/Tools

# Adjust SAMPLE and RUN parsing
FULL_SAMPLE=$(ls ${DAT}/*_1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p") # Extract full path
BASE_SAMPLE=$(basename $FULL_SAMPLE)                                    # Extract file name

# Extract SAMPLE and RUN
SAMPLE=$(echo $BASE_SAMPLE | sed 's/-run[234]//g' | sed 's/_1.fastq.gz//g') # Extract sample name
RUN=$(echo $BASE_SAMPLE | grep -o "\-run[234]" || echo "")                   # Extract run tag or leave empty for tagless

echo "Processing SAMPLE=${SAMPLE}, RUN=${RUN}"

# Navigate to the working directory
cd "${HOM}"

# Function to check command success
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

# Detect runs for a sample
detect_runs() {
    local sample=$1
    local runs=()

    # List all files for the given sample and extract run tags
    for file in ${DAT}/${sample}*; do
        run=$(basename "$file" | grep -o "\-run[234]")
        [[ -n "$run" ]] && runs+=("$run")
    done

    printf '%s ' "${runs[@]}"
}

# Merge BAM files for multiple run samples
merge_bam_files() {
    local sample=$1
    local merged_sample="${sample}_merged"
    local merge_files=()

    echo "Preparing to merge BAM files for SAMPLE=${sample}"

    # Detect all runs for the current sample
    for run in "" $(detect_runs "${sample}"); do
        local bam_file="${HOM}/mapping/${sample}${run}.RG.mapped.sorted.bam"
        if [ -f "${bam_file}" ]; then
            merge_files+=("${bam_file}")
            echo "Found BAM file for merging: ${bam_file}"
        fi
    done

    # Determine if merging is necessary
    if [ ${#merge_files[@]} -gt 1 ]; then
        echo "Merging ${#merge_files[@]} BAM files for SAMPLE=${sample}"
        samtools merge "${HOM}/mapping/${merged_sample}.RG.mapped.sorted.bam" "${merge_files[@]}"
        check_command "Samtools merge"
    fi    
        
    # Use merged name for downstream processing
    UPDATED_SAMPLE="${merged_sample}"
    grep -qx "$UPDATED_SAMPLE" "${HOM}/updated_samples_list.txt" 2>/dev/null \
        || echo "$UPDATED_SAMPLE" >> "${HOM}/updated_samples_list.txt"
    

    echo "BAM files merged (if applicable) for SAMPLE=${sample}. UPDATED_SAMPLE=${UPDATED_SAMPLE}"
}

# Perform merging step if necessary
    # detect runs for merging
RUNS=$(detect_runs "$SAMPLE")
RUNS_TRIM=$(echo "$RUNS" | tr -d '[:space:]')
HAS_MULTIPLE_RUNS=$([[ -n "$RUNS_TRIM" ]] && echo 1 || echo 0)

if [[ $HAS_MULTIPLE_RUNS -eq 1 && -z "$RUN" ]]; then
    # Perform merging for only one sample as representative
    merge_bam_files "${SAMPLE}"
elif [[ $HAS_MULTIPLE_RUNS -eq 0 && -z "$RUN" ]]; then
    # No run tag, proceed with deduplication and other steps
    UPDATED_SAMPLE="${SAMPLE}"
    grep -qx "$UPDATED_SAMPLE" "${HOM}/updated_samples_list.txt" 2>/dev/null \
        || echo "$UPDATED_SAMPLE" >> "${HOM}/updated_samples_list.txt"
    echo "Merging for ${UPDATED_SAMPLE} not needed."
fi

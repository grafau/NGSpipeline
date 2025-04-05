#!/bin/bash

#SBATCH --verbose
#SBATCH -J aDNA_read_processing
#SBATCH -p all
#SBATCH -o /path/to/log/a_read_processing_%A_%a.out
#SBATCH -e /path/to/log/a_read_processing_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --array=1-ARRAY_SIZE

# Activate conda environment
module purge
eval "$(conda shell.bash hook)"
conda activate conda_env

# Variables
DAT=/path/to/fastq
REF=/path/to/ref/reference.fasta
HOM=/path/to/main/working/directory
BIN=/path/to/bin

# Adjust SAMPLE and RUN parsing
FULL_SAMPLE=$(ls ${DAT}/*_1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
BASE_SAMPLE=$(basename $FULL_SAMPLE)

# Extract SAMPLE and RUN
SAMPLE=$(echo $BASE_SAMPLE | sed 's/-run[234]//g' | sed 's/_1.fastq.gz//g')
RUN=$(echo $BASE_SAMPLE | grep -o "\-run[234]" || echo "")

echo "Processing SAMPLE=${SAMPLE}, RUN=${RUN}"

# Navigate to the data directory
cd "${HOM}"

# Function to check command success
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

# Function to process until alignment and sorting
process_run() {
    local sample=$1
    local run=$2
    local rg="@RG\\tID:${sample}_${run}\\tSM:${sample}"

    # Adapter trimming
    AdapterRemoval --threads 16 --file1 ${DAT}/${sample}${run}_1.fastq.gz --file2 ${DAT}/${sample}${run}_2.fastq.gz --basename ${HOM}/trim/${sample}${run} --collapse --gzip
    check_command "AdapterRemoval"

    # Quality check - FASTQC
    fastqc -t 16 ${HOM}/trim/${sample}${run}.collapsed.gz ${HOM}/trim/${sample}${run}.pair1.truncated.gz ${HOM}/trim/${sample}${run}.pair2.truncated.gz -o ${HOM}/fastqc
    check_command "FASTQC"

    # Alignment and processing
    bwa aln -t 16 -l 1024 -f ${HOM}/mapping/${sample}${run}.collapsed.sai ${REF} ${HOM}/trim/${sample}${run}.collapsed.gz
    check_command "BWA aln ${sample}${run}"
    
    bwa samse -r $rg -f ${HOM}/mapping/${sample}${run}.RG.sam ${REF} ${HOM}/mapping/${sample}${run}.collapsed.sai ${HOM}/trim/${sample}${run}.collapsed.gz
    check_command "BWA samse ${sample}${run}"
    
    samtools flagstat ${HOM}/mapping/${sample}${run}.RG.sam > ${HOM}/mapping/${sample}${run}.flagstats.log
    check_command "Samtools flagstat ${sample}${run}"
    
    samtools view -@ 16 -F 4 -Sbh -o ${HOM}/mapping/${sample}${run}.RG.mapped.bam ${HOM}/mapping/${sample}${run}.RG.sam
    check_command "Samtools view ${sample}${run}"

    samtools sort -@ 16 -o ${HOM}/mapping/${sample}${run}.RG.mapped.sorted.bam ${HOM}/mapping/${sample}${run}.RG.mapped.bam
    check_command "Samtools sort ${sample}${run}"

    # Create a completion flag
    echo "Processing for ${sample}${run} completed."
}

# Process the current run
process_run $SAMPLE "$RUN"

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
        
        # Use merged name for downstream processing
        UPDATED_SAMPLE="${merged_sample}"
    else
        echo "Not all multiple BAM files found. Skipping merging for SAMPLE=${sample}."
        UPDATED_SAMPLE="${sample}"

        # Rename the single BAM file to follow the merged format for consistency
        mv "${merge_files[0]}" "${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam"
    fi

    echo "BAM files merged (if applicable) for SAMPLE=${sample}. UPDATED_SAMPLE=${UPDATED_SAMPLE}"
}

# Detect runs for a sample
detect_runs() {
    local sample=$1
    local runs=()

    # List all files for the given sample and extract run tags
    for file in ${DAT}/${sample}*; do
        run=$(basename "$file" | grep -o "\-run[234]")
        if [ -n "$run" ]; then
            runs+=("$run")
        fi
    done

    # Remove duplicates and return unique run tags
    echo "${runs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '
}

# Perform merging step if necessary
if [[ "$RUN" =~ ^-run[234]$ ]]; then
    # Run tag is -run[234], perform merging
    merge_bam_files "${SAMPLE}"
else
    # No run tag or other run tag, proceed with deduplication and other steps
    UPDATED_SAMPLE=$SAMPLE
fi

# Deduplication of mapped reads - DEDUP
dedup -i ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam -m -o ${HOM}/mapped/
check_command "Dedup"											

# Sorting and Indexing BAM files
samtools sort -@ 16 -o ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sorted_rmdup.bam
check_command "Sorting Samtools sort"
samtools index ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam
check_command "Indexing Samtools index"

# aDNA Authentication - AMBER
echo -e ${UPDATED_SAMPLE}'\t'${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam > ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv
${BIN}/AMBER --bamfiles ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv --output ${HOM}/aDNA_authentication/${UPDATED_SAMPLE} --errorbars --counts
check_command "AMBER"

# Haplotype calling - GATK
gatk --java-options -Xmx32G HaplotypeCaller -R ${REF} -I ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam -O ${HOM}/gvcf/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.g.vcf.gz -ERC GVCF
check_command "GATK"

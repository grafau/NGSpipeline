#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p himem
#SBATCH -J a_read_processing_part1
#SBATCH -t 0-330:00:00
#SBATCH -o log/a_read_processing_part1_%A_%a.log
#SBATCH -e err/a_read_processing_part1_%A_%a.err
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

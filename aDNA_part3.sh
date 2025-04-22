#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p all
#SBATCH -J a_read_processing_part3
#SBATCH -t 0-24:00:00
#SBATCH -o /data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA4/a_read_processing_part3_%A_%a.log
#SBATCH -e /data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA4/a_read_processing_part3_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE

# Module load & libraries
module purge
eval "$(conda shell.bash hook)"
conda activate ngs
module load samtools

# Variables
DAT=/data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA/fastq
REF=/data/users_area/yky10kg/GREENrice/Cons_Gen/datasets/ref/IRGSP-1.0_genome.fasta
HOM=/data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA4
BIN=/home/yky10kg/Tools/AMBER


# Define UPDATED_SAMPLE variable based on existing BAM files
UPDATED_SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HOM}/updated_samples_list.txt)
[[ -z "$UPDATED_SAMPLE" ]] && { echo "No sample for ID $SLURM_ARRAY_TASK_ID"; exit 0; }

echo "Using UPDATED_SAMPLE=${UPDATED_SAMPLE}"

# Navigate to the working directory
cd "${HOM}"

# Function to check command success
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

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

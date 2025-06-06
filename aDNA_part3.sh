#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p himem
#SBATCH -J a_read_processing_part3
#SBATCH -t 0-330:00:00
#SBATCH -o log/a_read_processing_part3_%A_%a.log
#SBATCH -e err/a_read_processing_part3_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE

# Module load & libraries
eval "$(conda shell.bash hook)"
conda activate ngs

# Variables
DAT=~/projects/rbgk/projects/greenrice/raw_fastq/embryo/unifIDlinks
REF=~/projects/rbgk/users_area/ykyungle/GREENrice/ref/IRGSP-1.0_genome.fasta
HOM=~/projects/rbgk/projects/greenrice/read_processing/historical/embryo
BIN=~/Tools

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
dedup -i ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam -m -o ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sorted_rmdup.bam
check_command "Dedup"											

# Sorting and Indexing BAM files
samtools sort -@ 16 -o ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sorted_rmdup.bam
check_command "Sorting Samtools sort"
samtools index ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam
check_command "Indexing Samtools index"

# aDNA Authentication - AMBER
echo -e ${UPDATED_SAMPLE}'\t'${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam > ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv
${BIN}/AMBER/AMBER --bamfiles ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv --output ${HOM}/aDNA_authentication/${UPDATED_SAMPLE} --errorbars --counts
check_command "AMBER"

# QC
## Mapped reads
echo "SAMPLE=${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
echo -e "\nMapped reads (samtools flagstat):" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
samtools flagstat ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt

## Read stats
echo -e "\nRead stats (samtools stats):" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
samtools stats ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam | grep ^SN | cut -f 2- >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt

# Haplotype calling - GATK
gatk --java-options -Xmx32G HaplotypeCaller -R ${REF} -I ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam -O ${HOM}/gvcf/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.g.vcf.gz -ERC GVCF
check_command "GATK"

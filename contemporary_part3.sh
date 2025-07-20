#!/bin/bash

#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p himem
#SBATCH -J c_read_processing_part3
#SBATCH -t 0-330:00:00
#SBATCH -o log/c_read_processing_part3_%A_%a.log
#SBATCH -e err/c_read_processing_part3_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE

# Module load & libraries
eval "$(conda shell.bash hook)"
conda activate ngs

# Variables
DAT=~/projects/rbgk/users_area/ykyungle/GREENrice/contemporary_test/fastq
REF=~/projects/rbgk/users_area/ykyungle/GREENrice/ref/IRGSP-1.0_genome.fasta
HOM=~/projects/rbgk/users_area/ykyungle/GREENrice/contemporary_test
BIN=/mnt/apps/users/ykyungle/conda/pkgs/picard-3.4.0-hdfd78af_0/share/picard-3.4.0-0

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

# Deduplication of mapped reads - Picard
java -Xmx16G -jar ${BIN}/picard.jar MarkDuplicates I=${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam O=${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.bam M=${HOM}/deduplicated/${UPDATED_SAMPLE}_MarkDup.metrics.txt REMOVE_DUPLICATES=true
check_command "Picard-rmdup ${UPDATED_SAMPLE}"

# Sorting and Indexing BAM files
samtools sort -@ 16 -o ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.bam
check_command "Samtools sort ${UPDATED_SAMPLE}"
samtools index ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam
check_command "Samtools index ${UPDATED_SAMPLE}"

# QC
## Mapped reads
echo "SAMPLE=${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
echo -e "\nMapped reads (samtools flagstat):" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
samtools flagstat ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt

## Read stats
echo -e "\nRead stats (samtools stats):" >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt
samtools stats ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam | grep ^SN | cut -f 2- >> ${HOM}/QC/quality_check_${UPDATED_SAMPLE}.txt

# Haplotype calling - GATK
gatk --java-options -Xmx32G HaplotypeCaller -R ${REF} -I ${HOM}/deduplicated/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.bam -O ${HOM}/gvcf/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.sorted.g.vcf.gz -ERC GVCF
check_command "GATK"

#!/bin/bash
#SBATCH --verbose
#SBATCH --mem=20G
#SBATCH -c 12
#SBATCH -p all
#SBATCH -J c_read_processing
#SBATCH -t 0-24:00:00
#SBATCH -o /path/to/log/c_read_processing_%A_%a.log
#SBATCH -e /path/to/log/c_read_processing_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE

#module load & libraries
module purge
eval "$(conda shell.bash hook)"
conda activate ngs

#Submit command: sbatch contemporary_read_processing.sh (or use Slurm_submit.sh for automated ARRAY_SIZE calculation)

#Print the task ID
cd "$SLURM_SUBMIT_DIR"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID ${SAMPLE}

#Variables
DAT=/data/users_area/yky10kg/GREENrice/Cons_Gen/datasets/farmers/fastq
REF=/data/users_area/yky10kg/GREENrice/Cons_Gen/datasets/ref
HOM=/data/users_area/yky10kg/GREENrice/Cons_Gen/datasets/farmers/trial
SAMPLE=$(ls ${DAT}/*_1.fastq.gz | rev | cut -d "/" -f 1 | rev | cut -f 1 -d "_" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Function to check command success
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

# Code
# Initial quality check - FastQC
fastqc ${DAT}/${SAMPLE}_1.fastq.gz ${DAT}/${SAMPLE}_2.fastq.gz -o ${HOM}/fastqc
check_command "FastQC"

# Adapter/quality trimming - Trimmomatic
trimmomatic PE -phred33 ${DAT}/${SAMPLE}_1.fastq.gz ${DAT}/${SAMPLE}_2.fastq.gz ${HOM}/trim/${SAMPLE}_1_paired.fastq.gz ${HOM}/trim/${SAMPLE}_1_unpaired.fastq.gz ${HOM}/trim/${SAMPLE}_2_paired.fastq.gz ${HOM}/trim/${SAMPLE}_2_unpaired.fastq.gz ILLUMINACLIP:/home/yky10kg/anaconda3/env/ngs/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
check_command "Trimmomatic"

# Mapping - bwa mem
bwa mem ${REF}/IRGSP-1.0_genome.fasta ${HOM}/trim/${SAMPLE}_1_paired.fastq.gz ${HOM}/trim/${SAMPLE}_2_paired.fastq.gz > ${HOM}/mapping/${SAMPLE}.sam
check_command "BWA mem"

bwa samse -r @RG\tID:${SAMPLE}\tSM:${SAMPLE} -f ${HOM}/mapping/${SAMPLE}.RG.sam ${REF}/IRGSP-1.0_genome.fasta ${HOM}/mapping/${SAMPLE}.collapsed.sai ${HOM}/trim/${SAMPLE}.collapsed.gz
check_command "BWA samse"

# Deduplication of mapped reads - samtools, picard
## Mapped
${BIN}/samtools view -bS -F 4 ${HOM}/mapping/${SAMPLE}.RG.sam > ${HOM}/mapped/${SAMPLE}.RG.mapped.bam
check_command "Samtools view"

${BIN}/samtools sort ${HOM}/mapped/${SAMPLE}.RG.mapped.bam -o ${HOM}/mapped/${SAMPLE}.RG.mapped.sort.bam
check_command "Samtools sort"

java -Xmx8G -jar ~/anaconda3/pkgs/picard-3.2.0-hdfd78af_0/share/picard-3.2.0-0/picard.jar MarkDuplicates I=${HOM}/mapped/${SAMPLE}.RG.mapped.sort.bam O=${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam M=${HOM}/mapped/${SAMPLE}.metrics.txt REMOVE_DUPLICATES=true
check_command "Picard-rmdup"

## Unmapped
${BIN}/samtools view -bS -f 4 ${HOM}/mapping/${SAMPLE}.RG.sam > ${HOM}/unmapped/${SAMPLE}.RG.unmapped.bam
check_command "Samtools view"

${BIN}/samtools sort ${HOM}/unmapped/${SAMPLE}.RG.unmapped.bam -o ${HOM}/unmapped/${SAMPLE}.RG.unmapped.sort.bam
check_command "Samtools sort"

java -Xmx3G -jar ~/anaconda3/pkgs/picard-3.2.0-hdfd78af_0/share/picard-3.2.0-0/picard.jar MarkDuplicates I=${HOM}/unmapped/${SAMPLE}.RG.unmapped.sort.bam O=${HOM}/unmapped/${SAMPLE}.RG.unmapped.sort.rmdup.bam M=${HOM}/unmapped/${SAMPLE}.metrics.txt REMOVE_DUPLICATES=true
check_command "Picard-rmdup"

# QC
## Raw fastq reads
echo "Raw fastq_1 reads:" > ${HOM}/QC/quality_check_${SAMPLE}.txt
zcat ${DAT}/${SAMPLE}_1.fastq.gz | wc -l | awk '{print $1/4}' >> ${HOM}/QC/quality_check_${SAMPLE}.txt
echo "Raw fastq_2 reads:" >> ${HOM}/QC/quality_check_${SAMPLE}.txt
zcat ${DAT}/${SAMPLE}_2.fastq.gz | wc -l | awk '{print $1/4}' >> ${HOM}/QC/quality_check_${SAMPLE}.txt

## Mapped reads
echo -e "\nMapped reads:" >> ${HOM}/QC/quality_check_${SAMPLE}.txt
${BIN}/samtools flagstat ${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam >> ${HOM}/QC/quality_check_${SAMPLE}.txt

## Unmapped reads
echo -e "\nUnmapped reads:" >> ${HOM}/QC/quality_check_${SAMPLE}.txt
${BIN}/samtools flagstat ${HOM}/unmapped/${SAMPLE}.RG.unmapped.sort.rmdup.bam >> ${HOM}/QC/quality_check_${SAMPLE}.txt

## Read stats
echo -e "\nRead stats:" >> ${HOM}/QC/quality_check_${SAMPLE}.txt
echo "SAMPLE=${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam" >> ${HOM}/QC/quality_check_${SAMPLE}.txt
${BIN}/samtools stats ${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam | grep ^SN | cut -f 2- >> ${HOM}/QC/quality_check_${SAMPLE}.txt

# Variant calling
## Index BAM files - samtools
${BIN}/samtools index ${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam
check_command "Samtools index"

## Haplotype calling - GATK
gatk --java-options -Xmx32G HaplotypeCaller -R ${REF}/IRGSP-1.0_genome.fasta -I ${HOM}/mapped/${SAMPLE}.RG.mapped.sort.rmdup.bam -O ${HOM}/gvcf/${SAMPLE}.RG.mapped.sort.rmdup.g.vcf.gz -ERC GVCF
check_command "GATK"

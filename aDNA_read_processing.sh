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
source /home/pho10kg/miniconda/etc/profile.d/conda.sh
conda activate /home/pho10kg/miniconda/envs/aDNA_Nov_env/

# Get ID for current array job
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Running analysis on ${SAMPLE}"

# Variables
DAT=/path/to/fastq
REF=/path/to/ref/reference.fasta
HOM=/path/to/main/working/directory
BIN=/path/to/bin
SAMPLE=$(ls ${DAT}/*_1.fastq.gz | rev | cut -d "/" -f 1 | rev | cut -f 1 -d "_" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Navigate to the data directory
cd "${HOM}"

# Function to check command success
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

# Function to process MAPPING for a single run
process_single_run() {
    local sample=$1
    local run=$2
    local rg="@RG\\tID:${sample}_${run}\\tSM:${sample}"

    bwa aln -t 16 -l 1024 -f ${HOM}/mapping/${sample}_${run}.collapsed.sai ${REF} ${HOM}/trim/${sample}.collapsed.gz
    check_command "BWA aln ${sample}_${run}"
    
    bwa samse -r $rg -f ${HOM}/mapping/${sample}_${run}.RG.sam ${REF} ${HOM}/mapping/${sample}_${run}.collapsed.sai ${HOM}/trim/${sample}.collapsed.gz
    check_command "BWA samse ${sample}_${run}"
    
    samtools flagstat ${HOM}/mapping/${sample}_${run}.RG.sam > ${HOM}/mapping/${sample}_${run}.flagstats.log
    check_command "Samtools flagstat ${sample}_${run}"
    
    samtools view -@ 16 -F 4 -Sbh -o ${HOM}/mapping/${sample}_${run}.RG.mapped.bam ${HOM}/mapping/${sample}_${run}.RG.sam
    check_command "Samtools view ${sample}_${run}"
    
    samtools sort -@ 16 -o ${HOM}/mapping/${sample}_${run}.RG.mapped.sorted.bam ${HOM}/mapping/${sample}_${run}.RG.mapped.bam
    check_command "Samtools sort ${sample}_${run}"
}


# Run pipeline
## Adapter trimming - ADAPTERREMOVAL
AdapterRemoval --file1 ${DAT}/${SAMPLE}_1.fastq.gz --file2 ${DAT}/${SAMPLE}_2.fastq.gz --basename ${HOM}/trim/${SAMPLE} --collapse --gzip --threads 16
check_command "AdapterRemoval" 

## Quality check - FASTQC
fastqc -t 16 ${HOM}/trim/${SAMPLE}.collapsed.gz ${HOM}/trim/${SAMPLE}.pair1.truncated.gz ${HOM}/trim/${SAMPLE}.pair2.truncated.gz -o ${HOM}/fastqc
check_command "FastQC"

## MAPPING
if [ -f "${DAT}/${SAMPLE}_run2_1.fastq.gz" ]; then
    # Process both batches
    	process_single_run $SAMPLE "run1"
    	process_single_run $SAMPLE "run2"

    # Merge BAM files
    	UPDATED_SAMPLE="${SAMPLE}_merged"
    	samtools merge ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam ${HOM}/mapping/${SAMPLE}_run1.RG.mapped.sorted.bam ${HOM}/mapping/${SAMPLE}_run2.RG.mapped.sorted.bam
    	check_command "Samtools merge"
else
    # Process single run
    	UPDATED_SAMPLE="${SAMPLE}"
    	process_single_batch $SAMPLE "run1"
    	mv ${HOM}/mapping/${SAMPLE}_run1.RG.mapped.sorted.bam ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam
fi

## Deduplication of mapped reads - DEDUP
dedup -i ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.bam -m -o ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.bam
check_command "Dedup"											

## aDNA Authentication - AMBER
samtools sort -@ 16 -o ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam ${HOM}/mapping/${UPDATED_SAMPLE}.RG.mapped.sorted.rmdup.bam
echo -e ${UPDATED_SAMPLE}'\t'${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam > ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv
${BIN}/AMBER --bamfiles ${HOM}/aDNA_authentication/${UPDATED_SAMPLE}.tsv --output ${HOM}/aDNA_authentication/${UPDATED_SAMPLE} --errorbars --counts

# Variant calling
## Index BAM files - samtools
samtools index ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam
check_command "Samtools index"

## Haplotype calling - GATK
gatk --java-options -Xmx32G HaplotypeCaller -R ${REF}/IRGSP-1.0_genome.fasta -I ${HOM}/mapped/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.bam -O ${HOM}/gvcf/${UPDATED_SAMPLE}.RG.mapped.sort.rmdup.g.vcf.gz -ERC GVCF
check_command "GATK"

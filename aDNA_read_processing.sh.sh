#!/bin/bash

#SBATCH --verbose
#SBATCH -J aDNA_read_processing
#SBATCH -p all
#SBATCH -o /data/users_area/pho10kg/output/arrayjobs_%A_%a.out
#SBATCH -e /data/users_area/pho10kg/errors/arrayjobs_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --array=1-ARRAY_SIZE


# Variables 
HOM="/data/users_area/pho10kg/seq_data"
REF="${HOM}/ref/Nipponbare_reference_genome.fasta"
SAMPLES=($(ls ${HOM}/raw_data/*_1.fastq.gz | sed 's|.*/||; s/_1.fastq.gz//'))
BIN="/home/pho10kg/bin"

#Number of samples
NUM_SAMPLES=${#SAMPLES[@]}

# Get ID for current array job

   SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Running analysis on ${SAMPLE}"

# Activate conda environment
source /home/pho10kg/miniconda/etc/profile.d/conda.sh
conda activate /home/pho10kg/miniconda/envs/aDNA_Nov_env/
 
# Navigate to the data directory
  cd "${HOM}"

# Run pipeline

 ## ADAPTERREMOVAL
	AdapterRemoval --file1 ${HOM}/raw_data/${SAMPLE}_1.fastq.gz --file2 ${HOM}/raw_data/${SAMPLE}_2.fastq.gz --basename ${HOM}/trimmed_merged/${SAMPLE} --collapse --gzip --threads 16

 ## FASTQC
	fastqc -t 16 -o ${HOM}/quality_control/ ${HOM}/trimmed_merged/${SAMPLE}.collapsed.gz ${HOM}/trimmed_merged/${SAMPLE}.pair1.truncated.gz ${HOM}/trimmed_merged/${SAMPLE}.pair2.truncated.gz	

 ## MAPPING
	bwa aln -t 16 -l 1024 -f ${HOM}/mapping/${SAMPLE}.collapsed.sai ${REF} ${HOM}/trimmed_merged/${SAMPLE}.collapsed.gz
									
	bwa samse -r @RG\\tID:${SAMPLE}\\tSM:${SAMPLE} -f ${HOM}/mapping/${SAMPLE}.sam ${REF} ${HOM}/mapping/${SAMPLE}.collapsed.sai ${HOM}/trimmed_merged/${SAMPLE}.collapsed.gz

	samtools flagstat ${HOM}/mapping/${SAMPLE}.sam > ${HOM}/mapping/${SAMPLE}.flagstats.log

	samtools view -@ 16 -F 4 -Sbh -o ${HOM}/mapping/${SAMPLE}.mapped.bam ${HOM}/mapping/${SAMPLE}.sam

	samtools sort -@ 16 -o ${HOM}/mapping/${SAMPLE}.mapped.sorted.bam ${HOM}/mapping/${SAMPLE}.mapped.bam

 ## DEDUP
	dedup -i ${HOM}/mapping/${SAMPLE}.mapped.sorted.bam -m -o ${HOM}/mapping/											

 ## AMBER

	samtools sort -@ 16 -o ${HOM}/mapping/${SAMPLE}.rmdup.sorted.bam ${HOM}/mapping/${SAMPLE}.mapped.sorted_rmdup.bam

	echo -e ${SAMPLE}'\t'${HOM}/mapping/${SAMPLE}.rmdup.sorted.bam > ${HOM}/aDNA_characteristics/paths/${SAMPLE}.tsv
	
	${BIN}/AMBER --bamfiles ${HOM}/aDNA_characteristics/paths/${SAMPLE}.tsv --output ${HOM}/aDNA_characteristics/${SAMPLE} --errorbars --counts




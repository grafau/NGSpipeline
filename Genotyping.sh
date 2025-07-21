#!/bin/bash
#SBATCH --verbose
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -p himem
#SBATCH -J Genotyping
#SBATCH -t 0-330:00:00
#SBATCH -o log/Genotyping_%j.log
#SBATCH -e log/Genotyping_%j.err

#module load & libraries
eval "$(conda shell.bash hook)"
conda activate ngs

#Submit command: sbatch ./Genotyping.sh

#Print the task ID
cd "$SLURM_SUBMIT_DIR"
echo "My SLURM_JOB_ID: " $SLURM_JOB_ID

#Variables
DAT=/path/to/g.vcf
REF=/path/to/ref
HOM=/path/to/working/directory


#Code
#Sample list for genotyping
ls ${DAT}/*.RG.mapped.sorted.rmdup.sorted.g.vcf.gz > ${HOM}/gvcf.list

#Combine GVCFs
gatk --java-options "-Xmx32G" CombineGVCFs -R ${REF}/IRGSP-1.0_genome.fasta --variant ${HOM}/gvcf.list -O ${HOM}/combined.g.vcf.gz

# Validate Variants
gatk --java-options "-Xmx32G" ValidateVariants -V ${HOM}/combined.g.vcf.gz

# Genotype GVCFs
gatk --java-options "-Xmx32G" GenotypeGVCFs -R ${REF}/IRGSP-1.0_genome.fasta -V ${HOM}/combined.g.vcf.gz -O ${HOM}/combined.vcf.gz

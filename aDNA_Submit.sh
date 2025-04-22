#!/bin/bash

# Common path and output directory
BASE=/data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA4
FASTQ_DIR=/data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA/fastq
mkdir -p $BASE/{fastqc,trim,mapping,mapped,unmapped,QC,aDNA_authentication,gvcf,vcf}

# Calculate the array size (For PART1&2)
ARRAY_SIZE=$(ls $FASTQ_DIR/*_1.fastq.gz | wc -l)

# Replace ARRAY_SIZE in the original scripts and generate modified scripts PART1&2
sed "s/ARRAY_SIZE/${ARRAY_SIZE}/" aDNA_part1.sh > $BASE/Array_aDNA_part1.sh
sed "s/ARRAY_SIZE/${ARRAY_SIZE}/" aDNA_part2.sh > $BASE/Array_aDNA_part2.sh

# Submit the job for part 1
PART1_JOB_ID=$(sbatch $BASE/Array_aDNA_part1.sh | awk '{print $4}')

# Submit the job for part 2 with dependency on part 1
PART2_JOB_ID=$(sbatch --dependency=afterok:$PART1_JOB_ID $BASE/Array_aDNA_part2.sh | awk '{print $4}')

# Create and submit PART3 wrapper job
cat > Submit_aDNA_part3_wrapper.sh <<EOF
#!/bin/bash
#SBATCH -J Submit_aDNA_part3
#SBATCH --dependency=afterok:${PART2_JOB_ID}
#SBATCH -p all
#SBATCH --mem=1G
#SBATCH -t 00:05:00

BASE=/data/users_area/yky10kg/GREENrice/Cons_Gen/pipeline_testing/aDNA4

# UPDATED_SAMPLE(after merging) into new array size
LINES=\$(wc -l < $BASE/updated_samples_list.txt)

# Update script
sed "s/ARRAY_SIZE/${LINES}/" $BASE/aDNA_part3.sh > $BASE/Array_aDNA_part3.sh
chmod +x $BASE/Array_aDNA_part3.sh

# Submit the job for Part3
PART3_JOB_ID=$(sbatch --dependency=afterok:$PART3_WRAPPER_JOB_ID $BASE/Array_aDNA_part3.sh | awk '{print $4}')
EOF

# Submit the wrapper job
chmod +x $BASE/Submit_aDNA_part3_wrapper.sh
PART3_WRAPPER_JOB_ID=$(sbatch --dependency=afterok:$PART2_JOB_ID $BASE/Submit_aDNA_part3_wrapper.sh | awk '{print $4}')


# Job Log
{
  echo "DATE                  : $(date)"
  echo "ARRAY_SIZE(FATQ)      : $ARRAY_SIZE"
  echo "SLURM_JOB_ID          : $SLURM_JOB_ID"
  echo "part1 JobID           : $PART1_JOB_ID"
  echo "part2 JobID           : $PART2_JOB_ID"
  echo "part3 wrapper         : $PART3_WRAPPER_JOB_ID"
  echo "UPDATED_ARRAY(merged) : $ARRAY_SIZE"
  echo "part3 JobID           : $PART3_JOB_ID"
  echo "-------------------------------------------------"
} >> $BASE/Slurm_Record.log

#!/bin/bash

# Common path and output directory
BASE=~/projects/rbgk/projects/greenrice/read_processing/historical/embryo
FASTQ_DIR=~/projects/rbgk/projects/greenrice/raw_fastq/embryo/unifIDlinks
mkdir -p $BASE/{log,err,fastqc,trim,mapping,deduplicated,QC,aDNA_authentication,gvcf}

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
WRAPPER=$BASE/Submit_aDNA_part3_wrapper.sh

cat > "$WRAPPER" <<'EOF'
#!/bin/bash
#SBATCH -J Submit_aDNA_part3
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH -t 00:05:00
#SBATCH --dependency=afterok:__PART2__

BASE=__BASE__

# UPDATED_SAMPLE(after merging) into new array size
LINES=$(wc -l < "$BASE/updated_samples_list.txt")
echo "UPDATED_ARRAY(merged) : $LINES" >> "$BASE/Slurm_Record.log"

# Update script
sed "s/ARRAY_SIZE/${LINES}/" "$BASE/aDNA_part3.sh" > "$BASE/Array_aDNA_part3.sh"
chmod +x "$BASE/Array_aDNA_part3.sh"

# Submit the job for Part3
PART3_JOB_ID=$(sbatch "$BASE/Array_aDNA_part3.sh" | awk '{print $4}')
echo "part3 JobID           : $PART3_JOB_ID" >> "$BASE/Slurm_Record.log" 
echo "-------------------------------------------------" >> "$BASE/Slurm_Record.log" 
EOF

# Submit the wrapper job
sed -e "s#__PART2__#${PART2_JOB_ID}#" \
    -e "s#__BASE__#${BASE}#"          \
    "$WRAPPER" > "${WRAPPER}.tmp"

mv "${WRAPPER}.tmp" "$WRAPPER"
chmod +x "$WRAPPER"

PART3_WRAPPER_JOB_ID=$(sbatch "$WRAPPER" | awk '{print $4}')
echo "part3 wrapper JobID   : $PART3_WRAPPER_JOB_ID"

# Job Log
{
  echo "DATE                  : $(date)"
  echo "ARRAY_SIZE(FATQ)      : $ARRAY_SIZE"
  echo "SLURM_JOB_ID          : $SLURM_JOB_ID"
  echo "part1 JobID           : $PART1_JOB_ID"
  echo "part2 JobID           : $PART2_JOB_ID"
  echo "part3 wrapper JobID   : $PART3_WRAPPER_JOB_ID"

} >> $BASE/Slurm_Record.log

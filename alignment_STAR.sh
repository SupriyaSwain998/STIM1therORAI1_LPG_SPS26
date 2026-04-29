#!/bin/bash
#SBATCH -t 4-00:00:00                   # Max runtime (4 days)
#SBATCH -c 10                            # Number of cores to use
#SBATCH --mem=40G                        # Memory allocation
#SBATCH --job-name=STAR_alignment_PE     # Job name
#SBATCH --array=0-23                     # Array job for 24 samples
#SBATCH -o alignment_output_%A_%a.log    # Standard output log
#SBATCH -e alignment_error_%A_%a.log     # Standard error log

# Load STAR
module load star

# Define paths
GENOME_INDEX_DIR="genome_index"
OUTPUT_DIR="03_Alignment"
INPUT_DIR="01_Trimming"

# Create output directory if not exists
mkdir -p ${OUTPUT_DIR}

# Define sample array (23 samples) (without KZLX9)
sample_array=(
"KZLX1" "KZLX2" "KZLX3" "KZLX4" "KZLX5" "KZLX6"
"KZLX7" "KZLX8" "KZLX10" "KZLX11" "KZLX12"
"KZLX13" "KZLX14" "KZLX15" "KZLX16" "KZLX17" "KZLX18"
"KZLX19" "KZLX20" "KZLX21" "KZLX22" "KZLX23" "KZLX24"
)

# Get current sample name
SAMPLE=${sample_array[$SLURM_ARRAY_TASK_ID]}

# Run STAR for paired-end reads
srun STAR --genomeDir $GENOME_INDEX_DIR \
    --runThreadN 10 \
    --readFilesIn ${INPUT_DIR}/${SAMPLE}_trimmed.R1.fastq.gz ${INPUT_DIR}/${SAMPLE}_trimmed.R2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}_ \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outSAMattributes Standard \
    --quantMode TranscriptomeSAM GeneCounts


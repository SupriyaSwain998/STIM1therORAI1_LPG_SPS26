# STIM1therORAI1_LPG_SPS26

RNA-Seq Analysis Pipeline
Overview

This repository contains a complete RNA-Seq analysis workflow for the project, including:

Read alignment using STAR
Generation of gene-level counts
Construction of expression matrix
Downstream analysis in R (normalization, differential expression, visualization, rescue status analysis, GO enrichment)
Workflow Summary
1. Read Alignment (STAR)

Alignment was performed using a SLURM-based batch script with STAR for paired-end reads.

Key features:

Parallel processing using SLURM array jobs
23 samples processed
Output includes:
Sorted BAM files
Gene counts (ReadsPerGene.out.tab)
Transcriptome alignments

Script:
alignment_STAR.sh

Important parameters:

Threads: 10
Memory: 40 GB
Input: trimmed FASTQ files
Output: BAM + gene counts
2. Gene Count Extraction

STAR generates gene-level counts directly using:

--quantMode GeneCounts

These count files were collected for all samples and used to build the expression matrix.

3. Expression Matrix Construction
Individual gene count files were merged
A matrix was created with:
Rows: genes
Columns: samples

This matrix serves as input for downstream analysis.

4. Downstream Analysis (R)

All downstream analyses were performed using the provided R script:

Script:
260428_analysis.R

Steps included:
Data import and preprocessing
Normalization of counts
Differential expression analysis
Filtering of significant genes
Rescue status analysis #Comparison of the treatments on the disease model
GO enrichment analysis (Clusterprofiler, Reactome)
Visualization (e.g., complex heatmaps, plots)

Directory Structure           
├── alignment_STAR.sh         # Alignment script
├── 20251023_merge.R and 20251023_merge.sh      #merging of alignment files to create an unified matrix
├── mergedReadCounts.csv       #readcount matrix
├── metadata.xlsx             #metadata of the samples
├── 260428_analysis.R         # Downstream analysis
├── 260428_output             # Output files obtained during the downstream analysis
└── README.md


Requirements
Software
STAR
R (with required packages for RNA-seq analysis, e.g. DESeq2 / edgeR / ggplot2)
HPC Environment
SLURM workload manager
Sufficient memory (≥40 GB recommended)
How to Run
Step 1: Alignment
Submit the STAR alignment job:
sbatch alignment_STAR.sh
Step 2: Generate Matrix
Collect ReadsPerGene.out.tab files
Merge into a count matrix (custom script or R)
Step 3: Run Analysis
source("260428_analysis.R") #Update the file paths inside the R script if needed
Notes
Ensure genome index is correctly built before running STAR
Sample names must match between FASTQ files and scripts
Adjust SLURM parameters based on your HPC system



Author
Supriya Priyadarshani Swain

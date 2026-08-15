# STIM1 / ORAI1 — LPG SPS26 RNA-Seq Analysis Pipeline

> **Project:** STIM1therORAI1_LPG_SPS26  
> **Author:** Supriya Priyadarshani Swain  
> **Last updated:** April 2026

---

## Overview

This repository contains a complete, end-to-end RNA-Seq analysis workflow, from raw FASTQ alignment through differential expression, rescue status classification, and pathway enrichment.

**Pipeline stages at a glance:**

| Stage | Tool / Environment | Output |
|---|---|---|
| Read alignment | STAR (SLURM HPC) | Sorted BAM + gene counts |
| Count extraction | STAR `--quantMode GeneCounts` | `ReadsPerGene.out.tab` per sample |
| Matrix construction | R / custom merge script | `mergedReadCounts.csv` |
| Downstream analysis | R (DESeq2, ClusterProfiler, …) | Normalized counts, DE results, plots |

---

## Repository Structure

```
STIM1therORAI1_LPG_SPS26/
├── alignment_STAR.sh            # SLURM batch script for STAR alignment
├── 20251023_merge.R             # R script — merges per-sample count files
├── 20251023_merge.sh            # Shell wrapper for the merge step
├── mergedReadCounts.csv         # Unified read-count matrix (genes × samples)
├── metadata.xlsx                # Sample metadata
├── 260428_analysis.R            # Full downstream analysis script
├── 260428_output/               # All outputs from downstream analysis
└── README.md
```

---

## Workflow Details

### 1 · Read Alignment (STAR)

Paired-end reads were aligned to the reference genome using **STAR**, submitted as a SLURM array job.

**Script:** `alignment_STAR.sh`

| Parameter | Value |
|---|---|
| Samples processed | 23 |
| Threads per job | 10 |
| Memory per job | 40 GB |
| Input | Trimmed FASTQ files |
| Output | Sorted BAM, `ReadsPerGene.out.tab`, transcriptome BAM |

Key STAR flags used:
- `--runMode alignReads`
- `--quantMode GeneCounts` — generates strand-aware gene-level counts directly
- `--outSAMtype BAM SortedByCoordinate`

---

### 2 · Gene Count Extraction

STAR's `--quantMode GeneCounts` produces one `ReadsPerGene.out.tab` file per sample containing unstranded, forward-stranded, and reverse-stranded counts. The appropriate column was selected based on library strandedness.

---

### 3 · Expression Matrix Construction

**Scripts:** `20251023_merge.R` + `20251023_merge.sh`

Individual count files were merged into a single matrix:
- **Rows:** Ensembl gene IDs
- **Columns:** Sample identifiers (matching metadata)
- **Output:** `mergedReadCounts.csv`

This matrix is the primary input for all downstream analyses.

---

### 4 · Downstream Analysis (R)

**Script:** `260428_analysis.R`

All statistical and biological analyses are performed in this single R script. Update file paths at the top of the script before running.

**Steps included:**

1. **Data import & preprocessing** — load count matrix and metadata; filter low-count genes
2. **Normalization** — library-size normalization (DESeq2 VST / rlog or edgeR TMM)
3. **Differential expression** — pairwise comparisons; adjusted *p*-value and log₂FC thresholds applied
4. **Significant gene filtering** — extraction of DE gene lists per contrast
5. **Rescue status analysis** — classification of genes by treatment effect on the disease model (treatment vs. disease vs. control comparisons)
6. **GO / pathway enrichment** — over-representation analysis via `clusterProfiler`; Reactome pathway analysis
7. **Visualization** — complex heatmaps (`ComplexHeatmap`), volcano plots, PCA, dot plots

Outputs are written to `260428_output/`.

---

## Requirements

### Software

| Tool | Purpose |
|---|---|
| [STAR](https://github.com/alexdobin/STAR) ≥ 2.7 | Read alignment |
| R ≥ 4.2 | Downstream analysis |
| DESeq2 | Normalization & differential expression |
| edgeR | Alternative normalization |
| clusterProfiler | GO & pathway enrichment |
| ReactomePA | Reactome enrichment |
| ComplexHeatmap | Heatmap visualization |
| ggplot2 | General plotting |

### HPC Environment

- **Scheduler:** SLURM workload manager
- **Memory:** ≥ 40 GB per alignment job
- **Reference genome index:** Must be pre-built with `STAR --runMode genomeGenerate` before running `alignment_STAR.sh`

---

## How to Run

### Step 1 — Align reads

```bash
sbatch alignment_STAR.sh
```

Submits a SLURM array job (one task per sample). Monitor with `squeue -u $USER`.

### Step 2 — Build the count matrix

```bash
bash 20251023_merge.sh
# or directly in R:
source("20251023_merge.R")
```

Collects all `ReadsPerGene.out.tab` files and writes `mergedReadCounts.csv`.

### Step 3 — Run downstream analysis

```R
source("260428_analysis.R")
```

> **Before running:** update the file path variables at the top of the script to match your local directory structure.

---

## Notes & Common Pitfalls

- **Genome index** — must be built for the same STAR version used for alignment; mismatches cause silent errors.
- **Sample name consistency** — sample identifiers must match exactly between FASTQ filenames, `metadata.xlsx`, and the count matrix columns.
- **Strandedness** — verify library strandedness (e.g., with `RSeQC infer_experiment.py`) and select the correct column from `ReadsPerGene.out.tab` before merging.
- **SLURM resources** — adjust `--mem`, `--cpus-per-task`, and `--array` range in `alignment_STAR.sh` to match your HPC cluster's policies.
- **R package versions** — Bioconductor packages are version-sensitive; use `BiocManager::install()` and record package versions in your session info (`sessionInfo()`).

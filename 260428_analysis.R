#hierarchial clustering
# =========================
# 1. Load libraries
# =========================
library(readxl)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)

# =========================
# 2. Read count data and metadata
# =========================
countData <- read.table("mergedReadCounts.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
coldata <- read_excel("metadata.xlsx")
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
remove_sample = "KZLX9"
countData <- countData[, colnames(countData) != remove_sample]
coldata   <- coldata[rownames(coldata) != remove_sample, ]


# Ensure alignment
coldata <- coldata[colnames(countData), , drop = FALSE]
all(colnames(countData) %in% coldata$sample)  # TRUE

# Keep conditions as unique factor
coldata$condition <- factor(trimws(coldata$condition))
levels(coldata$condition)  # e.g., "KIe", "KIn", "MOE", "sh", "WTe", "WTn"

# =========================
# 3. DESeq2 and VST
# =========================
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = coldata,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # filter low-count genes
vsd <- vst(dds)
stabilized_counts <- assay(vsd)

# =========================
# 4. Define consistent colors for conditions
# =========================
condition_levels <- levels(coldata$condition)
condition_colors <- c(
  "WT_empty" = "black",
  "WT_Nacl" = "darkgrey",
  "Stim1R304W/+_sh190" = "#FF7F00",
  "Stim1R304W/+_MOE" = "yellow",
  "Stim1R304W/+_empty" = "red",
  "Stim1R304W/+_Nacl" = "brown"
)


# =========================
# 5. Hierarchical clustering with Sample IDs at bottom
# =========================
dist_counts <- dist(t(stabilized_counts))

# Annotation for condition only (no color for SampleID)
annotation_col <- data.frame(condition = coldata$condition)
rownames(annotation_col) <- rownames(coldata)
png(
  filename = "260428_output/20260130_hierarchial_clustering.png",
  width = 10,
  height = 10,
  units = "in",
  res = 300
)
# Plot heatmap
pheatmap(
  dist_counts,
  annotation_col = annotation_col,
  annotation_colors = list(condition = condition_colors),  # Only condition colored
  annotation_names_col = FALSE,
  show_colnames = TRUE,
  labels_col = rownames(coldata),     # This will put sample IDs under columns
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  clustering_method = "complete")
dev.off()
# Save hierarchical clustering sample order
hc <- hclust(dist_counts, method = "complete")
hc_samples <- hc$labels[hc$order]

# =========================
# 6. PCA
# =========================
pca_res <- prcomp(t(stabilized_counts), scale. = FALSE)
pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- rownames(pca_df)
pca_df$condition <- coldata[rownames(pca_df), "condition"]

# Calculate percentage variance
percentVar <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Define pairs of PCs to plot
pc_pairs <- list(
  c(1, 2),
  c(1, 3),
  c(2, 3)
)

# Loop through PC pairs and save plots
for (pcs in pc_pairs) {
  x_pc <- pcs[1]
  y_pc <- pcs[2]
  
  p <- ggplot(
    pca_df,
    aes_string(x = paste0("PC", x_pc), y = paste0("PC", y_pc), color = "condition", label = "sample")
  ) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3) +
    scale_color_manual(values = condition_colors) +
    xlab(paste0("PC", x_pc, " (", round(percentVar[x_pc], 1), "%)")) +
    ylab(paste0("PC", y_pc, " (", round(percentVar[y_pc], 1), "%)")) +
    theme_classic() +
    ggtitle(paste0("PCA – Mouse RNA-seq: PC", x_pc, " vs PC", y_pc))
  
  # Save each PCA plot
  ggsave(
    filename = paste0("260428_output/20260130_PCA_PC", x_pc, "_PC", y_pc, ".png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

#Differential expression with out KZLX9

dds <- DESeqDataSetFromMatrix(
  countData = round(countData),
  colData = coldata,
  design = ~ condition
)

# ============================================================
# DESeq2 helper function
# ============================================================

run_deseq <- function(dds, ref, contrast) {
  dds$condition <- relevel(dds$condition, ref = ref)
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  res <- as.data.frame(res)
  res <- na.omit(res)
  res[, c("log2FoldChange", "padj")]
}

# ============================================================
# Common parameters
# ============================================================

log2FC_cutoff <- 0.5
padj_cutoff <- 0.05

# ============================================================
# ================= KIe – sh analysis =========================
# ============================================================
res_disease_e <- run_deseq(dds, "WT_empty", c("condition", "Stim1R304W/+_empty", "WT_empty"))
res_treated_e <- run_deseq(dds, "Stim1R304W/+_empty", c("condition", "Stim1R304W/+_sh190", "Stim1R304W/+_empty"))


head (res_disease_e)
head (res_treated_e)

# Filter significant disease genes (disease vs WT)
diseased_e <- res_disease_e[
  abs(res_disease_e$log2FoldChange) >= log2FC_cutoff &
    res_disease_e$padj < padj_cutoff, ]

# Find common genes between disease(WT) and treated datasets
common_genes_e <- intersect(rownames(diseased_e), rownames(res_treated_e))
diseased_matched_e <- diseased_e[common_genes_e, ]
treated_matched_e  <- res_treated_e[common_genes_e, ]

# Same order
treated_matched_e <- treated_matched_e[match(rownames(diseased_matched_e),
                                             rownames(treated_matched_e)), ]

# Log2FC definitions (MATCHING PROTEOMICS)
log2FC_disease_e <- diseased_matched_e$log2FoldChange              # disease vs WT  
log2FC_treated_e <- treated_matched_e$log2FoldChange                # treated vs disease
log2FC_treated_vs_WT_e <- log2FC_disease_e + log2FC_treated_e       # treated vs WT

# RNA-seq rescue metric: % reduction of disease effect (MATCHING PROTEOMICS)
rescue_metric_e <- 100 * (log2FC_disease_e - log2FC_treated_vs_WT_e) / log2FC_disease_e

# Create table (MATCHING PROTEOMICS STRUCTURE)
rescue_table_e <- data.frame(
  Gene               = rownames(diseased_matched_e),
  log2FC_disease     = log2FC_disease_e,
  log2FC_treated     = log2FC_treated_e,
  log2FC_treated_vs_WT = log2FC_treated_vs_WT_e,
  padj_disease       = diseased_matched_e$padj,
  padj_treated       = treated_matched_e$padj,
  rescue_metric      = rescue_metric_e
)

## ---- Fixed-threshold rescue categories ----
rescue_table_e$rescue_status_fixed <- cut(
  rescue_table_e$rescue_metric,
  breaks = c(-Inf, 0, 30, 80, 120, Inf),
  labels = c("Worsened", "Not rescued", "Partially rescued", "Rescued", "Over-corrected")
)
# Inspect (NO directionality filter)
cat("Number of significant disease genes (RNA):", nrow(diseased_e), "\n")
cat("Number of matched genes:", length(common_genes_e), "\n")
table(rescue_table_e$rescue_status_fixed)
head(rescue_table_e)

head (rescue_table_e)
library(org.Mm.eg.db)
library(clusterProfiler)

res_disease_e$GeneSymbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(res_disease_e),  # Your Ensembl IDs
  column = "SYMBOL",               # Want gene symbols
  keytype = "ENSEMBL",             # Input type
  multiVals = "first"              # If multiple symbols exist, take first
)
head (res_disease_e)


# Suppose your dataframe is called rescue table

rescue_table_e$GeneSymbol <- mapIds(org.Mm.eg.db,
                                    keys = rescue_table_e$Gene,
                                    column = "SYMBOL",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
write.csv(
  rescue_table_e,
  "260428_output/20260213_KIe_sh_try_rescue_table.csv",
  row.names = FALSE
)

# ============================================================
# ================= KIn – MOE analysis =======================
# ============================================================

res_disease_n <- run_deseq(dds, "WT_Nacl", c("condition", "Stim1R304W/+_Nacl", "WT_Nacl"))
res_treated_n <- run_deseq(dds, "Stim1R304W/+_Nacl", c("condition", "Stim1R304W/+_MOE", "Stim1R304W/+_Nacl"))
#res_sh = run_deseq(dds, "WT_empty", c("condition", "Stim1R304W/+_sh190", "WT_empty"))
#res_MOE = run_deseq(dds, "WT_Nacl", c("condition", "Stim1R304W/+_MOE", "WT_Nacl"))

diseased_n <- res_disease_n[
  abs(res_disease_n$log2FoldChange) >= log2FC_cutoff &
    res_disease_n$padj < padj_cutoff, ]
# Find common genes between disease(WT) and treated datasets
common_genes_n <- intersect(rownames(diseased_n), rownames(res_treated_n))
diseased_matched_n <- diseased_n[common_genes_n, ]
treated_matched_n  <- res_treated_n[common_genes_n, ]

# Same order (MATCHING PROTEOMICS)
treated_matched_n <- treated_matched_n[match(rownames(diseased_matched_n),
                                             rownames(treated_matched_n)), ]

# Log2FC definitions (MATCHING PROTEOMICS)
log2FC_disease_n <- diseased_matched_n$log2FoldChange              # disease vs WT  
log2FC_treated_n <- treated_matched_n$log2FoldChange                # treated vs disease
log2FC_treated_vs_WT_n <- log2FC_disease_n + log2FC_treated_n       # treated vs WT

# RNA-seq rescue metric: % reduction of disease effect (MATCHING PROTEOMICS)
rescue_metric_n <- 100 * (log2FC_disease_n - log2FC_treated_vs_WT_n) / log2FC_disease_n

# Create table (MATCHING PROTEOMICS STRUCTURE)
rescue_table_n <- data.frame(
  Gene               = rownames(diseased_matched_n),
  log2FC_disease     = log2FC_disease_n,
  log2FC_treated     = log2FC_treated_n,
  log2FC_treated_vs_WT = log2FC_treated_vs_WT_n,
  padj_disease       = diseased_matched_n$padj,
  padj_treated       = treated_matched_n$padj,
  rescue_metric      = rescue_metric_n
)

## ---- Fixed-threshold rescue categories ----
rescue_table_n$rescue_status_fixed <- cut(
  rescue_table_n$rescue_metric,
  breaks = c(-Inf, 0, 30, 80, 120, Inf),
  labels = c("Worsened", "Not rescued", "Partially rescued", "Rescued", "Over-corrected")
)
rescue_table_n$GeneSymbol <- mapIds(org.Mm.eg.db,
                                    keys = rescue_table_n$Gene,
                                    column = "SYMBOL",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
write.csv(
  rescue_table_n,
  "260428_output/20260213_KIn_MOE_try_rescue_table.csv",
  row.names = FALSE
)


# ============================================================
# ================= Treatment comparison =====================
# ============================================================

cat("\nKIe–sh fixed rescue summary:\n")
print(table(rescue_table_e$rescue_status_fixed))

cat("\nKIn–MOE fixed rescue summary:\n")
print(table(rescue_table_n$rescue_status_fixed))

cat("\nMedian rescue metric:\n")
cat("KIe–sh:", median(rescue_table_e$rescue_metric), "\n")
cat("KIn–MOE:", median(rescue_table_n$rescue_metric), "\n")

# ============================
# Plot distribution
# ============================

png(
  "260428_output/20260212_another_way_KI_sh_rescue_metric_distribution.png",
  width = 1200,
  height = 900,
  res = 150
)

hist(
  rescue_table_e$rescue_metric,
  breaks = 100,
  main = "Distribution of rescue metric_KIe_sh",
  xlab = "Rescue metric"
)

abline(v = 0, col = "red", lwd = 2)

dev.off()

png(
  "260428_output/20260212_another_way_KIn_MOE_rescue_metric_distribution.png",
  width = 1200,
  height = 900,
  res = 150
)

hist(
  rescue_table_n$rescue_metric,
  breaks = 100,
  main = "Distribution of rescue metric_KIn_MOE",
  xlab = "Rescue metric"
)

abline(v = 0, col = "red", lwd = 2)

dev.off()

#Comparison stack plot


library(dplyr)
library(ggplot2)

# ===============================
# Prepare data for comparison
# ===============================

# Filter only relevant statuses (exclude "Worsened")
rescue_table_sh_plot <- rescue_table_e %>%
  filter(rescue_status_fixed %in% c("Not rescued", "Partially rescued", "Rescued")) %>%
  select(rescue_status = rescue_status_fixed)

rescue_table_n_plot <- rescue_table_n %>%
  filter(rescue_status_fixed %in% c("Not rescued", "Partially rescued", "Rescued")) %>%
  select(rescue_status = rescue_status_fixed)

# Count percentages
rescue_summary_sh <- rescue_table_sh_plot %>%
  count(rescue_status) %>%
  mutate(percent = n / sum(n) * 100,
         Treatment = "Stim1R304W/+_sh190")
print (rescue_summary_sh)

rescue_summary_MOE <- rescue_table_n_plot %>%
  count(rescue_status) %>%
  mutate(percent = n / sum(n) * 100,
         Treatment = "Stim1R304W/+_MOE")

# Combine both treatments
combined_rescue <- bind_rows(rescue_summary_sh, rescue_summary_MOE)
unique(combined_rescue$Treatment)
combined_rescue$Treatment <- factor(combined_rescue$Treatment, 
                                    levels = c("Stim1R304W/+_sh190", "Stim1R304W/+_MOE"))
print (combined_rescue)
# ===============================
# Plot stacked bar
# ===============================

plot <- ggplot(combined_rescue, aes(x = Treatment, y = percent, fill = rescue_status)) +
  geom_bar(stat = "identity", width = 0.7) +
  #geom_text(aes(label = paste0(round(percent, 1), "%")), 
  #         position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("Not rescued" = "lightgrey", 
                               "Partially rescued" = "orange",
                               "Rescued" = "forestgreen")) +
  labs(x = NULL, y = "% of Dysregulated Genes", fill = "Rescue Status") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "right")


# Show plot
print(plot)
# Save
ggsave("260428_output/20260212_another_try_stack_plot_MOE_sh_rescue_pattern.png", 
       plot = plot, width = 7, height = 5, dpi = 300)

#Venn diagram
# ===============================
# Libraries
# ===============================
library(dplyr)
library(VennDiagram)

# ===============================
# Select rescued genes for each treatment
# ===============================

# For sh treatment: include both "Rescued" and "Partially rescued"
sh_genes_disease = res_disease_e$GeneSymbol
sh_genes <- rescue_table_e %>%
  filter(rescue_status_fixed %in% c("Rescued", "Partially rescued")) %>%
  pull(Gene)
bg_genes_sh = rescue_table_e %>% filter(rescue_status_fixed %in% c ("Rescued", "Partially rescued", "Not rescued", "Worsened", "Over-corrected")) %>% pull (Gene)
length(bg_genes_MOE)
# For MOE treatment: include both "Rescued" and "Partially rescued"
MOE_genes <- rescue_table_n %>%
  filter(rescue_status_fixed %in% c("Rescued", "Partially rescued")) %>%
  pull(Gene)
bg_genes_MOE = rescue_table_n %>% filter(rescue_status_fixed %in% c ("Rescued", "Partially rescued", "Not rescued", "Worsened", "Over-corrected")) %>% pull (Gene)

# ===============================
# Identify overlaps and uniques
# ===============================
only_sh <- setdiff(sh_genes, MOE_genes)
only_MOE <- setdiff(MOE_genes, sh_genes)
both <- intersect(sh_genes, MOE_genes)
print (both)
cat("Genes rescued only by sh:", length(only_sh), "\n")
cat("Genes rescued only by MOE:", length(only_MOE), "\n")
cat("Genes rescued by both treatments:", length(both), "\n")

# ===============================
# Create Venn list
# ===============================
venn_list <- list(
  "Stim1R304W/+_sh190" = sh_genes,
  "Stim1R304W/+_MOE" = MOE_genes
)

# ===============================
# Plot Venn diagram
# ===============================
library(VennDiagram)
library(grid)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("darkorange2", "gold1"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = c(-10, 10),
  cat.dist = c(0.05, 0.05),
  main = "Overlap of Rescued Genes"
)

# Draw the Venn
grid.newpage()
grid.draw(venn.plot)

# 🔑 Identify and move the intersection label ("6")
for (i in seq_along(venn.plot)) {
  if (inherits(venn.plot[[i]], "text") &&
      venn.plot[[i]]$label == "6") {
    
    venn.plot[[i]]$x <- unit(0.5, "npc")
    venn.plot[[i]]$y <- unit(0.5, "npc")
  }
}

# Redraw after moving
grid.newpage()
grid.draw(venn.plot)


# ===============================
# Save as PNG
# ===============================
png("260428_output/20260212_try_rescue_overlap_venn.png", width = 1000, height = 1000, res = 150)
grid::grid.draw(venn.plot)
dev.off()

#save the list
intersection_genes <- intersect(
  venn_list[["Stim1R304W/+_sh190"]],
  venn_list[["Stim1R304W/+_MOE"]]
)

write.csv(intersection_genes,
          file = "260428_output/20260213_intersection_rescued_genes_venn.csv",
          row.names = FALSE)



#volcano plot _significantly dysregulated genes


library(AnnotationDbi)
library(org.Mm.eg.db)  # mouse annotation
library(ggplot2)
library(dplyr)
library (ggrepel)
# ============================
# Function to convert Ensembl IDs to Gene Symbols
# ============================
convert_ensembl_to_symbol <- function(df) {
  df$GeneSymbol <- mapIds(
    org.Mm.eg.db,
    keys = rownames(df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  df <- df[!is.na(df$GeneSymbol), ]  # remove rows without gene symbol
  df
}

# ============================
# Convert Ensembl IDs
# ============================
res_disease_e <- convert_ensembl_to_symbol(res_disease_e)
res_treated_e <- convert_ensembl_to_symbol(res_treated_e)
#res_sh = convert_ensembl_to_symbol(res_sh)
#res_MOE = convert_ensembl_to_symbol(res_MOE)
# ============================
# Filter significant genes
# ============================
log2FC_cutoff <- 0.5
padj_cutoff <- 0.05

sig_disease <- res_disease_e %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < padj_cutoff)
sig_treated <- res_treated_e %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < padj_cutoff)

plot_volcano_labeled <- function(df, title) {
  df <- df %>%
    mutate(Significance = case_when(
      padj < padj_cutoff & log2FoldChange > log2FC_cutoff ~ "Up",
      padj < padj_cutoff & log2FoldChange < -log2FC_cutoff ~ "Down",
      TRUE ~ "NS"
    ))
  
  # Select top 10 up and down genes
  top_up <- df %>%
    filter(Significance == "Up") %>%
    arrange(padj) %>%
    head(10)
  
  top_down <- df %>%
    filter(Significance == "Down") %>%
    arrange(padj) %>%
    head(10)
  
  top_genes <- bind_rows(top_up, top_down)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_classic(base_size = 14) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)", color = "DE status") +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = top_genes,
      aes(label = GeneSymbol),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 20
    )
}

# ============================
# Volcano plots with labels
# ============================

volc_disease_labeled <- plot_volcano_labeled(res_disease_e, "KIe vs WTe (Disease effect)")
volc_treated_labeled <- plot_volcano_labeled(res_treated_e, "sh vs KIe (Treatment effect)")
volc_res_sh = plot_volcano_labeled(res_sh, "sh vs WTe (treatment effect with WT)")
volc_res_MOE = plot_volcano_labeled(res_MOE, "MOE vs WTe (treatment effect with WT)")
print(volc_res_sh)
print(volc_res_MOE)
print (volc_disease_labeled)
print (volc_treated_labeled)
# ============================
# Save plots
# ============================

ggsave("260428_output/20260202_volcano_disease_labeled.png", plot = volc_disease_labeled, width = 7, height = 6, dpi = 300)
ggsave("260428_output/20260202_volcano_treated_labeled.png", plot = volc_treated_labeled, width = 7, height = 6, dpi = 300)

res_disease_e$significance <- ifelse(
  res_disease_e$padj < 0.05,
  "significant",
  "non-significant"
)
res_treated_e$significance <- ifelse(
  res_treated_e$padj < 0.05,
  "significant",
  "non-significant"
)
head (res_treated_e)
write.csv(res_disease_e, file = "260428_output/20260202_DEGs_disease_KIe_WTe.csv", sep = "\t")
write.csv(res_treated_e, file = "260428_output/20260202_DEGs_treated_sh_KIe.csv", sep = "\t")


#KIn -MOE
library(AnnotationDbi)
library(org.Mm.eg.db)  # mouse annotation
library(ggplot2)
library(dplyr)
library (ggrepel)
# ============================
# Function to convert Ensembl IDs to Gene Symbols
# ============================
convert_ensembl_to_symbol <- function(df) {
  df$GeneSymbol <- mapIds(
    org.Mm.eg.db,
    keys = rownames(df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  df <- df[!is.na(df$GeneSymbol), ]  # remove rows without gene symbol
  df
}

# ============================
# Convert Ensembl IDs
# ============================
res_disease_n <- convert_ensembl_to_symbol(res_disease_n)
res_treated_n <- convert_ensembl_to_symbol(res_treated_n)

# ============================
# Filter significant genes
# ============================
log2FC_cutoff <- 0.5
padj_cutoff <- 0.05

sig_disease_n <- res_disease_n %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < padj_cutoff)

sig_treated_n <- res_treated_n %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < padj_cutoff)

plot_volcano_labeled <- function(df, title) {
  df <- df %>%
    mutate(Significance = case_when(
      padj < padj_cutoff & log2FoldChange > log2FC_cutoff ~ "Up",
      padj < padj_cutoff & log2FoldChange < -log2FC_cutoff ~ "Down",
      TRUE ~ "NS"
    ))
  
  # Select top 10 up and down genes
  top_up <- df %>%
    filter(Significance == "Up") %>%
    arrange(padj) %>%
    head(10)
  
  top_down <- df %>%
    filter(Significance == "Down") %>%
    arrange(padj) %>%
    head(10)
  
  top_genes <- bind_rows(top_up, top_down)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_classic(base_size = 14) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)", color = "DE status") +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = top_genes,
      aes(label = GeneSymbol),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 20
    )
}

# ============================
# Volcano plots with labels
# ============================

volc_disease_labeled_n <- plot_volcano_labeled(res_disease_n, "KIn vs WTn (Disease effect)")
volc_treated_labeled_n <- plot_volcano_labeled(res_treated_n, "MOE vs KIn (Treatment effect)")

print(volc_disease_labeled_n)
print(volc_treated_labeled_n)

head (res_treated_n)
# ============================
# Save plots
# ============================

ggsave("260428_output/20260202_volcano_disease_labeled_KIn_WTn.png", plot = volc_disease_labeled_n, width = 7, height = 6, dpi = 300)
ggsave("260428_output/20260202_volcano_treated_labeled_MOE_KIn.png", plot = volc_treated_labeled_n, width = 7, height = 6, dpi = 300)

res_disease_n$significance <- ifelse(
  res_disease_n$padj < 0.05,
  "significant",
  "non-significant"
)
res_treated_n$significance <- ifelse(
  res_treated_n$padj < 0.05,
  "significant",
  "non-significant"
)

write.csv(res_disease_n, file = "260428_output/20260202_DEGs_disease_KIn_WTn.csv", sep = "\t")
write.csv(res_treated_n, file = "260428_output/20260202_DEGs_treated_MOE_KIn.csv", sep = "\t")


# ------------------------------
# 0. Load libraries
# ------------------------------
library(DESeq2)
library(tidyverse)
library(ggpubr)

# ------------------------------
# 1. Define genes of interest
# ------------------------------
genes_of_interest <- c(
  "ENSMUSG00000026864",  # Hspa5
  "ENSMUSG00000020908",  # Myh3
  "ENSMUSG00000055775",  # Myh8
  "ENSMUSG00000049686",  # Orai1
  "ENSMUSG00000030987",  # Stim1
  "ENSMUSG00000030592",  # Ryr1
  "ENSMUSG00000030730",  # Serca1 / Atp2a1
  "ENSMUSG00000042045",   # Sln
  "ENSMUSG00000020484"
)

gene_names <- c(
  "Hspa5", "Myh3", "Myh8",
  "Orai1", "Stim1", "Ryr1", "Serca1", "Sln", "Xbp1"
)

# ------------------------------
# 2. Create DESeqDataSet and normalize
# ------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(countData),
  colData   = coldata,
  design    = ~ condition
)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# ------------------------------
# 3. Extract genes of interest
# ------------------------------
selected_counts <- norm_counts[genes_of_interest, ]
rownames(selected_counts) <- gene_names

# ------------------------------
# 4. Prepare long-format dataframe
# ------------------------------
coldata2 <- coldata %>%
  rownames_to_column("Sample") %>%
  rename(Condition = condition)

condition_order <- c(
  "WT_Nacl", "WT_empty",
  "Stim1R304W/+_Nacl", "Stim1R304W/+_empty",
  "Stim1R304W/+_MOE", "Stim1R304W/+_sh190"
)

long_df <- selected_counts %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(coldata2, by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = condition_order))

# ------------------------------
# 5. Mean ± SEM summary dataframe
# ------------------------------
summary_df <- long_df %>%
  group_by(Gene, Condition) %>%
  summarise(
    mean = mean(Expression),
    sem  = sd(Expression) / sqrt(n()),
    .groups = "drop"
  )

# ------------------------------
# 6. Define Tukey comparison pairs
# ------------------------------
my_pairs <- tribble(
  ~group1, ~group2,
  "WT_empty", "Stim1R304W/+_empty",
  "WT_empty", "Stim1R304W/+_sh190",
  "Stim1R304W/+_sh190", "Stim1R304W/+_empty",
  "WT_Nacl", "Stim1R304W/+_Nacl",
  "WT_Nacl", "Stim1R304W/+_MOE",
  "Stim1R304W/+_MOE", "Stim1R304W/+_Nacl"
) %>%
  mutate(pair = paste(pmin(group1, group2), pmax(group1, group2), sep = "_"))

# ------------------------------
# 7. Run one-way ANOVA + Tukey HSD
# ------------------------------
tukey_list <- list()

for (gene in gene_names) {
  
  df_gene <- long_df %>% filter(Gene == gene)
  fit <- aov(Expression ~ Condition, data = df_gene)
  tuk <- TukeyHSD(fit)
  
  tuk_df <- as.data.frame(tuk$Condition) %>%
    rownames_to_column("comparison") %>%
    separate(comparison, into = c("group1","group2"), sep = "-") %>%
    mutate(
      Gene = gene,
      pair = paste(pmin(group1, group2), pmax(group1, group2), sep = "_")
    ) %>%
    rename(p.adj = `p adj`)
  
  tukey_list[[gene]] <- tuk_df
}

tukey_results <- bind_rows(tukey_list) %>%
  semi_join(my_pairs, by = "pair")

# ------------------------------
# 8. Prepare Tukey plot dataframe
# ------------------------------
tukey_plot_df <- tukey_results %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  group_by(Gene) %>%
  mutate(
    y.position = max(long_df$Expression[long_df$Gene == unique(Gene)]) *
      seq(1.15, 1.45, length.out = n())
  ) %>%
  ungroup() %>%
  mutate(
    group1 = factor(group1, levels = condition_order),
    group2 = factor(group2, levels = condition_order)
  )

# ------------------------------
# 9. Define color scheme
# ------------------------------
condition_colors <- condition_colors

# ------------------------------
# 10. Plot function (mean ± SEM)
# ------------------------------
plot_gene <- function(gene_name) {
  
  p <- ggplot(long_df %>% filter(Gene == gene_name),
              aes(x = Condition, y = Expression, fill = Condition)) +
    
    geom_jitter(width = 0.12, size = 3, shape = 24) +
    
    geom_errorbar(
      data = summary_df %>% filter(Gene == gene_name),
      aes(x = Condition, y = mean, ymin = mean - sem, ymax = mean + sem),
      width = 0.25,
      size = 0.8,
      inherit.aes = FALSE
    ) +
    
    geom_crossbar(
      data = summary_df %>% filter(Gene == gene_name),
      aes(x = Condition, y = mean, ymin = mean, ymax = mean),
      width = 0.6,
      size = 0.4,
      color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = condition_colors) +
    theme_classic(base_size = 14) +
    labs(title = gene_name, y = "Normalized expression", x = NULL) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  tukey_gene <- tukey_plot_df %>% filter(Gene == gene_name)
  
  if (nrow(tukey_gene) > 0) {
    p <- p + stat_pvalue_manual(
      tukey_gene,
      label = "p.adj.signif",
      tip.length = 0.01,
      size = 5
    )
  }
  
  return(p)
}

# ------------------------------
# 11. Generate and save plots
# ------------------------------
plots <- lapply(gene_names, plot_gene)


for (i in seq_along(gene_names)) {
  ggsave(
    filename = paste0("", gene_names[i], "260428_output/_20260202_plot.png"),
    plot = plots[[i]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

#GO enrichment
length (sh_genes)
library(clusterProfiler)
library(org.Mm.eg.db)  # mouse annotation
library(dplyr)
head (sh_genes)
sh_genes_disease_entrez <- mapIds(org.Mm.eg.db, keys = sh_genes_disease, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
sh_genes_entrez <- mapIds(org.Mm.eg.db, keys = sh_genes, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
print (MOE_genes_entrez)
MOE_genes_entrez <- mapIds(org.Mm.eg.db, keys = MOE_genes, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
ego_rescued_sh <- enrichGO(
  gene = sh_genes_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",       # BP = Biological Process, you can also use "MF" or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
ego_rescued_MOE <- enrichGO(
  gene = MOE_genes_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
ego_disease_sh <- enrichGO(
  gene = sh_genes_disease_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",       # BP = Biological Process, you can also use "MF" or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2, readable = TRUE
)
head (ego_rescued_MOE)
p1 = barplot(ego_rescued_sh, showCategory = 35, title = "Enriched GO-BP terms by sh treatment", font.size = 8)
p2 = barplot(ego_rescued_MOE, showCategory = 35, title = "Enriched GO-BP terms by MOE treatment", font.size = 7)
p3 = barplot(ego_rescued_MOE, showCategory = 35, title = "Enriched GO-BP terms by MOE treatment_0.5", font.size = 7)
p4= barplot(ego_disease_sh, showCategory = 35, title = "Enriched GO-BP terms by sh disease", font.size = 7)
write.csv (ego_rescued_sh, file = "260428_output/20260213_try_sh_rescued_0.05_GO_BP.csv")
write.csv (ego_rescued_MOE, file = "260428_output/20260213_try_MOE_rescued_0.o5_GO_BP.csv")
print (p1)
print (p2)
print (p4)
write.csv (ego_disease_sh, file = "20260323_sh_disease_GO_BP.csv")
ggsave (filename = "260428_output/20260213_try_sh_rescued_0.05_GO_BP.png", plot = p1, width = 15, height = 10, dpi = 300)
ggsave (filename = "260428_output/20260202_try_MOE_rescued_0.05_GO_BP.png", plot = p2, width = 15, height = 10, dpi = 300)
ggsave (filename = "260428_output/20260202_MOE_rescued_0.5_GO_BP.png", plot = p3, width = 15, height = 10, dpi = 300)


run_GO_all_ont <- function(entrez_ids, label) {
  
  ego_BP <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  ego_MF <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  ego_CC <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  list(BP = ego_BP, MF = ego_MF, CC = ego_CC, label = label)
}


go_sh   <- run_GO_all_ont(sh_genes_entrez, "Genes rescued by sh")
go_MOE  <- run_GO_all_ont(MOE_genes_entrez, "Genes rescued by MOE")
library(dplyr)

combine_GO_results <- function(go_list_object) {
  
  bp_df <- as.data.frame(go_list_object$BP) %>% 
    mutate(Ontology = "BP")
  
  mf_df <- as.data.frame(go_list_object$MF) %>% 
    mutate(Ontology = "MF")
  
  cc_df <- as.data.frame(go_list_object$CC) %>% 
    mutate(Ontology = "CC")
  
  bind_rows(bp_df, mf_df, cc_df)
}
go_sh_df   <- combine_GO_results(go_sh)
go_MOE_df   <- combine_GO_results(go_MOE)

library(ggplot2)

plot_GO_combined <- function(go_df, title_text) {
  ggplot(go_df, aes(x = reorder(Description, Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~Ontology, , ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("BP" = "steelblue", 
                                 "MF" = "darkorange", 
                                 "CC" = "forestgreen")) +
    labs(title = title_text,
         x = "GO Term",
         y = "Gene Count") +
    theme_classic(base_size = 12)
}
p_sh_combined <- plot_GO_combined(go_sh_df, 
                                  "GO Enrichment (BP/MF/CC): Genes rescued by sh")

p_MOE_combined <- plot_GO_combined(go_MOE_df, 
                                   "GO Enrichment (BP/MF/CC): Genes rescued by MOE")


p_sh_combined
p_MOE_combined
ggsave("260428_output/20260213_GO_sh_combined.png", p_sh_combined, width = 12, height = 10, dpi = 300)
ggsave("260428_output/20260213_GO_MOE_combined.png", p_MOE_combined, width = 12, height = 10, dpi = 300)

write.csv(as.data.frame(go_sh_df), "260428_output/20260213_GO_sh_combined.csv", row.names = FALSE)
write.csv(as.data.frame(go_MOE_df), "260428_output/20260213_GO_MOE_combined.csv", row.names = FALSE)


#Reactome analysis
library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

sh_genes_entrez <- na.omit(sh_genes_entrez)
MOE_genes_entrez <- na.omit(MOE_genes_entrez)
reactome_sh <- enrichPathway(
  gene = sh_genes_entrez,
  organism = "mouse",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)
reactome_MOE <- enrichPathway(
  gene = MOE_genes_entrez,
  organism = "mouse",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)
reactome_sh_df  <- as.data.frame(reactome_sh)
reactome_MOE_df <- as.data.frame(reactome_MOE)
p_react_sh <- barplot(reactome_sh, 
                      showCategory = 30, 
                      title = "Reactome Pathways - sh treatment",
                      font.size = 8)

p_react_MOE <- barplot(reactome_MOE, 
                       showCategory = 30, 
                       title = "Reactome Pathways - MOE treatment",
                       font.size = 8)
print (p_react_MOE)
ggsave("260428_output/20260213_reactome_sh.png", plot = p_react_sh, width = 7, height = 6, dpi = 300)
ggsave("260428_output/20260213_reactome_MOE.png", plot = p_react_MOE, width = 7, height = 6, dpi = 300)

#GSEA 
#didnt work

#combined partial + complete
# For sh treatment
rescued_sh_combined <- rescue_table_e %>%
  filter(rescue_status_fixed %in% c("Rescued", "Partially rescued"))
# For MOE treatment
rescued_MOE_combined <- rescue_table_n %>%
  filter(rescue_status_fixed %in% c("Rescued", "Partially rescued"))
# sh
ranked_sh <- rescued_sh_combined$log2FC_treated
head (rescued_sh_combined)

names(ranked_sh) <- rescued_sh_combined$Gene
ranked_sh <- sort(ranked_sh, decreasing = TRUE)
# MOE
ranked_MOE <- rescued_MOE_combined$log2FC_treated
names(ranked_MOE) <- rescued_MOE_combined$Gene
ranked_MOE <- sort(ranked_MOE, decreasing = TRUE)
gsea_sh <- gseGO(
  geneList = ranked_sh,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",          # Biological Process
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# GSEA for MOE treatment
gsea_MOE <- gseGO(
  geneList = ranked_MOE,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)
ridgeplot(gsea_sh)
ridgeplot(gsea_MOE)

gsea_diff <- gseGO(
  geneList = ranked_diff,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

dotplot(gsea_diff, showCategory = 15, title = "GSEA - MOE vs sh (rescued + partially rescued)")

#upset plot


library(dplyr)

rescue_table_e <- rescue_table_e %>%
  mutate(
    rescue_status = case_when(
      rescue_status_fixed == "Rescued" ~ "Rescued",
      rescue_status_fixed == "Partially rescued" ~ "Partially rescued",
      rescue_status_fixed == "Not rescued" ~ "Not rescued",
      TRUE ~ NA_character_
    )
  )

rescue_table_n <- rescue_table_n %>%
  mutate(
    rescue_status = case_when(
      rescue_status_fixed == "Rescued" ~ "Rescued",
      rescue_status_fixed == "Partially rescued" ~ "Partially rescued",
      rescue_status_fixed == "Not rescued" ~ "Not rescued",
      TRUE ~ NA_character_
    )
  )
# sh treatment
sh_rescued    <- rescue_table_e %>% filter(rescue_status == "Rescued") %>% pull(Gene)
sh_partial    <- rescue_table_e %>% filter(rescue_status == "Partially rescued") %>% pull(Gene)
sh_not        <- rescue_table_e %>% filter(rescue_status == "Not rescued") %>% pull(Gene)

# MOE treatment
moe_rescued   <- rescue_table_n %>% filter(rescue_status == "Rescued") %>% pull(Gene)
moe_partial   <- rescue_table_n %>% filter(rescue_status == "Partially rescued") %>% pull(Gene)
moe_not       <- rescue_table_n %>% filter(rescue_status == "Not rescued") %>% pull(Gene)
upset_list <- list(
  "sh_Rescued"        = sh_rescued,
  "sh_Partially"     = sh_partial,
  "sh_Not_rescued"   = sh_not,
  "MOE_Rescued"      = moe_rescued,
  "MOE_Partially"    = moe_partial,
  "MOE_Not_rescued"  = moe_not
)
library(UpSetR)
png (filename = "260428_output/20260213_Upsetplot_sh_MOE_rescue_status.png", width = 3000, height = 2000, res = 300)

upset(
  fromList(upset_list),
  sets = names(upset_list),
  order.by = "freq",
  empty.intersections = "on",
  mainbar.y.label = "Gene overlap",
  sets.x.label = "Genes per category"
)
dev.off ()
saveRDS(upset_list, file = "UpSet_gene_lists_sh_vs_MOE.rds")
upset_matrix <- fromList(upset_list)

write.csv(
  upset_matrix,
  file = "260428_output/UpSet_binary_matrix_sh_vs_MOE.csv"
)





#complex heatmap
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)
head (stabilized_counts)
ca_genes <- read.csv("entire_pathway.csv", header = TRUE)
colnames(ca_genes) <- c("Ensembl", "Gene", "SubPathway", "Pathway")
colnames(ca_genes)
# Clean column names
genes_of_interest <- ca_genes$Ensembl
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(stabilized_counts)]
# Subset expression matrix
mat <- stabilized_counts[genes_of_interest, ]

# Log transform
mat <- log2(mat + 1)

# Scale per gene
mat <- t(scale(t(mat)))
desired_order <- c(
  "WT_empty",
  "Stim1R304W/+_empty",
  "Stim1R304W/+_sh190", 
  "WT_Nacl",
  "Stim1R304W/+_Nacl",
  "Stim1R304W/+_MOE"
)

coldata$condition <- factor(coldata$condition, levels = desired_order)
sample_order <- rownames(coldata[order(coldata$condition), ])

mat <- mat[, sample_order]
ca_genes <- ca_genes %>%
  filter(Ensembl %in% genes_of_interest)

# Ensure same order as matrix
ca_genes <- ca_genes[match(rownames(mat), ca_genes$Ensembl), ]
print (ca_genes)
ca_genes$Pathway <- factor(ca_genes$Pathway,
                           levels = c("Ca2+ handling", "Cell stress", "Mitochondria"))


ca_genes$SubPathway <- factor(ca_genes$SubPathway,
                              levels = c(
                                "Ca2+ extrusion",
                                "Sr refilling",
                                "EC coupling",
                                "Ca2+ signalling",
                                "Biogenesis",
                                "Migration",
                                "Fission/Fusion",
                                "UPR",
                                "Apoptosis",
                                "Regeneration"
                              ))
order_idx <- order(ca_genes$Pathway, ca_genes$SubPathway)

mat <- mat[order_idx, ]
ca_genes <- ca_genes[order_idx, ]
row_split <- ca_genes$Pathway
library(RColorBrewer)

row_ha <- rowAnnotation(
  Pathway = ca_genes$Pathway,
  SubPathway = ca_genes$SubPathway,
  
  col = list(
    Pathway = c(
      "Ca2+ handling" = "#1b9e77",
      "Mitochondria" = "#d95f02",
      "Cell stress" = "#7570b3"
    ),
    
    SubPathway = structure(
      colorRampPalette(brewer.pal(8, "Set3"))(length(levels(ca_genes$SubPathway))),
      names = levels(ca_genes$SubPathway)
    )
  )
)
col_ha <- HeatmapAnnotation(
  Condition = coldata[sample_order, "condition"],
  col = list(Condition = condition_colors)
)
rownames(mat) <- ca_genes$Gene
library(ComplexHeatmap)
library(circlize)

ht <- Heatmap(
  mat,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  top_annotation = col_ha,
  left_annotation = row_ha,
  row_split = row_split,
  
  cluster_rows = FALSE,          # <-- important to keep order
  cluster_row_slices = FALSE,    # <-- do not reorder within splits
  cluster_columns = FALSE,
  
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_side = "right",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 80,
  row_gap = unit(1, "mm"),
  column_title = "Conditions",
  row_title = "Calcium & Mitochondrial Pathways"
)
ht

png("260428_output/260323_heatmap_GOI.png",
    width = 3000, height = 1800, res = 300)

draw(ht)

dev.off()


#SPLIT THE HEATMAP
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(grid)

### -----------------------------
### MATRIX PREP
### -----------------------------

genes_of_interest <- ca_genes$Ensembl
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(stabilized_counts)]

mat <- stabilized_counts[genes_of_interest, ]
mat <- log2(mat + 1)
mat <- t(scale(t(mat)))

desired_order <- c(
  "WT_empty",
  "Stim1R304W/+_empty",
  "Stim1R304W/+_sh190", 
  "WT_Nacl",
  "Stim1R304W/+_Nacl",
  "Stim1R304W/+_MOE"
)

coldata$condition <- factor(coldata$condition, levels = desired_order)
sample_order <- rownames(coldata[order(coldata$condition), ])
mat <- mat[, sample_order]

ca_genes <- ca_genes %>%
  filter(Ensembl %in% genes_of_interest)

ca_genes <- ca_genes[match(rownames(mat), ca_genes$Ensembl), ]

# Factor levels
ca_genes$Pathway <- factor(ca_genes$Pathway,
                           levels = c("Ca2+ handling", "Mitochondria", "Cell stress"))

ca_genes$SubPathway <- factor(ca_genes$SubPathway,
                              levels = c(
                                "Ca2+ extrusion",
                                "Sr refilling",
                                "EC coupling",
                                "Ca2+ signalling",
                                "Biogenesis",
                                "Migration",
                                "Fission/Fusion",
                                "UPR",
                                "Apoptosis",
                                "Regeneration"
                              ))

# Order rows
order_idx <- order(ca_genes$Pathway, ca_genes$SubPathway)
mat <- mat[order_idx, ]
ca_genes <- ca_genes[order_idx, ]

rownames(mat) <- ca_genes$Gene
row_split <- ca_genes$Pathway

### -----------------------------
### ROW ANNOTATION (NO FIXED WIDTH)
### -----------------------------

row_ha <- rowAnnotation(
  Pathway = ca_genes$Pathway,
  SubPathway = ca_genes$SubPathway,
  col = list(
    Pathway = c(
      "Ca2+ handling" = "#1b9e77",
      "Mitochondria" = "#d95f02",
      "Cell stress" = "#7570b3"
    ),
    SubPathway = structure(
      colorRampPalette(brewer.pal(8, "Set3"))(length(levels(ca_genes$SubPathway))),
      names = levels(ca_genes$SubPathway)
    )
  )
)

### -----------------------------
### SPLIT GROUPS
### -----------------------------

group1 <- c("WT_empty",
            "Stim1R304W/+_empty",
            "Stim1R304W/+_sh190")

group2 <- c("WT_Nacl",
            "Stim1R304W/+_Nacl",
            "Stim1R304W/+_MOE")
deg_g1 <- read.csv("20260202_DEGs_disease_KIe_WTe.csv", header = TRUE)
deg_g2 <- read.csv("20260202_DEGs_disease_KIn_WTn.csv", header = TRUE)
samples_g1 <- rownames(coldata[coldata$condition %in% group1, ])
samples_g2 <- rownames(coldata[coldata$condition %in% group2, ])

samples_g1 <- samples_g1[order(factor(coldata[samples_g1, "condition"], levels = group1))]
samples_g2 <- samples_g2[order(factor(coldata[samples_g2, "condition"], levels = group2))]

mat_g1 <- mat[, samples_g1]
mat_g2 <- mat[, samples_g2]

# CRITICAL: enforce same row order
row_order_master <- rownames(mat)
mat_g1 <- mat_g1[row_order_master, ]
mat_g2 <- mat_g2[row_order_master, ]

row_split <- factor(row_split,
                    levels = c("Ca2+ handling", "Mitochondria", "Cell stress"))

### -----------------------------
### COLUMN ANNOTATIONS
### -----------------------------

col_ha_g1 <- HeatmapAnnotation(
  Condition = coldata[samples_g1, "condition"],
  col = list(Condition = condition_colors)
)

col_ha_g2 <- HeatmapAnnotation(
  Condition = coldata[samples_g2, "condition"],
  col = list(Condition = condition_colors)
)

### -----------------------------
### HEATMAP 1 (FULL)
### -----------------------------
head (deg_g1)
# For group 1
signif_g1 <- ifelse(rownames(mat_g1) %in% deg_g1$GeneSymbol[deg_g1$significance == "significant"],
                    "yes", "no")

# For group 2
signif_g2 <- ifelse(rownames(mat_g2) %in% deg_g2$GeneSymbol[deg_g2$significance == "significant"],
                    "yes", "no")
library(grid)

# Example: red and larger font for significant genes
row_gp_g1 <- grid::gpar(
  fontsize = ifelse(signif_g1 == "yes", 9, 8),  # larger for significant
  col = ifelse(signif_g1 == "yes", "darkgreen", "black")
)

row_gp_g2 <- grid::gpar(
  fontsize = ifelse(signif_g2 == "yes", 9, 8),
  col = ifelse(signif_g2 == "yes", "darkgreen", "black")
)
ht1 <- Heatmap(
  mat_g1,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  top_annotation = col_ha_g1,
  left_annotation = row_ha,
  row_split = row_split,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_names_side = "right",
  show_row_names = TRUE,
  row_names_gp = row_gp_g1,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Empty / sh190",
  column_names_rot = 80
)
ht1

### -----------------------------
### HEATMAP 2 (FULL)
### -----------------------------

ht2 <- Heatmap(
  mat_g2,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  top_annotation = col_ha_g2,
  left_annotation = row_ha,
  row_split = row_split,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = row_gp_g2,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "NaCl / MOE",
  column_names_rot = 80
)
ht2
### -----------------------------
### DRAW SEPARATELY
### -----------------------------

png("260428_output/heatmap_group1.png", width = 2000, height = 1800, res = 300)
draw(ht1)
dev.off()

png("260428_output/heatmap_group2.png", width = 2000, height = 1800, res = 300)
draw(ht2)
dev.off()


#different colours_up_down

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(grid)

### -----------------------------
### MATRIX PREP
### -----------------------------

# Keep only genes of interest in your stabilized_counts
genes_of_interest <- ca_genes$Ensembl
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(stabilized_counts)]

mat <- stabilized_counts[genes_of_interest, ]
mat <- log2(mat + 1)
mat <- t(scale(t(mat)))

# Desired order of conditions
desired_order <- c(
  "WT_empty",
  "Stim1R304W/+_empty",
  "Stim1R304W/+_sh190", 
  "WT_Nacl",
  "Stim1R304W/+_Nacl",
  "Stim1R304W/+_MOE"
)

coldata$condition <- factor(coldata$condition, levels = desired_order)
sample_order <- rownames(coldata[order(coldata$condition), ])
mat <- mat[, sample_order]

# Filter ca_genes and match row order
ca_genes <- ca_genes %>% filter(Ensembl %in% genes_of_interest)
ca_genes <- ca_genes[match(rownames(mat), ca_genes$Ensembl), ]

# Factor levels for pathway/subpathway
ca_genes$Pathway <- factor(ca_genes$Pathway, levels = c("Ca2+ handling", "Mitochondria", "Cell stress"))
ca_genes$SubPathway <- factor(ca_genes$SubPathway,
                              levels = c(
                                "Ca2+ extrusion","Sr refilling","EC coupling","Ca2+ signalling",
                                "Biogenesis","Migration","Fission/Fusion","UPR","Apoptosis","Regeneration"
                              ))

# Order rows by Pathway and SubPathway
order_idx <- order(ca_genes$Pathway, ca_genes$SubPathway)
mat <- mat[order_idx, ]
ca_genes <- ca_genes[order_idx, ]
rownames(mat) <- ca_genes$Gene
row_split <- ca_genes$Pathway

### -----------------------------
### ROW ANNOTATION
### -----------------------------

row_ha <- rowAnnotation(
  Pathway = ca_genes$Pathway,
  SubPathway = ca_genes$SubPathway,
  col = list(
    Pathway = c(
      "Ca2+ handling" = "#1b9e77",
      "Mitochondria" = "#d95f02",
      "Cell stress" = "#7570b3"
    ),
    SubPathway = structure(
      colorRampPalette(brewer.pal(8, "Set3"))(length(levels(ca_genes$SubPathway))),
      names = levels(ca_genes$SubPathway)
    )
  )
)

### -----------------------------
### SPLIT GROUPS AND READ DEG FILES
### -----------------------------

group1 <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190")
group2 <- c("WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")

deg_g1 <- read.csv("20260202_DEGs_disease_KIe_WTe.csv", header = TRUE)
deg_g2 <- read.csv("20260202_DEGs_disease_KIn_WTn.csv", header = TRUE)

samples_g1 <- rownames(coldata[coldata$condition %in% group1, ])
samples_g2 <- rownames(coldata[coldata$condition %in% group2, ])

samples_g1 <- samples_g1[order(factor(coldata[samples_g1, "condition"], levels = group1))]
samples_g2 <- samples_g2[order(factor(coldata[samples_g2, "condition"], levels = group2))]

mat_g1 <- mat[, samples_g1]
mat_g2 <- mat[, samples_g2]

# Ensure identical row order
row_order_master <- rownames(mat)
mat_g1 <- mat_g1[row_order_master, ]
mat_g2 <- mat_g2[row_order_master, ]

row_split <- factor(row_split, levels = c("Ca2+ handling", "Mitochondria", "Cell stress"))

### -----------------------------
### COLUMN ANNOTATIONS
### -----------------------------

col_ha_g1 <- HeatmapAnnotation(
  Condition = coldata[samples_g1, "condition"],
  col = list(Condition = condition_colors)
)

col_ha_g2 <- HeatmapAnnotation(
  Condition = coldata[samples_g2, "condition"],
  col = list(Condition = condition_colors)
)

### -----------------------------
### SIGNIFICANCE + DIRECTION COLORS
### -----------------------------

# Function to assign colors for row names
get_row_colors <- function(mat_group, deg_file) {
  sapply(rownames(mat_group), function(gene) {
    deg <- deg_file[deg_file$GeneSymbol == gene, ]
    if (nrow(deg) == 0 || deg$significance != "significant") {
      return("black")
    } else if (deg$log2FoldChange > 0) {
      return("red")   # upregulated
    } else {
      return("blue")  # downregulated
    }
  })
}

row_colors_g1 <- get_row_colors(mat_g1, deg_g1)
row_colors_g2 <- get_row_colors(mat_g2, deg_g2)

row_gp_g1 <- grid::gpar(
  fontsize = ifelse(row_colors_g1 != "black", 9, 8),
  col = row_colors_g1,
  fontface = "italic"   
)

row_gp_g2 <- grid::gpar(
  fontsize = ifelse(row_colors_g2 != "black", 9, 8),
  col = row_colors_g2,
  fontface = "italic"   

### -----------------------------
### HEATMAPS
### -----------------------------

ht1 <- Heatmap(
  mat_g1,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  top_annotation = col_ha_g1,
  left_annotation = row_ha,
  row_split = row_split,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_names_side = "right",
  show_row_names = TRUE,
  row_names_gp = row_gp_g1,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Empty / sh190",
  column_names_rot = 80
)

ht2 <- Heatmap(
  mat_g2,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  top_annotation = col_ha_g2,
  left_annotation = row_ha,
  row_split = row_split,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = row_gp_g2,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "NaCl / MOE",
  column_names_rot = 80
)
ht1
ht2
### -----------------------------
### DRAW SEPARATELY
### -----------------------------

png("260428_output/heatmap_sh_up_down.png", width = 2000, height = 1800, res = 300)
draw(ht1)
dev.off()

png("260428_output/heatmap_MOE_up_down.png", width = 2000, height = 1800, res = 300)
draw(ht2)
dev.off()

# ============================================================
# 20260811_revision_geneset_pipeline.R
#
# CLEAN REBUILD — supersedes 20260810_revision_geneset_heatmaps.R.
# Reviewer revision analysis for curated gene sets (Ca2+ response,
# muscle injury/regeneration, mTOR/Akt once ready) across the
# Stim1 R304W treatment groups (sh190 and MOE arms).
#
# Structure:
#   0. Libraries + paths
#   1. Load your previous DESeq2 workspace (skips re-running DESeq2),
#      or rebuild from raw counts if no saved workspace is found
#   2. Fresh, empty output folders for THIS run only — your existing
#      260810_output/ (gene lists + old RData) is left untouched
#   3. Shared helper functions, defined ONCE
#   4. Gene sets, read from Excel and APPENDED to one list (not
#      overwritten — this was the bug that ate your calcium list
#      when the muscle list was added)
#   5. One-call pipeline per gene set: heatmaps + rescue table +
#      rescue-status heatmap
#   6. Save a fresh workspace snapshot OUTSIDE the output folder
# ============================================================

# Absolute path — safe to re-run any number of times in the same session.
# (Relative "../..." paths drift further each time you re-run them.)
getwd ()
setwd("~/Desktop/IGBMC/ASO_ShRNA/STIM1therORAI1_LPG_SPS26/REVISION_1")

# ---- 0. Libraries ----
library(readxl)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)

# DESeq2 pulls in Bioconductor packages (S4Vectors, IRanges, BiocGenerics...)
# that load AFTER dplyr and silently mask several dplyr verbs — count,
# filter, rename, etc. This caused "Argument 'x' is not a vector: list"
# from count() further down. Re-point these names at the dplyr versions
# explicitly so it can't happen again regardless of package load order.
select   <- dplyr::select
filter   <- dplyr::filter
rename   <- dplyr::rename
count    <- dplyr::count
mutate   <- dplyr::mutate
group_by <- dplyr::group_by
arrange  <- dplyr::arrange
pull     <- dplyr::pull

# Gene-list input files — read from where they already are. NOT deleted.
calcium_xlsx <- "260810_output/260810_calciumresponse.xlsx"
muscle_xlsx  <- "260810_output/Muscle_tissue_regeneration/260810_muscletissueregeneration.xlsx"
mtorc_xlsx = "260811_output/mTORC_pathway/260811_mTORCgenes.xlsx"

# ALL new results go here — a clean folder, independent of 260810_output.
out_dir     <- "260811_output"
calcium_out <- file.path(out_dir, "Calcium_response")
muscle_out  <- file.path(out_dir, "Muscle_tissue_regeneration")
mtorc_out = file.path(out_dir, "mTORC_pathway")

# Workspace snapshots live here, never wiped by the "fresh folders" step.
workspace_dir <- "workspace_backups"
dir.create(workspace_dir, showWarnings = FALSE)

# ---- 1. Reuse previously computed DESeq2 results if available ----
# Skips re-running DESeq() (slow). Set to NULL to force a rebuild.
previous_workspace <- "260810_output/20260810_rescue_heatmap_workspace.RData"

if (!is.null(previous_workspace) && file.exists(previous_workspace)) {
  load(previous_workspace)
  message("Loaded dds / stabilized_counts / res_disease_* / res_treated_* from ", previous_workspace)
} else {
  message("No saved workspace found at that path — rebuilding from raw counts.")
  
  countData <- read.table("../mergedReadCounts.csv", header = TRUE, sep = ",",
                          check.names = FALSE, row.names = 1)
  coldata <- as.data.frame(read_excel("../metadata.xlsx"))
  rownames(coldata) <- coldata$sample
  
  remove_sample <- "KZLX9"
  countData <- countData[, colnames(countData) != remove_sample]
  coldata   <- coldata[rownames(coldata) != remove_sample, ]
  coldata   <- coldata[colnames(countData), , drop = FALSE]
  coldata$condition <- factor(trimws(coldata$condition))
  
  dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = coldata,
                                design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  vsd <- vst(dds)
  stabilized_counts <- assay(vsd)
  
  run_deseq <- function(dds, ref, contrast) {
    dds$condition <- relevel(dds$condition, ref = ref)
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds, contrast = contrast))
    na.omit(res[, c("log2FoldChange", "padj")])
  }
  
  res_disease_e <- run_deseq(dds, "WT_empty", c("condition", "Stim1R304W/+_empty", "WT_empty"))
  res_disease_n <- run_deseq(dds, "WT_Nacl",  c("condition", "Stim1R304W/+_Nacl",  "WT_Nacl"))
  res_treated_e <- run_deseq(dds, "Stim1R304W/+_empty", c("condition", "Stim1R304W/+_sh190", "Stim1R304W/+_empty"))
  res_treated_n <- run_deseq(dds, "Stim1R304W/+_Nacl",  c("condition", "Stim1R304W/+_MOE",    "Stim1R304W/+_Nacl"))
}

condition_colors <- c(
  "WT_empty"            = "black",
  "WT_Nacl"              = "darkgrey",
  "Stim1R304W/+_sh190"   = "#FF7F00",
  "Stim1R304W/+_MOE"     = "yellow",
  "Stim1R304W/+_empty"   = "red",
  "Stim1R304W/+_Nacl"    = "brown"
)

# ---- 2. Fresh, empty output folders for THIS run ----
# Only wipes 260811_output/ (created below) — never touches 260810_output/,
# your gene-list Excel files, or your raw data.
unlink(out_dir, recursive = TRUE)
dir.create(calcium_out, recursive = TRUE, showWarnings = FALSE)
dir.create(muscle_out,  recursive = TRUE, showWarnings = FALSE)

# ---- 3. Shared helper functions (defined ONCE, reused per gene set) ----
get_geneset_ensembl <- function(symbols) {
  ens <- mapIds(org.Mm.eg.db, keys = symbols, column = "ENSEMBL",
                keytype = "SYMBOL", multiVals = "first")
  ens <- ens[!is.na(ens)]
  data.frame(Symbol = names(ens), Ensembl = unname(ens))
}

# ---- 3. Shared helper functions (defined ONCE, reused per gene set) ----
# Maps gene symbols to mouse Ensembl IDs, robust to mixed casing.
# org.Mm.eg.db's SYMBOL keytype is case-sensitive and expects mouse-style
# casing (e.g. "Abcc2"), not human-style ALL-CAPS (e.g. "ABCC2") — a
# gene set copy-pasted from a human source (like MSigDB) will silently
# fail to map for every ALL-CAPS symbol unless corrected first.
# Strategy: try each symbol as-is, then retry unmatched ones converted to
# mouse-style case (first letter upper, rest lower) against SYMBOL, then
# against ALIAS as a last resort. Anything still unmapped is reported
# rather than silently dropped.
get_geneset_ensembl <- function(symbols) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[symbols != ""]
  
  to_mouse_case <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  
  final <- setNames(rep(NA_character_, length(symbols)), symbols)
  
  # Pass 1: symbols exactly as given
  ens1 <- suppressMessages(mapIds(org.Mm.eg.db, keys = symbols, column = "ENSEMBL",
                                  keytype = "SYMBOL", multiVals = "first"))
  final[names(ens1)] <- ifelse(is.na(final[names(ens1)]), ens1, final[names(ens1)])
  
  # Pass 2: unmatched, retried in mouse-style case against SYMBOL
  unmatched <- names(final)[is.na(final)]
  if (length(unmatched) > 0) {
    ens2 <- suppressMessages(mapIds(org.Mm.eg.db, keys = to_mouse_case(unmatched), column = "ENSEMBL",
                                    keytype = "SYMBOL", multiVals = "first"))
    names(ens2) <- unmatched  # keep the ORIGINAL symbol as the key
    final[unmatched] <- ifelse(is.na(final[unmatched]), ens2, final[unmatched])
  }
  
  # Pass 3: still unmatched, retried in mouse-style case against ALIAS
  unmatched <- names(final)[is.na(final)]
  if (length(unmatched) > 0) {
    ens3 <- suppressMessages(mapIds(org.Mm.eg.db, keys = to_mouse_case(unmatched), column = "ENSEMBL",
                                    keytype = "ALIAS", multiVals = "first"))
    names(ens3) <- unmatched
    final[unmatched] <- ifelse(is.na(final[unmatched]), ens3, final[unmatched])
  }
  
  still_missing <- names(final)[is.na(final)]
  if (length(still_missing) > 0) {
    message(length(still_missing), " of ", length(symbols),
            " symbols could not be mapped to a mouse Ensembl ID (SYMBOL or ALIAS, either casing): ",
            paste(head(still_missing, 20), collapse = ", "),
            if (length(still_missing) > 20) ", ..." else "")
  }
  
  final <- final[!is.na(final)]
  
  # Look up the CANONICAL mouse symbol for each Ensembl ID, so heatmap
  # labels always show proper mouse-style casing (e.g. "Hmox1") regardless
  # of how the gene was written in your input list (HMOX1, Hmox1, hmox1...).
  canonical_symbol <- suppressMessages(mapIds(org.Mm.eg.db, keys = unname(final),
                                              column = "SYMBOL", keytype = "ENSEMBL",
                                              multiVals = "first"))
  
  data.frame(Symbol = unname(canonical_symbol), Ensembl = unname(final), row.names = NULL)
}

# 3-panel dysregulation heatmap: combined / sh-only / MOE-only.
# Selects the top_n genes by disease-effect magnitude (not rescue status —
# see build_rescue_heatmap() below for the rescue-driven selection).
build_geneset_heatmaps <- function(geneset_name, symbols,
                                   stabilized_counts, coldata,
                                   res_disease_e, res_disease_n,
                                   condition_colors,
                                   top_n = 30,
                                   padj_cutoff = 0.05,
                                   row_axis_label,
                                   out_dir) {
  
  gs <- get_geneset_ensembl(symbols)
  gs <- gs[gs$Ensembl %in% rownames(stabilized_counts), ]
  
  stats_e <- res_disease_e[gs$Ensembl, ]
  colnames(stats_e) <- paste0(colnames(stats_e), "_e")
  stats_n <- res_disease_n[gs$Ensembl, ]
  colnames(stats_n) <- paste0(colnames(stats_n), "_n")
  gs <- cbind(gs, stats_e, stats_n)
  
  gs$sig_e <- !is.na(gs$padj_e) & gs$padj_e < padj_cutoff
  gs$sig_n <- !is.na(gs$padj_n) & gs$padj_n < padj_cutoff
  gs$max_abs_lfc <- pmax(abs(gs$log2FoldChange_e), abs(gs$log2FoldChange_n), na.rm = TRUE)
  
  gs_sig <- gs %>% filter(sig_e | sig_n) %>% arrange(desc(max_abs_lfc))
  gs_top <- head(gs_sig, top_n)
  
  write.csv(gs, file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_full_stats.csv")),
            row.names = FALSE)
  write.csv(gs_top, file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_top", top_n, "_stats.csv")),
            row.names = FALSE)
  
  if (nrow(gs_top) < 2) {
    message("Fewer than 2 significant genes for ", geneset_name, " — skipping heatmaps.")
    return(invisible(NULL))
  }
  
  mat <- stabilized_counts[gs_top$Ensembl, , drop = FALSE]
  mat <- t(scale(t(log2(mat + 1))))
  rownames(mat) <- gs_top$Symbol
  
  desired_order <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190",
                     "WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  coldata$condition <- factor(coldata$condition, levels = desired_order)
  
  group_sh  <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190")
  group_moe <- c("WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  
  samples_all <- rownames(coldata[order(coldata$condition), ])
  samples_sh  <- rownames(coldata[coldata$condition %in% group_sh, ])
  samples_sh  <- samples_sh[order(factor(coldata[samples_sh, "condition"], levels = group_sh))]
  samples_moe <- rownames(coldata[coldata$condition %in% group_moe, ])
  samples_moe <- samples_moe[order(factor(coldata[samples_moe, "condition"], levels = group_moe))]
  
  col_combined <- ifelse(gs_top$sig_e & gs_top$sig_n, "darkgreen", "black")
  col_sh  <- ifelse(!gs_top$sig_e, "black", ifelse(gs_top$log2FoldChange_e > 0, "red3", "blue3"))
  col_moe <- ifelse(!gs_top$sig_n, "black", ifelse(gs_top$log2FoldChange_n > 0, "red3", "blue3"))
  
  make_ht <- function(mat_full, samples, row_colors, title) {
    col_ha <- HeatmapAnnotation(
      Condition = coldata[samples, "condition"],
      col = list(Condition = condition_colors),
      annotation_name_gp = gpar(fontsize = 10)
    )
    Heatmap(
      mat_full[, samples],
      name = "z-score",
      col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
      top_annotation = col_ha,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 12, col = row_colors, fontface = "italic"),
      column_names_gp = gpar(fontsize = 11),
      column_names_rot = 80,
      column_title = title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      row_title = row_axis_label,
      row_title_gp = gpar(fontsize = 13, fontface = "bold")
    )
  }
  
  ht_combined <- make_ht(mat, samples_all, col_combined, "Empty/sh190 + Nacl/MOE (combined)")
  ht_sh       <- make_ht(mat, samples_sh,  col_sh,       "Empty / sh190")
  ht_moe      <- make_ht(mat, samples_moe, col_moe,      "Nacl / MOE")
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_heatmap_combined.pdf")), width = 9, height = 11)
  draw(ht_combined); dev.off()
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_heatmap_sh.pdf")), width = 7, height = 11)
  draw(ht_sh); dev.off()
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_heatmap_MOE.pdf")), width = 7, height = 11)
  draw(ht_moe); dev.off()
  
  invisible(list(table = gs_top, combined = ht_combined, sh = ht_sh, moe = ht_moe))
}

# Rescue table: did the treatment reverse the disease-dysregulated genes?
# rescue_metric: 0% = no effect, 100% = back to WT, <0% = worsened
# (same direction as disease, stronger), >120% = over-corrected (flipped
# past WT in the opposite direction).
build_geneset_rescue_table <- function(geneset_name, symbols,
                                       res_disease_e, res_treated_e,
                                       res_disease_n, res_treated_n,
                                       log2FC_cutoff = 0.5,
                                       padj_cutoff = 0.05,
                                       out_dir) {
  
  gs <- get_geneset_ensembl(symbols)
  
  compute_rescue <- function(res_disease, res_treated, gs, arm_label) {
    diseased <- res_disease[!is.na(res_disease$padj) &
                              res_disease$padj < padj_cutoff &
                              abs(res_disease$log2FoldChange) >= log2FC_cutoff, ]
    common <- intersect(gs$Ensembl, intersect(rownames(diseased), rownames(res_treated)))
    if (length(common) == 0) return(NULL)
    
    gs_sub <- gs[match(common, gs$Ensembl), ]
    d <- diseased[common, ]
    t <- res_treated[common, ]
    
    log2FC_disease <- d$log2FoldChange
    log2FC_treated <- t$log2FoldChange
    log2FC_treated_vs_WT <- log2FC_disease + log2FC_treated
    rescue_metric <- 100 * (log2FC_disease - log2FC_treated_vs_WT) / log2FC_disease
    
    out <- data.frame(
      Symbol = gs_sub$Symbol,
      Ensembl = gs_sub$Ensembl,
      log2FC_disease = log2FC_disease,
      log2FC_treated = log2FC_treated,
      log2FC_treated_vs_WT = log2FC_treated_vs_WT,
      padj_disease = d$padj,
      padj_treated = t$padj,
      rescue_metric = rescue_metric,
      Arm = arm_label
    )
    out$rescue_status <- cut(
      out$rescue_metric,
      breaks = c(-Inf, 0, 30, 80, 120, Inf),
      labels = c("Worsened", "Not rescued", "Partially rescued", "Rescued", "Over-corrected")
    )
    out
  }
  
  rescue_sh  <- compute_rescue(res_disease_e, res_treated_e, gs, "sh190")
  rescue_moe <- compute_rescue(res_disease_n, res_treated_n, gs, "MOE")
  rescue_all <- bind_rows(rescue_sh, rescue_moe)
  
  if (nrow(rescue_all) == 0) {
    message("No disease-dysregulated genes to score for ", geneset_name)
    return(invisible(NULL))
  }
  
  write.csv(rescue_all,
            file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_rescue_table.csv")),
            row.names = FALSE)
  
  summary_df <- rescue_all %>%
    filter(!is.na(rescue_status)) %>%
    count(Arm, rescue_status) %>%
    group_by(Arm) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()
  
  p <- ggplot(summary_df, aes(x = Arm, y = percent, fill = rescue_status)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c(
      "Worsened"           = "black",
      "Not rescued"        = "lightgrey",
      "Partially rescued"  = "orange",
      "Rescued"            = "forestgreen",
      "Over-corrected"     = "purple"
    )) +
    labs(title = paste0(gsub("_", " ", geneset_name), " — treatment rescue status"),
         x = NULL, y = "% of disease-dysregulated genes", fill = "Rescue status") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(size = 12, face = "bold"))
  
  ggsave(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_rescue_summary_barplot.pdf")),
         plot = p, width = 6, height = 5)
  
  invisible(list(table = rescue_all, summary = summary_df, plot = p))
}

# Rescue-STATUS heatmap: shows only genes rescued/partially rescued by
# sh190, MOE, or both — the gene selection you actually asked for
# ("show me the ones that were rescued").
build_rescue_heatmap <- function(geneset_name,
                                 rescue_csv,
                                 stabilized_counts,
                                 coldata,
                                 condition_colors,
                                 out_dir,
                                 rescue_cutoff = 30,
                                 row_axis_label = "Rescued genes") {
  
  rescue <- read.csv(rescue_csv, stringsAsFactors = FALSE)
  rescue$rescue_metric <- as.numeric(rescue$rescue_metric)
  
  rescue_selected <- rescue %>%
    filter(!is.na(rescue_metric)) %>%
    group_by(Symbol, Ensembl) %>%
    summarise(
      rescue_sh  = ifelse(any(Arm == "sh190"), max(rescue_metric[Arm == "sh190"], na.rm = TRUE), NA_real_),
      rescue_moe = ifelse(any(Arm == "MOE"),   max(rescue_metric[Arm == "MOE"],   na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      rescued_sh   = !is.na(rescue_sh)  & rescue_sh  >= rescue_cutoff,
      rescued_moe  = !is.na(rescue_moe) & rescue_moe >= rescue_cutoff,
      rescue_group = case_when(
        rescued_sh & rescued_moe  ~ "Both treatments",
        rescued_sh & !rescued_moe ~ "sh190 only",
        !rescued_sh & rescued_moe ~ "MOE only",
        TRUE                      ~ "Not rescued"
      )
    ) %>%
    filter(rescue_group != "Not rescued") %>%
    mutate(
      max_rescue  = pmax(rescue_sh, rescue_moe, na.rm = TRUE),
      group_order = case_when(
        rescue_group == "Both treatments" ~ 1,
        rescue_group == "sh190 only"      ~ 2,
        rescue_group == "MOE only"        ~ 3
      )
    ) %>%
    arrange(group_order, desc(max_rescue))
  
  message("\nGenes included in rescue heatmap: ", nrow(rescue_selected))
  message("Genes rescued by BOTH treatments: ",
          paste(rescue_selected %>% filter(rescue_group == "Both treatments") %>% pull(Symbol), collapse = ", "))
  message("Genes rescued by sh190 ONLY: ",
          paste(rescue_selected %>% filter(rescue_group == "sh190 only") %>% pull(Symbol), collapse = ", "))
  message("Genes rescued by MOE ONLY: ",
          paste(rescue_selected %>% filter(rescue_group == "MOE only") %>% pull(Symbol), collapse = ", "))
  
  write.csv(rescue_selected,
            file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_rescued_genes_for_heatmap.csv")),
            row.names = FALSE)
  
  rescue_selected <- rescue_selected %>% filter(Ensembl %in% rownames(stabilized_counts))
  if (nrow(rescue_selected) < 2) {
    message("Fewer than 2 rescued genes found in expression matrix.")
    return(invisible(NULL))
  }
  
  mat <- stabilized_counts[rescue_selected$Ensembl, , drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(scale(t(mat)))
  rownames(mat) <- rescue_selected$Symbol
  
  desired_order <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190",
                     "WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  coldata$condition <- factor(coldata$condition, levels = desired_order)
  
  group_sh  <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190")
  group_moe <- c("WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  
  samples_all <- rownames(coldata[order(coldata$condition), , drop = FALSE])
  samples_sh  <- rownames(coldata[coldata$condition %in% group_sh, , drop = FALSE])
  samples_sh  <- samples_sh[order(factor(coldata[samples_sh, "condition"], levels = group_sh))]
  samples_moe <- rownames(coldata[coldata$condition %in% group_moe, , drop = FALSE])
  samples_moe <- samples_moe[order(factor(coldata[samples_moe, "condition"], levels = group_moe))]
  
  row_ha <- rowAnnotation(
    Rescue = rescue_selected$rescue_group,
    col = list(Rescue = c(
      "Both treatments" = "darkgreen",
      "sh190 only"      = "royalblue3",
      "MOE only"        = "purple3"
    )),
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
    show_legend = TRUE
  )
  
  make_rescue_ht <- function(samples, title) {
    col_ha <- HeatmapAnnotation(
      Condition = coldata[samples, "condition"],
      col = list(Condition = condition_colors),
      annotation_name_gp = gpar(fontsize = 10)
    )
    Heatmap(
      mat[, samples, drop = FALSE],
      name = "Z-score",
      col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
      top_annotation = col_ha,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_dend = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 11, fontface = "italic"),
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 80,
      column_title = title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      row_title = row_axis_label,
      row_title_gp = gpar(fontsize = 13, fontface = "bold")
    )
  }
  
  ht_combined <- row_ha + make_rescue_ht(samples_all, "Genes rescued by either or both treatments")
  ht_sh       <- row_ha + make_rescue_ht(samples_sh,  "sh190 rescue")
  ht_moe      <- row_ha + make_rescue_ht(samples_moe, "MOE rescue")
  
  ht_height <- max(8, 0.25 * nrow(rescue_selected) + 4)
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_RESCUED_heatmap_combined.pdf")), width = 10, height = ht_height)
  draw(ht_combined, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_RESCUED_heatmap_sh190.pdf")), width = 8, height = ht_height)
  draw(ht_sh, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_RESCUED_heatmap_MOE.pdf")), width = 8, height = ht_height)
  draw(ht_moe, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
  
  invisible(list(genes = rescue_selected, combined = ht_combined, sh190 = ht_sh, moe = ht_moe))
}

# One-call wrapper: ties heatmaps + rescue table + rescue-status heatmap
# together for a single gene set. This is the piece that stops the
# copy-paste drift — add a new gene set by adding ONE call, not by
# duplicating three function definitions.
run_full_geneset_pipeline <- function(geneset_name, symbols, out_dir, row_axis_label,
                                      stabilized_counts, coldata,
                                      res_disease_e, res_disease_n,
                                      res_treated_e, res_treated_n,
                                      condition_colors,
                                      top_n = 30, rescue_cutoff = 30) {
  message("\n==== ", geneset_name, " ====")
  
  heatmaps <- build_geneset_heatmaps(
    geneset_name = geneset_name, symbols = symbols,
    stabilized_counts = stabilized_counts, coldata = coldata,
    res_disease_e = res_disease_e, res_disease_n = res_disease_n,
    condition_colors = condition_colors, top_n = top_n,
    row_axis_label = row_axis_label, out_dir = out_dir
  )
  
  rescue <- build_geneset_rescue_table(
    geneset_name = geneset_name, symbols = symbols,
    res_disease_e = res_disease_e, res_treated_e = res_treated_e,
    res_disease_n = res_disease_n, res_treated_n = res_treated_n,
    out_dir = out_dir
  )
  
  rescue_heatmap <- NULL
  if (!is.null(rescue)) {
    rescue_csv <- file.path(out_dir, paste0(Sys.Date(), "_", geneset_name, "_rescue_table.csv"))
    rescue_heatmap <- build_rescue_heatmap(
      geneset_name = geneset_name, rescue_csv = rescue_csv,
      stabilized_counts = stabilized_counts, coldata = coldata,
      condition_colors = condition_colors, out_dir = out_dir,
      rescue_cutoff = rescue_cutoff,
      row_axis_label = paste("Rescued", tolower(row_axis_label))
    )
  }
  
  list(heatmaps = heatmaps, rescue = rescue, rescue_heatmap = rescue_heatmap)
}


# ---- 4. Gene sets — APPENDED to one list, never overwritten ----
calcium_genes_raw <- read_excel(calcium_xlsx)
calcium_genes <- unique(trimws(as.character(calcium_genes_raw[[1]])))
calcium_genes <- calcium_genes[!is.na(calcium_genes) & calcium_genes != ""]

muscle_genes_raw <- read_excel(muscle_xlsx)
muscle_genes <- unique(trimws(as.character(muscle_genes_raw[[1]])))
muscle_genes <- muscle_genes[!is.na(muscle_genes) & muscle_genes != ""]

mtorc_genes_raw = read_excel(mtorc_xlsx)
mtorc_genes <- unique(trimws(as.character(mtorc_genes_raw[[1]])))
mtorc_genes <- mtorc_genes[!is.na(mtorc_genes) & mtorc_genes != ""]

geneset_list <- list()
geneset_list$Calcium_response           <- calcium_genes
geneset_list$Muscle_injury_regeneration <- muscle_genes
geneset_list$mTOR_Akt <-  mtorc_genes

# ---- 5. Run the pipeline per gene set ----
results_calcium <- run_full_geneset_pipeline(
  geneset_name      = "Calcium_response",
  symbols           = geneset_list$Calcium_response,
  out_dir           = calcium_out,
  row_axis_label    = "Genes regulated by Ca2+",
  stabilized_counts = stabilized_counts, coldata = coldata,
  res_disease_e = res_disease_e, res_disease_n = res_disease_n,
  res_treated_e = res_treated_e, res_treated_n = res_treated_n,
  condition_colors  = condition_colors,
  top_n = 30
)

results_muscle_regen <- run_full_geneset_pipeline(
  geneset_name      = "Muscle_injury_regeneration",
  symbols           = geneset_list$Muscle_injury_regeneration,
  out_dir           = muscle_out,
  row_axis_label    = "Genes involved in skeletal muscle regeneration",
  stabilized_counts = stabilized_counts, coldata = coldata,
  res_disease_e = res_disease_e, res_disease_n = res_disease_n,
  res_treated_e = res_treated_e, res_treated_n = res_treated_n,
  condition_colors  = condition_colors,
  top_n = 50
)

results_mtorc <- run_full_geneset_pipeline(
  geneset_name      = "mTORC pathway",
  symbols           = geneset_list$mTOR_Akt,
  out_dir           = mtorc_out,
  row_axis_label    = "Genes involved in mTORC pathway",
  stabilized_counts = stabilized_counts, coldata = coldata,
  res_disease_e = res_disease_e, res_disease_n = res_disease_n,
  res_treated_e = res_treated_e, res_treated_n = res_treated_n,
  condition_colors  = condition_colors,
  top_n = 30
)


# ---- 6. Save a fresh workspace snapshot OUTSIDE out_dir ----
save.image(file = file.path(workspace_dir, paste0(Sys.Date(), "_geneset_pipeline_workspace.RData")))
message("Workspace saved to ", workspace_dir)

# ============================================================
# ---- 7. Per-gene bar plots for rescued genes, with ANOVA + Tukey
# ============================================================
# Same style as your original genes_of_interest plot_gene() function:
# jittered points (triangles) + mean crossbar + SEM error bars, colored
# by condition, with Tukey-adjusted pairwise significance brackets for
# WT vs KO vs treatment within each arm (sh190 arm and MOE arm separately).
# Runs automatically over whichever genes came out of build_rescue_heatmap()
# — no gene list needs to be typed in by hand.
# NOTE: uses base ggplot2 only (no ggpubr) to sidestep the dplyr>=1.2.0
# version requirement that broke ggpubr's namespace load.

# Draws significance brackets manually (segment + tip + label), replacing
# ggpubr::stat_pvalue_manual so this doesn't depend on ggpubr at all.
add_significance_brackets <- function(p, tuk_df, condition_order, tip_length = 0.02) {
  if (nrow(tuk_df) == 0) return(p)
  
  tuk_df <- tuk_df %>%
    mutate(
      x1 = as.numeric(factor(as.character(group1), levels = condition_order)),
      x2 = as.numeric(factor(as.character(group2), levels = condition_order))
    )
  
  for (i in seq_len(nrow(tuk_df))) {
    row <- tuk_df[i, ]
    y <- row$y.position
    tip <- tip_length * y
    p <- p +
      annotate("segment", x = row$x1, xend = row$x2, y = y, yend = y) +
      annotate("segment", x = row$x1, xend = row$x1, y = y, yend = y - tip) +
      annotate("segment", x = row$x2, xend = row$x2, y = y, yend = y - tip) +
      annotate("text", x = (row$x1 + row$x2) / 2, y = y, label = row$p.adj.signif,
               vjust = -0.3, size = 5)
  }
  p
}

plot_rescued_gene_barplots <- function(rescue_selected, dds, coldata,
                                       condition_colors, out_dir) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # dds coming out of run_deseq()/vst() never had estimateSizeFactors()
  # applied to the object itself (only to internal copies) — set it here
  # so counts(dds, normalized=TRUE) has something to normalize by.
  if (is.null(DESeq2::sizeFactors(dds))) {
    dds <- DESeq2::estimateSizeFactors(dds)
  }
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  condition_order <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190",
                       "WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  
  coldata2 <- coldata %>%
    rownames_to_column("Sample") %>%
    rename(Condition = condition) %>%
    mutate(Condition = factor(Condition, levels = condition_order))
  
  # Only WT-vs-KO-vs-treatment comparisons WITHIN each arm — matches your
  # original design (not every pairwise combination across all 6 groups)
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
  
  plots <- list()
  
  for (i in seq_len(nrow(rescue_selected))) {
    gene_symbol <- rescue_selected$Symbol[i]
    ensembl_id  <- rescue_selected$Ensembl[i]
    
    if (!ensembl_id %in% rownames(norm_counts)) {
      message("Skipping ", gene_symbol, " — not found in normalized counts.")
      next
    }
    
    expr <- norm_counts[ensembl_id, ]
    df <- data.frame(Sample = names(expr), Expression = as.numeric(expr)) %>%
      left_join(coldata2, by = "Sample") %>%
      filter(!is.na(Condition))
    
    summary_df <- df %>%
      group_by(Condition) %>%
      summarise(mean = mean(Expression), sem = sd(Expression) / sqrt(n()), .groups = "drop")
    
    fit <- aov(Expression ~ Condition, data = df)
    tuk <- TukeyHSD(fit)
    tuk_df <- as.data.frame(tuk$Condition) %>%
      rownames_to_column("comparison") %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(pair = paste(pmin(group1, group2), pmax(group1, group2), sep = "_")) %>%
      rename(p.adj = `p adj`) %>%
      semi_join(my_pairs, by = "pair")
    
    if (nrow(tuk_df) > 0) {
      tuk_df <- tuk_df %>%
        mutate(
          p.adj.signif = case_when(
            p.adj < 0.0001 ~ "****",
            p.adj < 0.001  ~ "***",
            p.adj < 0.01   ~ "**",
            p.adj < 0.05   ~ "*",
            TRUE           ~ "ns"
          ),
          group1 = factor(group1, levels = condition_order),
          group2 = factor(group2, levels = condition_order),
          y.position = max(df$Expression) * seq(1.15, 1.45, length.out = n())
        )
    }
    
    p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_jitter(width = 0.12, size = 3, shape = 24) +
      geom_errorbar(data = summary_df,
                    aes(x = Condition, y = mean, ymin = mean - sem, ymax = mean + sem),
                    inherit.aes = FALSE, width = 0.25, linewidth = 0.8) +
      geom_crossbar(data = summary_df,
                    aes(x = Condition, y = mean, ymin = mean, ymax = mean),
                    inherit.aes = FALSE, width = 0.6, linewidth = 0.4, color = "black") +
      scale_fill_manual(values = condition_colors) +
      theme_classic(base_size = 14) +
      labs(title = gene_symbol, y = "Normalized expression", x = NULL) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    if (nrow(tuk_df) > 0) {
      p <- add_significance_brackets(p, tuk_df, condition_order)
    }
    
    ggsave(file.path(out_dir, paste0(Sys.Date(), "_", gene_symbol, "_barplot.png")),
           plot = p, width = 6, height = 5, dpi = 300)
    
    plots[[gene_symbol]] <- p
  }
  
  message("Saved ", length(plots), " gene bar plots to ", out_dir)
  invisible(plots)
}

# ---- Run for calcium rescued genes ----
if (!is.null(results_calcium$rescue_heatmap)) {
  barplots_calcium <- plot_rescued_gene_barplots(
    rescue_selected  = results_calcium$rescue_heatmap$genes,
    dds              = dds,
    coldata          = coldata,
    condition_colors = condition_colors,
    out_dir          = file.path(calcium_out, "gene_barplots")
  )
}

# ---- Run for muscle regeneration rescued genes ----
if (!is.null(results_muscle_regen$rescue_heatmap)) {
  barplots_muscle_regen <- plot_rescued_gene_barplots(
    rescue_selected  = results_muscle_regen$rescue_heatmap$genes,
    dds              = dds,
    coldata          = coldata,
    condition_colors = condition_colors,
    out_dir          = file.path(muscle_out, "gene_barplots")
  )
}

if (!is.null(results_mtorc$rescue_heatmap)) {
  barplots_mtorc <- plot_rescued_gene_barplots(
    rescue_selected  = results_mtorc$rescue_heatmap$genes,
    dds              = dds,
    coldata          = coldata,
    condition_colors = condition_colors,
    out_dir          = file.path(mtorc_out, "gene_barplots")
  )
}



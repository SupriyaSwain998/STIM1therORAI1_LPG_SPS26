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

# ============================================================
# ---- 8. Pathway-level GSEA: is each gene set coordinately
#          perturbed in disease, and does that signature resolve
#          after treatment?
# ============================================================
# Three contrasts per arm, all against the SAME gene set (no top_n
# filtering — this uses the entire ranked transcriptome, which is what
# GSEA is designed for):
#   1. Disease (KO vs WT)      — is the pathway perturbed as a set?
#   2. Treatment vs KO         — does treatment push the pathway in the
#                                 REVERSE direction of the disease effect?
#   3. Treatment vs WT         — after treatment, does the pathway look
#                                 statistically like WT again? (the actual
#                                 "cured" test — NES near 0 and non-sig
#                                 here is the strongest rescue evidence)
#
# Interpretation:
#   - Disease step: large |NES|, padj < 0.05  -> pathway is dysregulated
#   - Treatment vs WT step: NES shrinks toward 0 and/or padj > 0.05
#     compared to the Disease step -> pathway signature resolving
#   - Treatment vs WT step: NES stays large/significant, same sign as
#     Disease -> pathway NOT rescued
#   - Treatment vs WT step: NES large/significant, OPPOSITE sign ->
#     over-corrected at the pathway level

if (!requireNamespace("fgsea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("fgsea")
}
library(fgsea)

# Treatment vs WT — direct DESeq2 contrasts (not the log2FC-summing
# approximation used in the per-gene rescue table). run_deseq() was
# saved in your workspace RData, so this works without rebuilding dds.
res_treated_vs_WT_e <- run_deseq(dds, "WT_empty", c("condition", "Stim1R304W/+_sh190", "WT_empty"))
res_treated_vs_WT_n <- run_deseq(dds, "WT_Nacl",  c("condition", "Stim1R304W/+_MOE",    "WT_Nacl"))

run_pathway_gsea <- function(geneset_list, contrasts, out_dir,
                             minSize = 5, maxSize = 500) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Convert each gene set (Symbols) to Ensembl once, reusing the same
  # case-robust mapper used everywhere else in this script
  pathways_ensembl <- lapply(geneset_list, function(symbols) {
    get_geneset_ensembl(symbols)$Ensembl
  })
  
  all_results <- list()
  
  for (contrast_name in names(contrasts)) {
    res <- contrasts[[contrast_name]]
    ranks <- res$log2FoldChange
    names(ranks) <- rownames(res)
    ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)
    
    fgsea_res <- fgsea::fgsea(pathways = pathways_ensembl, stats = ranks,
                              minSize = minSize, maxSize = maxSize, eps = 0)
    fgsea_res$Contrast <- contrast_name
    all_results[[contrast_name]] <- as.data.frame(
      fgsea_res[, c("pathway", "NES", "pval", "padj", "size", "Contrast")]
    )
  }
  
  combined <- bind_rows(all_results)
  write.csv(combined, file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_results.csv")),
            row.names = FALSE)
  combined
}

# "Arm|Step" naming so both arms share identical step labels for faceting
gsea_contrasts <- list(
  "sh190 arm|Disease (KO vs WT)" = res_disease_e,
  "sh190 arm|Treatment vs KO"    = res_treated_e,
  "sh190 arm|Treatment vs WT"    = res_treated_vs_WT_e,
  "MOE arm|Disease (KO vs WT)"   = res_disease_n,
  "MOE arm|Treatment vs KO"      = res_treated_n,
  "MOE arm|Treatment vs WT"      = res_treated_vs_WT_n
)
# ============================================================
# Rebuild gene sets with mTOR symbols standardized to mouse case
# ============================================================

calcium_genes_raw <- read_excel(calcium_xlsx)
calcium_genes <- unique(trimws(as.character(calcium_genes_raw[[1]])))
calcium_genes <- calcium_genes[
  !is.na(calcium_genes) & calcium_genes != ""
]

muscle_genes_raw <- read_excel(muscle_xlsx)
muscle_genes <- unique(trimws(as.character(muscle_genes_raw[[1]])))
muscle_genes <- muscle_genes[
  !is.na(muscle_genes) & muscle_genes != ""
]

mtorc_genes_raw <- read_excel(mtorc_xlsx)
mtorc_genes <- unique(trimws(as.character(mtorc_genes_raw[[1]])))
mtorc_genes <- mtorc_genes[
  !is.na(mtorc_genes) & mtorc_genes != ""
]

# ------------------------------------------------------------
# Convert ALL-CAPS gene symbols to mouse-style symbols
#
# CAV3  -> Cav3
# CLTA  -> Clta
# MTOR  -> Mtor
# AKT1  -> Akt1
#
# Mixed-case symbols are left unchanged.
# ------------------------------------------------------------

to_mouse_case_explicit <- function(x) {
  
  is_all_caps <- x == toupper(x) &
    grepl("[A-Z]", x)
  
  x[is_all_caps] <- paste0(
    toupper(substr(x[is_all_caps], 1, 1)),
    tolower(substr(x[is_all_caps], 2, nchar(x[is_all_caps])))
  )
  
  x
}

mtorc_genes <- to_mouse_case_explicit(mtorc_genes)
mtorc_genes <- unique(mtorc_genes)

# ------------------------------------------------------------
# Build the complete gene-set list
# ------------------------------------------------------------

geneset_list <- list(
  Calcium_response           = calcium_genes,
  Muscle_injury_regeneration = muscle_genes,
  mTOR_Akt                   = mtorc_genes
)

# Check
cat("\nGene-set sizes:\n")
print(sapply(geneset_list, length))

cat("\nFirst mTOR genes after case correction:\n")
print(head(geneset_list$mTOR_Akt, 30))
pathway_gsea_results <- run_pathway_gsea(
  geneset_list = list(
    Calcium_response            = geneset_list$Calcium_response,
    Muscle_injury_regeneration  = geneset_list$Muscle_injury_regeneration,
    mTOR_Akt                    = geneset_list$mTOR_Akt
  ),
  contrasts = gsea_contrasts,
  out_dir   = out_dir
)

# ---- Summary plot: NES across Disease -> Treatment vs KO -> Treatment vs WT ----
gsea_plot_df <- pathway_gsea_results %>%
  separate(Contrast, into = c("Arm", "Step"), sep = "\\|") %>%
  mutate(
    Step = factor(Step, levels = c("Disease (KO vs WT)", "Treatment vs KO", "Treatment vs WT")),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  )

p_gsea <- ggplot(gsea_plot_df, aes(x = Step, y = NES, fill = significance)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  facet_grid(pathway ~ Arm) +
  scale_fill_manual(values = c("***" = "firebrick3", "**" = "orange",
                               "*" = "gold2", "ns" = "grey70")) +
  theme_classic(base_size = 13) +
  labs(y = "Normalized Enrichment Score (NES)", x = NULL, fill = "padj") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_NES_summary.pdf")),
       plot = p_gsea, width = 9, height = 8)

message("Pathway GSEA complete — see ", file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_NES_summary.pdf")))
EOF
echo "done"
Output



# ============================================================
# 8. STANDARD GENOME-WIDE GSEA
# ============================================================
#
# Goal:
#   Run unbiased GSEA across the ENTIRE ranked transcriptome.
#
# Comparisons:
#   1. Disease:       KO vs WT
#   2. Treatment vs KO
#   3. Treatment vs WT
#
# Arms:
#   sh190
#   MOE
#
# Gene sets:
#   Hallmark + GO Biological Process
#
# Output:
#   - Complete GSEA result table
#   - NES heatmap for all significantly enriched pathways
#   - Focused heatmap of pathways significant in disease
#     and/or altered by treatment
# ============================================================

if (!requireNamespace("fgsea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("fgsea")
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(fgsea)
library(msigdbr)


# ============================================================
# 8.1 Direct treatment-vs-WT contrasts
# ============================================================

res_treated_vs_WT_e <- run_deseq(
  dds,
  "WT_empty",
  c("condition", "Stim1R304W/+_sh190", "WT_empty")
)

res_treated_vs_WT_n <- run_deseq(
  dds,
  "WT_Nacl",
  c("condition", "Stim1R304W/+_MOE", "WT_Nacl")
)


# ============================================================
# 8.2 Build the six contrasts
# ============================================================

gsea_contrasts <- list(
  
  "sh190|Disease" =
    res_disease_e,
  
  "sh190|Treatment_vs_KO" =
    res_treated_e,
  
  "sh190|Treatment_vs_WT" =
    res_treated_vs_WT_e,
  
  "MOE|Disease" =
    res_disease_n,
  
  "MOE|Treatment_vs_KO" =
    res_treated_n,
  
  "MOE|Treatment_vs_WT" =
    res_treated_vs_WT_n
)


# ============================================================
# 8.3 Get mouse Hallmark gene sets
# ============================================================

hallmark_df <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

hallmark_sets <- split(
  hallmark_df$ensembl_gene,
  hallmark_df$gs_name
)


# ============================================================
# 8.4 Get mouse GO Biological Process gene sets
# ============================================================

gobp_df <- msigdbr(
  species = "Mus musculus",
  collection = "C5",
  subcollection = "GO:BP"
)

gobp_sets <- split(
  gobp_df$ensembl_gene,
  gobp_df$gs_name
)


# ============================================================
# 8.5 Combine gene sets
# ============================================================

gene_sets <- c(
  hallmark_sets,
  gobp_sets
)

# Remove duplicate genes within pathways
gene_sets <- lapply(
  gene_sets,
  unique
)

# Keep pathways of reasonable size
gene_sets <- gene_sets[
  lengths(gene_sets) >= 10 &
    lengths(gene_sets) <= 500
]

message(
  "Number of pathways used for GSEA: ",
  length(gene_sets)
)


# ============================================================
# 8.6 Run GSEA for every contrast
# ============================================================

all_gsea <- list()

for (contrast_name in names(gsea_contrasts)) {
  
  message("\nRunning GSEA: ", contrast_name)
  
  res <- gsea_contrasts[[contrast_name]]
  
  # ----------------------------------------------------------
  # Rank entire transcriptome by log2FC
  # ----------------------------------------------------------
  
  ranks <- res$log2FoldChange
  
  names(ranks) <- rownames(res)
  
  ranks <- ranks[
    !is.na(ranks) &
      is.finite(ranks)
  ]
  
  # Remove duplicated Ensembl IDs
  ranks <- ranks[
    !duplicated(names(ranks))
  ]
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  # ----------------------------------------------------------
  # Run fgsea
  # ----------------------------------------------------------
  
  fg <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )
  
  fg <- as.data.frame(fg)
  
  fg$Contrast <- contrast_name
  
  all_gsea[[contrast_name]] <- fg
  
}


# ============================================================
# 8.7 Combine all results
# ============================================================

gsea_results <- bind_rows(all_gsea)


# ============================================================
# 8.8 Clean pathway names
# ============================================================

gsea_results$Pathway <- gsea_results$pathway

gsea_results$Pathway_clean <- gsea_results$Pathway

gsea_results$Pathway_clean <- gsub(
  "^HALLMARK_",
  "",
  gsea_results$Pathway_clean
)

gsea_results$Pathway_clean <- gsub(
  "^GO_BP_",
  "",
  gsea_results$Pathway_clean
)

gsea_results$Pathway_clean <- gsub(
  "_",
  " ",
  gsea_results$Pathway_clean
)


# ============================================================
# 8.9 Save complete GSEA results
# ============================================================

# Remove list-columns such as leadingEdge before writing CSV
gsea_results_export <- gsea_results %>%
  dplyr::select(
    -dplyr::where(is.list)
  )

write.csv(
  gsea_results_export,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_all_pathways.csv"
    )
  ),
  row.names = FALSE
)

message(
  "GSEA results saved: ",
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_all_pathways.csv"
    )
  )
)


# ============================================================
# 8.10 Create pathway x contrast NES matrix
# ============================================================

nes_matrix <- gsea_results %>%
  select(
    Pathway_clean,
    Contrast,
    NES
  ) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from = Contrast,
    values_from = NES
  )

nes_mat <- as.data.frame(nes_matrix)

rownames(nes_mat) <- nes_mat$Pathway_clean

nes_mat$Pathway_clean <- NULL

nes_mat <- as.matrix(nes_mat)


# ============================================================
# 8.11 Keep pathways significantly enriched in at least one
#      contrast
# ============================================================

sig_pathways <- gsea_results %>%
  filter(padj < 0.05) %>%
  count(Pathway_clean, sort = TRUE)

sig_names <- sig_pathways$Pathway_clean

nes_sig <- nes_mat[
  rownames(nes_mat) %in% sig_names,
  ,
  drop = FALSE
]


# ============================================================
# 8.12 Order columns biologically
# ============================================================

desired_columns <- c(
  "sh190|Disease",
  "sh190|Treatment_vs_KO",
  "sh190|Treatment_vs_WT",
  "MOE|Disease",
  "MOE|Treatment_vs_KO",
  "MOE|Treatment_vs_WT"
)

desired_columns <- intersect(
  desired_columns,
  colnames(nes_sig)
)

nes_sig <- nes_sig[, desired_columns, drop = FALSE]


# ============================================================
# 8.13 Scale pathway NES for visualization
# ============================================================

nes_scaled <- t(
  scale(t(nes_sig))
)

# Remove pathways with undefined scaling
nes_scaled <- nes_scaled[
  complete.cases(nes_scaled),
  ,
  drop = FALSE
]


# ============================================================
# 8.14 Pathway NES heatmap
# ============================================================

library(ComplexHeatmap)
library(circlize)

ht <- Heatmap(
  nes_scaled,
  
  name = "NES\n(z-score)",
  
  col = colorRamp2(
    c(-2, 0, 2),
    c("navy", "white", "firebrick3")
  ),
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 7
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  column_names_rot = 45,
  
  row_title = "GSEA pathways",
  
  column_title =
    "Genome-wide pathway enrichment across disease and treatment",
  
  heatmap_legend_param = list(
    title = "NES\n(z-score)"
  )
)


pdf(
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_NES_heatmap.pdf"
    )
  ),
  width = 12,
  height = 14
)

draw(ht)

dev.off()


# ============================================================
# 8.15 Also save the unscaled NES matrix
# ============================================================

write.csv(
  nes_sig,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_NES_matrix.csv"
    )
  )
)


message(
  "\n============================================================\n",
  "GENOME-WIDE GSEA COMPLETE\n",
  "============================================================\n",
  "Results: ",
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_all_pathways.csv"
    )
  ),
  "\nHeatmap: ",
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_GENOME_WIDE_GSEA_NES_heatmap.pdf"
    )
  ),
  "\n============================================================"
)

# ============================================================
# Focused GSEA heatmap:
# Muscle regeneration + Calcium response + mTORC1
# ============================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ------------------------------------------------------------
# 1. Define pathways of interest
# ------------------------------------------------------------

pathways_of_interest <- c(
  "GOBP_SKELETAL_MUSCLE_TISSUE_REGENERATION",
  "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION",
  "HALLMARK_MTORC1_SIGNALING", "GOBP_CALCIUM_MEDIATED_SIGNALING"
)


# ------------------------------------------------------------
# 2. Extract these pathways from the complete GSEA results
# ------------------------------------------------------------

focused_gsea <- gsea_results %>%
  filter(Pathway %in% pathways_of_interest) %>%
  select(
    Pathway,
    Contrast,
    NES,
    pval,
    padj,
    size
  )


# Check
print(focused_gsea)


# ------------------------------------------------------------
# 3. Make pathway names publication-friendly
# ------------------------------------------------------------

focused_gsea <- focused_gsea %>%
  mutate(
    Pathway = case_when(
      Pathway == "GOBP_SKELETAL_MUSCLE_TISSUE_REGENERATION" ~
        "Skeletal muscle tissue regeneration",
      
      Pathway == "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION" ~
        "Cellular response to calcium ion",
      
      Pathway == "HALLMARK_MTORC1_SIGNALING" ~
        "mTORC1 signaling",
      
      Pathway == "GOBP_CALCIUM_MEDIATED_SIGNALING" ~
        "Caclcium mediated signaling",
      
      TRUE ~ Pathway
    )
  )


# ------------------------------------------------------------
# 4. Set the exact order of the six contrasts
# ------------------------------------------------------------

contrast_order <- c(
  "sh190|Disease",
  "sh190|Treatment_vs_KO",
  "sh190|Treatment_vs_WT",
  "MOE|Disease",
  "MOE|Treatment_vs_KO",
  "MOE|Treatment_vs_WT"
)

focused_gsea$Contrast <- factor(
  focused_gsea$Contrast,
  levels = contrast_order
)


# ------------------------------------------------------------
# 5. Create NES matrix
# ------------------------------------------------------------

nes_matrix <- focused_gsea %>%
  select(
    Pathway,
    Contrast,
    NES
  ) %>%
  pivot_wider(
    names_from = Contrast,
    values_from = NES
  )

nes_mat <- as.data.frame(nes_matrix)

rownames(nes_mat) <- nes_mat$Pathway

nes_mat$Pathway <- NULL

nes_mat <- as.matrix(nes_mat)

# Ensure exact column order
nes_mat <- nes_mat[
  ,
  contrast_order,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Create significance annotation
# ------------------------------------------------------------

sig_matrix <- focused_gsea %>%
  mutate(
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  select(
    Pathway,
    Contrast,
    significance
  ) %>%
  pivot_wider(
    names_from = Contrast,
    values_from = significance
  )

sig_mat <- as.data.frame(sig_matrix)

rownames(sig_mat) <- sig_mat$Pathway

sig_mat$Pathway <- NULL

sig_mat <- as.matrix(sig_mat)

sig_mat <- sig_mat[
  rownames(nes_mat),
  contrast_order,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Heatmap
# ------------------------------------------------------------

ht <- Heatmap(
  nes_mat,
  
  name = "NES",
  
  col = colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  ),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "left",
  
  row_names_gp = gpar(
    fontsize = 11
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  column_names_rot = 45,
  
  # ----------------------------------------------------------
  # White borders around individual tiles
  # ----------------------------------------------------------
  
  rect_gp = gpar(
    col = "white",
    lwd = 1.5
  ),
  
  # ----------------------------------------------------------
  # Significance stars
  # ----------------------------------------------------------
  
  cell_fun = function(
    j, i, x, y, width, height, fill
  ) {
    
    grid.text(
      sig_mat[i, j],
      x,
      y,
      gp = gpar(
        fontsize = 12,
        fontface = "bold"
      )
    )
    
  },
  
  # ----------------------------------------------------------
  # Separate sh190 and MOE
  # ----------------------------------------------------------
  
  column_split = factor(
    c(
      "sh190",
      "sh190",
      "sh190",
      "MOE",
      "MOE",
      "MOE"
    ),
    levels = c("sh190", "MOE")
  ),
  
  column_title = "Pathway-level GSEA across disease and treatment",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  heatmap_legend_param = list(
    title = "NES"
  )
)


# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

pdf(
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_focused_GSEA_3_pathways_heatmap.pdf"
    )
  ),
  width = 11,
  height = 4
)

draw(ht)

dev.off()


# Display in R
draw(ht)

calcium_leading_edge <- gsea_results %>%
  filter(
    Pathway == "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION"
  ) %>%
  select(
    Contrast,
    NES,
    pval,
    padj,
    leadingEdge
  )

calcium_leading_edge
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)

calcium_genes_long <- gsea_results %>%
  filter(
    Pathway == "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION"
  ) %>%
  select(
    Contrast,
    NES,
    pval,
    padj,
    leadingEdge
  ) %>%
  unnest_longer(leadingEdge) %>%
  rename(Ensembl = leadingEdge)

calcium_genes_long$Symbol <- mapIds(
  org.Mm.eg.db,
  keys = calcium_genes_long$Ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

calcium_genes_long %>%
  select(
    Contrast,
    Symbol,
    Ensembl,
    NES,
    pval,
    padj
  )
write.csv(
  calcium_genes_long,
  file.path(
    out_dir,
    paste0(Sys.Date(), "_calcium_response_leading_edge_genes.csv")
  ),
  row.names = FALSE
)

message(
  "Saved: ",
  file.path(
    out_dir,
    paste0(Sys.Date(), "_calcium_response_leading_edge_genes.csv")
  )
)
calcium_sh190_treat <- calcium_genes_long %>%
  filter(
    Contrast == "sh190|Treatment_vs_KO"
  )

calcium_sh190_treat
calcium_sh190_check <- res_treated_e %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Ensembl") %>%
  filter(
    Ensembl %in% calcium_sh190_treat$Ensembl
  ) %>%
  select(
    Ensembl,
    log2FoldChange,
    padj
  )

calcium_sh190_check$Symbol <- mapIds(
  org.Mm.eg.db,
  keys = calcium_sh190_check$Ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

calcium_sh190_check %>%
  arrange(desc(log2FoldChange))




# ============================================================
# KEGG CALCIUM SIGNALING PATHWAY — mm04020
# ============================================================

library(KEGGREST)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(tibble)

# ------------------------------------------------------------
# 1. Get all mouse genes belonging to KEGG Calcium signaling
# ------------------------------------------------------------

calcium_kegg_links <- keggLink(
  "mmu",
  "path:mmu04020"
)


calcium_kegg_entrez <- sub(
  "^mmu:",
  "",
  unname(calcium_kegg_links)
)

calcium_kegg_entrez <- unique(calcium_kegg_entrez)

message(
  "Number of KEGG Calcium signaling genes: ",
  length(calcium_kegg_entrez)
)
# ------------------------------------------------------------
# 2. Convert KEGG Entrez IDs -> Ensembl + Symbol
# ------------------------------------------------------------

calcium_kegg_annotation <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = calcium_kegg_entrez,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "ENTREZID"
)
calcium_kegg_annotation <- calcium_kegg_annotation %>%
  filter(
    !is.na(SYMBOL),
    !is.na(ENSEMBL)
  ) %>%
  distinct(ENTREZID, .keep_all = TRUE)

# ------------------------------------------------------------
# 3. Save the complete KEGG gene list
# ------------------------------------------------------------

write.csv(
  calcium_kegg_annotation,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_mm04020_genes.csv"
    )
  ),
  row.names = FALSE
)

print(calcium_kegg_annotation)

message(
  "Saved KEGG Calcium signaling gene list: ",
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_mm04020_genes.csv"
    )
  )
)
# ============================================================
# KEGG GENES PRESENT IN RNA-SEQ
# ============================================================

calcium_kegg_present <- calcium_kegg_annotation %>%
  filter(
    ENSEMBL %in% rownames(stabilized_counts)
  )

message(
  "KEGG Calcium signaling genes in RNA-seq: ",
  nrow(calcium_kegg_present),
  " / ",
  nrow(calcium_kegg_annotation)
)

# Genes missing from RNA-seq
calcium_kegg_missing <- calcium_kegg_annotation %>%
  filter(
    !ENSEMBL %in% rownames(stabilized_counts)
  )

message(
  "KEGG Calcium signaling genes absent from RNA-seq: ",
  nrow(calcium_kegg_missing)
)

write.csv(
  calcium_kegg_present,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_genes_in_RNAseq.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  calcium_kegg_missing,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_genes_missing_from_RNAseq.csv"
    )
  ),
  row.names = FALSE
)
# ============================================================
# KEGG CALCIUM SIGNALING — DESEQ2 RESULTS
# ============================================================

add_de_stats <- function(
    gene_df,
    res,
    contrast_name
) {
  
  tmp <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    select(
      ENSEMBL,
      log2FoldChange,
      padj
    )
  
  colnames(tmp)[2:3] <- c(
    paste0("log2FC_", contrast_name),
    paste0("padj_", contrast_name)
  )
  
  left_join(
    gene_df,
    tmp,
    by = "ENSEMBL"
  )
}

calcium_kegg_stats <- calcium_kegg_present

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_disease_e,
  "sh190_Disease"
)

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_treated_e,
  "sh190_Treatment_vs_KO"
)

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_treated_vs_WT_e,
  "sh190_Treatment_vs_WT"
)

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_disease_n,
  "MOE_Disease"
)

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_treated_n,
  "MOE_Treatment_vs_KO"
)

calcium_kegg_stats <- add_de_stats(
  calcium_kegg_stats,
  res_treated_vs_WT_n,
  "MOE_Treatment_vs_WT"
)

# Order by disease effect
calcium_kegg_stats <- calcium_kegg_stats %>%
  arrange(
    desc(abs(log2FC_sh190_Disease))
  )

write.csv(
  calcium_kegg_stats,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_DESeq2_statistics.csv"
    )
  ),
  row.names = FALSE
)

calcium_kegg_stats


# ============================================================
# KEGG CALCIUM SIGNALING — EXPRESSION HEATMAP
# ============================================================

calcium_kegg_mat <- stabilized_counts[
  calcium_kegg_present$ENSEMBL,
  ,
  drop = FALSE
]

# Log transform
calcium_kegg_mat <- log2(
  calcium_kegg_mat + 1
)

# Gene-wise z-score
calcium_kegg_mat <- t(
  scale(
    t(calcium_kegg_mat)
  )
)

# Replace Ensembl IDs with gene symbols
rownames(calcium_kegg_mat) <- calcium_kegg_present$SYMBOL

# ------------------------------------------------------------
# Sample order
# ------------------------------------------------------------

desired_order <- c(
  "WT_empty",
  "Stim1R304W/+_empty",
  "Stim1R304W/+_sh190",
  "WT_Nacl",
  "Stim1R304W/+_Nacl",
  "Stim1R304W/+_MOE"
)

coldata$condition <- factor(
  coldata$condition,
  levels = desired_order
)

samples_ordered <- rownames(
  coldata[
    order(coldata$condition),
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# Column annotation
# ------------------------------------------------------------

col_ha <- HeatmapAnnotation(
  Condition = coldata[
    samples_ordered,
    "condition"
  ],
  col = list(
    Condition = condition_colors
  ),
  annotation_name_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  )
)

# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

ht_calcium_kegg <- Heatmap(
  calcium_kegg_mat[
    ,
    samples_ordered,
    drop = FALSE
  ],
  
  name = "Z-score",
  
  col = colorRamp2(
    c(-2, 0, 2),
    c("navy", "white", "firebrick3")
  ),
  
  top_annotation = col_ha,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  
  row_names_side = "right",
  
  row_names_gp = gpar(
    fontsize = 8,
    fontface = "italic"
  ),
  
  column_names_gp = gpar(
    fontsize = 10
  ),
  
  column_names_rot = 80,
  
  row_title = "KEGG Calcium signaling genes",
  
  row_title_gp = gpar(
    fontsize = 13,
    fontface = "bold"
  ),
  
  column_title =
    "KEGG Calcium signaling pathway (mm04020)",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  )
)

pdf(
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_mm04020_heatmap.pdf"
    )
  ),
  width = 10,
  height = max(
    10,
    0.18 * nrow(calcium_kegg_mat) + 4
  )
)

draw(ht_calcium_kegg)

dev.off()

draw(ht_calcium_kegg)


# ============================================================
# SIGNIFICANT KEGG CALCIUM GENES IN DISEASE
# ============================================================

calcium_kegg_disease <- calcium_kegg_stats %>%
  mutate(
    sig_sh190 = !is.na(padj_sh190_Disease) &
      padj_sh190_Disease < 0.05,
    
    sig_MOE = !is.na(padj_MOE_Disease) &
      padj_MOE_Disease < 0.05,
    
    disease_absFC_sh190 =
      abs(log2FC_sh190_Disease),
    
    disease_absFC_MOE =
      abs(log2FC_MOE_Disease)
  ) %>%
  filter(
    sig_sh190 | sig_MOE
  ) %>%
  arrange(
    desc(
      pmax(
        disease_absFC_sh190,
        disease_absFC_MOE,
        na.rm = TRUE
      )
    )
  )

write.csv(
  calcium_kegg_disease,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_disease_DE_genes.csv"
    )
  ),
  row.names = FALSE
)

calcium_kegg_disease



# ============================================================
# KEGG CALCIUM SIGNALING — RESCUE ANALYSIS
# ============================================================

calculate_kegg_rescue <- function(
    gene_df,
    disease_res,
    treatment_res,
    arm_name,
    log2FC_cutoff = 0.5,
    padj_cutoff = 0.05
) {
  
  disease_df <- as.data.frame(disease_res) %>%
    rownames_to_column("ENSEMBL")
  
  treatment_df <- as.data.frame(treatment_res) %>%
    rownames_to_column("ENSEMBL")
  
  out <- gene_df %>%
    select(
      ENTREZID,
      SYMBOL,
      ENSEMBL
    ) %>%
    left_join(
      disease_df %>%
        select(
          ENSEMBL,
          log2FoldChange,
          padj
        ),
      by = "ENSEMBL"
    ) %>%
    rename(
      log2FC_disease = log2FoldChange,
      padj_disease = padj
    ) %>%
    left_join(
      treatment_df %>%
        select(
          ENSEMBL,
          log2FoldChange,
          padj
        ),
      by = "ENSEMBL"
    ) %>%
    rename(
      log2FC_treatment_vs_KO = log2FoldChange,
      padj_treatment_vs_KO = padj
    ) %>%
    mutate(
      Arm = arm_name
    ) %>%
    filter(
      !is.na(padj_disease),
      padj_disease < padj_cutoff,
      abs(log2FC_disease) >= log2FC_cutoff
    ) %>%
    mutate(
      
      # Reconstruct treatment vs WT
      log2FC_treatment_vs_WT =
        log2FC_disease +
        log2FC_treatment_vs_KO,
      
      # Same rescue definition as your existing pipeline
      rescue_metric =
        100 *
        (
          log2FC_disease -
            log2FC_treatment_vs_WT
        ) /
        log2FC_disease,
      
      rescue_status = cut(
        rescue_metric,
        breaks = c(
          -Inf,
          0,
          30,
          80,
          120,
          Inf
        ),
        labels = c(
          "Worsened",
          "Not rescued",
          "Partially rescued",
          "Rescued",
          "Over-corrected"
        )
      )
    )
  
  out
}

# sh190
calcium_kegg_rescue_sh190 <- calculate_kegg_rescue(
  calcium_kegg_present,
  res_disease_e,
  res_treated_e,
  "sh190"
)

# MOE
calcium_kegg_rescue_MOE <- calculate_kegg_rescue(
  calcium_kegg_present,
  res_disease_n,
  res_treated_n,
  "MOE"
)

calcium_kegg_rescue <- bind_rows(
  calcium_kegg_rescue_sh190,
  calcium_kegg_rescue_MOE
)

write.csv(
  calcium_kegg_rescue,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_rescue_table.csv"
    )
  ),
  row.names = FALSE
)

calcium_kegg_rescue



# ============================================================
# RESCUED KEGG CALCIUM SIGNALING GENES
# ============================================================

calcium_kegg_rescued <- calcium_kegg_rescue %>%
  filter(
    rescue_metric >= 30
  ) %>%
  arrange(
    desc(rescue_metric)
  )

write.csv(
  calcium_kegg_rescued,
  file.path(
    out_dir,
    paste0(
      Sys.Date(),
      "_KEGG_CALCIUM_SIGNALING_rescued_genes.csv"
    )
  ),
  row.names = FALSE
)

calcium_kegg_rescued

\


# ============================================================
# ---- 9. Pathway-annotated expression heatmap from a curated
#          Excel panel (Pathway, Genes columns)
# ============================================================
# Left annotation = Pathway/sub-pathway (from your Excel column, in the
# order it appears there), top annotation = genotype/condition, exactly
# like your original ca_genes ComplexHeatmap. Works for any two-column
# (Pathway, Genes) Excel file, not just the calcium signaling panel.
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
library(RColorBrewer)

build_pathway_annotated_heatmap <- function(xlsx_path,
                                            stabilized_counts,
                                            coldata,
                                            condition_colors,
                                            out_dir,
                                            heatmap_title = "Gene panel",
                                            out_filename = "pathway_heatmap",
                                            row_label_fontsize = 9) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  raw <- read_excel(xlsx_path)
  colnames(raw)[1:2] <- c("Pathway", "GeneSymbol")
  raw$Pathway    <- trimws(as.character(raw$Pathway))
  raw$GeneSymbol <- trimws(as.character(raw$GeneSymbol))
  raw <- raw[!is.na(raw$GeneSymbol) & raw$GeneSymbol != "", ]
  
  # Map each ORIGINAL symbol -> Ensembl, case-robust (as-is, then
  # mouse-cased SYMBOL, then mouse-cased ALIAS) — kept separate from
  # get_geneset_ensembl() so the Pathway label stays attached per gene.
  to_mouse_case <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  unique_symbols <- unique(raw$GeneSymbol)
  orig_to_ensembl <- sapply(unique_symbols, function(sym) {
    e <- suppressMessages(mapIds(org.Mm.eg.db, keys = sym, column = "ENSEMBL",
                                 keytype = "SYMBOL", multiVals = "first"))
    if (is.na(e)) e <- suppressMessages(mapIds(org.Mm.eg.db, keys = to_mouse_case(sym), column = "ENSEMBL",
                                               keytype = "SYMBOL", multiVals = "first"))
    if (is.na(e)) e <- suppressMessages(mapIds(org.Mm.eg.db, keys = to_mouse_case(sym), column = "ENSEMBL",
                                               keytype = "ALIAS", multiVals = "first"))
    unname(e)
  })
  
  raw$Ensembl <- orig_to_ensembl[raw$GeneSymbol]
  missing <- raw[is.na(raw$Ensembl), ]
  if (nrow(missing) > 0) {
    message(nrow(missing), " gene(s) could not be mapped and will be dropped: ",
            paste(missing$GeneSymbol, collapse = ", "))
  }
  raw <- raw[!is.na(raw$Ensembl), ]
  raw <- raw[raw$Ensembl %in% rownames(stabilized_counts), ]
  
  if (nrow(raw) < 2) {
    message("Fewer than 2 genes available in the expression matrix — skipping heatmap.")
    return(invisible(NULL))
  }
  
  raw$Symbol <- suppressMessages(mapIds(org.Mm.eg.db, keys = raw$Ensembl, column = "SYMBOL",
                                        keytype = "ENSEMBL", multiVals = "first"))
  
  mat <- stabilized_counts[raw$Ensembl, , drop = FALSE]
  mat <- t(scale(t(log2(mat + 1))))
  rownames(mat) <- raw$Symbol
  
  # Preserve pathway order as it appears in the Excel file (not alphabetical)
  raw$Pathway <- factor(raw$Pathway, levels = unique(raw$Pathway))
  row_split <- raw$Pathway
  
  desired_order <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190",
                     "WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  coldata$condition <- factor(coldata$condition, levels = desired_order)
  sample_order <- rownames(coldata[order(coldata$condition), ])
  mat <- mat[, sample_order]
  
  # Two-level column split: Arm (sh190 vs MOE) is the OUTER split, gets
  # the wider gap/border; Condition within each arm is the INNER split,
  # gets a thinner gap/border. ComplexHeatmap draws a border at every
  # split boundary when column_split has multiple columns like this.
  group_sh  <- c("WT_empty", "Stim1R304W/+_empty", "Stim1R304W/+_sh190")
  group_moe <- c("WT_Nacl", "Stim1R304W/+_Nacl", "Stim1R304W/+_MOE")
  sample_condition <- as.character(coldata[sample_order, "condition"])
  sample_arm <- ifelse(sample_condition %in% group_sh, "sh190 arm", "MOE arm")
  
  col_split_df <- data.frame(
    Arm       = factor(sample_arm, levels = c("sh190 arm", "MOE arm")),
    Condition = factor(sample_condition, levels = desired_order)
  )
  
  # column_gap needs ONE VALUE PER LEAF COLUMN SLICE (here: one per
  # condition present, in order) — not one per split level. Use a wide
  # gap only at the boundary where the arm actually changes.
  condition_levels_present <- levels(droplevels(col_split_df$Condition))
  arm_per_slice <- sapply(condition_levels_present, function(cond) {
    as.character(col_split_df$Arm[col_split_df$Condition == cond][1])
  })
  gap_values <- rep(1.2, length(condition_levels_present))
  if (length(gap_values) > 1) {
    for (i in seq_len(length(gap_values) - 1)) {
      if (arm_per_slice[i] != arm_per_slice[i + 1]) gap_values[i] <- 4
    }
  }
  column_gap_units <- do.call(unit.c, lapply(gap_values, unit, units = "mm"))
  
  n_pathways <- length(levels(raw$Pathway))
  pathway_colors <- setNames(
    colorRampPalette(brewer.pal(min(max(n_pathways, 3), 8), "Set2"))(n_pathways),
    levels(raw$Pathway)
  )
  row_ha <- rowAnnotation(
    Pathway = raw$Pathway,
    col = list(Pathway = pathway_colors),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    show_legend = TRUE
  )
  
  col_ha <- HeatmapAnnotation(
    Condition = coldata[sample_order, "condition"],
    col = list(Condition = condition_colors),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  ht <- Heatmap(
    mat,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    top_annotation = col_ha,
    left_annotation = row_ha,
    row_split = row_split,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    column_split = col_split_df,
    cluster_column_slices = FALSE,
    column_gap = column_gap_units,   # wide gap where arm changes, thin gap between conditions
    border = TRUE,
    column_title = NULL,                  # condition identity already shown by top annotation colors
    row_names_side = "right",
    row_names_gp = gpar(fontsize = row_label_fontsize, fontface = "italic"),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 80,
    row_title_gp = gpar(fontsize = 7, fontface = "bold")
  )
  
  pdf(file.path(out_dir, paste0(Sys.Date(), "_", out_filename, ".pdf")),
      width = 15, height = max(8, 0.22 * nrow(mat) + 4))
  draw(ht, column_title = heatmap_title, column_title_gp = gpar(fontsize = 13, fontface = "bold"),
       heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  message("Saved pathway heatmap: ", file.path(out_dir, paste0(Sys.Date(), "_", out_filename, ".pdf")))
  invisible(list(matrix = mat, annotation = raw, heatmap = ht))
}

# ---- Run for your curated calcium signaling panel ----
calcium_pathway_panel <- build_pathway_annotated_heatmap(
  xlsx_path         = "260811_output/260813_calcium_signalling_genes.xlsx",
  stabilized_counts = stabilized_counts,
  coldata           = coldata,
  condition_colors  = condition_colors,
  out_dir           = calcium_out,
  heatmap_title     = "Calcium signaling gene panel",
  out_filename      = "calcium_signaling_pathway_heatmap"
)
dir.create("workspace_backups", showWarnings = FALSE)

save.image(
  file = file.path(
    "workspace_backups",
    paste0(Sys.Date(), "_calcium_analysis_workspace.RData")
  )
)

message("Workspace saved successfully!")




# ============================================================
# 20260814_recover_wiped_outputs.R
#
# ONE-TIME RECOVERY after 260811_output/ got wiped by unlink().
# Nothing here re-runs DESeq2, GSEA, or any statistics — it only
# re-draws/re-writes from R objects that are still sitting in your
# current session (confirmed present via ls()).
#
# Run this top to bottom in the SAME R session where those objects
# still exist. Every block is wrapped in exists()/is.null() checks,
# so it's safe to run even if some objects aren't there.
# ============================================================

# ---- 0. Recreate folders (non-destructive — see the fixed pipeline script) ----
dir.create(out_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(calcium_out, recursive = TRUE, showWarnings = FALSE)
dir.create(muscle_out,  recursive = TRUE, showWarnings = FALSE)
dir.create(mtorc_out,   recursive = TRUE, showWarnings = FALSE)

resave_heatmap <- function(ht_obj, path, width = 9, height = 11, ...) {
  if (is.null(ht_obj)) { message("Skipping (NULL/not found): ", path); return(invisible(NULL)) }
  pdf(path, width = width, height = height)
  draw(ht_obj, ...)
  dev.off()
  message("Re-saved: ", path)
}

resave_csv <- function(df, path) {
  if (is.null(df)) { message("Skipping (NULL/not found): ", path); return(invisible(NULL)) }
  write.csv(df, path, row.names = FALSE)
  message("Re-saved: ", path)
}

# ============================================================
# ---- 1. Calcium_response (run_full_geneset_pipeline output) ----
# ============================================================
if (exists("results_calcium")) {
  resave_heatmap(results_calcium$heatmaps$combined,
                 file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_heatmap_combined.pdf")), 9, 11)
  resave_heatmap(results_calcium$heatmaps$sh,
                 file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_heatmap_sh.pdf")), 7, 11)
  resave_heatmap(results_calcium$heatmaps$moe,
                 file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_heatmap_MOE.pdf")), 7, 11)
  resave_csv(results_calcium$heatmaps$table,
             file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_top_stats.csv")))
  resave_csv(results_calcium$rescue$table,
             file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_rescue_table.csv")))
  if (!is.null(results_calcium$rescue_heatmap)) {
    resave_heatmap(results_calcium$rescue_heatmap$combined,
                   file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_RESCUED_heatmap_combined.pdf")), 10, 10)
    resave_heatmap(results_calcium$rescue_heatmap$sh190,
                   file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_RESCUED_heatmap_sh190.pdf")), 8, 10)
    resave_heatmap(results_calcium$rescue_heatmap$moe,
                   file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_RESCUED_heatmap_MOE.pdf")), 8, 10)
    resave_csv(results_calcium$rescue_heatmap$genes,
               file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_rescued_genes_for_heatmap.csv")))
  }
}
resave_csv(if (exists("rescue_calcium")) rescue_calcium$table else NULL,
           file.path(calcium_out, paste0(Sys.Date(), "_Calcium_response_rescue_table_standalone.csv")))

if (exists("barplots_calcium")) {
  message(length(barplots_calcium), " calcium bar plot objects found in memory — ",
          "re-run plot_rescued_gene_barplots() for calcium to regenerate the PNGs ",
          "(cheap, no stats recomputed).")
}

# ============================================================
# ---- 2. Muscle_injury_regeneration ----
# ============================================================
if (exists("results_muscle_regen")) {
  resave_heatmap(results_muscle_regen$heatmaps$combined,
                 file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_heatmap_combined.pdf")), 9, 11)
  resave_heatmap(results_muscle_regen$heatmaps$sh,
                 file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_heatmap_sh.pdf")), 7, 11)
  resave_heatmap(results_muscle_regen$heatmaps$moe,
                 file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_heatmap_MOE.pdf")), 7, 11)
  resave_csv(results_muscle_regen$heatmaps$table,
             file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_top_stats.csv")))
  resave_csv(results_muscle_regen$rescue$table,
             file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_rescue_table.csv")))
  if (!is.null(results_muscle_regen$rescue_heatmap)) {
    resave_heatmap(results_muscle_regen$rescue_heatmap$combined,
                   file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_RESCUED_heatmap_combined.pdf")), 10, 10)
    resave_heatmap(results_muscle_regen$rescue_heatmap$sh190,
                   file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_RESCUED_heatmap_sh190.pdf")), 8, 10)
    resave_heatmap(results_muscle_regen$rescue_heatmap$moe,
                   file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_RESCUED_heatmap_MOE.pdf")), 8, 10)
    resave_csv(results_muscle_regen$rescue_heatmap$genes,
               file.path(muscle_out, paste0(Sys.Date(), "_Muscle_injury_regeneration_rescued_genes_for_heatmap.csv")))
  }
}

# ============================================================
# ---- 3. mTOR_Akt / mTORC pathway ----
# ============================================================
if (exists("results_mtorc")) {
  resave_heatmap(results_mtorc$heatmaps$combined,
                 file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_heatmap_combined.pdf")), 9, 11)
  resave_heatmap(results_mtorc$heatmaps$sh,
                 file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_heatmap_sh.pdf")), 7, 11)
  resave_heatmap(results_mtorc$heatmaps$moe,
                 file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_heatmap_MOE.pdf")), 7, 11)
  resave_csv(results_mtorc$heatmaps$table,
             file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_top_stats.csv")))
  resave_csv(results_mtorc$rescue$table,
             file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_rescue_table.csv")))
  if (!is.null(results_mtorc$rescue_heatmap)) {
    resave_heatmap(results_mtorc$rescue_heatmap$combined,
                   file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_RESCUED_heatmap_combined.pdf")), 10, 10)
    resave_heatmap(results_mtorc$rescue_heatmap$sh190,
                   file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_RESCUED_heatmap_sh190.pdf")), 8, 10)
    resave_heatmap(results_mtorc$rescue_heatmap$moe,
                   file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_RESCUED_heatmap_MOE.pdf")), 8, 10)
    resave_csv(results_mtorc$rescue_heatmap$genes,
               file.path(mtorc_out, paste0(Sys.Date(), "_mTORC_pathway_rescued_genes_for_heatmap.csv")))
  }
}

# ============================================================
# ---- 4. Curated calcium signaling pathway panel (Excel-based) ----
# ============================================================
if (exists("calcium_pathway_panel") && !is.null(calcium_pathway_panel)) {
  resave_heatmap(calcium_pathway_panel$heatmap,
                 file.path(calcium_out, paste0(Sys.Date(), "_calcium_signaling_pathway_heatmap.pdf")),
                 width = 10, height = max(8, 0.22 * nrow(calcium_pathway_panel$matrix) + 4))
  resave_csv(calcium_pathway_panel$annotation,
             file.path(calcium_out, paste0(Sys.Date(), "_calcium_signaling_pathway_gene_annotation.csv")))
}

# ============================================================
# ---- 5. Custom 3-geneset pathway GSEA (fgsea, run_pathway_gsea) ----
# ============================================================
resave_csv(if (exists("pathway_gsea_results")) pathway_gsea_results else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_results.csv")))
if (exists("p_gsea")) {
  ggsave(file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_NES_summary.pdf")), plot = p_gsea, width = 9, height = 8)
  message("Re-saved: ", file.path(out_dir, paste0(Sys.Date(), "_pathway_GSEA_NES_summary.pdf")))
}

# ============================================================
# ---- 6. Genome-wide GSEA (msigdbr Hallmark + GO:BP) ----
# ============================================================
if (exists("gsea_results")) {
  gsea_results_export <- gsea_results %>% dplyr::select(-dplyr::where(is.list))
  resave_csv(gsea_results_export, file.path(out_dir, paste0(Sys.Date(), "_GENOME_WIDE_GSEA_all_pathways.csv")))
}
resave_csv(if (exists("nes_sig")) as.data.frame(nes_sig) else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_GENOME_WIDE_GSEA_NES_matrix.csv")))

# `ht` was reused across multiple heatmap blocks in your session, so it
# currently holds only the LAST one that ran — rebuild both explicitly
# from the underlying matrices (which ARE intact) rather than relying on `ht`.
if (exists("nes_scaled")) {
  ht_genome_wide <- Heatmap(
    nes_scaled, name = "NES\n(z-score)",
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    cluster_rows = TRUE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    row_title = "GSEA pathways",
    column_title = "Genome-wide pathway enrichment across disease and treatment",
    heatmap_legend_param = list(title = "NES\n(z-score)")
  )
  resave_heatmap(ht_genome_wide, file.path(out_dir, paste0(Sys.Date(), "_GENOME_WIDE_GSEA_NES_heatmap.pdf")),
                 width = 12, height = 14)
} else {
  message("nes_scaled not found — can't rebuild the genome-wide NES heatmap without re-running GSEA.")
}

# ============================================================
# ---- 7. Focused 4-pathway GSEA heatmap (muscle/calcium/mTORC1) ----
# ============================================================
resave_csv(if (exists("focused_gsea")) focused_gsea else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_focused_GSEA_pathways_table.csv")))

if (exists("nes_mat") && exists("sig_mat")) {
  ht_focused <- Heatmap(
    nes_mat, name = "NES",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left", row_names_gp = gpar(fontsize = 11),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"), column_names_rot = 45,
    rect_gp = gpar(col = "white", lwd = 1.5),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    },
    column_split = factor(c("sh190","sh190","sh190","MOE","MOE","MOE"), levels = c("sh190","MOE")),
    column_title = "Pathway-level GSEA across disease and treatment",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(title = "NES")
  )
  resave_heatmap(ht_focused, file.path(out_dir, paste0(Sys.Date(), "_focused_GSEA_pathways_heatmap.pdf")), 11, 4)
} else {
  message("nes_mat/sig_mat not found — can't rebuild the focused GSEA heatmap without re-deriving them from gsea_results.")
}

# ============================================================
# ---- 8. KEGG calcium signaling (mm04020) ----
# ============================================================
resave_csv(if (exists("calcium_kegg_annotation")) calcium_kegg_annotation else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_mm04020_genes.csv")))
resave_csv(if (exists("calcium_kegg_present")) calcium_kegg_present else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_genes_in_RNAseq.csv")))
resave_csv(if (exists("calcium_kegg_missing")) calcium_kegg_missing else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_genes_missing_from_RNAseq.csv")))
resave_csv(if (exists("calcium_kegg_stats")) calcium_kegg_stats else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_DESeq2_statistics.csv")))
resave_csv(if (exists("calcium_kegg_disease")) calcium_kegg_disease else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_disease_DE_genes.csv")))
resave_csv(if (exists("calcium_kegg_rescue")) calcium_kegg_rescue else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_rescue_table.csv")))
resave_csv(if (exists("calcium_kegg_rescued")) calcium_kegg_rescued else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_rescued_genes.csv")))
if (exists("ht_calcium_kegg")) {
  resave_heatmap(ht_calcium_kegg, file.path(out_dir, paste0(Sys.Date(), "_KEGG_CALCIUM_SIGNALING_mm04020_heatmap.pdf")),
                 width = 10, height = if (exists("calcium_kegg_mat")) max(10, 0.18 * nrow(calcium_kegg_mat) + 4) else 12)
}

# ============================================================
# ---- 9. Calcium leading-edge genes (per-contrast, earlier extraction) ----
# ============================================================
resave_csv(if (exists("calcium_genes_long")) calcium_genes_long else NULL,
           file.path(out_dir, paste0(Sys.Date(), "_calcium_response_leading_edge_genes.csv")))

message("\n=== Recovery pass complete. Check ", out_dir, " and its subfolders. ===\n")
message("Still outstanding (not lost — just never finished running):")
message(" - Section 10 (focused_leadingedge_stats + focused_pathways_leadingedge_heatmap): ",
        "re-run it now, everything it needs (gsea_results, res_disease_e/n, res_treated_e/n, ",
        "res_treated_vs_WT_e/n) is present.")
message(" - barplots_calcium / barplots_muscle_regen (per-gene PNGs): the plot objects exist ",
        "in memory but re-running plot_rescued_gene_barplots() is cheap and will regenerate ",
        "the PNG files directly.")



# ============================================================
# ---- 10. Genes behind the focused GSEA pathways, across
#           conditions
# ============================================================
# The NES summary tells you a pathway moved — this shows you WHICH genes
# drove that, and how each of them individually behaves across all 6
# genotypes/treatments. Pulls the leadingEdge gene list straight out of
# your genome-wide GSEA results (gsea_results, from the msigdbr/fgsea run)
# for the 4 focused pathways, unions the Disease-contrast leading edges
# per pathway (the genes actually driving the disease-state enrichment),
# and reuses the same heatmap style as section 9.

focused_pathways <- c(
  "GOBP_SKELETAL_MUSCLE_TISSUE_REGENERATION",
  "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION",
  "GOBP_CALCIUM_MEDIATED_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING"
)

focused_pathway_labels <- c(
  "GOBP_SKELETAL_MUSCLE_TISSUE_REGENERATION" = "Skeletal muscle tissue regeneration",
  "GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION"     = "Cellular response to calcium ion",
  "GOBP_CALCIUM_MEDIATED_SIGNALING"           = "Calcium mediated signaling",
  "HALLMARK_MTORC1_SIGNALING"                 = "mTORC1 signaling"
)

# Extract the leading-edge gene list per pathway from the Disease contrasts
# (both arms) — these are the genes actually driving each pathway's
# dysregulation signature, not the full (much larger) gene set membership.
focused_leadingedge_genes <- gsea_results %>%
  filter(Pathway %in% focused_pathways,
         Contrast %in% c("sh190|Disease", "MOE|Disease")) %>%
  select(Pathway, Contrast, leadingEdge) %>%
  tidyr::unnest_longer(leadingEdge) %>%
  rename(Ensembl = leadingEdge) %>%
  distinct(Pathway, Ensembl) %>%
  filter(Ensembl %in% rownames(stabilized_counts))

focused_leadingedge_genes$Symbol <- suppressMessages(mapIds(
  org.Mm.eg.db, keys = focused_leadingedge_genes$Ensembl,
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
))
focused_leadingedge_genes$Pathway <- recode(focused_leadingedge_genes$Pathway, !!!focused_pathway_labels)
focused_leadingedge_genes$Pathway <- factor(focused_leadingedge_genes$Pathway,
                                            levels = unname(focused_pathway_labels))
focused_leadingedge_genes <- focused_leadingedge_genes %>% arrange(Pathway)

# ---- Attach DE stats (log2FC, padj) across all 6 contrasts per gene ----
add_de_stats <- function(gene_df, res, contrast_label) {
  tmp <- as.data.frame(res) %>%
    rownames_to_column("Ensembl") %>%
    select(Ensembl, log2FoldChange, padj)
  colnames(tmp)[2:3] <- c(paste0("log2FC_", contrast_label), paste0("padj_", contrast_label))
  left_join(gene_df, tmp, by = "Ensembl")
}

focused_leadingedge_stats <- focused_leadingedge_genes %>%
  add_de_stats(res_disease_e,        "sh190_Disease") %>%
  add_de_stats(res_treated_e,        "sh190_Treatment_vs_KO") %>%
  add_de_stats(res_treated_vs_WT_e,  "sh190_Treatment_vs_WT") %>%
  add_de_stats(res_disease_n,        "MOE_Disease") %>%
  add_de_stats(res_treated_n,        "MOE_Treatment_vs_KO") %>%
  add_de_stats(res_treated_vs_WT_n,  "MOE_Treatment_vs_WT")

write.csv(focused_leadingedge_stats,
          file.path(out_dir, paste0(Sys.Date(), "_focused_pathways_leadingedge_genes_DEstats.csv")),
          row.names = FALSE)
message("Saved leading-edge gene DE stats: ",
        file.path(out_dir, paste0(Sys.Date(), "_focused_pathways_leadingedge_genes_DEstats.csv")))
load ("workspace_backups/2026-08-21_geneset_pipeline_workspace.RData")



downstream_genes_panel <- build_pathway_annotated_heatmap(
  xlsx_path         = "260811_output/260827_downstream_genes.xlsx",
  stabilized_counts = stabilized_counts,
  coldata           = coldata,
  condition_colors  = condition_colors,
  out_dir           = out_dir,
  heatmap_title     = "Downstream muscle regeneration & mTORC signaling genes",
  out_filename      = "downstream_genes_heatmap"
)
barplots_downstream_genes <- plot_rescued_gene_barplots(
  rescue_selected  = downstream_genes_panel$annotation,
  dds              = dds,
  coldata          = coldata,
  condition_colors = condition_colors,
  out_dir          = file.path(out_dir, "downstream_genes_barplots")
)

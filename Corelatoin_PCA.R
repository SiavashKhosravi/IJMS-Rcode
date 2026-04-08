
## ============================================================
## CHOOSE PATHS (EDIT THIS)
## ============================================================
BASE_DIR <- "C:/Users/khosravs/OneDrive - Boehringer Ingelheim/Desktop/RNA PAPER/RNA data/Transcriptom_ Analysis/Corelation_PCA"  # <-- change to your folder

INPUT_DIR   <- file.path(BASE_DIR, "input")
OUTPUT_DIR  <- file.path(BASE_DIR, "output")
INPUT_RDATA <- file.path(INPUT_DIR, "exprMatrix.Rdata")

## Output subfolders
fig_dir     <- file.path(OUTPUT_DIR, "figures_qc")
results_dir <- file.path(OUTPUT_DIR, "qc_results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## Load the SummarizedExperiment object
if (!file.exists(INPUT_RDATA)) stop("Input file not found: ", INPUT_RDATA)
load(INPUT_RDATA)  # should load object named: exprMatrix

## ============================================================
## BLOCK 1 — Setup & Package Check (optimized; no DESeq2/mapping)
## ============================================================

## --- TUNABLE PARAMETERS ---
base_font_size <- 12
fig_dpi <- 600

## Column in colData that holds time points / groups
group_col <- "MFGroup"  # e.g., Day0, Day1, Day3, Day7, Day10

## Desired timepoint order & colors (editable)
time_levels <- c("Day0", "Day1", "Day3", "Day7", "Day10")
time_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
names(time_colors) <- time_levels

## --- Install helper (single definition; robust) ---
options(stringsAsFactors = FALSE)
set.seed(42)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("Installing package: %s", p))
      # Save you from stale locks
      lockdir <- file.path(.libPaths()[1], paste0("00LOCK-", p))
      if (dir.exists(lockdir)) unlink(lockdir, recursive = TRUE, force = TRUE)
      
      if (bioc) {
        BiocManager::install(p, update = FALSE, ask = FALSE)
      } else {
        install.packages(p, dependencies = TRUE)
      }
    }
  }
}

## --- Packages actually needed for Blocks 1–5 ---
cran_pkgs <- c(
  "ggplot2", "dplyr", "tibble", "tidyr", "readr",
  "patchwork", "scales", "matrixStats", "ggrepel",
  "data.table", "glue", "rlang", "svglite"
)

bioc_pkgs <- c(
  "SummarizedExperiment",
  "ComplexHeatmap",
  "circlize"
)

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tibble); library(tidyr); library(readr)
  library(patchwork); library(scales); library(matrixStats); library(ggrepel)
  library(data.table); library(glue); library(rlang); library(svglite)
  library(SummarizedExperiment)
  library(ComplexHeatmap); library(circlize)
})

## --- Publication theme (unchanged look) ---
theme_pub <- function(base_size = base_font_size) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}
theme_set(theme_pub())

## --- Helpers for consistent timepoint scales ---
scale_color_time <- function() scale_color_manual(values = time_colors, drop = FALSE)
scale_fill_time  <- function() scale_fill_manual(values = time_colors, drop = FALSE)

## --- Save ggplot helper: now saves PNG + PDF + SVG ---
save_plot <- function(plot, filename,
                      width = 6.5, height = 5.0, units = "in",
                      dpi = fig_dpi,
                      devices = c("png", "pdf", "svg")) {
  for (device in devices) {
    out_path <- file.path(fig_dir, sprintf("%s.%s", filename, device))
    if (device == "png") {
      ggsave(out_path, plot = plot, width = width, height = height, units = units,
             dpi = dpi, limitsize = FALSE)
    } else if (device == "pdf") {
      ggsave(out_path, plot = plot, width = width, height = height, units = units,
             limitsize = FALSE)
    } else if (device == "svg") {
      # uses svglite for high quality if available
      ggsave(out_path, plot = plot, width = width, height = height, units = units,
             limitsize = FALSE, device = "svg")
    } else {
      stop("Unsupported device: ", device)
    }
    message(glue("Saved: {out_path}"))
  }
  invisible(TRUE)
}

## --- Object existence check ---
if (!exists("exprMatrix")) stop("Object 'exprMatrix' not found after loading Rdata.")
if (!is(exprMatrix, "SummarizedExperiment")) stop("'exprMatrix' is not a SummarizedExperiment.")

## --- Session info ---
sess_file <- file.path(results_dir, "sessionInfo.txt")
capture.output(utils::sessionInfo(), file = sess_file)
message(glue("Session info saved to: {sess_file}"))

## ============================================================
## BLOCK 2 — Sanity checks & overview (kept; small robustness)
## ============================================================

get_col <- function(se, col) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!col %in% colnames(cd)) {
    stop(glue("Column '{col}' not found in colData. Available: {paste(colnames(cd), collapse=', ')}"))
  }
  cd[[col]]
}

overview_txt       <- file.path(results_dir, "object_overview.txt")
assay_presence_txt <- file.path(results_dir, "assay_presence.txt")
sample_head_csv    <- file.path(results_dir, "sample_metadata_head.csv")
group_counts_csv   <- file.path(results_dir, "group_counts.csv")

ngenes   <- nrow(exprMatrix)
nsamples <- ncol(exprMatrix)

assay_names <- SummarizedExperiment::assayNames(exprMatrix)
expected_assays <- c("counts", "tpm", "normCounts", "vstTransform")
assay_presence  <- setNames(expected_assays %in% assay_names, expected_assays)

rn_genes   <- rownames(exprMatrix)
cn_samples <- colnames(exprMatrix)
cd         <- SummarizedExperiment::colData(exprMatrix)
consistent_cols <- identical(cn_samples, rownames(cd))

dup_genes   <- any(duplicated(rn_genes))
dup_samples <- any(duplicated(cn_samples))

any_na_assay <- function(se, nm) {
  if (!(nm %in% SummarizedExperiment::assayNames(se))) return(NA)
  anyNA(SummarizedExperiment::assay(se, nm))
}
na_counts     <- any_na_assay(exprMatrix, "counts")
na_tpm        <- any_na_assay(exprMatrix, "tpm")
na_normCounts <- any_na_assay(exprMatrix, "normCounts")
na_vst        <- any_na_assay(exprMatrix, "vstTransform")

ensembl_like <- grepl("^ENSSSCG", rn_genes)
prop_ens <- mean(ensembl_like) * 100

## Group / timepoint handling
group_vec <- get_col(exprMatrix, group_col)
unique_groups <- unique(as.character(group_vec))
lvl_order <- c(time_levels, setdiff(sort(unique_groups), time_levels))
group_fac <- factor(as.character(group_vec), levels = lvl_order)

SummarizedExperiment::colData(exprMatrix)[[group_col]] <- group_fac

tmp_df <- tibble(sample = cn_samples, grp = group_fac)
group_df <- tmp_df %>%
  dplyr::count(grp, name = "n_samples") %>%
  dplyr::mutate(grp = factor(grp, levels = lvl_order)) %>%
  dplyr::arrange(grp)
names(group_df)[names(group_df) == "grp"] <- group_col
readr::write_csv(group_df, group_counts_csv)

## Minimal library checks if 'counts' exists
libsize_summary <- NULL
zero_lib <- NULL
if ("counts" %in% assay_names) {
  counts_mat <- SummarizedExperiment::assay(exprMatrix, "counts")
  lib_sizes <- colSums(counts_mat, na.rm = TRUE)
  zero_lib <- which(lib_sizes == 0)
  libsize_summary <- summary(lib_sizes)
}

## Write head of sample metadata (avoid duplicate "sampleId")
cd <- SummarizedExperiment::colData(exprMatrix)
cd_df <- as.data.frame(cd)
rn_var <- if ("sampleId" %in% colnames(cd_df)) ".sampleId" else "sampleId"
cd_df %>%
  tibble::rownames_to_column(var = rn_var) %>%
  head(20) %>%
  readr::write_csv(sample_head_csv)

assay_lines <- c(
  glue("Assays present: {paste(assay_names, collapse=', ')}"),
  glue("Expected assays presence:"),
  paste0(" - counts : ", ifelse(is.na(assay_presence[["counts"]]), "NA", assay_presence[["counts"]])),
  paste0(" - tpm : ", ifelse(is.na(assay_presence[["tpm"]]), "NA", assay_presence[["tpm"]])),
  paste0(" - normCounts : ", ifelse(is.na(assay_presence[["normCounts"]]), "NA", assay_presence[["normCounts"]])),
  paste0(" - vstTransform : ", ifelse(is.na(assay_presence[["vstTransform"]]), "NA", assay_presence[["vstTransform"]])),
  "",
  glue("Any NA in assays (NA means assay missing):"),
  paste0(" - counts : ", na_counts),
  paste0(" - tpm : ", na_tpm),
  paste0(" - normCounts : ", na_normCounts),
  paste0(" - vstTransform : ", na_vst)
)
writeLines(assay_lines, con = assay_presence_txt)

overview_lines <- c(
  glue("Date: {Sys.time()}"),
  glue("Dimensions: genes = {ngenes}, samples = {nsamples}"),
  glue("colData columns: {ncol(cd)}"),
  glue("First 10 assay names: {paste(head(assay_names, 10), collapse=', ')}"),
  glue("Column names aligned with colData rownames: {consistent_cols}"),
  glue("Duplicated gene IDs: {dup_genes}"),
  glue("Duplicated sample IDs: {dup_samples}"),
  glue("Proportion of gene IDs that look like Ensembl Sus scrofa (ENSSSCG*): {sprintf('%.1f', prop_ens)}%"),
  glue("Grouping column used: '{group_col}'"),
  glue("Defined timepoint order: {paste(time_levels, collapse=' -> ')}"),
  if (!is.null(libsize_summary)) c("Library size summary (counts per sample):", capture.output(print(libsize_summary)))
  else "Library size summary: counts assay not present.",
  if (!is.null(zero_lib) && length(zero_lib) > 0)
    glue("WARNING: Samples with zero library size: {paste(names(zero_lib), collapse=', ')}")
  else "No zero-library samples detected."
)
writeLines(unlist(overview_lines), con = overview_txt)

message(glue("SUMMARY: {ngenes} genes x {nsamples} samples"))
message(glue("Assays: {paste(assay_names, collapse=', ')}"))
if (!all(assay_presence, na.rm = TRUE)) message("NOTE: Not all expected assays are present. See assay_presence.txt.")
if (!consistent_cols) warning("colnames(exprMatrix) do not match rownames(colData).")
if (dup_genes) warning("Duplicated gene IDs detected.")
if (dup_samples) warning("Duplicated sample IDs detected.")
if (!is.null(zero_lib) && length(zero_lib) > 0) warning(glue("Found {length(zero_lib)} samples with zero library size."))

message(glue("Wrote: {overview_txt}"))
message(glue("Wrote: {assay_presence_txt}"))
message(glue("Wrote: {sample_head_csv}"))
message(glue("Wrote: {group_counts_csv}"))

## ============================================================
## BLOCK 3 — Per-sample QC metrics & plots (fixes + same visuals)
## ============================================================

tpm_min_thresh <- 1
use_log10_libsize <- TRUE
point_size <- 1.4
point_alpha <- 0.7
jitter_width <- 0.15
violin_alpha <- 0.8
label_outliers <- TRUE
nmads_threshold <- 3
rank_label_n <- 3

label_si <- function(...) scales::label_number(scale_cut = scales::cut_short_scale(), ...)

assay_names <- SummarizedExperiment::assayNames(exprMatrix)
if (!("counts" %in% assay_names)) stop("Assay 'counts' is required for Block 3.")
has_tpm <- "tpm" %in% assay_names

counts_mat <- SummarizedExperiment::assay(exprMatrix, "counts")
tpm_mat <- if (has_tpm) SummarizedExperiment::assay(exprMatrix, "tpm") else NULL

sample_ids <- colnames(exprMatrix)
cd <- as.data.frame(SummarizedExperiment::colData(exprMatrix))
stopifnot(group_col %in% colnames(cd))
group_vec <- cd[[group_col]]

lib_size <- colSums(counts_mat, na.rm = TRUE)
detected_counts <- colSums(counts_mat > 0, na.rm = TRUE)
pct_zero_counts <- colMeans(counts_mat == 0, na.rm = TRUE) * 100

detected_tpm <- if (has_tpm) {
  colSums(tpm_mat >= tpm_min_thresh, na.rm = TRUE)
} else rep(NA_integer_, length(sample_ids))

sample_qc <- tibble::tibble(
  sample_id = sample_ids,
  !!group_col := factor(group_vec, levels = levels(cd[[group_col]])),
  lib_size = as.numeric(lib_size),
  detected_counts = as.integer(detected_counts),
  detected_tpm = as.integer(detected_tpm),
  pct_zero_counts = as.numeric(pct_zero_counts)
)

## Robust z via MAD (fixed condition; same behavior intended)
robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  md  <- mad(x, constant = 1.4826, na.rm = TRUE)
  if (is.na(md) || md == 0) return(rep(0, length(x)))
  abs(x - med) / md
}

sample_qc <- sample_qc %>%
  mutate(
    z_libsize  = robust_z(lib_size),
    z_detected = robust_z(detected_counts),
    z_pctzero  = robust_z(pct_zero_counts),
    out_libsize  = z_libsize  > nmads_threshold,
    out_detected = z_detected > nmads_threshold,
    out_pctzero  = z_pctzero  > nmads_threshold
  )

out_csv_metrics <- file.path(results_dir, "sample_qc_metrics.csv")
readr::write_csv(sample_qc, out_csv_metrics)

group_sym <- rlang::sym(group_col)
group_stats <- sample_qc %>%
  group_by(!!group_sym) %>%
  summarise(
    n = n(),
    lib_size_median = median(lib_size, na.rm = TRUE),
    lib_size_IQR    = IQR(lib_size, na.rm = TRUE),
    detected_median = median(detected_counts, na.rm = TRUE),
    detected_IQR    = IQR(detected_counts, na.rm = TRUE),
    pctzero_median  = median(pct_zero_counts, na.rm = TRUE),
    pctzero_IQR     = IQR(pct_zero_counts, na.rm = TRUE),
    .groups = "drop"
  )
out_csv_group <- file.path(results_dir, "sample_qc_group_stats.csv")
readr::write_csv(group_stats, out_csv_group)

## Plot 1: Violin + jitter (same visuals; fixed patchwork syntax)
p_lib <- ggplot(sample_qc, aes(x = !!group_sym, y = lib_size, fill = !!group_sym)) +
  geom_violin(alpha = violin_alpha, color = NA, trim = FALSE) +
  geom_jitter(aes(color = !!group_sym), width = jitter_width, size = point_size, alpha = point_alpha) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white") +
  scale_fill_time() + scale_color_time() +
  labs(x = group_col, y = "Library size (reads)", title = "Library size by timepoint") +
  theme_pub()

if (use_log10_libsize) {
  p_lib <- p_lib + scale_y_continuous(trans = "log10", labels = label_si())
} else {
  p_lib <- p_lib + scale_y_continuous(labels = label_si())
}

p_det <- ggplot(sample_qc, aes(x = !!group_sym, y = detected_counts, fill = !!group_sym)) +
  geom_violin(alpha = violin_alpha, color = NA, trim = FALSE) +
  geom_jitter(aes(color = !!group_sym), width = jitter_width, size = point_size, alpha = point_alpha) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white") +
  scale_fill_time() + scale_color_time() +
  labs(x = group_col, y = "Detected genes (counts > 0)", title = "Detected genes by timepoint") +
  theme_pub()

p_zero <- ggplot(sample_qc, aes(x = !!group_sym, y = pct_zero_counts, fill = !!group_sym)) +
  geom_violin(alpha = violin_alpha, color = NA, trim = FALSE) +
  geom_jitter(aes(color = !!group_sym), width = jitter_width, size = point_size, alpha = point_alpha) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white") +
  scale_fill_time() + scale_color_time() +
  labs(x = group_col, y = "% zero counts per sample", title = "Sparsity (% zeros) by timepoint") +
  theme_pub()

if (label_outliers) {
  p_lib <- p_lib +
    ggrepel::geom_text_repel(
      data = subset(sample_qc, out_libsize),
      aes(label = sample_id, y = lib_size),
      size = 3, max.overlaps = 15, min.segment.length = 0
    )
  p_det <- p_det +
    ggrepel::geom_text_repel(
      data = subset(sample_qc, out_detected),
      aes(label = sample_id, y = detected_counts),
      size = 3, max.overlaps = 15, min.segment.length = 0
    )
  p_zero <- p_zero +
    ggrepel::geom_text_repel(
      data = subset(sample_qc, out_pctzero),
      aes(label = sample_id, y = pct_zero_counts),
      size = 3, max.overlaps = 15, min.segment.length = 0
    )
}

p_violin <- (p_lib | p_det | p_zero) + plot_layout(guides = "collect")

save_plot(p_violin, "qc_metrics_violin_18_11", width = 11, height = 4.2)

## Plot 2: Ranked scatter (same visuals; fixed patchwork syntax)
rank_df <- sample_qc %>%
  arrange(desc(lib_size)) %>% mutate(rank_lib  = row_number()) %>%
  arrange(desc(detected_counts)) %>% mutate(rank_det  = row_number()) %>%
  arrange(pct_zero_counts) %>% mutate(rank_zero = row_number())

p_rank_lib <- ggplot(rank_df, aes(x = rank_lib, y = lib_size, color = !!group_sym)) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_time() +
  labs(x = "Rank by library size (desc)", y = "Library size (reads)", title = "Ranked: Library size") +
  theme_pub()

if (use_log10_libsize) {
  p_rank_lib <- p_rank_lib + scale_y_continuous(trans = "log10", labels = label_si())
} else {
  p_rank_lib <- p_rank_lib + scale_y_continuous(labels = label_si())
}

p_rank_det <- ggplot(rank_df, aes(x = rank_det, y = detected_counts, color = !!group_sym)) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_time() +
  labs(x = "Rank by detected genes (desc)", y = "Detected genes (counts > 0)", title = "Ranked: Detected genes") +
  theme_pub()

p_rank_zero <- ggplot(rank_df, aes(x = rank_zero, y = pct_zero_counts, color = !!group_sym)) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_time() +
  labs(x = "Rank by % zeros (asc)", y = "% zero counts", title = "Ranked: % zeros") +
  theme_pub()

if (label_outliers && rank_label_n > 0) {
  lab_lib  <- rank_df %>% slice_max(order_by = lib_size, n = rank_label_n, with_ties = FALSE)
  lab_det  <- rank_df %>% slice_max(order_by = detected_counts, n = rank_label_n, with_ties = FALSE)
  lab_zero <- rank_df %>% slice_min(order_by = pct_zero_counts, n = rank_label_n, with_ties = FALSE)
  
  p_rank_lib  <- p_rank_lib  + ggrepel::geom_text_repel(data = lab_lib,  aes(label = sample_id), size = 3, max.overlaps = 20)
  p_rank_det  <- p_rank_det  + ggrepel::geom_text_repel(data = lab_det,  aes(label = sample_id), size = 3, max.overlaps = 20)
  p_rank_zero <- p_rank_zero + ggrepel::geom_text_repel(data = lab_zero, aes(label = sample_id), size = 3, max.overlaps = 20)
}

p_ranked <- (p_rank_lib | p_rank_det | p_rank_zero) + plot_layout(guides = "collect")

save_plot(p_ranked, "qc_metrics_ranked_18_11", width = 11, height = 4.2)

message(glue("Saved per-sample QC: {out_csv_metrics}"))
message(glue("Saved per-group stats: {out_csv_group}"))

## ============================================================
## BLOCK 4 — Expression distributions & RLE plots (same visuals)
## ============================================================

use_vst_for_density <- TRUE
log_base <- 2
pseudocount <- 1
facet_density <- FALSE
alpha_density <- 0.25
line_width_density <- 0.6
violin_alpha_fill <- 0.85
trim_violins <- FALSE
rle_show_points <- FALSE

label_si <- function(...) scales::label_number(scale_cut = scales::cut_short_scale(), ...)

assay_names <- SummarizedExperiment::assayNames(exprMatrix)
has_vst  <- "vstTransform" %in% assay_names
has_norm <- "normCounts" %in% assay_names

expr_for_density <- NULL
expr_label <- NULL

if (use_vst_for_density && has_vst) {
  expr_for_density <- SummarizedExperiment::assay(exprMatrix, "vstTransform")
  expr_label <- "VST"
} else if (has_norm) {
  nc <- SummarizedExperiment::assay(exprMatrix, "normCounts")
  expr_for_density <- log(nc + pseudocount, base = log_base)
  expr_label <- glue("log{log_base}(normCounts + {pseudocount})")
} else if ("counts" %in% assay_names) {
  counts <- SummarizedExperiment::assay(exprMatrix, "counts")
  lib_size <- colSums(counts, na.rm = TRUE)
  cpm <- sweep(counts, 2, lib_size, FUN = "/") * 1e6
  expr_for_density <- log(cpm + pseudocount, base = log_base)
  expr_label <- glue("log{log_base}(CPM + {pseudocount})")
} else {
  stop("No suitable assay found for distribution plots (need vstTransform or normCounts or counts).")
}

sample_ids <- colnames(exprMatrix)
cd <- as.data.frame(SummarizedExperiment::colData(exprMatrix))
stopifnot(group_col %in% colnames(cd))
group_vec <- factor(cd[[group_col]], levels = levels(cd[[group_col]]))

set.seed(42)
ngenes <- nrow(expr_for_density)
max_genes_for_density <- 25000
sel_idx <- if (ngenes > max_genes_for_density) sample(ngenes, max_genes_for_density) else seq_len(ngenes)
expr_sub <- expr_for_density[sel_idx, , drop = FALSE]

expr_mean   <- colMeans(expr_sub, na.rm = TRUE)
expr_median <- matrixStats::colMedians(expr_sub, na.rm = TRUE)
expr_sd     <- matrixStats::colSds(expr_sub, na.rm = TRUE)

summary_df <- tibble(
  sample_id = sample_ids,
  !!group_col := group_vec,
  expr_mean   = as.numeric(expr_mean),
  expr_median = as.numeric(expr_median),
  expr_sd     = as.numeric(expr_sd),
  assay_used  = expr_label
)
readr::write_csv(summary_df, file.path(results_dir, "expression_summary_stats.csv"))

## Build long_df (use data.table melt for speed; same downstream ggplot)
dt <- data.table::as.data.table(expr_sub)
data.table::setnames(dt, sample_ids)
dt[, gene := sel_idx]
long_df <- data.table::melt(dt, id.vars = "gene", variable.name = "sample_id", value.name = "expr")
long_df <- as_tibble(long_df) %>%
  left_join(tibble(sample_id = sample_ids, !!group_col := group_vec), by = "sample_id")

p_density <- ggplot(long_df, aes(x = expr, color = !!sym(group_col), group = sample_id)) +
  geom_density(linewidth = line_width_density, alpha = alpha_density, adjust = 1) +
  scale_color_manual(values = time_colors, drop = FALSE) +
  labs(
    x = expr_label, y = "Density",
    title = "Expression distributions across samples",
    subtitle = "Each curve = one sample; colored by timepoint"
  ) +
  theme_pub()

if (facet_density) {
  p_density <- p_density + facet_wrap(vars(!!sym(group_col)), scales = "free_y")
}

save_plot(p_density, "expr_density_panel_18_11", width = 8.5, height = 5.5)

p_violin_expr <- ggplot(long_df, aes(x = !!sym(group_col), y = expr, fill = !!sym(group_col))) +
  geom_violin(alpha = violin_alpha_fill, color = NA, trim = trim_violins) +
  scale_fill_manual(values = time_colors, drop = FALSE) +
  labs(x = group_col, y = expr_label, title = "Expression distribution by timepoint") +
  theme_pub()

save_plot(p_violin_expr, "expr_violin_panel_18_11", width = 6.5, height = 5.0)

## RLE computation
if (has_norm) {
  rle_mat <- log(SummarizedExperiment::assay(exprMatrix, "normCounts") + pseudocount, base = log_base)
  rle_label <- glue("RLE on log{log_base}(normCounts + {pseudocount})")
} else if (has_vst) {
  rle_mat <- SummarizedExperiment::assay(exprMatrix, "vstTransform")
  rle_label <- "RLE on VST"
} else if ("counts" %in% assay_names) {
  counts <- SummarizedExperiment::assay(exprMatrix, "counts")
  lib_size <- colSums(counts, na.rm = TRUE)
  cpm <- sweep(counts, 2, lib_size, FUN = "/") * 1e6
  rle_mat <- log(cpm + pseudocount, base = log_base)
  rle_label <- glue("RLE on log{log_base}(CPM + {pseudocount})")
} else {
  stop("No suitable assay available for RLE computation.")
}

gene_medians <- matrixStats::rowMedians(rle_mat, na.rm = TRUE)
rle_vals <- rle_mat - gene_medians

dt_rle <- data.table::as.data.table(rle_vals)
data.table::setnames(dt_rle, sample_ids)
dt_rle[, gene := rownames(exprMatrix)]
rle_long <- data.table::melt(dt_rle, id.vars = "gene", variable.name = "sample_id", value.name = "rle")
rle_long <- as_tibble(rle_long) %>%
  left_join(tibble(sample_id = sample_ids, !!group_col := group_vec), by = "sample_id")

p_rle_box <- ggplot(rle_long, aes(x = sample_id, y = rle)) +
  geom_boxplot(outlier.size = 0.4, width = 0.5, color = "grey30", fill = "grey85") +
  labs(
    x = "Sample",
    y = "Relative Log Expression (centered per gene)",
    title = "RLE plot (per sample)",
    subtitle = rle_label
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

if (rle_show_points) {
  p_rle_box <- p_rle_box + geom_jitter(size = 0.2, alpha = 0.1, width = 0.2)
}

save_plot(p_rle_box, "rle_boxplot_18_11", width = 11, height = 5.5)

p_rle_violin <- ggplot(rle_long, aes(x = !!sym(group_col), y = rle, fill = !!sym(group_col))) +
  geom_violin(alpha = 0.85, color = NA, trim = FALSE) +
  scale_fill_manual(values = time_colors, drop = FALSE) +
  labs(
    x = group_col,
    y = "Relative Log Expression",
    title = "RLE by timepoint",
    subtitle = rle_label
  ) +
  theme_pub()

save_plot(p_rle_violin, "rle_violin_by_group_18_11", width = 7.0, height = 5.2)

message("Saved: expr_density_panel, expr_violin_panel, rle_boxplot, rle_violin_by_group (png/pdf/svg)")
message("Saved expression summary stats: qc_results/expression_summary_stats.csv")

## ============================================================
## BLOCK 5 — PCA & Sample–Sample Correlation Heatmap (same visuals)
## ============================================================

use_vst_for_pca <- TRUE
log_base <- 2
pseudocount <- 1

n_top_var <- 5000
center_data <- TRUE
scale_data <- FALSE
label_samples <- FALSE
ellipse_level <- 0.95
point_size <- 2.0
point_alpha <- 0.9
show_scree <- TRUE

cor_method <- "pearson"
cluster_distance <- "correlation"
cluster_method <- "complete"
heat_lims <- c(0.8, 1.0)
show_corr_numbers <- FALSE
heatmap_width_in <- 8.5
heatmap_height_in <- 8.5

assay_names <- SummarizedExperiment::assayNames(exprMatrix)
has_vst <- "vstTransform" %in% assay_names
has_norm <- "normCounts" %in% assay_names
has_counts <- "counts" %in% assay_names

if (use_vst_for_pca && has_vst) {
  expr_mat <- SummarizedExperiment::assay(exprMatrix, "vstTransform")
  expr_label <- "VST"
} else if (has_norm) {
  nc <- SummarizedExperiment::assay(exprMatrix, "normCounts")
  expr_mat <- log(nc + pseudocount, base = log_base)
  expr_label <- glue("log{log_base}(normCounts + {pseudocount})")
} else if (has_counts) {
  counts <- SummarizedExperiment::assay(exprMatrix, "counts")
  lib_size <- colSums(counts, na.rm = TRUE)
  cpm <- sweep(counts, 2, lib_size, FUN = "/") * 1e6
  expr_mat <- log(cpm + pseudocount, base = log_base)
  expr_label <- glue("log{log_base}(CPM + {pseudocount})")
} else {
  stop("No suitable assay found for PCA/heatmap (need vstTransform or normCounts or counts).")
}

sample_ids <- colnames(exprMatrix)
cd <- as.data.frame(SummarizedExperiment::colData(exprMatrix))
stopifnot(group_col %in% colnames(cd))
group_vec <- factor(cd[[group_col]], levels = levels(cd[[group_col]]))
group_sym <- rlang::sym(group_col)

gene_vars <- matrixStats::rowVars(expr_mat, na.rm = TRUE)
if (is.na(n_top_var) || n_top_var <= 0 || n_top_var >= nrow(expr_mat)) {
  sel_idx <- seq_len(nrow(expr_mat))
} else {
  sel_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_top_var)]
}
expr_sel <- expr_mat[sel_idx, , drop = FALSE]

pca <- prcomp(t(expr_sel), center = center_data, scale. = scale_data)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as_tibble(pca$x, rownames = "sample_id") %>%
  left_join(tibble(sample_id = sample_ids, !!group_col := group_vec), by = "sample_id")

loadings <- as_tibble(pca$rotation, rownames = "gene_id")
topN <- 50
top_pc1 <- loadings %>% arrange(desc(abs(PC1))) %>% slice_head(n = topN)
top_pc2 <- loadings %>% arrange(desc(abs(PC2))) %>% slice_head(n = topN)

readr::write_csv(scores,  file.path(results_dir, "pca_scores.csv"))
readr::write_csv(top_pc1, file.path(results_dir, "pca_loadings_top_PC1.csv"))
readr::write_csv(top_pc2, file.path(results_dir, "pca_loadings_top_PC2.csv"))

p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = !!group_sym)) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_manual(values = time_colors, drop = FALSE) +
  labs(
    title = "PCA (PC1 vs PC2)",
    subtitle = glue("{expr_label} — top {length(sel_idx)} variable genes; center={center_data}, scale={scale_data}"),
    x = glue("PC1 ({round(100 * var_explained[1], 2)}% var)"),
    y = glue("PC2 ({round(100 * var_explained[2], 2)}% var)")
  ) +
  theme_pub()

p_pca <- p_pca + stat_ellipse(aes(group = !!group_sym), type = "norm",
                              level = ellipse_level, linewidth = 0.6, alpha = 0.4)

## Centroid labels (kept as in your script)
centroids_df <- scores %>%
  dplyr::group_by(!!group_sym) %>%
  dplyr::summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(group_val = !!group_sym)

p_pca <- p_pca +
  geom_point(
    data = centroids_df,
    aes(x = PC1, y = PC2, fill = group_val),
    shape = 21, size = 4.5, color = "black", stroke = 0.6, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = time_colors, drop = FALSE, guide = "none") +
  ggrepel::geom_label_repel(
    data = centroids_df,
    aes(x = PC1, y = PC2, label = as.character(group_val), color = group_val),
    size = 3.5, fontface = "bold", label.size = 0.15,
    box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf,
    inherit.aes = FALSE, show.legend = FALSE
  )

if (label_samples) {
  p_pca <- p_pca +
    ggrepel::geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 30)
}

p_scree <- ggplot(
  tibble(PC = factor(seq_along(var_explained)), pct = 100 * var_explained),
  aes(x = PC, y = pct, group = 1)
) +
  geom_col(fill = "grey70") +
  geom_point(size = 1.5) +
  geom_line() +
  labs(x = "Principal Component", y = "Explained variance (%)", title = "PCA Scree plot") +
  theme_pub()

save_plot(p_pca, "pca_scatter_pc1_pc2_18_11", width = 7.0, height = 5.6)
if (show_scree) save_plot(p_scree, "pca_scree_18_11", width = 6.0, height = 4.0)

## Correlation heatmap (ComplexHeatmap) + save to png/pdf/svg
cor_mat <- cor(expr_sel, method = cor_method, use = "pairwise.complete.obs")
rownames(cor_mat) <- colnames(cor_mat) <- sample_ids

ha_col <- HeatmapAnnotation(
  Group = group_vec,
  col = list(Group = structure(time_colors[levels(group_vec)], names = levels(group_vec))),
  annotation_legend_param = list(title = group_col)
)

col_fun <- circlize::colorRamp2(
  c(heat_lims[1], mean(heat_lims), heat_lims[2]),
  c("#313695", "white", "#A50026")
)

if (cluster_distance == "correlation") {
  d_rows <- as.dist(1 - cor_mat)
  d_cols <- d_rows
} else {
  d_rows <- dist(cor_mat, method = cluster_distance)
  d_cols <- dist(t(cor_mat), method = cluster_distance)
}
row_clust <- hclust(d_rows, method = cluster_method)
col_clust <- hclust(d_cols, method = cluster_method)

cell_fun <- NULL
if (isTRUE(show_corr_numbers)) {
  cell_fun <- function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = grid::gpar(fontsize = 7))
  }
}

ht <- Heatmap(
  cor_mat,
  name = glue("corr ({cor_method})"),
  col = col_fun,
  cluster_rows = row_clust,
  cluster_columns = col_clust,
  top_annotation = ha_col,
  show_row_names = FALSE,
  show_column_names = FALSE,
  rect_gp = grid::gpar(col = NA),
  cell_fun = cell_fun,
  heatmap_legend_param = list(
    title = glue("corr ({cor_method})"),
    at = seq(heat_lims[1], heat_lims[2], by = 0.05)
  )
)

save_ht <- function(ht, filename,
                    width = heatmap_width_in, height = heatmap_height_in,
                    units = "in", dpi = fig_dpi,
                    devices = c("png", "pdf", "svg")) {
  for (device in devices) {
    out_path <- file.path(fig_dir, sprintf("%s.%s", filename, device))
    if (device == "png") {
      png(filename = out_path, width = width, height = height, units = units, res = dpi)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
    } else if (device == "pdf") {
      pdf(file = out_path, width = width, height = height, onefile = FALSE)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
    } else if (device == "svg") {
      svglite::svglite(file = out_path, width = width, height = height)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
    } else {
      stop("Unsupported device: ", device)
    }
    message(glue("Saved: {out_path}"))
  }
  invisible(TRUE)
}

save_ht(ht, "sample_correlation_heatmap_18_11")

readr::write_csv(
  as_tibble(cor_mat, rownames = "sample_id"),
  file.path(results_dir, "sample_correlation_matrix.csv")
)

message(glue("PCA done on {nrow(expr_sel)} genes ({expr_label})."))
message("Saved: pca_scatter_pc1_pc2, pca_scree (if enabled) (png/pdf/svg)")
message("Saved: sample_correlation_heatmap (png/pdf/svg) + sample_correlation_matrix.csv")

## ===============================
## END (Blocks 1–5)
## ===============================
``

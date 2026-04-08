
## ============================================================
## DESeq2-only DGE script (no QC, no plots, no CSV)
## Input : exprMatrix.Rdata containing a SummarizedExperiment
## Output: RDS files under <base_path>/DESeq2_result/
## ============================================================

## -----------------------------
## 0) USER SETTINGS (edit these)
## -----------------------------
base_path <- "C:/Users/khosravs/OneDrive - Boehringer Ingelheim/Desktop/RNA PAPER/RNA data/Transcriptom_ Analysis/DESeq2_DGE"
out_dir_name <- "DESeq2_result"

# File name of your input object
expr_rdata_filename <- "exprMatrix.Rdata"

# Grouping variable in colData(exprMatrix)
group_col <- "MFGroup"
time_levels <- c("Day0", "Day1", "Day3", "Day7", "Day10")
ref_level <- "Day0"

# Paired/block design
use_pair_block <- TRUE
pair_block_col <- "EyeID"

# Filtering
filter_min_count <- 10
filter_min_n <- 3

# Significance threshold in results()
alpha_fdr <- 0.05

# lfcThreshold affects the hypothesis:
# 0 = standard DE (LFC != 0)
# 1 = thresholded DE (|LFC| > 1)
lfc_threshold_test <- 1

# Optional: shrink LFCs and save as RDS (no plots)
do_shrink <- TRUE
shrink_type <- "apeglm"   # recommended

set.seed(42)
options(stringsAsFactors = FALSE)

## -----------------------------
## 1) Package setup (robust)
## -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

req_bioc <- c("DESeq2", "SummarizedExperiment")
for (p in req_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

if (isTRUE(do_shrink) && !requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm", update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
})

## -----------------------------
## 2) Paths & output folder
## -----------------------------
out_dir <- file.path(base_path, out_dir_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Write sessionInfo (useful for troubleshooting)
sess_file <- file.path(out_dir, "sessionInfo.txt")
capture.output(utils::sessionInfo(), file = sess_file)

message("Output folder: ", out_dir)

## -----------------------------
## 3) Load exprMatrix.Rdata + input checks
## -----------------------------
expr_rdata_path <- file.path(base_path, expr_rdata_filename)

if (!file.exists(expr_rdata_path)) {
  stop("Cannot find input file: ", expr_rdata_path,
       "\nFix: check base_path and expr_rdata_filename.")
}

loaded_names <- load(expr_rdata_path)

# If 'exprMatrix' doesn't exist, try to detect SummarizedExperiment and assign
if (!exists("exprMatrix", envir = .GlobalEnv)) {
  
  if (length(loaded_names) == 1 && exists(loaded_names[1], envir = .GlobalEnv)) {
    obj <- get(loaded_names[1], envir = .GlobalEnv)
    
    if (is(obj, "SummarizedExperiment")) {
      exprMatrix <- obj
      message("Loaded object '", loaded_names[1], "' and assigned it to exprMatrix.")
    } else {
      stop("Loaded object '", loaded_names[1], "' is not a SummarizedExperiment.\nClass: ",
           paste(class(obj), collapse = ", "))
    }
    
  } else {
    se_hits <- loaded_names[sapply(loaded_names, function(nm) {
      exists(nm, envir = .GlobalEnv) &&
        is(get(nm, envir = .GlobalEnv), "SummarizedExperiment")
    })]
    
    if (length(se_hits) == 1) {
      exprMatrix <- get(se_hits[1], envir = .GlobalEnv)
      message("Detected SummarizedExperiment '", se_hits[1], "' and assigned it to exprMatrix.")
    } else if (length(se_hits) > 1) {
      stop("Multiple SummarizedExperiment objects found in the Rdata: ",
           paste(se_hits, collapse = ", "),
           "\nFix: rename the intended object to 'exprMatrix' before saving, or specify which to use.")
    } else {
      stop("No SummarizedExperiment object found in the loaded Rdata.\nObjects loaded: ",
           paste(loaded_names, collapse = ", "))
    }
  }
}

# Sanity check
stopifnot(is(exprMatrix, "SummarizedExperiment"))

assay_names <- SummarizedExperiment::assayNames(exprMatrix)
if (!("counts" %in% assay_names)) {
  stop("Assay 'counts' not found in exprMatrix.\nAvailable assays: ",
       paste(assay_names, collapse = ", "))
}

counts <- SummarizedExperiment::assay(exprMatrix, "counts")
cd <- as.data.frame(SummarizedExperiment::colData(exprMatrix))

# Alignment check
if (!identical(colnames(counts), rownames(cd))) {
  stop("Sample ID mismatch:\n  colnames(counts) != rownames(colData)\n",
       "Fix: ensure rownames(colData(exprMatrix)) exactly match colnames(assay(exprMatrix,'counts')).")
}

# Duplicate checks
if (any(duplicated(rownames(counts)))) stop("Duplicated gene IDs in rownames(counts). Resolve before DGE.")
if (any(duplicated(colnames(counts)))) stop("Duplicated sample IDs in colnames(counts). Resolve before DGE.")

message(sprintf("Loaded exprMatrix successfully: %s genes x %s samples.", nrow(counts), ncol(counts)))

## -----------------------------
## 4) Prepare metadata factors
## -----------------------------
if (!(group_col %in% colnames(cd))) {
  stop("group_col '", group_col, "' not found in colData. Available columns: ",
       paste(colnames(cd), collapse = ", "))
}

cd[[group_col]] <- factor(as.character(cd[[group_col]]), levels = time_levels)
if (any(is.na(cd[[group_col]]))) {
  bad <- rownames(cd)[is.na(cd[[group_col]])]
  stop("Some samples have group values not in time_levels. Offending samples: ",
       paste(bad, collapse = ", "))
}

cd[[group_col]] <- relevel(cd[[group_col]], ref = ref_level)

if (isTRUE(use_pair_block)) {
  if (!(pair_block_col %in% colnames(cd))) {
    stop("pair_block_col '", pair_block_col, "' not found in colData. Available columns: ",
         paste(colnames(cd), collapse = ", "))
  }
  cd[[pair_block_col]] <- factor(as.character(cd[[pair_block_col]]))
}

## -----------------------------
## 5) Filter genes
## -----------------------------
keep <- rowSums(counts >= filter_min_count, na.rm = TRUE) >= filter_min_n
message("Filtering kept ", sum(keep), " / ", nrow(counts), " genes.")
counts_f <- counts[keep, , drop = FALSE]

## -----------------------------
## 6) Design formula and DESeq fit (Wald)
## -----------------------------
if (isTRUE(use_pair_block)) {
  design_formula <- as.formula(paste0("~ ", pair_block_col, " + ", group_col))
} else {
  design_formula <- as.formula(paste0("~ ", group_col))
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData   = cd,
  design    = design_formula
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

## -----------------------------
## 7) Contrasts (baseline & adjacent) + results (no CSV)
## -----------------------------
groups <- levels(cd[[group_col]])

baseline_contrasts <- lapply(setdiff(groups, ref_level), function(g) c(group_col, g, ref_level))
baseline_names <- paste0(setdiff(groups, ref_level), "_vs_", ref_level)

adjacent_contrasts <- lapply(seq_len(length(groups) - 1), function(i) c(group_col, groups[i + 1], groups[i]))
adjacent_names <- paste0(groups[-1], "_vs_", groups[-length(groups)])

contrast_list <- c(baseline_contrasts, adjacent_contrasts)
names(contrast_list) <- c(baseline_names, adjacent_names)

# Helper to find coefficient name for shrinkage (avoids hard-coded names)
find_coef_name <- function(dds_obj, gcol, level, ref) {
  rn <- resultsNames(dds_obj)
  
  candidates <- c(
    paste0(gcol, "_", level, "_vs_", ref),
    paste0(gcol, level, "_vs_", ref),
    paste0(gcol, level, "vs", ref)
  )
  hit <- candidates[candidates %in% rn]
  if (length(hit) > 0) return(hit[1])
  
  pat <- paste0("^", gcol, ".*", level, ".*vs.*", ref, "$")
  hit2 <- grep(pat, rn, value = TRUE)
  if (length(hit2) > 0) return(hit2[1])
  
  return(NA_character_)
}

res_unshrunken <- vector("list", length(contrast_list))
names(res_unshrunken) <- names(contrast_list)

res_shrunken <- NULL
if (isTRUE(do_shrink)) {
  res_shrunken <- vector("list", length(contrast_list))
  names(res_shrunken) <- names(contrast_list)
}

for (nm in names(contrast_list)) {
  con <- contrast_list[[nm]]
  
  res_unshrunken[[nm]] <- results(
    dds,
    contrast     = con,
    alpha        = alpha_fdr,
    lfcThreshold = lfc_threshold_test
  )
  
  if (isTRUE(do_shrink)) {
    lvl <- con[2]; ref <- con[3]
    coef_nm <- find_coef_name(dds, con[1], lvl, ref)
    
    if (!is.na(coef_nm) && identical(shrink_type, "apeglm")) {
      res_shrunken[[nm]] <- lfcShrink(dds, coef = coef_nm, type = "apeglm")
    } else {
      warning("Shrinkage skipped for ", nm, " (coef not found or shrink_type unsupported).")
      res_shrunken[[nm]] <- NA
    }
  }
}

## -----------------------------
## 8) LRT (time effect) — PITFALL-FREE
## -----------------------------
# Reduced model:
# - if paired: ~ <pair_block_col>
# - if not paired: ~ 1
reduced_formula <- if (isTRUE(use_pair_block)) {
  as.formula(paste0("~ ", pair_block_col))
} else {
  ~ 1
}

dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_formula)
res_lrt <- results(dds_lrt, alpha = alpha_fdr)

## -----------------------------
## 9) Save outputs (RDS only)
## -----------------------------
saveRDS(dds, file.path(out_dir, "dds_wald_v2.rds"))
saveRDS(dds_lrt, file.path(out_dir, "dds_lrt_v2.rds"))

saveRDS(res_unshrunken, file.path(out_dir, "results_unshrunken_by_contrast_v2.rds"))
if (isTRUE(do_shrink)) {
  saveRDS(res_shrunken, file.path(out_dir, "results_shrunken_by_contrast_2.rds"))
}
saveRDS(res_lrt, file.path(out_dir, "results_LRT_time_effect_v2.rds"))

saveRDS(cd, file.path(out_dir, "colData_used_v2.rds"))
saveRDS(list(
  group_col = group_col,
  time_levels = time_levels,
  ref_level = ref_level,
  use_pair_block = use_pair_block,
  pair_block_col = pair_block_col,
  filter_min_count = filter_min_count,
  filter_min_n = filter_min_n,
  alpha_fdr = alpha_fdr,
  lfc_threshold_test = lfc_threshold_test,
  do_shrink = do_shrink,
  shrink_type = shrink_type,
  design = design_formula,
  reduced = reduced_formula,
  resultsNames = resultsNames(dds)
), file.path(out_dir, "run_settings_v2.rds"))

# Combined list object like your previous pipeline
all_DESeq2_results <- list(
  dds       = dds,
  dds_lrt   = dds_lrt,
  metadata  = cd,
  contrasts = contrast_list
)
saveRDS(all_DESeq2_results, file.path(out_dir, "all_DESeq2_results_v2.rds"))

message("DONE. Saved DESeq2 results to: ", out_dir)
message("Key file: ", file.path(out_dir, "all_DESeq2_results.rds"))


## ============================
## Sanity check: all_DESeq2_results.rds
## ============================

out_dir <- "C:/Users/khosravs/OneDrive - Boehringer Ingelheim/Desktop/RNA PAPER/RNA data/Transcriptom_ Analysis/DESeq2_DGE/DESeq2_result"
f <- file.path(out_dir, "all_DESeq2_results.rds")

stopifnot(file.exists(f))
x <- readRDS(f)

cat("\n[1] Top-level structure\n")
stopifnot(is.list(x))
req_names <- c("dds", "dds_lrt", "metadata", "contrasts")
missing <- setdiff(req_names, names(x))
if (length(missing) > 0) stop("Missing entries in all_DESeq2_results: ", paste(missing, collapse = ", "))
cat("OK: contains ->", paste(req_names, collapse = ", "), "\n")

dds <- x$dds
dds_lrt <- x$dds_lrt
md <- x$metadata
cons <- x$contrasts

cat("\n[2] Class checks\n")
cat("dds class     :", paste(class(dds), collapse = ", "), "\n")
cat("dds_lrt class :", paste(class(dds_lrt), collapse = ", "), "\n")
stopifnot(inherits(dds, "DESeqDataSet"))
stopifnot(inherits(dds_lrt, "DESeqDataSet"))
stopifnot(is.data.frame(md))
stopifnot(is.list(cons) && length(cons) > 0)

cat("\n[3] Dimensions & alignment\n")
cat("dds dims (genes x samples):", nrow(dds), "x", ncol(dds), "\n")
cat("metadata rows:", nrow(md), "\n")

# colData alignment
stopifnot(identical(colnames(counts(dds)), rownames(colData(dds))))
# metadata rownames should match sample IDs
if (!identical(rownames(md), rownames(colData(dds)))) {
  warning("metadata rownames != dds colData rownames (sample IDs). Downstream merges by rowname may fail.")
} else {
  cat("OK: metadata rownames match dds sample IDs\n")
}

cat("\n[4] Model/design sanity\n")
cat("Design formula (dds):\n")
print(design(dds))

cat("\nresultsNames(dds):\n")
print(resultsNames(dds))

cat("\n[5] Contrasts sanity\n")
stopifnot(is.character(names(cons)))
stopifnot(all(sapply(cons, length) == 3))

# Confirm contrast variable exists in colData
contrast_vars <- unique(sapply(cons, `[`, 1))
if (length(contrast_vars) != 1) {
  warning("More than one contrast variable detected: ", paste(contrast_vars, collapse = ", "))
} else {
  v <- contrast_vars[1]
  if (!(v %in% colnames(colData(dds)))) stop("Contrast variable '", v, "' not found in colData(dds).")
  lvls <- levels(colData(dds)[[v]])
  cat("Contrast variable:", v, "\n")
  cat("Levels:", paste(lvls, collapse = ", "), "\n")
  
  # Validate each contrast uses known levels
  bad_cons <- names(cons)[!sapply(cons, function(cc) all(cc[2:3] %in% lvls))]
  if (length(bad_cons) > 0) {
    stop("Some contrasts reference levels not in colData levels: ", paste(bad_cons, collapse = ", "))
  } else {
    cat("OK: All contrasts reference valid levels.\n")
  }
}

cat("\n[6] Check that differential testing produced non-empty output\n")
# Spot-check 1–2 contrasts that should exist (first two)
nm_check <- names(cons)[seq_len(min(2, length(cons)))]
cat("Will spot-check contrasts:", paste(nm_check, collapse = ", "), "\n")

# Load result lists saved by script (recommended sanity check)
res_un_path <- file.path(out_dir, "results_unshrunken_by_contrast_v2.rds")
stopifnot(file.exists(res_un_path))
res_un <- readRDS(res_un_path)

stopifnot(is.list(res_un))
stopifnot(all(nm_check %in% names(res_un)))

for (nm in nm_check) {
  r <- res_un[[nm]]
  cat("\n--", nm, "--\n")
  cat("Class:", paste(class(r), collapse = ", "), "\n")
  cat("Rows:", nrow(r), "\n")
  cat("Non-NA pvalues:", sum(!is.na(r$pvalue)), "\n")
  cat("Non-NA padj   :", sum(!is.na(r$padj)), "\n")
  if (sum(!is.na(r$padj)) > 0) {
    cat("Significant (padj<=0.05):", sum(r$padj <= 0.05, na.rm = TRUE), "\n")
  }
  # Basic numeric sanity
  if (all(is.na(r$log2FoldChange))) warning("All log2FoldChange are NA for contrast ", nm)
}

cat("\n[7] LRT sanity\n")
res_lrt_path <- file.path(out_dir, "results_LRT_time_effect_v2.rds")
stopifnot(file.exists(res_lrt_path))
rL <- readRDS(res_lrt_path)

cat("LRT results rows:", nrow(rL), "\n")
cat("Non-NA padj:", sum(!is.na(rL$padj)), "\n")
cat("Significant (padj<=0.05):", sum(rL$padj <= 0.05, na.rm = TRUE), "\n")

cat("\n[8] Quick dispersion/size-factor checks\n")
sf <- sizeFactors(dds)
cat("Size factors summary:\n")
print(summary(sf))
if (any(!is.finite(sf)) || any(sf <= 0)) warning("Non-finite or non-positive size factors detected!")

cat("\n[9] Save quick snapshot for downstream reproducibility\n")
# Not required, but useful:
qc_snapshot <- list(
  out_dir = out_dir,
  n_genes = nrow(dds),
  n_samples = ncol(dds),
  sample_ids = rownames(colData(dds)),
  design = design(dds),
  resultsNames = resultsNames(dds),
  contrast_names = names(cons)
)
saveRDS(qc_snapshot, file.path(out_dir, "sanity_snapshot_v2.rds"))
cat("Saved sanity_snapshot.rds\n")



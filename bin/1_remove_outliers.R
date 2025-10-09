#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 1_remove_outliers.R
# Project: Bacillus subtilis QTL mapping
# Purpose: Remove outliers from continuous phenotypic traits
#          using Median Absolute Deviation (MAD). Outliers -> NA.
#
# Input TSV must contain at least:
#   - progeny_strain            (genotype / strain ID)
#   - one or more trait columns (e.g., spore_optical_density, spore_area)
#   - any other columns are preserved unchanged (e.g., measurement_date, Batch)
#
# Example (auto-detect numeric traits except 'progeny_strain'):
#   Rscript bin/1_remove_outliers.R \
#     --input data/progeny_phenotypes.tsv \
#     --output results/progeny_phenotypes_clean.tsv \
#     --mad_cutoff 6
#
# Example (explicit traits):
#   Rscript bin/1_remove_outliers.R \
#     --input data/progeny_phenotypes.tsv \
#     --output results/progeny_phenotypes_clean.tsv \
#     --mad_cutoff 6 \
#     --traits spore_optical_density,spore_area
#
# Author: Hari B. Chhetri
# License: MIT
# --------------------------------------------------------------------

suppressPackageStartupMessages({ library(optparse) })

# ---- CLI options ----
opt_list <- list(
  make_option("--input",  type="character", help="Input TSV file (required)"),
  make_option("--output", type="character", default="progeny_phenotypes_clean.tsv",
              help="Output TSV file [default: %default]"),
  make_option("--mad_cutoff", type="double", default=6,
              help="MAD cutoff (|x - median| / MAD > cutoff) [default: %default]"),
  make_option("--traits", type="character", default="",
              help="Comma-separated trait column names to clean. If omitted, all numeric columns except 'progeny_strain' are cleaned.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input) || !file.exists(opt$input))
  stop("Please provide --input <file.tsv> (file must exist).")

# ---- Helper: MAD-based outlier removal ----
remove_outliers_mad <- function(x, cutoff = 6, b = 1.4826) {
  if (!is.numeric(x)) return(x)
  med  <- median(x, na.rm = TRUE)
  madv <- mad(x, constant = b, na.rm = TRUE)
  if (is.na(madv) || madv == 0) return(x)  # no variability; nothing to flag
  z <- abs((x - med) / madv)
  x[z > cutoff] <- NA
  round(x, 4)
}

# ---- Read data ----
d <- read.table(opt$input, header = TRUE, sep = "\t",
                check.names = FALSE, stringsAsFactors = FALSE)

# ---- Validate columns ----
if (!"progeny_strain" %in% names(d))
  stop("Missing required column: 'progeny_strain'")

# ---- Determine traits to clean ----
if (nzchar(opt$traits)) {
  trait_cols <- trimws(strsplit(opt$traits, ",")[[1]])
  missing <- setdiff(trait_cols, names(d))
  if (length(missing))
    stop("Trait column(s) not found: ", paste(missing, collapse = ", "))
} else {
  # auto-detect numeric columns, excluding the ID column
  numeric_cols <- names(d)[sapply(d, is.numeric)]
  trait_cols <- setdiff(numeric_cols, "progeny_strain")
}

if (length(trait_cols) == 0)
  stop("No trait columns to process. Pass --traits or include numeric trait columns in the TSV.")

cat("Trait columns to clean (MAD): ", paste(trait_cols, collapse = ", "), "\n", sep="")

# ---- Clean each trait & track summary ----
summary_counts <- data.frame(
  Trait = trait_cols,
  Total = NA_integer_, NA_before = NA_integer_, NA_after = NA_integer_,
  Outliers_flagged = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(trait_cols)) {
  cn <- trait_cols[i]
  if (!is.numeric(d[[cn]]))
    d[[cn]] <- suppressWarnings(as.numeric(d[[cn]]))

  n_total   <- length(d[[cn]])
  na_before <- sum(is.na(d[[cn]]))
  cleaned   <- remove_outliers_mad(d[[cn]], cutoff = opt$mad_cutoff)
  na_after  <- sum(is.na(cleaned))
  flagged   <- max(na_after - na_before, 0)

  d[[cn]] <- cleaned

  summary_counts$Total[i]            <- n_total
  summary_counts$NA_before[i]        <- na_before
  summary_counts$NA_after[i]         <- na_after
  summary_counts$Outliers_flagged[i] <- flagged
}

# ---- Write cleaned TSV ----
write.table(d, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved cleaned file: %s\n", opt$output))

# ---- Print summary ----
cat("\nSummary (per trait):\n")
print(summary_counts, row.names = FALSE)
cat("\n")


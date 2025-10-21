#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 1_remove_outliers.R
# Purpose: MAD-based outlier removal for continuous trait columns.
# Input TSV must contain: progeny_strain; may also contain: batch (lowercase).
# Other columns are preserved; numeric traits are auto-detected unless specified.
#
# Usage:
#   Rscript bin/1_remove_outliers.R \
#     --input data/progeny_phenotypes.tsv \
#     --output results/progeny_phenotypes_clean.tsv \
#     --mad_cutoff 6
#
# Options:
#   --traits <comma,separated,names>   # optional explicit trait list
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--input",  type="character", help="Input TSV (required)"),
  make_option("--output", type="character", help="Output TSV (required)"),
  make_option("--mad_cutoff", type="double", default=6,
              help="MAD cutoff (abs devs > cutoff become NA) [default: %default]"),
  make_option("--traits", type="character", default="",
              help="Comma-separated trait columns. If omitted, use all numeric columns except 'progeny_strain' and 'batch'.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input)  || !file.exists(opt$input))  stop("Provide --input <file.tsv> that exists.")
if (is.null(opt$output)) stop("Provide --output <file.tsv>.")

# ---------- IO (base R) ----------
dt <- read.delim(opt$input, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

if (!"progeny_strain" %in% names(dt)) stop("Missing required column 'progeny_strain'.")
has_batch <- "batch" %in% names(dt)

# ---------- trait selection ----------
if (nzchar(opt$traits)) {
  trait_cols <- trimws(strsplit(opt$traits, ",")[[1]])
  missing <- setdiff(trait_cols, names(dt))
  if (length(missing)) stop("Trait column(s) not found: ", paste(missing, collapse=", "))
} else {
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  trait_cols <- setdiff(numeric_cols, c("progeny_strain", if (has_batch) "batch"))
}
if (!length(trait_cols)) stop("No numeric trait columns found to clean.")

cat("Trait columns to clean (MAD): ", paste(trait_cols, collapse=", "), "\n", sep="")

# ---------- MAD outlier function ----------
outliersMAD <- function(x, MADcutoff=6, replace=NA_real_, bConstant=1.4826) {
  x_num <- suppressWarnings(as.numeric(x))
  # leave non-numeric entirely untouched
  if (all(is.na(x_num)) && !all(is.na(x))) return(x)
  med <- stats::median(x_num, na.rm=TRUE)
  madv <- stats::mad(x_num, constant=bConstant, na.rm=TRUE)
  if (!is.finite(madv) || madv == 0) return(round(x_num, 6))  # nothing to do
  absMADAway <- abs((x_num - med) / madv)
  x_num[absMADAway > MADcutoff] <- replace
  round(x_num, 6)
}

# ---------- clean & summarize ----------
summary_rows <- list()

for (tr in trait_cols) {
  orig <- dt[[tr]]
  if (!is.numeric(orig)) orig <- suppressWarnings(as.numeric(orig))
  na_before <- sum(is.na(orig))
  cleaned <- outliersMAD(orig, MADcutoff = opt$mad_cutoff)
  na_after  <- sum(is.na(cleaned))
  flagged   <- (na_after - na_before)
  dt[[tr]]  <- cleaned

  summary_rows[[tr]] <- data.frame(
    Trait = tr,
    Total = length(orig),
    NA_before = na_before,
    NA_after  = na_after,
    Outliers_flagged = max(flagged, 0),
    stringsAsFactors = FALSE
  )
}

# ---------- write output ----------
dir.create(dirname(opt$output), showWarnings=FALSE, recursive=TRUE)
write.table(dt, opt$output, sep="\t", quote=FALSE, row.names=FALSE)
cat("Saved cleaned file: ", opt$output, "\n", sep="")

# ---------- print summary ----------
smry <- do.call(rbind, summary_rows)
print(smry, row.names=FALSE)


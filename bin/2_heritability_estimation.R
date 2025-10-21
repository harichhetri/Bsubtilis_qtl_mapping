#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 2_heritability_estimation.R  (nlme version, stable under Rscript)
# Purpose: Estimate variance components and heritability (H2) per trait
#          using random intercept(s): progeny_strain (always) and batch
#          (optional if present). No BLUPs here.
#
# Usage:
#   Rscript bin/2_heritability_estimation.R INPUT_TSV OUTPREFIX [TRAITS_CSV]
#   - INPUT_TSV  : e.g., results/progeny_phenotypes_clean.tsv
#   - OUTPREFIX  : e.g., results/heritability
#   - TRAITS_CSV : optional comma-separated trait names; if omitted, all
#                  numeric columns except 'progeny_strain' and 'batch' are used.
#
# Output:
#   <OUTPREFIX>_variance.tsv with columns:
#     trait, N, var_G, var_batch, var_E, H2
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 2_heritability_estimation.R INPUT_TSV OUTPREFIX [TRAITS_CSV]")
}
input_tsv  <- args[1]
outprefix  <- args[2]
traits_csv <- if (length(args) >= 3) args[3] else ""

if (!file.exists(input_tsv)) stop("Input file not found: ", input_tsv)

# --- Read (base R) ---
dt <- read.delim(input_tsv, header = TRUE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)

if (!"progeny_strain" %in% names(dt)) stop("Missing required column 'progeny_strain'.")
has_batch <- "batch" %in% names(dt)

# --- Determine trait columns ---
if (nzchar(traits_csv)) {
  trait_cols <- trimws(strsplit(traits_csv, ",")[[1]])
  missing <- setdiff(trait_cols, names(dt))
  if (length(missing)) stop("Trait column(s) not found: ", paste(missing, collapse = ", "))
} else {
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  trait_cols <- setdiff(numeric_cols, c("progeny_strain", if (has_batch) "batch"))
}
if (!length(trait_cols)) stop("No trait columns to model.")
cat("Traits to model:", paste(trait_cols, collapse = ", "), "\n")
if (has_batch) cat("Including 'batch' as a random effect.\n")

# --- Types ---
dt$progeny_strain <- as.factor(dt$progeny_strain)
if (has_batch && !is.factor(dt$batch)) dt$batch <- as.factor(dt$batch)

# --- Fit per trait (nlme-only) ---
min_n  <- 20
digits <- 6
rows   <- list()

for (tr in trait_cols) {
  y <- dt[[tr]]
  if (!is.numeric(y)) y <- suppressWarnings(as.numeric(y))
  n_obs <- sum(!is.na(y))
  if (n_obs < min_n) {
    warning(sprintf("Skipping '%s': only %d non-missing (< %d).", tr, n_obs, min_n))
    next
  }

  if (has_batch) {
    sub <- data.frame(trait = y, progeny_strain = dt$progeny_strain, batch = dt$batch)
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) < min_n) {
      warning(sprintf("Skipping '%s': insufficient complete cases after filtering.", tr))
      next
    }
    fit <- try(nlme::lme(trait ~ 1,
                         random = list(progeny_strain = ~1, batch = ~1),
                         data = sub, method = "REML",
                         control = nlme::lmeControl(msMaxIter = 200, opt = "optim")),
               silent = TRUE)
  } else {
    sub <- data.frame(trait = y, progeny_strain = dt$progeny_strain)
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) < min_n) {
      warning(sprintf("Skipping '%s': insufficient complete cases after filtering.", tr))
      next
    }
    fit <- try(nlme::lme(trait ~ 1,
                         random = ~1 | progeny_strain,
                         data = sub, method = "REML",
                         control = nlme::lmeControl(msMaxIter = 200, opt = "optim")),
               silent = TRUE)
  }

  if (inherits(fit, "try-error")) {
    warning(sprintf("Model failed for '%s': %s", tr, as.character(fit)))
    next
  }

  # Variance components from nlme::VarCorr (StdDev → square to get variance)
  vc <- nlme::VarCorr(fit)
  sd_E <- as.numeric(vc[nrow(vc), "StdDev"]); var_E <- sd_E^2

  var_G <- NA_real_
  idx_G <- which(grepl("^progeny_strain", rownames(vc)))
  if (length(idx_G)) {
    sd_G <- as.numeric(vc[idx_G[1], "StdDev"]); var_G <- sd_G^2
  }

  var_batch <- 0
  if (has_batch) {
    idx_B <- which(grepl("^batch", rownames(vc)))
    if (length(idx_B)) {
      sd_B <- as.numeric(vc[idx_B[1], "StdDev"]); var_batch <- sd_B^2
    }
  }

  denom <- var_G + var_batch + var_E
  H2 <- if (is.finite(denom) && denom > 0) var_G / denom else NA_real_

  rows[[tr]] <- data.frame(
    trait     = tr,
    N         = nrow(sub),
    var_G     = round(var_G,     digits),
    var_batch = round(var_batch, digits),
    var_E     = round(var_E,     digits),
    H2        = round(H2,        digits),
    stringsAsFactors = FALSE
  )
}

out_var <- paste0(outprefix, "_variance.tsv")
res <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.table(res, out_var, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Variance components written:", out_var, "\n")


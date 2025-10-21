#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 3_random_effect_model.R
# Purpose: Estimate BLUPs and mean-adjusted BLUPs per trait, writing
#          a separate file for each trait.
#
# Model:
#   trait ~ 1 + (1|progeny_strain) [+ (1|batch) if present]
#
# Inputs (TSV):
#   - required: progeny_strain
#   - numeric trait columns (e.g., spore_area, optical_density)
#   - optional: batch (lowercase; used as random intercept if present)
#
# Output files (per trait):
#   If --outprefix ends with "/" or is a directory:
#       <outprefix>/<trait>_blups.tsv
#   Else if basename(outprefix) == "blups" (case-insensitive):
#       <dirname(outprefix)>/blups_<trait>.tsv
#   Else:
#       <dirname(outprefix)>/<basename>_<trait>_blups.tsv
#
# Example:
#   Rscript --vanilla bin/3_random_effect_model.R \
#     --input results/progeny_phenotypes_clean.tsv \
#     --outprefix results/blups/ \
#     --traits spore_area,optical_density
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
})

# ---- CLI ----
opt_list <- list(
  make_option("--input",     type="character", help="Input cleaned TSV (from script 1) [required]"),
  make_option("--outprefix", type="character", default="results/blups",
              help="Output prefix or directory [default: %default]"),
  make_option("--traits",    type="character", default="",
              help="Comma-separated trait columns; if omitted, use all numeric columns except 'progeny_strain' and 'batch'"),
  make_option("--min_n",     type="integer", default=10,
              help="Minimum non-missing observations per trait to fit [default: %default]"),
  make_option("--digits",    type="integer", default=6,
              help="Rounding for outputs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("Please provide --input <file.tsv> that exists.")
}

# ---- read (base R) ----
dt <- read.delim(opt$input, header = TRUE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)
if (!"progeny_strain" %in% names(dt)) stop("Missing required column 'progeny_strain'.")
has_batch <- "batch" %in% names(dt)

# ---- choose traits ----
if (nzchar(opt$traits)) {
  trait_cols <- trimws(strsplit(opt$traits, ",")[[1]])
  missing <- setdiff(trait_cols, names(dt))
  if (length(missing)) stop("Trait column(s) not found: ", paste(missing, collapse = ", "))
} else {
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  trait_cols <- setdiff(numeric_cols, c("progeny_strain", if (has_batch) "batch"))
}
if (!length(trait_cols)) stop("No trait columns to model.")

cat("Traits to model: ", paste(trait_cols, collapse = ", "), "\n", sep = "")
if (has_batch) cat("Including 'batch' as a random effect.\n")

# ---- ensure factors (+ trim whitespace) ----
dt$progeny_strain <- as.factor(trimws(dt$progeny_strain))
if (has_batch) dt$batch <- as.factor(trimws(dt$batch))

# ---- helpers ----
extract_blups <- function(fit) {
  # Standard path
  re_list <- try(ranef(fit, drop = TRUE), silent = TRUE)
  if (!inherits(re_list, "try-error")) {
    u <- re_list[["progeny_strain"]]
    if (!is.null(u) && is.data.frame(u) && nrow(u) > 0) {
      return(u[, 1, drop = FALSE])  # intercept BLUPs
    }
  }
  # Fallback via coef() - fixef()
  cf <- try(coef(fit)$progeny_strain, silent = TRUE)
  fx <- try(fixef(fit)[["(Intercept)"]], silent = TRUE)
  if (!inherits(cf, "try-error") && !is.null(cf) &&
      !inherits(fx, "try-error") && is.finite(fx)) {
    v <- cf[, "(Intercept)"] - fx
    u <- data.frame(`(Intercept)` = v, check.names = FALSE)
    rownames(u) <- rownames(cf)
    return(u)
  }
  return(NULL)
}

sanitize <- function(x) {
  x <- gsub("[^[:alnum:]]+", "_", x)
  gsub("_+$", "", gsub("^_+", "", x))
}

split_outprefix <- function(outprefix) {
  if (grepl("/$", outprefix) || dir.exists(outprefix)) {
    out_dir <- if (grepl("/$", outprefix)) outprefix else paste0(outprefix, "/")
    base <- ""
  } else {
    out_dir <- if (dirname(outprefix) %in% c(".", "")) "" else paste0(dirname(outprefix), "/")
    base <- basename(outprefix)
  }
  list(dir = out_dir, base = base)
}

# ---- fit per trait and write per-trait file ----
for (tr in trait_cols) {
  y <- dt[[tr]]
  if (!is.numeric(y)) y <- suppressWarnings(as.numeric(y))
  n_obs <- sum(!is.na(y))
  if (n_obs < opt$min_n) {
    warning(sprintf("Skipping '%s': only %d non-missing (< %d).", tr, n_obs, opt$min_n))
    next
  }

  if (has_batch) {
    sub <- data.frame(trait = y, progeny_strain = dt$progeny_strain, batch = dt$batch)
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) < opt$min_n) {
      warning(sprintf("Skipping '%s': insufficient complete cases after filtering.", tr))
      next
    }
    fml <- as.formula("trait ~ 1 + (1|progeny_strain) + (1|batch)")
  } else {
    sub <- data.frame(trait = y, progeny_strain = dt$progeny_strain)
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) < opt$min_n) {
      warning(sprintf("Skipping '%s': insufficient complete cases after filtering.", tr))
      next
    }
    fml <- as.formula("trait ~ 1 + (1|progeny_strain)")
  }

  fit <- try(lmer(fml, data = sub, REML = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning(sprintf("Model failed for '%s': %s", tr, as.character(fit)))
    next
  }

  u <- extract_blups(fit)
  if (is.null(u)) {
    warning(sprintf("Could not extract BLUPs for '%s'.", tr))
    next
  }

  overall_mean <- mean(y, na.rm = TRUE)
  digits <- opt$digits

  out <- data.frame(
    progeny_strain = rownames(u),
    trait          = tr,
    BLUP           = round(u[, 1], digits),
    OverallMean    = round(overall_mean, digits),
    MeanAdj_BLUP   = round(overall_mean + u[, 1], digits),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Build clean per-trait filename
  trait_safe <- sanitize(tr)
  parts <- split_outprefix(opt$outprefix)
  out_dir <- parts$dir
  base    <- parts$base
  if (nzchar(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!nzchar(base)) {
    out_file <- paste0(out_dir, trait_safe, "_blups.tsv")
  } else if (tolower(base) == "blups") {
    out_file <- paste0(out_dir, "blups_", trait_safe, ".tsv")
  } else {
    out_file <- paste0(out_dir, base, "_", trait_safe, "_blups.tsv")
  }

  write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("BLUPs written: ", out_file, "\n", sep = "")
}


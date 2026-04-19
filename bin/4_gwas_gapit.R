#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 4_gwas_gapit.R
# Purpose: Run GWAS using GAPIT (MLM) on phenotype + genotype data
#
# Usage:
#   Rscript 4_gwas_gapit.R \
#     --pheno phenotypes.tsv \
#     --geno genotype.hmp.txt \
#     --outdir gapit_results \
#     [--trait trait_name]  # required only for multi-column phenotype
#
# Inputs:
#   --pheno : phenotype file
#             (either GAPIT format: Taxa|Trait OR wide format)
#   --geno  : HapMap genotype file
#   --trait : required if phenotype has multiple traits
#
# Output:
#   GAPIT output files written to --outdir
# --------------------------------------------------------------------

# ---- Fix CRAN mirror for CLI ----
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# -----------------------------
# CLI options
# -----------------------------
opt_list <- list(
  make_option("--pheno", type="character", help="Phenotype file (required)"),
  make_option("--geno",  type="character", help="Genotype HapMap file (required)"),
  make_option("--trait", type="character", default=NULL,
              help="Trait column (required for multi-column phenotype)"),
  make_option("--outdir", type="character", default="gapit_results",
              help="Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$pheno) || !file.exists(opt$pheno)) {
  stop("Provide --pheno <file> that exists.")
}
if (is.null(opt$geno) || !file.exists(opt$geno)) {
  stop("Provide --geno <file> that exists.")
}

# -----------------------------
# Load GAPIT
# -----------------------------
cat("Loading GAPIT...\n")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# -----------------------------
# Read phenotype
# -----------------------------
cat("Reading phenotype...\n")
df <- read.table(opt$pheno, header = TRUE, check.names = FALSE)

# ---- AUTO-DETECT FORMAT ----
if (ncol(df) == 2 && all(c("Taxa", "Trait") %in% names(df))) {
  
  cat("Detected GAPIT-ready phenotype format (Taxa | Trait)\n")
  myY <- df
  
} else {
  
  cat("Detected multi-column phenotype format\n")
  
  if (is.null(opt$trait)) {
    stop("Provide --trait for multi-column phenotype file.")
  }
  
  if (!("progeny_strain" %in% names(df))) {
    stop("Column 'progeny_strain' not found.")
  }
  
  if (!(opt$trait %in% names(df))) {
    stop("Trait not found: ", opt$trait)
  }
  
  cat("Using trait:", opt$trait, "\n")
  
  myY <- df[, c("progeny_strain", opt$trait)]
  colnames(myY) <- c("Taxa", "Trait")
}

# ---- Clean phenotype ----
myY[,2] <- suppressWarnings(as.numeric(myY[,2]))
myY <- myY[complete.cases(myY), ]

if (nrow(myY) == 0) {
  stop("No valid phenotype values after filtering.")
}

cat("Using", nrow(myY), "samples for GWAS\n")

# -----------------------------
# Read genotype
# -----------------------------
cat("Reading genotype...\n")
myG <- fread(
  opt$geno,
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

# Fix marker IDs (keep original behavior)
myG[,1] <- c("rs", 1:(nrow(myG)-1))

# -----------------------------
# Prepare output directory
# -----------------------------
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$outdir)

# -----------------------------
# Run GAPIT
# -----------------------------
cat("Running GAPIT MLM...\n")

myGAPIT_MLM <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 0,
  Random.model = FALSE,
  SNP.MAF = 0.03,
  SNP.impute = "Major",
  Major.allele.zero = TRUE,
  Multiple_analysis = FALSE,
  model = "MLM"
)

cat("GAPIT analysis completed.\n")
#!/usr/bin/env Rscript
# --------------------------------------------------------------------
# Script: 5_plot_variant_map.R
# Purpose: Plot SNP presence/absence heatmap from VCF file
#
# Usage:
#   Rscript 5_plot_variant_map.R \
#     --vcf input.vcf \
#     --out b_subtilis_variant_map.png
#
# Inputs:
#   --vcf : input VCF file (required)
#   --out : output PNG file [default: variant_map.png]
#
# Output:
#   PNG heatmap of variant presence/absence
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(vcfR)
})

# -----------------------------
# CLI options
# -----------------------------
opt_list <- list(
  make_option("--vcf", type="character", help="Input VCF file (required)"),
  make_option("--out", type="character", default="variant_map.png",
              help="Output PNG file [default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$vcf) || !file.exists(opt$vcf)) {
  stop("Provide --vcf <file.vcf> that exists.")
}

# -----------------------------
# READ VCF
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(opt$vcf)

# -----------------------------
# EXTRACT GENOTYPE
# -----------------------------
cat("Extracting genotype matrix...\n")
gt <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# convert to binary presence/absence
gt_bin <- ifelse(gt > 0, 1, 0)
gt_bin[is.na(gt)] <- NA

# transpose: strains x SNPs
mat <- t(gt_bin)

cat("Matrix dimensions:", nrow(mat), "strains x", ncol(mat), "SNPs\n")

# -----------------------------
# PLOT
# -----------------------------
cat("Generating plot...\n")

png(
  filename = opt$out,
  width = 2400,
  height = 1200,
  res = 200
)

par(mar = c(3, 3, 2, 2))

cols <- colorRampPalette(c("gray20", "dodgerblue"))(200)

image(
  t(mat[nrow(mat):1, ]),
  col = cols,
  axes = FALSE
)

# -----------------------------
# LEGEND
# -----------------------------
legend(
  "topright",
  inset = 0.01,
  legend = c(
    paste0("Strains: 1 (top) → ", nrow(mat), " (bottom)"),
    paste0("SNPs: 1 (left) → ", round(ncol(mat)/1000), "K (right)")
  ),
  bty = "n",
  cex = 0.9,
  text.col = "white"
)

dev.off()

cat("Saved plot:", opt$out, "\n")
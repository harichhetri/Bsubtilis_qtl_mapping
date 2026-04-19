# Bacillus subtilis QTL Mapping

This repository contains the example analysis code used for association mapping of spore germination–related traits in *Bacillus subtilis*, including:

- Phenotype QC (MAD-distance based outlier removal)
- Mixed-model variance components and heritability estimation
- BLUP estimation (random effects)
- QTL mapping using GAPIT MLM
- Variant map visualization from VCF data

---

## 🚀 Quick start

### Set up the environment

```bash
# Create environment
conda env create -f envs/environment.yml

# Activate (conda or micromamba)
conda activate bsubtilis-qtl-r44
# OR
micromamba activate bsubtilis-qtl-r44

## Run the full pipeline
Run the following commands from the repository root directory (`Bsubtilis_qtl_mapping`):
```
```bash
Rscript 1_remove_outliers.R \
  --input phenotypes_raw.tsv \
  --output phenotypes_outlier_removed.tsv \
  --mad_cutoff 6

Rscript 2_heritability_estimation.R \
  phenotypes_outlier_removed.tsv \
  heritability

Rscript 3_random_effect_model.R \
  --input phenotypes_outlier_removed.tsv \
  --outprefix blups \
  --traits spore_area,optical_density

Rscript 4_gwas_gapit.R \
  --pheno pheno_data_qtl_mapping.txt \
  --geno genotype_data_qtl_mapping_hmp.txt \
  --outdir gapit_results

Rscript 5_plot_variant_map.R \
  --vcf genomic_data_variant_map_10k_snps_200_strains.vcf \
  --out variant_map.png

## 📁 Input files

All required input files are included in this repository and can be used to run the pipeline directly.

- phenotypes_raw.tsv  
  Raw phenotype data  

- phenotypes_outlier_removed.tsv  
  Cleaned phenotype data after outlier removal  

- pheno_data_qtl_mapping.txt  
  GAPIT-ready phenotype file (Taxa, Trait format)  

- genotype_data_qtl_mapping_hmp.txt  
  HapMap genotype file used for GWAS  

- genomic_data_variant_map_10k_snps_200_strains.vcf  
  VCF file used for variant visualization  
  
```bash
# Create environment
conda env create -f envs/environment.yml

# Activate (conda or micromamba)
conda activate bsubtilis-qtl-r44
# OR
micromamba activate bsubtilis-qtl-r44
```
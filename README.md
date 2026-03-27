# Whole-Genome Sequencing Reveals Genomic Diversity and Population Structure of Saudi Arabian Mango Germplasm

**Species:** *Mangifera indica* L.
**Dataset:** 60 accessions | ~840,000 LD-pruned SNPs | 20 chromosomes
**Sampling site:** Jazan region, Saudi Arabia

> Raw sequencing data: NCBI SRA [Accession numbers to be provided upon acceptance]

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Repository Structure](#repository-structure)
4. [Pipeline](#pipeline)
   - [Step 1 — Variant Calling & Filtering (GATK)](#step-1--variant-calling--filtering-gatk)
   - [Step 2 — Sample QC: Missingness](#step-2--sample-qc-missingness)
   - [Step 3 — Sample QC: Heterozygosity](#step-3--sample-qc-heterozygosity)
   - [Step 4 — LD Pruning](#step-4--ld-pruning)
   - [Step 5 — Principal Component Analysis (PCA)](#step-5--principal-component-analysis-pca)
   - [Step 6 — Ancestry Inference (ADMIXTURE)](#step-6--ancestry-inference-admixture)
   - [Step 7 — Kinship Estimation (KING)](#step-7--kinship-estimation-king)
   - [Step 8 — Identity by Descent (IBD)](#step-8--identity-by-descent-ibd)
   - [Step 9 — Pairwise IBS Tree](#step-9--pairwise-ibs-tree)
   - [Step 10 — Phylogenetic Tree (IQ-TREE)](#step-10--phylogenetic-tree-iq-tree)
   - [Step 11 — SNP Fingerprinting Panel](#step-11--snp-fingerprinting-panel)
   - [Step 12 — Population Statistics (Ho, He, FST)](#step-12--population-statistics-ho-he-fst)
   - [Step 13 — Nucleotide Diversity (π)](#step-13--nucleotide-diversity-π)
   - [Step 14 — LD Decay](#step-14--ld-decay)
   - [Step 15 — Sampling Map](#step-15--sampling-map)
5. [Key Results](#key-results)
6. [Citation](#citation)

---

## Overview

This repository contains the complete analysis pipeline for population genomics of 60 mango (*Mangifera indica*) accessions collected from the Jazan region, Saudi Arabia. The pipeline covers SNP quality control, population structure analysis (PCA, ADMIXTURE, phylogenetics), relatedness estimation, and the development of a 150-SNP fingerprinting panel for cultivar identification.

---

## Requirements

### External Tools

| Tool | Version | Purpose |
|------|---------|---------|
| [GATK](https://gatk.broadinstitute.org) | ≥ 4.0 | Variant calling & filtering |
| [PLINK](https://www.cog-genomics.org/plink/) | 1.9 | QC, PCA, IBD, heterozygosity |
| [KING](https://www.kingrelatedness.com) | ≥ 2.2 | Kinship estimation |
| [ADMIXTURE](https://dalexander.github.io/admixture/) | 1.3 | Ancestry inference |
| [IQ-TREE](http://www.iqtree.org) | ≥ 2.0 | Maximum likelihood phylogeny |
| [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay) | ≥ 3.4 | LD decay analysis |

### Python Packages

```bash
pip install pandas numpy matplotlib seaborn scikit-allel geopandas scikit-bio ete3
```

### R Packages

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap",
                   "ape", "dendextend", "RColorBrewer"))
```

---

## Repository Structure

```
.
├── Map/
│   └── sampling_map.py           # Sampling location map (Jazan, Saudi Arabia)
├── Metrics/
│   ├── nucleotide_diversity.py   # Per-cluster nucleotide diversity (π)
│   └── ld_decay.py               # LD decay curve
├── PCA/
│   └── vis.py                    # PCA scatter plot with passport ellipses
├── Admixture/
│   ├── vis.r                     # ADMIXTURE stacked bar plot
│   └── make_tables.py            # Ho / He / FST summary tables
├── Kinship/
│   └── heat.r                    # KING kinship heatmap
├── IBD/
│   └── heatmap.py                # IBD PI_HAT heatmap
├── pIBS/
│   ├── vis.py                    # Neighbor-joining tree from IBS matrix
│   └── vs.r                      # Dendrogram from IBS distance matrix
├── FingerPrint/
│   ├── gen.py                    # Generate allele-call CSV for 150-SNP panel
│   └── heatmapfinger.py          # Genotype heatmap
├── Tree/
│   ├── vcpTofaasta.py            # VCF → FASTA alignment for IQ-TREE
│   └── vis.py                    # Phylogenetic dendrogram visualization
└── QC/
    ├── hetro.py                  # Heterozygosity bar plot
    └── bar.py                    # Missingness bar plot
```

---

## Pipeline

### Step 1 — Variant Calling & Filtering (GATK)

SNPs were called using GATK HaplotypeCaller in GVCF mode per sample, then jointly genotyped with `GenomicsDBImport` and `GenotypeGVCFs`.

Hard-filtering thresholds applied to the raw SNP callset:

```bash
gatk VariantFiltration \
  -R reference.fa \
  -V raw_snps.vcf.gz \
  --filter-expression "MQ < 37.0"    --filter-name "MQ37" \
  --filter-expression "QD < 24.0"    --filter-name "QD24" \
  --filter-expression "FS > 60.0"    --filter-name "FS60" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
  -O filtered_snps.vcf.gz
```

**Threshold rationale:**
- **MQ > 37**: Slightly relaxed from the GATK-recommended 40 to retain variants in repeat-rich regions of the mango genome while maintaining mapping error probability < 2 × 10⁻⁴.
- **QD ≥ 24**: Conservative threshold (vs. GATK default ≥ 2) to ensure high-confidence calls at 15–20× sequencing depth.

---

### Step 2 — Sample QC: Missingness

```bash
# Calculate per-sample missingness
plink --vcf filtered_snps.vcf.gz \
      --keep-allele-order \
      --missing \
      --out QC/mango_SNPsOnly_no_M1toM9_missing

# Exclude samples with > 20% missingness
plink --vcf filtered_snps.vcf.gz \
      --mind 0.2 \
      --make-bed \
      --out QC/mango_after_mind_filter
```

**Visualise:**
```bash
# Run from QC/ directory (with .imiss file present)
python QC/bar.py
# Output: mango_SNPsOnly_no_M1toM9_missingness_barplot.png
```

---

### Step 3 — Sample QC: Heterozygosity

```bash
plink --bfile QC/mango_after_mind_filter \
      --het \
      --out QC/mango_QC_v5_het
```

**Visualise:**
```bash
# Run from QC/ directory (with .het and .fam files present)
python QC/hetro.py
# Output: heterozygosity_barplot_accessions.png
# Outlier threshold: F > +0.4 or F < -0.4
```

---

### Step 4 — LD Pruning

After missingness and heterozygosity QC, LD pruning was applied to produce the final dataset used for all downstream analyses.

```bash
# Step 4a — identify SNPs to keep after LD pruning
plink --bfile QC/mango_after_mind_filter \
      --indep-pairwise 50 10 0.1 \
      --out Raw-After-QC/mango_QC_v5_LD

# Step 4b — extract pruned SNP set
plink --bfile QC/mango_after_mind_filter \
      --extract Raw-After-QC/mango_QC_v5_LD.prune.in \
      --make-bed \
      --out Raw-After-QC/mango_QC_v5_LD_final
```

**Parameters:** window = 50 SNPs, step = 10, r² threshold = 0.10
**Result:** ~840,000 SNPs retained across 60 samples and 20 chromosomes.

> All subsequent steps use: `Raw-After-QC/mango_QC_v5_LD_final.{bed,bim,fam}`

---

### Step 5 — Principal Component Analysis (PCA)

```bash
plink --bfile Raw-After-QC/mango_QC_v5_LD_final \
      --pca 20 \
      --out PCA/mango_QC_v5_pca
```

**Visualise** (requires `PCA/passport_cleaned_sorted.csv`):
```bash
python PCA/vis.py
# Output: PCA/Mango_PCA_with_Passport_Ellipses.png (300 DPI)
```

Passport metadata columns: `SentCode`, `LocalNames`, `countryoforigin`
Confidence ellipses drawn for countries with ≥ 3 samples.

---

### Step 6 — Ancestry Inference (ADMIXTURE)

```bash
# Run ADMIXTURE for K = 2 to 10 with 10-fold cross-validation
for K in 2 3 4 5 6 7 8 9 10; do
    admixture --cv=10 Raw-After-QC/mango_QC_v5_LD_final.bed $K \
              | tee Admixture/log_K${K}.out
done

# Collect CV errors
grep "CV error" Admixture/log_K*.out > Admixture/cv_errors.txt
```

**K = 3** was selected based on minimum cross-validation error.

**Visualise:**
```bash
# Run from Admixture/ directory
Rscript Admixture/vis.r
# Output: admixture_K3_nogap.png / .pdf (400 DPI)
```

---

### Step 7 — Kinship Estimation (KING)

```bash
king --bfile Raw-After-QC/mango_QC_v5_LD_final \
     --kinship \
     --prefix Kinship/mango_QC_v5_kinship
```

**Visualise:**
```bash
# Run from Kinship/ directory
Rscript Kinship/heat.r
# Output: kinship_relationship_clustermap.png (300 DPI)
```

Kinship categories used:

| Kinship coefficient | Relationship |
|---------------------|-------------|
| > 0.354 | Duplicate / MZ twin |
| 0.177 – 0.354 | 1st-degree |
| 0.088 – 0.177 | 2nd-degree |
| 0.044 – 0.088 | 3rd-degree |
| < 0.044 | Unrelated |

---

### Step 8 — Identity by Descent (IBD)

```bash
plink --bfile Raw-After-QC/mango_QC_v5_LD_final \
      --genome \
      --out IBD/mango_QC_v5_ibd
```

**Visualise:**
```bash
# Run from IBD/ directory
python IBD/heatmap.py
# Output: ibd_heatmap_accession_final_publication.png (600 DPI)
# Duplicates flagged at PI_HAT > 0.95
```

---

### Step 9 — Pairwise IBS Tree

```bash
# Compute pairwise IBS matrix (homozygous SNPs only)
plink --bfile Raw-After-QC/mango_QC_v5_LD_final \
      --distance ibs \
      --out pIBS/mango_QC_v5_pIBS
```

**Visualise:**
```bash
# Neighbor-joining tree (Python / ete3)
python pIBS/vis.py
# Output: unrooted_tree.png (300 DPI)

# Dendrogram (R)
Rscript pIBS/vs.r
# Output: IBS_dendrogram.png
```

---

### Step 10 — Phylogenetic Tree (IQ-TREE)

Convert LD-pruned VCF to FASTA alignment:

```bash
python Tree/vcpTofaasta.py
# Input:  Tree/mango_LD_final.vcf
# Output: Tree/mango_snps_alignment.fasta
```

Run IQ-TREE with GTR+G model and 1000 ultrafast bootstrap replicates:

```bash
iqtree -s Tree/mango_snps_alignment.fasta \
       -m GTR+G \
       -bb 1000 \
       -pre Tree/mango_tree \
       -nt AUTO
# Output: Tree/mango_tree.treefile (Newick format)
```

**Visualise:**
```bash
python Tree/vis.py
# Output: phylogenetic dendrogram PNG
```

---

### Step 11 — SNP Fingerprinting Panel

MAF and PIC were calculated for all SNPs. The top 150 high-information SNPs were selected by:
- MAF > 0.25
- PIC > 0.25
- Heterozygosity filter: P(AB) > 0.30
- Even distribution across all 20 chromosomes

```bash
# Extract genotypes for the top 150 SNPs
plink --bfile Raw-After-QC/mango_QC_v5_LD_final \
      --extract FingerPrint/top150_mango_SNP_list.txt \
      --recode A \
      --out FingerPrint/top150_mango_genotypes

# Generate allele-call matrix
python FingerPrint/gen.py
# Output: FingerPrint/top150_mango_genotypes_alleles.csv

# Generate heatmap
python FingerPrint/heatmapfinger.py
# Output: FingerPrint/top150_genotype_heatmap_final_legend_bold_fixed.png (300 DPI)
```

---

### Step 12 — Population Statistics (Ho, He, FST)

```bash
# Per-sample heterozygosity (also used for QC above)
plink --bfile Raw-After-QC/mango_QC_v5_LD_final \
      --het \
      --out Raw-After-QC/Hetro/mango_QC_v5_het

# Generate Ho, He, and pairwise FST tables (K=3 clusters)
# Run from Admixture/ directory
python Admixture/make_tables.py
```

**Outputs:** `table_Ho_summary_k3.tsv`, `table_He_summary_k3.tsv`, `table_FST_matrix_k3.tsv`

**Summary (K = 3 ADMIXTURE clusters):**

| Cluster | N | Mean Ho | Mean He | Weighted FST vs C1 | Weighted FST vs C2 |
|---------|---|---------|---------|--------------------|--------------------|
| Cluster 1 | 19 | 0.305 | 0.300 | — | 0.048 |
| Cluster 2 | 25 | 0.312 | 0.300 | 0.048 | — |
| Cluster 3 | 16 | 0.286 | 0.300 | 0.089 | 0.065 |

---

### Step 13 — Nucleotide Diversity (π)

Computed from the LD-pruned VCF using `scikit-allel`:

```bash
python Metrics/nucleotide_diversity.py
# Input:  Tree/mango_LD_final.vcf
#         Admixture/admixture_cluster_k3.txt
# Output: Metrics/nucleotide_diversity_summary.tsv
#         Metrics/nucleotide_diversity_perchrom.tsv
#         Metrics/nucleotide_diversity_perchrom.png  (300 DPI)
```

**Genome-wide π:**

| Group | N | Mean π |
|-------|---|--------|
| All samples | 60 | 5.87 × 10⁻⁴ |
| Cluster 1 | 19 | 5.58 × 10⁻⁴ |
| Cluster 2 | 25 | 5.78 × 10⁻⁴ |
| Cluster 3 | 16 | 5.41 × 10⁻⁴ |

---

### Step 14 — LD Decay

LD decay was computed from the full (non-pruned) SNP dataset using PopLDdecay:

```bash
PopLDdecay -InVCF mango_LD_final.vcf \
           -OutStat Metrics/pop.ld_decay \
           -MaxDist 1000 \
           -MAF 0.05

# Visualise
python Metrics/ld_decay.py
# Output: Metrics/ld_decay.png       (300 DPI)
#         Metrics/ld_decay_summary.tsv
```

**Result:** LD half-decay distance = **39.4 kb**; r² declines to 0.1 at **142.5 kb**.

---

### Step 15 — Sampling Map

```bash
python Map/sampling_map.py
# Output: Map/sampling_map.png  (300 DPI)
# Shows Saudi Arabia with Jazan sampling site marked
```

---

## Key Results

| Metric | Value |
|--------|-------|
| Total samples | 60 accessions |
| SNPs after QC & LD pruning | ~840,000 |
| Chromosomes | 20 |
| ADMIXTURE K (optimal) | K = 3 |
| Genome-wide π | 5.87 × 10⁻⁴ per bp |
| LD half-decay distance | 39.4 kb |
| Fingerprint panel size | 150 SNPs |
| Pairwise FST (max) | 0.089 (Cluster 1 vs Cluster 3) |

---

## Citation

> [Authors]. Whole-Genome Sequencing Reveals Genomic Diversity and Population
> Structure of Saudi Arabian Mango Germplasm. *[Journal]*, [Year].
> doi: [DOI]

---

## Data Availability

Raw sequencing reads are deposited in the NCBI Sequence Read Archive (SRA).
**Accession numbers:** [to be provided upon acceptance]

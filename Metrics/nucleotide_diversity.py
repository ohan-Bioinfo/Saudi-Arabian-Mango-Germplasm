#!/usr/bin/env python3
"""
reviewer/Metrics/nucleotide_diversity.py

Compute nucleotide diversity (pi) per ADMIXTURE cluster and per chromosome.

Inputs (relative to QC-5 root):
    Tree/mango_LD_final.vcf
    Admixture/admixture_cluster_k3.txt

Outputs (written to reviewer/Metrics/):
    nucleotide_diversity_summary.tsv        genome-wide summary per group
    nucleotide_diversity_perchrom.tsv       per-chromosome pi values
    nucleotide_diversity_perchrom.png       grouped bar chart (300 DPI)

Dependencies:
    scikit-allel, pandas, numpy, matplotlib
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import allel
except ImportError:
    sys.exit("ERROR: scikit-allel is required.  Install: pip install scikit-allel")

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BASE_DIR     = os.path.join(SCRIPT_DIR, "..", "..")
VCF_PATH     = os.path.join(BASE_DIR, "Tree",      "mango_LD_final.vcf")
CLUSTER_PATH = os.path.join(BASE_DIR, "Admixture", "admixture_cluster_k3.txt")
OUT_SUMMARY  = os.path.join(SCRIPT_DIR, "nucleotide_diversity_summary.tsv")
OUT_CHROM    = os.path.join(SCRIPT_DIR, "nucleotide_diversity_perchrom.tsv")
OUT_PNG      = os.path.join(SCRIPT_DIR, "nucleotide_diversity_perchrom.png")

# ── Load VCF ───────────────────────────────────────────────────────────────
print("Loading VCF …")
# Read sample names from header first
headers     = allel.read_vcf_headers(VCF_PATH)
samples_raw = list(headers.samples)

callset = allel.read_vcf(VCF_PATH,
                         fields=["CHROM", "POS", "calldata/GT"],
                         alt_number=1)

# VCF sample IDs are "M55_M55" format; strip suffix to get "M55"
samples_ids = [s.split("_")[0] for s in samples_raw]
chrom_arr   = callset["variants/CHROM"]
pos_arr     = callset["variants/POS"]
gt          = allel.GenotypeArray(callset["calldata/GT"])

chroms = sorted(set(chrom_arr), key=lambda x: int(x))
print(f"  Variants : {len(pos_arr):,}")
print(f"  Samples  : {len(samples_ids)}")
print(f"  Chromosomes: {len(chroms)}")

# ── Load cluster assignments ───────────────────────────────────────────────
clusters_df = pd.read_csv(CLUSTER_PATH, sep="\t")
clusters_df.columns = clusters_df.columns.str.strip()
cluster_map   = dict(zip(clusters_df["IID"], clusters_df["Cluster"]))
cluster_names = sorted(clusters_df["Cluster"].unique())

cluster_indices = {}
for cname in cluster_names:
    idxs = [i for i, sid in enumerate(samples_ids)
             if cluster_map.get(sid) == cname]
    cluster_indices[cname] = idxs
    print(f"  {cname}: {len(idxs)} samples")

# ── Per-chromosome pi ──────────────────────────────────────────────────────
records = []

for chrom in chroms:
    mask  = chrom_arr == chrom
    pos_c = pos_arr[mask]
    gt_c  = gt.compress(mask, axis=0)

    if len(pos_c) < 2:
        continue

    start, stop = int(pos_c[0]), int(pos_c[-1])

    # All samples
    ac_all = gt_c.count_alleles()
    pi_all = allel.sequence_diversity(pos_c, ac_all, start=start, stop=stop)

    row = {"Chromosome": chrom, "N_SNPs": int(mask.sum()), "Pi_All": pi_all}

    for cname, cidx in cluster_indices.items():
        if len(cidx) < 2:
            row[f"Pi_{cname}"] = np.nan
            continue
        gt_cl = gt_c.take(cidx, axis=1)
        ac_cl = gt_cl.count_alleles()
        row[f"Pi_{cname}"] = allel.sequence_diversity(
            pos_c, ac_cl, start=start, stop=stop)

    records.append(row)
    cluster_pi_str = "  ".join(
        f"{c}={row[f'Pi_{c}']:.5f}" for c in cluster_names)
    print(f"  Chr {chrom:>2s}  SNPs={mask.sum():>6,}  "
          f"Pi_All={pi_all:.5f}  {cluster_pi_str}")

df = pd.DataFrame(records)
df.to_csv(OUT_CHROM, sep="\t", index=False, float_format="%.8f")
print(f"\nSaved per-chromosome table: {OUT_CHROM}")

# ── Genome-wide weighted-mean summary ─────────────────────────────────────
pi_cols = ["Pi_All"] + [f"Pi_{c}" for c in cluster_names]
summary_rows = []
for col in pi_cols:
    group    = col.replace("Pi_", "")
    n_smp    = (len(samples_ids) if group == "All"
                else len(cluster_indices.get(group, [])))
    valid    = df[col].dropna()
    weights  = df.loc[valid.index, "N_SNPs"]
    pi_mean  = float(np.average(valid, weights=weights))
    pi_sd    = float(np.sqrt(np.average((valid - pi_mean) ** 2, weights=weights)))
    summary_rows.append({
        "Group":   group,
        "N":       n_smp,
        "Mean_Pi": round(pi_mean, 8),
        "SD_Pi":   round(pi_sd,   8),
        "Min_Pi":  round(float(valid.min()), 8),
        "Max_Pi":  round(float(valid.max()), 8),
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_SUMMARY, sep="\t", index=False)
print(f"Saved summary table: {OUT_SUMMARY}")
print("\n=== Nucleotide Diversity (pi) Summary ===")
print(summary_df.to_string(index=False))

# ── Per-chromosome bar chart ───────────────────────────────────────────────
bar_colors  = ["#444444", "#66c2a5", "#fc8d62", "#8da0cb"]
plot_cols   = pi_cols
plot_labels = ["All samples"] + cluster_names

n_groups  = len(plot_cols)
bar_width = 0.18
x         = np.arange(len(df))

fig, ax = plt.subplots(figsize=(18, 6))

for i, (col, label, color) in enumerate(zip(plot_cols, plot_labels, bar_colors)):
    offsets = x + (i - n_groups / 2.0 + 0.5) * bar_width
    ax.bar(offsets, df[col], width=bar_width,
           color=color, alpha=0.85, edgecolor="black",
           linewidth=0.4, label=label)

ax.set_xticks(x)
ax.set_xticklabels([f"Chr {c}" for c in df["Chromosome"]],
                   rotation=45, ha="right", fontsize=9, fontweight="bold")
ax.set_ylabel("Nucleotide diversity (π)", fontsize=12, fontweight="bold")
ax.set_xlabel("Chromosome", fontsize=12, fontweight="bold")
ax.set_title(
    "Per-chromosome nucleotide diversity (π) by ADMIXTURE cluster (K=3)\n"
    "60 Mangifera indica accessions | ~699,000 LD-pruned SNPs",
    fontsize=13, fontweight="bold"
)
ax.legend(fontsize=9, framealpha=0.9)
ax.yaxis.grid(True, linestyle="--", alpha=0.45)
ax.set_axisbelow(True)
ax.tick_params(axis="y", labelsize=9)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved figure: {OUT_PNG}")

#!/usr/bin/env python3
"""
Generate population genetics summary tables for Mango K=3 ADMIXTURE clusters.
Produces:
  - Table 1: Observed Heterozygosity (Ho) per cluster
  - Table 2: Expected Heterozygosity (He) per cluster
  - Table 3: Pairwise FST matrix
"""
import pandas as pd
import numpy as np

# ── Load per-sample Ho and He ─────────────────────────────────────────────────
ho_df = pd.read_csv("table_observed_heterozygosity_Ho_k3.tsv", sep="\t")
he_df = pd.read_csv("table_expected_heterozygosity_He_k3.tsv", sep="\t")

# ── Load pairwise FST summary ─────────────────────────────────────────────────
fst_df = pd.read_csv("table_pairwise_fst_k3.tsv", sep="\t")

# ─────────────────────────────────────────────────────────────────────────────
# Table 1 – Observed Heterozygosity (Ho) by cluster
# ─────────────────────────────────────────────────────────────────────────────
cluster_order = ["Cluster1", "Cluster2", "Cluster3"]

ho_summary = (
    ho_df.groupby("Cluster")["Ho"]
    .agg(N="count", Mean="mean", SD="std", Min="min", Max="max")
    .reindex(cluster_order)
    .round(6)
    .reset_index()
)
ho_summary.columns = ["Cluster", "N", "Mean_Ho", "SD_Ho", "Min_Ho", "Max_Ho"]
ho_summary.to_csv("table_Ho_summary_k3.tsv", sep="\t", index=False)
print("=== Table 1: Observed Heterozygosity (Ho) ===")
print(ho_summary.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Table 2 – Expected Heterozygosity (He) by cluster
# ─────────────────────────────────────────────────────────────────────────────
he_summary = (
    he_df.groupby("Cluster")["He"]
    .agg(N="count", Mean="mean", SD="std", Min="min", Max="max")
    .reindex(cluster_order)
    .round(6)
    .reset_index()
)
he_summary.columns = ["Cluster", "N", "Mean_He", "SD_He", "Min_He", "Max_He"]
he_summary.to_csv("table_He_summary_k3.tsv", sep="\t", index=False)
print("\n=== Table 2: Expected Heterozygosity (He) ===")
print(he_summary.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Table 3a – Pairwise FST (long format, both Mean and Weighted)
# ─────────────────────────────────────────────────────────────────────────────
fst_long = fst_df[["Cluster_A", "Cluster_B", "Mean_FST", "Weighted_FST"]].copy()
fst_long["Mean_FST"]     = fst_long["Mean_FST"].round(4)
fst_long["Weighted_FST"] = fst_long["Weighted_FST"].round(4)
fst_long.to_csv("table_FST_long_k3.tsv", sep="\t", index=False)
print("\n=== Table 3a: Pairwise FST (long format) ===")
print(fst_long.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Table 3b – Pairwise FST matrix (Weighted FST, lower-triangle)
# ─────────────────────────────────────────────────────────────────────────────
mat = pd.DataFrame(np.nan, index=cluster_order, columns=cluster_order)
for _, row in fst_df.iterrows():
    a, b = row["Cluster_A"], row["Cluster_B"]
    mat.loc[a, b] = round(row["Weighted_FST"], 4)
    mat.loc[b, a] = round(row["Weighted_FST"], 4)

# Diagonal = 0
for c in cluster_order:
    mat.loc[c, c] = 0.0

mat.to_csv("table_FST_matrix_k3.tsv", sep="\t")
print("\n=== Table 3b: Pairwise Weighted FST matrix ===")
print(mat.to_string())

# ─────────────────────────────────────────────────────────────────────────────
# Combined Table – Ho, He, N side by side
# ─────────────────────────────────────────────────────────────────────────────
combined = ho_summary[["Cluster", "N", "Mean_Ho", "SD_Ho"]].merge(
    he_summary[["Cluster", "Mean_He", "SD_He"]], on="Cluster"
)
combined.to_csv("table_HoHe_combined_k3.tsv", sep="\t", index=False)
print("\n=== Combined Table: Ho & He per cluster ===")
print(combined.to_string(index=False))

print("\nAll tables saved.")

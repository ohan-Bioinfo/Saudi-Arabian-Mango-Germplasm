import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# === Step 1: Load the .raw file ===
df = pd.read_csv("mango_pure_homo.raw", delim_whitespace=True)

# === Step 2: Drop non-genotype columns ===
df = df.drop(columns=["FID", "PAT", "MAT", "SEX", "PHENOTYPE"])
sample_ids = df["IID"]
df = df.drop(columns=["IID"])
df.index = sample_ids  # Set Accession IDs as index

# === Step 3: Convert genotypes to numeric ===
geno = df.apply(pd.to_numeric, errors="coerce")

# === Step 4: Compute pairwise IBS (0% difference, homozygous only) ===
n = geno.shape[0]
ibs_matrix = np.full((n, n), np.nan)

for i in range(n):
    for j in range(i, n):
        row1 = geno.iloc[i, :]
        row2 = geno.iloc[j, :]
        valid = (~row1.isna()) & (~row2.isna())
        total = valid.sum()
        if total == 0:
            sim = np.nan
        else:
            matches = (row1[valid] == row2[valid]).sum()
            sim = matches / total
        ibs_matrix[i, j] = sim
        ibs_matrix[j, i] = sim

# === Step 5: Convert to DataFrame with accession IDs ===
ibs_df = pd.DataFrame(ibs_matrix, index=sample_ids, columns=sample_ids)

# === Step 6: Plot heatmap without dendrograms ===
sns.set(style="white")
plt.figure(figsize=(14, 12))

sns.heatmap(
    ibs_df,
    cmap="YlOrRd",
    linewidths=0.5,
    square=True,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "IBS (Homozygous SNPs only)"}
)

plt.title("Pairwise IBS Heatmap (Only Homozygous SNPs)", fontsize=14, weight="bold")
plt.xlabel("Accession ID", fontsize=12, weight="bold")
plt.ylabel("Accession ID", fontsize=12, weight="bold")

# === Bold & increase font size of sample IDs ===
plt.xticks(rotation=90, fontsize=9, fontweight="bold")
plt.yticks(fontsize=9, fontweight="bold")

plt.tight_layout()
plt.savefig("IBS_heatmap_homozygous_nodendro_boldIDs.png", dpi=600)
plt.show()

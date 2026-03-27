import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# === Load IBD genome file ===
df = pd.read_csv("mango_QC_v5_ibd.genome", delim_whitespace=True)

# === Prepare unique accession IDs ===
accessions = sorted(set(df["IID1"]).union(set(df["IID2"])))
n = len(accessions)
matrix = pd.DataFrame(0, index=accessions, columns=accessions, dtype=float)

# === Fill PI_HAT values ===
for _, row in df.iterrows():
    i, j = row["IID1"], row["IID2"]
    matrix.loc[i, j] = row["PI_HAT"]
    matrix.loc[j, i] = row["PI_HAT"]

# === Set diagonal to distinguishable value ===
for acc in accessions:
    matrix.loc[acc, acc] = 1.01

# === Custom perceptual colormap ===
colors = ["#F7FCFD", "#66C2A4", "#2CA25F", "#006D2C", "#000000"]  # from light to dark
cmap = LinearSegmentedColormap.from_list("pi_hat_gradient", colors, N=256)

# === Plot ===
plt.figure(figsize=(24, 22))
sns.heatmap(
    matrix,
    cmap=cmap,
    square=True,
    linewidths=0.6,
    linecolor='gray',
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "PI_HAT (IBD Estimate)", "shrink": 0.8}
)

# === Titles and Labels ===
plt.title(f"IBD Heatmap — PI_HAT Between Accessions (n = {n})", fontsize=24, fontweight='bold', pad=20)
plt.xlabel("Accession ID", fontsize=18, fontweight='bold', labelpad=12)
plt.ylabel("Accession ID", fontsize=18, fontweight='bold', labelpad=12)
plt.xticks(rotation=90, fontsize=14, fontweight='bold')  # Increased x-axis font
plt.yticks(fontsize=14, fontweight='bold')               # Increased y-axis font

# === Layout and Save ===
plt.tight_layout()
plt.savefig("ibd_heatmap_accession_final_publication.png", dpi=600)
plt.show()

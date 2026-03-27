import pandas as pd
import matplotlib.pyplot as plt

# === Load PLINK .het file and .fam file ===
het = pd.read_csv("mango_QC_v5_het.het", delim_whitespace=True)
fam = pd.read_csv("mango_QC_v5_LD_final.fam", delim_whitespace=True, header=None)
fam.columns = ['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype']

# === Merge het and fam to get accession names ===
het = het.merge(fam[['FID', 'IID']], on=['FID', 'IID'])

# === Mark outliers ===
het['Outlier'] = het['F'].apply(lambda x: 'Outlier' if x > 0.4 or x < -0.4 else 'Normal')

# === Sort by F value and reset index ===
het_sorted = het.sort_values(by='F').reset_index(drop=True)

# === Set color per sample ===
colors = ['red' if val == 'Outlier' else 'steelblue' for val in het_sorted['Outlier']]

# === Plot ===
plt.figure(figsize=(16, 6))
bars = plt.bar(het_sorted['IID'], het_sorted['F'], color=colors, edgecolor='black')

# === Threshold lines ===
plt.axhline(0.4, linestyle='--', color='gray', linewidth=1)
plt.axhline(-0.4, linestyle='--', color='gray', linewidth=1)

# === Labels & style ===
plt.title("QC-5: Heterozygosity (F Coefficient) Per Sample", fontsize=14, fontweight='bold')
plt.xlabel("Accession ID", fontsize=12, fontweight='bold')
plt.ylabel("F (Inbreeding Coefficient)", fontsize=12, fontweight='bold')
plt.xticks(rotation=90, fontsize=9, fontweight='bold')
plt.yticks(fontsize=10)
plt.tight_layout()

# === Save ===
plt.savefig("heterozygosity_barplot_accessions.png", dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# === Load .imiss and .fam ===
imiss = pd.read_csv("mango_SNPsOnly_no_M1toM9_missing.imiss", delim_whitespace=True)
fam = pd.read_csv("mango_SNPsOnly_no_M1toM9.fam", delim_whitespace=True, header=None)
fam.columns = ['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype']

# === Merge accession info ===
data = imiss.merge(fam[['FID', 'IID']], on=['FID', 'IID'])

# === Define outlier threshold ===
threshold = 0.2
data['Outlier'] = data['F_MISS'].apply(lambda x: 'Outlier' if x > threshold else 'Normal')

# === Sort by missing rate and assign color ===
data_sorted = data.sort_values(by='F_MISS', ascending=False).reset_index(drop=True)
colors = ['red' if o == 'Outlier' else 'lightgray' for o in data_sorted['Outlier']]

# === Plot ===
plt.figure(figsize=(18, 6))
plt.bar(data_sorted['IID'], data_sorted['F_MISS'], color=colors, edgecolor='black')

# === Labels & Styling ===
plt.title("Sample Missingness Rate (Per Individual)", fontsize=14, fontweight='bold')
plt.xlabel("Accession ID", fontsize=12, fontweight='bold')
plt.ylabel("Missing Rate", fontsize=12, fontweight='bold')
plt.xticks(rotation=90, fontsize=9, fontweight='bold')
plt.yticks(fontsize=10)
plt.axhline(threshold, linestyle='--', color='gray', linewidth=1, label='Outlier Threshold (0.2)')
plt.legend()
plt.tight_layout()



plt.savefig("mango_SNPsOnly_no_M1toM9_missingness_barplot.png", dpi=300)
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# === Load genotype label matrix ===
geno_labels = pd.read_csv("geno_labels_top150.csv", index_col=0)

# === Define genotype color map ===
color_dict = {
    "A/A": "#87CEFA",    # Light Blue
    "T/T": "#FFD700",    # Gold
    "G/G": "#32CD32",    # Green
    "C/C": "#FF4500",    # Red-Orange
    "Hetero": "#A9A9A9", # Gray
    "NA": "#FFFFFF"      # White
}

# === Map genotypes to integers for plotting ===
value_to_int = {label: i for i, label in enumerate(color_dict)}
geno_numeric = geno_labels.applymap(lambda x: value_to_int.get(x, len(color_dict)))
geno_array = geno_numeric.values

# === Plot setup ===
fig, ax = plt.subplots(figsize=(28, 20))
cmap = ListedColormap([color_dict[label] for label in color_dict])
im = ax.imshow(geno_array, aspect='auto', cmap=cmap)

# === Format X-axis: SNP labels on top ===
ax.set_xticks(np.arange(geno_array.shape[1]))
ax.set_xticklabels(geno_labels.columns, rotation=90, fontsize=9, fontweight='bold')
ax.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)

# === Format Y-axis: Accession IDs on left ===
ax.set_yticks(np.arange(geno_array.shape[0]))
ax.set_yticklabels(geno_labels.index, fontsize=15, fontweight='bold')

# === Titles ===
ax.set_xlabel("SNPs", fontsize=18, fontweight='bold')
ax.set_ylabel("Accession ID", fontsize=18, fontweight='bold')
plt.title("Genotype Heatmap (Top 150 SNPs)", fontsize=20, fontweight='bold')

# === Build Manual Bold Legend ===
legend_elements = [
    Patch(facecolor=color_dict[label], edgecolor='black', label=label)
    for label in color_dict
]

# Draw the legend manually in an external axis
from matplotlib.legend import Legend

legend_ax = fig.add_axes([0.91, 0.3, 0.05, 0.4])  # [left, bottom, width, height]
legend_ax.axis('off')  # Hide the box

custom_legend = Legend(
    legend_ax,
    handles=legend_elements,
    labels=[f"{label}" for label in color_dict],
    title="Genotype",
    title_fontsize=15,
    fontsize=14,
    loc='center'
)

# Apply bold to each label manually
for text in custom_legend.get_texts():
    text.set_fontweight('bold')
custom_legend.get_title().set_fontweight('bold')

legend_ax.add_artist(custom_legend)

# === Save ===
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Leave space for external legend
plt.savefig("top150_genotype_heatmap_final_legend_bold_fixed.png", dpi=300)
plt.close()

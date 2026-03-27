import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# === Load data ===
ids = pd.read_csv("pIBS/mango_QC_v5_pIBS.mibs.id", delim_whitespace=True, header=None)
ibs = pd.read_csv("pIBS/mango_QC_v5_pIBS.mibs", delim_whitespace=True, header=None)
labels = ids[1].tolist()

# === Convert IBS to distance ===
distance_matrix = 1 - ibs.values
condensed_dist = squareform(distance_matrix, checks=False)
Z = linkage(condensed_dist, method="average")

# === Plot vertical dendrogram (top-down) ===
plt.figure(figsize=(10, 18))  # taller for vertical tree
dendrogram(
    Z,
    labels=labels,
    orientation='left',  # vertical orientation
    leaf_font_size=12,
    color_threshold=0.5
)
plt.title("Phylogenetic Tree", fontsize=16, weight='bold')
plt.xlabel("Genetic Distance", fontsize=13)
plt.tight_layout()
plt.savefig("Tree/mango_phylo_tree_vertical.png", dpi=600)
plt.show()

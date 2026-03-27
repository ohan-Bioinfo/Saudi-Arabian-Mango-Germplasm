import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree, TreeStyle, TextFace

# === Step 1: Load IBS matrix and accession IDs ===
ibs = np.loadtxt("mango_QC_v5_pIBS.mibs")
ids = pd.read_csv("mango_QC_v5_pIBS.mibs.id", sep="\s+", header=None)
accessions = ids[1].tolist()

# === Step 2: Convert IBS to distance matrix ===
distance = 1 - ibs

# === Step 3: Create DistanceMatrix and build Neighbor-Joining tree ===
dm = DistanceMatrix(distance, accessions)
tree = nj(dm)

# === Step 4: Save Newick tree to file ===
with open("unrooted_tree.nwk", "w") as f:
    f.write(str(tree))

# === Step 5: Load tree using ete3 and render unrooted (circular) layout ===
t = Tree("unrooted_tree.nwk")

# Configure tree style
ts = TreeStyle()
ts.mode = "c"  # circular mode (unrooted look)
ts.show_leaf_name = True
ts.title.add_face(TextFace("Unrooted Genetic Distance Tree", fsize=14, bold=True), column=0)

# Save to high-res PNG
t.render("unrooted_tree.png", w=1000, dpi=300, tree_style=ts)

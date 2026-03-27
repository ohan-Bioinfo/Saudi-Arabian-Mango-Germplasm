#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import math
from collections import Counter

# --------------------------
# Config (tweak as you like)
# --------------------------
PC_AXES = [(1, 2)]            # add (2,3) if you also want PC2 vs PC3
DRAW_ELLIPSES = True          # draw country confidence ellipses when n>=3
ANNOTATE = False              # set True to annotate Named cultivars with IID
DPI = 300
OUT_BASENAME = "PCA/pca_by_passport"

# --------------------------
# Load data
# --------------------------
pca = pd.read_csv("PCA/mango_QC_v5_pca.eigenvec", sep=r"\s+", header=None)
# First two columns are FID, IID; rest are PCs
pca_cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca.shape[1]-1)]
pca.columns = pca_cols
pca["IID"] = pca["IID"].astype(str)

# Variance explained for axis labels
try:
    eigenval = np.loadtxt("PCA/mango_QC_v5_pca.eigenval")
    var_exp = (eigenval / eigenval.sum()) * 100.0
except Exception:
    var_exp = None

passport = pd.read_csv("PCA/passport_cleaned_sorted.csv")
passport.columns = passport.columns.str.strip()
passport = passport.rename(columns={"Sent Code": "IID"})
passport["IID"] = passport["IID"].astype(str)

# Clean whitespace & simple normalizations
passport["Local Names"] = passport["Local Names"].astype(str).str.strip()
passport["country of origin"] = passport["country of origin"].astype(str).str.strip()
passport.loc[passport["country of origin"].isin(["", "nan", "NaN"]), "country of origin"] = "Unknown"

# Merge
df = pca.merge(passport, on="IID", how="left")
df["country of origin"] = df["country of origin"].fillna("Unknown")
# Define Type (marker)
df["Type"] = np.where(df["Local Names"].str.lower().str.replace(r"\s+", "", regex=True) == "new_variety",
                      "New variety", "Named cultivar")

# --------------------------
# Palette & markers
# --------------------------
countries = df["country of origin"].tolist()
counts = Counter(countries)

# Put KSA first in legend, then by frequency descending
unique_countries = list({ "KSA": None, **{c:None for c,_ in sorted(counts.items(), key=lambda x: -x[1])} }.keys())
if "KSA" not in counts:
    unique_countries = [c for c,_ in sorted(counts.items(), key=lambda x: -x[1])]

# Build a color map using tab20 (cycled if needed)
tab20 = plt.get_cmap("tab20").colors
colors = []
while len(colors) < len(unique_countries):
    colors.extend(tab20)
country2color = {c: colors[i] for i, c in enumerate(unique_countries)}

type2marker = {"Named cultivar": "s", "New variety": "o"}  # square vs circle

# --------------------------
# Ellipse helper
# --------------------------
def draw_confidence_ellipse(ax, x, y, facecolor, edgecolor, alpha=0.12, lw=1.2, conf=0.95):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if np.linalg.matrix_rank(cov) < 2:
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    # Chi-square quantile for 2 dof at ~95% ≈ 5.991; radius scale = sqrt(q)
    q = 5.991 if conf == 0.95 else 2.4477**2
    width, height = 2 * np.sqrt(vals * q)
    ellipse = Ellipse((np.mean(x), np.mean(y)), width, height, theta,
                      facecolor=facecolor, edgecolor=edgecolor, lw=lw, alpha=alpha)
    ax.add_patch(ellipse)

# --------------------------
# Plotting
# --------------------------
for (pcx, pcy) in PC_AXES:
    xcol, ycol = f"PC{pcx}", f"PC{pcy}"
    fig, ax = plt.subplots(figsize=(10, 8), dpi=DPI)

    # Draw ellipses first (under points), only for groups with >=3
    if DRAW_ELLIPSES:
        for c in unique_countries:
            sub = df[df["country of origin"] == c]
            if len(sub) >= 3:
                draw_confidence_ellipse(ax, sub[xcol], sub[ycol],
                                        facecolor=country2color[c],
                                        edgecolor=country2color[c],
                                        alpha=0.10, lw=1.0, conf=0.95)

    # Scatter points
    # Color = country; Marker = type
    for t, marker in type2marker.items():
        sub = df[df["Type"] == t]
        ax.scatter(sub[xcol], sub[ycol],
                   c=[country2color[c] for c in sub["country of origin"]],
                   marker=marker, s=80, linewidths=0.6, edgecolors="black", alpha=0.95, label=t)

    # Optional annotations (keep tidy)
    if ANNOTATE:
        ann = df[df["Type"] == "Named cultivar"]  # annotate only named cultivars
        for _, r in ann.iterrows():
            ax.text(r[xcol], r[ycol], r["IID"], fontsize=8, ha="left", va="bottom")

    # Legends: one for country (colors), one for Type (markers)
    country_handles = [Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=country2color[c], markeredgecolor="black",
                              label=f"{c} (n={counts[c]})", markersize=8)
                       for c in unique_countries]
    type_handles = [Line2D([0], [0], marker=m, color='black', linestyle='',
                           label=t, markersize=8) for t, m in type2marker.items()]

    leg1 = ax.legend(handles=country_handles, title="Country of origin", loc="upper right",
                     bbox_to_anchor=(1.02, 1), frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=type_handles, title="Sample type", loc="lower right",
              bbox_to_anchor=(1.02, 0), frameon=True)

    # Axis labels (with % variance if available)
    if var_exp is not None and len(var_exp) >= max(pcx, pcy):
        ax.set_xlabel(f"PC{pcx} ({var_exp[pcx-1]:.2f}% var.)")
        ax.set_ylabel(f"PC{pcy} ({var_exp[pcy-1]:.2f}% var.)")
    else:
        ax.set_xlabel(f"PC{pcx}")
        ax.set_ylabel(f"PC{pcy}")

    ax.set_title("Mango PCA by Passport: color=Country, marker=Type", pad=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linewidth=0.4, alpha=0.4)

    plt.tight_layout()
    out = f"{OUT_BASENAME}_PC{pcx}_PC{pcy}.png"
    plt.savefig(out, dpi=DPI)
    print(f"Saved: {out}")

# Also export the color key for reference
pd.DataFrame({
    "country": unique_countries,
    "color_rgb": [country2color[c] for c in unique_countries]
}).to_csv("PCA/country_color_map.csv", index=False)
print("Saved: PCA/country_color_map.csv")

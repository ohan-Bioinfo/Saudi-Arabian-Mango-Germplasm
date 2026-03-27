#!/usr/bin/env python3
"""
reviewer/Map/sampling_map.py

Map of Saudi Arabia with Jazan region highlighted as the sampling site
for all 60 mango accessions.

Output:
    reviewer/Map/sampling_map.png   (300 DPI)

Dependencies:
    geopandas, matplotlib, numpy
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import geopandas as gpd
from shapely.geometry import Point

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PNG    = os.path.join(SCRIPT_DIR, "sampling_map.png")

# ── Jazan coordinates (city centre) ───────────────────────────────────────
JAZAN_LON, JAZAN_LAT = 42.558, 16.889

# ── Load world shapefile ───────────────────────────────────────────────────
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

# Separate Saudi Arabia and neighbours
ksa        = world[world["name"] == "Saudi Arabia"]
neighbours = world[world["name"].isin([
    "Yemen", "Oman", "United Arab Emirates", "Qatar",
    "Bahrain", "Kuwait", "Iraq", "Jordan", "Egypt", "Sudan", "Eritrea",
])]
other      = world[~world.index.isin(ksa.index) & ~world.index.isin(neighbours.index)]

# ── Plot ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 9))
fig.patch.set_facecolor("white")

# Ocean background
ax.set_facecolor("#d6e8f5")

# Other countries (muted)
other.plot(ax=ax, color="#e0e0e0", edgecolor="#bbbbbb", linewidth=0.3)

# Neighbouring countries
neighbours.plot(ax=ax, color="#c8d8b0", edgecolor="#999999", linewidth=0.4)

# Saudi Arabia — highlighted
ksa.plot(ax=ax, color="#f5c87a", edgecolor="#666666", linewidth=0.8)

# Jazan marker
ax.scatter(JAZAN_LON, JAZAN_LAT,
           s=300, marker="*", color="#d62728",
           edgecolors="black", linewidth=0.8, zorder=8)

# Jazan label with arrow
ax.annotate(
    "Jazan\n(sampling site\nn = 60)",
    xy=(JAZAN_LON, JAZAN_LAT),
    xytext=(JAZAN_LON - 8.5, JAZAN_LAT + 4.5),
    arrowprops=dict(arrowstyle="->", color="#222222",
                    lw=1.3, connectionstyle="arc3,rad=-0.2"),
    fontsize=12, fontweight="bold", color="#222222",
    bbox=dict(boxstyle="round,pad=0.4", fc="white",
              ec="#d62728", alpha=0.95, linewidth=1.2)
)

# "Saudi Arabia" label inside the country
ax.text(44.5, 24.5, "Saudi Arabia",
        fontsize=13, fontweight="bold", color="#555533",
        ha="center", va="center", style="italic")

# Map extent — Arabian Peninsula + context
ax.set_xlim(25, 65)
ax.set_ylim(8, 38)

ax.set_xlabel("Longitude (°E)", fontsize=11, fontweight="bold")
ax.set_ylabel("Latitude (°N)", fontsize=11, fontweight="bold")
ax.set_title(
    "Sampling location of mango accessions\nJazan region, Saudi Arabia",
    fontsize=13, fontweight="bold", pad=10
)
ax.tick_params(labelsize=9)

# Legend
legend_handles = [
    mpatches.Patch(facecolor="#f5c87a", edgecolor="#666666", label="Saudi Arabia"),
    mpatches.Patch(facecolor="#c8d8b0", edgecolor="#999999", label="Neighbouring countries"),
    plt.scatter([], [], s=200, marker="*", color="#d62728",
                edgecolors="black", linewidth=0.7, label="Jazan (sampling site)"),
]
ax.legend(handles=legend_handles, fontsize=10, loc="upper right",
          framealpha=0.95, edgecolor="#aaaaaa")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_PNG}")

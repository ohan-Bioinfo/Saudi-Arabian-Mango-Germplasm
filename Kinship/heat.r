# === Install and load required package ===
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# === Load kinship file ===
kin <- read.table("mango_QC_v5_kinship.kin0", header = FALSE)
colnames(kin) <- c("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP")

# === Define categorization function ===
categorize_kin <- function(k) {
  if (k > 0.354) return("Duplicates")
  else if (k > 0.177) return("1st-degree")
  else if (k > 0.0884) return("2nd-degree")
  else if (k > 0.0442) return("3rd-degree")
  else return("Unrelated")
}

# === Apply categories ===
kin$Category <- sapply(kin$KINSHIP, categorize_kin)

# === Get unique sample IDs ===
samples <- sort(unique(c(kin$IID1, kin$IID2)))

# === Create category matrix ===
cat_matrix <- matrix("NA", nrow = length(samples), ncol = length(samples),
                     dimnames = list(samples, samples))

for (i in 1:nrow(kin)) {
  id1 <- kin$IID1[i]
  id2 <- kin$IID2[i]
  cat <- kin$Category[i]
  cat_matrix[id1, id2] <- cat
  cat_matrix[id2, id1] <- cat
}
diag(cat_matrix) <- "Duplicates"  # self-self

# === Convert to integer matrix for pheatmap ===
relationship_levels <- c("Duplicates", "1st-degree", "2nd-degree", "3rd-degree", "Unrelated")
rel_to_int <- setNames(seq_along(relationship_levels), relationship_levels)

cat_int_matrix <- matrix(rel_to_int[cat_matrix], nrow = nrow(cat_matrix), dimnames = dimnames(cat_matrix))

# === Color palette ===
color_palette <- c(
  "Duplicates" = "#FF0000",
  "1st-degree" = "#FF9900",
  "2nd-degree" = "#FFD700",
  "3rd-degree" = "#87CEFA",
  "Unrelated"  = "#D3D3D3"
)

# === Output to PNG ===
png("kinship_relationship_clustermap.png", width = 2400, height = 2200, res = 300)

pheatmap(
  cat_int_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = unname(color_palette[relationship_levels]),
  legend_breaks = rel_to_int,
  legend_labels = names(rel_to_int),
  border_color = "grey90",
  main = "Kinship Relationship Clustermap",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  fontface_row = "bold",
  fontface_col = "bold"
)

dev.off()

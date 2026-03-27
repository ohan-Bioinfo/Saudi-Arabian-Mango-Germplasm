# === Load libraries ===
library(ggplot2)
library(dplyr)
library(tidyr)

# === Load ADMIXTURE Q matrix and FAM file
q <- read.table("mango_QC_v5_LD_final.3.Q")
fam <- read.table("mango_QC_v5_LD_final.fam")

# === Combine Q with sample IDs
colnames(q) <- paste0("Cluster", 1:3)
q$IID <- fam$V2

# === Sort by dominant ancestry
q$max_cluster <- apply(q[, 1:3], 1, which.max)
q$sort_key <- apply(q[, 1:3], 1, max)
q_sorted <- q[order(q$max_cluster, -q$sort_key), ]

# === Long format for plotting
q_long <- q_sorted %>%
  select(-sort_key) %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion")

q_long$IID <- factor(q_long$IID, levels = unique(q_sorted$IID))

# === Build plot with no gaps between bars
p <- ggplot(q_long, aes(x = IID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +  # full width = no gap
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "ADMIXTURE Plot (K = 3)",
    x = "Accession ID",
    y = "Ancestry Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(10, 20, 10, 20)
  )

# === Save image (PNG and PDF)
ggsave("admixture_K3_nogap.png", plot = p, width = 10, height = 5, dpi = 400)
ggsave("admixture_K3_nogap.pdf", plot = p, width = 10, height = 5)

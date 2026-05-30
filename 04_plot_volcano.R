###############################
# 04_plot_volcano.R
# Volcano plot
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

result <- read.csv(file.path(table_dir, paste0(comparison_name, "_all_genes.csv")))

result$state <- "Not statistically significant"

result$state[
  result$padj < padj_cutoff & result$log2FoldChange > lfc_cutoff
] <- paste0("Higher in ", target_group)

result$state[
  result$padj < padj_cutoff & result$log2FoldChange < -lfc_cutoff
] <- paste0("Higher in ", reference_group)

label_genes <- result$gene[
  result$padj < 0.01 &
    abs(result$log2FoldChange) > lfc_cutoff
]

result$point_label <- ifelse(result$gene %in% label_genes, result$gene, "")

p <- ggplot(result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = state), size = 1.8, alpha = 0.8) +
  geom_text_repel(aes(label = point_label), max.overlaps = 25, size = 3) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  labs(
    title = paste0(target_group, " minus ", reference_group),
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Expression pattern"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(plot_dir, paste0("02_", comparison_name, "_volcano.pdf")),
  plot = p,
  width = 8,
  height = 6,
  dpi = plot_dpi
)

cat("Volcano plot saved\\n")

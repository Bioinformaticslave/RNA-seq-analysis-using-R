###############################
# 05_run_gsea.R
# GSEA using fgsea and MSigDB mouse C7 gene sets
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(fgsea)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

result <- read.csv(file.path(table_dir, paste0(comparison_name, "_all_genes.csv")))

load(msigdb_mouse_c7)  # expected object: Mm.c7
pathways <- Mm.c7

result$entrez <- mapIds(
  org.Mm.eg.db,
  keys = result$gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gsea_data <- result[!is.na(result$entrez), ]

if ("stat" %in% colnames(gsea_data)) {
  ranks <- setNames(gsea_data$stat, gsea_data$entrez)
} else {
  ranks <- setNames(gsea_data$log2FoldChange, gsea_data$entrez)
}

ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

fgsea_res <- fgseaMultilevel(pathways, stats = ranks)
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

write.csv(as.data.frame(fgsea_res),
          file.path(table_dir, paste0(comparison_name, "_fgsea_results.csv")),
          row.names = FALSE)

plot_data <- as.data.frame(fgsea_res) %>%
  arrange(padj) %>%
  slice_head(n = 30)

p <- ggplot(plot_data, aes(x = reorder(pathway, padj), y = -log10(padj))) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Adjusted p-values of GSEA pathways",
    subtitle = "Dashed line indicates adjusted p = 0.05",
    x = "Pathway",
    y = expression(-log[10](adjusted~p))
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(plot_dir, paste0("03_", comparison_name, "_gsea_padj_plot.pdf")),
  plot = p,
  width = 10,
  height = 8,
  dpi = plot_dpi
)

cat("GSEA complete\\n")

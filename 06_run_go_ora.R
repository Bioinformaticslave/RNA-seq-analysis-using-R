###############################
# 06_run_go_ora.R
# GO Biological Process ORA
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

result_all <- read.csv(file.path(table_dir, paste0(comparison_name, "_all_genes.csv")))

gene_vector <- result_all$gene[
  result_all$padj < padj_cutoff &
    abs(result_all$log2FoldChange) > lfc_cutoff
]

all_genes <- result_all$gene
all_genes <- all_genes[!is.na(all_genes)]

ego <- enrichGO(
  gene          = gene_vector,
  universe      = all_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

write.csv(as.data.frame(ego),
          file.path(table_dir, paste0(comparison_name, "_GO_ORA_BP.csv")),
          row.names = FALSE)

p <- dotplot(ego, showCategory = 20) +
  ggtitle(paste0("GO Biological Process ORA: ", comparison_name)) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(plot_dir, paste0("04_", comparison_name, "_GO_ORA_dotplot.pdf")),
  plot = p,
  width = 8,
  height = 10,
  dpi = plot_dpi
)

cat("GO ORA complete\\n")

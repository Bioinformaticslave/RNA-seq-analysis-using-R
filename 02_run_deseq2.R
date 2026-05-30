###############################
# 02_run_deseq2.R
# Run DESeq2 and extract target vs reference results
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
})

dds <- readRDS(file.path(object_dir, "dds_unfitted.rds"))
dds <- DESeq(dds)

res <- results(dds, contrast = c("group", target_group, reference_group))

result_all <- as.data.frame(res)
result_all$gene <- rownames(result_all)

result_all <- result_all[
  !is.na(result_all$log2FoldChange) &
    !is.na(result_all$padj),
]

result_sig <- result_all[
  result_all$padj < padj_cutoff &
    abs(result_all$log2FoldChange) > lfc_cutoff,
]

result_all <- result_all[order(result_all$padj), ]
result_sig <- result_sig[order(result_sig$padj), ]

saveRDS(dds, file.path(object_dir, "dds_fitted.rds"))
saveRDS(res, file.path(object_dir, paste0(comparison_name, "_DESeqResults.rds")))

write.csv(result_all, file.path(table_dir, paste0(comparison_name, "_all_genes.csv")), row.names = FALSE)
write.csv(result_sig, file.path(table_dir, paste0(comparison_name, "_significant_DEGs.csv")), row.names = FALSE)

cat("DESeq2 complete\\n")
cat("All analyzable genes:", nrow(result_all), "\\n")
cat("Significant DEGs:", nrow(result_sig), "\\n")

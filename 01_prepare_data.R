###############################
# 01_prepare_data.R
# Load counts/metadata and create DESeq2 object
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
})

load(countdata_rdata)   # expected object: countdata

metadata_DP <- read.csv(metadata_DP_csv, row.names = 1, stringsAsFactors = FALSE)
metadata_SP <- read.csv(metadata_SP_csv, row.names = 1, stringsAsFactors = FALSE)

metadata_DP$source <- "DP"
metadata_SP$source <- "SP"
metadata_all <- rbind(metadata_DP, metadata_SP)

counts_all <- countdata[, rownames(metadata_all)]

keep_genes <- rowSums(counts_all) > 10
counts_all <- counts_all[keep_genes, ]

metadata_all$group <- factor(metadata_all$group)

dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData   = metadata_all,
  design    = ~ group
)

saveRDS(dds, file.path(object_dir, "dds_unfitted.rds"))
saveRDS(metadata_all, file.path(object_dir, "metadata_all.rds"))
saveRDS(counts_all, file.path(object_dir, "counts_all_filtered.rds"))

cat("Prepared DESeq2 object\\n")
print(dds)

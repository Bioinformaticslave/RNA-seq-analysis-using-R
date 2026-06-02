###############################
# 03_plot_pca.R
# PCA plot for all samples
###############################

source("R/00_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS(file.path(object_dir, "dds_fitted.rds"))

# perform vsd transformation for more convenient later analysis
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(object_dir, "vsd.rds"))

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of immune-cell transcriptomes",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Group"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
# modify the Principal Component level accordingly for interest

# save the file accordingly
ggsave(
  filename = file.path(plot_dir, "01_PCA_all_samples.pdf"),
  plot = p,
  width = 8,
  height = 6,
  dpi = plot_dpi
)

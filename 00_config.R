###############################
# 00_config.R
# Project configuration
###############################

project_dir <- normalizePath(".", mustWork = FALSE)

data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
plot_dir    <- file.path(results_dir, "plots")
table_dir   <- file.path(results_dir, "tables")
object_dir  <- file.path(results_dir, "objects")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)

# Edit these paths for your Linux system.
countdata_rdata <- "/filepath/countdata" # Rows = genes, Columns = samples, Values = raw read counts
metadata_DP_csv <- "/filepath/metadata_1.csv"
metadata_SP_csv <- "/filepath/metadata_2.csv"
key_gene_file   <- "/filepath/key_genes_to_display.txt"
msigdb_mouse_c7 <- "/filepath/mouse_c7_v5p2.rdata"

# DESeq2 contrast is target minus reference.
# Positive log2FC = higher in target.
# Negative log2FC = higher in reference.
target_group    <- "A" # A can be any grouof interest
reference_group <- "B" # Same as B
comparison_name <- paste0(target_group, "_minus_", reference_group)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

candidate_genes <- c("Gene1", "Gene2")

plot_dpi <- 300

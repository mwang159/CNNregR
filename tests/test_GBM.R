library(CNNregR)

options(stringsAsFactors = FALSE)

data("pseudoBulk_GBM", package = "CNNregR")
data("ref_GBM", package = "CNNregR")
data("gene_annotation", package = "CNNregR")

work_dir <- file.path(tempdir(), "CNNregR_GBM_test")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
bulk_file <- file.path(work_dir, "bulk.csv")
ref_file <- file.path(work_dir, "ref.csv")
output_dir <- file.path(work_dir, "CNNreg_output")


###############################################################################
### pseudoBulk_GBM$pseudoBulk is a Sample X Gene count matrix
bulk_counts <- pseudoBulk_GBM$pseudoBulk

###############################################################################
### Scenario 1: use one single scRNA-seq reference
## ref_GBM is a list of count matrices, one for each cell type.
sc_counts_list <- list(ref_GBM)

cat("\n=== Test: Check CNNreg Installation ===\n")
check_cnnreg_install()

cat("\n=== Test: Validate Input Data ===\n")
validation <- validate_input_data(bulk_counts, sc_counts_list)

cat("\n=== Test: Select Marker Genes ===\n")
gene_result <- select_genes(
  bulk_counts     = bulk_counts,
  sc_counts_list  = sc_counts_list,
  gene_annotation = gene_annotation,
  chromosomes     = paste0("chr", 1:22),
  filter_biotype  = TRUE
)
genes <- gene_result$genes
cat("Selected", length(genes), "genes\n")


cat("\n=== Test: Preprocess Bulk ===\n")
bulk_df <- preprocess_bulk(
  bulk_counts = bulk_counts,
  genes = genes,
  output_file = bulk_file
)
cat("Bulk dimensions:", nrow(bulk_df), "×", ncol(bulk_df), "\n")


cat("\n=== Test: Obtain suggested number of cluster ===\n")
k_est <- estimate_clusterNumber(
  sc_counts = sc_counts_list[[1]],
  bulk_df = bulk_df,
  genes = bulk_df$Gene,
  quantile_norm = 0.99,
  seed = 1235
)

cat("\n=== Test: Preprocess Reference ===\n")
ref_df <- preprocess_reference(
    sc_counts_list = sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,  # Use all genes from bulk_df
    cluster = "k-means",
    n_clusters = k_est$k,
    output_file = ref_file
  )



cat("\n=== Test 6: Run CNNreg ===\n")
if (check_cnnreg_install()) {
  result <- run_cnnreg(
    bulk_file = bulk_file,
    ref_file = ref_file,
    output_dir = output_dir,
    n_celltypes = length(sc_counts_list[[1]]),
    epochs = 300,  # Short run for testing
    prefix = "test"
  )

  
  if (file.exists(result$proportions_file)) {
    props <- read_cnnreg_proportions(result$proportions_file)
    print(props)
  }
} else {
  cat("Skipping CNNreg test (not installed)\n")
}

cat("\n=== All Tests Complete ===\n")


unlink(work_dir, recursive = TRUE)
                                                        

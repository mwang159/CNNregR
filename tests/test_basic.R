library(CNNregR)

# Simulate example data for testing
set.seed(123)

# Create bulk RNA-seq data (10 samples × 1000 genes)
bulk_counts <- matrix(
  rpois(10 * 1000, lambda = 50),
  nrow = 10,
  ncol = 1000
)
colnames(bulk_counts) <- paste0("Gene_", 1:1000)
rownames(bulk_counts) <- paste0("Sample_", 1:10)

# Create scRNA-seq data with cell-type-specific patterns (3 cell types)
# Make first 100 genes differentially expressed to create distinct patterns
sc_counts_list <- list(
  CellType_A = cbind(
    matrix(rpois(200 * 100, lambda = 100), nrow = 200, ncol = 100),  # High in A
    matrix(rpois(200 * 900, lambda = 30), nrow = 200, ncol = 900)
  ),
  CellType_B = cbind(
    matrix(rpois(300 * 100, lambda = 30), nrow = 300, ncol = 100),
    matrix(rpois(300 * 900, lambda = 100), nrow = 300, ncol = 900)   # High in B
  ),
  CellType_C = matrix(rpois(150 * 1000, lambda = 60), nrow = 150, ncol = 1000)  # Medium everywhere
)

# Add gene names
for (ct in names(sc_counts_list)) {
  colnames(sc_counts_list[[ct]]) <- paste0("Gene_", 1:1000)
}

# Test 1: Check CNNreg installation
cat("\n=== Test 1: Check CNNreg Installation ===\n")
check_cnnreg_install()

# Test 2: Validate input data
cat("\n=== Test 2: Validate Input Data ===\n")
validation <- validate_input_data(bulk_counts, sc_counts_list)

# Test 3: Select marker genes
cat("\n=== Test 3: Select Marker Genes ===\n")
gene_result <- select_genes(
  bulk_counts = bulk_counts,
  sc_counts_list = sc_counts_list,
  filter_biotype = FALSE  # Skip biotype filtering for test data
)
genes <- gene_result$genes
size_factors <- gene_result$size_factors
cat("Selected", length(genes), "genes\n")

# Test 4: Preprocess bulk
cat("\n=== Test 4: Preprocess Bulk ===\n")
bulk_df <- preprocess_bulk(
  bulk_counts = bulk_counts,
  genes = genes,
  size_factors = size_factors,
  output_file = "test_bulk.csv"
)
cat("Bulk dimensions:", nrow(bulk_df), "×", ncol(bulk_df), "\n")

# Test 5: Preprocess reference
cat("\n=== Test 5: Preprocess Reference ===\n")
ref_df <- preprocess_reference(
    sc_counts_list = sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,  # Use all genes from bulk_df
    n_clusters = 10,
    output_file = "test_ref.csv"
  )
cat("Reference dimensions:", nrow(ref_df), "×", ncol(ref_df), "\n")

# Test 6: Run CNNreg (if installed)
cat("\n=== Test 6: Run CNNreg ===\n")
if (check_cnnreg_install()) {
  result <- run_cnnreg(
    bulk_file = "test_bulk.csv",
    ref_file = "test_ref.csv",
    output_dir = "test_output/",
    n_celltypes = 3,
    epochs = 100,  # Short run for testing
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

# Cleanup
unlink("test_bulk.csv")
unlink("test_ref.csv")
unlink("test_output", recursive = TRUE)

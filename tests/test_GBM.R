library(tidyverse)

options(stringsAsFactors=F)


###############################################################################
### read pseudoBulk 
### ll.bk$pseudoBulk is Sample X Gene count matrix, with colnames and rownames
ll.bk        <- readRDS("./data/pseudoBulk_GBM.RDS")
bulk_counts  <- ll.bk$pseudoBulk

###############################################################################
### senario 1:  read one single scRNAseq ref
sc        <- readRDS("./data/ref_GBM.RDS")
celltypes <- names(sc)  ## a list of count matrices (Sample X Gene) with colnames and rownames, each for one cell type 
sc_counts_list <- list(sc)

# Test: Check CNNreg installation
cat("\n=== Test: Check CNNreg Installation ===\n")
check_cnnreg_install()

# Test: Validate input data
cat("\n=== Test: Validate Input Data ===\n")
validation <- validate_input_data(bulk_counts, sc_counts_list)

# Test: Select genes
gene_annotation     <- read_delim("./data/gene_annotation.tsv", delim="\t")
cat("\n=== Test: Select Marker Genes ===\n")
gene_result <- select_genes(
  bulk_counts     = bulk_counts,
  sc_counts_list  = sc_counts_list,
  gene_annotation = gene_annotation,
  chromosomes     = paste0("chr", 1:22),
  filter_biotype  = T
)
genes <- gene_result$genes
cat("Selected", length(genes), "genes\n")


# Test: Preprocess bulk
cat("\n=== Test: Preprocess Bulk ===\n")
bulk_df <- preprocess_bulk(
  bulk_counts = bulk_counts,
  genes = genes,
  output_file = "./test/bulk.csv"
)
cat("Bulk dimensions:", nrow(bulk_df), "×", ncol(bulk_df), "\n")


# Test: Obtain suggested number of clusters
cat("\n=== Test: Obtain suggested number of cluster ===\n")
k_est  <-  estimate_clusterNumber(sc_counts_list[[1]],
                                 bulk_df,
                                 genes,
                                 quantile_norm = 0.99,
                                 seed = 1235)

# Test: Preprocess reference
cat("\n=== Test: Preprocess Reference ===\n")
ref_df <- preprocess_reference(
    sc_counts_list = sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,  # Use all genes from bulk_df
    cluster = "k-means",
    n_clusters = k_est$k,
    output_file = "./test/ref.csv" 
  )



# Test: Run CNNreg (if installed)
cat("\n=== Test 6: Run CNNreg ===\n")
if (check_cnnreg_install()) {
  result <- run_cnnreg(
    bulk_file   = "./test/bulk.csv",
    ref_file    = "./test/ref.csv",
    output_dir  = "./test/CNNreg_output/",
    n_celltypes = 7,
    epochs = 10000,  # Short run for testing
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
                                                        

library(testthat)
library(CNNregR)

make_test_data <- function() {
  set.seed(123)

  bulk_counts <- matrix(
    rpois(10 * 1000, lambda = 50),
    nrow = 10,
    ncol = 1000,
    dimnames = list(
      paste0("Sample_", 1:10),
      paste0("Gene_", 1:1000)
    )
  )

  sc_reference <- list(
    CellType_A = cbind(
      matrix(rpois(200 * 100, lambda = 100), nrow = 200, ncol = 100),
      matrix(rpois(200 * 900, lambda = 30), nrow = 200, ncol = 900)
    ),
    CellType_B = cbind(
      matrix(rpois(300 * 100, lambda = 30), nrow = 300, ncol = 100),
      matrix(rpois(300 * 900, lambda = 100), nrow = 300, ncol = 900)
    ),
    CellType_C = matrix(
      rpois(150 * 1000, lambda = 60),
      nrow = 150,
      ncol = 1000
    )
  )

  for (cell_type in names(sc_reference)) {
    colnames(sc_reference[[cell_type]]) <- paste0("Gene_", 1:1000)
  }

  list(
    bulk_counts = bulk_counts,
    sc_counts_list = list(sc_reference)
  )
}

test_that("synthetic input data pass validation", {
  input <- make_test_data()

  validation <- validate_input_data(
    input$bulk_counts,
    input$sc_counts_list
  )

  expect_true(validation$valid)
  expect_length(input$sc_counts_list, 1)
  expect_length(input$sc_counts_list[[1]], 3)
})

test_that("synthetic data can be preprocessed into a reference", {
  input <- make_test_data()
  output_dir <- tempfile("CNNregR-preprocess-")
  dir.create(output_dir)

  gene_result <- select_genes(
    bulk_counts = input$bulk_counts,
    sc_counts_list = input$sc_counts_list,
    filter_biotype = FALSE
  )

  expect_gt(length(gene_result$genes), 0)

  bulk_df <- preprocess_bulk(
    bulk_counts = input$bulk_counts,
    genes = gene_result$genes,
    output_file = file.path(output_dir, "bulk.csv")
  )

  expect_true(file.exists(file.path(output_dir, "bulk.csv")))
  expect_gt(nrow(bulk_df), 0)

  ref_df <- preprocess_reference(
    sc_counts_list = input$sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,
    cluster = "k-means",
    n_clusters = 10,
    output_file = file.path(output_dir, "reference.csv")
  )

  expect_true(file.exists(file.path(output_dir, "reference.csv")))
  expect_equal(nrow(ref_df), 10)
})

test_that("CNNreg can run the synthetic workflow", {
  skip_if_not(check_cnnreg_install())

  input <- make_test_data()
  output_dir <- tempfile("CNNregR-cnnreg-")
  dir.create(output_dir)

  gene_result <- select_genes(
    bulk_counts = input$bulk_counts,
    sc_counts_list = input$sc_counts_list,
    filter_biotype = FALSE
  )
  bulk_df <- preprocess_bulk(
    bulk_counts = input$bulk_counts,
    genes = gene_result$genes,
    output_file = file.path(output_dir, "bulk.csv")
  )
  preprocess_reference(
    sc_counts_list = input$sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,
    cluster = "k-means",
    n_clusters = 10,
    output_file = file.path(output_dir, "reference.csv")
  )

  result <- run_cnnreg(
    bulk_file = file.path(output_dir, "bulk.csv"),
    ref_file = file.path(output_dir, "reference.csv"),
    output_dir = file.path(output_dir, "cnnreg_output"),
    n_celltypes = length(input$sc_counts_list[[1]]),
    epochs = 100,
    prefix = "test"
  )

  expect_equal(result$exit_code, 0)
  expect_true(file.exists(result$proportions_file))
})

#' Run CNNreg Deconvolution via Command Line
#'
#' @description
#' Calls the CNNreg Python CLI tool from R using system() command.
#' Requires CNNreg to be installed in the Python environment.
#'
#' @param bulk_file Path to preprocessed bulk RNA-seq CSV file
#' @param ref_file Path to preprocessed reference CSV file
#' @param output_dir Output directory for results
#' @param n_celltypes Number of cell types (kernel size)
#' @param epochs Maximum training epochs (default: 50000)
#' @param prefix Output file prefix (default: "deconv")
#' @param mode CNNreg mode: "train", "evaluate", "predict", "explain" (default: "train")
#' @param python_path Path to Python executable with CNNreg installed (default: "python")
#' @param additional_args Additional arguments to pass to CNNreg (character vector)
#'
#' @return List with:
#'   \item{exit_code}{System exit code (0 = success)}
#'   \item{output_dir}{Path to output directory}
#'   \item{command}{Full command executed}
#'   \item{proportions_file}{Path to predicted cell proportions CSV}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_cnnreg(
#'   bulk_file = "bulk.csv",
#'   ref_file = "ref.csv",
#'   output_dir = "output/",
#'   n_celltypes = 7,
#'   epochs = 10000,
#'   prefix = "deconv"
#' )
#' }
run_cnnreg <- function(bulk_file,
                       ref_file,
                       output_dir,
                       n_celltypes,
                       epochs = 50000,
                       prefix = "deconv",
                       mode = "train",
                       python_path = "python",
                       additional_args = NULL) {
  
  # Check if cnnreg is available
  check_cmd <- sprintf("%s -m CNNreg.cli --help", python_path)
  test_result <- system(check_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  if (test_result != 0) {
    # Try direct cnnreg command
    check_cmd_2 <- "cnnreg --help"
    test_result_2 <- system(check_cmd_2, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (test_result_2 != 0) {
      stop("CNNreg not found. Please install: pip install CNNreg\n",
           "Or specify correct python_path with CNNreg installed.")
    }
    cnnreg_cmd <- "cnnreg"
  } else {
    cnnreg_cmd <- sprintf("%s -m CNNreg.cli", python_path)
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # Build command
  cmd <- sprintf(
    "%s -M %s -bulk %s -ref %s -o %s -C %d -EP %d -pre %s",
    cnnreg_cmd,
    mode,
    bulk_file,
    ref_file,
    output_dir,
    n_celltypes,
    epochs,
    prefix
  )
  
  # Add additional arguments
  if (!is.null(additional_args)) {
    cmd <- paste(cmd, paste(additional_args, collapse = " "))
  }
  
  # Display command
  message("Running CNNreg:")
  message(cmd)
  message("")
  
  # Execute
  exit_code <- system(cmd)
  
  if (exit_code != 0) {
    warning("CNNreg exited with code ", exit_code)
  } else {
    message("\nCNNreg completed successfully!")
  }
  
  # Construct output file path
  prop_file <- file.path(output_dir, paste0("Prop_predicted_", prefix, ".csv"))
  
  return(list(
    exit_code = exit_code,
    output_dir = output_dir,
    command = cmd,
    proportions_file = prop_file
  ))
}


#' High-Level Deconvolution Wrapper
#'
#' @description
#' Complete workflow: preprocess data and run CNNreg deconvolution using
#' edgeR TMM normalization and k-means clustering
#'
#' @param bulk_counts Matrix or data.frame of bulk RNA-seq counts (samples × genes)
#' @param sc_counts_list Named list of scRNA-seq count matrices per cell type
#' @param output_dir Output directory for all results
#' @param prefix Output file prefix (default: "deconv")
#' @param epochs Maximum training epochs (default: 50000)
#' @param n_clusters Number of k-means clusters (references) per cell type (default: 5)
#' @param filter_biotype Logical, whether to filter by gene biotype (default: TRUE)
#' @param gene_biotypes Gene biotypes to include if filter_biotype=TRUE (default: protein_coding, lincRNA, processed_pseudogene)
#' @param organism Organism name for BioMart (default: "hsapiens")
#' @param gene_annotation Optional gene annotation data.frame
#' @param chromosomes Character vector of chromosome names to include (default: NULL)
#' @param exclude_genes Character vector of genes to exclude (default: MALAT1)
#' @param bulk_min_expr Minimum median expression in bulk (default: 1)
#' @param bulk_min_median_quantile Minimum quantile for bulk median expression filtering (default: 0.1)
#' @param bulk_max_median_quantile Maximum quantile for bulk median expression filtering (default: 0.999)
#' @param bulk_min_cv_quantile Minimum quantile for bulk CV filtering (default: 0.25)
#' @param bulk_max_cv_quantile Maximum quantile for bulk CV filtering (default: 1.0)
#' @param bulk_min_cv_final Minimum CV threshold after final normalization (default: 0.0001)
#' @param ref_min_detection Minimum fraction of cells expressing gene per cell type (default: 0.2)
#' @param ref_max_housekeeping_factor Exclude genes with scRNA-seq expression > factor × median (default: 50)
#' @param ref_var_quantile scRNA-seq variance quantile threshold (default: 0.5)
#' @param python_path Path to Python with CNNreg (default: "python")
#' @param keep_intermediates Keep intermediate CSV files (default: TRUE)
#'
#' @return List with:
#'   \item{genes}{Selected marker genes}
#'   \item{size_factors}{TMM normalization factors}
#'   \item{bulk_file}{Path to bulk CSV}
#'   \item{ref_file}{Path to reference CSV}
#'   \item{proportions}{Data.frame of predicted cell proportions}
#'   \item{cnnreg_result}{Output from run_cnnreg()}
#'
#' @export
#'
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' # Human samples with biotype filtering
#' result <- deconvolve(
#'   bulk_counts = bulk_data,
#'   sc_counts_list = list(Astro = astro_counts, Oligo = oligo_counts),
#'   output_dir = "results/",
#'   prefix = "deconv",
#'   filter_biotype = TRUE,
#'   organism = "hsapiens"
#' )
#' }
deconvolve <- function(bulk_counts,
                       sc_counts_list,
                       output_dir,
                       prefix = "deconv",
                       epochs = 50000,
                       n_clusters = 5,
                       filter_biotype = TRUE,
                       gene_biotypes = c("protein_coding", "lincRNA", "processed_pseudogene"),
                       organism = "hsapiens",
                       gene_annotation = NULL,
                       chromosomes = NULL,
                       exclude_genes = c("MALAT1"),
                       bulk_min_expr = 1,
                       bulk_min_median_quantile = 0.1,
                       bulk_max_median_quantile = 0.999,
                       bulk_min_cv_quantile = 0.25,
                       bulk_max_cv_quantile = 1.0,
                       bulk_min_cv_final = 0.0001,
                       ref_min_detection = 0.2,
                       ref_max_housekeeping_factor = 50,
                       ref_var_quantile = 0.5,
                       python_path = "python",
                       keep_intermediates = TRUE) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("=== CNNreg Deconvolution Workflow (edgeR TMM + K-means) ===")
  if (filter_biotype) {
    message("Organism: ", organism)
    message("Gene biotypes: ", paste(gene_biotypes, collapse = ", "))
  } else {
    message("Biotype filtering: disabled")
  }
  message("")
  
  # Step 1: Select marker genes with TMM normalization
  message("Step 1/4: Selecting marker genes with edgeR TMM...")
  gene_result <- select_genes(
    bulk_counts = bulk_counts,
    sc_counts_list = sc_counts_list,
    filter_biotype = filter_biotype,
    gene_biotypes = gene_biotypes,
    organism = organism,
    gene_annotation = gene_annotation,
    chromosomes = chromosomes,
    exclude_genes = exclude_genes,
    bulk_min_expr = bulk_min_expr,
    bulk_min_median_quantile = bulk_min_median_quantile,
    bulk_max_median_quantile = bulk_max_median_quantile,
    bulk_min_cv_quantile = bulk_min_cv_quantile,
    bulk_max_cv_quantile = bulk_max_cv_quantile,
    ref_min_detection = ref_min_detection,
    ref_max_housekeeping_factor = ref_max_housekeeping_factor,
    ref_var_quantile = ref_var_quantile
  )
  genes <- gene_result$genes
  size_factors <- gene_result$size_factors
  message("")
  
  # Step 2: Preprocess bulk with size factors
  message("Step 2/4: Preprocessing bulk RNA-seq with TMM normalization...")
  bulk_file <- file.path(output_dir, paste0("bulk_", prefix, ".csv"))
  bulk_df <- preprocess_bulk(
    bulk_counts = bulk_counts,
    genes = genes,
    size_factors = size_factors,
    bulk_min_cv_final = bulk_min_cv_final,
    output_file = bulk_file
  )
  message("")
  
  # Step 3: Preprocess reference with k-means
  message("Step 3/4: Generating references with k-means clustering...")
  ref_file <- file.path(output_dir, paste0("ref_", prefix, ".csv"))
  ref_df <- preprocess_reference(
    sc_counts_list = sc_counts_list,
    bulk_df = bulk_df,
    genes = bulk_df$Gene,  # Use genes after CV filter
    n_clusters = n_clusters,
    output_file = ref_file
  )
  message("")
  
  # Step 4: Run CNNreg
  message("Step 4/4: Running CNNreg deconvolution...")
  n_celltypes <- length(sc_counts_list)
  
  cnnreg_result <- run_cnnreg(
    bulk_file = bulk_file,
    ref_file = ref_file,
    output_dir = output_dir,
    n_celltypes = n_celltypes,
    epochs = epochs,
    prefix = prefix,
    python_path = python_path
  )
  
  # Load results
  proportions <- NULL
  if (file.exists(cnnreg_result$proportions_file)) {
    proportions <- readr::read_csv(cnnreg_result$proportions_file, show_col_types = FALSE)
    message("\nPredicted proportions:")
    print(head(proportions))
  } else {
    warning("Proportions file not found: ", cnnreg_result$proportions_file)
  }
  
  # Clean up intermediates
  if (!keep_intermediates) {
    file.remove(bulk_file, ref_file)
    message("\nRemoved intermediate files")
  }
  
  message("\n=== Deconvolution Complete ===")
  
  return(list(
    genes = genes,
    size_factors = size_factors,
    bulk_file = bulk_file,
    ref_file = ref_file,
    proportions = proportions,
    cnnreg_result = cnnreg_result
  ))
}

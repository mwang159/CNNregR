#' Get Gene Annotation from BioMart
#'
#' @description
#' Retrieves gene annotations (biotype, chromosome) from Ensembl BioMart
#' for specified organism
#'
#' @param genes Character vector of gene symbols
#' @param organism Organism name for BioMart. Options:
#'   - "hsapiens" (human, default)
#'   - "mmusculus" (mouse)
#'   Or full dataset name like "hsapiens_gene_ensembl"
#' @param attributes Attributes to retrieve (default: gene symbol, biotype, chromosome)
#' @param host BioMart host (default: "https://www.ensembl.org")
#'
#' @return Data.frame with gene annotations
#' @export
#'
#' @examples
#' \dontrun{
#' # Human genes
#' anno <- get_gene_annotation(
#'   genes = c("TP53", "GAPDH", "ACTB"),
#'   organism = "hsapiens"
#' )
#' 
#' # Mouse genes
#' anno_mouse <- get_gene_annotation(
#'   genes = c("Tp53", "Gapdh"),
#'   organism = "mmusculus"
#' )
#' }
get_gene_annotation <- function(genes,
                                 organism = "hsapiens",
                                 attributes = c("hgnc_symbol", "gene_biotype", "chromosome_name"),
                                 host = "https://www.ensembl.org") {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' required. Install with: BiocManager::install('biomaRt')")
  }
  
  # Map short names to full dataset names
  dataset_map <- list(
    "hsapiens" = "hsapiens_gene_ensembl",
    "mmusculus" = "mmusculus_gene_ensembl"
  )
  
  # Get dataset name
  if (organism %in% names(dataset_map)) {
    dataset <- dataset_map[[organism]]
  } else if (grepl("_gene_ensembl$", organism)) {
    dataset <- organism
  } else {
    stop("Unknown organism: ", organism, 
         "\nSupported: ", paste(names(dataset_map), collapse = ", "),
         "\nOr provide full dataset name (e.g., 'hsapiens_gene_ensembl')")
  }
  
  # Adjust gene symbol attribute based on organism
  symbol_attr <- switch(organism,
    "hsapiens" = "hgnc_symbol",
    "mmusculus" = "mgi_symbol",
    "external_gene_name"  # Default
  )
  
  # Update attributes
  if (attributes[1] == "hgnc_symbol") {
    attributes[1] <- symbol_attr
  }
  
  message("Querying BioMart for ", organism, "...")
  message("  Dataset: ", dataset)
  message("  Genes: ", length(genes))
  
  tryCatch({
    # Connect to BioMart
    ensembl <- biomaRt::useMart("ensembl", dataset = dataset, host = host)
    
    # Query annotations
    anno <- biomaRt::getBM(
      attributes = attributes,
      filters = symbol_attr,
      values = genes,
      mart = ensembl
    )
    
    # Rename columns
    colnames(anno)[1] <- "GeneName"
    if (length(attributes) >= 2) colnames(anno)[2] <- "GeneBiotype"
    if (length(attributes) >= 3) colnames(anno)[3] <- "Chromosome"
    
    # Remove duplicates
    anno <- anno[!duplicated(anno$GeneName), ]
    
    message("Retrieved annotations for ", nrow(anno), " genes")
    
    return(anno)
    
  }, error = function(e) {
    warning("BioMart query failed: ", e$message)
    message("Returning NULL. You can provide gene_annotation manually.")
    return(NULL)
  })
}


#' Check CNNreg Installation
#'
#' @description
#' Verifies that CNNreg Python package is installed and accessible
#'
#' @param python_path Path to Python executable (default: "python")
#'
#' @return Logical: TRUE if CNNreg is available, FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' if (check_cnnreg_install()) {
#'   message("CNNreg is ready!")
#' } else {
#'   message("Install CNNreg: pip install CNNreg")
#' }
#' }
check_cnnreg_install <- function(python_path = "python") {
  # Try python module
  check_cmd <- sprintf("%s -c 'import CNNreg; print(CNNreg.__version__)'", python_path)
  result <- suppressWarnings(system(check_cmd, intern = TRUE, ignore.stderr = TRUE))
  
  if (length(result) > 0) {
    message("CNNreg version: ", result[1])
    return(TRUE)
  }
  
  # Try cnnreg command
  result2 <- suppressWarnings(system("cnnreg --help", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if (result2 == 0) {
    message("CNNreg CLI found")
    return(TRUE)
  }
  
  message("CNNreg not found. Install with: pip install CNNreg")
  return(FALSE)
}


#' Read CNNreg Output Proportions
#'
#' @description
#' Helper function to read cell proportion predictions from CNNreg output
#'
#' @param file Path to Prop_predicted_*.csv file
#' @param rename_samples Optional named vector to rename samples
#'
#' @return Data.frame with cell proportions
#' @export
#'
#' @examples
#' \dontrun{
#' props <- read_cnnreg_proportions("output/Prop_predicted_GBM.csv")
#' }
read_cnnreg_proportions <- function(file, rename_samples = NULL) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  if (!is.null(rename_samples)) {
    if ("Sample" %in% colnames(df)) {
      df$Sample <- rename_samples[df$Sample]
    }
  }
  
  return(df)
}


#' Validate Input Data Format
#'
#' @description
#' Checks that input data meets CNNreg requirements
#'
#' @param bulk_counts Bulk RNA-seq count matrix
#' @param sc_counts_list List of scRNA-seq count matrices
#'
#' @return List with validation results and warnings
#' @export
validate_input_data <- function(bulk_counts, sc_counts_list) {
  issues <- list()
  
  # Check bulk
  if (!is.matrix(bulk_counts) && !is.data.frame(bulk_counts)) {
    issues <- c(issues, "bulk_counts must be matrix or data.frame")
  }
  
  if (any(is.na(bulk_counts))) {
    issues <- c(issues, "bulk_counts contains NA values")
  }
  
  if (any(bulk_counts < 0)) {
    issues <- c(issues, "bulk_counts contains negative values")
  }
  
  # Check scRNA-seq
  if (!is.list(sc_counts_list) || is.null(names(sc_counts_list))) {
    issues <- c(issues, "sc_counts_list must be a named list")
  }
  
  for (ct in names(sc_counts_list)) {
    mat <- sc_counts_list[[ct]]
    if (!is.matrix(mat) && !is.data.frame(mat)) {
      issues <- c(issues, paste0("sc_counts_list[[", ct, "]] must be matrix or data.frame"))
    }
    if (any(is.na(mat))) {
      issues <- c(issues, paste0("sc_counts_list[[", ct, "]] contains NA values"))
    }
    if (any(mat < 0)) {
      issues <- c(issues, paste0("sc_counts_list[[", ct, "]] contains negative values"))
    }
  }
  
  # Check overlap
  bulk_genes <- colnames(bulk_counts)
  sc_genes <- unique(unlist(lapply(sc_counts_list, colnames)))
  overlap <- intersect(bulk_genes, sc_genes)
  
  if (length(overlap) < 100) {
    issues <- c(issues, sprintf("Only %d genes overlap between bulk and scRNA-seq (need >100)", 
                                length(overlap)))
  }
  
  if (length(issues) > 0) {
    message("Validation issues found:")
    for (issue in issues) {
      message("  - ", issue)
    }
    return(list(valid = FALSE, issues = issues))
  }
  
  message("Input data validation passed!")
  message("  Bulk samples: ", nrow(bulk_counts))
  message("  Bulk genes: ", ncol(bulk_counts))
  message("  Cell types: ", length(sc_counts_list))
  message("  Total scRNA-seq cells: ", sum(sapply(sc_counts_list, nrow)))
  message("  Gene overlap: ", length(overlap))
  
  return(list(valid = TRUE, issues = NULL))
}


#' Format Seurat Object for CNNregR
#'
#' @description
#' Extracts count matrices from Seurat object by cell type
#'
#' @param seurat_obj Seurat object with cell type annotations
#' @param celltype_col Column name in metadata containing cell type labels
#' @param assay Assay to use (default: "RNA")
#' @param layer Layer to use: "counts" or "data" (default: "counts")
#'
#' @return Named list of count matrices per cell type
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires Seurat installed
#' sc_counts_list <- format_seurat_object(
#'   seurat_obj = sc_data,
#'   celltype_col = "celltype"
#' )
#' }
format_seurat_object <- function(seurat_obj, 
                                  celltype_col, 
                                  assay = "RNA",
                                  layer = "counts") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' required. Install with: install.packages('Seurat')")
  }
  
  # Get cell type labels
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", celltype_col, "' not found in Seurat metadata")
  }
  
  celltypes <- unique(seurat_obj@meta.data[[celltype_col]])
  celltypes <- celltypes[!is.na(celltypes)]
  
  # Extract counts per cell type
  counts_list <- list()
  for (ct in celltypes) {
    cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == ct]
    
    if (layer == "counts") {
      mat <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")[, cells])
    } else {
      mat <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data")[, cells])
    }
    
    counts_list[[ct]] <- t(mat)  # Transpose to cells × genes
  }
  
  message("Extracted ", length(counts_list), " cell types from Seurat object:")
  for (ct in names(counts_list)) {
    message("  ", ct, ": ", nrow(counts_list[[ct]]), " cells")
  }
  
  return(counts_list)
}

#' CNNregR: R Interface for CNNreg Cell Type Deconvolution
#'
#' @description
#' CNNregR provides R functions for preprocessing scRNA-seq and bulk RNA-seq
#' data for CNNreg deconvolution analysis using edgeR TMM normalization and
#' k-means clustering. It includes:
#' 
#' - Marker gene selection with edgeR TMM normalization and discrimination scoring
#' - Optional gene biotype filtering via BioMart for multiple organisms
#' - K-means clustering-based reference generation from scRNA-seq
#' - Data formatting for CNNreg compatibility
#' - Wrapper functions to run CNNreg from R
#' 
#' @section Main Functions:
#' 
#' **High-level workflow:**
#' \itemize{
#'   \item \code{\link{deconvolve}}: Complete preprocessing + deconvolution pipeline
#' }
#' 
#' **Preprocessing:**
#' \itemize{
#'   \item \code{\link{select_genes}}: Select discriminative genes with edgeR TMM
#'   \item \code{\link{preprocess_bulk}}: Format bulk RNA-seq data
#'   \item \code{\link{preprocess_reference}}: Generate reference CSE profiles via k-means
#' }
#' 
#' **Deconvolution:**
#' \itemize{
#'   \item \code{\link{run_cnnreg}}: Call CNNreg CLI from R
#' }
#' 
#' **Utilities:**
#' \itemize{
#'   \item \code{\link{check_cnnreg_install}}: Verify CNNreg installation
#'   \item \code{\link{validate_input_data}}: Check data format
#'   \item \code{\link{format_seurat_object}}: Convert Seurat to CNNregR format
#'   \item \code{\link{get_gene_annotation}}: Query BioMart for gene annotations
#' }
#'
#' @name CNNregR-package
#' @aliases CNNregR
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom matrixStats colMedians
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter pull
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom stats cor median quantile var kmeans
#' @importFrom utils head
#' @importFrom magrittr %>%
## usethis namespace: end

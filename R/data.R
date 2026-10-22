#' Example GBM pseudobulk expression data
#'
#' A pseudobulk RNA-seq expression matrix used in CNNregR examples
#' and tests.
#'
#' @format A list of length 2. First element is matrix containing pseudobulk expression data with samples in rows and genes in columns. The second element is a matrix containing the corresponding cell-type proportions for each sample.
#' @source Internal example dataset.
#' @usage data(pseudoBulk_GBM)
"pseudoBulk_GBM"

#' Example GBM reference expression data
#'
#' A cell-type reference expression dataset used by CNNregR examples
#' and tests.
#'
#' @format A list of length 7. Each element is a matrix containing expression data for a specific cell type.
#' @source Internal example dataset.
#' @usage data(ref_GBM)
"ref_GBM"

#' Example gene annotation data
#'
#' A gene annotation dataset used by CNNregR examples
#' and tests.
#'
#' @format A data.frame with columns for gene ID, gene name, chromosome, start position, end position, and biotype.
#' @source Internal example dataset.
#' @usage data(gene_annotation)
"gene_annotation"

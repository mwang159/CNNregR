#' Select Marker Genes for Deconvolution (edgeR-based)
#'
#' @description
#' Selects discriminative genes suitable for cell type deconvolution using
#' edgeR TMM normalization and coefficient of variation filtering:
#' 1. Optional gene biotype filtering via BioMart
#' 2. Excludes MT and ribosomal genes
#' 3. TMM normalization with edgeR
#' 4. CV and expression filtering
#' 5. Variance-based discrimination from scRNA-seq
#'
#' @param bulk_counts Matrix or data.frame of bulk RNA-seq counts (samples × genes)
#' @param sc_counts_list Named list of scRNA-seq count matrices per cell type (cells × genes)
#' @param filter_biotype Logical, whether to filter by gene biotype (default: TRUE)
#' @param gene_biotypes Gene biotypes to include if filter_biotype=TRUE 
#'   (default: c("protein_coding", "lincRNA", "processed_pseudogene"))
#' @param organism Organism name for BioMart queries if filter_biotype=TRUE
#'   Options: "hsapiens" (default), "mmusculus", or full dataset name like "hsapiens_gene_ensembl"
#' @param gene_annotation Optional data.frame with gene annotations (GeneName, GeneBiotype, Chromosome).
#'   If NULL and filter_biotype=TRUE, will query BioMart.
#' @param chromosomes Character vector of chromosome names to include (e.g., c("1", "2", "X")).
#'   If NULL (default), all chromosomes are included. Only used if gene_annotation has Chromosome column.
#' @param exclude_genes Character vector of genes to exclude (default: MALAT1).
#' @param bulk_min_expr Minimum median expression in bulk (default: 1)
#' @param bulk_min_median_quantile Minimum quantile for bulk median expression filtering (default: 0.1)
#' @param bulk_max_median_quantile Maximum quantile for bulk median expression filtering (default: 0.999)
#' @param bulk_min_cv_quantile Minimum quantile for bulk coefficient of variation filtering (default: 0.25)
#' @param bulk_max_cv_quantile Maximum quantile for bulk coefficient of variation filtering (default: 1.0)
#' @param ref_min_detection Minimum fraction of cells expressing gene per cell type in scRNA-seq (default: 0.2)
#' @param ref_max_housekeeping_factor Exclude genes with scRNA-seq expression > factor × median (default: 50)
#' @param ref_var_quantile scRNA-seq variance quantile threshold for discrimination (default: 0.5) 
#'
#' @return List with:
#'   \item{genes}{Character vector of selected gene names}
#'   \item{gene_stats}{Data.frame with gene statistics}
#'   \item{gene_annotation}{Data.frame with gene annotations (if queried)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Human genes with biotype filtering
#' result <- select_genes(
#'   bulk_counts = bulk_data,
#'   sc_counts_list = list(Astro = astro_counts, Oligo = oligo_counts),
#'   filter_biotype = TRUE,
#'   organism = "hsapiens"
#' )
#' 
#' # Skip biotype filtering (any organism)
#' result <- select_genes(
#'   bulk_counts = bulk_data,
#'   sc_counts_list = sc_list,
#'   filter_biotype = FALSE
#' )
#' 
#' # Custom biotypes for mouse
#' result <- select_genes(
#'   bulk_counts = bulk_data,
#'   sc_counts_list = sc_list,
#'   filter_biotype = TRUE,
#'   gene_biotypes = c("protein_coding", "lincRNA"),
#'   organism = "mmusculus"
#' )
#' }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median quantile var kmeans
select_genes <- function(bulk_counts,
                                 sc_counts_list,
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
                                 ref_min_detection = 0.2,
                                 ref_max_housekeeping_factor = 50,
                                 ref_var_quantile = 0.5) {
  
  # Step 1: Calculate TMM size factors using edgeR
  message("Calculating TMM normalization factors...")
  tmm <- edgeR::DGEList(t(bulk_counts))
  tmm <- edgeR::calcNormFactors(tmm)

  
  # Step 2: Filter bulk genes by annotation
  genes_bulk <- colnames(bulk_counts)
  anno_result <- NULL
  
  if (filter_biotype) {
    message("Filtering by gene biotype...")
    
    # Get or query gene annotation
    if (is.null(gene_annotation)) {
      message("Querying BioMart for ", organism, " gene annotations...")
      gene_annotation <- get_gene_annotation(
        genes = genes_bulk,
        organism = organism
      )
      anno_result <- gene_annotation
    }
    
    if (!is.null(gene_annotation)) {
      # Filter by biotype and optionally by chromosome
      if ("Chromosome" %in% colnames(gene_annotation) && !is.null(chromosomes)) {
        genes_incl <- gene_annotation %>%
          dplyr::filter(
            .data$Chromosome %in% chromosomes,
            .data$GeneBiotype %in% gene_biotypes
          ) %>%
          dplyr::pull(.data$GeneName)
        message("  Filtering by chromosomes: ", paste(chromosomes, collapse = ", "))
      } else {
        genes_incl <- gene_annotation %>%
          dplyr::filter(.data$GeneBiotype %in% gene_biotypes) %>%
          dplyr::pull(.data$GeneName)
      }
      
      genes_bulk <- intersect(genes_bulk, genes_incl)
      message("  Kept ", length(genes_bulk), " genes with biotypes: ", 
              paste(gene_biotypes, collapse = ", "))
    } else {
      warning("BioMart query failed. Proceeding without biotype filtering.")
      filter_biotype <- FALSE
    }
  }
  
  # Exclude MALAT1, MT genes, ribosomal genes
  genes_bulk <- setdiff(genes_bulk, exclude_genes)
  genes_bulk <- genes_bulk[!grepl("^MT-|^MT|^mt-|^RPL|^RPS|^Rpl|^Rps", genes_bulk)]
  
  message("After MT/ribosomal exclusion: ", length(genes_bulk), " genes")
  
  bulk_sub <- bulk_counts[, genes_bulk]
  
  # Step 3: Calculate normalized expression and CV
  effective_lib_size <- tmm$samples$lib.size * tmm$samples$norm.factors
  bulk_norm      <- sweep(bulk_sub, 1, effective_lib_size, "/") * 1e6
  median_bulk <- matrixStats::colMedians(as.matrix(bulk_sub))
  names(median_bulk) <- genes_bulk
  
  # Filter by minimum expression
  keep_genes <- genes_bulk[median_bulk >= bulk_min_expr]
  bulk_sub <- bulk_sub[, keep_genes]
  bulk_norm <- bulk_norm[, keep_genes]
  median_bulk <- median_bulk[keep_genes]
  
  # Calculate CV
  cv_bulk <- apply(bulk_norm, 2, function(x) {
    var(x) / (mean(x) + 1)
  })
  
  # Step 4: Apply quantile filters
  genes_bulk <- colnames(bulk_sub)[
    median_bulk >= quantile(median_bulk, probs = bulk_min_median_quantile) &
    median_bulk <= quantile(median_bulk, probs = bulk_max_median_quantile) &
    cv_bulk >= quantile(cv_bulk, probs = bulk_min_cv_quantile) &
    cv_bulk <= quantile(cv_bulk, probs = bulk_max_cv_quantile)
  ]
  
  message("After expression/CV filters: ", length(genes_bulk), " genes")
  
  # Step 5: Filter scRNA-seq genes
  celltypes <- names(sc_counts_list[[1]])
  genes     <- colnames(sc_counts_list[[1]][[1]])
  # Select cell-type-specific genes
  for(ii in 1:length(sc_counts_list)){
      sc_counts <- sc_counts_list[[ii]]
      gene_spe  <- gene_excl <- vector("list", length(celltypes))
      names(gene_spe) <- names(gene_excl) <- celltypes
      for (ct in celltypes) {
          x <- sc_counts[[ct]]
          detection_rate <- colSums(x > 0) / nrow(x)
          genes_detected <- colnames(x)[detection_rate >= ref_min_detection]
          sum_expr <- colSums(x[, genes_detected])
          med_expr <- median(sum_expr)
          gene_excl[[ct]] <- genes_detected[sum_expr >= ref_max_housekeeping_factor * med_expr]
          gene_spe[[ct]] <- setdiff(genes_detected, gene_excl[[ct]])
      }
      gene_sc  <- setdiff(Reduce("union", gene_spe), Reduce("union", gene_excl))
      genes    <- intersect(genes, gene_sc)
  }
  # Intersect bulk and scRNA genes
  genes <- intersect(genes, colnames(genes_bulk))
  

  
  message("After scRNA-seq filtering: ", length(genes), " genes")
  
  # Step 6: Calculate discrimination score from scRNA-seq
  gene_sel <- list(length(sc_counts_list))
  for(ii in 1:length(sc_counts_list)){
      sum_celltype <- lapply(celltypes, function(ct) {
        colSums(sc_counts_list[[ct]][, genes])
      })
  
      norm_celltype <- lapply(sum_celltype, function(x) {
        log2(1 + 1e6 * x / sum(x))
      })
  
      expr_matrix <- do.call(rbind, norm_celltype)
  
      dis_score <- apply(expr_matrix, 2, function(x) {
        var(x)
      })
  
      # Keep top variance genes
      threshold <- quantile(dis_score, ref_var_quantile)
      gene_sel[[ii]] <- genes[dis_score >= threshold]
  }
  genes <- Reduce("union", gene_sel)
  
  # Gene statistics
  gene_stats <- data.frame(
    gene = genes,
    median_bulk = median_bulk[genes],
    cv_bulk = cv_bulk[genes],
    stringsAsFactors = FALSE
  )
  
  message("\n=== Final Gene Selection ===")
  message("Selected ", length(genes), " marker genes")
  message("  Median expression: ", round(min(median_bulk[genes]), 2), " - ", round(max(median_bulk[genes]), 2))
  message("  CV: ", round(min(cv_bulk[genes]), 3), " - ", round(max(cv_bulk[genes]), 3))
  
  result <- list(
    genes = genes,
    gene_stats = gene_stats
  )
  
  if (!is.null(anno_result)) {
    result$gene_annotation <- anno_result
  }
  
  return(result)
}


#' Preprocess Bulk RNA-seq Data (edgeR TMM)
#'
#' @description
#' Normalizes bulk RNA-seq data using TMM size factors and quantile normalization
#'
#' @param bulk_counts Matrix or data.frame of bulk RNA-seq counts (samples × genes)
#' @param genes Character vector of genes to include
#' @param quantile_norm Quantile for normalization (default: 0.99)
#' @param bulk_min_cv_final Minimum CV threshold after final normalization (default: 0.0001)
#' @param output_file Path to save CSV file (optional)
#'
#' @return Data.frame with format: Gene, Sample1, Sample2, ... (transposed)
#' @export
#'
#' @examples
#' \dontrun{
#' result <- select_genes(bulk_data, sc_list)
#' bulk_df <- preprocess_bulk(
#'   bulk_counts = bulk_data,
#'   genes = result$genes,
#'   output_file = "bulk_processed.csv"
#' )
#' }
preprocess_bulk <- function(bulk_counts,
                             genes,
                             quantile_norm = 0.99,
                             bulk_min_cv_final = 0.0001,
                             output_file = NULL) {
  
  #Calculate TMM size factors using edgeR
  message("Calculating TMM normalization factors...")
  tmm <- edgeR::DGEList(t(bulk_counts))
  tmm <- edgeR::calcNormFactors(tmm)

  # Subset genes
  expr <- bulk_norm[, genes] + 0.00001
  
  # Transpose to genes × samples
  expr <- t(expr)
  
  # Quantile normalization across samples
  quant <- apply(expr, 2, function(x) {
    quantile(x, quantile_norm)
  })
  expr_norm <- expr / median(quant)
  
  # Additional CV filtering
  cv <- apply(expr_norm, 1, function(x) {
    var(x) / mean(x)
  })
  
  keep_genes <- genes[cv >= bulk_min_cv_final]
  expr_norm <- expr_norm[keep_genes, ]
  
  # Format: Gene column + sample columns (transposed)
  df_out <- data.frame(Gene = keep_genes, expr_norm, check.names = FALSE)
  
  # Write output
  if (!is.null(output_file)) {
    readr::write_csv(df_out, output_file)
    message("Bulk data saved to: ", output_file)
  }
  
  message("Final genes after CV filter: ", length(keep_genes))
  
  return(df_out)
}


#' Preprocess Reference scRNA-seq Data (K-means Clustering)
#'
#' @description
#' Creates reference Cell-type Specific Expression (CSE) profiles using
#' k-means clustering. Generates reference samples by clustering cells within
#' each cell type and matching expression scales to bulk samples.
#'
#' @param sc_counts_list Named list of scRNA-seq count matrices per cell type (cells × genes)
#' @param bulk_df Preprocessed bulk data.frame (output from preprocess_bulk)
#' @param genes Character vector of genes to include. Should be bulk_df$Gene from preprocess_bulk output
#'   to ensure genes match between bulk and reference after final CV filtering.
#' @param cluster should use the reference ("byRef") or divide cells by k-mean clustering ("k-mean")
#' @param n_clusters if we do k-mean clustering, the number of clusters per cell type (default: 1)
#' @param quantile_norm Quantile for normalization (default: 0.99)
#' @param seed Random seed for k-means (default: 1235)
#' @param output_file Path to save CSV file (optional)
#'
#' @return Data.frame with format: Sample, Gene1_CellType1, Gene1_CellType2, ...
#' @export
#'
#' @examples
#' \dontrun{
#' # Complete workflow
#' gene_result <- select_genes(bulk_data, sc_list)
#' bulk_df <- preprocess_bulk(bulk_data, gene_result$genes, gene_result$size_factors)
#' 
#' # Use genes from bulk_df output (after final CV filtering)
#' ref_df <- preprocess_reference(
#'   sc_counts_list = sc_list,
#'   bulk_df = bulk_df,
#'   genes = bulk_df$Gene,
#'   n_clusters = 5,
#'   output_file = "ref_processed.csv"
#' )
#' }
preprocess_reference <- function(sc_counts_list,
                                  bulk_df,
                                  genes,
                                  cluster="byRef",
                                  n_clusters = 1,
                                  quantile_norm = 0.99,
                                  seed = 1235,
                                  output_file = NULL) {
  
  celltypes <- names(sc_counts_list[[1]])
  mm.bk       <- as.matrix(bulk_df[,2:ncol(bulk_df)])
  quant.bk    <- quantile(apply(mm.bk, 2, function(x){quantile(x, 0.75)}), 0.5)
  quant.mx    <- quantile(apply(mm.bk, 2, function(x){quantile(x, 1)}), 0.5)
  ll.multi    <- list()
  if(cluster == "byRef"){ N.cl <- length(sc_counts_list)}

  ## to prepare the reference by input sn/scRNAseq datasets
  if(cluster == "byRef"){      
        for(ii in 1:length(sc_counts_list)){
            ll.sc <- sc_counts_list[[ii]]
            ll.sort <- vector("list", length(celltypes))
            names(ll.sort) <- celltypes
            for(ss in celltypes){
                mm             <- ll.sc[[ss]][, genes, drop=F]
                ll.sort[[ss]]  <- colMeans(mm)
                quant          <- quantile(ll.sort[[ss]], quantile_norm)
                ll.sort[[ss]]  <- ll.sort[[ss]]/quant
                quant.sc       <- quantile(ll.sort[[ss]], 0.75)
                s.factor       <- quant.bk/quant.sc
                print(s.factor)
                s.factor       <- min(s.factor, 5)
                s.factor       <- max(s.factor, 0.2)
                idx.large      <- which(ll.sort[[ss]] > 1)
                idx.small      <- which(ll.sort[[ss]] <= 1)
                ll.sort[[ss]]  <- s.factor * ll.sort[[ss]]
                tmp            <- ll.sort[[ss]]
                if(length(idx.large)>=2){
                    tmp.idx        <- ll.sort[[ss]][idx.large]
                    tmp.min        <- min(tmp[idx.large])
                    tmp.max        <- max(tmp[idx.large])
                    vv             <- quant.mx*(tmp.idx-tmp.min)/(tmp.max-tmp.min)
                    ll.sort[[ss]][idx.large] <- s.factor + vv
                }}
            cc <- c()
            for(gene in genes){
                cc <- c(cc, sapply(celltypes, function(x){ ll.sort[[x]][gene] })) }
            ll.multi[[ii]] <- cc
      }
  }
  # K-means clustering for each cell type
  set.seed(seed)
  if( cluster=="k-mean" && length(N.cl) == 1 && !is.na(N.cl) && N.cl > 0){
      message("Performing k-means clustering...")
      cl_list <- list()
      for (ct in celltypes) {
        message("  Clustering ", ct, "...")
        sc_matrix <- sc_counts_list[[ct]][, genes]
        cl_list[[ct]] <- kmeans(log2(1 + sc_matrix), centers = n_clusters, nstart = 10)
      }
  
     # Calculate bulk quantiles for scaling
      bulk_matrix <- as.matrix(bulk_df[, 2:ncol(bulk_df)])
      quant_bulk_75 <- quantile(apply(bulk_matrix, 2, function(x) {
        quantile(x, 0.75)
      }), 0.5)
      quant_bulk_max <- quantile(apply(bulk_matrix, 2, function(x) {
        quantile(x, 1)
      }), 0.5)
  
      message("Bulk quantiles: 75th=", round(quant_bulk_75, 3), ", max=", round(quant_bulk_max, 3))
  
      ## Generate simulated references from clusters
      for (i in 1:n_clusters) {
        message("Generating reference ", i, "/", n_clusters)
    
        sim_expr <- list()
        for (ct in celltypes) {
          # Get cells in cluster i
          cluster_idx <- which(cl_list[[ct]]$cluster == i)
          cluster_expr <- sc_counts_list[[ct]][cluster_idx, genes, drop = FALSE]
      
          # Average expression in cluster
          mean_expr <- colMeans(cluster_expr)
      
          # Normalize by quantile
          quant <- quantile(mean_expr, quantile_norm)
          mean_expr_norm <- mean_expr / quant
      
          # Calculate scaling factor to match bulk
          quant_sc_75 <- quantile(mean_expr_norm, 0.75)
          s_factor <- quant_bulk_75 / quant_sc_75
      
          # Constrain scaling factor
          s_factor <- min(s_factor, 5)
          s_factor <- max(s_factor, 0.2)
      
          message("  ", ct, ": scale factor = ", round(s_factor, 3))
      
          # Apply scaling
          mean_expr_scaled <- s_factor * mean_expr_norm
      
          # Special scaling for highly expressed genes
          idx_large <- which(mean_expr_norm > 1)
          if (length(idx_large) >= 2) {
            tmp_idx <- mean_expr_scaled[idx_large]
            tmp_min <- min(mean_expr_scaled[idx_large])
            tmp_max <- max(mean_expr_scaled[idx_large])
        
            # Scale to bulk max
            vv <- quant_bulk_max * (tmp_idx - tmp_min) / (tmp_max - tmp_min)
               mean_expr_scaled[idx_large] <- s_factor + vv
             }
      
           sim_expr[[ct]] <- mean_expr_scaled
         }
    
         # Create wide format: gene_celltype columns
        sim_vec <- c()
        for (gene in genes) {
          for (ct in celltypes) {
            sim_vec <- c(sim_vec, sim_expr[[ct]][gene])
          }
        }
        ll.multi[[ii]] <- sim_vec
      }
  }

  
  # Create data.frame
  expr_matrix <- do.call(rbind, ll.multi)
  df_out <- data.frame(Sample = paste0("sim_", 1:n_clusters), expr_matrix, check.names = FALSE)
  
  # Column names: Gene_CellType
  gene_cell <- as.vector(sapply(genes, function(g) {
    paste0(g, "_", celltypes)
  }))
  colnames(df_out) <- c("Sample", gene_cell)
  
  # Write output
  if (!is.null(output_file)) {
    readr::write_csv(df_out, output_file)
    message("Reference data saved to: ", output_file)
  }
  
  return(df_out)
}

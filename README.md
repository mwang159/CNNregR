# CNNregR

R interface for CNNreg cell type deconvolution of bulk RNA-seq using scRNA-seq reference data.

## Overview

CNNregR provides R functions to:

1. **Preprocess scRNA-seq and bulk RNA-seq data** with edgeR TMM normalization
2. **Select discriminative marker genes** using variance and CV-based scoring
3. **Generate reference CSE profiles** via k-means clustering
4. **Run CNNreg deconvolution** directly from R

## Key Features

- **edgeR TMM normalization** for robust library size adjustment
- **K-means clustering** for generating diverse reference profiles
- **CV and variance filtering** for optimal marker gene selection
- **BioMart integration** for organism-specific gene biotype filtering

## Installation

### Prerequisites

CNNregR requires the [CNNreg](https://github.com/mwang159/CNNreg) Python package to be installed:

```bash
pip install CNNreg
```

See the [CNNreg](https://github.com/mwang159/CNNreg) for more instructions. 

### Install CNNregR

```r
# Install from GitHub (recommended)
devtools::install_github("mwang159/CNNregR")
```

## Quick Start

### Complete Workflow (Recommended)

```r
library(CNNregR)

# Check CNNreg installation
check_cnnreg_install()

# Wrapper function for the whole process
result <- deconvolve(
  bulk_counts = my_bulk_data,
  sc_counts_list = list(
    Astrocyte = astro_counts,
    Oligodendrocyte = oligo_counts,
    Neuron = neuron_counts,
    Microglia = micro_counts
  ),
  output_dir = "results/deconvolution/",
  prefix = "deconv",
  epochs = 10000,
  n_clusters = 5,
  filter_biotype = TRUE,
  organism = "hsapiens"
)

# View results
print(result$genes)           # Selected marker genes
print(result$size_factors)    # TMM normalization factors
head(result$proportions)      # Predicted cell proportions
```

### Step-by-Step Workflow

```r
library(CNNregR)

# Step 1: Get gene annotations (optional, auto-queried if needed)
anno <- get_gene_annotation(
  genes = colnames(bulk_data),
  organism = "hsapiens"
)

# Step 2: Select genes with TMM normalization
gene_result <- select_genes(
  bulk_counts = bulk_data,
  sc_counts_list = sc_list,
  filter_biotype = TRUE,
  gene_biotypes = c("protein_coding", "lincRNA", "processed_pseudogene"),
  organism = "hsapiens",
  bulk_min_median_quantile = 0.1,
  bulk_max_median_quantile = 0.999,
  bulk_min_cv_quantile = 0.25,
  ref_var_quantile = 0.5
)

# Step 3: Preprocess bulk RNA-seq with TMM
bulk_df <- preprocess_bulk(
  bulk_counts = bulk_data,
  genes = gene_result$genes,
  size_factors = gene_result$size_factors,
  output_file = "bulk_processed.csv"
)

# Step 4: Generate reference CSE via k-means
ref_df <- preprocess_reference(
  sc_counts_list = sc_list,
  bulk_df = bulk_df,
  genes = bulk_df$Gene,
  n_clusters = 5,
  output_file = "ref_processed.csv"
)

# Step 5: Run CNNreg deconvolution
result <- run_cnnreg(
  bulk_file = "bulk_processed.csv",
  ref_file = "ref_processed.csv",
  output_dir = "output/",
  n_celltypes = length(sc_list),
  epochs = 10000,
  prefix = "my_deconv"
)

# Read proportions
props <- read_cnnreg_proportions(result$proportions_file)
```

## Input Data Format

### Bulk RNA-seq
- **Matrix or data.frame**: rows = samples, columns = genes
- **Raw counts** (not normalized)
- Example:
```r
         Gene1  Gene2  Gene3
Sample1    450    120    890
Sample2    520    135    920
Sample3    380    110    750
```

### scRNA-seq Reference
- **Named list** of matrices, one per cell type
- Each matrix: rows = cells, columns = genes
- **Raw counts** (not normalized)
- Example:
```r
sc_counts_list <- list(
  Astrocyte = astro_matrix,        # 500 cells × 20000 genes
  Oligodendrocyte = oligo_matrix,  # 800 cells × 20000 genes
  Neuron = neuron_matrix           # 1200 cells × 20000 genes
)
```

### From Seurat Object

```r
# Extract counts from Seurat object
sc_counts_list <- format_seurat_object(
  seurat_obj = my_seurat,
  celltype_col = "celltype",
  assay = "RNA",
  slot = "counts"
)
```

## Key Functions

| Function | Description |
|----------|-------------|
| `deconvolve()` | **Complete workflow**: TMM normalize + k-means + deconvolve |
| `get_gene_annotation()` | Query BioMart for gene biotype/chromosome annotations |
| `select_genes()` | Select genes using edgeR TMM + variance + CV |
| `preprocess_bulk()` | Normalize bulk with TMM size factors |
| `preprocess_reference()` | Generate reference via k-means clustering |
| `run_cnnreg()` | Call CNNreg CLI from R |
| `check_cnnreg_install()` | Verify CNNreg Python package |
| `validate_input_data()` | Check data format and quality |

## Output Files

| File | Description |
|------|-------------|
| `bulk_*.csv` | Preprocessed bulk data (Gene × Samples) |
| `ref_*.csv` | Reference CSE (Sample × Gene_CellType) |
| `Prop_predicted_*.csv` | **Final cell proportions** |
| `Prop_predicted_*_epoch_N.csv` | Checkpoint proportions |
| `cellprop_model.pt` | Trained PyTorch model |

## Citation

If you use CNNregR in your research, please cite:

## License

MIT License

## Authors

- Xue Wang (wang.xue@mayo.edu)
- Yuanhang Liu (liu.yuanhang@mayo.edu)

## Links

- Python Package: [CNNreg on PyPI](https://pypi.org/project/CNNreg/)
- CNNreg: [CNNreg on GitHub](https://github.com/mwang159/CNNreg)
- CNNregR: [mwang159/CNNregR](https://github.com/mwang159/CNNregR)
- Issues: [Report bugs](https://github.com/mwang159/CNNregR/issues)
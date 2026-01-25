# Data Directory

This directory contains the data files used in the ferroDG project.

## Directory Structure

```
data/
├── GSE233363/              # Original single-cell RNA-seq data
│   └── Seurat_combined_with_Celltype_Article.rds
└── processed/              # Processed data files
    ├── combined_processed.rds
    └── neuro_with_ferroptosis_scores.rds
```

## Data Sources

### GSE233363
- **Description**: Single-cell RNA sequencing data of mouse dentate gyrus
- **Organism**: Mus musculus (mouse)
- **Tissue**: Hippocampus - Dentate Gyrus
- **Age Groups**: Young, Middle, Old
- **Cell Types**: Multiple cell types including neurogenic lineage cells (qNSC, nIPC, Neuroblast, GC)

### Data Access
The data can be downloaded from the Gene Expression Omnibus (GEO):
- **GEO Accession**: GSE233363
- **URL**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363

## Data Preparation

### Required Files
To run the analysis pipeline, you need to download the following file:

1. **Seurat_combined_with_Celltype_Article.rds**
   - Place in: `data/GSE233363/`
   - This is the main Seurat object containing all cells and metadata

### Download Instructions

#### Option 1: Using GEOquery in R
```r
# Install GEOquery if not already installed
if (!require("GEOquery")) {
  install.packages("BiocManager")
  BiocManager::install("GEOquery")
}

# Download the data
library(GEOquery)
# Note: Check the GEO page for the correct file to download
```

#### Option 2: Manual Download
1. Visit the GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363
2. Download the Seurat object file
3. Place it in the `data/GSE233363/` directory

### Data Processing

The analysis pipeline will automatically process the data:

1. **Stage 1** ([`R/01_data_preprocessing.R`](../R/01_data_preprocessing.R)):
   - Load the Seurat object
   - Collapse age groups (Young/Old)
   - Subset neurogenic lineage cells
   - Save processed objects

2. **Stage 2** ([`R/02_ferroptosis_analysis.R`](../R/02_ferroptosis_analysis.R)):
   - Calculate ferroptosis module scores
   - Extract regulator genes
   - Save updated objects

## Processed Data Files

### combined_processed.rds
- **Description**: Full Seurat object with age grouping
- **Location**: `data/processed/combined_processed.rds`
- **Content**: All cells with Age_collapsed metadata

### neuro_with_ferroptosis_scores.rds
- **Description**: Neurogenic lineage subset with ferroptosis scores
- **Location**: `data/processed/neuro_with_ferroptosis_scores.rds`
- **Content**: 
  - Neurogenic lineage cells only (qNSC, nIPC, Neuroblast, GC)
  - Ferroptosis promoter score
  - Ferroptosis inhibitor score
  - All original metadata

## Data Format

### Seurat Object Structure
```r
# Access metadata
neuro@meta.data

# Access expression data
neuro[["RNA"]]@counts
neuro[["RNA"]]@data

# Access dimensionality reduction
neuro[["umap"]]@cell.embeddings
```

### Key Metadata Columns
- `Celltype_Article`: Cell type annotation
- `timepoint`: Original age timepoint
- `Age_collapsed`: Collapsed age group (Young/Old)
- `FerroPromoter1`: Ferroptosis promoter module score
- `FerroInhibitor1`: Ferroptosis inhibitor module score

## Data Size Considerations

- **Original Seurat object**: ~500 MB - 1 GB (depending on compression)
- **Processed neurogenic subset**: ~100-200 MB
- **Memory requirements**: Minimum 8 GB RAM recommended

## Data Privacy and Usage

- This data is publicly available from GEO
- Please cite the original study when using this data
- Follow GEO's data usage policies

## Troubleshooting

### File Not Found Error
If you encounter "Seurat object file not found" error:
1. Check that the file is in the correct directory: `data/GSE233363/`
2. Verify the file name matches exactly: `Seurat_combined_with_Celltype_Article.rds`
3. Ensure the file is not corrupted

### Memory Issues
If you encounter memory issues:
1. Close other applications
2. Increase R memory limit (if applicable)
3. Consider using a machine with more RAM

### Data Format Issues
If the Seurat object format is incompatible:
1. Check your Seurat version (requires >= 4.0.0)
2. Try re-downloading the data
3. Contact the original data authors

## Additional Resources

- [Seurat Documentation](https://satijalab.org/seurat/)
- [GEO Data Submission Guidelines](https://www.ncbi.nlm.nih.gov/geo/info/submission.html)
- [Single-cell RNA-seq Best Practices](https://www.10xgenomics.com/resources/analysis-guides)

## Contact

For data-related questions or issues, please:
1. Check the original GEO page for data documentation
2. Submit an issue on GitHub
3. Contact the project maintainers

# Methodology

This document describes the methodology used in the ferroDG project for analyzing ferroptosis activity in the dentate gyrus neurogenic lineage.

## Overview

The ferroDG project analyzes single-cell RNA sequencing data to investigate ferroptosis activity across the neurogenic lineage (qNSC → nIPC → Neuroblast → GC) in young and old mice.

## Data Source

### Dataset
- **GEO Accession**: GSE233363
- **Organism**: Mus musculus (mouse)
- **Tissue**: Hippocampus - Dentate Gyrus
- **Technology**: Single-cell RNA sequencing

### Age Groups
The original dataset contains multiple timepoints. For this analysis, we collapse them into two main groups:
- **Young**: Young, Young1, Young2
- **Old**: Old, Old1, Old2

Middle-aged samples are excluded to focus on the comparison between young and old animals.

## Analysis Pipeline

### Stage 1: Data Preprocessing

#### 1.1 Data Loading
- Load the Seurat object containing all cells
- Set the default assay to "RNA"

#### 1.2 Age Grouping
```r
# Collapse timepoints into Young/Old
combined$Age_collapsed <- NA_character_
combined$Age_collapsed[combined$Age %in% c("Young", "Young1", "Young2")] <- "Young"
combined$Age_collapsed[combined$Age %in% c("Old", "Old1", "Old2")] <- "Old"
```

#### 1.3 Neurogenic Lineage Subset
Extract cells belonging to the neurogenic lineage:
- **qNSC**: Quiescent neural stem cells
- **nIPC**: Neural intermediate progenitor cells
- **Neuroblast**: Immature neurons
- **GC**: Granule cells

#### 1.4 Data Filtering
- Remove cells with NA age values
- Keep only Young and Old age groups
- Ensure proper factor ordering

### Stage 2: Ferroptosis Analysis

#### 2.1 Ferroptosis Gene Sets

We use three categories of ferroptosis-related genes:

##### Promoter Genes (25 genes)
Genes that promote ferroptosis:
- **Iron metabolism**: TFRC, NCOA4, HMOX1, ACO1, IREB2
- **Lipid peroxidation**: ACSL4, LPCAT3, ALOX5, ALOX12, ALOX15
- **Other regulators**: PTGS2, CHAC1, NOX1, NOX4, etc.

##### Inhibitor Genes (40 genes)
Genes that inhibit ferroptosis:
- **Antioxidant system**: GPX4, SLC7A11, GCLC, GCLM, GSS
- **Iron storage**: FTH1, NFS1
- **Nrf2 pathway**: NFE2L2, KEAP1
- **Other protective factors**: AKR1C1, AKR1C2, AKR1C3, etc.

##### Regulator Genes (4 genes)
Key ferroptosis regulators:
- NQO1, VDAC2, TP53, RAS

#### 2.2 Gene Symbol Conversion
Human gene symbols are converted to mouse format:
```r
human_to_mouse_symbol <- function(genes) {
  genes <- tolower(genes)
  paste0(toupper(substr(genes, 1, 1)), substr(genes, 2, nchar(genes)))
}
```

#### 2.3 Gene Set Intersection
Only genes present in the dataset are retained for analysis.

#### 2.4 Module Score Calculation
Using Seurat's `AddModuleScore` function:

```r
# Promoter score
neuro <- AddModuleScore(
  object   = neuro,
  features = list(ferroptosis_genes_use$Ferro_Promoter),
  name     = "FerroPromoter",
  assay    = "RNA"
)

# Inhibitor score
neuro <- AddModuleScore(
  object   = neuro,
  features = list(ferroptosis_genes_use$Ferro_Inhibitor),
  name     = "FerroInhibitor",
  assay    = "RNA"
)
```

The module score represents the average expression of genes in the set, normalized for control gene expression.

### Stage 3: Visualization

#### 3.1 UMAP Visualization
- Display cell distribution in UMAP space
- Color by cell type
- Split by age group (Young vs Old)
- Include cell counts in labels

#### 3.2 Violin Plots
Show the distribution of ferroptosis scores:
- X-axis: Cell type
- Split by age group
- Y-axis: Module score or gene expression

#### 3.3 Line Plots
Display mean values across age groups:
- X-axis: Age group (Young → Old)
- Y-axis: Mean score or expression
- Lines connect cell types

#### 3.4 Trajectory Plots
Show changes along the neurogenic lineage:
- X-axis: Cell type (qNSC → nIPC → Neuroblast → GC)
- Y-axis: Median score or expression
- Separate lines for Young and Old

## Statistical Considerations

### Module Score Calculation
- Uses Seurat's `AddModuleScore` algorithm
- Normalizes for control gene expression
- Random seed set to 1 for reproducibility

### Data Filtering
- Cells with NA age values are excluded
- Only Young and Old age groups are analyzed
- Middle-aged samples are not included

### Multiple Testing
For differential expression analysis (if performed):
- Consider adjusting p-values for multiple testing
- Use appropriate thresholds (e.g., padj < 0.05)

## Reproducibility

### Random Seed
```r
set.seed(1)
```

### Version Control
- R scripts are version-controlled
- Gene lists are stored as plain text files
- All parameters are documented

### Output Format
All figures are saved in both:
- **PNG**: Raster format, 300 DPI
- **PDF**: Vector format, using cairo_pdf for font compatibility

## Limitations

### Data Limitations
- Cross-sectional data (not longitudinal)
- Limited sample size per age group
- Batch effects may exist

### Methodological Limitations
- Module scores are relative, not absolute
- Gene set overlap may affect interpretation
- Human-to-mouse gene conversion is approximate

### Biological Limitations
- Ferroptosis activity is inferred from gene expression
- Post-transcriptional regulation is not captured
- Protein-level activity may differ from mRNA

## Future Directions

### Potential Extensions
1. **Pseudotime Analysis**: Infer differentiation trajectories
2. **Spatial Analysis**: Incorporate spatial transcriptomics
3. **Integration with Other Modalities**: ATAC-seq, proteomics
4. **Machine Learning**: Predict ferroptosis susceptibility

### Additional Analyses
1. **Pathway Enrichment**: Identify enriched pathways
2. **Gene Regulatory Networks**: Infer regulatory relationships
3. **Cell-Cell Communication**: Analyze ligand-receptor interactions
4. **Metabolic Analysis**: Investigate metabolic changes

## References

### Key Papers on Ferroptosis
1. Stockwell BR, et al. (2017). Ferroptosis: A Regulated Cell Death Nexus Linking Metabolism, Redox Biology, and Disease. Cell.
2. Dixon SJ, et al. (2012). Ferroptosis: An Iron-Dependent Form of Nonapoptotic Cell Death. Cell.

### Seurat Documentation
- Seurat: https://satijalab.org/seurat/
- AddModuleScore: https://satijalab.org/seurat/reference/addmodulescore

### Original Data
- GSE233363: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363

## Contact

For questions about the methodology, please:
1. Check the main README.md
2. Review the R scripts with inline comments
3. Submit an issue on GitHub

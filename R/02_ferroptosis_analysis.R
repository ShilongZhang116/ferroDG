## =========================================================
## Script: Ferroptosis Analysis in Neurogenic Lineage
## Purpose: Calculate ferroptosis module scores and analyze gene expression
## Input: Processed Seurat object from 01_data_preprocessing.R
## Output: Seurat object with ferroptosis scores
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## =========================================================
## 0. Set working directory and paths
## =========================================================
project_dir <- file.path(dirname(getwd()))
data_dir <- file.path(project_dir, "data", "processed")
output_dir <- file.path(project_dir, "data", "processed")

message("Project directory: ", project_dir)
message("Data directory: ", data_dir)

## =========================================================
## 1. Load processed Seurat object
## =========================================================
neuro_file <- file.path(data_dir, "neuro_processed.rds")

if (!file.exists(neuro_file)) {
  stop("Processed Seurat object not found: ", neuro_file,
       "\nPlease run 01_data_preprocessing.R first.")
}

neuro <- readRDS(neuro_file)
DefaultAssay(neuro) <- "RNA"

message("Loaded neurogenic lineage subset with ", ncol(neuro), " cells")

## =========================================================
## 2. Define ferroptosis gene sets
## =========================================================

# Human gene symbols (will be converted to mouse format)
ferroptosis_genes <- list(
  Ferro_Promoter = c(
    "TFRC","NCOA4","HMOX1","PHKG2","ACO1","IREB2",
    "ACSL4","LPCAT3","ALOX5","ALOX12","ALOX15",
    "DPP4","ACSF2","ZEB1","PTGS2","PEBP1",
    "CARS","CHAC1","NOX1","ABCC1","NOX4",
    "GLS2","ATP5G3","VDAC3","SAT1"
  ),
  Ferro_Inhibitor = c(
    "FTH1","NFS1","MT1G","STEAP3","HSBP1","FANCD2","TERF1",
    "GPX4","AKR1C1","AKR1C2","AKR1C3","HMGCR","SQLE","CISD1","ACSL3",
    "CS","CBS","FDFT1","FADS2","ACACA",
    "SLC7A11","GCLC","GCLM","GSS","NFE2L2","KEAP1","SLC38A1",
    "G6PD","PGD","SLC1A5","GOT1",
    "CD44","CRYAB","EMC2","HSPB1","AIFM2","RPL8"
  ),
  Ferro_Regulator = c("NQO1","VDAC2","TP53","RAS")
)

# Convert human to mouse gene symbols
human_to_mouse_symbol <- function(genes) {
  genes <- tolower(genes)
  paste0(toupper(substr(genes, 1, 1)), substr(genes, 2, nchar(genes)))
}

ferroptosis_genes_mouse <- lapply(ferroptosis_genes, human_to_mouse_symbol)

# Intersect with genes present in the dataset
genes_present <- rownames(neuro)
ferroptosis_genes_use <- lapply(ferroptosis_genes_mouse, function(g) intersect(g, genes_present))

message("Ferroptosis genes retained after intersection:")
print(sapply(ferroptosis_genes_use, length))
print(ferroptosis_genes_use)

## =========================================================
## 3. Calculate module scores
## =========================================================
set.seed(1)  # For reproducibility

# Promoter score
if (!any(grepl("^FerroPromoter", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Promoter),
    name     = "FerroPromoter",
    assay    = "RNA"
  )
  message("Calculated Ferroptosis Promoter module score")
}

# Inhibitor score
if (!any(grepl("^FerroInhibitor", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Inhibitor),
    name     = "FerroInhibitor",
    assay    = "RNA"
  )
  message("Calculated Ferroptosis Inhibitor module score")
}

## =========================================================
## 4. Get column names for scores
## =========================================================
promoter_col <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]

message("Promoter score column: ", promoter_col)
message("Inhibitor score column: ", inhibitor_col)

## =========================================================
## 5. Extract regulator genes
## =========================================================
regulator_genes <- unique(ferroptosis_genes_use$Ferro_Regulator)
regulator_genes <- regulator_genes[regulator_genes %in% rownames(neuro)]

message("Regulator genes available: ", paste(regulator_genes, collapse = ", "))

## =========================================================
## 6. Summary statistics for scores
## =========================================================
message("\n=== Ferroptosis Score Summary ===")

# Promoter score summary
promoter_stats <- neuro@meta.data[[promoter_col]]
message("Promoter score:")
message("  Mean: ", round(mean(promoter_stats, na.rm = TRUE), 4))
message("  Median: ", round(median(promoter_stats, na.rm = TRUE), 4))
message("  SD: ", round(sd(promoter_stats, na.rm = TRUE), 4))

# Inhibitor score summary
inhibitor_stats <- neuro@meta.data[[inhibitor_col]]
message("Inhibitor score:")
message("  Mean: ", round(mean(inhibitor_stats, na.rm = TRUE), 4))
message("  Median: ", round(median(inhibitor_stats, na.rm = TRUE), 4))
message("  SD: ", round(sd(inhibitor_stats, na.rm = TRUE), 4))

# Scores by cell type and age
message("\nPromoter score by cell type and age:")
print(neuro@meta.data %>%
        group_by(Celltype_Article, Age_collapsed) %>%
        summarise(
          mean_promoter = mean(.data[[promoter_col]], na.rm = TRUE),
          .groups = "drop"
        ))

message("\nInhibitor score by cell type and age:")
print(neuro@meta.data %>%
        group_by(Celltype_Article, Age_collapsed) %>%
        summarise(
          mean_inhibitor = mean(.data[[inhibitor_col]], na.rm = TRUE),
          .groups = "drop"
        ))

## =========================================================
## 7. Save updated Seurat object
## =========================================================
saveRDS(neuro, file.path(output_dir, "neuro_with_ferroptosis_scores.rds"))
message("\nSaved: ", file.path(output_dir, "neuro_with_ferroptosis_scores.rds"))

## =========================================================
## 8. Save gene lists
## =========================================================
gene_list_dir <- file.path(project_dir, "gene_lists")
dir.create(gene_list_dir, showWarnings = FALSE, recursive = TRUE)

# Save promoter genes
writeLines(ferroptosis_genes_use$Ferro_Promoter,
           file.path(gene_list_dir, "ferroptosis_promoter.txt"))
message("Saved: ", file.path(gene_list_dir, "ferroptosis_promoter.txt"))

# Save inhibitor genes
writeLines(ferroptosis_genes_use$Ferro_Inhibitor,
           file.path(gene_list_dir, "ferroptosis_inhibitor.txt"))
message("Saved: ", file.path(gene_list_dir, "ferroptosis_inhibitor.txt"))

# Save regulator genes
writeLines(regulator_genes,
           file.path(gene_list_dir, "ferroptosis_regulator.txt"))
message("Saved: ", file.path(gene_list_dir, "ferroptosis_regulator.txt"))

message("\nFerroptosis analysis completed successfully!")

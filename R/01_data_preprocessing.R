## =========================================================
## Script: Data Preprocessing for Ferroptosis Analysis
## Purpose: Load and preprocess single-cell RNA-seq data
## Input: Seurat object (GSE233363)
## Output: Processed Seurat object with age grouping
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## =========================================================
## 0. Set working directory and paths
## =========================================================
# Note: Update this path to your project directory
project_dir <- file.path(dirname(getwd()))
data_dir <- file.path(project_dir, "data", "GSE233363")
output_dir <- file.path(project_dir, "data", "processed")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Project directory: ", project_dir)
message("Data directory: ", data_dir)
message("Output directory: ", output_dir)

## =========================================================
## 1. Load Seurat object
## =========================================================
seurat_file <- file.path(data_dir, "Seurat_combined_with_Celltype_Article.rds")

if (!file.exists(seurat_file)) {
  stop("Seurat object file not found: ", seurat_file, 
       "\nPlease download the data from GSE233363 and place it in the data directory.")
}

combined <- readRDS(seurat_file)
DefaultAssay(combined) <- "RNA"

message("Loaded Seurat object with ", ncol(combined), " cells and ", nrow(combined), " genes")

## =========================================================
## 2. Age grouping
## =========================================================
# Collapse timepoints into Young/Old groups
combined$Age <- combined$timepoint

combined$Age_collapsed <- NA_character_
combined$Age_collapsed[combined$Age %in% c("Young", "Young1", "Young2")] <- "Young"
combined$Age_collapsed[combined$Age %in% c("Old", "Old1", "Old2")] <- "Old"

combined$Age_collapsed <- factor(combined$Age_collapsed, levels = c("Young", "Old"))

message("Age distribution:")
print(table(combined$Age, useNA = "ifany"))
message("Collapsed age groups:")
print(table(combined$Age_collapsed, useNA = "ifany"))

## =========================================================
## 3. Define neurogenic lineage cell types
## =========================================================
neuro_ct <- c("qNSC", "nIPC", "Neuroblast", "GC")

message("Neurogenic lineage cell types: ", paste(neuro_ct, collapse = ", "))

## =========================================================
## 4. Subset neurogenic lineage cells
## =========================================================
neuro <- subset(combined, subset = Celltype_Article %in% neuro_ct)

# Keep only Young/Old cells
neuro <- subset(neuro, subset = !is.na(Age_collapsed) & Age_collapsed %in% c("Young", "Old"))
neuro$Age_collapsed <- factor(neuro$Age_collapsed, levels = c("Young", "Old"))
neuro$Celltype_Article <- factor(neuro$Celltype_Article, levels = neuro_ct)

message("Neurogenic lineage subset:")
print(table(neuro$Celltype_Article, neuro$Age_collapsed, useNA = "ifany"))

## =========================================================
## 5. Save processed objects
## =========================================================
# Save full combined object
saveRDS(combined, file.path(output_dir, "combined_processed.rds"))
message("Saved: ", file.path(output_dir, "combined_processed.rds"))

# Save neurogenic lineage subset
saveRDS(neuro, file.path(output_dir, "neuro_processed.rds"))
message("Saved: ", file.path(output_dir, "neuro_processed.rds"))

## =========================================================
## 6. Summary statistics
## =========================================================
message("\n=== Summary Statistics ===")
message("Total cells (combined): ", ncol(combined))
message("Total genes (combined): ", nrow(combined))
message("Neurogenic lineage cells: ", ncol(neuro))
message("Unique cell types: ", length(unique(neuro$Celltype_Article)))
message("Age groups: ", paste(levels(neuro$Age_collapsed), collapse = ", "))

message("\nData preprocessing completed successfully!")

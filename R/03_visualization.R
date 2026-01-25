## =========================================================
## Script: Visualization for Ferroptosis Analysis
## Purpose: Generate publication-quality figures
## Input: Seurat object with ferroptosis scores from 02_ferroptosis_analysis.R
## Output: Figures (PNG + PDF) in figures/ directory
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(grid)
})

## =========================================================
## 0. Set working directory and paths
## =========================================================
project_dir <- file.path(dirname(getwd()))
data_dir <- file.path(project_dir, "data", "processed")
fig_dir <- file.path(project_dir, "figures", "stage-2")

# Create figures directory
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("Project directory: ", project_dir)
message("Figures directory: ", fig_dir)

## =========================================================
## Helper functions
## =========================================================

# Save plot to PNG and PDF
save_plot_local <- function(p, filename, width = 8, height = 5, dpi = 300) {
  out_png <- file.path(fig_dir, paste0(filename, ".png"))
  out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  
  ggsave(filename = out_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(filename = out_pdf, plot = p, width = width, height = height, device = cairo_pdf)
  
  message("[Saved] ", out_png)
  message("[Saved] ", out_pdf)
}

# Shorten title for plots
short_title <- function(x) {
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

## =========================================================
## 1. Load Seurat object with ferroptosis scores
## =========================================================
neuro_file <- file.path(data_dir, "neuro_with_ferroptosis_scores.rds")

if (!file.exists(neuro_file)) {
  stop("Seurat object with ferroptosis scores not found: ", neuro_file,
       "\nPlease run 02_ferroptosis_analysis.R first.")
}

neuro <- readRDS(neuro_file)
DefaultAssay(neuro) <- "RNA"

message("Loaded neurogenic lineage subset with ", ncol(neuro), " cells")

## =========================================================
## 2. Define colors
## =========================================================

# Cell type colors
celltype_colors <- c(
  "qNSC"       = "#FDAE6B",
  "nIPC"       = "#3182BD",
  "Neuroblast" = "#08519C",
  "GC"         = "#9ECAE1"
)

# Age colors
age_colors <- c("Young" = "#4DBBD5", "Old" = "#D73027")

# Neurogenic lineage order
neuro_order <- c("qNSC", "nIPC", "Neuroblast", "GC")

## =========================================================
## 3. Get ferroptosis score column names
## =========================================================
promoter_col <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]

## =========================================================
## 4. Load gene lists
## =========================================================
gene_list_dir <- file.path(project_dir, "gene_lists")

regulator_genes <- readLines(file.path(gene_list_dir, "ferroptosis_regulator.txt"))
regulator_genes <- regulator_genes[regulator_genes %in% rownames(neuro)]

message("Regulator genes: ", paste(regulator_genes, collapse = ", "))

## =========================================================
## 5. Figure 1: UMAP of neurogenic lineage split by age
## =========================================================

# Add cell counts to age labels
age_n <- neuro@meta.data %>%
  dplyr::count(Age_collapsed, name = "n_cells") %>%
  dplyr::mutate(Age_strip = paste0(as.character(Age_collapsed), " (n=", n_cells, ")"))

neuro$Age_collapsed_strip <- as.character(neuro$Age_collapsed)
neuro$Age_collapsed_strip[neuro$Age_collapsed_strip == "Young"] <-
  age_n$Age_strip[age_n$Age_collapsed == "Young"]
neuro$Age_collapsed_strip[neuro$Age_collapsed_strip == "Old"] <-
  age_n$Age_strip[age_n$Age_collapsed == "Old"]

neuro$Age_collapsed_strip <- factor(
  neuro$Age_collapsed_strip,
  levels = age_n$Age_strip[match(c("Young", "Old"),
                                 as.character(age_n$Age_collapsed))]
)

# Add cell counts to cell type labels
ct_age_n <- neuro@meta.data %>%
  dplyr::count(Celltype_Article, Age_collapsed, name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from  = Age_collapsed,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    legend_label = paste0(
      Celltype_Article,
      " (Y=", Young, ", O=", Old, ")"
    )
  )

neuro$Celltype_Article_legend <- as.character(neuro$Celltype_Article)
for (i in seq_len(nrow(ct_age_n))) {
  ct <- as.character(ct_age_n$Celltype_Article[i])
  lb <- as.character(ct_age_n$legend_label[i])
  neuro$Celltype_Article_legend[neuro$Celltype_Article_legend == ct] <- lb
}

neuro$Celltype_Article_legend <- factor(
  neuro$Celltype_Article_legend,
  levels = ct_age_n$legend_label[match(neuro_order,
                                       ct_age_n$Celltype_Article)]
)

celltype_colors_legend <- celltype_colors[neuro_order]
names(celltype_colors_legend) <- ct_age_n$legend_label[
  match(neuro_order, ct_age_n$Celltype_Article)
]

p_fig1 <- DimPlot(
  neuro,
  reduction = "umap",
  group.by  = "Celltype_Article_legend",
  split.by  = "Age_collapsed_strip",
  label     = FALSE
) +
  scale_color_manual(values = celltype_colors_legend, name = "Cell type") +
  ggtitle("Figure 1. Neurogenic lineage UMAP split by age") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

save_plot_local(p_fig1, "Fig1_UMAP_Neuro_Celltype_SplitByAge", width = 10, height = 5)

## =========================================================
## 6. Figure 2: Violin plots of ferroptosis scores
## =========================================================

# 2A: Promoter score
p_fig2A <- VlnPlot(
  neuro,
  features = promoter_col,
  group.by = "Celltype_Article",
  split.by = "Age_collapsed",
  pt.size  = 0
) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Figure 2A. Promoter module score by cell type and age") +
  theme(plot.title = element_text(hjust = 0.5))

save_plot_local(p_fig2A, "Fig2A_Violin_PromoterScore_Celltype_byAge", width = 8, height = 5)

# 2B: Inhibitor score
p_fig2B <- VlnPlot(
  neuro,
  features = inhibitor_col,
  group.by = "Celltype_Article",
  split.by = "Age_collapsed",
  pt.size  = 0
) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Figure 2B. Inhibitor module score by cell type and age") +
  theme(plot.title = element_text(hjust = 0.5))

save_plot_local(p_fig2B, "Fig2B_Violin_InhibitorScore_Celltype_byAge", width = 8, height = 5)

# 2C: Regulator genes
if (length(regulator_genes) > 0) {
  p_fig2C <- VlnPlot(
    neuro,
    features = regulator_genes,
    group.by = "Celltype_Article",
    split.by = "Age_collapsed",
    pt.size  = 0,
    ncol     = min(2, length(regulator_genes))
  ) +
    scale_fill_manual(values = age_colors, drop = TRUE) +
    theme(plot.title = element_text(hjust = 0.5))
  
  save_plot_local(p_fig2C, "Fig2C_Violin_RegulatorGenes_Celltype_byAge", width = 10, height = 5)
}

## =========================================================
## 7. Figure 3: Line plots of ferroptosis scores
## =========================================================

plot_line_summary <- function(df_long, value_col,
                              title_text, ylab_text,
                              out_name,
                              width = 6, height = 5) {
  df2 <- df_long %>%
    select(all_of(c(value_col, "Age_collapsed", "Celltype_Article"))) %>%
    rename(Value = all_of(value_col)) %>%
    mutate(Age_collapsed = factor(Age_collapsed, levels = c("Young", "Old"))) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      n_cells    = n(),
      .groups    = "drop"
    )
  
  p <- ggplot(
    df2,
    aes(
      x     = Age_collapsed,
      y     = mean_value,
      group = Celltype_Article,
      color = Celltype_Article
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = celltype_colors, name = "Cell type") +
    theme_bw() +
    labs(x = "Age group", y = ylab_text, title = short_title(title_text)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  save_plot_local(p, out_name, width = width, height = height)
  invisible(df2)
}

# 3A: Promoter score
df_promoter <- FetchData(neuro, vars = c(promoter_col, "Age_collapsed", "Celltype_Article"))
colnames(df_promoter)[1] <- promoter_col
plot_line_summary(
  df_long    = df_promoter,
  value_col  = promoter_col,
  title_text = "Figure 3A. Promoter module score change from Young to Old",
  ylab_text  = "Mean promoter module score",
  out_name   = "Fig3A_Line_PromoterScore_Celltype_Young_vs_Old",
  width      = 6, height = 5
)

# 3B: Inhibitor score
df_inhib <- FetchData(neuro, vars = c(inhibitor_col, "Age_collapsed", "Celltype_Article"))
colnames(df_inhib)[1] <- inhibitor_col
plot_line_summary(
  df_long    = df_inhib,
  value_col  = inhibitor_col,
  title_text = "Figure 3B. Inhibitor module score change from Young to Old",
  ylab_text  = "Mean inhibitor module score",
  out_name   = "Fig3B_Line_InhibitorScore_Celltype_Young_vs_Old",
  width      = 6, height = 5
)

# 3C: Regulator genes
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    df_g <- FetchData(neuro, vars = c(g, "Age_collapsed", "Celltype_Article"))
    colnames(df_g)[1] <- g
    plot_line_summary(
      df_long    = df_g,
      value_col  = g,
      title_text = paste0("Figure 3C. Regulator gene ", g, " expression change from Young to Old"),
      ylab_text  = paste0("Mean expression of ", g),
      out_name   = paste0("Fig3C_Line_RegulatorGene_", g, "_Celltype_Young_vs_Old"),
      width      = 6, height = 5
    )
  }
}

## =========================================================
## 8. Figure 4: Neurogenic trajectory plots
## =========================================================

plot_traj_median <- function(df_long, value_col,
                             title_text, ylab_text,
                             out_name,
                             width = 6.5, height = 4.5) {
  dfA <- df_long %>%
    filter(Celltype_Article %in% neuro_order) %>%
    mutate(
      Celltype_Article = factor(Celltype_Article, levels = neuro_order),
      Age_collapsed    = factor(Age_collapsed, levels = c("Young", "Old"))
    ) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(median_value = median(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(dfA, aes(x = Celltype_Article, y = median_value, color = Age_collapsed, group = Age_collapsed)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = age_colors, name = "Age") +
    theme_bw() +
    labs(title = short_title(title_text), x = "Neurogenic lineage", y = ylab_text) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold"),
          panel.grid = element_blank())
  
  save_plot_local(p, out_name, width = width, height = height)
  invisible(dfA)
}

# 4A: Promoter score trajectory
df_promoter_meta <- neuro@meta.data %>%
  dplyr::select(Celltype_Article, Age_collapsed, all_of(promoter_col))

plot_traj_median(
  df_long    = df_promoter_meta,
  value_col  = promoter_col,
  title_text = "Figure 4A. Neurogenic trajectory vs ferroptosis promoter activity",
  ylab_text  = "Median promoter module score",
  out_name   = "Fig4A_Trajectory_PromoterScore",
  width      = 6.5, height = 4.5
)

# 4B: Inhibitor score trajectory
df_inhibitor_meta <- neuro@meta.data %>%
  dplyr::select(Celltype_Article, Age_collapsed, all_of(inhibitor_col))

plot_traj_median(
  df_long    = df_inhibitor_meta,
  value_col  = inhibitor_col,
  title_text = "Figure 4B. Neurogenic trajectory vs ferroptosis inhibitor activity",
  ylab_text  = "Median inhibitor module score",
  out_name   = "Fig4B_Trajectory_InhibitorScore",
  width      = 6.5, height = 4.5
)

# 4C: Regulator genes trajectory
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    if (!g %in% rownames(neuro)) next
    
    df_g_meta <- neuro@meta.data %>%
      dplyr::select(Celltype_Article, Age_collapsed)
    
    expr_g <- FetchData(neuro, vars = g)[, 1]
    df_g_meta[[g]] <- expr_g
    
    plot_traj_median(
      df_long    = df_g_meta,
      value_col  = g,
      title_text = paste0("Figure 4C. Neurogenic trajectory vs regulator gene ", g),
      ylab_text  = paste0("Median expression of ", g),
      out_name   = paste0("Fig4C_Trajectory_Regulator_", g),
      width      = 6.5, height = 4.5
    )
  }
}

message("\nVisualization completed successfully!")
message("All figures saved to: ", fig_dir)

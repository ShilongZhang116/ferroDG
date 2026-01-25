## =========================================================
## Script: Ferroptosis × Neurogenesis (qNSC / nIPC / Neuroblast / GC)
## Focus: Young vs Old
## Input:  Seurat_combined_with_Celltype_Article.rds
## Output: Figures saved to ./figs (PNG + PDF)
## Notes:
##   - Age is collapsed into {Young, Old}; other timepoints are removed.
##   - Ferroptosis gene sets are defined as Promoter / Inhibitor / Regulator.
##   - Regulator set is analyzed at single-gene level (module size too small).
##   - All ggplot figures are saved via save_plot_local().
##   - Heatmaps (pheatmap) are saved via dedicated helpers (PDF + PNG).
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
## 0. Paths & I/O
## =========================================================
setwd("H:/Projects/mmFerroptosis/stage-2")

data_dir <- "data/GSE233363"
fig_dir  <- file.path("figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- helper: save ggplot to PNG + PDF ----
save_plot_local <- function(p, filename, width = 8, height = 5, dpi = 300) {
  out_png <- file.path(fig_dir, paste0(filename, ".png"))
  out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  
  ggsave(filename = out_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(filename = out_pdf, plot = p, width = width, height = height, device = cairo_pdf)
  
  message("[Saved] ", out_png)
  message("[Saved] ", out_pdf)
}

## ---- helper: safe title (short) ----
short_title <- function(x) {
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

## =========================================================
## 1. Load Seurat object & collapse age groups
## =========================================================
combined <- readRDS(file.path(data_dir, "Seurat_combined_with_Celltype_Article.rds"))
DefaultAssay(combined) <- "RNA"

combined$Age <- combined$timepoint

## Collapse to Young/Old only; everything else is NA then removed
combined$Age_collapsed <- NA_character_
combined$Age_collapsed[combined$Age %in% c("Young", "Young1", "Young2")] <- "Young"
combined$Age_collapsed[combined$Age %in% c("Old", "Old1", "Old2")]       <- "Old"
combined$Age_collapsed <- factor(combined$Age_collapsed, levels = c("Young", "Old"))

message("== Age table (combined) ==")
print(table(combined$Age, useNA = "ifany"))
message("== Age_collapsed table (combined) ==")
print(table(combined$Age_collapsed, useNA = "ifany"))

## =========================================================
## 2. Subset neurogenic lineage cells (qNSC → nIPC → NB → GC)
## =========================================================
neuro_ct <- c("qNSC", "nIPC", "Neuroblast", "GC")

neuro <- subset(combined, subset = Celltype_Article %in% neuro_ct)
neuro <- subset(neuro, subset = !is.na(Age_collapsed) & Age_collapsed %in% c("Young", "Old"))

neuro$Age_collapsed    <- factor(neuro$Age_collapsed, levels = c("Young", "Old"))
neuro$Celltype_Article <- factor(neuro$Celltype_Article, levels = neuro_ct)

message("== Celltype × Age (neuro) ==")
print(table(neuro$Celltype_Article, neuro$Age_collapsed, useNA = "ifany"))

## =========================================================
## 3. Color system
## =========================================================
celltype_colors_full <- c(
  "qNSC"            = "#FDAE6B",
  "nIPC"            = "#3182BD",
  "Neuroblast"      = "#08519C",
  "GC"              = "#9ECAE1",
  "Pyr"             = "#006D2C",
  "Cajal-Retzius"   = "#74C476",
  "OPC"             = "#E6550D",
  "Oligodendrocyte" = "#FDAE6B",
  "Microglia"       = "#C51B8A",
  "PVM"             = "#F768A1",
  "Astrocyte"       = "#D73027",
  "Endothelial"     = "#542788",
  "Pericyte"        = "#8073AC",
  "SMC"             = "#636363",
  "VLMC"            = "#969696"
)

age_colors_B <- c("Young" = "#4DBBD5", "Old" = "#D73027")
age_colors   <- age_colors_B

## =========================================================
## Figure 0: Global UMAP of all cell types (no age split)
## - Color by Celltype_Article
## - Legend shows Celltype with total cell counts
## Output:
##   - Fig0_UMAP_AllCelltypes_WithCounts.(png/pdf)
## =========================================================

## ---- 0) make sure we are using the full object (not neuro subset) ----
all_cells <- combined
DefaultAssay(all_cells) <- "RNA"

## ---- 1) cell counts per cell type (global, no age split) ----
ct_n_all <- all_cells@meta.data %>%
  dplyr::count(Celltype_Article, name = "n_cells") %>%
  dplyr::mutate(
    legend_label = paste0(Celltype_Article, " (n=", n_cells, ")")
  )

## ---- 2) create legend-mapped celltype variable ----
all_cells$Celltype_Article_legend <- as.character(all_cells$Celltype_Article)

for (i in seq_len(nrow(ct_n_all))) {
  ct <- as.character(ct_n_all$Celltype_Article[i])
  lb <- as.character(ct_n_all$legend_label[i])
  all_cells$Celltype_Article_legend[all_cells$Celltype_Article_legend == ct] <- lb
}

## keep a stable order in legend (use original Celltype order if exists)
celltype_levels <- intersect(names(celltype_colors_full),
                             ct_n_all$Celltype_Article)

legend_levels <- ct_n_all$legend_label[
  match(celltype_levels, ct_n_all$Celltype_Article)
]

all_cells$Celltype_Article_legend <- factor(
  all_cells$Celltype_Article_legend,
  levels = legend_levels
)

## ---- 3) remap colors to legend labels ----
celltype_colors_legend <- celltype_colors_full[celltype_levels]
names(celltype_colors_legend) <- legend_levels

## ---- 4) plot UMAP (no age split) ----
p_fig0_umap <- DimPlot(
  all_cells,
  reduction = "umap",
  group.by  = "Celltype_Article_legend",
  label     = FALSE
) +
  scale_color_manual(values = celltype_colors_legend, name = "Cell type") +
  ggtitle("Figure 0. UMAP of all cell types") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 9)
  )

save_plot_local(
  p_fig0_umap,
  "Fig0_UMAP_AllCelltypes_WithCounts",
  width  = 7.5,
  height = 6
)


## =========================================================
## Figure 1: UMAP of neurogenic lineage split by age (Young vs Old)
## Add: cell counts (n=) per Celltype within each Age facet
## Output:
##   - Fig1_UMAP_Neuro_Celltype_SplitByAge.(png/pdf)
## =========================================================

## ---- 1) facet strip: total cell numbers per age ----
age_n <- neuro@meta.data %>%
  dplyr::count(Age_collapsed, name = "n_cells") %>%
  dplyr::mutate(
    Age_strip = paste0(as.character(Age_collapsed), " (n=", n_cells, ")")
  )

neuro$Age_collapsed_strip <- as.character(neuro$Age_collapsed)
neuro$Age_collapsed_strip[neuro$Age_collapsed_strip == "Young"] <-
  age_n$Age_strip[age_n$Age_collapsed == "Young"]
neuro$Age_collapsed_strip[neuro$Age_collapsed_strip == "Old"]   <-
  age_n$Age_strip[age_n$Age_collapsed == "Old"]

neuro$Age_collapsed_strip <- factor(
  neuro$Age_collapsed_strip,
  levels = age_n$Age_strip[match(c("Young", "Old"),
                                 as.character(age_n$Age_collapsed))]
)

## ---- 2) legend: per Celltype, split Young / Old counts ----
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

## 映射到新的 legend 变量
neuro$Celltype_Article_legend <- as.character(neuro$Celltype_Article)
for (i in seq_len(nrow(ct_age_n))) {
  ct <- as.character(ct_age_n$Celltype_Article[i])
  lb <- as.character(ct_age_n$legend_label[i])
  neuro$Celltype_Article_legend[neuro$Celltype_Article_legend == ct] <- lb
}

neuro$Celltype_Article_legend <- factor(
  neuro$Celltype_Article_legend,
  levels = ct_age_n$legend_label[match(neuro_ct,
                                       ct_age_n$Celltype_Article)]
)

## ---- 3) colors (保持原 Celltype 配色，仅重命名) ----
celltype_colors_legend <- celltype_colors_full[neuro_ct]
names(celltype_colors_legend) <- ct_age_n$legend_label[
  match(neuro_ct, ct_age_n$Celltype_Article)
]

## ---- 4) UMAP ----
p_fig1_umap <- DimPlot(
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

save_plot_local(
  p_fig1_umap,
  "Fig1_UMAP_Neuro_Celltype_SplitByAge",
  width = 10,
  height = 5
)

## =========================================================
## 4. Ferroptosis gene sets (human symbols → mouse-like symbols)
## =========================================================
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
  ## Regulator: analyzed at single-gene level (often too few genes remain)
  Ferro_Regulator = c("NQO1","VDAC2","TP53","RAS")
)

human_to_mouse_symbol <- function(genes) {
  genes <- tolower(genes)
  paste0(toupper(substr(genes, 1, 1)), substr(genes, 2, nchar(genes)))
}

ferroptosis_genes_mouse <- lapply(ferroptosis_genes, human_to_mouse_symbol)

genes_present <- rownames(neuro)
ferroptosis_genes_use <- lapply(ferroptosis_genes_mouse, function(g) intersect(g, genes_present))

message("== Ferroptosis genes retained after intersection (neuro) ==")
print(sapply(ferroptosis_genes_use, length))
print(ferroptosis_genes_use)

regulator_genes <- unique(ferroptosis_genes_use$Ferro_Regulator)
regulator_genes <- regulator_genes[regulator_genes %in% rownames(neuro)]

## =========================================================
## 5. Compute module scores (Promoter / Inhibitor)
## =========================================================
set.seed(1)

if (!any(grepl("^FerroPromoter", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Promoter),
    name     = "FerroPromoter",
    assay    = "RNA"
  )
}
promoter_col <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]

if (!any(grepl("^FerroInhibitor", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Inhibitor),
    name     = "FerroInhibitor",
    assay    = "RNA"
  )
}
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]

## =========================================================
## Figure 2: Violin plots
## Panels:
##   A) Promoter module score (by Celltype × Age)
##   B) Inhibitor module score (by Celltype × Age)
##   C) Regulator genes (single-gene expression; facetted)
## Outputs:
##   - Fig2A_Violin_PromoterScore_Celltype_byAge.(png/pdf)
##   - Fig2B_Violin_InhibitorScore_Celltype_byAge.(png/pdf)
##   - Fig2C_Violin_RegulatorGenes_Celltype_byAge.(png/pdf)   (if any)
## =========================================================
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

regulator_genes <- unique(ferroptosis_genes_use$Ferro_Regulator)
regulator_genes <- regulator_genes[regulator_genes %in% rownames(neuro)]

if (length(regulator_genes) > 0) {
  p_vln_regulator <- VlnPlot(
    neuro,
    features = regulator_genes,
    group.by = "Celltype_Article",
    split.by = "Age_collapsed",
    pt.size  = 0,
    ncol     = min(2, length(regulator_genes))
  ) +
    scale_fill_manual(values = age_colors, drop = TRUE) +
    theme(plot.title = element_text(hjust = 0.5))
  
  save_plot_local(p_vln_regulator, "Fig2C_Violin_RegulatorGenes_Celltype_byAge", width = 10, height = 5)
} else {
  message("[WARN] No regulator genes retained after intersection in neuro. Skip regulator violin.")
}

## =========================================================
## Figure 3: Summary line plots (mean) across age for each cell type
## Panels:
##   A) Promoter module score
##   B) Inhibitor module score
##   C) Regulator genes (single-gene; one plot per gene)
## Outputs:
##   - Fig3A_Line_PromoterScore_Celltype_Young_vs_Old.(png/pdf)
##   - Fig3B_Line_InhibitorScore_Celltype_Young_vs_Old.(png/pdf)
##   - Fig3C_Line_RegulatorGene_<Gene>_Celltype_Young_vs_Old.(png/pdf) (each gene)
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
    scale_color_manual(values = celltype_colors_full[neuro_ct], name = "Cell type") +
    theme_bw() +
    labs(x = "Age group", y = ylab_text, title = short_title(title_text)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  save_plot_local(p, out_name, width = width, height = height)
  invisible(df2)
}

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
} else {
  message("[WARN] No regulator genes available; skip Figure 3C.")
}

## =========================================================
## Figure 4: Neurogenic trajectory vs ferroptosis activity (median along lineage)
## Panels:
##   A) Promoter module score
##   B) Inhibitor module score
##   C) Regulator genes (single-gene; one plot per gene)
## Outputs:
##   - Fig4A_Trajectory_PromoterScore.(png/pdf)
##   - Fig4B_Trajectory_InhibitorScore.(png/pdf)
##   - Fig4C_Trajectory_Regulator_<Gene>.(png/pdf) (each gene)
## =========================================================
neuro_order <- c("qNSC", "nIPC", "Neuroblast", "GC")

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

df_promoter_meta <- neuro@meta.data %>% select(Celltype_Article, Age_collapsed, all_of(promoter_col))
plot_traj_median(
  df_long    = df_promoter_meta,
  value_col  = promoter_col,
  title_text = "Figure 4A. Neurogenic trajectory vs promoter activity",
  ylab_text  = "Median promoter module score",
  out_name   = "Fig4A_Trajectory_PromoterScore",
  width      = 6.5, height = 4.5
)

df_inhib_meta <- neuro@meta.data %>% select(Celltype_Article, Age_collapsed, all_of(inhibitor_col))
plot_traj_median(
  df_long    = df_inhib_meta,
  value_col  = inhibitor_col,
  title_text = "Figure 4B. Neurogenic trajectory vs inhibitor activity",
  ylab_text  = "Median inhibitor module score",
  out_name   = "Fig4B_Trajectory_InhibitorScore",
  width      = 6.5, height = 4.5
)

if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    if (!g %in% rownames(neuro)) next
    df_g_meta <- neuro@meta.data %>% select(Celltype_Article, Age_collapsed)
    df_g_meta[[g]] <- FetchData(neuro, vars = g)[, 1]
    plot_traj_median(
      df_long    = df_g_meta,
      value_col  = g,
      title_text = paste0("Figure 4C. Neurogenic trajectory vs regulator ", g),
      ylab_text  = paste0("Median expression of ", g),
      out_name   = paste0("Fig4C_Trajectory_Regulator_", g),
      width      = 6.5, height = 4.5
    )
  }
} else {
  message("[WARN] No regulator genes available; skip Figure 4C.")
}

## =========================================================
## Figure 5: Differential expression volcano plots (Old vs Young) per cell type
## Outputs:
##   - Fig5_Volcano_<Celltype>_Old_vs_Young.(png/pdf)
##   - Table: Fig5_DE_<Celltype>_Old_vs_Young.csv
## =========================================================
run_de_and_volcano <- function(
    obj,
    celltype,
    group_col = "Age_collapsed",
    ident_col = "Celltype_Article",
    group1 = "Old",
    group2 = "Young",
    assay = "RNA",
    slot = "data",
    min.pct = 0.1,
    logfc.threshold = 0.0,
    test.use = "wilcox",
    fc_cutoff = 0.25,
    padj_cutoff = 0.05,
    pcap = 50,
    top_label_n = 15
) {
  DefaultAssay(obj) <- assay
  
  sub <- subset(obj, subset = !!as.name(ident_col) == celltype)
  sub[[group_col]] <- factor(sub[[group_col]][, 1], levels = c(group2, group1))
  sub <- subset(sub, subset = !is.na(!!as.name(group_col)) & !!as.name(group_col) %in% c(group1, group2))
  
  Idents(sub) <- sub[[group_col]][, 1]
  
  de <- FindMarkers(
    sub,
    ident.1 = group1,
    ident.2 = group2,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    assay = assay,
    slot  = slot
  )
  
  de$gene <- rownames(de)
  
  fc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(de) ~ "avg_log2FC",
    "avg_logFC"  %in% colnames(de) ~ "avg_logFC",
    "log2FC"     %in% colnames(de) ~ "log2FC",
    TRUE ~ NA_character_
  )
  if (is.na(fc_col)) stop("Cannot find logFC column in FindMarkers result.")
  
  padj_col <- if ("p_val_adj" %in% colnames(de)) "p_val_adj" else if ("p_adj" %in% colnames(de)) "p_adj" else NA_character_
  if (is.na(padj_col)) stop("Cannot find adjusted p-value column in FindMarkers result.")
  
  de <- de %>%
    mutate(
      avg_log2FC = .data[[fc_col]],
      p_adj      = .data[[padj_col]],
      neglog10p  = -log10(p_adj + 1e-300),
      neglog10p_capped = pmin(neglog10p, pcap),
      sig = ifelse(!is.na(p_adj) & p_adj < padj_cutoff & abs(avg_log2FC) >= fc_cutoff,
                   "Significant", "NS")
    )
  
  de_sig <- de %>%
    filter(sig == "Significant") %>%
    arrange(p_adj, desc(abs(avg_log2FC))) %>%
    head(top_label_n)
  
  p <- ggplot(de, aes(x = avg_log2FC, y = neglog10p_capped)) +
    geom_point(aes(color = sig), alpha = 0.8, size = 1.6) +
    scale_color_manual(values = c("NS" = "grey80", "Significant" = "#B2182B"), name = "") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = de_sig %>%
        mutate(neglog10p_capped_jitter = neglog10p_capped + runif(n(), -0.6, 0.6)),
      aes(x = avg_log2FC, y = neglog10p_capped_jitter, label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2,
      segment.size = 0.2,
      min.segment.length = 0,
      segment.color = "grey50",
      color = "black"
    ) +
    theme_bw() +
    labs(
      title = paste0("Figure 5. Volcano: ", celltype, " (Old vs Young)"),
      x = "log2FC (Old vs Young)",
      y = "-log10(adj. P) (capped)"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  list(de = de, de_sig = de_sig, plot = p)
}

neuro_types <- neuro_ct
volcano_results <- list()

for (ct in neuro_types) {
  res <- run_de_and_volcano(obj = neuro, celltype = ct)
  volcano_results[[ct]] <- res
  
  save_plot_local(res$plot, filename = paste0("Fig5_Volcano_", ct, "_Old_vs_Young"), width = 6.5, height = 5)
  
  out_csv <- file.path(fig_dir, paste0("Fig5_DE_", ct, "_Old_vs_Young.csv"))
  write.csv(res$de, out_csv, row.names = FALSE)
  message("[Saved] ", out_csv)
}

## =========================================================
## Figure 6: Correlation heatmaps (ferroptosis genes × neurogenesis markers)
## Outputs:
##   - Fig6A_CorrHeatmap_All.(png/pdf)
##   - Fig6B_CorrHeatmap_Young.(png/pdf)
##   - Fig6C_CorrHeatmap_Old.(png/pdf)
## =========================================================
build_corr_sub <- function(obj, ferro_genes, neuro_genes) {
  expr_df <- FetchData(obj, vars = c(ferro_genes, neuro_genes))
  
  sd_vec   <- apply(expr_df, 2, sd, na.rm = TRUE)
  expr_df2 <- expr_df[, sd_vec > 0, drop = FALSE]
  
  ferro_use2 <- intersect(ferro_genes, colnames(expr_df2))
  neuro_use2 <- intersect(neuro_genes, colnames(expr_df2))
  
  if (length(ferro_use2) == 0 || length(neuro_use2) == 0) {
    stop("After filtering sd>0, no overlapping ferro or neuro markers remain.")
  }
  
  corr_mat <- cor(expr_df2, method = "spearman", use = "pairwise.complete.obs")
  corr_sub <- corr_mat[ferro_use2, neuro_use2, drop = FALSE]
  
  row_keep <- apply(corr_sub, 1, function(x) any(!is.na(x)))
  col_keep <- apply(corr_sub, 2, function(x) any(!is.na(x)))
  corr_sub <- corr_sub[row_keep, col_keep, drop = FALSE]
  
  if (nrow(corr_sub) == 0 || ncol(corr_sub) == 0) {
    stop("corr_sub became empty after removing all-NA rows/cols.")
  }
  corr_sub
}

save_pheatmap_mat <- function(mat, base_name, main_title,
                              width = 6.5, height = 10,
                              cluster_rows = TRUE, cluster_cols = TRUE) {
  fig_pdf <- file.path(fig_dir, paste0(base_name, ".pdf"))
  fig_png <- file.path(fig_dir, paste0(base_name, ".png"))
  
  hm_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  ph <- pheatmap::pheatmap(
    mat,
    color        = hm_cols,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    border_color = NA,
    main         = main_title,
    silent       = TRUE
  )
  
  pdf(fig_pdf, width = width, height = height)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  dev.off()
  
  png(fig_png, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  dev.off()
  
  message("[Saved] ", fig_pdf)
  message("[Saved] ", fig_png)
  
  invisible(ph)
}

ferro_genes_all <- unique(unlist(ferroptosis_genes_use))
ferro_genes_all <- intersect(ferro_genes_all, rownames(neuro))

# neuro_markers_raw <- c("Sox2","Hes1","Hopx","Ascl1","Dcx","Neurod1","Prox1","Calb2","Nestin","Nes")
neuro_markers_raw <- c("Sox2","Hopx","Ascl1","Dcx","Neurod1","Prox1","Calb2","Nestin","Nes")
neuro_markers <- intersect(unique(neuro_markers_raw), rownames(neuro))
if (length(neuro_markers) == 0) stop("No specified neurogenesis markers found in neuro.")

corr_all <- build_corr_sub(neuro, ferro_genes_all, neuro_markers)
save_pheatmap_mat(
  mat        = corr_all,
  base_name  = "Fig6A_CorrHeatmap_All",
  main_title = "Figure 6A. Correlation (All cells): ferroptosis genes vs neurogenesis markers",
  width      = 6.5, height = 10
)

neuro_young <- subset(neuro, subset = Age_collapsed == "Young")
corr_young  <- build_corr_sub(neuro_young, ferro_genes_all, neuro_markers)
save_pheatmap_mat(
  mat        = corr_young,
  base_name  = "Fig6B_CorrHeatmap_Young",
  main_title = "Figure 6B. Correlation (Young): ferroptosis genes vs neurogenesis markers",
  width      = 6.5, height = 10
)

neuro_old <- subset(neuro, subset = Age_collapsed == "Old")
corr_old  <- build_corr_sub(neuro_old, ferro_genes_all, neuro_markers)
save_pheatmap_mat(
  mat        = corr_old,
  base_name  = "Fig6C_CorrHeatmap_Old",
  main_title = "Figure 6C. Correlation (Old): ferroptosis genes vs neurogenesis markers",
  width      = 6.5, height = 10
)

## =========================================================
## Figure 7: Fixed-order correlation heatmaps (All defines row/col order)
## Goal: Use an identical color scale across Fig7A/Fig7B/Fig7C
## Outputs:
##   - Fig7A_FixedOrder_CorrHeatmap_All.(png/pdf)
##   - Fig7B_FixedOrder_CorrHeatmap_Young.(png/pdf)
##   - Fig7C_FixedOrder_CorrHeatmap_Old.(png/pdf)
## =========================================================

library(pheatmap)
library(grid)

## ---------------- helper: save pheatmap (pdf + png) ----------------
save_pheatmap_obj <- function(ph_obj, base_name, width = 6.5, height = 10) {
  fig_pdf <- file.path(fig_dir, paste0(base_name, ".pdf"))
  fig_png <- file.path(fig_dir, paste0(base_name, ".png"))
  
  pdf(fig_pdf, width = width, height = height)
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable)
  dev.off()
  
  png(fig_png, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable)
  dev.off()
  
  message("[Saved] ", fig_pdf)
  message("[Saved] ", fig_png)
}

## ---------------- (1) unified colormap + unified breaks ----------------
# 颜色（蓝-白-红）
hm_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# 方案A（推荐、最审稿友好）：固定相关系数范围 [-1, 1]
global_min <- -0.5
global_max <-  0.5

# 若你更想按数据自适应范围，可改为：
# all_values  <- c(corr_all, corr_young, corr_old)
# global_min  <- max(-1, min(all_values, na.rm = TRUE))
# global_max  <- min( 1, max(all_values, na.rm = TRUE))

# breaks 数量必须比 colors 多 1
hm_breaks <- seq(global_min, global_max, length.out = length(hm_cols) + 1)

## ---------------- (2) Fig7A: All cells (defines clustering order) ----------------
ph_all_fixed <- pheatmap::pheatmap(
  corr_all,
  color        = hm_cols,
  breaks       = hm_breaks,      # <- key: fixed color scale
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main         = "Figure 7A. Fixed-order correlation (All cells)",
  silent       = TRUE
)

# 固定行列顺序（来自 All 的聚类）
row_order <- ph_all_fixed$tree_row$order
col_order <- ph_all_fixed$tree_col$order
ferro_order_fixed <- rownames(corr_all)[row_order]
neuro_order_fixed <- colnames(corr_all)[col_order]

save_pheatmap_obj(ph_all_fixed, "Fig7A_FixedOrder_CorrHeatmap_All", width = 6.5, height = 10)

## ---------------- (3) helper: fixed-order heatmap (Young/Old follow All order) ----------------
plot_fixed_order_heatmap <- function(mat, base_name, main_title,
                                     ferro_fixed, neuro_fixed,
                                     width = 6.5, height = 10,
                                     hm_cols, hm_breaks) {
  # 用 All 的行列顺序构造一个完整矩阵；缺失用 NA 填充，保证三图结构可比
  mat2 <- matrix(
    NA_real_,
    nrow = length(ferro_fixed),
    ncol = length(neuro_fixed),
    dimnames = list(ferro_fixed, neuro_fixed)
  )
  
  common_r <- intersect(rownames(mat), ferro_fixed)
  common_c <- intersect(colnames(mat), neuro_fixed)
  mat2[common_r, common_c] <- mat[common_r, common_c]
  
  ph <- pheatmap::pheatmap(
    mat2,
    color        = hm_cols,
    breaks       = hm_breaks,    # <- key: fixed color scale
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    main         = main_title,
    silent       = TRUE
  )
  
  save_pheatmap_obj(ph, base_name, width = width, height = height)
  invisible(mat2)
}

## ---------------- (4) Fig7B: Young ----------------
plot_fixed_order_heatmap(
  mat         = corr_young,
  base_name   = "Fig7B_FixedOrder_CorrHeatmap_Young",
  main_title  = "Figure 7B. Fixed-order correlation (Young)",
  ferro_fixed = ferro_order_fixed,
  neuro_fixed = neuro_order_fixed,
  width       = 6.5, height = 10,
  hm_cols     = hm_cols,
  hm_breaks   = hm_breaks
)

## ---------------- (5) Fig7C: Old ----------------
plot_fixed_order_heatmap(
  mat         = corr_old,
  base_name   = "Fig7C_FixedOrder_CorrHeatmap_Old",
  main_title  = "Figure 7C. Fixed-order correlation (Old)",
  ferro_fixed = ferro_order_fixed,
  neuro_fixed = neuro_order_fixed,
  width       = 6.5, height = 10,
  hm_cols     = hm_cols,
  hm_breaks   = hm_breaks
)

message("All figures completed (with a unified color scale).")

message("All figures completed.")


## =========================================================
## Script: Ferroptosis × Neurogenesis (qNSC / nIPC / NB / GC)
## Focus: Young vs Old
## Output: Fig4 (Promoter score / Inhibitor score / Regulator single-gene)
## =========================================================

library(Seurat)
library(dplyr)
library(ggplot2)

## ---------------- 0. Paths ----------------
setwd("H:/Projects/mmFerroptosis/stage-2")

data_dir <- "data/GSE233363"
fig_dir  <- file.path("figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------- helper: save plot (png + pdf) ----------------
save_plot_local <- function(p, filename,
                            width = 8, height = 5, dpi = 300) {
  ## PNG
  out_png <- file.path(fig_dir, paste0(filename, ".png"))
  ggsave(
    filename = out_png,
    plot     = p,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = "white"
  )
  
  ## PDF (vector, no dpi)
  out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(
    filename = out_pdf,
    plot     = p,
    width    = width,
    height   = height,
    device   = cairo_pdf   # strongly recommended to avoid font/transparency issues
  )
  
  message("[Saved] ", out_png)
  message("[Saved] ", out_pdf)
}

## ---------------- 0. Load object ----------------
combined <- readRDS(file.path(data_dir, "Seurat_combined_with_Celltype_Article.rds"))
DefaultAssay(combined) <- "RNA"

## ---------------- 0.1 Age collapsed (STRICT: Young vs Old only) ----------------
combined$Age <- combined$timepoint

## Collapse original timepoints into Young/Old; set all others to NA to prevent Middle from being plotted later as NA
combined$Age_collapsed <- NA_character_
combined$Age_collapsed[combined$Age %in% c("Young", "Young1", "Young2")] <- "Young"
combined$Age_collapsed[combined$Age %in% c("Old", "Old1", "Old2")]       <- "Old"

combined$Age_collapsed <- factor(combined$Age_collapsed, levels = c("Young", "Old"))

message("== Age table (combined) ==")
print(table(combined$Age, useNA = "ifany"))
message("== Age_collapsed table (combined) ==")
print(table(combined$Age_collapsed, useNA = "ifany"))

## ---------------- 1. Subset neurogenic lineage ----------------
neuro_ct <- c("qNSC", "nIPC", "Neuroblast", "GC")

neuro <- subset(combined, subset = Celltype_Article %in% neuro_ct)

## Keep only Young/Old and remove NA values so split.by does not create an NA group
neuro <- subset(neuro, subset = !is.na(Age_collapsed) & Age_collapsed %in% c("Young", "Old"))
neuro$Age_collapsed <- factor(neuro$Age_collapsed, levels = c("Young", "Old"))

neuro$Celltype_Article <- factor(neuro$Celltype_Article, levels = neuro_ct)

message("== Celltype × Age (neuro) ==")
print(table(neuro$Celltype_Article, neuro$Age_collapsed, useNA = "ifany"))

## Optional: quick UMAP check
p_umap <- DimPlot(
  neuro,
  reduction = "umap",
  group.by  = "Celltype_Article",
  split.by  = "Age_collapsed"
)
save_plot_local(p_umap, "UMAP_Neuro_Celltype_splitByAge", width = 10, height = 5)

## ---------------- 2. Ferroptosis gene sets ----------------
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
  ## Regulator: analyzed only at the single-gene level later; note that RAS is not a gene symbol, but keeping it is fine because it will be removed after intersection
  Ferro_Regulator = c("NQO1","VDAC2","TP53","RAS")
)

## human -> mouse symbol
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

## ---------------- 3. Colors ----------------
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

age_colors_A <- c("Young" = "#1F78B4", "Middle" = "#A6CEE3", "Old" = "#E31A1C")
age_colors_B <- c("Young" = "#4DBBD5", "Middle" = "#EFC000", "Old" = "#D73027")
age_colors_C <- c("Young" = "#3288BD", "Middle" = "#ABDDDE", "Old" = "#D53E4F")
age_palettes <- list("A_simple" = age_colors_A, "B_neuro" = age_colors_B, "C_nature" = age_colors_C)

## Use only the Young/Old two-color subset from the selected palette
age_colors <- age_palettes[["B_neuro"]][c("Young", "Old")]

## ---------------- 4. Fig4: Promoter & Inhibitor module scores ----------------
set.seed(1)  # ensure AddModuleScore is reproducible

## 4.1 Promoter score
if (!any(grepl("^FerroPromoter", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Promoter),
    name     = "FerroPromoter",
    assay    = "RNA"
  )
}
promoter_col <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]

p_vln_promoter <- VlnPlot(
  neuro,
  features = promoter_col,
  group.by = "Celltype_Article",
  split.by = "Age_collapsed",
  pt.size  = 0
) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Ferroptosis Promoter module score by cell type and age") +
  theme(plot.title = element_text(hjust = 0.5))

save_plot_local(p_vln_promoter, "Fig4_FerroPromoterScore_Violin_Celltype_byAge", width = 8, height = 5)

## 4.2 Inhibitor score
if (!any(grepl("^FerroInhibitor", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Inhibitor),
    name     = "FerroInhibitor",
    assay    = "RNA"
  )
}
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]

p_vln_inhibitor <- VlnPlot(
  neuro,
  features = inhibitor_col,
  group.by = "Celltype_Article",
  split.by = "Age_collapsed",
  pt.size  = 0
) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Ferroptosis Inhibitor module score by cell type and age") +
  theme(plot.title = element_text(hjust = 0.5))

save_plot_local(p_vln_inhibitor, "Fig4_FerroInhibitorScore_Violin_Celltype_byAge", width = 8, height = 5)

## ---------------- 5. Fig4: Regulator single-gene expression (violin) ----------------
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
  
  save_plot_local(p_vln_regulator, "Fig4_FerroRegulatorGenes_Violin_Celltype_byAge", width = 10, height = 5)
} else {
  message("[WARN] No regulator genes retained after intersection in neuro. Skip regulator violin.")
}

message("Done: Fig4 (Promoter/Inhibitor scores + Regulator genes).")


## =========================================================
## Fig 5: Summary line plots (mean ± SD) by Age × Celltype
## - Promoter score (module)
## - Inhibitor score (module)
## - Regulator genes (single-gene expression; each gene one plot)
## =========================================================


## ---------------- 5.1 helper: make line plot for a numeric column (NO error bar) ----------------
plot_line_summary <- function(df_long, value_col,
                              title_text,
                              ylab_text,
                              out_name,
                              width = 6, height = 5) {
  df2 <- df_long %>%
    select(all_of(c(value_col, "Age_collapsed", "Celltype_Article"))) %>%
    rename(Value = all_of(value_col)) %>%
    mutate(
      Age_collapsed = factor(Age_collapsed, levels = c("Young", "Old"))
    ) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      sd_value   = sd(Value, na.rm = TRUE),
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
    labs(
      x = "Age group",
      y = ylab_text,
      title = title_text
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  save_plot_local(p, out_name, width = width, height = height)
  return(invisible(df2))
}

## ---------------- 5.2 Promoter score line plot ----------------
# promoter_col was already generated above in Fig4
ferro_df_promoter <- FetchData(
  neuro,
  vars = c(promoter_col, "Age_collapsed", "Celltype_Article")
)
colnames(ferro_df_promoter)[1] <- promoter_col

plot_line_summary(
  df_long    = ferro_df_promoter,
  value_col  = promoter_col,
  title_text = "Promoter module score change from Young to Old",
  ylab_text  = "Mean promoter module score",
  out_name   = "Fig5_FerroPromoterScore_Line_Celltype_Young_vs_Old",
  width      = 6,
  height     = 5
)

## ---------------- 5.3 Inhibitor score line plot ----------------
# inhibitor_col was already generated above in Fig4
ferro_df_inhibitor <- FetchData(
  neuro,
  vars = c(inhibitor_col, "Age_collapsed", "Celltype_Article")
)
colnames(ferro_df_inhibitor)[1] <- inhibitor_col

plot_line_summary(
  df_long    = ferro_df_inhibitor,
  value_col  = inhibitor_col,
  title_text = "Inhibitor module score change from Young to Old",
  ylab_text  = "Mean inhibitor module score",
  out_name   = "Fig5_FerroInhibitorScore_Line_Celltype_Young_vs_Old",
  width      = 6,
  height     = 5
)

## ---------------- 5.4 Regulator genes: single-gene line plots ----------------
# Draw one Fig5 panel per regulator gene to avoid overcrowded line plots
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    df_g <- FetchData(
      neuro,
      vars = c(g, "Age_collapsed", "Celltype_Article")
    )
    colnames(df_g)[1] <- g
    
    plot_line_summary(
      df_long    = df_g,
      value_col  = g,
      title_text = paste0("Regulator gene ", g, " expression change from Young to Old"),
      ylab_text  = paste0("Mean expression of ", g),
      out_name   = paste0("Fig5_RegulatorGene_", g, "_Line_Celltype_Young_vs_Old"),
      width      = 6,
      height     = 5
    )
  }
} else {
  message("[WARN] No regulator genes available for Fig5 single-gene line plots.")
}

message("Done: Fig5 (Promoter/Inhibitor line plots + Regulator gene line plots).")


## ===================================================================
## Fig10. Neurogenic trajectory vs ferroptosis activity (line plot)
## - Promoter module score
## - Inhibitor module score
## - Regulator genes (single-gene expression; one plot per gene)
## ===================================================================

## ---------------- 10.0 neurogenic order ----------------
## If you already have neuro_order, remove this block; otherwise use the default order
if (!exists("neuro_order")) {
  neuro_order <- c("qNSC", "nIPC", "Neuroblast", "GC")
}

## ---------------- 10.1 helper: trajectory line plot (median by Celltype × Age) ----------------
plot_traj_median <- function(df_long, value_col,
                             title_text,
                             ylab_text,
                             out_name,
                             width = 6.5, height = 4.5) {
  dfA <- df_long %>%
    dplyr::filter(Celltype_Article %in% neuro_order) %>%
    mutate(
      Celltype_Article = factor(Celltype_Article, levels = neuro_order),
      Age_collapsed    = factor(Age_collapsed, levels = c("Young", "Old"))
    ) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(
      median_value = median(.data[[value_col]], na.rm = TRUE),
      .groups      = "drop"
    )
  
  p <- ggplot(
    dfA,
    aes(
      x     = Celltype_Article,
      y     = median_value,
      color = Age_collapsed,
      group = Age_collapsed
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = age_colors, name = "Age") +
    theme_bw() +
    labs(
      title = title_text,
      x     = "Neurogenic lineage",
      y     = ylab_text
    ) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      axis.title  = element_text(face = "bold"),
      panel.grid  = element_blank()
    )
  
  save_plot_local(p, out_name, width = width, height = height)
  return(invisible(dfA))
}

## ---------------- 10.2 Promoter score trajectory ----------------
df_promoter_meta <- neuro@meta.data %>%
  dplyr::select(Celltype_Article, Age_collapsed, all_of(promoter_col))

plot_traj_median(
  df_long    = df_promoter_meta,
  value_col  = promoter_col,
  title_text = "Fig10. Neurogenic trajectory vs ferroptosis promoter activity",
  ylab_text  = "Median promoter module score",
  out_name   = "Fig10_NeurogenicTrajectory_FerroPromoterScore",
  width      = 6.5,
  height     = 4.5
)

## ---------------- 10.3 Inhibitor score trajectory ----------------
df_inhibitor_meta <- neuro@meta.data %>%
  dplyr::select(Celltype_Article, Age_collapsed, all_of(inhibitor_col))

plot_traj_median(
  df_long    = df_inhibitor_meta,
  value_col  = inhibitor_col,
  title_text = "Fig10. Neurogenic trajectory vs ferroptosis inhibitor activity",
  ylab_text  = "Median inhibitor module score",
  out_name   = "Fig10_NeurogenicTrajectory_FerroInhibitorScore",
  width      = 6.5,
  height     = 4.5
)

## ---------------- 10.4 Regulator genes trajectory (single-gene) ----------------
## One Fig10 plot per regulator gene
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    if (!g %in% rownames(neuro)) next
    
    df_g_meta <- neuro@meta.data %>%
      dplyr::select(Celltype_Article, Age_collapsed)
    
    ## Pull gene expression from the RNA assay to avoid missing expression values in meta.data
    expr_g <- FetchData(neuro, vars = g)[, 1]
    df_g_meta[[g]] <- expr_g
    
    plot_traj_median(
      df_long    = df_g_meta,
      value_col  = g,
      title_text = paste0("Fig10. Neurogenic trajectory vs regulator gene ", g),
      ylab_text  = paste0("Median expression of ", g),
      out_name   = paste0("Fig10_NeurogenicTrajectory_Regulator_", g),
      width      = 6.5,
      height     = 4.5
    )
  }
} else {
  message("[WARN] No regulator genes available for Fig10 single-gene trajectories.")
}

message("Done: Fig10 (Promoter/Inhibitor trajectories + Regulator gene trajectories).")


## ---------------- 10.4 Volcano plot ----------------
library(ggrepel)

## =========================================================
## Volcano helper: DE (Old vs Young) within one cell type
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
    logfc.threshold = 0.0,  # do not pre-filter; the volcano plot uses its own threshold lines
    test.use = "wilcox",
    fc_cutoff = 0.25,
    padj_cutoff = 0.05,
    pcap = 50,              # upper cap for -log10(p_adj) to avoid extreme values dominating the scale
    top_label_n = 15        # maximum number of genes to label per plot (ranked by p_adj and |FC|)
) {
  DefaultAssay(obj) <- assay
  
  ## Subset to the current cell type only
  sub <- subset(obj, subset = !!as.name(ident_col) == celltype)
  
  ## Ensure the grouping column is a factor with only Young/Old
  sub[[group_col]] <- factor(sub[[group_col]][,1], levels = c(group2, group1))
  sub <- subset(sub, subset = !is.na(!!as.name(group_col)) & !!as.name(group_col) %in% c(group1, group2))
  
  ## Differential expression: Old vs Young
  Idents(sub) <- sub[[group_col]][,1]
  
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
  
  ## Handle FC column names across different Seurat versions
  fc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(de) ~ "avg_log2FC",
    "avg_logFC"  %in% colnames(de) ~ "avg_logFC",
    "log2FC"     %in% colnames(de) ~ "log2FC",
    TRUE ~ NA_character_
  )
  if (is.na(fc_col)) stop("Cannot find logFC column in FindMarkers result.")
  
  ## Handle p-value column names across different versions
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
  
  ## Select genes to label by taking the strongest significant hits
  de_sig <- de %>%
    filter(sig == "Significant") %>%
    arrange(p_adj, desc(abs(avg_log2FC))) %>%
    head(top_label_n)
  
  ## Volcano plot
  p <- ggplot(de, aes(x = avg_log2FC, y = neglog10p_capped)) +
    geom_point(aes(color = sig), alpha = 0.8, size = 1.6) +
    scale_color_manual(
      values = c("NS" = "grey80", "Significant" = "#B2182B"),
      name = ""
    ) +
    geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed", color = "grey60", linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed", color = "grey60", linewidth = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = de_sig %>%
        mutate(
          neglog10p_capped_jitter =
            neglog10p_capped + runif(n(), min = -0.6, max = 0.6)
        ),
      aes(
        x = avg_log2FC,
        y = neglog10p_capped_jitter,
        label = gene
      ),
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
      title = paste0("Volcano: ", celltype, " (Old vs Young)"),
      x = "log2FC (Old vs Young)",
      y = "-log10(adj. P) (capped)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  list(de = de, de_sig = de_sig, plot = p)
}

## =========================================================
## Run volcano for 4 neurogenic cell types
## =========================================================

neuro_types <- c("qNSC", "nIPC", "Neuroblast", "GC")

volcano_results <- list()

for (ct in neuro_types) {
  res <- run_de_and_volcano(
    obj = neuro,
    celltype = ct,
    group_col = "Age_collapsed",
    ident_col = "Celltype_Article",
    group1 = "Old",
    group2 = "Young",
    assay = "RNA",
    slot = "data",
    min.pct = 0.1,
    logfc.threshold = 0.0,
    fc_cutoff = 0.25,
    padj_cutoff = 0.05,
    pcap = 50,
    top_label_n = 15
  )
  
  volcano_results[[ct]] <- res
  
  save_plot_local(
    res$plot,
    filename = paste0("Volcano_", ct, "_Old_vs_Young"),
    width = 6.5,
    height = 5
  )
}

message("Done: Volcano plots for qNSC/nIPC/Neuroblast/GC (Old vs Young).")


## Save DE results for each cell type
for (ct in names(volcano_results)) {
  de_tbl <- volcano_results[[ct]]$de
  out_csv <- file.path(fig_dir, paste0("DE_", ct, "_Old_vs_Young.csv"))
  write.csv(de_tbl, out_csv, row.names = FALSE)
  message("[Saved] ", out_csv)
}

## ===================================================================
## Fig13. Correlation heatmap: ferroptosis genes × neurogenesis markers
## All / Young-only / Old-only
## ===================================================================
library(pheatmap)
library(grid)

## ---------------- 13.0 Ferroptosis genes (combine 3 categories) ----------------
## Use the previously obtained ferroptosis_genes_use object after intersection
## Merge Promoter + Inhibitor + Regulator gene sets and remove duplicates
ferro_genes_all <- unique(unlist(ferroptosis_genes_use))
ferro_genes_all <- intersect(ferro_genes_all, rownames(neuro))

if (length(ferro_genes_all) == 0) {
  stop("No ferroptosis genes found in rownames(neuro). Check ferroptosis_genes_use.")
}

## ---------------- 13.1 Neurogenesis markers (+ Nestin/Nes) ----------------
## The mouse Nestin gene symbol is typically Nes; also keep compatibility with direct use of Nestin
neuro_markers_raw <- c("Sox2","Hes1","Hopx","Ascl1","Dcx","Neurod1","Prox1","Calb2","Nestin","Nes")
neuro_markers <- unique(neuro_markers_raw)
neuro_markers <- intersect(neuro_markers, rownames(neuro))

if (length(neuro_markers) == 0) {
  stop("None of the specified neurogenesis markers are present in rownames(neuro).")
}

message("Neuro markers retained: ", paste(neuro_markers, collapse = ", "))

## ---------------- 13.2 Helper: build correlation submatrix ----------------
build_corr_sub <- function(obj, ferro_genes, neuro_genes) {
  ## 1) Expression matrix (by cell)
  expr_df <- FetchData(obj, vars = c(ferro_genes, neuro_genes))
  
  ## 2) Remove variables with sd == 0 to avoid all-NA correlations
  sd_vec   <- apply(expr_df, 2, sd, na.rm = TRUE)
  expr_df2 <- expr_df[, sd_vec > 0, drop = FALSE]
  
  ferro_use2 <- intersect(ferro_genes, colnames(expr_df2))
  neuro_use2 <- intersect(neuro_genes, colnames(expr_df2))
  
  if (length(ferro_use2) == 0 || length(neuro_use2) == 0) {
    stop("After filtering sd>0, no overlapping ferro or neurogenic markers remain.")
  }
  
  ## 3) Spearman correlation
  corr_mat <- cor(expr_df2, method = "spearman", use = "pairwise.complete.obs")
  
  ## 4) ferro × neuro submatrix
  corr_sub <- corr_mat[ferro_use2, neuro_use2, drop = FALSE]
  
  ## 5) Remove rows/columns that are entirely NA
  row_keep <- apply(corr_sub, 1, function(x) any(!is.na(x)))
  col_keep <- apply(corr_sub, 2, function(x) any(!is.na(x)))
  corr_sub <- corr_sub[row_keep, col_keep, drop = FALSE]
  
  if (nrow(corr_sub) == 0 || ncol(corr_sub) == 0) {
    stop("corr_sub became empty after removing all-NA rows/cols.")
  }
  
  return(corr_sub)
}

## ---------------- 13.3 Helper: save pheatmap to pdf+png ----------------
save_pheatmap_pdf_png <- function(mat, base_name,
                                  main_title,
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
  
  ## PDF
  pdf(fig_pdf, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
  
  ## PNG
  png(fig_png, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
  
  message("Saved: ", fig_pdf, " and ", fig_png)
  invisible(ph)
}

## ---------------- 13.4 (1) All cells (Young + Old combined) ----------------
corr_all <- build_corr_sub(neuro, ferro_genes_all, neuro_markers)

save_pheatmap_pdf_png(
  mat        = corr_all,
  base_name  = "Fig13_Corr_Ferro_vs_NeuroMarkers_All",
  main_title = "Fig13. Correlation: ferroptosis genes vs neurogenesis markers (All cells)",
  width      = 6.5,
  height     = 10
)

## ---------------- 13.5 (2) Young-only ----------------
neuro_young <- subset(neuro, subset = Age_collapsed == "Young")
corr_young  <- build_corr_sub(neuro_young, ferro_genes_all, neuro_markers)

save_pheatmap_pdf_png(
  mat        = corr_young,
  base_name  = "Fig13_Corr_Ferro_vs_NeuroMarkers_Young",
  main_title = "Fig13. Correlation: ferroptosis genes vs neurogenesis markers (Young)",
  width      = 6.5,
  height     = 10
)

## ---------------- 13.6 (3) Old-only ----------------
neuro_old <- subset(neuro, subset = Age_collapsed == "Old")
corr_old  <- build_corr_sub(neuro_old, ferro_genes_all, neuro_markers)

save_pheatmap_pdf_png(
  mat        = corr_old,
  base_name  = "Fig13_Corr_Ferro_vs_NeuroMarkers_Old",
  main_title = "Fig13. Correlation: ferroptosis genes vs neurogenesis markers (Old)",
  width      = 6.5,
  height     = 10
)

message("Done: Fig13 correlation heatmaps (All / Young / Old).")


## ===================================================================
## Fig13 (Fixed-order). Correlation heatmap: ferroptosis genes × neurogenesis markers
## All defines clustering order; Young/Old follow the same row/col order
## Output names: *_FixedOrder_*
## ===================================================================
library(pheatmap)
library(grid)

## ---------------- 13F.0 Ferro genes + markers (reuse from your pipeline) ----------------
ferro_genes_all <- unique(unlist(ferroptosis_genes_use))
ferro_genes_all <- intersect(ferro_genes_all, rownames(neuro))
if (length(ferro_genes_all) == 0) stop("No ferroptosis genes found in rownames(neuro).")

neuro_markers_raw <- c("Sox2","Hes1","Hopx","Ascl1","Dcx","Neurod1","Prox1","Calb2","Nestin","Nes")
neuro_markers <- unique(neuro_markers_raw)
neuro_markers <- intersect(neuro_markers, rownames(neuro))
if (length(neuro_markers) == 0) stop("No neurogenesis markers found in rownames(neuro).")

message("Neuro markers retained: ", paste(neuro_markers, collapse = ", "))

## ---------------- 13F.1 Helper: build correlation submatrix ----------------
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

## ---------------- 13F.2 Helper: save pheatmap to pdf+png ----------------
save_pheatmap_pdf_png <- function(ph_obj, base_name, width = 6.5, height = 10) {
  fig_pdf <- file.path(fig_dir, paste0(base_name, ".pdf"))
  fig_png <- file.path(fig_dir, paste0(base_name, ".png"))
  
  pdf(fig_pdf, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph_obj$gtable)
  dev.off()
  
  png(fig_png, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(ph_obj$gtable)
  dev.off()
  
  message("Saved: ", fig_pdf, " and ", fig_png)
}

## ---------------- 13F.3 (1) All: compute corr + cluster to get fixed order ----------------
corr_all <- build_corr_sub(neuro, ferro_genes_all, neuro_markers)

hm_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

ph_all <- pheatmap::pheatmap(
  corr_all,
  color        = hm_cols,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main         = "Fig13 (Fixed order). Correlation (All cells)",
  silent       = TRUE
)

## Extract a fixed order from the clustering tree of All
row_order <- ph_all$tree_row$order
col_order <- ph_all$tree_col$order
ferro_order_fixed <- rownames(corr_all)[row_order]
neuro_order_fixed <- colnames(corr_all)[col_order]

## Save All with clustering
save_pheatmap_pdf_png(
  ph_all,
  base_name = "Fig13_FixedOrder_Corr_Ferro_vs_NeuroMarkers_All",
  width = 6.5, height = 10
)

## ---------------- 13F.4 Helper: fixed-order pheatmap (no reclustering) ----------------
plot_fixed_order_heatmap <- function(mat, base_name, main_title,
                                     ferro_fixed, neuro_fixed,
                                     width = 6.5, height = 10) {
  ## Align to the fixed order: fill missing rows/columns with NA to keep the three plots structurally consistent
  mat2 <- matrix(NA_real_,
                 nrow = length(ferro_fixed),
                 ncol = length(neuro_fixed),
                 dimnames = list(ferro_fixed, neuro_fixed))
  
  common_r <- intersect(rownames(mat), ferro_fixed)
  common_c <- intersect(colnames(mat), neuro_fixed)
  mat2[common_r, common_c] <- mat[common_r, common_c]
  
  ph <- pheatmap::pheatmap(
    mat2,
    color        = hm_cols,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    main         = main_title,
    silent       = TRUE
  )
  
  save_pheatmap_pdf_png(ph, base_name = base_name, width = width, height = height)
  invisible(mat2)
}

## ---------------- 13F.5 (2) Young-only: fixed order ----------------
neuro_young <- subset(neuro, subset = Age_collapsed == "Young")
corr_young  <- build_corr_sub(neuro_young, ferro_genes_all, neuro_markers)

plot_fixed_order_heatmap(
  mat         = corr_young,
  base_name   = "Fig13_FixedOrder_Corr_Ferro_vs_NeuroMarkers_Young",
  main_title  = "Fig13 (Fixed order). Correlation (Young)",
  ferro_fixed = ferro_order_fixed,
  neuro_fixed = neuro_order_fixed,
  width       = 6.5,
  height      = 10
)

## ---------------- 13F.6 (3) Old-only: fixed order ----------------
neuro_old <- subset(neuro, subset = Age_collapsed == "Old")
corr_old  <- build_corr_sub(neuro_old, ferro_genes_all, neuro_markers)

plot_fixed_order_heatmap(
  mat         = corr_old,
  base_name   = "Fig13_FixedOrder_Corr_Ferro_vs_NeuroMarkers_Old",
  main_title  = "Fig13 (Fixed order). Correlation (Old)",
  ferro_fixed = ferro_order_fixed,
  neuro_fixed = neuro_order_fixed,
  width       = 6.5,
  height      = 10
)

message("Done: Fig13 fixed-order heatmaps (All / Young / Old).")


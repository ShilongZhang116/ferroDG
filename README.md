# ferroDG: Ferroptosis in Dentate Gyrus Neurogenesis

## Overview

ferroDG is an R-based analysis workflow for studying ferroptosis activity during
mouse dentate gyrus neurogenesis. It uses single-cell RNA-seq data to compare
Young vs Old groups across the neurogenic lineage (qNSC → nIPC → Neuroblast → GC).

## Key Features

- Data preprocessing for GSE233363 single-cell data
- Neurogenic lineage cell type identification
- Ferroptosis module scoring (promoter and inhibitor gene sets)
- Visualization: UMAP, violin plots, line plots, volcano plots, heatmaps
- Trajectory-style summaries across the neurogenic lineage

## Project Structure

```
ferroDG/
├── README.md                 # Chinese README
├── README_EN.md              # English README
├── LICENSE                   # License
├── .gitignore                # Git ignore rules
├── R/                        # Main R scripts
│   ├── 01_data_preprocessing.R
│   ├── 02_ferroptosis_analysis.R
│   └── 03_visualization.R
├── analysis/                 # Archived stage-based notebooks/scripts
│   ├── stage-1/
│   └── stage-2/
├── data/                     # Data directory
│   └── README.md
├── figures/                  # Figure outputs
│   ├── stage-1/
│   └── stage-2/
├── docs/                     # Documentation
│   ├── methodology.md
│   └── results.md
└── gene_lists/               # Gene lists used for scoring
    ├── ferroptosis_promoter.txt
    ├── ferroptosis_inhibitor.txt
    └── ferroptosis_regulator.txt
```

## Requirements

- R >= 4.0.0

### Core R packages

```r
Seurat
dplyr
tidyr
ggplot2
ggrepel
pheatmap
```

## Quick Start

### 1) Clone the repo

```bash
git clone https://github.com/shilongzhang/ferroDG.git
cd ferroDG
```

### 2) Prepare data

Place the Seurat object here:

```
data/
└── GSE233363/
    └── Seurat_combined_with_Celltype_Article.rds
```

### 3) Run the pipeline

```r
source("R/01_data_preprocessing.R")
source("R/02_ferroptosis_analysis.R")
source("R/03_visualization.R")
```

Or from the command line:

```bash
Rscript R/01_data_preprocessing.R
Rscript R/02_ferroptosis_analysis.R
Rscript R/03_visualization.R
```

### 4) Outputs

Figures are saved under `figures/` as PNG and PDF.

## Analysis Stages

### Stage 1: Initial analysis

- Global UMAPs
- Neurogenic lineage UMAPs
- Cell composition analysis
- Ferroptosis scoring and basic visualization

### Stage 2: Deep analysis

- Promoter/inhibitor module scores
- Regulator gene expression
- Differential expression (Old vs Young)
- Trajectory summaries
- Correlation analysis

## Gene Sets

### Promoters
TFRC, NCOA4, HMOX1, ACO1, IREB2, ACSL4, LPCAT3, ALOX5, ALOX12, ALOX15, PTGS2,
CHAC1, NOX1, NOX4, and others.

### Inhibitors
GPX4, SLC7A11, GCLC, GCLM, GSS, FTH1, NFS1, NFE2L2, KEAP1, and others.

### Regulators
NQO1, VDAC2, TP53, RAS.

## Citation

```bibtex
@software{ferroDG,
  title = { ... },
  author = { ... },
  year = {2026},
  url = {https://github.com/shilongzhang/ferroDG},
}
```

## Data Source

- GSE233363: single-cell RNA-seq from mouse dentate gyrus (Young vs Old)

## License

MIT. See `LICENSE`.

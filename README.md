# PCI Time-Course RNA-seq & Proteomics Analysis

R scripts for analyzing a mouse time-course RNA-seq experiment comparing PCI-treated and control samples at 0h, 4h, 8h, and 12h, with optional integration of proteomics data.

---

## Repository Contents

| File | Description |
|------|-------------|
| `pcsk.R` | Main analysis script: data loading, PCA, DESeq2 DEA, heatmaps, boxplots, proteomics integration |
| `pathway_pcsk.R` | GO/pathway enrichment dot plot visualization |

---

## Experiment Overview

- **Species:** Mouse (*Mus musculus*)
- **Design:** Time-course (0h, 4h, 8h, 12h) × Treatment (Control vs PCI)
- **Data types:** Bulk RNA-seq (count matrix) ± mass spectrometry proteomics

---

## Dependencies

Install Bioconductor packages once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("PCAtools", "Biobase", "DESeq2", "ComplexHeatmap"))
install.packages(c("ggplot2", "ggrepel", "ggpubr", "pheatmap",
                   "dplyr", "tidyr", "tibble", "readxl", "openxlsx",
                   "stringr", "circlize", "rstatix", "VennDiagram"))
```

---

## Script Structure (`pcsk.R`)

| Section | Content |
|---------|---------|
| 1 | Package installation & library loading |
| 2 | Data loading & preprocessing (count matrix + metadata) |
| 3 | PCA |
| 4 | Group construction (Time × Condition labels) |
| 5 | ComplexHeatmap for candidate gene sets |
| 6 | Sample-to-sample correlation heatmap (pheatmap) |
| 7 | DESeq2 DEA (stratified by time point) |
| 8 | Volcano plots |
| 9 | Violin / Boxplot functions (final version with Wilcoxon + significance brackets) |
| 10 | Gene lists for plotting |
| 11 | Global sample correlation analysis |
| 12 | Cross-group comparison: 8h PCI vs 12h CTL |
| 13 | Proteomics + RNA-seq integration (heatmap + scatter) |
| 14 | GO Biological Process bar plots from proteomics annotation |

---



---

## Citation / Contact
Proprotein convertase activity regulates cumulus-oocyte-complex matrix integrity and cumulus cell migration during ovulation via a GDF9-dependent mechanism

Caroline E. Kratka, et al. 

Question please email: ruixu.huang0405@gmail.com

Relevant data is in GSE331136

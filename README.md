# Differential Gene Expression Analysis with DESeq2 (Airway Dataset)

**Question:** Which genes change their expression when human airway smooth muscle cells are treated with dexamethasone, a corticosteroid used to treat asthma?

This project performs a complete bulk RNA-seq differential expression workflow in R using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), from raw counts to quality control, statistical testing, gene annotation and visualisation.

## Dataset

- **Source:** the [`airway`](https://bioconductor.org/packages/release/data/experiment/html/airway.html) Bioconductor package (Himes et al., 2014, *PLoS ONE*; GEO: GSE52778)
- **Samples:** 8 RNA-seq samples from 4 human airway smooth muscle cell lines, each untreated and treated with dexamethasone (paired design)
- **Input:** raw gene-level read counts (~63,000 Ensembl genes)

## Methods

1. **Data preparation:** extracted the count matrix and sample metadata from the `airway` package
2. **Filtering:** removed genes with fewer than 10 total reads across samples
3. **Normalisation and testing:** DESeq2 (median-of-ratios normalisation, negative binomial model, Wald test), with "untreated" as the reference level
4. **Quality control:** variance-stabilising transformation (VST), PCA and sample-to-sample distance heatmap
5. **Annotation:** mapped Ensembl IDs to gene symbols with `org.Hs.eg.db`
6. **Visualisation:** MA plot and volcano plot (thresholds: adjusted p < 0.05 and |log2 fold change| > 1)

## Results

**Treatment has a clear effect.** On the PCA plot, treated and untreated samples separate into two distinct groups based on gene expression alone.

<img src="pca_plot.png" width="550">

**872 genes changed significantly.** Of about 18,400 genes tested, 872 passed the thresholds: 475 were up-regulated and 397 down-regulated after treatment.

<img src="volcano_plot.png" width="600">

**The results agree with known biology.** Several top genes, including *PER1*, *DUSP1*, *ZBTB16* and *MAOA*, are known glucocorticoid-responsive genes, which supports the validity of the analysis.

**Limitation.** The sample-distance heatmap shows that part of the variation comes from the cell line (donor), not only the treatment. A design that accounts for this, `design = ~ cellLine + dexamethasone`, would be the natural next step.

<img src="sample_distance_heatmap.png" width="500">

## Files

| File | Description |
|---|---|
| `script_for_airwaypackage.R` | Extracts `counts_data.csv` and `sample_info.csv` from the `airway` package |
| `DESeq2_project.R` | Main analysis: builds the DESeqDataSet, filters genes, runs DESeq2 and summarises results (MA plot) |
| `sample_qc_and_volcano.R` | PCA, sample-distance heatmap, gene annotation and volcano plot |
| `deseq2_results_annotated.csv` | Full results table with gene symbols |
| `Rplot_DESeq2.pdf` | MA plot |

## How to run

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "airway", "org.Hs.eg.db", "AnnotationDbi"))
install.packages(c("tidyverse", "pheatmap", "ggrepel", "RColorBrewer"))

source("script_for_airwaypackage.R")   # generates counts_data.csv and sample_info.csv
source("DESeq2_project.R")             # runs the DESeq2 analysis
source("sample_qc_and_volcano.R")      # QC plots, annotation and volcano plot
```

## Skills demonstrated

Bulk RNA-seq analysis · DESeq2 · normalisation and variance-stabilising transformation · PCA · multiple-testing correction · gene annotation · data visualisation with ggplot2 · reproducible R workflows

## Related project

[Data wrangling and visualisation of breast cancer gene expression (GSE183947)](https://github.com/Sayali-k28/gse183947-data-wrangling-visualization)

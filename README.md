# ATManalysis
# Nerve-associated macrophages control adipose homeostasis across lifespan and restrain age-related inflammation

This repository contains the full analysis code and notebooks used to generate the single-cell RNA sequencing results presented in the study:  
**"Nerve-associated macrophages control adipose homeostasis across lifespan and restrain age-related inflammation."**

## Overview

The repository includes code for:

- Single-cell RNA-sequencing preprocessing, dimensionality reduction, clustering, and cell type annotation
- Differential gene expression analysis between conditions
- Gene set enrichment analysis (GSEA)
- Gene set enrichment scoring (GSES)
- Differential abundance testing using scProportional Test
- RNA velocity analysis

All files are provided to enable full reproducibility and transparency of the analyses performed in this study.

---

## Files

- **Notebook SC1** — R notebook for single-cell RNA-sequencing preprocessing, quality control, dimensionality reduction, clustering, and cell type identification.
- **Notebook SC2** — R notebook for gene set enrichment analysis (GSEA) and single-cell gene set enrichment scoring (GSES) using clusterProfiler, msigdbr, and escape.
- **Notebook SC3** — R notebook for differential abundance testing using scProportional Test.
- **Notebook SC4** — Python notebook for RNA velocity inference using the steady state model (scVelo).

---

## Software and Packages

### R Packages
- Seurat v4.3.0
- clusterProfiler v4.11.0
- msigdbr v7.5.1
- escape v1.10.0
- scProportionTest

### Python Packages
- scVelo v0.2.5
- Velocyto

---

## Analysis Workflow

### 1️⃣ Quality Control and Preprocessing
- Cells expressing fewer than 363 genes or with mitochondrial gene content exceeding 10% were excluded.
- Data were normalized using Log-Normalization (scale factor = 10,000).
- Top variable genes were identified using Seurat’s `mean.var.plot` method.

### 2️⃣ Dimensionality Reduction & Clustering
- PCA was performed and the top 20 principal components were retained.
- Shared Nearest Neighbor (SNN) graph constructed with default parameters.
- Louvain clustering was performed at resolution 0.5.

### 3️⃣ Cell Type Annotation
- Cluster-specific markers identified using Wilcoxon Rank Sum test (min.pct = 0.1, logfc.threshold = 0.1).
- Clusters 9, 12, and 14 were excluded based on absence of Ptprc (CD45), Itgam (CD11b), and Adgre1 (F4/80).

### 4️⃣ Differential Gene Expression
- Wilcoxon Rank Sum test used to identify differentially expressed genes between Old and Young conditions (min.cell = 25, min.pct = 0.01, logfc.threshold = log(2)).

### 5️⃣ Enrichment Analysis
- Performed using clusterProfiler with one-sided hypergeometric test.
- Benjamini-Hochberg correction used for multiple testing.
- Gene Ontology (GO) terms (C5: GO:BP, GO:CC, GO:MF) and Hallmark gene sets (category H) were retrieved from msigdbr for *Mus musculus*.
- p-value and q-value cutoffs were set to 0.05.

### 6️⃣ Gene Set Enrichment Scoring (GSES)
- Calculated using the `enrichIt` function in escape (method = UCell).
- Gene sets analyzed included:
  - SenMayo (SAUL_SEN_MAYO, GSEA: MM16098)
  - NAMs (LEVEAU_NERVE_ASSOCIATED_MACROPHAGES, Table S3)

### 7️⃣ RNA Velocity
- Velocyto used to generate spliced/unspliced matrices.
- scVelo (steady-state model) used for velocity inference.
- Genes filtered with a minimum shared count of 20; top 5000 highly variable genes used.
- Velocities projected onto UMAP embeddings for visualization.

---

## Contact

For questions or additional information, please contact Elsie Gonzalez-Hurtado at elsie.gonzalez@icloud.com.

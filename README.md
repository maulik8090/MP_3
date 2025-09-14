# MP_3: Single-Nucleus RNA-seq Analysis of Mouse Kidney Samples

## 📘 Project Overview

This repository contains the analysis pipeline for single-nucleus RNA sequencing (snRNA-seq) data from mouse kidney samples. The study investigates the transcriptional differences across age and genotype using Seurat-based processing and differential gene expression analysis.

---

## 🧪 Sample Description

The dataset includes 12 mouse kidney samples grouped as follows:

- **3-month Wildtype (WT):** Smpl1, Smpl2, Smpl9  
- **3-month Cpt1a Knockout (KO):** Smpl3, Smpl4, Smpl5  
- **24-month Wildtype (WT):** Smpl6, Smpl7, Smpl8  
- **24-month Cpt1a Knockout (KO):** Smpl10, Smpl11, Smpl12  

All samples were processed and merged into a single Seurat object (`ahfbase`) for downstream analysis.

---

## 🧬 Analysis Workflow

### 1. **Preprocessing**
- Read and merge `.h5` or `.h5Seurat` files using `Read10X_h5()` or `LoadH5Seurat()`.
- Create Seurat objects for each sample.
- Annotate sample metadata and merge into a unified object.

### 2. **Quality Control**
- Filter cells based on mitochondrial, hemoglobin, ribosomal, and lncRNA content.
- Remove predicted doublets using `DoubletFinder`.

### 3. **Normalization and Feature Selection**
- Normalize data using `LogNormalize`.
- Identify highly variable genes.
- Scale data and perform PCA.

### 4. **Batch Correction**
- Apply Harmony integration using `RunHarmony()` on the `Sample` metadata.

### 5. **Clustering and Visualization**
- Perform clustering using `FindNeighbors()` and `FindClusters()`.
- Visualize clusters using UMAP (`RunUMAP()` and `DimPlot()`).
- Annotate kidney-specific cell types using marker genes.

### 6. **Differential Gene Expression (DGE)**
- Compare **24-month KO vs WT** within each kidney cell type using `FindMarkers()`.
- Filter significant genes (adjusted p-value < 0.05).
- Save DGE results per cell type.

### 7. **Ontology Enrichment**
- Use top upregulated and downregulated genes for pathway analysis.
- Enrichment performed using:
  - `clusterProfiler` (GO Biological Process)
  - `EnrichR` with databases:
    - KEGG_2019_Mouse
    - GO_Molecular_Function_2023
    - GO_Cellular_Component_2023
    - GO_Biological_Process_2023
    - WikiPathways_2019_Mouse

---

## 📦 Packages Used

- [Seurat](https://satijalab.org/seurat/)
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
- [EnrichR](https://cran.r-project.org/web/packages/enrichR/)
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- [Harmony](https://github.com/immunogenomics/harmony)

---

## 📁 Repository Contents

- `auto_git_push.R` – Script to automate Git updates.
- `DGE_*.csv` – Differential expression results per cell type.
- `GO_BP_*.csv` – GO enrichment results.
- `Enrichr_*.csv` – EnrichR pathway results.
- `README.md` – Project documentation.

---

## 📞 Contact

For questions or contributions, please contact [Maulik](https://github.com/maulik8090).


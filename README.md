
# 🧬 Breast Cancer Transcriptome Insights : A deep dive into Tumor and Normal Gene Expression" 
### A Deep Dive into Tumor and Normal Gene Expression

This project explores differential gene expression (DGE) between **tumor** and **normal** breast tissues using microarray data from the **GSE15852** dataset. It highlights significant genes, enriched biological processes, and pathways involved in breast cancer using R-based tools.

---

## 📌 Objectives

- Identify differentially expressed genes (DEGs) between tumor and normal breast tissues
- Visualize gene expression patterns using PCA, heatmaps, and volcano plots
- Perform GO and KEGG pathway enrichment analyses
- Conduct GSEA (Gene Set Enrichment Analysis)
- Generate interpretable biological insights for breast cancer

---

## 🧰 Tools & Packages Used

| Type               | Packages / Tools                        |
|--------------------|-----------------------------------------|
| Data Source        | GEO (GSE15852)                          |
| DGE Analysis       | `limma`, `GEOquery`                     |
| Annotation         | `org.Hs.eg.db`, `biomaRt`              |
| Enrichment         | `clusterProfiler`, `enrichplot`        |
| Visualization      | `pheatmap`, `EnhancedVolcano`, `ggplot2`|

---

## 📁 Directory Structure

breast-cancer-dge/
├── data/ # Filtered expression matrix (optional)
├── results/ # CSVs for DEGs and enrichment
├── plots/ # PNGs for volcano, heatmap, PCA, etc.
├── scripts/ # R scripts for each step
│ ├── dge_analysis.R
│ ├── enrichment_analysis.R
│ └── qc_plots.R
├── breast_cancer_dge.Rmd # R Markdown main report
├── README.md # This file
├── sessionInfo.txt # R environment details
└── .gitignore

yaml
Copy
Edit

---



### 🧬 Step-by-Step Workflow

#### 1. 📥 Data Acquisition
- Dataset: **GSE15852** from NCBI GEO
- Platform: **GPL570** (Affymetrix Human Genome U133 Plus 2.0 Array)
- Samples: 43 tumor & 43 matched normal breast tissues
- Tool used: `GEOquery` (R/Bioconductor)

#### 2. 🧹 Data Preprocessing
- Extracted expression matrix from GEO object
- Checked and handled missing values
- Applied log2 transformation (if needed)
- Defined sample groups for comparison (Tumor vs Normal)

#### 3. 📊 Quality Control & Exploratory Data Analysis
- **PCA Plot**: Visualized sample clustering
- **Heatmap**: Top 50 most variable genes
- **Boxplots** (optional): Checked normalization across samples

#### 4. 📈 Differential Gene Expression (DGE) Analysis
- Tool: `limma`
- Designed model matrix and fitted linear model
- Used empirical Bayes method to calculate statistics
- Selected DEGs with:
  - |log2FoldChange| ≥ 1  
  - Adjusted p-value < 0.05
- Output: `results/deg_results.csv`

#### 5. 🧠 Functional Enrichment Analysis
- Tool: `clusterProfiler`
- **GO Enrichment (Biological Process)**:
  - Identified overrepresented functions in DEGs
- **KEGG Pathway Analysis**:
  - Mapped DEGs to biological pathways
- Outputs:
  - `results/go_enrichment.csv`
  - `results/kegg_enrichment.csv`

#### 6. 🔍 Gene Set Enrichment Analysis (GSEA)
- Tool: `gseGO()` from `clusterProfiler`
- Used ranked gene list (all genes)
- Identified coordinated gene expression patterns
- Visuals:
  - `plots/gsea_GO_top1.png`
  - Dotplots of enriched pathways

#### 7. 🔗 Network and Pathway Visualization
- Tool: `cnetplot()` from `enrichplot`
- Network plot of gene-function relationships
- Helped visualize how DEGs contribute to enriched biological functions

---

### 🧾 Summary of Tools and Packages

| Step                     | R Package(s) Used                           |
|--------------------------|---------------------------------------------|
| Data Retrieval           | `GEOquery`                                  |
| DGE Analysis             | `limma`                                     |
| Annotation               | `org.Hs.eg.db`, `biomaRt`                   |
| Enrichment Analysis      | `clusterProfiler`, `enrichplot`             |
| Visualization            | `ggplot2`, `pheatmap`, `EnhancedVolcano`   |

---
## 🚀 How to Reproduce

1. Clone this repository:
   ```bash
   git clone https://github.com/Asthayadav09/Breast-cancer--normal-and-tumor/new/main?filename=README.md

Install required packages:

r
Copy
Edit
install.packages("BiocManager")
BiocManager::install(c(
  "limma", "GEOquery", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "pheatmap", "EnhancedVolcano"
))


📊 Outputs
Output Type	File Example
DEGs Table	results/deg_results.csv
GO Enrichment	results/go_enrichment.csv
KEGG Pathways	results/kegg_enrichment.csv
Volcano Plot	plots/volcano_plot.png
PCA Plot	plots/pca_plot.png
GSEA Plot	plots/gsea_GO_top1.png
Network Plot	plots/go_network_cnetplot.png


📚 Dataset Reference
GSE15852: NCBI GEO

Platform: GPL570

Samples: 43 tumor + 43 matched normal breast tissue

🤝 Acknowledgements
Special thanks to:

The NCBI GEO team for providing open-access data

R/Bioconductor for open-source bioinformatics tools

🧠 Author
Astha Yadav
📧 asthayadav0904@gmail.com

📌 License
This project is open-source and available under the MIT License.

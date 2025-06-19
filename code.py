# Install Bioconductor and required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db", "annotate", "pheatmap", "EnhancedVolcano"))

# Load libraries
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(annotate)
library(pheatmap)
library(EnhancedVolcano)

# Download data from GEO
gse_list <- getGEO("GSE15852", GSEMatrix = TRUE)
gse <- gse_list[[1]]

# Expression data and phenotype data
expr <- exprs(gse)
pheno <- pData(gse)

# Check relevant metadata columns
colnames(pheno)
unique(pheno$`source_name_ch1`)  # Helps find sample condition

# Assign groups
pheno$group <- ifelse(grepl("tumor", pheno$`source_name_ch1`, ignore.case = TRUE), "Tumor",
                      ifelse(grepl("normal", pheno$`source_name_ch1`, ignore.case = TRUE), "Normal", NA))

# Keep only samples labeled Tumor or Normal
keep <- !is.na(pheno$group)
expr <- expr[, keep]
pheno <- pheno[keep, ]
pheno$group <- factor(pheno$group)

# Confirm group counts
table(pheno$group)

# Create design matrix
design <- model.matrix(~ 0 + pheno$group)
colnames(design) <- levels(pheno$group)

# Fit linear model
fit <- lmFit(expr, design)

# Define contrast: Tumor vs Normal
contrast.matrix <- makeContrasts(TumorVsNormal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract top results
results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

# Add gene symbols from annotation DB
results$Symbol <- getSYMBOL(rownames(results), "hgu133plus2.db")
results <- results[!is.na(results$Symbol), ]  # Remove NAs

EnhancedVolcano(results,
                lab = results$Symbol,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Breast Tumor vs Normal (GSE15852)",
                subtitle = "Volcano Plot of Differential Expression")
# Install if needed
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot: Breast Tumor vs Normal",
       x = "Principal Component 1",
       y = "Principal Component 2")


# Install if not already done
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)

# Convert gene symbols to ENTREZ IDs
gene_df <- bitr(results$Symbol, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Merge ENTREZ IDs with results
results$ENTREZID <- gene_df$ENTREZID[match(results$Symbol, gene_df$SYMBOL)]

# Create named vector for GSEA
geneList <- results$logFC
names(geneList) <- results$ENTREZID
geneList <- sort(na.omit(geneList), decreasing = TRUE)

gsea_go <- gseGO(geneList = geneList,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

gseaplot2(gsea_go, geneSetID = 1, title = gsea_go@result$Description[1])

# Use significantly DEGs for enrichment (adj.P.Val < 0.05, logFC > 1)
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
entrez_sig <- bitr(sig_genes$Symbol, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO BP
go_enrich <- enrichGO(gene = entrez_sig$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05)

# KEGG
kegg_enrich <- enrichKEGG(gene = entrez_sig$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

# Barplot
barplot(go_enrich, showCategory = 10, title = "Top GO BP Terms")
barplot(kegg_enrich, showCategory = 10, title = "Top KEGG Pathways")

# Dotplot
dotplot(go_enrich, showCategory = 10, title = "GO Enrichment Dotplot")
dotplot(kegg_enrich, showCategory = 10, title = "KEGG Dotplot")

# Network plot for GO
cnetplot(go_enrich, showCategory = 5, foldChange = geneList, circular = FALSE)

# For KEGG
cnetplot(kegg_enrich, showCategory = 5, foldChange = geneList)

ggsave("GSEA_TopTerm.png", plot = last_plot(), width = 8, height = 6)


write.csv(results, "Breast_Tumor_vs_Normal_DEG_results.csv")


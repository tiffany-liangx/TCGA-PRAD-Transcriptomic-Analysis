
# Create a hierarchical clustering of TCGA-PRAD dataset.

Top 100 genes are filtered based on fold change and FDR threshold resulting in 88 genes DE genes for hierarchical clustering: 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

#colnames(colData(SECoding))
limma_res = limma_pipeline(
  tcga_data=SECoding,
  condition_variable="sample_type",
  reference_group="Solid Tissue Normal")

#outputs:

# 1. DE analysis- all genes
allgenelist<- limma_res$allGenes
all_gene_names<- allgenelist$gene_name

# 2. normalized expression data
expr <- as.data.frame(limma_res$voomObj$E)
expr$gene_id <- rownames(expr)

# 3. top 100 DE genes
top100genes<-limma_res$topGenes


#Set thresholds for FDR and fold-change
fdr_threshold <- 0.05
lfc_threshold <- 1

# filter for subset genes that pass both thresholds
sig_fc <- top100genes[top100genes$adj.P.Val < fdr_threshold &  abs(top100genes$logFC) > lfc_threshold, ]
sorted_fc <- sig_fc[order(-sig_fc$logFC), ]

#extract normalized expression data from filtered top 100 genes 
sig_ids <- rownames(sorted_fc)
expr_sig <- expr[sig_ids,]


```


```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#set up data input for hclust. clusters will be labeled with gene names. 
w_df <- sorted_fc[, c("gene_id", "gene_name")]
sig_genes_expr <- merge(w_df, expr, by = "gene_id")
rownames(sig_genes_expr)<-sig_genes_expr$gene_name
sig_genes_expr<-sig_genes_expr[,3:106]


h<-hclust(dist(sig_genes_expr))
plot(h)
```

After filtering, the genes that passed both the FDR (< 0.05) and fold-change (|logFC| > 1) thresholds were ranked by decreasing log fold change, and their expression values were pulled across all samples. Hierarchical clustering was then performed using Euclidean distance, which groups genes based on how similar their overall expression patterns are. Genes that merge at lower branch heights tend to show very similar expression behavior, while those joining at higher heights are more distinct.

What stands out is that the dendrogram forms a few major branches, suggesting that these significant genes naturally group into broader expression modules. This implies that the changes between tumor and normal tissue are not happening one gene at a time, but rather through coordinated shifts in sets of genes. These clusters likely reflect shared biological programs, such as pathways involved in proliferation, metabolism, or loss of normal tissue function.

The remaining genes are then divided based on fold change, where a positive logFC indicates upregulation and a negative logFC represents downregulation.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#Set thresholds for FDR and fold-change
fdr_threshold <- 0.05
lfc_threshold <- 1

# filter for subset genes that pass both thresholds
sig_fc <- top100genes[top100genes$adj.P.Val < fdr_threshold &  abs(top100genes$logFC) > lfc_threshold, ]
sorted_fc <- sig_fc[order(-sig_fc$logFC), ]
head(sorted_fc)

pos_logfc<-sorted_fc[sorted_fc$logFC>0, ] #upregulation
neg_logfc<-sorted_fc[sorted_fc$logFC<0, ] #down-regulation 


```


# Annotate each list with GO Annotations and KEGG.

Used reference link https://yulab-smu.top/biomedical-knowledge-mining-book/021-go.html accessed via vignette("clusterProfiler") in R for GO and KEGG enrichment analysis. 

GO enrichment requires Entrez IDs as input for the enrichGO() function. The code first loads the necessary annotation and enrichment packages, then retrieves the SYMBOL-to-ENTREZID mapping from org.Hs.eg.db for all genes in the dataset. The goal is to convert the gene symbols from the upregulated and downregulated gene lists into their corresponding Entrez IDs to ensure proper annotation during enrichment analysis. Next, enrichGO() performs an over-representation analysis separately for the upregulated and downregulated gene sets across all Gene Ontology categories (Biological Process, Molecular Function, and Cellular Component), using Benjamini–Hochberg (BH) correction to control for multiple testing. The results are visualized using dot plots that display the top 10 enriched GO terms for each group. This approach shifts the interpretation from individual gene-level differences to pathway-level insights to identify broader biological themes that distinguish tumor from normal samples.

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
library(org.Hs.eg.db)
library(data.table)
library(clusterProfiler)
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- unlist(as.list(x[mapped_genes]))
go <- data.table(SYMBOL=xx, ENTREZID=names(xx))
go <- go[!duplicated(go$SYMBOL),]
go <- go[go$SYMBOL %in% all_gene_names,]
head(go)
dim(go)


pos_entrez <- go[SYMBOL %in% pos_logfc$gene_name]$ENTREZID
neg_entrez <- go[SYMBOL %in% neg_logfc$gene_name]$ENTREZID
length(pos_entrez)
length(neg_entrez)
all_entrez <- go$ENTREZID

head(pos_entrez)
head(all_entrez)


library(enrichplot)
ego_poslogfc <- enrichGO(
  gene          = pos_entrez,
  universe      = all_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff = 1,
  readable      = TRUE
)
dotplot(ego_poslogfc, showCategory = 10, font.size = 10)+ ggtitle("up-regulated genes- GO")

ego_neglogfc <- enrichGO(
  gene          = neg_entrez,
  universe      = all_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff = 1,
  readable      = TRUE
)

dotplot(ego_neglogfc, showCategory = 10, font.size = 10)+ ggtitle("down-regulated genes- GO")



```

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

ekegg_poslogfc <- enrichKEGG(
  gene          = pos_entrez,
  universe      = all_entrez,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1
)

dotplot(ekegg_poslogfc, showCategory = 10, font.size = 10)+ ggtitle("up-regulated genes- KEGG")

ekegg_neglogfc <- enrichKEGG(
  gene          = neg_entrez,
  universe      = all_entrez,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1
)

dotplot(ekegg_neglogfc, showCategory = 10, font.size = 10)+ ggtitle("down-regulated genes- KEGG")
```
The enrichment analysis suggests that tumor-associated DEGs are strongly enriched in developmental and differentiation-related programs. GO terms for upregulated genes include pattern specification, regionalization, and neuronal/glial differentiation indicating that tumor cells may be reactivating embryonic lineage programs. For instance, pattern specification has the highest count from go enrichment analysis. It refers to developmental processes that establish defined spatial identities that guide cell differentiation [1]. This is consistent with evidence that cancer cells reactivate embryonic signaling pathways, such as Wnt signaling, which are normally inactive in differentiated tissues [2]. In prostate cancer, loss or epigenetic silencing of the tumor suppressor APC disrupts β-catenin degradation and activating β-catenin mutations have shown to further drive early tumor development [2]. Together, the dotplot results support the idea that tumors adopt developmental and stem-like transcriptional states. In this dataset, GO enrichment appears more statistically robust than KEGG, as reflected by consistently lower q-values. KEGG enrichment for upregulated genes is weaker and less statistically robust, with higher adjusted p-values and small gene counts. Based on this, upregulated tumor transcriptome appears more characterized by development and cell identity, rather than metabolic pathways. 

For the downregulated genes, GO enrichment points to processes related to membrane organization and vesicular trafficking, including receptor-mediated endocytosis and membrane raft assembly. This suggests that tumor tissue may have altered receptor internalization and disrupted membrane dynamics compared to normal tissue. Supporting this, Boyu et al. reported increased clathrin expression and reduced caveolin-1 (Cav-1) in prostate cancer, with enhanced clathrin-mediated EGFR recycling correlating with increased tumor aggressiveness [3]. Although Cav-1 can function as a tumor suppressor in certain contexts, it has also been associated with malignant progression depending on the cellular environment [3]. In addition, KEGG enrichment highlights reduced detoxification and metabolic pathways, such as cytochrome P450, xenobiotic metabolism, and glutathione metabolism. These findings suggest that tumors may exhibit both altered receptor signaling regulation and diminished normal metabolic capacity relative to normal tissue.

References:
[1] https://amigo.geneontology.org/amigo/term/GO:0007389 

[2] da Silva, Henrique B et al. “Dissecting Major Signaling Pathways throughout the Development of Prostate Cancer.” Prostate cancer vol. 2013 (2013): 920612. doi:10.1155/2013/920612

[3] Xie, Boyu et al. “Opposite changes in the expression of clathrin and caveolin-1 in normal and cancerous human prostate tissue: putative clathrin-mediated recycling of EGFR.” Histochemistry and cell biology vol. 159,6 (2023): 489-500. doi:10.1007/s00418-023-02183-8


# GSEA
Next, I created a ranked gene list from the full differential expression results and to prepare it for Gene Set Enrichment Analysis (GSEA) by mapping gene symbols to Entrez IDs and removing duplicates or missing values. It then runs gseGO() and gseKEGG() to identify enriched biological processes and pathways using the entire ranked list rather than a subset of significant genes. The goal is to capture coordinated, pathway-level shifts between tumor and normal samples.

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
sorted_allgenelist <- allgenelist[order(-allgenelist$logFC), ]
sorted_allgenelist$ENSEMBL <- rownames(sorted_allgenelist)
head(sorted_allgenelist)

#create a working data frame of ranked genelist as input for GSEA:
rank_df <- data.frame(
  logFC   = as.numeric(sorted_allgenelist$logFC),
  ENSEMBL = sorted_allgenelist$ENSEMBL,
  SYMBOL  = as.character(sorted_allgenelist$gene_name),
  stringsAsFactors = FALSE)

rank_df <- rank_df[!duplicated(rank_df$SYMBOL),]
rank_df$ENTREZID <- go$ENTREZID[match(rank_df$SYMBOL, go$SYMBOL)]


#keep only rows with an Entrez ID and numeric logFC
rank_df <- rank_df[!is.na(rank_df$ENTREZID) & !is.na(rank_df$logFC), ]
head(rank_df)
geneList <- rank_df$logFC
names(geneList) <- rank_df$ENTREZID


library(fgsea); library(qusage)
gsea_go <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  by            = "fgsea",
  verbose       = FALSE
)

head(as.data.frame(gsea_go), 10)

dotplot(gsea_go, showCategory = 10) +
  ggtitle("GSEA GO: Top 10")



gsea_kegg <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(gsea_kegg)

dotplot(gsea_kegg, showCategory = 10) +
  ggtitle("GSEA KEGG: Top 10")
```

The GSEA results provide a broader, pathway-level view of how tumor samples differ from normal tissue without relying on strict cutoffs for significance. In the GO analysis, the strongest signals are centered around mitosis and chromosome segregation, which is consistent with increased proliferation and altered chromosomal regulation as key features of the tumor transcriptome. Overall points to coordinated activation of cell cycle programs across many genes. The KEGG analysis supports a similar theme but at a more pathway-oriented level. The enriched pathways include cell cycle, pathways in cancer, Calcium signaling, and Ribosome, may reflect on the enhanced growth signaling, structural remodeling, and increased protein synthesis associated with a highly proliferative and transcriptionally rewired tumor state. In result, the GSEA summarizes that tumors are driven by coordinated activation of proliferation and cancer-related pathways, reflecting consistent biological shifts rather than random expression changes.

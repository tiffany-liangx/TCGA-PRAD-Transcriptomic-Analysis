
## 1. load in matched tumor- normal tissue sample from TCGA-PRAD

For this project, I will explore the TCGA-PRAD dataset to compare differentially expressed genes between primary tumor and solid tissue normal samples using the limma-voom workflow (using the pipline from this link: https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html#3_RNASeq_Normalization. 

The TCGA-PRAD dataset, part of The Cancer Genome Atlas, provides a comprehensive molecular profile of prostate adenocarcinoma, a common malignancy affecting approximately one in nine men. Most diagnoses occur after age 66, and the disease disproportionately impacts African American men [1]. The analysis begins by querying TCGA-PRAD RNA-seq transcriptome data processed with the STAR pipeline to obtain raw gene-level counts. Then, the metadata from the query are examined to summarize the sample composition revealing 501 primary tumor samples, 52 solid tissue normal samples, and one metastatic sample. To enable a biologically controlled comparison, the goal is to identifies patients who contribute both a primary tumor (TP) and a solid tissue normal (NT) sample based on the TCGA barcodes. Patient matching is based on the first 12 characters of each barcode to ensure that tumor samples are directly compared with their paired normal tissue from the same individual. Although it reduces the possible samples, this approach can will also reduce inter-patient variability. Using this strategy, approximately 52 matched tumor–normal pairs are identified and to be downloaded into R. Lastly, filtering is applied to the dataset to retain only protein-coding genes to reduce noise from low-confidence or noncoding features. 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='hide', cache=FALSE}
#import libraries:
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(SummarizedExperiment)
library(gplots)
library(tibble)
library(NormalyzerDE)
library(DT)
library(tibble)
library(tidyverse)
library(DESeq2) 



#project summary
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

#query all PRAD STAR-counts (gene expression quantification)
TargetSamples <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")


CaseInfo <- getResults(TargetSamples)
#as_tibble(head(CaseInfo))
summary(factor(CaseInfo$sample_type))

#Subset samples (Tumor TP vs Normal NT)
dataPrimary_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "TP")
dataNormal_Target  <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "NT")


# match barcodes based on sample type
TP_barcode <- substr(dataPrimary_Target, 1, 12)
NT_barcode <- substr(dataNormal_Target, 1, 12)
matched_barcode <- intersect(TP_barcode, NT_barcode)
length(matched_barcode) #resulted in 52 matched pairs!

#extract matching patient ID
matched_TP <- dataPrimary_Target[TP_barcode %in% matched_barcode]
matched_NT <- dataNormal_Target[NT_barcode %in% matched_barcode]

length(matched_TP)# two of the patients have 2 tumor samples
length(matched_NT)
matched_TP <- matched_TP[!duplicated(substr(matched_TP, 1, 12))] # keep exactly one TP and one NT per patient

TargetSamples_matched <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(matched_TP, matched_NT)
)


TSM<-getResults(TargetSamples_matched)
summary(factor(TSM$sample_type))

GDCdownload(TargetSamples_matched)
data <- GDCprepare(TargetSamples_matched)

#normalization:
SECoding <- data[rowData(data)$gene_type == "protein_coding", ]
counts <- assay(SECoding) #raw matrix counts

```

 

## 2. limma-voom workflow DE + FC analysis
Goal: Normalized data and generate list of differentially expressed genes with limma_pipeline wrapper function. Fold change analysis to determine top 3 upregulated & downregulated genes.

Next, it will use the wrapper function, limma_pipline() given from the reference link to normalize, filter, and generate a list of differently expressed genes using bioconductor packages like limma, edgeR, and SummarizedExperiment. Within the function, the condition variable (sample type) is extracted from the sample metadata and converted into a factor. The reference group is set to “Solid Tissue Normal,” which establishes normal tissue as the baseline when constructing the design matrix for it to be modeled statistically. Before differential expression analysis, genes with very low expression are filtered out using filterByExpr() from the edgeR package. This step removes genes with insufficient read counts, which improves statistical power and reduces noise. The filtered count data are stored in a DGEList object, enabling access to edgeR’s normalization tools. Normalization is then performed using the trimmed mean of M-values (TMM) method, which corrects for differences in library size and composition across samples. To make RNA-seq data compatible with limma’s linear modeling framework, the normalized counts are transformed using the voom method. Voom estimates the mean–variance relationship and produces log2-transformed expression values with associated precision weights. Linear models are then fitted to each gene using lmFit, and empirical Bayes moderation is applied with eBayes to stabilize variance estimates across genes. Finally, topTable extracts the top 100 differentially expressed genes ranked by p-value, summarizing the most statistically significant expression changes between tumor and normal tissues.

Function: 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

limma_pipeline = function(
    tcga_data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes

    )
  )
}

```

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#colnames(colData(SECoding))
limma_res = limma_pipeline(
  tcga_data=SECoding,
  condition_variable="sample_type",
  reference_group="Solid Tissue Normal")

#outputs:
#limma_res$topGenes
#limma_res$fit
top100genes<-limma_res$topGenes# used for Fold change and PCA analysis
dim(top100genes)
allgenes <- topTable(limma_res$fit, coef=ncol(design), number=Inf, sort.by="p") # used later for volcano plots

```

Fold change analysis:
The goal of the fold-change analysis is to further evaluate magnitude of change in differentially expressed genes between tumor vs normal samples rather than only statistical significance amongst the top 100 differentially expressed genes. An adjusted p-value (FDR) threshold of 0.05 is applied to control for multiple testing and reduce false positives, while a log2 fold-change threshold greater than 1 is used to focus on genes with at least a two-fold difference in expression between tumor and normal tissue samples. As a result, 12 genes are filtered out from the initial list, leaving 88 of the top 100 genes have a strong effect in fold change and statistical support. In result, DLX1, ZIC2, and NKX2-3 genes have the highest positive fold change indicating that they are likely to be strongly upregulated in tumor tissue. CRYAB, ID4, and AOX1 show the largest negative log fold changes indicating strong downregulation in tumors. 
As for clinical relevance, DLX1 has been shown to be upregulated in prostate cancer clinical samples and to promote tumor cell growth, migration, and colony formation through activation of the β-catenin/TCF signaling pathway [2]. In contrast, CRYAB appears to play a tumor-suppressive role [3]. As high expressions of CRYAB negatively correlated with the advanced clinicopathological characteristics and metastasis status of PCa patients [3]. 

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

#Set thresholds for FDR and fold-change
fdr_threshold <- 0.05
lfc_threshold <- 1


# filter for subset genes that pass both thresholds
sig_fc <- top100genes[top100genes$adj.P.Val < fdr_threshold &  abs(top100genes$logFC) > lfc_threshold, ]

dim(sig_fc)

sorted_fc <- sig_fc[order(-sig_fc$logFC), ]
head(sorted_fc$gene_name)[1:3]
tail(sorted_fc$gene_name)[1:3]
```



| Code| printed output |
|---|---|
|dim(sig_fc)|## [1] 88 16|
|head(sorted_fc$gene_name)[1:3]|## [1] "DLX1"   "ZIC2"   "NKX2-3"|
|tail(sorted_fc$gene_name)[1:3]|## [1] "CRYAB" "ID4"   "AOX1"|

## 3. PCA and volcano plots

Import plot_PCA wrapper function from the assignment's reference link. It produces a PCA plot given from the voom object outputted in the limma_pipline function. It defines PC1 and PC2 based on the largest variance to plot for tumor vs normal sample types. In results, there's a clear subset distinction of RNA expression between the two sample types. 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
```

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
dim(limma_res$voomObj)#all genes from 104 samples 
res_pca <-plot_PCA(limma_res$voomObj, "sample_type")

```
**Results/PCA.png**

The differential gene expression between normal and tumor prostate tissue using all analyzed genes in a volcano plot. The resulted plot demonstrates clear transcriptional differences between the two sample types. Genes that are significantly upregulated in tumors are shown in red, while downregulated genes are shown in blue, using solid normal tissue as the reference. The symmetrical shape of the plot indicates balanced up- and down-regulation, and the clustering of genes at high −log10(adjusted p-values) highlights a substantial number of statistically significant differentially expressed genes.

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

Sc <- factor(limma_res$voomObj$design[,2])
groups <- make.names(c("normal","tumor"))
levels(Sc) <- groups

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(allgenes$logFC, -log10(allgenes$adj.P.Val), main=paste(groups[1], "vs", groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(allgenes, adj.P.Val<0.05 & abs(logFC) >= 0),
     points(logFC, -log10(adj.P.Val), pch=20, col=(sign(logFC) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)


```
**Results/volcanoplot.png**


## References:

Maié, Tiago, and Martin Manolov. Analysis of Cancer Genome Atlas in R, 20 Nov. 2020, www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html#2_TCGA_data. 

[1] Wang T, Lewis B, Elaimy A, et al. Adenocarcinoma of the Prostate: Future Directions for Translational Science. In: Bott SRJ, Ng KL, editors. Prostate Cancer [Internet]. Brisbane (AU): Exon Publications; 2021 May 27. Chapter 7. Available from: https://www.ncbi.nlm.nih.gov/books/NBK571329/ doi: 10.36255/exonpublications.prostatecancer.translationalscience.2021

[2] Liang, Ming, et al. “DLX1, a Binding Protein of Beta-Catenin, Promoted the Growth and Migration of Prostate Cancer Cells.” Experimental Cell Research [United States], vol. 363, no. 1, 2018, pp. 26–32, https://doi.org/10.1016/j.yexcr.2018.01.007.

[3] Oh-Hohenhorst, Su Jung, and Tobias Lange. “Role of Metastasis-Related microRNAs in Prostate Cancer Progression and Treatment.” Cancers [Basel], vol. 13, no. 17, 2021, p. 4492, https://doi.org/10.3390/cancers13174492.


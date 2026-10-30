# TCGA- PRAD Transcriptomic Analysis
**R | RNA-seq | Differential Gene Expression | Pathway Enrichment | Gene CorrelationNetwork Analysis | WGCNA**

In this analysis , RNA-seq data from the TCGA-PRAD dataset were analyzed to identify genes and pathways that differentiate prostate tumor tissue from matched normal samples. The analysis integrates several complementary approaches, including differential expression analysis, pathway enrichment, gene set enrichment analysis (GSEA), and gene co-expression network analysis using WGCNA. This workflow moves from identifying individual differentially expressed genes, to interpreting pathway-level biological changes, and finally linking gene expression patterns to clinical traits through co-expression modules. One strength of this study design is the use of patient-matched tumor and normal samples, which helps reduce inter-patient variability and improves the biological relevance of the comparisons. However, restricting the analysis to matched samples reduces the dataset to approximately 52 pairs, which may limit statistical power in some enrichment analyses such as KEGG over-representation tests. In this context, correlation and network-based approaches, such as WGCNA, can be particularly useful for identifying coordinated transcriptional programs that may not be detected through single-gene analyses alone.

## Research Questions

- Which genes are differentially expressed between prostate tumor and normal tissue?
- Which biological pathways are associated with prostate cancer?
- Which genes emerge as network hubs, and how do they compare with known protein interactions?
- Which co-expression modules and genes are associated with tumor status and Gleason score?

## Workflow

1. **Data Acquisition and Preprocessing**
   - Retrieved TCGA-PRAD RNA-seq gene counts and sample metadata.
   - Identified approximately 52 matched tumor-normal patient pairs using TCGA barcodes.
   - Filtered the dataset to retain protein-coding genes.

2. **Differential Gene Expression Analysis**
   - Normalized count data using TMM and applied voom transformation.
   - Used `edgeR`, `limma`, and `SummarizedExperiment` for differential expression analysis.
   - Identified significantly differentially expressed genes using FDR < 0.05 and |log2FC| > 1 thresholds.
   - Visualized expression patterns using PCA, volcano plots, and hierarchical clustering.

3. **Functional Enrichment Analysis**
   - Separated differentially expressed genes into upregulated and downregulated groups.
   - Performed Gene Ontology (GO) and KEGG over-representation analyses.
   - Applied Benjamini-Hochberg correction to account for multiple testing.
   - Used Gene Set Enrichment Analysis (GSEA) to identify coordinated pathway-level expression changes.

4. **Gene Correlation Network Analysis**
   - Constructed a Pearson correlation network using a correlation threshold of r >= 0.7.
   - Used `igraph` to evaluate network connectivity and betweenness centrality.
   - Compared co-expression relationships with high-confidence protein-protein interactions from STRING (combined score >= 700).

5. **Weighted Gene Co-expression Network Analysis (WGCNA)**
   - Constructed a weighted gene co-expression network using a soft-thresholding power of 20.
   - Identified co-expression modules and calculated module eigengenes.
   - Assessed module associations with tumor status, age, Gleason score, and race.
   - Investigated candidate genes in the MEpink module based on their associations with tumor status and Gleason score.

## Methods

| Analysis Stage | Methods and Tools | Purpose |
|---|---|---|
| **Data Acquisition & Preprocessing** | TCGA-PRAD, STAR gene-level counts, TCGA barcodes, protein-coding gene filtering | Retrieved RNA-seq data and identified approximately 52 matched tumor-normal patient pairs. |
| **Differential Gene Expression** | `edgeR`, `limma`, `SummarizedExperiment`, `filterByExpr()`, TMM normalization, `voom()`, `lmFit()`, `eBayes()`, `topTable()` | Normalized expression data and identified genes with differential expression between tumor and normal tissue. |
| **Gene Filtering & Visualization** | FDR < 0.05, \|log2FC\| > 1, PCA, volcano plots, hierarchical clustering, Euclidean distance | Selected 88 candidate genes meeting the reported thresholds and visualized expression patterns. |
| **Functional Enrichment** | GO, KEGG, `org.Hs.eg.db`, `enrichGO()`, Benjamini-Hochberg correction | Identified enriched biological processes and pathways among upregulated and downregulated genes. |
| **Gene Set Enrichment Analysis** | GSEA, full ranked gene list | Evaluated coordinated pathway-level expression changes across the ranked genes. |
| **Gene Correlation Network** | Pearson correlation (r >= 0.7), `igraph`, STRING (combined score >= 700) | Examined gene co-expression relationships, network connectivity, candidate hub genes, and overlap with known protein-protein interactions. |
| **WGCNA** | Median absolute deviation (MAD) filtering, `pickSoftThreshold()`, soft-thresholding power β = 20, `blockwiseModules()` | Identified co-expression modules and evaluated their associations with tumor status, Gleason score, age, and race. |

## Key Findings

### 1. Differential Gene Expression
- Identified key tumor-associated genes, including **DLX1, ZIC2, and NKX2-3** as strongly upregulated and **CRYAB, ID4, and AOX1** as downregulated in prostate tumor tissue.
- These genes represent potential candidates for further investigation in prostate cancer biomarker discovery and disease progression.

### 2. Pathway Enrichment
- Upregulated genes were enriched in developmental, differentiation-related, and proliferative processes.
- Downregulated genes were associated with membrane organization, receptor-mediated processes, and metabolic pathways.
- GSEA highlighted coordinated changes in mitosis, chromosome segregation, and cell-cycle-related pathways, consistent with altered tumor growth programs.

### 3. Gene Correlation Network Analysis
- Identified highly connected genes, including **APOBEC3C** and **GSTP1**, as network hubs in the correlation-based network.
- Comparison with high-confidence STRING protein-protein interactions identified overlapping gene pairs associated with oxidative stress response, protein quality control, and calcium-dependent signaling.

### 4. WGCNA
- Identified the **MEpink** co-expression module as positively associated with tumor status and moderately correlated with Gleason score (*r* = 0.33).
- Highlighted **HJURP, MKI67, RRM2, NEIL3, MYBL2, and CDC25C** as candidate genes associated with tumor presence and disease severity.
- **HJURP** emerged as a candidate of interest due to its association with tumor status and Gleason score, consistent with previously reported roles in cell proliferation.

Overall, this systems-level perspective highlights potential biomarkers and regulatory hubs, such as DLX1, GSTP1, and HJURP that may have diagnostic or prognostic value if their expression correlates with disease severity or treatment response.

## Data Source & References: 

The Cancer Genome Atlas (TCGA), Prostate Adenocarcinoma (PRAD).

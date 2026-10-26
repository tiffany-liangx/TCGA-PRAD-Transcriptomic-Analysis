# Transcriptomic-Analysis-of-Prostate-Cancer-Using-TCGA-PRAD
**R | RNA-seq | Differential Gene Expression | Pathway Enrichment | Gene CorrelationNetwork Analysis | WGCNA**

The TCGA-PRAD dataset, part of The Cancer Genome Atlas, provides a comprehensive molecular profile of prostate adenocarcinoma, a common malignancy affecting approximately one in nine men. Most diagnoses occur after age 66, and the disease disproportionately affects African American men [1]. Advances in screening and molecular profiling have improved early detection and risk stratification, increasing interest in identifying biomarkers that can distinguish tumor from normal tissue and help inform diagnosis, prognosis, and treatment decisions [1]. In this analysis , RNA-seq data from the TCGA-PRAD dataset were analyzed to identify genes and pathways that differentiate prostate tumor tissue from matched normal samples. The analysis integrates several complementary approaches, including differential expression analysis, pathway enrichment, gene set enrichment analysis (GSEA), and gene co-expression network analysis using WGCNA. This workflow moves from identifying individual differentially expressed genes, to interpreting pathway-level biological changes, and finally linking gene expression patterns to clinical traits through co-expression modules. One strength of this study design is the use of patient-matched tumor and normal samples, which helps reduce inter-patient variability and improves the biological relevance of the comparisons. However, restricting the analysis to matched samples reduces the dataset to approximately 52 pairs, which may limit statistical power in some enrichment analyses, such as KEGG over-representation tests. In this context, correlation and network-based approaches, such as WGCNA, can be particularly useful for identifying coordinated transcriptional programs that may not be detected through single-gene analyses alone.

## Research Questions

- Which genes are differentially expressed between prostate tumor and normal tissue?
- Which biological pathways and functional processes are associated with tumor development?
- Which co-expression modules are associated with tumor status and Gleason score?
- Which genes emerge as candidates for further investigation based on their expression patterns and network relationships?

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
- These findings provide additional context for understanding relationships among tumor-associated genes.

### 4. WGCNA and Clinical Associations
- Identified the **MEpink** co-expression module as positively associated with tumor status and moderately correlated with Gleason score (*r* = 0.33).
- Highlighted **HJURP, MKI67, RRM2, NEIL3, MYBL2, and CDC25C** as candidate genes associated with tumor presence and disease severity.
- **HJURP** emerged as a candidate of interest due to its association with tumor status and Gleason score, consistent with previously reported roles in cell proliferation.

## Limitations and Future Directions
One limitation of this analysis relates to the WGCNA network construction. The dendrogram suggests a relatively sparse network that may be limited by the selected high soft-thresholding power (20) used to approximate scale-free topology. The analysis can be further optimized by exploring with lower soft thresholds powers such as 10, to produce a more connected network and to potentially identify more co-expression modules. Another interesting direction would be to further examine the **MEblack** module, which showed a moderate correlation with the Black/African American clinical trait (r = 0.25). Because prostate cancer disproportionately affects African American men, investigating genes within this module may help identify transcriptional patterns linked to population-specific differences in disease biology [1]. Future studies could further validate and expand these findings by testing the gene signatures in independent patient cohorts and experimentally confirming expression changes using qPCR, Western blotting, or immunohistochemistry. Functional experiments, such as CRISPR or siRNA knockdown, could also help clarify the roles of key genes like HJURP or DLX1 in prostate cancer progression. 


**Overall:** While many of the genes identified in this analysis have already been linked to prostate cancer biology, combining differential expression analysis with pathway enrichment, co-expression networks, and module–trait relationships provides a broader view of how these genes function within coordinated transcriptional programs. This systems-level perspective highlights potential biomarkers and regulatory hubs, such as DLX1, GSTP1, and HJURP, that may have diagnostic or prognostic value if their expression correlates with disease severity or treatment response.


## Repository Structure

```text
tcga-prad-transcriptomic-analysis/
├── README.md
├── scripts/
│   ├── 01_data_preprocessing.R
│   ├── 02_differential_expression.R
│   ├── 03_enrichment_GSEA.R
│   ├── 04_network_analysis.R
│   └── 05_WGCNA.R
├── results/
│   ├── figures/
│   └── tables/
└── data/
    └── README.md
```


## Data Source & References: 

The Cancer Genome Atlas (TCGA), Prostate Adenocarcinoma (PRAD).

[1] Wang T, Lewis B, Elaimy A, et al. Adenocarcinoma of the Prostate: Future Directions for Translational Science. In: Bott SRJ, Ng KL, editors. Prostate Cancer [Internet]. Brisbane (AU): Exon Publications; 2021 May 27. Chapter 7. Available from: https://www.ncbi.nlm.nih.gov/books/NBK571329/ doi: 10.36255/exonpublications.prostatecancer.translationalscience.2021

[2] Liang, Ming, et al. “DLX1, a Binding Protein of Beta-Catenin, Promoted the Growth and Migration of Prostate Cancer Cells.” Experimental Cell Research [United States], vol. 363, no. 1, 2018, pp. 26–32, https://doi.org/10.1016/j.yexcr.2018.01.007.

[3] Oh-Hohenhorst, Su Jung, and Tobias Lange. “Role of Metastasis-Related microRNAs in Prostate Cancer Progression and Treatment.” Cancers [Basel], vol. 13, no. 17, 2021, p. 4492, https://doi.org/10.3390/cancers13174492.

[4] Jiang, Z., Zhang, Y., Chen, X. et al. Inactivation of the Wnt/β-catenin signaling pathway underlies inhibitory role of microRNA-129-5p in epithelial–mesenchymal transition and angiogenesis of prostate cancer by targeting ZIC2. Cancer Cell Int 19, 271 (2019). https://doi.org/10.1186/s12935-019-0977-9

[5] Wu, L., Quan, W., Yue, G. et al. Identification of a novel six autophagy-related genes signature for the prognostic and a miRNA-related autophagy predictor for anti-PD-1 therapy responses in prostate cancer. BMC Cancer 21, 15 (2021). https://doi.org/10.1186/s12885-020-07725-0

[6] Carey, Jason P W et al. “Inhibitor of differentiation 4 (Id4) is a potential tumor suppressor in prostate cancer.” BMC cancer vol. 9 173. 7 Jun. 2009, doi:10.1186/1471-2407-9-173

[7] Li, Weiqiang et al. “Genome-wide Scan Identifies Role for AOX1 in Prostate Cancer Survival.” European urology vol. 74,6 (2018): 710-719. doi:10.1016/j.eururo.2018.06.021

[8] https://amigo.geneontology.org/amigo/term/GO:0007389

[9] da Silva, Henrique B et al. “Dissecting Major Signaling Pathways throughout the Development of Prostate Cancer.” Prostate cancer vol. 2013 (2013): 920612. doi:10.1155/2013/920612

[10] Xie, Boyu et al. “Opposite changes in the expression of clathrin and caveolin-1 in normal and cancerous human prostate tissue: putative clathrin-mediated recycling of EGFR.” Histochemistry and cell biology vol. 159,6 (2023): 489-500. doi:10.1007/s00418-023-02183-8

[11] Pang, Z., Wang, J., Xu, Y., Ji, B., Ren, M., & Ding, B. (2026). APOBEC3C Suppresses Prostate Cancer by Regulating Key Molecules Involved in Cellular Inflammation, Cell Cycle Arrest, and DNA Damage Response. Cancers, 18(1), 170. https://doi.org/10.3390/cancers18010170 

[12] Martignano, Filippo et al. “GSTP1 Methylation and Protein Expression in Prostate Cancer: Diagnostic Implications.” Disease markers vol. 2016 (2016): 4358292. doi:10.1155/2016/4358292

[13] Lai, W., Zhu, W., Xiao, C., Li, X., Wang, Y., Han, Y., Zheng, J., Li, Y., Li, M., & Wen, X. (2021). HJURP promotes proliferation in prostate cancer cells through increasing CDKN1A degradation via the GSK3β/JNK signaling pathway. Cell Death & Disease, 12(6), Article 583. https://doi.org/10.1038/s41419-021-03870-x


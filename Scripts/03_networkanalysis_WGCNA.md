# Create a correlation network for TCGA-PRAD analysis and compare to STRINGdb

Computation is reduced to expression matrix of top 100 genes. Top 100 genes are filtered based on fold change and FDR threshold resulting in 88 genes DE genes. 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

#import additional libraries:

library(WGCNA)
library(CorLevelPlot)
library(gridExtra)
library(igraph)
library(STRINGdb)

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

#convert rownames to gene_names
w_df <- sorted_fc[, c("gene_id", "gene_name")]
sig_genes_expr <- merge(w_df, expr, by = "gene_id")
rownames(sig_genes_expr)<-sig_genes_expr$gene_name
sig_genes_expr<-sig_genes_expr[,3:106]


```


Next, a gene–gene correlation network was constructed using Pearson correlation to examine relationships among the significant genes. Because the analysis already focuses on a subset of differentially expressed genes, the network is expected to be relatively dense. The main objective is to determine whether these genes form coordinated transcriptional programs and to identify highly connected hub genes within this significant gene set. Pearson correlations were calculated for each gene, and an adjacency matrix was generated using a correlation threshold of 0.7, while self-correlations were removed. In the resulting network, nodes represent genes and edges represent co-expression relationships. The edge and node tables were exported for potential visualization in Cytoscape. For exploratory analysis, the network was reconstructed in igraph, where basic network statistics were calculated for each gene, including in-degree, out-degree, and betweenness centrality. Degree measures the number of connections a gene has within the network, while betweenness reflects how important a gene is for linking different regions of the network. In result, genes were then ranked by in-degree and out-degree to identify the most connected hub genes.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}


#transform table so that it matches the input needed for the corrlation analysis 
cor.matrix<- t(sig_genes_expr)

#computes pairwise Pearson correlation of significant genes.
gene_cor <- cor(cor.matrix, method = "pearson")

#define correlation threshold and create an adjacency matrix
threshold <- 0.7
adj_matrix <- abs(gene_cor)>=threshold
diag(adj_matrix) <- 0

#prep adjacency matrix for igraph 
gene_network <- graph_from_adjacency_matrix(
  adjmatrix = adj_matrix,
  mode = "undirected",
  weighted = FALSE,
  diag = FALSE)
gene_network

#extract a list of edges and nodes from the gene network:
edge_df <- as.data.frame(get.edgelist(gene_network))
colnames(edge_df) <- c("source", "target")

head(edge_df,3)
write.table(edge_df, file = "edges.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

node_df <- data.frame(
  id = colnames(cor.matrix),
  stringsAsFactors = FALSE)

head(node_df,3)
write.table(node_df, file = "nodes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#reconstruct the network in igraph using the edge and node tables
net <- graph_from_data_frame(edge_df, directed=TRUE, vertices=node_df)
plot(net)


#calculate network statistics for each gene.
summary_net <- data.frame(
  node = V(net)$name,
  in_degree = degree(net, mode = "in"), #number of incoming edges
  out_degree = degree(net, mode = "out"), #number of outgoing edges
  betweenness = betweenness(net, directed = TRUE),
  stringsAsFactors = FALSE
)

incoming_hubs  <- summary_net[order(summary_net$in_degree, decreasing = TRUE), ]
outgoing_hubs <- summary_net[order(summary_net$out_degree, decreasing = TRUE), ]

head(incoming_hubs, 10)   # top 10 incoming hubs
head(outgoing_hubs, 10)  # top 10 outgoing hubs


```

Genes with high in-degree receive many connections and often represent convergence points of shared co-expression patterns, acting as downstream effectors or markers of coordinated biological programs. In this network, APOBEC3C emerged as the top incoming hub. Pang et al. reported that APOBEC3C is downregulated in prostate cancer and that reduced expression is associated with advanced tumor stage, higher Gleason scores, and poorer survival. APOBEC3C has been linked to an activated anti-tumor immune microenvironment, including increased CD8⁺ T cell infiltration and reduced M2 macrophages. It also furthers activates downstream signaling pathways involved in DNA damage response and inflammatory signaling, highlighting important characteristics of a hub gene receiving a lot of incoming influences of the network.  On the otherhand, genes with high out-degree connect to many other genes and may act as regulatory or signaling hubs that influence multiple downstream pathways. In this network, GSTP1 ranked as the top outgoing hub. GSTP1 functions in detoxification and interacts with signaling pathways that regulate proliferation, differentiation, and apoptosis. Filippo et al. demonstrated that GSTP1 is frequently hypermethylated in prostate cancer, leading to reduced expression compared with normal prostate tissue. Because also GSTP1 belongs to the GSTs family, a group of enzymes,  is consistent with its high out-degree in the network. 

References: 
[1] Pang, Z., Wang, J., Xu, Y., Ji, B., Ren, M., & Ding, B. (2026). APOBEC3C Suppresses Prostate Cancer by Regulating Key Molecules Involved in Cellular Inflammation, Cell Cycle Arrest, and DNA Damage Response. Cancers, 18(1), 170. https://doi.org/10.3390/cancers18010170
[2] Martignano, Filippo et al. “GSTP1 Methylation and Protein Expression in Prostate Cancer: Diagnostic Implications.” Disease markers vol. 2016 (2016): 4358292. doi:10.1155/2016/4358292





```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
library(STRINGdb)

#create a STRINGdb object for human PPI data
string_db = STRINGdb$new(version="11.5",species=9606, score_threshold=0,input_directory="")

#input gene set for STRING.
genes <- unique(c(edge_df$source, edge_df$target))
gene_df <- data.frame(gene = genes, stringsAsFactors = FALSE)

#map the gene symbols from the correlation network to STRING protein ID
mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
head(mapped)

# keep only high-confidence edges
score_cutoff <- 700

#get all pairwise interactions among the mapped STRING IDs
ppi <- string_db$get_interactions(mapped$STRING_id)
ppi <- ppi[ppi$combined_score >= score_cutoff, ]

#convert STRING protein IDs back to gene symbols
id2gene <- setNames(mapped$gene, mapped$STRING_id)

#create STRING edge table 
string_edges <- data.frame(
  source = id2gene[ppi$from],
  target = id2gene[ppi$to],
  score  = ppi$combined_score,
  stringsAsFactors = FALSE)
string_edges <- unique(string_edges[, c("source","target","score")])

edges_corr <- edge_df[, c("source","target")]
edges_corr <- unique(edges_corr)

#helper function to create an undirected key for each edge.
make_key <- function(df) {
  apply(df[, c("source","target")], 1, function(x) paste(sort(x), collapse="--"))}

#generate unique undirected edge identifiers
corr_keys   <- unique(make_key(edges_corr))
string_keys <- unique(make_key(string_edges))

#determine genes that are represented in both network
overlap_keys <- intersect(corr_keys, string_keys)
overlap_keys

#summarize the comparison of degree of overlap from both networks
comparison_summary <- data.frame(
  corr_edges = length(corr_keys),
  string_edges = length(string_keys),
  overlap_edges = length(overlap_keys),
  overlap_percent.frac_of_string = (length(overlap_keys) / length(string_keys))*100
)

comparison_summary

```
The next step was to compare the correlation-based gene network with a known protein–protein interaction (PPI) network from the STRING database to see whether strongly co-expressed gene pairs are also supported by known biological interactions. First, a STRING database object for human genes was created, and the genes from the correlation network were used as the input set. These genes were then mapped to STRING protein identifiers so that known interactions could be retrieved. To keep the analysis reliable, the STRING interactions were filtered to include only high-confidence edges with a combined score of at least 700. To compare the two networks, both the correlation edges and the STRING interactions were converted into undirected gene-pair identifiers, which allowed overlapping interactions to be identified regardless of gene order. The analysis revealed three overlapping gene pairs: GSTM4–GSTP1, CRYAB–HSPB8, and ANXA2–S100A6. In result, the correlation network contained 2010 edges, while the STRING PPI network contained 6 high-confidence interactions among the same gene set. Of these, 3 interactions overlapped, meaning that these gene pairs are not only strongly co-expressed in the dataset but also known to interact at the protein level.

At first, the small overlap was somewhat surprising, but it makes sense because correlation networks capture co-expression relationships, which can reflect shared regulation, pathway membership, or indirect functional associations. In contrast, the STRING database is more conservative and focuses on experimentally supported or predicted protein interactions, which often represent direct physical binding. After entering the proteins into the STRING database to explore the interactions further, the results showed that GSTM4–GSTP1 belongs to the glutathione S-transferase family involved in detoxification and oxidative stress responses. CRYAB–HSPB8 participates in protein quality control and cellular stress-response pathways, while ANXA2–S100A6 is associated with calcium-dependent signaling and cytoskeletal regulation. These interactions capture pathways related to cellular stress responses, metabolism, and tumor-associated signaling.





Notes copied from STRING database: 

(Reference link: https://string-db.org/cgi/network?taskId=beUvvxeMZCU6&sessionId=bppfZt18b1b0)

GSTM4	Glutathione S-transferase Mu 4; Conjugation of reduced glutathione to a wide number of exogenous and endogenous hydrophobic electrophiles. Active on 1-chloro- 2,4-dinitrobenzene. (218 aa)

GSTP1	Glutathione S-transferase P; Conjugation of reduced glutathione to a wide number of exogenous and endogenous hydrophobic electrophiles. Regulates negatively CDK5 activity via p25/p35 translocation to prevent neurodegeneration. (210 aa)

CRYAB	Alpha-crystallin B chain; May contribute to the transparency and refractive index of the lens. Has chaperone-like activity, preventing aggregation of various proteins under a wide range of stress conditions. (175 aa)

HSPB8	Heat shock protein beta-8; Displays temperature-dependent chaperone activity; Belongs to the small heat shock protein (HSP20) family. (196 aa)

ANXA2	Annexin A2; Calcium-regulated membrane-binding protein whose affinity for calcium is greatly enhanced by anionic phospholipids. It binds two calcium ions with high affinity. May be involved in heat-stress response. Inhibits PCSK9-enhanced LDLR degradation, probably reduces PCSK9 protein levels via a translational mechanism but also competes with LDLR for binding with PCSK9 ; Belongs to the annexin family. (357 aa)

S100A6	Protein S100-A6; May function as calcium sensor and modulator, contributing to cellular calcium signaling. May function by interacting with other proteins, such as TPR-containing proteins, and indirectly play a role in many physiological processes such as the reorganization of the actin cytoskeleton and in cell motility. Binds 2 calcium ions. Calcium binding is cooperative. (90 aa)


Additional references:
“STRINGdb: Stringdb (an R Interface to Http://String-Db.Org).” RDocumentation, www.rdocumentation.org/packages/STRINGdb/versions/1.8.1/topics/STRINGdb. Accessed 3 Mar. 2026. 


# WGCNA
Before diving in, I want to explicitly shout out the references and resources for writing the code for segment that I used to guide this analysis. They were incredibly helpful to better develop my understanding of WGCNA and wanted to make sure it's accredit to properly for this assignment:

- Vy-P-Tran. “Functionally_enigmatic_genes/Scripts.r at Master · Vy-P-Tran/Functionally_enigmatic_genes.” GitHub, github.com/vy-p-tran/Functionally_enigmatic_genes/blob/master/Scripts.R. Accessed 3 Mar. 2026. 

- kpatel427. “YouTubeTutorials/WGCNA.R at Main · KPATEL427/Youtubetutorials.” GitHub, github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R. Accessed 3 Mar. 2026. 

--

To prepare the data for constructing a WGCNA network and identifying related gene modules, we first need to determine the appropriate soft-thresholding power for network construction. The voom-normalized expression matrix is first converted into a data frame, and gene IDs are matched with gene symbols so that the expression matrix is easier to interpret. Quality control is then performed using goodSamplesGenes() to detect problematic genes or samples. In this analysis, all genes passed filtering, while one sample was flagged as an outlier and removed. Hierarchical clustering of the remaining samples confirmed that no additional obvious outliers were present.

Next, genes are ranked by median absolute deviation (MAD) so that the most variable genes are retained. This filtering step helps reduce noise and focuses the analysis on genes that are more likely to contribute meaningful biological signal. The expression matrix is then transposed into the format required by WGCNA, with samples as rows and genes as columns.

Finally, the pickSoftThreshold() function is used to test a range of candidate powers and evaluate how well each one approximates scale-free topology, a common property of biological networks. Two diagnostic plots are generated: one showing the scale-free topology fit and another showing mean connectivity across the tested powers. In result, soft-thresholding power of 20 was selected, indicating that this value provides the best balance between achieving a scale-free network structure and maintaining sufficient connectivity for downstream module detection.


```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# set gene names as row names for expression data frame
expr <- as.data.frame(limma_res$voomObj$E)
w_df <-cbind(gene_name= limma_res$voomObj$genes$gene_name, gene_id=limma_res$voomObj$genes$gene_id)
expr$gene_id<- rownames(expr)
all.expr <- merge(w_df, expr, by = "gene_id")
rownames(expr)<- all.expr$gene_name
expr<-expr[1:104]


#detect and exclude outlier samples: 
gsg<- goodSamplesGenes(t(expr))
summary(gsg)
table(gsg$goodGenes)#shows all gene pass
table(gsg$goodSamples) #shows all 1 sample outlier 


PRAD <- expr[, gsg$goodSamples ==TRUE]
dim(PRAD)

#hierarchical clustering 
htree <- hclust(dist(t(PRAD)),method="average")
plot(htree) # no visual outliers

#Convert data frame to numeric matrix and preserve row and column names
PRAD1 <- as.matrix(PRAD)
dim(PRAD1)
colnames(PRAD1) = colnames(PRAD)
rownames(PRAD1) = rownames(PRAD)


#filter based on variance to 10K genes. transposed. 
PRADdata <- t(PRAD1[order(apply(PRAD1,1,mad), decreasing = T)[1:10000],])
dim(PRADdata)
##head(PRADdata)

# determine a set of soft-thresholding powers to test
power <- c(c(1:10), seq( from =12, to=50, by=2))
sft<- pickSoftThreshold(PRADdata, powerVector=power, networkType= "unsigned", verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="red")

## result, power = 20
soft_power <- 20
```

blockwiseModules() is used to build the WGCNA network using the previously selected soft-thresholding power. This step constructs a gene co-expression network from the expression matrix and identifies groups of genes, called modules, that show similar expression patterns across samples. The network is based on an unsigned topological overlap matrix (TOM), meaning genes are grouped according to the strength of their co-expression relationships regardless of whether the correlation is positive or negative. The output of this step includes module eigengenes, which summarize the overall expression pattern of each module. These eigengenes serve as representative profiles for the modules and can later be correlated with clinical traits.

The dendrogram plot shows how genes cluster based on their co-expression similarity, with the colored bar underneath indicating module membership. In this network, a large portion of genes falls into the grey module, which usually represents genes that were not strongly assigned to any specific co-expression module. However, several distinct colored modules are also present, including larger blue and yellow modules and smaller clusters such as red and turquoise.

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#use WGCNA correlation function 
cor <- WGCNA::cor


# create the WGCNA network
bwnet <- blockwiseModules(PRADdata,
                          TOMType = "unsigned",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          randomSeed = 1234, minModuleSize = 30,
                          reassignThreshold = 0,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "PRADTOM",
                          verbose = 3)

# convert labels to colors names
mergedColors_PRAD <- labels2colors(bwnet$colors)
moduleLabels_PRAD <- bwnet$colors

module_eigengenes <- bwnet$MEs
MEs_PRAD <- orderMEs(moduleEigengenes(PRADdata, mergedColors_PRAD)$eigengenes)


# plot the dendrogram and the module colors underneath
par(mfrow = c(1,1))

plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors_PRAD[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cluster dendrogram for PRAD")




```



Prepare phenotype trait table using reference: 
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}


pheno<-as.data.frame(colData(SECoding))
data.barcode<- colnames(PRAD1)

pData_filtered <- pheno[pheno$barcode %in% data.barcode, ]
pData_filtered<- as.data.frame(pData_filtered[, c("sample_type", "age_at_index","gleason_score", "race")])


# create traits file - binarize categorical variables
traits <- pData_filtered %>%
  mutate(
    tumor = ifelse(grepl("Tumor", sample_type), 1, 0),
    age = as.numeric(age_at_index),
    gleason = as.numeric(gleason_score)
  ) %>%
  select(tumor, age, gleason)

table(pData_filtered$race)
pData_filtered$race <- factor(pData_filtered$race, levels = c("Unknown", "black or african american", "white"))

race.binarized <- binarizeCategoricalColumns(
  pData_filtered$race,
  includePairwise = FALSE,
  includeLevelVsAll = TRUE,
  minCount = 1
)
colnames(race.binarized)

#finalized table: 
traits<- cbind(traits,race.binarized)
head(traits,3)

```



Identify which modules are associated with tumor status, gleason score, age, and race
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
nGenes = ncol(PRADdata)
nSamples = nrow(PRADdata)

model.trait.cor<- cor(module_eigengenes, traits, use='p')
model.trait.cor.pvals<- corPvalueStudent(model.trait.cor, nSamples)



# visualize module-trait relationship as a heatmap
heatmap.data <- cbind(MEs_PRAD, traits)

CorLevelPlot(
  heatmap.data,
  x = colnames(traits),
  y = colnames(MEs_PRAD),
  col = c("blue", "skyblue", "white", "pink", "red")
)
```


The heatmap shows the correlation between WGCNA module eigengenes and clinical traits including tumor status, age, Gleason score, and race. Several modules show strong associations with tumor status, particularly MEgreen, MEgrey, and MEturquoise to be negatively correlated with tumor samples. In contrast, MEpink and MEbrown show strong positive correlations with tumor presence, suggesting these modules may be associated with tumor-related gene expression programs. MEpink also shows a moderate positive correlation with Gleason score (0.33), indicating potential involvement in tumor aggressiveness. Age and race show relatively weak correlations with most modules, suggesting limited influence on the detected co-expression patterns.

Determine how strongly each gene belongs to each module by calculating module membership for each gene:
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

modNames_PRAD <- substring(colnames(MEs_PRAD), 3)

geneModuleMembership_PRAD <- as.data.frame(cor(PRADdata, MEs_PRAD, use = "p"))
MMPvalue_PRAD <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_PRAD), nSamples))

colnames(geneModuleMembership_PRAD) <- paste0("MM_", modNames_PRAD)
colnames(MMPvalue_PRAD) <- paste0("p.MM_", modNames_PRAD)

#Gene info table
geneInfoPRAD <- data.frame(
  Genename = colnames(PRADdata),
  moduleColor = mergedColors_PRAD,
  geneModuleMembership_PRAD,
  MMPvalue_PRAD
)
geneInfoPRAD <- geneInfoPRAD[order(geneInfoPRAD$moduleColor), ]

head(geneInfoPRAD, 2)
```

Examine the genes in MEpink to find the overlap between tumor to which gleason score:
```{r}

pink_genes <- colnames(PRADdata)[mergedColors_PRAD == "pink"]
head(pink_genes)
length(pink_genes)

#determine corelation of each gene with tumor and gleason score 
cor_tumor <- as.data.frame(cor(PRADdata, traits$tumor, use = "p"))
cor_gleason <- as.data.frame(cor(PRADdata, traits$gleason, use = "p"))

colnames(cor_tumor) <- "cor_tumor"
colnames(cor_gleason) <- "cor_gleason"

pink_gene_cor_traits <- data.frame(
  gene = colnames(PRADdata),
  moduleColor = mergedColors_PRAD,
  cor_tumor,
  cor_gleason
)


pink_geneInfo <- pink_gene_cor_traits[pink_gene_cor_traits$moduleColor == "pink", ]

head(pink_geneInfo)

#rank pink-module genes by combining tumor + Gleason association
top_genes <- pink_geneInfo[order(-(pink_geneInfo$cor_tumor + pink_geneInfo$cor_gleason)), ][1:10, ]
top_genes

#scatter plot of gene significance for tumor vs Gleason
ggplot(pink_geneInfo, aes(x = cor_tumor, y = cor_gleason)) +
  geom_point(color = "deeppink3", size = 3) +
  geom_text(data = top_genes, aes(label = gene), vjust = -0.6, size = 3) +
  theme_minimal() +
  labs(
    title = "Top MEpink genes associated with tumor and Gleason",
    x = "Gene significance for tumor",
    y = "Gene significance for Gleason"
  )

```

I was mainly interested in identifying which genes in the MEpink module were most strongly associated with tumor status and Gleason score. To do this, the correlation between each gene’s expression and these clinical traits was calculated and combined into a table containing the gene name, module color, and the two correlation values. The genes were then ranked based on their combined association with tumor presence and disease severity. A scatter plot was generated to visualize the relationship between gene significance for tumor (x-axis) and gene significance for Gleason score (y-axis), with the top genes labeled.
One of the top genes identified was HJURP (Holliday Junction Recognition Protein). Previous work by Lai et al. reported that HJURP is significantly upregulated in prostate cancer and that higher mRNA and protein levels are associated with poorer patient outcomes, serving as an independent prognostic marker for survival. Functionally, HJURP promotes prostate cancer cell proliferation by enhancing the ubiquitination and degradation of the cell cycle inhibitor CDKN1A (p21) through the GSK3β/JNK signaling pathway. This supports the observation that HJURP appears among the top MEpink genes, suggesting it may play an important role in tumor progression 

Reference:
Lai, W., Zhu, W., Xiao, C., Li, X., Wang, Y., Han, Y., Zheng, J., Li, Y., Li, M., & Wen, X. (2021). HJURP promotes proliferation in prostate cancer cells through increasing CDKN1A degradation via the GSK3β/JNK signaling pathway. Cell Death & Disease, 12(6), Article 583. https://doi.org/10.1038/s41419-021-03870-x


Confirm whether if the pink modules become more active in higher-grade tumors by determining which gleason score shows highest MEpink activity:
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

pink_gleason <- data.frame(
  MEpink = MEs_PRAD$MEpink,
  gleason = traits$gleason
)

head(pink_gleason)

pink_by_gleason <- aggregate(MEpink ~ gleason, data = pink_gleason, FUN = mean)

pink_by_gleason <- pink_by_gleason[order(-pink_by_gleason$MEpink), ]
pink_by_gleason 

```

Results does confirm that, and indicates that pink modules become progressively more active as tumor gets more aggressive.



Export pink module network to cytoscape: 
```{r}

dataExpr_mod <- PRADdata[, pink_genes]

adjMat <- adjacency(dataExpr_mod, power = soft_power, type = "signed")
diag(adjMat) <- 0

cyt <- exportNetworkToCytoscape(
  adjMat,
  edgeFile = paste0("Cytoscape_edges_pink_genes.txt"),
  nodeFile = paste0("Cytoscape_nodes_pink_genes.txt"),
  weighted = TRUE,
  threshold = 0.05,
  nodeNames = colnames(dataExpr_mod),
  nodeAttr = mergedColors_PRAD["pink"],
  includeColNames = TRUE
)


head(cyt$edgeData)
head(cyt$nodeData)


#create an igraph object from cytoscape edge and node tables
net_2nd <- graph_from_data_frame(
  d = cyt$edgeData,
  directed = FALSE,
  vertices = cyt$nodeData
)
plot(net_2nd)
```



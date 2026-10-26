

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

#Load in functions: 
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
  allGenes= topTable(fit, coef=ncol(design), number=Inf, sort.by="p")
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes, # the 100 most differentially expressed genes
      allGenes=allGenes
    )
  )
}



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

# ------------------------------------------------------------------
#Goal: create a hierarchical clustering of TCGA-PRAD dataset.
##Top 100 genes are filtered based on fold change and FDR threshold resulting in 88 genes DE genes for hierarchical clustering: 

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

#set up data input for hclust. clusters will be labeled with gene names. 
w_df <- sorted_fc[, c("gene_id", "gene_name")]
sig_genes_expr <- merge(w_df, expr, by = "gene_id")
rownames(sig_genes_expr)<-sig_genes_expr$gene_name
sig_genes_expr<-sig_genes_expr[,3:106]


h<-hclust(dist(sig_genes_expr))
plot(h)
#---------------------------------------------------
##The remaining genes are then divided based on fold change, where a positive logFC indicates upregulation and a negative logFC represents downregulation.
#Set thresholds for FDR and fold-change
fdr_threshold <- 0.05
lfc_threshold <- 1

# filter for subset genes that pass both thresholds
sig_fc <- top100genes[top100genes$adj.P.Val < fdr_threshold &  abs(top100genes$logFC) > lfc_threshold, ]
sorted_fc <- sig_fc[order(-sig_fc$logFC), ]
head(sorted_fc)

pos_logfc<-sorted_fc[sorted_fc$logFC>0, ] #upregulation
neg_logfc<-sorted_fc[sorted_fc$logFC<0, ] #down-regulation 

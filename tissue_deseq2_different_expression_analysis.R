library(DESeq2)
library(tidyr)
library(dplyr)
ps_rawCounts <- read.csv("H:/prna_transcript_analysis/hcc_tissue/ps_output.txt",sep = '\t',skip = 1)
ps_rawCounts$Geneid = gsub('^','ps_',ps_rawCounts$Geneid)

pas_rawCounts <- read.csv("H:/prna_transcript_analysis/hcc_tissue/pas_output.txt",sep = '\t',skip = 1)
pas_rawCounts$Geneid = gsub('^','pas_',pas_rawCounts$Geneid)

prna = rbind(ps_rawCounts,pas_rawCounts)
prna = prna %>% column_to_rownames('Geneid')
rawCounts = prna

rccm = colnames(rawCounts)
colnames(rawCounts) = gsub(".sort.bam", "", rccm)
head(rawCounts)

sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])


sampleData <- read.csv("H:/prna_transcript_analysis/hcc_tissue/SraRunTable.txt",sep = ',',row.names = 1)


keep <- c("Tissue", "PATIENT_ID")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "PATIENT_ID")
sampleData$PATIENT_ID <- factor(sampleData$PATIENT_ID)
sampleData = sampleData[which(rownames(sampleData) %in%  colnames(rawCounts)),]

all(colnames(rawCounts) == rownames(sampleData))

sampleData$tissueType = gsub("\\\\.*", "",sampleData$tissueType)
#sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$binary = factor(rep(c("normal", "cancer"), times = c(16, 109)))
#sampleData$tissueType <- factor(sampleData$tissueType, levels=c("Normal tissue","Fibrosis-low tissue","Fibrosis-high tissue","Cirrhosis tissue","Dysplastic nodule-low tissue","Dysplastic nodule-high tissue","Gastric cancer"))
library(DESeq2)
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ binary)
#dim(deseq2Data)
#dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
#deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
#deseq2Data <- DESeq(deseq2Data)
#deseq2Results <- results(deseq2Data, contrast=c("tissueType", "Normal tissue", "Fibrosis-high tissue"))
#result = deseq2Results[which(deseq2Results$log2FoldChange > 2 & deseq2Results$padj < 0.05),]

#plotMA(deseq2Results)
#res <- results(deseq2Data,alpha = 0.05)
#res <- res[!is.na(res$padj), ]
#res <- res[res$padj <= 0.05, ]
#table(res$padj <= 0.05)
#res <- res[order(res$padj), ]
#resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
library(EnhancedVolcano)
dds <- DESeq(deseq2Data, betaPrior=FALSE)
res1 <- results(dds)
res1
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05)



up = as.data.frame(res1) %>% filter(log2FoldChange > 1 & pvalue < 0.05) 
down = as.data.frame(res1) %>% filter(log2FoldChange < -1 & pvalue < 0.05) 

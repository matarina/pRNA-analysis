rm(list = ls())
##from epd promoter id to ensemble gene id 
epd_id = read.csv("H:/prna_transcript_analysis/epd_promoter/4_cell_line/hepg2_1.count",sep = '\t')
ensembl_id = read.csv("H:/prna_transcript_analysis/epd_promoter/epd_annotaion/promoter_ensembl.txt",sep = '\t',header = FALSE)
colnames(ensembl_id) = c("Geneid","ensembl")
G_list = as.data.frame(G_list[,-2])
colnames(G_list)[1] = c("ensembl")
house_gene = merge(ensembl_id,G_list,by="ensembl")
house_gene = as.data.frame(house_gene[,-1])
house_gene = house_gene[!duplicated(house_gene), ]
load("H:/chrome downlaod/Housekeeping_GenesHuman.RData")
###ensembl biomart geneid converter
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- Housekeeping_Genes$Ensembl

G_list <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_gene_id","ensembl_transcript_id"),values=genes,mart= mart)



####

library(clusterProfiler)
library(org.Hs.eg.db)
#go富集分析流程
go <- enrichGO(merge$ensemble, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05,  qvalueCutoff = 0.2,keyType = 'ENSEMBL')
dotplot(go)
## 查看基因id
keytypes(org.Hs.eg.db) 

####g:profiler2
library(gprofiler2)
gost =gost(merge$ensemble,organism = "hsapiens",user_threshold = 0.05)
gostplot(gost)

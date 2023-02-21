###source code for go and gsea enrichment analysis



library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(ReactomePA)
library(tidyverse)
library(data.table)

library(enrichplot)
library(ggridges)

f <- read.table("41420_2020_328_MOESM5_ESM.csv",header=T,sep=',')
f = filter(f,logfc_RNA..OVX_vs_sham.>2)
x = f[,4]
x
#hg<-bitr(x,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Mm.eg.db")###gene id transfer to entrezid ,24.87% of input gene IDs are fail to map...
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
# query biomart
results <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),##try "hgnc_symbol",
                 filters ="ensembl_gene_id", values = x,
                 mart = mart)
results

##considering the low map rate, try ensemble trans from symbol id to entrezid,12.6% failed
#keytypes(org.Mmu.eg.db)
results$entrezgene_id = as.character(results$entrezgene_id)
keytypes()
go <- enrichGO(results$entrezgene_id,org.Mm.eg.db,keyType = "ENTREZID", ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2)
write.csv(go,file="H:/go.csv") ##save go results
barplot(go, split="ONTOLOGY",cex.main=0.05)+ facet_grid(ONTOLOGY~.,scale="free") ## bar plot visualization






### reshape data to gsea analysis
colnames(results)[1] = "geneID"
subset_f = subset(f,select = c("geneID","logfc_RNA..OVX_vs_sham."))
subset_f = distinct(subset_f, geneID, .keep_all= TRUE)
merge = merge(results,subset_f,by="geneID")
merge = na.omit(merge)
merge = distinct(merge,entrezgene_id,.keep_all= TRUE)
merge$entrezgene_id = as.character(merge$entrezgene_id)

subset_merge = subset(merge,select=c(entrezgene_id,logfc_RNA..OVX_vs_sham.))
subset_merge_sorted = subset_merge[order(subset_merge$logfc_RNA..OVX_vs_sham.,decreasing = TRUE),]
geneList = subset_merge_sorted[,2]
names(geneList) = as.character(subset_merge_sorted[,1])
geneList



#GSEA analysis——GO
Go_gseresult <- gseGO(geneList, 'org.Mm.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA analysis——KEGG
KEGG_gseresult <- gseKEGG(geneList, organism = "mmu", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA analysis——Reactome
Go_Reactomeresult <- gsePathway(geneList,organism ="mouse", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#save gsea results
write.table (Go_gseresult, file ="H:/Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="H:/KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="H:/Go_Reactomeresult.csv", sep =",", row.names =TRUE)
###ridge plot visualization for top 10 terms
ridgeplot(Go_gseresult,10)
ridgeplot(KEGG_gseresult,10)
ridgeplot(Go_Reactomeresult,10)

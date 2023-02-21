rm(list=ls())
library(tidyverse)
library(tibble)
clinical  = read.csv('H:/prna_transcript_analysis/plasma/plasma_clinical_metadata.txt',sep = ',')
ps_output = read.csv('H:/prna_transcript_analysis/plasma/ps_output.txt',sep = '\t',skip = 1)
ps_output$Geneid <- gsub('^', 'ps_', ps_output$Geneid)
pas_output = read.csv('H:/prna_transcript_analysis/plasma/pas_output.txt',sep = '\t',skip = 1)
pas_output$Geneid <- gsub('^', 'pas_', pas_output$Geneid)

##map gene id and gene type 
map = read.csv("H:/prna_transcript_analysis/plasma/map_id_type.txt",sep = '\t')
map = map %>% distinct()
map3 = map %>% filter(X.1 == "protein_coding")
#map4 = map %>% filter(X.1 == "lincRNA")

mrna_output = read.csv('H:/prna_transcript_analysis/plasma/genecode_output.txt',sep = '\t',skip = 1)
mrna_output = mrna_output[which(mrna_output$Geneid %in% map3$X),]


expr = rbind(ps_output,pas_output,mrna_output)

rownames(expr) = expr$Geneid
expr = expr[,-1]

cancer = expr[,1:35]
normal = expr[,36:65]
expr_cancer = data.frame(edgeR::cpm(cancer))
expr_normal = data.frame(edgeR::cpm(normal))

filter <- apply(expr_cancer, 1, function(x) length(x[x>15])>=12)
filter_cancer = expr_cancer[filter,]
length(grep('pas_',rownames(filter_cancer)))
#ps 1500 pas 2166 mrna mRNA 9258
filter <- apply(expr_normal, 1, function(x) length(x[x>15])>=10)
filter_normal = expr_normal[filter,]
length(grep('ENSG',rownames(filter_normal)))
#ps 1730 pas 3150 mRNA 8733

## get coefficient of variation
keep_expr_cancer = c(rownames(filter_cancer)[grep('ps_|pas',rownames(filter_cancer))])
keep_expr_normal = c(rownames(filter_normal)[grep('ps_|pas',rownames(filter_normal))])
expr_prna = unique(keep_expr_cancer,keep_expr_normal)
t_expr_cancer = as.data.frame(t(expr_cancer[expr_prna,]))
cv_cancer = sapply(t_expr_cancer, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)

t_expr_normal = as.data.frame(t(expr_normal[expr_prna,]))
cv_normal = sapply(t_expr_normal, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)

ps_cv_cancer = data.frame(cv_cancer[grep('ps',names(cv_cancer))])
colnames(ps_cv_cancer) = 'cv_cancer'
ps_cv_normal = data.frame(cv_normal[grep('ps',names(cv_normal))])
colnames(ps_cv_normal) = 'cv_normal'
ps_cv = merge(ps_cv_cancer,ps_cv_normal,by="row.names")
library(ggpubr)
ggplot(ps_cv, aes(x = cv_normal, y = cv_cancer)) +geom_smooth(method = "lm", color="#DE616A")+ theme_classic()+
  geom_point()+stat_cor(method="pearson",size = 8)+theme(text = element_text(size = 20)) 

pas_cv_cancer = data.frame(cv_cancer[grep('pas',names(cv_cancer))])
colnames(pas_cv_cancer) = 'cv_cancer'
pas_cv_normal = data.frame(cv_normal[grep('pas',names(cv_normal))])
colnames(pas_cv_normal) = 'cv_normal'
pas_cv = merge(pas_cv_cancer,pas_cv_normal,by="row.names")
ggplot(pas_cv, aes(x = cv_normal, y = cv_cancer)) +geom_smooth(method = "lm", color="#DE616A")+ theme_classic()+
  geom_point()+stat_cor(method="pearson",size = 8)+theme(text = element_text(size = 20)) 


##heatmap of biomarkers
expr_prna = column_to_rownames(merge(expr_cancer,expr_normal,by = 'row.names'),var = 'Row.names')
biomarker_expr =expr_prna[rfeImp$feature,]

library("heatmaply")
library(pheatmap)
library(ComplexHeatmap)
library("tinyarray")
library(circlize)

colnames(biomarker_expr) = c(rep(c('cancer','normal'),times=c(35,30)))
biomarker_expr = as.matrix(biomarker_expr)
Group = factor(rep(c("cancer","normal"),times = c(35,30)))
Group = factor(Group,levels = c("cancer","normal"))

#heatmaply(scale(biomarker_expr),Colv = FALSE,Rowv = FALSE)
#draw_heatmap(biomarker_expr,Group,cluster_cols = F)
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#2E5F6B", "white", "#C3573C"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2E5F6B", "#C3573C")),
                       labels = c("Normal","Cancer"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Heatmap(t(scale(t(biomarker_expr))),name = " ",
        col = col_fun,
        top_annotation = top_annotation,row_dend_side = 'right', row_names_side = 'left',
        column_split = Group,
        show_heatmap_legend = T,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL)









## here deliver expressed prna to volcano plot
expressed_prna = c(rownames(filter_cancer),rownames(filter_normal))
uniq_expressed_prna = unique(expressed_prna)



library(edgeR)
cpm_expr = data.frame(edgeR::cpm(expr[,2:ncol(expr)]))
rownames(cpm_expr) = expr$Geneid

keep =  c("SRR10822565","SRR10822566","SRR10822568","SRR10822569","SRR10822570","SRR10822571","SRR10822572","SRR10822573","SRR10822574","SRR10822579","SRR10822580","SRR10822591")
keep = gsub('$',".sort.deduplicated.bam",keep)



select_cpm_expr = cpm_expr %>% select(keep)
table(rowMeans(select_cpm_expr[grep('^pas',rownames(select_cpm_expr)),])>10)
table(rowMeans(select_cpm_expr[grep('^ps',rownames(select_cpm_expr)),])>10)
table(rowMeans(select_cpm_expr[which(rownames(select_cpm_expr) %in% map3$X),])>10)
table(rowMeans(select_cpm_expr[which(rownames(select_cpm_expr) %in% map4$X),])>10)


#pas 2319  ps 1608  protein 8261 linc 301


expressed_mrna = mrna_output2[which(rownames(mrna_output2) %in% map2$X),]

expressed_mrna$mean = rowMeans(expressed_mrna)
table(expressed_mrna$mean >= 5)



expressed_lncrna= mrna_output2[which(rownames(mrna_output2) %in% map3$X),]

expressed_lncrna$mean = rowMeans(expressed_lncrna)
table(expressed_lncrna$mean >= 5)

#pas 3036 ps 1851 protein 11053 lincrna 404





#filter data by alignment rate
align = read.csv('alignrate.txt',header = FALSE)
a = seq(1,195,3)
b = seq(3,195,3)
align2 = data.frame(a1 = align[a,], a2 = align[b,])
align2 = separate(align2,a2,into = c('rate','dump'),sep = '%')
align2$rate = as.numeric(align2$rate)
align2 = filter(align2,rate > 70)
keep_align = align2$a1
keep_align




######plasma different expressed prna between normal and hcc patient volvano plot
library(DESeq2)
library(tidyr)
library(dplyr)
library(tibble)
ps_rawCounts = read.csv('H:/prna_transcript_analysis/plasma/ps_output.txt',sep = '\t',skip = 1)
ps_rawCounts$Geneid = gsub('^','ps_',ps_rawCounts$Geneid)

pas_rawCounts = read.csv('H:/prna_transcript_analysis/plasma/pas_output.txt',sep = '\t',skip = 1)
pas_rawCounts$Geneid = gsub('^','pas_',pas_rawCounts$Geneid)
prna = rbind(ps_rawCounts,pas_rawCounts)
prna = prna %>% column_to_rownames('Geneid')


rawCounts = prna
rccm = colnames(rawCounts)
colnames(rawCounts) = gsub(".sort.deduplicated.bam", "", rccm)
head(rawCounts)


sampleData <- read.csv("H:/prna_transcript_analysis/plasma/plasma_clinical_metadata.txt",sep = ',')
sampleData = sampleData %>% select('Run','source_name') %>% column_to_rownames('Run')

sampleData$source_name <- factor(sampleData$source_name)
head(sampleData)
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

rawCounts = rawCounts[uniq_expressed_prna,]
sampleData$binary = factor(rep(c( "cancer","normal"), times = c(35,30)))


library(DESeq2)
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ binary)

library(EnhancedVolcano)
dds <- DESeq(deseq2Data, betaPrior=FALSE)
res1 <- results(dds)
res1
EnhancedVolcano(res1,
                lab =NA, legendPosition = 'right',
                x = 'log2FoldChange',
                y = 'padj',parseLabels = FALSE,caption = 'log2FC cutoff: |2|; p-value cutoff: 0.05',
                pCutoff = 0.05)

#The default cut-off for log2FC is >|2|; 
#select up/down regulated pRNA
res2 = as.data.frame(res1) %>% filter(log2FoldChange > 1 & pvalue < 0.05) 
res3 = as.data.frame(res1) %>% filter(log2FoldChange < -1 & pvalue < 0.05) 
different_pRNA_index = c(rownames(res2),rownames(res3))
saveRDS(different_pRNA_index,'H:/prna_transcript_analysis/plasma/different_pRNA_index.rds')




table(rownames(res3) %in% rownames(up))
table(rownames(res2) %in% rownames(down))



library(eulerr)

vd <- euler(c(plasma_down = 461, tissue_down = 621, "plasma_down&tissue_down" = 17))
plot(vd,
     fills = list(fill = c("#fbb4ae", "#ccebc5"), alpha = 0.6),#"#ccebc5"
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE)




####################
##################plasma sample pca plot 

rm(list=ls())
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(ggthemes)
library(tidyr)
library(dplyr)
library(tibble)

data1= read.csv('H:/prna_transcript_analysis/plasma/pas_output.txt',sep = '\t',skip = 1)
data2= read.csv('H:/prna_transcript_analysis/plasma/ps_output.txt',sep = '\t',skip = 1)
data1$Geneid = gsub('^','pas_',data1$Geneid)
data2$Geneid = gsub('^','ps_',data2$Geneid)
data = rbind(data1,data2)
data = data %>% column_to_rownames('Geneid')

coldata = data.frame(sample = colnames(data),condition = rep(c('cancer','normal'),times=c(35,28)))

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
#if error in antisense
#vsd = varianceStabilizingTransformation(dds,blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
pca = plotPCA(vsd, intgroup=c("condition"),returnData = TRUE)
ggplot(pca, aes(x=PC1, y=PC2,color=condition))+ geom_point(size=2.5)+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                      panel.grid.minor=element_blank(),
                                                                                      axis.line= element_line(colour = "black"),
                                                                                      legend.title = element_text(size=15),
                                                                                      legend.text = element_text(size=15),
                                                                                      axis.title.x = element_text(size = 15),
                                                                                      axis.title.y = element_text(size = 16),
                                                                                      plot.title = element_text(hjust = 0.4))+
  labs(title = "pRNA PCA")
#+geom_text_repel (aes(PC1,PC2, label=group))  ###add text lable on every dot
###669 473



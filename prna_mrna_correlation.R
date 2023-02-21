rm(list=ls())
setwd("H:/")
SRR8131644prna_expression = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR8131644antisense.count",sep = '\t',skip = 1)
SRR8131645prna_expression = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR8131645antisense.count",sep = '\t',skip = 1)

SRR8131644mrna_expression = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR8131644mrna.count",sep = '\t',skip = 1)
SRR8131645mrna_expression = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR8131645mrna.count",sep = '\t',skip = 1)
promoter_ensembl = read.csv("H:/prna_transcript_analysis/epd_promoter/epd_annotaion/promoter_ensembl.txt",header = FALSE,sep = '\t')
colnames(promoter_ensembl) = c("Geneid","ensembl")

SRR8131644prna_expression$average = (SRR8131644prna_expression$SRR8131644.sort.bam+SRR8131645prna_expression$SRR8131645.sort.bam)/2
SRR8131644mrna_expression$mean = (SRR8131644mrna_expression$SRR8131644.sort.bam+SRR8131645mrna_expression$SRR8131645.sort.bam)/2
prna_transid  =  merge(SRR8131644prna_expression,promoter_ensembl,by='Geneid')                                          
colnames(SRR8131644mrna_expression)[1] = "ensembl"
library(tidyr)
SRR8131644mrna_expression = separate(SRR8131644mrna_expression,ensembl,into = c("ensembl","version"),sep = "\\.")
prna_mrna = merge(prna_transid,SRR8131644mrna_expression,by='ensembl')
library(ggplot2)
library(ggpubr)
prna_mrna = prna_mrna[,c("average","mean")]
colnames(prna_mrna) = c("antisense_pRNA_count","mRNA_count")
 ggplot(prna_mrna, aes(x=antisense_pRNA_count, y=mRNA_count)) +geom_point(size=2, shape=19,color = "#5AC3DE")+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,size = 9)+theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 20), axis.title.y = element_text(angle = 90, vjust = 0.5), 
        axis.text = element_text(size = 15))+xlim(NA, 1000)+ylim(NA,1e+05)
ggsave("plot.png",p, path = "H:/",width = 6,height = 6,)




colnames(prna_mrna)

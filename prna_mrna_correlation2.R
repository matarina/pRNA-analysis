rm(list=ls())
setwd("H:/")
library(tidyverse)
ps = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/ps_output.txt",sep = '\t')
ps$FEATURE_ID = gsub('^','ps_',ps$FEATURE_ID)
pas = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/pas_output.txt",sep = '\t')
pas$FEATURE_ID = gsub('^','pas_',pas$FEATURE_ID)
mrna = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/mrna_output.txt",sep = '\t')
id = read.csv('H:/prna_transcript_analysis/bp_target_decay_cdf/promoter_ensembl.txt',sep = '\t',header = FALSE)
colnames(id) = c('epd','ensembl')
total = rbind(ps,pas,mrna)[,1:3]

total = total %>% filter(!grepl('ps_ERCC-\\|pas_ERCC-',FEATURE_ID)) %>% column_to_rownames(var = 'FEATURE_ID')
total2 = data.frame(edgeR::cpm(total))
total2$mean = rowMeans(total2)
total2 = filter(total2,mean>=0.2)
eps = total2[grep('ps_',rownames(total2)),]
epas = total2[grep('pas_',rownames(total2)),]
emrna = total2[grep('ENSG',rownames(total2)),]
eps  = eps %>% rownames_to_column() 
eps$rowname = gsub('ps_','',eps$rowname)
colnames(eps)[1] = 'epd'
eps = merge(eps,id,by = 'epd')

epas  = epas %>% rownames_to_column() 
epas$rowname = gsub('pas_','',epas$rowname)
colnames(epas)[1] = 'epd'
epas = merge(epas,id,by = 'epd')

emrna  = emrna %>% rownames_to_column() 
colnames(emrna)[1] = 'ensembl'
emrna$ensembl = gsub('\\..*$','',emrna$ensembl)
emrna2 = merge(emrna,id,by ='ensembl')
colnames(epas)[4] = 'pas_mean'
colnames(eps)[4] = 'ps_mean'
colnames(emrna2)[4] = 'mrna_mean'

b = merge(emrna2,eps,by = 'epd')
b$group = 'others'
library(biomaRt)
# define biomart object

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
load('H:/prna_transcript_analysis/bp_target_decay_cdf/Housekeeping_GenesHuman.RData')
value = Housekeeping_Genes$Ensembl
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                 filters = "ensembl_transcript_id", values = value,
                 mart = mart)
results
table(b$ensembl.x %in% results$ensembl_gene_id)

for ( i in 1:nrow(b)) {
  if (b$ensembl.x[i] %in% results$ensembl_gene_id) {  
    b$group[i] = 'house_keeping'}
}

library(ggplot2)
library(ggpubr)
ggplot(b, aes(x=mrna_mean, y=ps_mean,color = group)) + geom_point()+ theme_classic()+
  labs(x = "mrna")+labs(y = "ps")+theme(text = element_text(size = 20)) +
  scale_x_continuous(limits = c(NA, 2000))+scale_y_continuous(limits = c(NA,75))



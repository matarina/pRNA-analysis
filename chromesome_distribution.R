library(dplyr)
library(tidyr)
library(ggplot2)
### statistics of mrna chromesome distribution
mrna = read.csv("H:/ENCFF010XLY.count",header = TRUE,skip = 1,sep = '\t')
mrna = separate(mrna,Chr,into = c("chr","dump"),sep = ';')
mrna = mrna[,c(1,2)]
mrna_expr = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/latest_annotation/mrna_expr.txt",sep = '\t')
colnames(mrna_expr)[1] = "Geneid"
mrna = inner_join(mrna,mrna_expr,by='Geneid')

mrna$mean_mrna = rowMeans(subset(mrna, select = c(3:65)))
mrna_mask = filter(mrna,mean_mrna>0.5)
table_mrna = table(mrna_mask$chr)
table_mrna = as.data.frame(table_mrna)
colnames(table_mrna)[2] = "mrna"

### statistics of psrna chromesome distribution
psrna = read.csv("H:/ps_ENCFF010XLY.count",header = TRUE,skip = 1,sep = '\t')
psrna = psrna[,c(1,2)]
psrna_expr = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/latest_annotation/ps_expr.tsv",sep = '\t')
colnames(psrna_expr)[1] = "Geneid"
psrna = inner_join(psrna,psrna_expr,by='Geneid')

psrna$mean_psrna = rowMeans(subset(psrna, select = c(3:65)))
psrna_mask = filter(psrna,mean_psrna>0.5)
table_psrna = table(psrna_mask$Chr)
table_psrna = as.data.frame(table_psrna)
colnames(table_psrna)[2] = "psrna"


### statistics of pasrna chromesome distribution
pasrna = read.csv("H:/pas_ENCFF010XLY.count",header = TRUE,skip = 1,sep = '\t')
pasrna = pasrna[,c(1,2)]
pasrna_expr = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/latest_annotation/pas_expr.tsv",sep = '\t')
colnames(pasrna_expr)[1] = "Geneid"
pasrna = inner_join(pasrna,pasrna_expr,by='Geneid')

pasrna$mean_pasrna = rowMeans(subset(pasrna, select = c(3:65)))
pasrna_mask = filter(pasrna,mean_pasrna>0.5)
table_pasrna = table(pasrna_mask$Chr)
table_pasrna = as.data.frame(table_pasrna)
colnames(table_pasrna)[2] = "pasrna"



library(tidyverse)

#put all data frames into list
df_list <- list(table_mrna,table_psrna,table_pasrna)

#merge all data frames in list

chromesome_distribution = df_list %>% reduce(full_join, by='Var1')

library(ggplot2)
library(ggpubr)
### dot plot
chrome_ps_pas = ggplot(chromesome_distribution, aes(x=psrna, y=pasrna)) + geom_point(color="#44D4EB")+ theme_classic()+geom_smooth(method = "lm", color="#DE616A")+
  labs(x = "psRNA_count")+labs(y = "pasRNA_count")+stat_cor(method="pearson",size = 4)+theme(text = element_text(size = 13))     
#549 477
chrome_ps_mrna = ggplot(chromesome_distribution, aes(x=psrna, y=mrna)) + geom_point(color="#44D4EB")+ theme_classic()+geom_smooth(method = "lm", color="#DE616A")+
  labs(x = "psRNA_count")+labs(y = "mRNA_count")+stat_cor(method="pearson",size = 4)+theme(text = element_text(size = 13))  

chrome_pas_mrna =  ggplot(chromesome_distribution, aes(x=pasrna, y=mrna)) + geom_point(color="#44D4EB")+ theme_classic()+geom_smooth(method = "lm", color="#DE616A")+
  labs(x = "pasRNA_count")+labs(y = "mRNA_count")+stat_cor(method="pearson",size = 4)+theme(text = element_text(size = 13)) 
chrome_distr = ggarrange(chrome_ps_pas,chrome_ps_mrna,chrome_pas_mrna,nrow = 1)

ggsave(filename = 'H:/fig_project/chrome_distr.pdf',plot = chrome_distr,width = 8.27,height =2.69,units = 'in')

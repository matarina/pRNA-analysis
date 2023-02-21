library(dplyr)
library(ggplot2)
library(ggsci)
library(yyplot)
#fig1A

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
df = pd.DataFrame([
  ['allTSS',34698, 'Total'],
  ['distal TSS', 16455, 'Total'],
  ['psRNA', 13449, 'distal TSS'],
  ['pasRNA',13658, 'distal TSS']],
  columns=['country', 'pop', 'continent'])

fig = px.sunburst(df, names='country', values='pop',parents='continent', title='',labels = 'pop',
                  color_discrete_sequence=["#79B4B7","#D7E9F7"])
fig.update_traces(textinfo="label+percent parent")
fig.update_layout( font=dict(
  family="arial",
  size=18))

fig.show()
fig.write_image("H:/prna_transcript_analysis/annotationnewplot/fig/fig1a.pdf")


#fig2A
library(dplyr)
library(tibble)
library(reshape2)
ps = read.csv("H:/prna_transcript_analysis/annotationnewplot/pca/ps_output.txt",sep = '\t',header = TRUE,row.names = 1,skip = 1)#from pca/featurecount/p(a)s/ directory
pas = read.csv("H:/prna_transcript_analysis/annotationnewplot/pca/pas_output.txt",sep = '\t',header = TRUE,row.names = 1,skip = 1)
rownames(ps) = gsub('^','ps_',rownames(ps))
rownames(pas) = gsub('^','pas_',rownames(pas))
output = rbind(ps,pas)
ps_length = read.csv('H:/prna_transcript_analysis/annotationnewplot/pca/ps_length.count',sep = '\t',skip = 1)
pas_length = read.csv('H:/prna_transcript_analysis/annotationnewplot/pca/pas_length.count',sep = '\t',skip = 1)
tpm <- function(raw_counts, gene_lengths) {
  x <- raw_counts*1e3 / gene_lengths
  return(t(t(x)*1e6 / colSums(x)))
  
}
length = c(ps_length$Length,pas_length$Length)
tpm_output = tpm(output,length)



colddata = read.csv("H:/prna_transcript_analysis/annotationnewplot/pca/coldata.csv")
#colddata$X = gsub("\\..*","",colddata$X)
mix = as.data.frame(t(merge(t(tpm_output),colddata, by.x = 'row.names',by.y = 'X')))
colnames(mix) = mix[rownames(mix) %in% c("condition"),]
encode_id = t(mix['cell',])
mix = mix[!rownames(mix) %in% c('Row.names',"type","cell","condition"),]
mix_ps = mix[grep('ps',rownames(mix)),]
mix_pas = mix[grep('pas',rownames(mix)),]
####cell line specific expressed prna bar plot
mix_ps =  mix_ps[,order(names(mix_ps))]
mix_pas =  mix_pas[,order(names(mix_pas))]
stat1 = apply(apply(mix_ps,2,function(x) as.numeric(x)),2,function(x) table(x>=5))##here set cpm cutoff
stat2 = apply(apply(mix_pas,2,function(x) as.numeric(x)),2,function(x) table(x>=5))
### replace stat1/stat2 to plot bar scheme
melt_stat1 = melt(stat1) %>% filter(Var1 == 'TRUE')
melt_stat1$group = gsub("\\..*","",melt_stat1$Var2)

plot1 = merge(melt_stat1,encode_id,by.x ='group' , by.y = 'row.names' ) %>% arrange(desc(value)) %>% mutate(Var2= factor(Var2, levels=Var2))

#library(ggplot2)
#plot_melt_stat1 = melt_stat1 %>% arrange(desc(value)) %>% mutate(Var2= factor(Var2, levels=Var2))
mypal <- c('#BFA2DB','#0E5E6F','#DD5353','#FFD4D4','#829460','#F32424','#FF8D29','#ADDDD0','#D36B00','#319DA0',
                    '#7895B2','#A1E3D8','#2C3639','#4B5D67','#E3ACF9','#CCD6A6','#90A17D','#57CC99',
                    '#A25B5B','#ECB390','#285430','#91D8E4','#557153','#BFA2DB','#0E5E6F','#DD5353','#FFD4D4','#829460','#F32424','#FF8D29','#ADDDD0','#D36B00','#319DA0',
                    '#7895B2','#A1E3D8','#2C3639','#4B5D67','#E3ACF9','#CCD6A6','#90A17D','#57CC99',
                    '#A25B5B','#ECB390','#285430')
fig2a = ggplot(plot1, aes(x = cell, y = value, fill = group)) +
                       geom_col(stat = "identity", position = position_dodge2(preserve="single"),) +
                       labs(x = "Group", y = "Value", fill = "Var2") +labs(x = "Samples")+
  labs(y = "Expressed psRNA Count")+theme_classic() +scale_fill_manual(values = mypal) + 
  theme(panel.background = element_blank(),
        legend.text = element_text(size = 7),
        #legend.position="bottom",legend.direction = "vertical",
        #axis.text.x=element_blank(),
        legend.title=element_blank(),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.ticks.x.bottom = element_line(linewidth = 1),
        axis.ticks.length.x.bottom = unit(1, "mm"),
        axis.ticks.length.y.left  = unit(1, "mm"),
        axis.line.y.left = element_line(linewidth = 1),
        axis.ticks.y.left = element_line(linewidth = 1),
        axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 7),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin=unit(c(1, 1, 1, 1),'cm'))
fig2a
# p = ggplot(data=plot_melt_stat1, aes(x=Var2,weight = value,fill=group)) +geom_bar(position='dodge')+ 
#   labs(x = "Samples")+labs(y = "Expressed psRNA Count")+theme_classic() +scale_fill_manual(values = mypal)+
#   theme(panel.background = element_blank(),
#         legend.text = element_text(size = 45),
#         axis.text.x=element_blank(),
#         legend.title=element_blank(),
#         axis.line.x.bottom = element_line(linewidth = 2),
#         axis.line.y.left = element_line(linewidth = 2),
#         axis.text.y  = element_text(size = 45),
#         axis.title.x = element_text(size = 55),
#         axis.title.y = element_text(size = 55),
#         plot.margin=unit(c(1, 1, 1, 1),'cm'))  

#p = p + scale_y_continuous(expand = c(0,0),
#                           limits = c(0,5000))




#fig2b

melt_stat2 = melt(stat2) %>% filter(Var1 == 'TRUE')
melt_stat2$group = gsub("\\..*","",melt_stat2$Var2)

plot2 = merge(melt_stat2,encode_id,by.x ='group' , by.y = 'row.names' ) %>% arrange(desc(value)) %>% mutate(Var2= factor(Var2, levels=Var2))

#library(ggplot2)
#plot_melt_stat2 = melt_stat2 %>% arrange(desc(value)) %>% mutate(Var2= factor(Var2, levels=Var2))
mypal <- c('#BFA2DB','#0E5E6F','#DD5353','#FFD4D4','#829460','#F32424','#FF8D29','#ADDDD0','#D36B00','#319DA0',
                    '#7895B2','#A1E3D8','#2C3639','#4B5D67','#E3ACF9','#CCD6A6','#90A17D','#57CC99',
                    '#A25B5B','#ECB390','#285430','#91D8E4','#557153','#BFA2DB','#0E5E6F','#DD5353','#FFD4D4','#829460','#F32424','#FF8D29','#ADDDD0','#D36B00','#319DA0',
                    '#7895B2','#A1E3D8','#2C3639','#4B5D67','#E3ACF9','#CCD6A6','#90A17D','#57CC99',
                    '#A25B5B','#ECB390','#285430')
fig2b = ggplot(plot2, aes(x = cell, y = value, fill = group)) +
                      geom_col(stat = "identity", position = position_dodge2(preserve="single",),) +
  labs(x = "Group", y = "Value", fill = "Var2") +labs(x = "Samples")+
                      labs(y = "Expressed pasRNA Count")+theme_classic() +scale_fill_manual(values = mypal) + 
  theme(panel.background = element_blank(),
        legend.text = element_text(size = 7),
        #legend.position="bottom",legend.direction = "vertical",
        #axis.text.x=element_blank(),
        legend.title=element_blank(),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.ticks.x.bottom = element_line(linewidth = 1),
        axis.ticks.length.x.bottom = unit(1, "mm"),
        axis.ticks.length.y.left  = unit(1, "mm"),
        axis.line.y.left = element_line(linewidth = 1),
        axis.ticks.y.left = element_line(linewidth = 1),
        axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 7),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin=unit(c(1, 1, 1, 1),'cm'))

fig2b
fig2a_b = ggarrange(fig2a,fig2b,nrow = 2 ,align = 'v',common.legend = TRUE,legend = 'right')
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig2a_b.pdf',plot = fig2a_b,width = 8.27,height =5.69,units = 'in')

##fig2c
rm(list=ls())
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(ggthemes)
library(tidyr)
library(dplyr)
sense = read.csv("H:/prna_transcript_analysis/annotationnewplot/pca/ps_output.txt",sep = '\t',row.names = 1,skip = 1)
coldata = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/sense/single_coldata.csv")
#coldata = coldata[grep(pattern = 'iPSC|foreskin_fiber|lymphoblasts', x = coldata$cell),]
coldata$condition = gsub("-","_",coldata$condition)
sense = sense[,which(colnames(sense) %in% coldata$X)]
coldata = coldata[order(coldata$X),]
sense = sense[,order(colnames(sense))]
#colnames(sense) = coldata$cell
##delete abnormal data
sense = sense[,-c(2,8,22)]
coldata = coldata[-c(2,8,22),]
dds <- DESeqDataSetFromMatrix(countData = sense,
                              colData = coldata,
                              design = ~ cell)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
#if error in antisense
#vsd = varianceStabilizingTransformation(dds,blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
pca = plotPCA(vsd, intgroup=c("cell"),returnData = TRUE)
mypal <- c('#219F94','#C65D7B','#325288','#C3B091','#FC997C')
fig2c = ggplot(pca, aes(x=PC1, y=PC2,color=group))+ geom_point(size=0.5)+ stat_ellipse(level = 0.999, lwd = 0.4)+
  scale_color_manual(values= mypal)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 
fig2c

ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig2c.pdf',plot = fig2c,width = 4,height =2.6,units = 'in')


#fig2d
rm(list=ls())
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(ggthemes)
library(tidyr)
library(dplyr)
sense = read.csv("H:/prna_transcript_analysis/annotationnewplot/pca/pas_output.txt",sep = '\t',row.names = 1,skip = 1)
coldata = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/sense/single_coldata.csv")
#coldata = coldata[grep(pattern = 'iPSC|foreskin_fiber|lymphoblasts', x = coldata$cell),]
coldata$condition = gsub("-","_",coldata$condition)
sense = sense[,which(colnames(sense) %in% coldata$X)]
coldata = coldata[order(coldata$X),]
sense = sense[,order(colnames(sense))]
#colnames(sense) = coldata$cell
##delete abnormal data
sense = sense[,-c(2,8,22)]
coldata = coldata[-c(2,8,22),]
dds <- DESeqDataSetFromMatrix(countData = sense,
                              colData = coldata,
                              design = ~ cell)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
#if error in antisense
#vsd = varianceStabilizingTransformation(dds,blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
pca = plotPCA(vsd, intgroup=c("cell"),returnData = TRUE)
mypal <- c('#219F94','#C65D7B','#325288','#C3B091','#FC997C')
fig2d = ggplot(pca, aes(x=PC1, y=PC2,color=group))+ geom_point(size=0.5)+ stat_ellipse(level = 0.999, lwd = 0.4)+
  scale_color_manual(values= mypal)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 
fig2d

ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig2d.pdf',plot = fig2d,width = 4,height =2.6,units = 'in')

#fig2c_d = ggarrange(fig2c,fig2d,align = 'h')
#ggsave(filename = 'H:/fig_project/fig2c_d.pdf',plot = fig2c_d,width = 8.27,height =5.69,units = 'in')

##fig2e
sense = tpm_output[grep('ps_',rownames(tpm_output)),]
coldata = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/sense/single_coldata.csv")
#coldata$X = gsub("\\..*","",coldata$X)
#coldata = coldata[grep(pattern = 'iPSC|foreskin_fiber|lymphoblasts', x = coldata$cell),]
sense = sense[,which(colnames(sense) %in% coldata$X)]
coldata = coldata[order(coldata$X),]
sense = as.data.frame(sense[,order(colnames(sense))])
colnames(sense) = coldata$cell

stem_cell = sense[,grep("stem_cell",colnames(sense))]
stem_cell$mean_stem_cell = rowMeans(stem_cell)
stem_cell = filter(stem_cell,mean_stem_cell>=5)

foreskin_fiber = sense[,grep("foreskin_fiber",colnames(sense))]
foreskin_fiber$mean_foreskin_fiber = rowMeans(foreskin_fiber)
foreskin_fiber = filter(foreskin_fiber,mean_foreskin_fiber>=5)

prostate = sense[,c(35,36)]

prostate = filter(prostate,prostate>=5)

breast = sense[,grep("breast",colnames(sense))]
breast$mean_breast = rowMeans(breast)
breast = filter(breast,mean_breast>=5)

lymphoblasts = sense[,grep("lymphoblasts",colnames(sense))]
lymphoblasts$mean_lymphoblasts = rowMeans(lymphoblasts)
lymphoblasts = filter(lymphoblasts,mean_lymphoblasts>=5)

x = list(a = rownames(stem_cell),
         b = rownames(foreskin_fiber),
         c = rownames(prostate),
         d = rownames(lymphoblasts),
         e = rownames(breast))

##plot 5 set venn diagram
library("extrafont")
library(VennDiagram)
library(ggVennDiagram)
#font_import()
#loadfonts() 
venn.diagram(
  x,filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig2e.png',
  category.names = c("stem_cell" , "foreskin_fiber " , "prostate", "lymphoblasts","breast"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#FF8787", "#FFD372", "#C6EBC5", "#BAD7E9","#579BB1"),
  # Numbers
  cex = 1.7,
  #fontface = "arial",
  cat.fontfamily ="Arial", ##category font
  fontfamily ="Arial",  ##number font
  # Set names
  cat.cex = 1.9,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(1),
  cat.dist = c(0.055, 0.081, -0.1, -0.1,0.1),
  height = 18,
  width = 18, resolution = 300, imagetype = "png",
  units = "cm"
  )



##fig2f
sense = tpm_output[grep('pas_',rownames(tpm_output)),]
coldata = read.csv("H:/prna_transcript_analysis/epd_promoter/encode_data_pca/sense/single_coldata.csv")
#coldata$X = gsub("\\..*","",coldata$X)
#coldata = coldata[grep(pattern = 'iPSC|foreskin_fiber|lymphoblasts', x = coldata$cell),]
sense = sense[,which(colnames(sense) %in% coldata$X)]
coldata = coldata[order(coldata$X),]
sense = as.data.frame(sense[,order(colnames(sense))])
colnames(sense) = coldata$cell

stem_cell = sense[,grep("stem_cell",colnames(sense))]
stem_cell$mean_stem_cell = rowMeans(stem_cell)
stem_cell = filter(stem_cell,mean_stem_cell>=5)

foreskin_fiber = sense[,grep("foreskin_fiber",colnames(sense))]
foreskin_fiber$mean_foreskin_fiber = rowMeans(foreskin_fiber)
foreskin_fiber = filter(foreskin_fiber,mean_foreskin_fiber>=5)

prostate = sense[,c(35,36)]

prostate = filter(prostate,prostate>=5)

breast = sense[,grep("breast",colnames(sense))]
breast$mean_breast = rowMeans(breast)
breast = filter(breast,mean_breast>=5)

lymphoblasts = sense[,grep("lymphoblasts",colnames(sense))]
lymphoblasts$mean_lymphoblasts = rowMeans(lymphoblasts)
lymphoblasts = filter(lymphoblasts,mean_lymphoblasts>=5)

x = list(a = rownames(stem_cell),
         b = rownames(foreskin_fiber),
         c = rownames(prostate),
         d = rownames(lymphoblasts),
         e = rownames(breast))

##plot 5 set venn diagram
library("extrafont")
library(ggVennDiagram)
#font_import()
#loadfonts() 
venn.diagram(
  x,filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig2f.png',
  category.names = c("stem_cell" , "foreskin_fiber " , "prostate", "lymphoblasts","breast"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#FF8787", "#FFD372", "#C6EBC5", "#BAD7E9","#579BB1"),
  # Numbers
  cex = 1.7,
  #fontface = "arial",
  cat.fontfamily ="Arial", ##category font
  fontfamily ="Arial",  ##number font
  # Set names
  cat.cex = 1.9,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(1),
  cat.dist = c(0.055, 0.081, -0.1, -0.1,0.1),
  height = 18,
  width = 18, resolution = 300, imagetype = "png",
  units = "cm"
)

##correlation
rm(list=ls())
setwd("H:/")
library(tidyverse)
ps = read.csv("H:/prna_transcript_analysis/annotationnewplot/hepg2/ps_output.txt",sep = '\t',skip = 1)
colnames(ps)[1] = 'FEATURE_ID'
ps$FEATURE_ID = gsub('^','ps_',ps$FEATURE_ID)
pas = read.csv("H:/prna_transcript_analysis/annotationnewplot/hepg2/pas_output.txt",sep = '\t',skip = 1)
colnames(pas)[1] = 'FEATURE_ID'
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
mypal = c('#03C988','#325288')
ps_mrna = ggplot(b, aes(x=mrna_mean, y=ps_mean,color = group)) + geom_point(size = 1)+ theme_classic()+
  labs(x = "mRNA")+labs(y = "psRNA")+scale_x_continuous(limits = c(NA, 1500))+scale_y_continuous(limits = c(NA,75))+
  scale_color_manual(values= mypal)+theme_classic()+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 



ps_mrna
#correlation2 
b = merge(emrna2,epas,by = 'epd')
b$group = 'others'
library(biomaRt)
# define biomart object

#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
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

mypal = c('#C65D7B','#325288')
pas_mrna = ggplot(b, aes(x=mrna_mean, y=pas_mean,color = group)) + geom_point(size = 1)+ theme_classic()+
  labs(x = "mRNA")+labs(y = "pasRNA")+scale_x_continuous(limits = c(NA, 1500))+scale_y_continuous(limits = c(NA,75))+
  scale_color_manual(values= mypal)+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 

pas_mrna
corr = ggarrange(ps_mrna,pas_mrna,align = 'v')
corr
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/corr.pdf',plot = corr,width = 8.27,height =2.69,units = 'in')

#fig3d

ps = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/ps.bed',sep = '\t',header = FALSE)
pas = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/pas.bed',sep = '\t',header = FALSE)
psnew = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/psnew.bed',sep = '\t',header = FALSE)
pasnew = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/pasnew.bed',sep = '\t',header = FALSE)
ps_diff = ps[!(ps$V4 %in% psnew$V4),]$V4
pas_diff = ps[!(pas$V4 %in% pasnew$V4),]$V4

bp3target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/ps_bp3_merip.bed",header = FALSE,sep = "\t")
bp1target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/ps_bp1_merip.bed",header = FALSE,sep = "\t")
target = rbind(bp3target,bp1target)
decaydata = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/ps_result.txt")
decaydata = decaydata[-grep("ERCC-",decaydata$X),]

decaydata = filter(decaydata,sigma2<0.066)
decaydata = decaydata[!(decaydata$X %in% pas_diff),]
bp13 = decaydata[(decaydata$X %in% bp3target$V4),]
bp13$halflife = log(2)/bp13$alpha_hepg2
mean(bp13$halflife)
bp13$group = rep("target \nmean : 208.99min",nrow(bp13))
without_bp13 = decaydata[!(decaydata$X %in% bp3target$V4),]###删除在另一个表中含有同名的行
without_bp13$halflife = log(2)/without_bp13$alpha_hepg2
mean(without_bp13$halflife)
without_bp13$group = rep("untarget \nmean :139.96min",nrow(without_bp13))
bp13_all = bind_rows(bp13,without_bp13)
ks.test(bp13$halflife,without_bp13$halflife, alternative = "l")

ps_bp = ggplot(bp13_all, aes( x=halflife,color = group))+stat_ecdf(size=0.7)+
  annotate("text", x = 300.05, y = 0.073,size = 4, label = "KS-test P-value:1.212e-4")+
  scale_color_manual(values = c("#cf3550","#048e9c"))+
  scale_x_continuous(expand = c(0, 5),limits = c(0,495))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Cumulative fraction")+labs(x = " halflife")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
        legend.title=element_blank(),legend.position = c(0.64, 0.4),
          legend.text=element_text(size=13),
          axis.text=element_text(size=13),
          axis.title=element_text(size=13),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.line = element_line(colour = "black"))+
          theme(plot.margin=unit(c(0.5, 1.5, 0.5, 0.5),'cm'))
ps_bp



##fig3e
ps = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/ps.bed',sep = '\t',header = FALSE)
pas = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/pas.bed',sep = '\t',header = FALSE)
psnew = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/psnew.bed',sep = '\t',header = FALSE)
pasnew = read.csv('H:/prna_transcript_analysis/annotationnewplot/annot_diff_genes/pasnew.bed',sep = '\t',header = FALSE)
ps_diff = ps[!(ps$V4 %in% psnew$V4),]$V4
pas_diff = ps[!(pas$V4 %in% pasnew$V4),]$V4

bp3target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/pas_bp3_merip.bed",header = FALSE,sep = "\t")
bp1target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/pas_bp1_merip.bed",header = FALSE,sep = "\t")
target = rbind(bp3target,bp1target)
decaydata = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/pas_result.txt")
decaydata = decaydata[-grep("ERCC-",decaydata$X),]
decaydata = decaydata[!(decaydata$X %in% pas_diff),]
decaydata = filter(decaydata,sigma2<0.066)
bp13 = decaydata[(decaydata$X %in% bp3target$V4),]
bp13$halflife = log(2)/bp13$alpha_hepg2
mean(bp13$halflife)
bp13$group = rep("target \nmean :177.21min",nrow(bp13))
without_bp13 = decaydata[!(decaydata$X %in% bp3target$V4),]###删除在另一个表中含有同名的行
without_bp13$halflife = log(2)/without_bp13$alpha_hepg2
mean(without_bp13$halflife)
without_bp13$group = rep("untarget \nmean :124.31min",nrow(without_bp13))
bp13_all = bind_rows(bp13,without_bp13)
ks.test(bp13$halflife,without_bp13$halflife, alternative = "l")

pas_bp = ggplot(bp13_all, aes( x=halflife,color = group))+stat_ecdf(size=0.7)+
  annotate("text", x = 300.05, y = 0.073,size = 4, label = "KS-test P-value:3.3e-10")+
  scale_color_manual(values = c("#cf3550","#048e9c"))+
  scale_x_continuous(expand = c(0, 5),limits = c(0,495))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Cumulative fraction")+labs(x = " halflife")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
        legend.title=element_blank(),legend.position = c(0.64, 0.4),
        legend.text=element_text(size=13),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"))+
  theme(plot.margin=unit(c(0.5, 1.5, 0.5, 0.5),'cm'))
pas_bp
prna_bp = ggarrange(ps_bp,pas_bp,align = 'v')
prna_bp
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig3de.pdf',plot = prna_bp,width = 8.27,height =3.69,units = 'in')
##fig3c
createDGE <- function(countMatrix = NA, annotationFile = "annotation.gtf", sampleInfo = "targets.txt", strand = 0, paired = c("YES", "NO"), nthreads = 1) {
  ## Create a DGE list by counting reads on bam files that are
  ## listed in the sampleInfo dataframe. With no arguments specified
  ## the function reads the file targets.txt and  annotation.gtf in the
  ## current directory and assumes that the data is unstranded.
  ##
  ## Args:
  ## countMatrix: Character vector length 1. Name and path to available count matrix.
  ## The file should have column header and the first column should be gene names.
  ## annotationFile: Character vector length 1. The complete path to the gtf file used
  ## for analysis. Annotation.gtf in the current folder will be used as default
  ## sampleInfo: Character vector length 1. Name of the file containing
  ## information on the bam files to be analysed.
  ## This file needs to have a column named "sample_name" containing
  ## filenames and path to all input files
  ## strand: Numeric vector of length 1. Is the data stranded or not.
  ## 0 = unstranded,
  ## 1 = stranded with first read in direction of annotation,
  ## 2 = stranded with first read opposite of annotation (Typical for Illumina)
  ## paired: Character vector length 1. Is the data paired end (YES), or not (NO)
  ## nthreads: Numeric vector length 1. Number of threads used for counting reads.
  wd <- getwd()
  if (is.na(countMatrix)) {
    pe = match.arg(paired, c("YES", "NO"))
    sampleInfo <- read.table(sampleInfo, header = TRUE)
    bamFiles = sampleInfo$sample_name
    fc <- featureCounts(files = bamFiles, annot.ext = annotationFile,
                        isPaired = pe, nthreads = nthreads,
                        isGTFAnnotationFile = TRUE, strandSpecific = strand)
    dge <- DGEList(counts = fc$counts, genes = fc.annotation)
    dge
  } else {
    cm <- read.table(paste0(wd,"/",countMatrix), header = TRUE)
    dge <- DGEList(counts = cm[,-1], genes = cm[,1])
    dge
  }
}

#### 
apply_expression_cutoff <- function(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10){
  ## Takes a DGEList object as an input and outputs the filtered DGEList object
  ## The samples must be named after the time point "t" followed by a number, eg "t0", "t1", "t2"
  ## filter on cpm of t0, t1, t2
  ## filter on sum of raw count
  ## recalculates library size
  ## prints to screen the number of genes filtered out and how many genes are left
  
  obj <- DGEList
  all_rows <- nrow(obj)
  keep <- rowSums(cpm(obj)) > CountCutoff
  obj <- obj[keep, , keep.lib.sizes = FALSE]
  countfilter <- nrow(obj)
  
  if(all(c("t0", "t1","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t1"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t0", "t1") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t1"]) > earlyCPM_cutoff 
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t0","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t1","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t1"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else{stop("There are not enough early time points ... exiting ...")}
  
  cat("There were ", all_rows, " genes in the input. \n", "After filtering for sum of row counts (", CountCutoff, "), there are ",countfilter, " left. \n", "After filtering for CPM at t=0 (", earlyCPM_cutoff,"), there are ", cpmfilter, " genes left. \n", sep="")
  return(obj)
}

#####


trim_late_time_points <- function(DGeList =obj, CPMcutoff= 0.5){
  #replaces late time points with cpm lower than threshold with NA
  #replaces all time points later than a NA time point with NA
  #if there is a cpm increase of more than 20 % between two time points, turn into NA
  #stores the "trimmed cpm" data in obj$trimmed
  q <- ifelse(cpm(obj)<CPMcutoff,NA,cpm(obj))
  q <- apply(q,1, function(x){
    for(i in 1:(length(x)-1)){
      if(((x[i+1]-x[i])/x[i]) > 0.2 & !is.na(x[i]) & !is.na(x[i+1])){
        x[i+1] <-NA
      }
      if(is.na(x[i])) {
        x[i+1] <- NA
      }
    }
    return(x)
  })
  
  obj$cpm_trimmed <- t(q)
  obj
  
  
  
  # q <- apply(q,1, function(x){
  #   for(i in 1:(length(x))){
  #     if(i < length(x)){
  #       if(((x[i+1]-x[i])/x[i]) > 0.2 & !is.na(x[i]) & !is.na(x[i+1])){
  #         x[i+1] <-NA
  #       }
  #     }
  #     if(is.na(x[i])) {
  #       x[i+1] <- NA
  #     }
  #   }
  #   return(x)
  # })
  # 
  # obj$cpm_trimmed <- t(q)
  # obj
  
}
#Use non transformed data set fit nls

calculate_normalization_factors3 <- function(DGEList = obj, method = c("mean", "median", "peak")){
  
  #get the data, only keep the complete rows for calculation of correction coefficients
  data <- DGEList$cpm_trimmed
  data <- data[complete.cases(data),]
  
  ##create a time point vector that fits the main data
  t <- sapply(strsplit(colnames(data),"t"), function(x){x[!x ==""]})
  
  ##loop over the genes to use for normalization and collect normalization factors
  ##normalization factors are coefficients that transform the curve into a perfect logarithmic decay
  corr_coeffs <- as.data.frame(data[0,])
  for(i in 1:nrow(data)){
    # adjust the length of the time point vector to the actual data
    names(corr_coeffs) <- t
    t <- as.numeric(t)
    realt <- suppressWarnings(t+data[i,]-data[i,])
    
    #fit
    values <- data[i,]
    fit = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= data[i,1],b=0.1)), error=function(e){print(paste("fitting failed row:" ,i))})
    
  }
  for(i in 1:nrow(data)){
    # adjust the length of the time point vector to the actual data
    names(corr_coeffs) <- t
    t <- as.numeric(t)
    realt <- suppressWarnings(t+data[i,]-data[i,])
    
    #fit
    values <- data[i,]
    fit = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= data[i,1],b=0.1)), error=function(e){print(paste("fitting failed row:" ,i))})
    
    ifelse(grepl("fitting failed",fit),{corr_coeffs <- rbind(corr_coeffs, values*NA)}, {corr_coeffs <- rbind(corr_coeffs,fitted(fit)/data[i,])})
    
    
    
  }
  
  # write new normalization factors in the object
  corr_coeffs_no_NA <- corr_coeffs[complete.cases(corr_coeffs),]
  if(method == "mean"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, mean)
  }
  if(method == "median"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, median)
  }
  if(method == "peak"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, function(z){
      den <- density(z)
      den$x[which.max(den$y)]
    } #closes function(z)
    )#close apply
  }#closes if
  
  print(paste("The calculated normalization factors were : "))
  print(DGEList$samples$norm.factors)
  
  ###PLots the distribution of the normalization factors for each time point
  par(mfrow = c(3,2), main = method,oma = c(0, 0, 2, 0))
  for(i in 1:length(DGEList$samples$norm.factors)){
    plot(density(corr_coeffs_no_NA[,i]), main = paste(" t =", colnames(corr_coeffs_no_NA)[i]))
  }
  mtext(paste("Distribution of correction coefficients"), outer = TRUE, cex = 1.5)
  
  DGEList$cpm_trimmed_normalized <- t(t(DGEList$cpm_trimmed)*DGEList$samples$norm.factors)
  
  return(DGEList)
}#closes function


calculate_half_life <- function(DGEList = obj){
  
  data <- DGEList$cpm_trimmed_normalized
  ##create a time point vector that fits the main data
  t <- sapply(strsplit(colnames(data),"t"), function(x){x[!x ==""]})
  
  results <- as.data.frame(data[0,])
  
  
  
  for(i in 1:nrow(data)){
    # adjust the length of the time point vector to the actual data
    t <- as.numeric(t)
    realt <- suppressWarnings(t+data[i,]-data[i,])
    
    #logtranform data
    log2data <- apply(data,2,log2)
    
    #fit_nls
    values <- data[i,]
    fit_nls = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= data[i,1],b=0.1)), error=function(e){paste("nls fitting failed row:" ,i)})
    
    ifelse(grepl("fitting failed",fit_nls),{b <- NA ; hl_nls <- NA},{ b <- coef(fit_nls)["b"] ;hl_nls <- log(2)/b })
    #fit_lm
    values <- log2data[i,]
    fit_lm <- tryCatch(lm(formula =  values ~ realt), error=function(e){paste("lm fitting failed row:" ,i)})
    
    ifelse(grepl("fitting failed",fit_lm),{hl_lm <- NA},{hl_lm <- -1/coef(fit_lm)[[2]]})
    
    results <- rbind(results, c(hl_nls, hl_lm))
  }#end for loop
  names(results) <- c("hl_nls","hl_lm")
  results <- cbind(DGEList$genes, results)
  return(results)
}#end function
## calculate mrna vs ps /pas halflife
library(dplyr)
setwd("H:/")
mrna = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/mrna_output.txt",sep = '\t')
mrna = mrna[-grep("ERCC-",mrna$FEATURE_ID),]


mrna_rep1 = mrna[,c(1,2,4,6,8)]
colnames(mrna_rep1) = c("Geneid","t0","t1","t3","t6")
mrna_rep2 = mrna[,c(1,3,5,7,9)]
colnames(mrna_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = mrna_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
mrna_rep1_halflife = obj


cm = mrna_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
mrna_rep2_halflife = obj


mrna = merge(mrna_rep1_halflife,mrna_rep2_halflife,by = "genes")
mrna$mean_nls = rowMeans(subset(mrna,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

mrna_filtered = filter(mrna,mean_nls>0)
mean(mrna_filtered$mean_nls,trim = 0.2)

#mean mrna halflife(trim 0.2)  8.843557h = 530.6134min



#calculate total psRNA halflife 
library(dplyr)
setwd("H:/")
ps = read.csv("H:/prna_transcript_analysis/annotationnewplot/hepg2/ps_output.txt",skip = 1,sep = '\t')
colnames(ps)[1] = 'FEATURE_ID'
ps = ps[-grep("ERCC-",ps$FEATURE_ID),]


ps_rep1 = ps[,c(1,2,4,6,8)]
colnames(ps_rep1) = c("Geneid","t0","t1","t3","t6")
ps_rep2 = ps[,c(1,3,5,7,9)]
colnames(ps_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = ps_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
ps_rep1_halflife = obj


cm = ps_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
ps_rep2_halflife = obj


ps = merge(ps_rep1_halflife,ps_rep2_halflife,by = "genes")
ps$mean_nls = rowMeans(subset(ps,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

ps_filtered = filter(ps,mean_nls>0)
mean(ps_filtered$mean_nls,trim = 0.2)
#ps total halflife 4.022274h = 241.3364min



#calculate total pasRNA halflife
library(dplyr)
setwd("H:/")
pas = read.csv("H:/prna_transcript_analysis/annotationnewplot/hepg2/pas_output.txt",skip = 1,sep = '\t')
colnames(pas)[1] = 'FEATURE_ID'
pas = pas[-grep("ERCC-",pas$FEATURE_ID),]

pas_rep1 = pas[,c(1,2,4,6,8)]
colnames(pas_rep1) = c("Geneid","t0","t1","t3","t6")
pas_rep2 = pas[,c(1,3,5,7,9)]
colnames(pas_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = pas_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
pas_rep1_halflife = obj


cm = pas_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
pas_rep2_halflife = obj


pas = merge(pas_rep1_halflife,pas_rep2_halflife,by = "genes")
pas$mean_nls = rowMeans(subset(pas,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

pas_filtered = filter(pas,mean_nls>0)
mean(pas_filtered$mean_nls,trim = 0.2)
mean(pas_filtered$mean_nls,trim = 0.2)
pas_filtered$min = pas_filtered$mean_nls*60
pas_filtered$group = rep("pasRNA\nmean :252.20min",857)

mean(ps_filtered$mean_nls,trim = 0.2)
ps_filtered$group = rep("psRNA\nmean :241.34min",523)
ps_filtered$min = ps_filtered$mean_nls*60

mean(mrna_filtered$mean_nls,trim = 0.2)
mrna_filtered$group = rep("mRNA\nmean :530.61min",6871)
mrna_filtered$min = mrna_filtered$mean_nls*60
ks.test(mrna_filtered$min,pas_filtered$min, alternative = "l")
ks.test(mrna_filtered$min,ps_filtered$min, alternative = "l")
mrna_vs_ps = rbind(mrna_filtered,pas_filtered)
mrna_vs = rbind(mrna_vs_ps,ps_filtered)
library(ggplot2)
prna_mrna = ggplot(mrna_vs, aes( x=min,color = group))+stat_ecdf(size=0.7)+
  annotate("text", x = 300.05, y = 0.073,size = 4, label = "KS-test P-value:2.2e-16")+
  scale_color_manual(values = c("#13274f","#ce1141","#e7a801"))+
  scale_x_continuous(expand = c(0, 5),limits = c(0,495))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Cumulative fraction")+labs(x = " halflife")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
        legend.title=element_blank(),legend.position = c(0.64, 0.4),
        legend.text=element_text(size=13),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"))+
  theme(plot.margin=unit(c(0.5, 1.5, 0.5, 0.5),'cm'))
prna_mrna = prna_mrna+ scale_x_continuous(limits = c(0,600))
prna_mrna
prna_decay  = ggarrange(ps_bp,pas_bp,align = 'v')
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig3c.pdf',plot = prna_mrna,width = 4.27,height =3.69,units = 'in')
ggsave(filename = 'H:/fig_project/fig3d_e.pdf',plot = prna_decay,width = 8.27,height =3.69,units = 'in')



halflife1 = data.frame(a = c('a','b','c'),value = c(530.61,131.11,144.57))
halflife1 = ggplot(halflife1,aes(x = a,y = value ,fill = a))+geom_bar(stat="identity")+scale_fill_manual(values = c("#13274f","#ce1141","#e7a801"))+theme_classic()
halflife2 = data.frame(a = c('a','b'),value = c(208.99,139.96))
halflife2 = ggplot(halflife2,aes(x = a,y = value ,fill = a))+geom_bar(stat="identity")+scale_fill_manual(values = c("#cf3550","#048e9c"))+theme_classic()
halflife3 = data.frame(a = c('a','b'),value = c(177.21,124.31))
halflife3 = ggplot(halflife3,aes(x = a,y = value ,fill = a))+geom_bar(stat="identity")+scale_fill_manual(values = c("#cf3550","#048e9c"))+theme_classic()
halflife = ggarrange(halflife1,halflife2,halflife3,ncol = 3)
halflife
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/fig3cde_halflife.pdf',plot = halflife,width = 8.27,height =3.69,units = 'in')


##fig3b
library(cliProfiler)
library(rtracklayer)
library(ChIPpeakAnno)
## Extract all the exon annotation
test_anno <- rtracklayer::import.gff3('H:/prna_transcript_analysis/annotationnewplot/decay/ps.gtf')
test_anno <- test_anno[test_anno$type == "exon"]


peak<-toGRanges('H:/prna_transcript_analysis/annotationnewplot/decay/ps_peaks.bed',format=c('BED'))
## Run the windowProfile
window_profile <- windowProfile(peak, test_anno)

profile1 = window_profile[[2]]+scale_y_continuous(expand = c(0,0))+geom_density(adjust = 1/5 ,linewidth = 1.2,color="#cf3550")+
  labs(x = "psRNA profile")+labs(y = "Density of peaks")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.line.y.left = element_line(linewidth = 1),
        axis.ticks.length.x.bottom = unit(1, "mm"),
        axis.ticks.length.y.left  = unit(1, "mm"),
        axis.ticks.y.left = element_line(linewidth = 1),
        axis.ticks.x.bottom  = element_line(linewidth = 1),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 



test_anno <- rtracklayer::import.gff3('H:/prna_transcript_analysis/annotationnewplot/decay/pas.gtf')
test_anno <- test_anno[test_anno$type == "exon"]

peak2 <-toGRanges('H:/prna_transcript_analysis/annotationnewplot/decay/pas_peaks.bed',format=c('BED'))
## Run the windowProfile
window_profile2 <- windowProfile(peak2, test_anno)
profile2 = window_profile2[[2]]+scale_y_continuous(expand = c(0,0))+geom_density(adjust = 1/5 ,linewidth = 1.2,color="#cf3550")+
  labs(x = "pasRNA profile")+labs(y = "Density of peaks")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.line.y.left = element_line(linewidth = 1),
        axis.ticks.length.x.bottom = unit(1, "mm"),
        axis.ticks.length.y.left  = unit(1, "mm"),
        axis.ticks.y.left = element_line(linewidth = 1),
        axis.ticks.x.bottom  = element_line(linewidth = 1),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/profile_pas.pdf',plot = profile2,width = 4.27,height =3.69,units = 'in')
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig/profile_ps.pdf',plot = profile1,width = 4.27,height =3.69,units = 'in')


#fig5a plasma expressed prna analysis
specie <- c(rep("mRNA" , 2) , rep("psRNA" , 2) , rep("pasRNA" , 2))
condition <- rep(c("cancer" , "normal") , 3)
value <- as.numeric(c("9258","8733","2166","3150","1500","1730"))
data <- data.frame(specie,condition,value)

# Grouped
plasma_prna_stat = ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("#153462","#FEC260"))+
  #  geom_text(aes(label = len), size = 15,color="white",vjust= 1.81, position = "stack")+
  labs(x = "TSS")+labs(y = "Numbers")+theme_classic() +
  theme(panel.background = element_blank(),legend.position="bottom",
        legend.text = element_text(size = 13),
        legend.title=element_blank(),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.line.y.left = element_line(linewidth = 1),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        plot.margin=unit(c(1, 1, 1, 1),'cm'))
plasma_prna_stat



#fig5b 
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
keyvals <- ifelse(
  res1$log2FoldChange < -2.5, '#0081B4',
  ifelse(res1$log2FoldChange > 2.5, '#ff6361',
         'black'))
         keyvals[is.na(keyvals)] <- '#B6EADA'
           names(keyvals)[keyvals == '#ff6361'] <- 'high'
           names(keyvals)[keyvals == '#B6EADA'] <- 'mid'
           names(keyvals)[keyvals == '#0081B4'] <- 'low'
plasma_volcano = EnhancedVolcano(res1,legendPosition = 'none',
                gridlines.minor=FALSE, selectLab = rownames(res1)[which(names(keyvals) %in% c('high', 'low'))],
                gridlines.major=FALSE,colCustom = keyvals,
                lab =NA,title = NULL,subtitle = NULL,#legendLabels = NULL,
                caption = '',hline = NULL,vline = NULL,
                x = 'log2FoldChange',
                y = 'padj',parseLabels = FALSE,
                border = 'full',
                borderWidth = 1,pointSize = 2.0,
                axisLabSize = 6,
                col=c("#E69F00","#86A3B8","#BDCDD6","#ff6361"),
                pCutoff = 0.05)


plasma_volcano

##fig5c

#rm(list=ls())
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
### remove a abnormal data
data = data[,-64]

coldata = data.frame(sample = colnames(data),condition = rep(c('cancer','normal'),times=c(35,29)))

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
plasma_pca = ggplot(pca, aes(x=PC1, y=PC2,color=condition))+ geom_point(size=1)+ 
  scale_color_manual(values= c("#153462","#FEC260"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.line.y.left = element_line(linewidth = 1),
        axis.ticks.length.x.bottom = unit(1, "mm"),
        axis.ticks.length.y.left  = unit(1, "mm"),
        axis.ticks.y.left = element_line(linewidth = 1),
        axis.ticks.x.bottom  = element_line(linewidth = 1),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 
plasma_pca
plasma = ggarrange(plasma_prna_stat,plasma_pca,plasma_volcano,ncol = 3)
plasma
ggsave(filename = 'H:/fig_project/plasma_pca.pdf',plot = plasma_pca,width = 4.27,height =3.69,units = 'in')
ggsave(filename = 'H:/fig_project/plasma_volcano.pdf',plot = plasma_volcano,width = 4.27,height =3.69,units = 'in')


#fig6a plasma mechine learning
library(extrafont)
loadfonts(device = "win")
rfefit = readRDS('H:/prna_transcript_analysis/plasma/plasma_rfefit.RDS')
rfepl <- ggplot(data=rfeAccuracy, aes(x=Variables, y=Accuracy))+
  geom_point(color="grey", size=2, shape=19)+ xlim(NA, 30)+
  geom_line(color="black", linetype=1, size=1)+
  # geom_text(aes(label=label), nudge_y=0.002)+
  annotate(geom="point",
           x=rfeAccuracy[grep(max(rfeAccuracy$Accuracy), rfeAccuracy$Accuracy), ]$Variables,
           y=max(rfeAccuracy$Accuracy), color="#d74a49", size=4)+
  labs(x="Features (Numbers)",
       y="Accuracy (Bootstrap)")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(size=.5, color="black"),
        axis.title = element_text(color='black', size=13),
        axis.text = element_text(color='black', size=13),
        text = element_text(size=13, color="black"))
rfepl
ggsave('H:/fig_project/plasma_feature_number.pdf',rfepl,width = 4.27,height =3.69,units = 'in')


#fig6b
library("heatmaply")
library(pheatmap)
library(ComplexHeatmap)
library("tinyarray")
library(circlize)
library(caret)
expr_prna = column_to_rownames(merge(expr_cancer,expr_normal,by = 'row.names'),var = 'Row.names')
rfefit = readRDS('H:/prna_transcript_analysis/plasma/plasma_rfefit.RDS')
rfeImp <- varImp(rfefit) %>% 
  rownames_to_column("feature") %>%
  filter(feature%in%rfefit$optVariables) %>%
  arrange(desc(Overall))
biomarker_expr =expr_prna[rfeImp$feature,]



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
png("H:/fig_project/plasma_heatmap.png",width=13,height=10,units="cm",res=300)
Heatmap(t(scale(t(biomarker_expr))),name = " ",
        col = col_fun,
        top_annotation = top_annotation,row_dend_side = 'right', row_names_side = 'left',
        column_split = Group,
        show_heatmap_legend = T,show_column_dend = F,
        show_column_names = F,
        show_row_names = T,row_names_gp = gpar(fontsize = 15),
        column_title = NULL)
dev.off()


#fig6c 
roc_pl <- ggplot(data=roc, aes(x=fpr, y=tpr))+
  geom_path(color="#2E5F6B", size=1)+
  geom_abline(intercept=0, slope=1, color="grey", size=1, linetype=2)+
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensivity or Recall)")+
  annotate("text", x=.75, y=.25, label=paste("AUC =", auc),
           size=5, family="serif")+
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(size=.5, color="black"),
        axis.title = element_text(color='black', size=13),
        axis.text = element_text(color='black', size=13),
        text = element_text(size=13, color="black"))

roc_pl
ggsave('H:/fig_project/roc.pdf',roc_pl,width = 4.27,height =3.69,units = 'in')


#fig6c 
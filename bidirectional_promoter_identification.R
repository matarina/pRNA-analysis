bidirection = read.csv("H:/prna_transcript_analysis/annotationnewplot/bidirectional_promoter/closed.text",header = FALSE,sep = '\t')
keep_bi = bidirection[abs(bidirection$V13)<= 1000,]
keep_bi$TSS_distance = -keep_bi$V13
tss_distange = ggplot(keep_bi, aes(TSS_distance)) +geom_histogram(color="black", fill="#91D8E4",) + theme_classic()#+ ylim(0, 400)+scale_x_continuous(limits = c(NA, 2000))
ggsave(filename = 'H:/prna_transcript_analysis/annotationnewplot/fig_project/tss_distance.pdf',plot = tss_distange,width = 3.27,height =2.69,units = 'in')


#get bidirecitonal bed file
bidirection = read.csv("H:/prna_transcript_analysis/annotationnewplot/bidirectional_promoter/closed.text",header = FALSE,sep = '\t')
bidirection  <- bidirection[bidirection$V8 > bidirection$V2, ]
bidirection = filter(bidirection,abs(bidirection$V13) <= 1000)
bidirection$name  = paste0(bidirection$V4,bidirection$V10)
bidirection$score = rep(1,nrow(bidirection))
bidirection$strand = rep('.',nrow(bidirection))
bi_promoter = bidirection[,c(1,2,9,14,15,16)]
write.table(bi_promoter,'H:/prna_transcript_analysis/annotationnewplot/bidirectional_promoter/bi_promoter.bed',sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)



#bi prna expression analysis in hepg2_2 cell line without stranded 
bi_exp = read.csv("H:/prna_transcript_analysis/annotationnewplot/bidirectional_promoter/bi_prna_count.txt",sep = '\t',skip = 1)
bi_exp$mean = (bi_exp$SRR8131644.sort.bam+bi_exp$SRR8131645.sort.bam)/2
library(dplyr)
bi_exp  = filter(bi_exp,mean>4)
bi_prna_masked  = read.csv("H:/bi_prna_masked.bed",header = FALSE,sep = '\t')

import pandas as pd
import plotly.express as px
# This dataframe has 244 lines, but 4 distinct values for `day`
data = {'group': ['expressed bidirectional promoter pairs', 'others'],
  'value': [351,1117]}
df = pd.DataFrame(data)
fig = px.pie(df, values='value', names='group')
fig['data'][0].update({'textinfo' : 'text+value+percent'})
fig.show()
fig.write_image("H:/prna_transcript_analysis/annotationnewplot/fig/bi_expr_pie.pdf")




map_id = read.csv('H:/prna_transcript_analysis/annotationnewplot/bidirectional_promoter/promoter_ensembl.txt',sep = '\t',header = FALSE)
maped = merge(bidirection,map_id,by.x = 'V4',by.y = 'V1')

maped2 = merge(maped,map_id,by.x = 'V10',by.y = 'V1')
enrich = c(maped2$V2,maped2$V2.y)
library(clusterProfiler)
library(org.Hs.eg.db)
ggo <- enrichGO(gene         = enrich,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

go = barplot(ggo)+theme_classic()
go
ggsave('H:/prna_transcript_analysis/annotationnewplot/fig/bi_go.pdf',plot = go,width = 4.27,height =3.69,units = 'in')


ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)





















bi_exp_prna_masked = bi_prna_masked[which(bi_exp$Geneid %in% bi_prna_masked$V4),]
write.table(bi_exp_prna_masked,"H:/biprna_expressed_masked.bed",sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
#179 expresed
#1307 feature
#84 cds overlap
#316 lnc overlap
frequency = read.csv("H:/closed.bed",header = FALSE,sep = '\t')
frequency$TSS_distance = abs(frequency$V13)



tss_distange = ggplot(frequency, aes(TSS_distance)) +geom_histogram(color="black", fill="#91D8E4",) + theme_classic()+ ylim(0, 400)+scale_x_continuous(limits = c(NA, 2000))
ggsave(filename = 'H:/fig_project/tss_distance.pdf',plot = tss_distange,width = 3.27,height =2.69,units = 'in')
WhatIHave <- c("expressed_pRNA_non_CpG" = 9,
               "expressed_pRNA_CpG" = 170)
data <- data.frame("group" = names(WhatIHave), WhatIHave)
vals <- paste(data$WhatIHave, sep = "")
library(plotly)
cpg_prna <- plot_ly(data, labels = ~group, values = ~WhatIHave, type = 'pie',textinfo = "text + values", text = vals, textfont = list( size = 15)) %>%
  layout(title = "",          
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),legend = list(font = list(size = 24)))
kaleido(cpg_prna, "H:/fig_project/cpg_prna.pdf", width = 4, height = 4)



df <- data.frame(
  group = c("expressed bi-pRNA","all bi-pRNA"),
  value = c(179,1128)
)
bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pie
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
cpg_prna = pie + scale_fill_manual(values = c('#91D8E4','#ECB390')) + blank_theme+
  theme(axis.text.x=element_blank()) +geom_text(aes(label = value), position = position_stack(vjust = 0.4),size=6)
 # geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
   #             label = percent(value/100)), size=5)
cpg_prna
ggsave(filename = 'H:/fig_project/cpg_prna.pdf',plot = cpg_prna,width = 6.27,height =2.69,units = 'in')

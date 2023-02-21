a = read.csv("H:/prna_transcript_analysis/epd_promoter/promoter_feature/expressed_ccaat.txt",sep = ' ',header = FALSE)
a = a[,c(2,4)]
colnames(a) = c("position","expressed_pRNA_CCAAT_box")
b = read.csv("H:/prna_transcript_analysis/epd_promoter/promoter_feature/unexpressed_ccaat.txt",sep = ' ',header = FALSE)
b = b[,c(2,4)]
colnames(b) = c("position","unexpressed_pRNA_CCAAT_box")
mer = merge(a,b,by="position")
library(ggplot2)
library(reshape2)
meltmer = na.omit(melt(mer,id.vars = "position"))
colnames(meltmer)[2] = "Group"
ggplot(meltmer, aes(x = position, y = value)) + 
  geom_line(aes(color = Group, linetype = Group)) + 
  scale_color_manual(values = c("#E25555", "steelblue"))+labs(x = "Position")+labs(y = "Frequency")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text=element_text(size=19),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text = element_text(size = 13),legend.position = "bottom",legend.title = element_text(size = 20),panel.background = element_blank())

library(tidyr)
a = read.table(text = gsub("\\s+", " ", readLines("H:/prna_transcript_analysis/epd_promoter/promoter_feature/expressed_ccaat.txt")))
colnames(a) = c("position","expressed_pRNA_CCAAT_box")
b = read.table(text = gsub("\\s+", " ", readLines("H:/prna_transcript_analysis/epd_promoter/promoter_feature/unexpressed_ccaat.txt")))
colnames(b) = c("position","unexpressed_pRNA_CCAAT_box")
mer = merge(a,b,by="position")
library(ggplot2)
library(reshape2)
meltmer = na.omit(melt(mer,id.vars = "position"))
colnames(meltmer)[2] = "Group"
ggplot(meltmer, aes(x = position, y = value)) + xlim(-1000, 100)+
  geom_line(aes(color = Group, linetype = Group)) + 
  scale_color_manual(values = c("#E25555", "steelblue"))+labs(x = "Position")+labs(y = "Frequency")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text=element_text(size=19),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text = element_text(size = 13),legend.position = "bottom",legend.title = element_text(size = 20),panel.background = element_blank())


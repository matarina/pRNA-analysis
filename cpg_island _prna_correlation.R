hepg2 = read.csv("H:/prna_transcript_analysis/new_prna/sense_hepg2_2.txt",sep = '\t',skip = 1)
h2h = read.csv("H:/prna_transcript_analysis/new_prna/unmask_promoter_with_cpg.bed",sep = '\t',header = FALSE)
library(dplyr)
hepg2$mean_expressed = (hepg2$SRR8131644.sort.bam+hepg2$SRR8131645.sort.bam)/2
hepg2 = filter(hepg2,mean_expressed>4)
hepg2  = hepg2[,c(1,10)]
hepg2 = hepg2[-grep("ERCC-",hepg2$Geneid),]
colnames(h2h)[4]="Geneid"
meg = merge(h2h,hepg2,by="Geneid")
cpg_promoter = read.csv("H:/cpg_island.txt",sep = '\t',header = FALSE)
unmask_mrna_promoter = read.csv("H:/prna_transcript_analysis/new_prna/mrna_promoter_unmask.bed",sep = '\t',header = FALSE)
unmask_prna_cpg = read.csv("H:/prna_transcript_analysis/new_prna/unmask_prna_with_cpg.bed",sep = '\t',header = FALSE)
unmask_prna_without_cpg = unmask_mrna_promoter[-which(unmask_mrna_promoter$V4 %in% unmask_prna_cpg$V4),] 
colnames(unmask_prna_cpg)[4]="Geneid"
unmask_prna_cpg_expressed = merge(hepg2,unmask_prna_cpg,by="Geneid")
unmask_prna_without_cpg_expressed = merge(hepg2,unmask_prna_without_cpg,by="Geneid")
unexpressed_pRNA = unmask_mrna_promoter[-which(unmask_mrna_promoter$V4 %in% hepg2$Geneid),]
table(unexpressed_pRNA$V4 %in% unmask_prna_cpg$Geneid)
library(tidyr)
separate(data = cpg_promoter, col = V6, into = c("left", "cpg_length"), sep = "\\|")
cpg_promoter = separate(data = cpg_promoter, col = V6, into = c("left", "cpg_length"), sep = "\\|")
colnames(cpg_promoter)[4] = "Geneid"
merge = merge(hepg2,cpg_promoter,by="Geneid")
library(ggplot2)
ggplot(merge, aes(x=SRR8131644.sort.bam, y=cpg_length)) +
  geom_point(size=2, shape=23)
cpg_promoter$cpg_length = as.numeric(cpg_promoter$cpg_length)
hepg22 = read.csv("H:/single_cell_sepcific_tpm.csv",sep = ',',skip = 1,header = FALSE)
colnames(hepg22)=c("Geneid","FPKM")
merge2 = merge(hepg22,cpg_promoter,by="Geneid")
filtermerge2 = filter(filtermerge2,FPKM<5000)
## stacked barplot
library(ggplot2)
library(tidyverse)
library(grid)
# create a dataset
specie <- c(rep("total_promoter" , 2) , rep("expressed_pRNA" , 2) , rep("unexpressed_pRNA" , 2) )
condition <- rep(c("CpG" , "without_CpG" ) , 3)
value <- c(11949,4506,1215,189,10734,4317)
data <- data.frame(specie,condition,value)
library(cartography)
mypal <- carto.pal(pal1 = "pastel.pal", n1 = 13,middle = TRUE, transparency = FALSE)
# Stacked
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = mypal)+
  geom_text(aes(label = value), size = 5, hjust = 0.5, vjust = 2,color="white",position =     "fill")+
  labs(x = "Group")+labs(y = "Percentage")+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size = 20),panel.background = element_blank())+
  theme_classic() +
  theme(axis.text = element_text(size = 15))     
        

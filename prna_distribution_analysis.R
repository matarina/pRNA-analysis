pseudogene_annotation = read.csv("H:/prna_transcript_analysis/prna_count/subdivision_genecode_annotation/pseudogene_annotation.csv",header = FALSE)
mrna_annotation = read.csv("H:/prna_transcript_analysis/prna_count/subdivision_genecode_annotation/mrna_annotation.csv",header = FALSE)
nodecay_annotation = read.csv("H:/prna_transcript_analysis/prna_count/subdivision_genecode_annotation/nonsense_mediated_decay_annotation.csv",header = FALSE)
lncrna_annotation = read.csv("H:/prna_transcript_analysis/prna_count/subdivision_genecode_annotation/lncrna_annotation.csv",header = FALSE)
processed_transcripts_annotation = read.csv("H:/prna_transcript_analysis/prna_count/subdivision_genecode_annotation/processed_transcript_annotation.csv",header = FALSE)
meripid = read.csv("H:/prna_transcript_analysis/prna_count/m6a_prna_distribution/drach_true_cut_uniq.csv",header = FALSE)

table(meripid$V1 %in% lncrna_annotation$V1)
table(meripid$V1 %in% mrna_annotation$V1)
table(meripid$V1 %in% pseudogene_annotation$V1)
table(meripid$V1 %in% processed_transcripts_annotation$V1)
table(meripid$V1 %in% nodecay_annotation$V1)







table(is.na(X293T_rep1_mrna$V1))

length(unique(X293T_rep1_mrna$V8))
p <- ggplot(data = staked_barplot, aes(x = cell_line, y = count_number)) +geom_col(aes(fill = group), width = 0.7)+
  geom_text(aes( label = count_number, group =group), color = "white",vjust=2.3)
p
intersect_Hela_mrna = intersect(Hela_rep1_mrna,Hela_rep2_mrna)
intersect_HepG2_mrna = intersect(HepG2_rep1_mrna,HepG2_rep2_mrna)
intersect_RPE_mrna = intersect(RPE_rep1_mrna,RPE_rep2_mrna)

x293Trep1 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/293T_0h_1.tab",sep="\t",header = TRUE)
x293Trep2 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/293T_0h_2.tab",sep="\t",header = TRUE)
helarep1 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/Hela_0h_1.tab",sep="\t",header = TRUE)
helarep2 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/Hela_0h_2.tab",sep="\t",header = TRUE)
hepg2rep1 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/HepG2_0h_rep1.tab",sep="\t",header = TRUE)
hepg2rep2 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/HepG2_0h_rep2.tab",sep="\t",header = TRUE)
rperep1 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/RPE_0h_A.tab",sep="\t",header = TRUE)
rperep2 = read.table("H:/prna_transcript_analysis/prna_count/0h_prna_expression_tab/RPE_0h_B.tab",sep="\t",header = TRUE)
x293Trep1 = x293Trep1[x293Trep1$TPM> 1,]
x293Trep2 = x293Trep2[x293Trep2$TPM> 1,]
helarep1 = helarep1[helarep1$TPM> 1,]
helarep2 = helarep2[helarep2$TPM> 1,]
hepg2rep1 = hepg2rep1[hepg2rep1$TPM> 1,]
hepg2rep2 = hepg2rep2[hepg2rep2$TPM> 1,]
rperep1 = rperep1[rperep1$TPM> 1,]
rperep2 = rperep2[rperep2$TPM> 1,]
x293Trep1_mrna = na.omit(pseu_annotation[match(x293Trep1$Gene.ID,pseu_annotation$V1),])
length(x293Trep1_mrna)
length(unique(x293Trep1_mrna))
x293Trep2_mrna = na.omit(mrna_annotation[match(x293Trep2$Gene.ID,mrna_annotation$V4),])
length(x293Trep2_mrna$V8)
length(unique(x293Trep2_mrna$V8))
helarep1_mrna = na.omit(mrna_annotation[match(helarep1$Gene.ID,mrna_annotation$V4),])
length(helarep1_mrna$V8)
length(unique(helarep1_mrna$V8)) 
helarep2_mrna = na.omit(mrna_annotation[match(helarep2$Gene.ID,mrna_annotation$V4),])
length(helarep2_mrna$V8)
length(unique(helarep2_mrna$V8)) 
hepg2rep1_mrna = na.omit(mrna_annotation[match(hepg2rep1$Gene.ID,mrna_annotation$V4),])
length(hepg2rep1_mrna$V8)
length(unique(hepg2rep1_mrna$V8)) 
hepg2rep2_mrna = na.omit(mrna_annotation[match(hepg2rep2$Gene.ID,mrna_annotation$V4),])
length(hepg2rep2_mrna$V8)
length(unique(hepg2rep2_mrna$V8)) 
rperep1_mrna = na.omit(mrna_annotation[match(rperep1$Gene.ID,mrna_annotation$V4),])
length(rperep1_mrna$V8)
length(unique(rperep1_mrna$V8)) 
rperep2_mrna = na.omit(mrna_annotation[match(rperep2$Gene.ID,mrna_annotation$V4),])
length(rperep2_mrna$V8)
length(unique(rperep2_mrna$V8)) 


x293Trep1_lncrna = na.omit(lncrna_annotation[match(x293Trep1$Gene.ID,lncrna_annotation$V4),])
length(x293Trep1_lncrna$V8)
length(unique(x293Trep1_lncrna$V8))
x293Trep2_lncrna = na.omit(lncrna_annotation[match(x293Trep2$Gene.ID,lncrna_annotation$V4),])
length(x293Trep2_lncrna$V8)
length(unique(x293Trep2_lncrna$V8))
helarep1_lncrna = na.omit(lncrna_annotation[match(helarep1$Gene.ID,lncrna_annotation$V4),])
length(helarep1_lncrna$V8)
length(unique(helarep1_lncrna$V8)) 
helarep2_lncrna = na.omit(lncrna_annotation[match(helarep2$Gene.ID,lncrna_annotation$V4),])
length(helarep2_lncrna$V8)
length(unique(helarep2_lncrna$V8)) 
hepg2rep1_lncrna = na.omit(lncrna_annotation[match(hepg2rep1$Gene.ID,lncrna_annotation$V4),])
length(hepg2rep1_lncrna$V8)
length(unique(hepg2rep1_lncrna$V8)) 
hepg2rep2_lncrna = na.omit(lncrna_annotation[match(hepg2rep2$Gene.ID,lncrna_annotation$V4),])
length(hepg2rep2_lncrna$V8)
length(unique(hepg2rep2_lncrna$V8)) 
rperep1_lncrna = na.omit(lncrna_annotation[match(rperep1$Gene.ID,lncrna_annotation$V4),])
length(rperep1_lncrna$V8)
length(unique(rperep1_lncrna$V8)) 
rperep2_lncrna = na.omit(lncrna_annotation[match(rperep2$Gene.ID,lncrna_annotation$V4),])
length(rperep2_lncrna$V8)
length(unique(rperep2_lncrna$V8))


stacked_barplot = read.csv("H:/prna_transcript_analysis/prna_count/stacked_barplot.csv")
library(ggplot2)
p <- ggplot(data = stacked_barplot, aes(x = cell_line, y = count_number)) +geom_col(aes(fill = group), width = 0.7)+geom_text(aes( label = count_number, group =group), color = "white", position = position_stack(vjust = 0.5))
p
p


#甲基化分布
transcripts_peak = read.csv("H:/merip_prna/extractid.csv")
mrna_peak = na.omit(mrna_annotation[match(transcripts_peak$name,mrna_annotation$V4),])
length(mrna_peak$V1)
lncrna_peak = na.omit(lncrna_annotation[match(transcripts_peak$name,lncrna_annotation$V4),])
length(lncrna_peak$v1)
library(ggplot2)
stacked_barplot = read.csv("H:/merip_prna/merip_stack_column.csv")


##pseugene分布
x293Trep1_mrna = na.omit(pseu_annotation[match(x293Trep1$Gene.ID,pseu_annotation$V1),])
length(x293Trep1_mrna)
length(unique(x293Trep1_mrna))
x293Trep2_mrna = na.omit(pseu_annotation[match(x293Trep2$Gene.ID,pseu_annotation$V1),])
length(x293Trep2_mrna)
length(unique(x293Trep2_mrna))
helarep1_mrna = na.omit(pseu_annotation[match(helarep1$Gene.ID,pseu_annotation$V1),])
length(helarep1_mrna)
length(unique(helarep1_mrna))
helarep2_mrna = na.omit(pseu_annotation[match(helarep2$Gene.ID,pseu_annotation$V1),])
length(helarep2_mrna)
length(unique(helarep2_mrna))
hepg2rep1_mrna = na.omit(pseu_annotation[match(hepg2rep1$Gene.ID,pseu_annotation$V1),])
length(hepg2rep1_mrna)
length(unique(hepg2rep1_mrna))
hepg2rep2_mrna = na.omit(pseu_annotation[match(hepg2rep2$Gene.ID,pseu_annotation$V1),])
length(hepg2rep2_mrna)
length(unique(hepg2rep2_mrna))
rperep1_mrna = na.omit(pseu_annotation[match(rperep1$Gene.ID,pseu_annotation$V1),])
length(rperep1_mrna)
length(unique(rperep1_mrna))
rperep2_mrna = na.omit(pseu_annotation[match(rperep2$Gene.ID,pseu_annotation$V1),])
length(rperep2_mrna)
length(unique(rperep2_mrna))      


##processed_transcript
process_annotation = read.table("H:/prna_transcript_analysis/prna_count/nonsense_anno.txt")
x293Trep1_mrna = na.omit(process_annotation[match(x293Trep1$Gene.ID,process_annotation$V1),])
length(x293Trep1_mrna)
length(unique(x293Trep1_mrna))
x293Trep2_mrna = na.omit(process_annotation[match(x293Trep2$Gene.ID,process_annotation$V1),])
length(x293Trep2_mrna)
length(unique(x293Trep2_mrna))
helarep1_mrna = na.omit(process_annotation[match(helarep1$Gene.ID,process_annotation$V1),])
length(helarep1_mrna)
length(unique(helarep1_mrna))
helarep2_mrna = na.omit(process_annotation[match(helarep2$Gene.ID,process_annotation$V1),])
length(helarep2_mrna)
length(unique(helarep2_mrna))
hepg2rep1_mrna = na.omit(process_annotation[match(hepg2rep1$Gene.ID,process_annotation$V1),])
length(hepg2rep1_mrna)
length(unique(hepg2rep1_mrna))
hepg2rep2_mrna = na.omit(process_annotation[match(hepg2rep2$Gene.ID,process_annotation$V1),])
length(hepg2rep2_mrna)
length(unique(hepg2rep2_mrna))
rperep1_mrna = na.omit(process_annotation[match(rperep1$Gene.ID,process_annotation$V1),])
length(rperep1_mrna)
length(unique(rperep1_mrna))
rperep2_mrna = na.omit(process_annotation[match(rperep2$Gene.ID,process_annotation$V1),])
length(rperep2_mrna)
length(unique(rperep2_mrna))      


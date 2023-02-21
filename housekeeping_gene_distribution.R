setwd("H:/")
SRR5062105_prna = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR5062105_prna.csv",sep = '\t',skip = 1)
SRR5062106_prna = read.csv("H:/prna_transcript_analysis/epd_promoter/mrna_relation/SRR5062106_prna.csv",sep = '\t' ,skip = 1)
promoter_ensembl = read.csv("H:/prna_transcript_analysis/epd_promoter/epd_annotaion/promoter_ensembl.txt",sep = '\t',header = FALSE)
colnames(promoter_ensembl) = c("Geneid","ensmbl")

SRR5062105_prna$average  = (SRR5062105_prna$SRR5062105.sorted.bam+SRR5062106_prna$SRR5062106.sorted.bam)/2


###reshape and clean data to linux rnadecay package
output = read.csv("H:/pas_cpm.txt",sep = '\t',row.names = 1)
library(dplyr)

output =  dplyr::select(output, -c('SRR8131651.sort.bam'))
output$SRR8131651.sort.bam = output$SRR8131650.sort.bam
rpms = output
### latest shape data 
#rpms = read.csv("H:/prna_transcript_analysis/new_prna/antisense_hepg2_2.txt",sep = '\t',skip = 1,row.names = 1)
colnames(rpms) = c("hepg2_00_r1","hepg2_00_r2","hepg2_60_r1","hepg2_60_r2","hepg2_180_r1","hepg2_180_r2","hepg2_360_r1","hepg2_360_r2")
rpms = filter(rpms,hepg2_00_r1>0,hepg2_00_r2>0)
#rpms[] = lapply(rpms, function(x) x*10^6/sum(x))
sham = rpms
colnames(sham)=c("sham_00_r1","sham_00_r2","sham_60_r1","sham_60_r2","sham_180_r1","sham_180_r2","sham_360_r1","sham_360_r2")
hepg2_linux = as.data.frame(cbind(rpms,sham))
rm (sham)
write.csv(na.omit(hepg2_linux),"H:/pas_hepg2_linux.txt")




rm(list=ls())
library(ggplot2)
library(tidyverse)


### 
### compare IGF2BP1 target pRNAs eClip data CDF plot

bp3target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/pas_bp3_merip.bed",header = FALSE,sep = "\t")
decaydata = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/pas_result.txt")
decaydata = decaydata[-grep("ERCC-",decaydata$X),]
decaydata = filter(decaydata,sigma2<0.066)
bp1 = decaydata[(decaydata$X %in% bp3target$V4),]
bp1$halflife = log(2)/bp1$alpha_hepg2
bp1$group = rep("target \nmean :303.17min",198)
without_bp1 = decaydata[!(decaydata$X %in% bp3target$V4),]###删除在另一个表中含有同名的行
without_bp1$halflife = log(2)/without_bp1$alpha_hepg2
without_bp1$group = rep("untarget \nmean :242.05min",1379)
bp1andwithout = bind_rows(bp1,without_bp1)
ks.test(bp1$halflife,without_bp1$halflife, alternative = "l")

p = ggplot(bp1andwithout,
       aes(
           x=halflife,
         color = group
       ))+
  stat_ecdf( # ggplot2中的经验累积分布函数
    size=1   # 线条粗细
  )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
          legend.title=element_text(size=18),legend.position = c(0.6, 0.4),
          legend.text=element_text(size=18),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          annotate("text", x = 340.05, y = 0.073,size = 5, label = "KS-test P-value:3.192e-10")+
          scale_x_continuous(expand = c(0, 5),limits = c(0,495)) + scale_y_continuous(expand = c(0, 0))+#+ theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
          theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text=element_text(size=14),
              axis.title=element_text(size=18,face="bold"))+labs(y = "Cumulative fraction")+labs(x = " halflife")+theme(plot.margin=unit(c(2, 2, 2, 2),'cm'))
p#761 705










### compare mrna and psRNA pasRNA halflife


###source from RNA_decay_master.R object 
mean(pas_filtered$mean_nls,trim = 0.2)
pas_filtered$min = pas_filtered$mean_nls*60
pas_filtered$group = rep("pasRNA\nmean :252.20min",986)

mean(ps_filtered$mean_nls,trim = 0.2)
ps_filtered$group = rep("psRNA\nmean :241.34min",787)
ps_filtered$min = ps_filtered$mean_nls*60

mean(mrna_filtered$mean_nls,trim = 0.2)
mrna_filtered$group = rep("mRNA\nmean :530.61min",6871)
mrna_filtered$min = mrna_filtered$mean_nls*60
ks.test(mrna_filtered$min,pas_filtered$min, alternative = "l")
ks.test(mrna_filtered$min,ps_filtered$min, alternative = "l")
mrna_vs_ps = rbind(mrna_filtered,pas_filtered)
mrna_vs = rbind(mrna_vs_ps,ps_filtered)
library(ggplot2)
p = ggplot(mrna_vs,
           aes(
             x=min,
             color = group
           ))+
  stat_ecdf( # ggplot2中的经验累积分布函数
    size=1   # 线条粗细
  )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
          legend.title=element_text(size=18),legend.position = c(0.75, 0.4),
          legend.text=element_text(size=18),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 450.05, y = 0.073,size = 5, label = "KS-test P-value:2.2e-16")
p  = p+ scale_x_continuous(limits = c(0,600))
p =  p + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text=element_text(size=14),
             axis.title=element_text(size=18,face="bold"))+labs(y = "Cumulative fraction")+labs(x = " halflife")
p#682 537



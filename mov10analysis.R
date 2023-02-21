## ggplot box 
library(ggplot2)
library(ggpubr)
library(edgeR)
library(tidyverse)
mapid = read.csv('H:/prna_transcript_analysis/bp_target_decay_cdf/map_id_type.txt',sep = '\t',skip = 4)
expr = read.csv('H:/mov10_decay/genes_count.csv',row.names = 1)

expr = data.frame(edgeR::cpm(expr))
expr[grep('ENSG00000155363',rownames(expr)),]

p = data.frame(group = c('nc','nc','sh','sh'),RPM = as.numeric(c('60.77505','51.58195','13.83542','15.22056')))
ggplot(p, aes(x=as.factor(group), y=RPM,fill=group)) + stat_boxplot(geom = "errorbar", width = 0.15)+ theme(axis.text=element_text(size=20),
                                                                                                            axis.title=element_text(size=20),
                                                                                                            panel.grid.major = element_blank(), 
                                                                                                            panel.grid.minor = element_blank(),
                                                                                                            panel.background = element_blank(),
                                                                                                            axis.line = element_line(colour = "black"))+
  geom_boxplot(show.legend = FALSE) + scale_fill_manual(values=c("#AE3E3E", "#90FCFF"))+xlab('')

##mRNA and lncRNA lifetime determine 
decay  = expr[,c(10,5,7,6,9,8,16,11,13,12,15,14)]
colnames(decay) = c('nc_0h_r1','nc_0h_r2','nc_2h_r1','nc_2h_r2','nc_4h_r1','nc_4h_r2','sh_0h_r1','sh_0h_r2','sh_2h_r1','sh_2h_r2','sh_4h_r1','sh_4h_r2')
mrnadecay = decay[which(rownames(decay) %in% distinct(mapid[which(mapid[,2] == 'protein_coding'),] )[,1]),]

for (i in 1:nrow(mrnadecay)){
  mrnadecay$nc1_hl[i] = 2/(log2(mrnadecay$nc_2h_r1[i]/mrnadecay$nc_0h_r1[i])/-120+log2(mrnadecay$nc_4h_r1[i]/mrnadecay$nc_0h_r1[i])/-360)
  mrnadecay$nc2_hl[i] = 2/(log2(mrnadecay$nc_2h_r2[i]/mrnadecay$nc_0h_r2[i])/-120+log2(mrnadecay$nc_4h_r2[i]/mrnadecay$nc_0h_r2[i])/-360)
  mrnadecay$sh1_hl[i] = 2/(log2(mrnadecay$sh_2h_r1[i]/mrnadecay$sh_0h_r1[i])/-120+log2(mrnadecay$sh_4h_r1[i]/mrnadecay$sh_0h_r1[i])/-360)
  mrnadecay$sh2_hl[i] = 2/(log2(mrnadecay$sh_2h_r2[i]/mrnadecay$sh_0h_r2[i])/-120+log2(mrnadecay$sh_4h_r2[i]/mrnadecay$sh_0h_r2[i])/-360)
  }
halflifetab_2 = mrnadecay2

for (i in 1:nrow(halflifetab_2)){
  halflifetab_2$foldchange[i]<-log2(mean(halflifetab_2$sh1_hl[i]+halflifetab_2$sh2_hl[i])/
                                      mean(halflifetab_2$nc1_hl[i]+halflifetab_2$nc2_hl[i]))
  halflifetab_2$pvalue[i]<-t.test(c(halflifetab_2$sh1_hl[i],halflifetab_2$sh2_hl[i]),
                                  c(halflifetab_2$nc1_hl[i],halflifetab_2$nc2_hl[i]))$p.value
}

library(dplyr)
library(limma)
library(ggrepel)

mrnadecay2  = mrnadecay[,c((ncol(mrnadecay)-3):ncol(mrnadecay))]
rm <- apply(mrnadecay2, 1, function(x){
  sum(x == 0) >= 2
})
mrnadecay2 <- na.omit(mrnadecay2[!rm,])
mrnadecay2 = mrnadecay2 %>% filter(nc1_hl>0,nc1_hl>0,sh2_hl>0,sh1_hl>0)
mrnadecay2 <- na.omit(log2(mrnadecay2 + 1)) 

group_list <- factor(c(rep("A",2), rep("B",2)))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(mrnadecay2)
fit <- lmFit(mrnadecay2, design)
fit <- eBayes(fit, trend=TRUE)
result_limma <- topTable(fit, coef=2,n=Inf)
sum(result_limma$P.Value < 0.05 & result_limma$logFC > -1)
genes = result_limma 
genes = rownames_to_column(genes,var = 'Gene')
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes$Gene <- sub("[.][0-9]*","",genes$Gene)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),values=genes$Gene,mart= mart)
genes = merge(genes,G_list,by.x="Gene",by.y="ensembl_gene_id")
genes = halflifetab_2
colnames(genes)[5:6] = c('logFC','P.Value')
genes$logFC = log2(genes$logFC)
genes$Significant <- ifelse(genes$P.Value < 0.05 & abs(genes$logFC) >= 1, 
                            ifelse(genes$logFC > 1, "Up", "Down"), "Stable")
gg = ggplot(
  # 数据、映射、颜色
  genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("#354259","#C2DED1", "#CDC2AE")) +
  # 注释
 # geom_text_repel(
  #  data = subset(genes, P.Value < 0.05 & abs(genes$logFC) >= 1),
  #  aes(label = ),
  #  size = 5,
  #  box.padding = unit(0.35, "lines"),
  #  point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # 图例
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.background = element_blank(),
        axis.line= element_line(colour = "black"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.4))
gg
ggsave("plot3.png",gg,width = 8, height = 8, dpi = 300 )
























coldata = data.frame(row.names = colnames(mrnadecay),condition =  c('nc_0h','nc_0h','nc_2h','nc_2h','nc_4h','nc_4h','sh_0h','sh_0h','sh_2h','sh_2h','sh_4h','sh_4h'))
all(colnames(mrnadecay) == rownames(coldata))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=mrnadecay, colData=coldata, design= ~ condition)


keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
#if error in antisense
#vsd = varianceStabilizingTransformation(dds,blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
pca = plotPCA(vsd, intgroup=c("condition"),returnData = TRUE)
ggplot(pca, aes(x=PC1, y=PC2,color=group))+ geom_point(size=2.5)+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                  panel.grid.minor=element_blank(),
                                                                                  axis.line= element_line(colour = "black"),
                                                                                  legend.title = element_text(size=15),
                                                                                  legend.text = element_text(size=15),
                                                                                  axis.title.x = element_text(size = 15),
                                                                                  axis.title.y = element_text(size = 16),
                                                                                  plot.title = element_text(hjust = 0.4))


###reshape and clean data to linux rnadecay package
expr = read.csv('H:/output.txt',skip = 1,sep = '\t',row.names = 1)
expr = data.frame(edgeR::cpm(expr))
expr  = expr[,c(10,5,7,6,9,8,16,11,13,12,15,14)]
### latest shape data 
#rpms = read.csv("H:/prna_transcript_analysis/new_prna/antisense_hepg2_2.txt",sep = '\t',skip = 1,row.names = 1)
colnames(expr) = c("nc_00_r1","nc_00_r2","nc_120_r1","nc_120_r2","nc_240_r1","nc_240_r2","sh_00_r1","sh_00_r2","sh_120_r1","sh_120_r2","sh_240_r1","sh_240_r2")
library(dplyr)
rpms = expr 
rpms = filter(rpms,nc_00_r1>0,nc_00_r2>0,sh_00_r1>0,sh_00_r2>0,)
#rpms[] = lapply(rpms, function(x) x*10^6/sum(x))

write.csv(na.omit(rpms),"H:/decay_linux.txt")




###cdf decay plot
result = read.csv('H:/results.txt',sep = ',')
decaydata = filter(result,sigma2<0.066)
decaydata = result
decaydata$nc_hl = log(2)/decaydata$alpha_nc
decaydata$sh_hl = log(2)/decaydata$alpha_sh
#pdata = data.frame(group = rep(c('nc','sh'),each = 7343),halflife = c(decaydata$nc_hl,decaydata$sh_hl))
colnames(decaydata)[6] = 'pvalue'
decaydata$log2FoldChange = as.numeric(foldchange(decaydata$sh_hl,decaydata$nc_hl))
decay = decaydata
decaydata = decay
decaydata$xx = -log10(decaydata$pvalue)
for (i in 1:nrow(decaydata)){
  if (decaydata[i,16] > 0){decaydata[i,16] = log2(decaydata[i,16])}
  else
  {decaydata[i,16] = -log2(abs(decaydata[i,16]))}
}


EnhancedVolcano(decaydata,
                lab =NA, legendPosition = 'right',
                x = 'log2FoldChange',
                y = 'pvalue',parseLabels = FALSE,caption = 'log2FC cutoff: |2|; p-value cutoff: 0.05',
                pCutoff = 0.065)


library(tidyverse)
library(dplyr)
nc1_peak = read.csv('H:/nc1_peaks.csv')
nc2_peak = read.csv('H:/nc2_peaks.csv')
mov10_clip  =read.csv('H:/mov10_clip.csv')
nc1_peak = nc1_peak %>% distinct(geneID,.keep_all = TRUE) %>% filter(log2FC>2,pvalue<0.05,fdr<0.05)
nc2_peak = nc2_peak %>% distinct(geneID,.keep_all = TRUE) %>% filter(log2FC>2,pvalue<0.05,fdr<0.05)
m6a_mov10  = intersect(mov10_clip$Target.gene.ID , gsub('\\..*','',intersect(nc1_peak$geneID ,nc2_peak$geneID)))
result = read.csv('H:/limit_result.txt',sep = ',')
decaydata = filter(result,sigma2<0.066)
decaydata$nc_hl = log(2)/decaydata$alpha_nc
decaydata$sh_hl = log(2)/decaydata$alpha_sh

decaydata$ensembl = gsub('\\..*','',decaydata$X)

target = decaydata[(decaydata$ensembl %in% m6a_mov10),]
target$group = rep('target',1517)
mean(untarget$nc_hl)
untarget = decaydata[!(decaydata$ensembl %in% m6a_mov10),]
untarget$group = rep("untarget",5826)
pdata = rbind(target,untarget)
pdata$log_hf = log(pdata$nc_hl)
ks.test(untarget$nc_hl,target$nc_hl, alternative = "l")
p = ggplot(pdata,
           aes(
             x=nc_hl,
             color = group
           ))+
  stat_ecdf( # ggplot2中的经验累积分布函数
    size=1   # 线条粗细
  )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),##空白背景
          legend.title=element_text(size=18),legend.position = c(0.6, 0.4),
          legend.text=element_text(size=18),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #annotate("text", x = 200.05, y = 0.073,size = 5, label = "KS-test P-value:3.758e-07")+
  scale_x_continuous(expand = c(0, 5),limits = c(0,489)) +
  scale_y_continuous(expand = c(0, 0))+#+ theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"))+labs(y = "Cumulative fraction")+labs(x = " halflife")+theme(plot.margin=unit(c(2, 2, 2, 2),'cm'))
p#7


library(gtools)
#volcanoplot 

results = t.test(x, y)
decaydata$pvalue = t.test(decaydata$nc_hl,decaydata$sh_hl)
library(EnhancedVolcano)
 EnhancedVolcano(decaydata,
                  lab =NA, legendPosition = 'right',
                  x = 'log2FoldChange',
                  y = 'padj',parseLabels = FALSE,caption = 'log2FC cutoff: |2|; p-value cutoff: 0.05',
                  pCutoff = 0.05)
 
 
 library(RCurl)
 
 library(ggpubr)
 decaydata = decaydata[order(decaydata$foldchange),]
 up = decaydata[c(1:20,7324:7343),]
 up$ensembl = gsub('\\..*$','',up$X)
 library('biomaRt')
 mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = up$ensembl
 G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","ensembl_gene_id_version", "hgnc_symbol"),values=genes,mart= mart)
 G_list[c(20,30,40),3] = c('pseudogene','pseudogene2','uncategorized_gene')
 colnames(G_list)[1] = 'ensembl'
 merge = merge(up,G_list,by= 'ensembl')
merge$hgnc_symbol = as.factor(merge$hgnc_symbol)
 ggdotchart(merge, x = "hgnc_symbol", y = "log2Foldchange",
            color = "red",                            
            sorting = "ascending",                        
            add = "segments",                             
            xlab="", ylab = 'fold_change',
            rotate = TRUE,
            group = "log2Foldchange", 

            
 )
 ggdotchart(merge, x = "hgnc_symbol", y = "log2Foldchange",position = p
            color = "hgnc_symbol", #点得颜色分组依据
            palette = "npg", rotate = TRUE,
            sorting = "ascending",#根据y轴从大到小排序
            add = "segments")#添加棒棒糖的棍子
class(merge$log2Foldchange)




###guitar
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
stBedFiles <- list(('H:/nc1_peaks.bed'),('H:/sh1_peaks.bed'))
# Build Guitar Coordinates
txdb_file <- system.file("extdata", "mm10_toy.sqlite",
                         package="Guitar")
txdb <- loadDb(txdb_file)
# Guitar Plot
GuitarPlot(txTxdb = txdb,
           stBedFiles = stBedFiles,
           headOrtail = TRUE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = c("nc1","sh1"))           


## venn plot
library(VennDiagram)

# Generate 3 sets of 200 words
set1 <-  gsub('\\..*','',intersect(nc1_peak$geneID ,nc2_peak$geneID))
set2 <- mov10_clip$Target.gene.ID
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)

venn.diagram(
  x = list(set1, set2),
  category.names = c( "m6A genes","MOV10 Clip genes" ),
  imagetype="png" ,
  height = 1980 , 
  width = 1980 , 
  resolution = 300,
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  lwd = 4,
  lty = 'solid',
  col=c("#440154ff", '#21908dff'),
  cat.pos = c(-27, 27)
)










###m6a metaplotR density profile 
m6a.dist <- read.delim ("H:/m6adisk/m6a.dist.measures.txt", header = T)
# assign the regions to new dataframes
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
# rescale 5'UTR and 3'UTR
library("scales")
library(plotly)
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
# Combine and plot
## Histogram
m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
############
m6a.dist <- read.delim ("H:/chrome downlaod/sh/sh_m6a.dist.measures.txt", header = T)
# assign the regions to new dataframes
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
# rescale 5'UTR and 3'UTR
library("scales")
library(plotly)
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
# Combine and plot
## Histogram
sh_m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)



metagene.cord <- c(m6a.metagene.coord, sh_m6a.metagene.coord)
group <- c(rep("nc", length(m6a.metagene.coord)), 
         rep("sh", length(sh_m6a.metagene.coord))) 
df <- data.frame(metagene.cord, group)
gg= ggplot(df) + geom_density(aes(x = metagene.cord, color = group),adjust = 4, size = 1.1) + xlim(0, 3) +
  geom_vline(xintercept = 1:2, col = "grey")+
  theme(axis.text.x=element_blank(), 
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
   labs(x = "5UTR              CDS                 3UTR")+scale_color_manual(values=c("#495579", "#CE7777"))
 


 ggsave("plot.png",gg,width = 8, height = 8, dpi = 300 )






p <- qplot(m6a.metagene.coord, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()
qplot(m6a.metagene.coord, geom="freqpoly") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
qplot(m6a.metagene.coord, geom="density") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
df = data.frame(m6a.metagene.coord)
write.csv(df,'H:/df.csv')
gg = ggplot(df) + geom_density(aes(x = m6a.metagene.coord),adjust = 4, colour = mod,size = 0.7, outline.type = "upper") + xlim(0, 3) + theme(axis.text.x=element_blank(),
                                                                                       axis.title=element_text(size=20),
                                                                                       panel.grid.major = element_blank(), 
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank(),
                                                                                       axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = 1:2, col = "grey")
gg
ggplotly(gg,source = "B")




nc <- t(read.delim ("H:/m6adisk/m6a.dist.measures.txt", header = T))
sh <- t(read.delim ("H:/chrome downlaod/sh/sh_m6a.dist.measures.txt", header = T))


metagene.cord <- c(nc, sh)
mod <- c(rep("nc", length(nc)), 
         rep("sh", length(sh))) 
df <- data.frame(metagene.cord, mod)

ggplot(df) + geom_density(aes(x = metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")




































##### ###########RNADECAY MASTER 


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


























###
###
### compare IGFBPbp1 IGFBPbp3 clip target pRNA halflife data 



#calculate psRNA/pasRNA bp1 clip target halflife
library(dplyr)
setwd("H:/prna_transcript_analysis/bp_target_decay_cdf/")
pas = read.csv("H:/mov10_decay/output.txt",sep = '\t',skip = 1)


bp1_rep1 = pas[,c(1,6,7,9)]
colnames(bp1_rep1) = c("Geneid","t0","t2","t4")


library(edgeR)
cm = bp1_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
nc1_hl = obj

nc2 = pas[,c(1,11,8,10)]
colnames(nc2) = c("Geneid","t0","t2","t4")
cm = nc2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
nc2_hl = obj


sh1 = pas[,c(1,17,14,16)]
colnames(sh1) = c("Geneid","t0","t2","t4")
cm = sh1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
sh1_hl = obj


sh2 = pas[,c(1,12,13,15)]
colnames(sh2) = c("Geneid","t0","t2","t4")
cm = sh2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
sh2_hl = obj


## preprocess ,remove NA value
library(ggplot2)
library(dplyr)
library(tidyverse)
colnames(nc1_hl)[3] = "nc1_hl"
colnames(nc2_hl)[3] = "nc2_hl"
colnames(sh1_hl)[3] = "sh1_hl"
colnames(sh2_hl)[3] = "sh2_hl"
df_list = list(nc1_hl,nc2_hl,sh1_hl,sh2_hl)
data = Reduce(function(x, y) merge(x, y, by="genes"), df_list)
data = data[,c(1,3,5,7,9)]
data[is.na(data)] = 0
data = column_to_rownames(data,var = 'genes')

rm <- apply(data, 1, function(x){
  sum(x == 0) >= 2
})
df <- data[!rm,]

library(limma)
library(ggrepel)
df1 <- na.omit(log2(df + 1))
group_list <- factor(c(rep("A",2), rep("B",2)))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(df1)

fit <- lmFit(df1, design)
fit <- eBayes(fit, trend=TRUE)
result_limma <- topTable(fit, coef=2,n=Inf)
sum(result_limma$P.Value < 0.05 & result_limma$logFC < -1)
genes = result_limma 
genes = rownames_to_column(genes,var = 'Gene')

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
x <- genes$Gene
genes$Gene <- sub("[.][0-9]*","",genes$Gene)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),values=genes$Gene,mart= mart)
genes = merge(genes,G_list,by.x="Gene",by.y="ensembl_gene_id")









genes$Significant <- ifelse(genes$P.Value < 0.05 & abs(genes$logFC) >= 0.5, 
                            ifelse(genes$logFC > 0.5, "Up", "Down"), "Stable")
gg = ggplot(
  # 数据、映射、颜色
  genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("#354259","#C2DED1", "#CDC2AE")) +
  # 注释
  geom_text_repel(
    data = subset(genes, P.Value < 0.05 & abs(genes$logFC) >= 1),
    aes(label = hgnc_symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # 图例
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.background = element_blank(),
        axis.line= element_line(colour = "black"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.4))
ggsave("plot3.png",g,width = 8, height = 8, dpi = 300 )

library(pathfindR)

detach("package:org.Hs.eg.db", unload=TRUE)
library(dplyr)
class(genes)
enrich_data = genes[,c(11,2,5)]
colnames(enrich_data) = c("Gene.symbol","logFC","adj.P.Val")
enrich_data$Gene.symbol
enrich_data = enrich_data[!(enrich_data$Gene.symbol == ""),]
output_df <- run_pathfindR(enrich_data, gene_sets = "GO-MF")
g = enrichment_chart(result_df = output_df, 
                     top_terms = 10)
result_df = output_df
g = enrichment_chart(result_df = output_df, 
                     top_terms = 15)
g = g + scale_colour_gradient(low = "#65647C",high = "#FF8787",na.value = "grey50")  +theme(text = element_text(size = 20))                  
ggsave("plot3.png",g,width = 8, height = 8, dpi = 300 )

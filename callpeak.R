setwd('~/biodata/shTRA2A/')
library("exomePeak2")
library('GenomicFeatures')
txdb <- makeTxDbFromGFF('~/biodata/annotation/GENCODE_NONCODE_ann.gtf',format = 'gtf',organism='Homo sapiens')
exomePeak2(bam_ip = 'markdup/shNC-IP.markdup.bam',
  bam_input ='markdup/shNC-input.markdup.bam',
  txdb=txdb,
  #gff_dir = '~/biodata/annotation/gencode.v35.annotation.sorted.gtf',
  genome = 'hg38',
  p_cutoff = 0.01,
  peak_calling_mode='exon', #我也不好说是不是内存的问题，找个时间再试一下吧
  log2FC_cutoff = 1,
  parallel = TRUE,
  correct_GC_bg = TRUE,
  save_dir = "exomepeak/shNC_genome")
library()
library(ggplot2)
install.packages("languageserver")
require(ggplot2)
data(diamonds)
set.seed(42)
small <- diamonds[sample(nrow(diamonds), 1000), ]
head(small)
p <- ggplot(data = small, mapping = aes(x = carat, y = price))

p + geom_point()
a = rep("1:2", 15)

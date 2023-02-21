library(GenomicFeatures)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters <- promoters(genes(txdb), upstream = 1000, downstream = 0)
rtracklayer::export(promoters,"H:/promoters_ucsc_hg38_knowgene.bed","bed")

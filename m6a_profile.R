library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(ggplot2)
peak <- readPeakFile("H:/peaks.bed")
peak
peak2 <- readPeakFile("H:/macs2_peaks.narrowPeak")
peak2
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=500)
tagMatrix <- getTagMatrix(peak2, windows=promoter)
#tagHeatmap(tagMatrix, xlim=c(-1000, 100), color="red")
plotAvgProf(tagMatrix, xlim=c(-1000, 500),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",conf = 0.95)+
  scale_color_manual(values=c("#00AFBB","#FC4E07")) 















plotPeakProf2(peak = peak2, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)

five_UTR_body <- getTagMatrix(peak = peak2, 
                              TxDb = txdb,
                              upstream = rel(0.2),
                              downstream = rel(0.2), 
                              type = "body",
                              by = "3UTR",
                              weightCol = "V5",
                              nbin = 50)

plotPeakProf(tagMatrix = five_UTR_body, conf = 0.95)

genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")
matrix_actual_extension <- getTagMatrix(peak2,windows = genebody, nbin = 800,
                                        upstream = 1000,downstream = 1000)
plotPeakProf(matrix_actual_extension,conf = 0.95)

tagMatrixList <- lapply("H:/merip_prna/call_peak_file/macs2_peaks.narrowPeak", getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-1000, 500), conf=0.95,resample=500, facet="row")

library(rtracklayer)
library(cliProfiler)
###github r package to get granges file from bed file 
library(bedr)
setwd("H:/")
myfile = 'peak.bed'
my_granges <- bed_to_granges(myfile)
head(my_granges)

meta <- metaGeneProfile(object = test, annotation = "H:/prna_transcript_analysis/gencode_reference/gencode.v40.annotation.gff3/gencode.v40.annotation.gff3")
meta <- metaGeneProfile(object = test, annotation = test_gff3)
meta[[1]]
shown_gff3 <- rtracklayer::import.gff3(test_gff3)
shown_gff3
my_gff3 <- rtracklayer::import.gff3("H:/prna_transcript_analysis/gencode_reference/gencode.v40.annotation.gff3/gencode.v40.annotation.gff3")
my_gff3
library(ggplot2)
## For example if user want to have a new name for the plot
meta[[2]] + ggtitle("Meta Profile 2")
testpath <- system.file("extdata", package = "cliProfiler")
## loading the test GRanges object
test <- readRDS(file.path(testpath, "test.rds"))
class(test)
stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
                               package="Guitar"),
                   system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed6.bed",
                               package="Guitar"))
# Build Guitar Coordinates
txdb_file <- system.file("extdata", "mm10_toy.sqlite",
                         package="Guitar")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(Guitar)
# Guitar Plot


GuitarPlot(txTxdb = txdb,
           stBedFiles = "H:/peaks.bed",
           headOrtail = TRUE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = c("BED12"))






##cliprofiler plot prna region coverage 
library(cliProfiler)
library(rtracklayer)
library(ChIPpeakAnno)
## Extract all the exon annotation
test_anno <- rtracklayer::import.gff3('H:/prna_transcript_analysis/annotation/ps.gtf')
test_anno <- test_anno[test_anno$type == "exon"]


peak<-toGRanges('H:/prna_transcript_analysis/annotation/ps_peaks.bed',format=c('BED'))
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



test_anno <- rtracklayer::import.gff3('H:/prna_transcript_analysis/annotation/pas.gtf')
test_anno <- test_anno[test_anno$type == "exon"]

peak2 <-toGRanges('H:/prna_transcript_analysis/annotation/pas_peaks.bed',format=c('BED'))
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
ggsave(filename = 'H:/fig_project/profile2.pdf',plot = profile2,width = 4.27,height =3.69,units = 'in')
ggsave(filename = 'H:/fig_project/profile1.pdf',plot = profile1,width = 4.27,height =3.69,units = 'in')

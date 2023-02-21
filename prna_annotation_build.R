#load library
library("biomaRt")

#use ensembl mart
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#store the filters
filters <- listFilters(ensembl)

filters$name#look for filters with the string refseq
grep('refseq', filters$name, ignore.case=T, value=T)
#[1] "with_refseq_peptide"            "with_refseq_peptide_predicted"  "with_ox_refseq_mrna"           
#[4] "with_ox_refseq_mrna_predicted"  "with_ox_refseq_ncrna"           "with_ox_refseq_ncrna_predicted"
#[7] "refseq_mrna"                    "refseq_mrna_predicted"          "refseq_ncrna"                  
#[10] "refseq_ncrna_predicted"         "refseq_peptide"                 "refseq_peptide_predicted"

#store the attributes we can fetch
attributes <- listAttributes(ensembl)

#check out the first 10 attributes
head(attributes,10)
#                    name           description
#1        ensembl_gene_id       Ensembl Gene ID
#2  ensembl_transcript_id Ensembl Transcript ID
#3     ensembl_peptide_id    Ensembl Protein ID
#4        ensembl_exon_id       Ensembl Exon ID
#5            description           Description
#6        chromosome_name       Chromosome Name
#7         start_position       Gene Start (bp)
#8           end_position         Gene End (bp)
#9                 strand                Strand
#10                  band                  Band

#create vector of chromosomes
my_chr <- c(1:22,'X','Y')
my_chr
#fetch refseqs
my_refseq <- getBM(attributes='refseq_mrna',
                   filters = 'chromosome_name',
                   values = my_chr,
                   mart = ensembl)

#how many entries 2013 November 28th
length(my_refseq$refseq_mrna)
#[1] 35423
my_refseq = my_refseq$refseq_mrna

#check out the first few
head(my_refseq)
#[1] "NM_001084392" "NM_001355"    "NR_036221"    "NM_138450"    "NM_022840"   
#[6] "NM_173505"

#I only want entries starting with a NM (curated mRNA)
my_refseq <- my_refseq[grep(pattern="^NM",x=my_refseq,perl=T)]
length(my_refseq)
#[1] 33584

#check to see if only NM
head(my_refseq)
#[1] "NM_001084392" "NM_001355"    "NM_138450"    "NM_022840"    "NM_173505"   
#[6] "NM_033453"

#build attribute vector
my_attribute <- c('refseq_mrna',
                  'chromosome_name',
                  'transcript_start',
                  'transcript_end',
                  'strand')

#fetch refseqs and their chromosomal locations
my_refseq_loci <- getBM(attributes=my_attribute,
                        filters = c('refseq_mrna', 'chromosome_name'),
                        values = list(refseq_mrna=my_refseq, chromosome_name=my_chr),
                        mart = ensembl)
                        
dim(my_refseq_loci)
#[1] 33657     5

#how many refseq ids are listed in multiple places?
table(duplicated(my_refseq_loci$refseq_mrna))

#FALSE  TRUE 
#33584    73

#get rid of the duplicated entry
my_refseq_loci <- my_refseq_loci[!duplicated(my_refseq_loci$refseq_mrna),]
dim(my_refseq_loci)
#[1] 33584     5

#convert the strand into '-' and '+'
my_refseq_loci$strand <- gsub(pattern='-1', replacement='-', my_refseq_loci$strand)
my_refseq_loci$strand <- gsub(pattern='1', replacement='+', my_refseq_loci$strand)

#add a 'chr' into the chromosome_name
my_refseq_loci$chromosome_name <- gsub(pattern="^",
                                       replacement='chr',
                                       my_refseq_loci$chromosome_name)

#we now have a data frame of all human mRNA refseqs and their chromosomal locations
head(my_refseq_loci)
#   refseq_mrna chromosome_name transcript_start transcript_end strand
#1 NM_001084392           chr22         24313554       24316773      -
#2    NM_001355           chr22         24313554       24322019      -
#3    NM_138450           chr13         50202435       50208008      +
#4    NM_022840           chr18          2538452        2571485      -
#5    NM_033453           chr20          3190006        3204516      +
#6 NM_001267623           chr20          3190171        3204516      +
my_refseq_loci$score = rep(1, 62356)
my_refseq_loci = my_refseq_loci[,c(2,3,4,1,6,5)]
write.table(my_refseq_loci,"H:/ncbi_ref_mrna.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t')



###get most upstream tss from epd database downloaded all tss bed file 

library(tidyr)

raw_tss = read.table("H:/chrome downlaod/human_epdnew_Btnsl.bed",sep = '\t',header = FALSE)
raw_tss = separate(raw_tss,V4,into = c("ID","sub"),sep = '_')
minus_tss = raw_tss[grep("\\-",raw_tss$V6),]
plus_tss = raw_tss[grep("\\+",raw_tss$V6),]


library(dplyr)

distal_plus = plus_tss %>% 
  group_by(ID) %>% 
  mutate(min_tss=min(V3)) %>% 
  ungroup()

distal_plus= distal_plus[which(distal_plus$V3 == distal_plus$min_tss),]
distal_plus$idv = paste(distal_plus$ID,distal_plus$sub,sep = "_")
distal_plus= distal_plus[,c(1,2,3,9,6,7)]

txnip
library(dplyr)

distal_minus = minus_tss %>% 
  group_by(ID) %>% 
  mutate(min_tss=max(V3)) %>% 
  ungroup()

distal_minus= distal_minus[which(distal_minus$V3 == distal_minus$min_tss),]
distal_minus$idv = paste(distal_minus$ID,distal_minus$sub,sep = "_")
distal_minus= distal_minus[,c(1,2,3,9,6,7)]


most_distal = rbind(distal_minus,distal_plus)
write.table(most_distal,"H:/most_distal.bed",sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)

###



a = read.csv("H:/chrome downlaod/human_epdnew_fFfAl.bed",sep = '\t',header = FALSE)

table(most_distal$idv %in% a$V4)
b = separate(most_distal,idv,into = c("ID","sub"),sep = '_')



###get intersect partion length distribution
overlap = read.table("H:/overlap.txt",sep = ';',fill=TRUE)
overlap$max = apply(overlap, 1,function(x) max(x, na.rm=T))
ggplot(a, aes(length)) +
  geom_bar(fill = "#0073C2FF") 

### bedtools filter prna

library(tidyverse)
psrna = read.csv("H:/pasrna_3.bed",sep = '\t',header = FALSE)
psrna$length = psrna$V3 - psrna$V2
psrna = filter(psrna,length>200)## remove 145 rows less than 200bp(pasrna remove 1455 less than 200 bp)
ggplot(psrna, aes(length)) +
  geom_bar(fill = "#0073C2FF") 
psrna = psrna[,-7]
write.table(psrna,"H:/psrna_ff.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t')

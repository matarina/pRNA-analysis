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
pas = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/pas_output.txt",sep = '\t')
pas = pas[-grep("ERCC-",pas$FEATURE_ID),]
bp1target = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/m6a_igfbp1bp_target_decay/pas_bp3_merip.bed",header = FALSE,sep = "\t")
bp1 = pas[which(pas$FEATURE_ID %in% bp1target$V4),]
non_bp1 =  pas[!(pas$FEATURE_ID %in% bp1target$V4),]

bp1_rep1 = bp1[,c(1,2,4,6,8)]
colnames(bp1_rep1) = c("Geneid","t0","t1","t3","t6")
non_bp1_rep1 = non_bp1[,c(1,2,4,6,8)]
colnames(non_bp1_rep1) = c("Geneid","t0","t1","t3","t6")

bp1_rep2 = bp1[,c(1,3,5,7,9)]
colnames(bp1_rep2) = c("Geneid","t0","t1","t3","t6")
non_bp1_rep2 = non_bp1[,c(1,3,5,7,9)]
colnames(non_bp1_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = bp1_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
bp1_rep1_halflife = obj


cm = non_bp1_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
non_bp1_rep1_halflife = obj

cm = bp1_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
bp1_rep2_halflife = obj


cm = non_bp1_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
non_bp1_rep2_halflife = obj


bp1 = merge(bp1_rep1_halflife,bp1_rep2_halflife,by = "genes")
non_bp1 = merge(non_bp1_rep1_halflife,non_bp1_rep2_halflife,by = "genes")

bp1$mean_nls = rowMeans(subset(bp1,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)
non_bp1$mean_nls = rowMeans(subset(non_bp1,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

bp1_filtered = filter(bp1,mean_nls>0)
non_bp1_filtered = filter(non_bp1,mean_nls>0)
mean(bp1_filtered$mean_nls,trim = 0.2)*60
mean(non_bp1_filtered$mean_nls,trim = 0.2)*60


#ps bp1 vs unbp1 target halflife with m6a modification:273.6971  248.9378
#ps bp3 vs unbp3 target halflife with m6a modification:351.884 230.7303
#pas bp1 vs unbp1 target halflife with m6a modification:314.3783 242.4283
#pas bp3 vs unbp3 target halflife with m6a modification:303.1724 242.0461







###################################total halflife of mrna and prna are calcuclated by this script
#mrna  count Command:featureCounts -T 18 -s 2 -t exon -g gene_id
# -a /home/ma/genomereference/genecode/protein_ercc.gtf -o ./mrna_quant/SRR8131644.count SRR8131644.sort.bam
library(dplyr)
setwd("H:/")
mrna = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/mrna_output.txt",sep = '\t')
mrna = mrna[-grep("ERCC-",mrna$FEATURE_ID),]


mrna_rep1 = mrna[,c(1,2,4,6,8)]
colnames(mrna_rep1) = c("Geneid","t0","t1","t3","t6")
mrna_rep2 = mrna[,c(1,3,5,7,9)]
colnames(mrna_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = mrna_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
mrna_rep1_halflife = obj


cm = mrna_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
mrna_rep2_halflife = obj


mrna = merge(mrna_rep1_halflife,mrna_rep2_halflife,by = "genes")
mrna$mean_nls = rowMeans(subset(mrna,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

mrna_filtered = filter(mrna,mean_nls>0)
mean(mrna_filtered$mean_nls,trim = 0.2)

#mean mrna halflife(trim 0.2)  8.843557h = 530.6134min



#calculate total psRNA halflife 
library(dplyr)
setwd("H:/")
ps = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/ps_output.txt",sep = '\t')
ps = ps[-grep("ERCC-",ps$FEATURE_ID),]


ps_rep1 = ps[,c(1,2,4,6,8)]
colnames(ps_rep1) = c("Geneid","t0","t1","t3","t6")
ps_rep2 = ps[,c(1,3,5,7,9)]
colnames(ps_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = ps_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
ps_rep1_halflife = obj


cm = ps_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
ps_rep2_halflife = obj


ps = merge(ps_rep1_halflife,ps_rep2_halflife,by = "genes")
ps$mean_nls = rowMeans(subset(ps,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

ps_filtered = filter(ps,mean_nls>0)
mean(ps_filtered$mean_nls,trim = 0.2)
#ps total halflife 4.022274h = 241.3364min



#calculate total pasRNA halflife
library(dplyr)
setwd("H:/")
pas = read.csv("H:/prna_transcript_analysis/bp_target_decay_cdf/pas_output.txt",sep = '\t')
pas = pas[-grep("ERCC-",pas$FEATURE_ID),]


pas_rep1 = pas[,c(1,2,4,6,8)]
colnames(pas_rep1) = c("Geneid","t0","t1","t3","t6")
pas_rep2 = pas[,c(1,3,5,7,9)]
colnames(pas_rep2) = c("Geneid","t0","t1","t3","t6")


library(edgeR)
cm = pas_rep1
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
pas_rep1_halflife = obj


cm = pas_rep2
dge <- DGEList(counts = cm[,-1], genes = cm[,1])
obj = dge
obj = apply_expression_cutoff(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10)
obj = trim_late_time_points(DGeList = obj, CPMcutoff= 0.5)
obj = calculate_normalization_factors3(DGEList = obj, method = c("mean"))
obj = calculate_half_life(DGEList = obj)
pas_rep2_halflife = obj


pas = merge(pas_rep1_halflife,pas_rep2_halflife,by = "genes")
pas$mean_nls = rowMeans(subset(pas,select = c("hl_nls.x","hl_nls.y")),na.rm = TRUE)

pas_filtered = filter(pas,mean_nls>0)
mean(pas_filtered$mean_nls,trim = 0.2)
#pas total halflife 4.203391h = 252.2035min





#

#mean mrna halflife(trim 0.2)  8.843557 (without m6a modification)
#halflife  pas bp1=4.722062=283.32min non_bp1=4.088213=245.29min bp3=4.765821=285.95min non_bp3=3.943781=236.63min
#halflife   ps bp1=5.429605=325.78min non_bp1=3.930975=235.85min  bp3=5.666296=339.98min non_bp3=3.847075=230.82min



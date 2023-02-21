##deseq2 get diffrernt expressed prnas
library(DESeq2)
library(tidyr)
library(dplyr)
ps_rawCounts <- read.csv("H:/prna_transcript_analysis/hcc_tissue/mongolia/ps_gene.csv")
ps_rawCounts$geneid = gsub('^','ps_',ps_rawCounts$geneid)
pas_rawCounts <- read.csv("H:/prna_transcript_analysis/hcc_tissue/mongolia/pas_gene.csv")
pas_rawCounts$geneid = gsub('^','pas_',pas_rawCounts$geneid)
prna = rbind(ps_rawCounts,pas_rawCounts)
prna = prna %>% column_to_rownames('geneid')
rawCounts = prna


patient = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/patient_sample_metadata.txt',sep = '\t')
patient = patient %>% drop_na(RNASeq_T)
a = patient
a = a[,-5]
patient = patient[,-6]
colnames(a)[5] = 'RNASeq_T'
c = rbind(a,patient)
sra = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/SraRunTable.txt') %>% select(Run,Sample.Name)
srx = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/sra_result.csv') %>% select(Experiment.Accession,Experiment.Title) %>% 
      separate(Experiment.Title,into = c('a','b','c','d'),sep = ':|;')
srr = merge(sra,srx,by.x = 'Sample.Name',by.y = 'a')
srr = srr %>% select(Run,b) %>% separate(b,into = c('b1','b2','b3'),sep = ' \\[|\\]') %>% separate(b1,into = c('c1','c2','c3','c4'),sep = ' ')
srr = srr %>%  filter(Run != 'SRR14511012'& Run != 'SRR14511013')
meta = merge(c,srr,by.x='RNASeq_T',by.y = 'b2')
write.csv(meta,'H:/prna_transcript_analysis/hcc_tissue/mongolia/omic_pipe/mongolia_meta.csv',quote = FALSE,row.names = FALSE)















## mechine learning predict biomarker
rm(list=ls())
library(dplyr)
library(tidyverse)
library(caret)
library(pROC)
patient = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/patient_sample_metadata.txt',sep = '\t')
patient = patient %>% drop_na(RNASeq_T)
sra = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/SraRunTable.txt')
sra = sra %>% select(Run,Sample.Name)
sra = sra[-c(141,142),]
srx = read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/sra_result.csv')
srx  = srx %>% select(Experiment.Accession,Experiment.Title) %>% separate(Experiment.Title,into = c('a','b','c','d'),sep = ':|;')
srr = merge(sra,srx,by.x = 'Sample.Name',by.y = 'a')
srr = srr %>% select(Run,b) %>% separate(b,into = c('b1','b2','b3'),sep = ' \\[|\\]') %>% separate(b1,into = c('c1','c2','c3','c4'),sep = ' ')
meta = merge(patient,srr,by.x='Patient',by.y = 'c3')


data1= read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/pas_gene.csv')
data2= read.csv('H:/prna_transcript_analysis/hcc_tissue/mongolia/ps_gene.csv')
data1$geneid = gsub('^','pas_',data1$geneid)
data2$geneid = gsub('^','ps_',data2$geneid)
data = rbind(data1,data2)
rownames(data) = data$geneid
different_pRNA_index = readRDS('H:/prna_transcript_analysis/hcc_tissue/mongolia/different_pRNA_index.rds')
data = data %>% filter(geneid %in% different_pRNA_index)
counts_table = data[,-1]

counts_table = counts_table %>% select(-c(SRR10969228,SRR10969234, SRR10969260, SRR10969262, SRR10969263 ,SRR10969265, SRR10969266, SRR10969288,
                         SRR10969289, SRR10969291, SRR10969292))


counts_to_logcpm <- function(counts_table) {
  vapply(as.data.frame(counts_table), 
         function(x) log10(x * 1e+06/sum(x) + 1),
         c(rep(1.0,nrow(counts_table)))) %>%
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(counts_table)) %>% 
    magrittr::set_rownames(rownames(counts_table))
}

logcpm_table = counts_table %>% counts_to_logcpm() %>% base::t() %>% base::as.data.frame() %>% drop_na()
descrCorr = cor(logcpm_table )
highCorr = findCorrelation(descrCorr, 0.7)
newdata2 =logcpm_table[, -highCorr]
meta2 = meta %>% filter(c4 == 'tumor') %>% select(Run,obesity) %>% column_to_rownames(var = 'Run')
meta_data = merge(meta2,newdata2,by = 'row.names')
meta_data = meta_data %>% column_to_rownames(var = 'Row.names') %>% drop_na()
meta_data$obesity = gsub('0','obesity',meta_data$obesity)
meta_data$obesity = gsub('1','nonobesity',meta_data$obesity)
meta_data$obesity = as.factor(meta_data$obesity)


library(randomForest)
library(doParallel)
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

seeds <- vector(mode = "list", length = 31)#length is = (n_repeats*nresampling)+1 10*3+1
for(i in 1:30) seeds[[i]]<- sample.int(n=10000, 3392) #
seeds[[31]]<-sample.int(1000, 1)
rfeControl <- rfeControl(functions=rfFuncs, seeds = seeds,
                         method="cv",
                         number=10)
rfebreaks <- seq(2, ncol(meta_data), 2)## step size 2 to gradient select variables
##
cl <- makeCluster(12, type="SOCK",outfile="")
registerDoParallel(cl)

rfefit <- rfe(meta_data[,-1],
              meta_data$obesity,
              sizes=rfebreaks,
              rfeControl=rfeControl)
saveRDS(rfefit,'H:/prna_transcript_analysis/hcc_tissue/mongolia/obesity_rfefit.rds')
rfefit = readRDS('H:/prna_transcript_analysis/hcc_tissue/mongolia/rfefit.rds')
stopCluster(cl)
rfeAccuracy <- tibble(Accuracy = signif(rfefit$results$Accuracy, 3),
                      Variables = rfefit$results$Variables) %>%
  mutate(label=ifelse(Accuracy == max(Accuracy),
                      paste0("(Features=", Variables, ", Accuracy=", Accuracy,")"), NA))

rfepl <- ggplot(data=rfeAccuracy, aes(x=Variables, y=Accuracy))+
  geom_point(color="grey", size=3, shape=19)+ xlim(NA, 30)+
  geom_line(color="black", linetype=1, linewidth=1)+
  # geom_text(aes(label=label), nudge_y=0.002)+
  annotate(geom="point",
           x=rfeAccuracy[grep(max(rfeAccuracy$Accuracy), rfeAccuracy$Accuracy), ]$Variables,
           y=max(rfeAccuracy$Accuracy), color="red", size=4)+
  labs(x="Features (Numbers)",
       y="Accuracy (Bootstrap)")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent"),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(size=.5, color="black"),
        axis.title = element_text(color='black', size=12, face="bold"),
        axis.text = element_text(color='black', size=10),
        text = element_text(size=8, color="black", family="serif"))
rfepl

rfeImp <- varImp(rfefit) %>% 
  rownames_to_column("feature") %>%
  filter(feature%in%rfefit$optVariables) %>%
  arrange(desc(Overall))
caret::featurePlot(meta_data[,rfeImp$feature], meta_data$obesity, plot="box",par.settings = list( box.umbrella=list(col= c("red", "green")), 
                                                                                              box.dot=list(col= c("red", "green")), 
                                                                                              box.rectangle = list(col= c("red", "green")) 
))
library(ggpubr)

 ggdotchart(rfeImp, x = "feature", y = "Overall",
           color = "red",   size = 3,                            
           sorting = "ascending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           group = "Overall", 
           dot.size = 5
           
)

 
 
 newdata3 <- meta_data %>% dplyr::select(c("obesity",rfeImp$feature))
 set.seed(1234)
 inTrain = createDataPartition(newdata3$obesity, p = .75, list = FALSE,times = 1)
 new_data_train = newdata3[inTrain,]
 new_data_test = newdata3[-inTrain,]
 newdata3_preProcess <- preProcess(newdata3, method = c("center", "scale"))
 
 
 
 trainTransformed <- predict(newdata3_preProcess, new_data_train)
 testTransformed <- predict(newdata3_preProcess, new_data_test)
 dat_table <- as.matrix(trainTransformed %>% dplyr::select(-obesity) %>% data.frame())
 plot(density(dat_table))
 
 
 
 kfold=5
 times=10
 
 if(dim(trainTransformed)[1] > 10 * 8){
   kfold <- 10
 }else{
   kfold <- kfold
 }
 set.seed(12)
 myControl <- trainControl(method = "repeatedcv", 
                           number = kfold,
                           repeats = times,
                           search = "grid",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,                           
                           verboseIter = TRUE)
 
 unregister_dopar <- function() {
   env <- foreach:::.foreachGlobals
   rm(list=ls(name=env), pos=env)
 }
 unregister_dopar()
 
 fit <- train(trainTransformed[,-1],trainTransformed$obesity,
              method = "rf",
              trControl = myControl,
              tuneLength = 15,
              metric = "ROC",
              verbose = FALSE) 
 fit_Acc <- tibble(Feartures = fit$results$mtry,
                   ROC = signif(fit$results$ROC, 3),
                   Sens = signif(fit$results$Sens, 3),
                   Spec = signif(fit$results$Spec, 3)) %>%
   mutate(label=ifelse(ROC == max(ROC), 
                       paste0("(Features=", Feartures, "\n ROC=", ROC,")"), NA))
 fit_final <- fit
 fit_Acc_final <- fit_Acc
 testData_final <- testTransformed    
 
 fit_Acc_pl <- ggplot(data=fit_Acc_final, aes(x=Feartures, y=ROC))+
   geom_point(color="grey", size=3, shape=19)+
   geom_line(color="black", linetype=1, size=1)+
   geom_text(aes(label=label), nudge_y=-0.001)+
   annotate(geom="point",
            x=fit_Acc_final[grep(max(fit_Acc_final$ROC),
                                 fit_Acc_final$ROC), ]$Feartures,
            y=max(fit_Acc_final$ROC), color="red", size=4) +
   labs(x="Randomly Selected Predictors",
        y="ROC (Repeated Cross-Validation)")+
   theme_bw()+
   theme(panel.background = element_rect(fill="transparent"),
         axis.ticks.length = unit(0.4, "lines"),
         axis.ticks = element_line(color='black'),
         axis.line = element_line(size=.5, color="black"),
         axis.title = element_text(color='black', size=12, face="bold"),
         axis.text = element_text(color='black', size=10),
         text = element_text(size=8, color="black", family="serif"))
 fit_Acc_pl
 pred_raw <- predict(fit_final, newdata = testData_final, type = "raw")
 print(confusionMatrix(pred_raw, testData_final$obesity))
 pred_prob <- predict(fit_final, newdata = testData_final, type = "prob")  
 rocobj <- roc(testData_final$obesity, pred_prob[, 1])
 auc <- round(auc(testData_final$obesity, pred_prob[, 1]), 3)
 roc <- tibble(tpr=rocobj$sensitivities,
               fpr=1 - rocobj$specificities)
 roc_pl <- ggplot(data=roc, aes(x=fpr, y=tpr))+
   geom_path(color="red", size=1)+
   geom_abline(intercept=0, slope=1, color="grey", size=1, linetype=2)+
   labs(x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensivity or Recall)")+
   annotate("text", x=.75, y=.25, label=paste("AUC =", auc),
            size=5, family="serif")+
   coord_cartesian(xlim=c(0, 1), ylim=c(0, 1))+
   theme_bw()+
   theme(panel.background = element_rect(fill="transparent"),
         axis.ticks.length = unit(0.4, "lines"),
         axis.ticks = element_line(color="black"),
         axis.line = element_line(size=.5, color="black"),
         axis.title = element_text(color='black', size=12, face="bold"),
         axis.text = element_text(color='black', size=10),
         text = element_text(size=8, color="black", family="serif"))
 
 roc_pl
 
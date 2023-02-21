rm(list = ls())
##biomakrer analysis
library(magrittr)
library(caret)
library(tibble)
library(tidyverse)
library(ggplot2)
library(pROC)
data1= read.csv('H:/prna_transcript_analysis/plasma/pas_output.txt',sep = '\t',skip = 1)
data2= read.csv('H:/prna_transcript_analysis/plasma/ps_output.txt',sep = '\t',skip = 1)
data1$Geneid = gsub('^','pas_',data1$Geneid)
data2$Geneid = gsub('^','ps_',data2$Geneid)
data = rbind(data1,data2)
different_pRNA_index = readRDS('H:/prna_transcript_analysis/plasma/different_pRNA_index.rds')
data = data[which(data$Geneid %in% different_pRNA_index),]

rownames(data) = data$Geneid
counts_table = data[,-1]

counts_to_logcpm <- function(counts_table) {
  vapply(as.data.frame(counts_table), 
         function(x) log10(x * 1e+06/sum(x) + 1),
         c(rep(1.0,nrow(counts_table)))) %>%
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(counts_table)) %>% 
    magrittr::set_rownames(rownames(counts_table))
}

logcpm_table = counts_table %>% counts_to_logcpm() %>% base::t() %>% base::as.data.frame()



vali_data1= read.csv('H:/prna_transcript_analysis/vali/pas_output.txt',sep = '\t',skip = 1)
vali_data2= read.csv('H:/prna_transcript_analysis/vali/ps_output.txt',sep = '\t',skip = 1)
vali_data1$Geneid = gsub('^','pas_',vali_data1$Geneid)
vali_data2$Geneid = gsub('^','ps_',vali_data2$Geneid)
vali_data = rbind(vali_data1,vali_data2)
vali_data = vali_data[which(vali_data$Geneid %in% different_pRNA_index),]

rownames(vali_data) = vali_data$Geneid
counts_table2 = vali_data[,-1]
counts_to_logcpm2 <- function(counts_table2) {
  vapply(as.data.frame(counts_table2), 
         function(x) log10(x * 1e+06/sum(x) + 1),
         c(rep(1.0,nrow(counts_table2)))) %>%
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(counts_table2)) %>% 
    magrittr::set_rownames(rownames(counts_table2))
}
logcpm_table2 = counts_table2 %>% counts_to_logcpm() %>% base::t() %>% base::as.data.frame()

logcpm_table2 = logcpm_table2[-which(rownames(logcpm_table2) %in% rowkeep),]
rowkeep = c('SRR14506729','SRR14506731','SRR14506732','SRR14506843','SRR14506845'
            ,'SRR14506751','SRR14506848','SRR14506860','SRR14506880','SRR14506883','SRR14506698','SRR14506789')
rowkeep = gsub('$','.sort.deduplicated.bam',rowkeep)
 rowkeep %in% rownames(logcpm_table2)

      
#zerologcpm_table=nearZeroVar(logcpm_table)
#newdata1=logcpm_table[,-zerologcpm_table]
descrCorr = cor(logcpm_table )
highCorr = findCorrelation(descrCorr, 0.9)
newdata2 =logcpm_table[, -highCorr]
#comboInfo = findLinearCombos(newdata2)
#length(comboInfo$remove)
#newdata2=newdata2[, -comboInfo$remove]

group = as.factor(rep(c('cancer','normal'),times=c(35,30)))
newdata2$group = as.factor(rep(c('cancer','normal'),times=c(35,30)))
set.seed(1)

library(doParallel)
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
set.seed(123)
seeds <- vector(mode = "list", length = 11)#length is = (n_repeats*nresampling)+1
for(i in 1:10) seeds[[i]]<- sample.int(n=10000, 1433) #(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)

seeds[[11]]<-sample.int(1000, 1)





rfeControl <- rfeControl(functions=rfFuncs, seeds = seeds,
                         method="cv",
                         number=10)
rfebreaks <- seq(2, ncol(newdata2), 2)## step size 2 to gradient select variables

cl <- makeCluster(12, type="SOCK",outfile="")
registerDoParallel(cl)

rfefit <- rfe(newdata2[-ncol(newdata2)],
              newdata2$group,
              sizes=rfebreaks,
              rfeControl=rfeControl)
stopCluster(cl)


rfeAccuracy <- tibble(Accuracy = signif(rfefit$results$Accuracy, 3),
                      Variables = rfefit$results$Variables) %>%
  mutate(label=ifelse(Accuracy == max(Accuracy),
                      paste0("(Features=", Variables, ", Accuracy=", Accuracy,")"), NA))

rfepl <- ggplot(data=rfeAccuracy, aes(x=Variables, y=Accuracy))+
  geom_point(color="grey", size=3, shape=19)+ xlim(NA, 30)+
  geom_line(color="black", linetype=1, size=1)+
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
caret::featurePlot(newdata2[,rfeImp$feature], newdata2$group, plot="box",par.settings = list( box.umbrella=list(col= c("red", "green")), 
                                                                                              box.dot=list(col= c("red", "green")), 
                                                                                              box.rectangle = list(col= c("red", "green")) 
))


library(ggpubr)

ggdotchart(rfeImp, x = "feature", y = "Overall",
           color = "red",   size = 5,                            
           sorting = "ascending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           group = "Overall", 
           dot.size = 9
           
)



#### independent test data 
newdata3 <- newdata2 %>% dplyr::select(c("group",rfeImp$feature))
newdata3_preProcess <- preProcess(newdata3, method = c("center", "scale"))
trainTransformed <- predict(newdata3_preProcess, newdata3)
dat_table <- as.matrix(Transformed %>% dplyr::select(-group) %>% data.frame())
plot(density(dat_table))

#logcpm_table2 = logcpm_table2[1:40,]
logcpm_table2$group = rep(c("cancer","normal"),times = c(22,36))
logcpm_table3 = logcpm_table2%>% dplyr::select(c("group",rfeImp$feature))
logcpm_table3_preProcess <- preProcess(logcpm_table3, method = c("center", "scale"))
testTransformed <- predict(logcpm_table3_preProcess, logcpm_table3)
test_dat_table <- as.matrix(testTransformed %>% dplyr::select(-group) %>% data.frame())



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

fit <- train(trainTransformed[,-1],trainTransformed$group,
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

pred_raw <- predict(fit_final, newdata = testTransformed, type = "raw")
print(confusionMatrix(pred_raw, as.factor(testData_final$group)))
pred_prob <- predict(fit_final, newdata = testTransformed, type = "prob")  
rocobj <- roc(testTransformed$group, pred_prob[, 1])
auc <- round(auc(testTransformed$group, pred_prob[, 1]), 3)
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























########self devided test
newdata3 <- newdata2 %>% dplyr::select(c("group",rfeImp$feature))
set.seed(1234)
inTrain = createDataPartition(newdata3$group, p = .75, list = FALSE,times = 1)
new_data_train = newdata3[inTrain,]
new_data_test = newdata3[-inTrain,]
newdata3_preProcess <- preProcess(newdata3, method = c("center", "scale"))



trainTransformed <- predict(newdata3_preProcess, new_data_train)
testTransformed <- predict(newdata3_preProcess, new_data_test)
dat_table <- as.matrix(trainTransformed %>% dplyr::select(-group) %>% data.frame())
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

fit <- train(trainTransformed[,-1],trainTransformed$group,
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
print(confusionMatrix(pred_raw, testData_final$group))
pred_prob <- predict(fit_final, newdata = testData_final, type = "prob")  
rocobj <- roc(testData_final$group, pred_prob[, 1])
auc <- round(auc(testData_final$group, pred_prob[, 1]), 3)
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






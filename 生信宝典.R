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
group = as.factor(rep(c('cancer','normal'),times=c(35,28)))


set.seed(1)
train_index = createDataPartition(group,p=0.75,list = F)
train_data = logcpm_table[train_index,]
train_data_group = group[train_index]
test_data = logcpm_table[-train_index,]
test_data_group = group[-train_index]
dim(test_data)

library(Boruta)
set.seed(2)
boruta <- Boruta(x=train_data, y=train_data_group, pValue=0.01, mcAdj=T, 
                 maxRuns=300)
boruta



boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)
head(boruta.variable.imp)
library(ImageGP)
sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90)

boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")
caret::featurePlot(train_data[,boruta.finalVarsWithTentative$Item], train_data_group, plot="box")
generateTestVariableSet <- function(num_toal_variable){
  max_power <- ceiling(log10(num_toal_variable))
  tmp_subset <- c(unlist(sapply(1:max_power, function(x) (1:10)^x, simplify = F)), ceiling(max_power/3))
  #return(tmp_subset)
  base::unique(sort(tmp_subset[tmp_subset<num_toal_variable]))
}
boruta_train_data <- train_data[, boruta.finalVarsWithTentative$Item]
boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
trControl <- trainControl(method="repeatedcv", number=10, repeats=5)
set.seed(3)
tuneGrid <- expand.grid(mtry=boruta_mtry)
borutaConfirmed_rf_default <- train(x=boruta_train_data, y=train_data_group, method="rf", 
                                    tuneGrid = tuneGrid, # 
                                    metric="Accuracy", #metric='Kappa'
                                    trControl=trControl)
borutaConfirmed_rf_default
plot(borutaConfirmed_rf_default)
dotPlot(varImp(borutaConfirmed_rf_default))
borutaConfirmed_rf_default_finalmodel <- borutaConfirmed_rf_default$finalModel   
prediction_prob <- predict(borutaConfirmed_rf_default_finalmodel, newdata=test_data, type="prob")
library(pROC)
roc_curve <- roc(test_data_group, prediction_prob[,1])
roc_curve
roc <- roc(test_data_group, factor(prediction_prob, ordered=T))
plot(roc)







trControl <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
rf_random <- train(x=train_data, y=train_data_group, method="rf", 
                   tuneLength = 15, # 随机15个参数值或参数值组合
                   metric="Accuracy", #metric='Kappa'
                   trControl=trControl)
print(rf_random)
tunegrid <- expand.grid(mtry=c(556))
for (ntree in c(500,700, 800, 1000)) {
  set.seed(1)
  fit <- train(x=train_data, y=train_data_group, method="rf", 
               metric="Accuracy", tuneGrid=tunegrid, 
               trControl=trControl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

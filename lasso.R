library("glmnet") 
library("survival")
library(dplyr)
library(tidyverse)
#
ps = read.csv('H:/ps_output.txt',sep = '\t')
pas = read.csv('H:/pas_output.txt',sep = '\t')
ps$strand = rep('ps_',nrow(ps))
pas$strand = rep('pas_',nrow(pas))
ps$Geneid = paste0(ps$strand,ps$Geneid)
pas$Geneid = paste0(pas$strand,pas$Geneid)
ps = subset(ps,select = -c(strand))
pas = subset(pas,select = -c(strand))
expr = rbind(ps,pas)
rownames(expr) = expr$Geneid
expr = expr[,-1]
exprSet = expr


#logCPM or logTPM data transform
exprSet =log2(edgeR::cpm(exprSet)+1)
k = apply(exprSet,1, function(x){sum(x>0)>0.5*ncol(exprSet)});table(k)
exprSet = exprSet[k,]

#merge expression data and clinical data
expr_cpm = exprSet
nrow(expr_cpm)
sample = read.csv('H:/SraRunTable.txt',sep = ',')
sample = sample[,c('Run','Sample.Name')]

id =  read.csv("H:/sample_id.txt",sep = '\t',header = FALSE)
colnames(id) = c('Sample.Name','V2')
id = separate(id,V2,into = c("v2","v3","sample"),sep = '(_|])')
id = id[grep('(TG|DH013)',id$sample),]
id = id[,-c(2,3)]

clinical = read.csv('H:/clinical_meta_data.txt',sep = '\t')
clinical_fileter = merge(clinical,id,by='sample')
cli = merge(clinical_fileter,sample,by='Sample.Name')
cli_stat = subset(cli,select=c('vital','os','Run'))



colnames(expr_cpm) <- gsub("\\.sort\\.bam", "", colnames(expr_cpm))

t_expr = data.frame(t(expr_cpm))
t_expr$Run = row.names(t_expr)
rt = merge(cli_stat,t_expr,by='Run')
rt =rt[,c(1,3,2,4:ncol(rt))]

colnames(rt)[2:3] = c('futime','fustat')
rt$fustat = gsub('Alive','0',rt$fustat)
rt$fustat = gsub('Dead','1',rt$fustat)

rt$fustat = as.numeric(rt$fustat)
rt$futime = round(rt$futime/12,digits = 0)
rownames(rt) = rt$Run
rt = rt %>% select(!Run)
rt = rt[,!apply(rt==0,2,all)]
rt<- filter(rt, !futime == 0)
set.seed(123)   #设定随机种子

x=as.matrix(rt[,c(3:ncol(rt))])

y=data.matrix(Surv(rt$futime,rt$fustat))

fit=glmnet(x, y, family = "cox")
fit <- glmnet(x, y, family = 'cox', type.measure = "deviance", nfolds = 10)
print(fit)
plot(fit, xvar = "lambda", label = TRUE)


cvfit <- cv.glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10)

plot(cvfit)

lambda.1se <- cvfit$lambda.1se
lambda.1se
model_lasso_1se <- glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10,lambda = lambda.1se)
gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]#as.numeric后"."会转化为0
gene_1se

coef = coef(fit, s = cvfit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef




FinalGeneExp = rt[,lassoGene] 
myFun = function(x){crossprod(as.numeric(x),actCoef)} 
riskScore = apply(FinalGeneExp,1,myFun) 
outCol = c("futime", "fustat", lassoGene) 
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low")) 
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
library(ggpubr)   #使用ggpubr包绘制散点图

p <- ggboxplot(dat, x = "fustat", y = "riskScore", 
               color = "fustat", palette = "jco", 
               add = "jitter") 
p <- p + stat_compare_means()   #  Add p-value 
p   #得出预测结果

pred <- prediction(dat$riskScore, dat$fustat) 
auc_min <- performance(pred,"auc")@y.values[[1]] 
perf <- performance(pred,"tpr","fpr") 
performance(pred,"auc")   # shows calculated AUC for model 
plot(perf,colorize=FALSE, col="red")   #绘制ROC曲线 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线 
library(glmnet) 
library(caret)

FinalGeneExp = rt[,lassoGene] 
myFun = function(x){crossprod(as.numeric(x),actCoef)} 
riskScore = apply(FinalGeneExp,1,myFun) 
outCol = c("futime", "fustat", lassoGene) 
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low")) 
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
pred <- prediction(dat$riskScore, dat$fustat) 
perf <- performance(pred,"tpr","fpr") 
performance(pred,"auc")   # shows calculated AUC for model 
plot(perf,colorize=FALSE, col="red")   #绘制ROC曲线 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )



##survival Kaplan-Meier plot
library(survminer)
rt = rt[,c(1:3)]
colnames(rt)[2:3] = c('time','status')
#rt$time = rt$time*30
dat$Run = rownames(dat)
sur_dat = merge(rt,dat,by = 'Run')
sur_dat$risk = gsub('low','1',sur_dat$risk)
sur_dat$risk = gsub('high','2',sur_dat$risk)
sur_dat$status  = gsub('1','2',sur_dat$status)
sur_dat$status  = gsub('0','1',sur_dat$status)
sur_dat$status = as.numeric(sur_dat$status)
sur_dat$time = as.numeric(sur_dat$time)
sur_dat$risk = as.numeric(sur_dat$risk)
ggsurvplot(fit)
data(lung) # 加载lung数据集
View(lung) # 查看数据集
attach(lung) 

dat$time = dat$futime
Surv(time,sex) 
fit <- survfit(Surv(time,status) ~ ph.karno,  # 创建生存对象 
               data = dat) 
ggsurvplot(fit, data = dat)
dat$Run

library(data.table)
data = fread('H:/chrome downlaod/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena',header = T,stringsAsFactors = F)
head(data)
library(tibble)
library(dplyr)
library(forcats)
data = data[12309,]
colnames(data)
data = column_to_rownames(data,'sample')

tdata = data.frame(t(data))
tdata = rownames_to_column(tdata,var = 'sample')
colnames(tdata) = c('sample','OGFR')
tdata = tdata[-1,]
tdata$sample = gsub('\\.','-',tdata$sample)
clinical = read.csv('H:/chrome downlaod/Survival_SupplementalTable_S1_20171025_xena_sp',sep = '\t')
subset_clinical = clinical %>%  select('sample','cancer.type.abbreviation')
merge  = merge(subset_clinical,tdata,by = 'sample')
colnames(merge)[2] = 'cancer'
merge$OGFR = as.numeric(merge$OGFR)
merge = merge %>% mutate( cancer= fct_reorder(cancer, OGFR, .fun='median',.desc = TRUE))

library(ggplot2)
library(RColorBrewer)


ggplot(merge, aes(cancer,OGFR,fill=cancer)) +  geom_jitter(width =0.2,shape = 21,size=0.05)+
  geom_boxplot(alpha=0.7)+theme_classic()+labs(y = "OGFR Expression",x = 'Cancer Type')

subset_clinical2 = clinical %>%  select('vital_status','OS.time','sample')
survival_merge = merge(merge,subset_clinical2,by='sample')
survival_merge$mouth=round(survival_merge$OS.time/30,2)

library("survival")
library("survminer")

OV = filter(survival_merge,cancer == 'OV')
OV$level = ifelse(OV$OGFR>median(OV$OGFR),'high','low')
survData = Surv(time=OV$mouth, 
                event=OV$vital_status=='Dead') 
KMfit <- survfit(survData ~ OV$level)
png("OVplot.png",width=8,height=8,units="in",res=400)
myplot = ggsurvplot(KMfit, 
           data = OV,  
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, # 添加风险表
           ncensor.plot = FALSE, #？？图
           xlab = "Follow up time(m) in OV", # 指定x轴标签
           break.x.by = 10, # 设置x轴刻度间距
           palette = c("#E7B800", "#2E9FDF"),
           #legend = c(0.8,0.75), # 指定图例位置
           legend.title = 'OGFR', # 设置图例标题
           #legend.labs = c("old", "young"), # 指定图例分组标签
) 
print(myplot)
dev.off()























































aaa = read.csv('H:/chrome downlaod/POSTAR3_CLIPdb_module_RBP_binding_sites_sub_info.csv')
colnames(aaa)[1:2] = c('geneid','ensembl_id')
aaa$geneid
library(xlsx)
bbb = read.csv('H:/41420_2020_328_MOESM5_ESM.csv')

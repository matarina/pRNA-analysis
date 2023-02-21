library(dplyr)
df <- data.frame(group=c("all TSS","psRNA","pasRNA"),
                 len=c(16455,15442,12825))
df$group <- factor(df$group, levels=c("all TSS", "psRNA", "pasRNA"))
mypal <- carto.pal(pal1 = "pastel.pal", n1 = 3,middle = TRUE, transparency = FALSE)
mypal
p<-ggplot(data=df, aes(x=group, y=len,fill=group)) +
  geom_bar(stat="identity")+
 scale_fill_manual(values = c("#DA9191","#6978DB","#69D282"))+
  geom_text(aes(label = len), size = 5,color="white",vjust= 1.81, position = "stack")+
  labs(x = "TSS")+labs(y = "Numbers")+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size = 20),panel.background = element_blank())+
  theme_classic() +
  theme(axis.text = element_text(size = 20),axis.title.y = element_text(size = 20),legend.title=element_text(size=20),legend.text=element_text(size=20),axis.title.x = element_text(size = 20),axis.text.x=element_blank())   
p
library(ggplot2)
library(cartography)
mypal <- carto.pal(pal1 = "pastel.pal", n1 = 3,middle = TRUE, transparency = FALSE)
##595 482
mypal


a = read.csv("H:/pas.bed",header = FALSE,sep = '\t')
length(unique(a$V4))

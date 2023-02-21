setwd('H:/')


df <- c("pasRNA_IGF2BP3" = 103, "pasRNA_IGF2BP3_m6A" = 290)
data <- data.frame("group" = names(df), df)
vals <- paste(data$df, sep = "")
pal <- c("#68B8E1","#EC6161")
library(plotly)
p = plot_ly(data, labels = ~group, values = ~df, type = 'pie',textinfo = "text + values", text = vals,
             marker = list(colors = pal),textfont = list( size = 15)) %>%
  layout(title = "",legend = list(orientation = "h",xanchor = "center",  x = 0.5),    
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),legend = list(font = list(size = 24)))
p
#527 535

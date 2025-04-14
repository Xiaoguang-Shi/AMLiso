library(ggplot2)
library(ggalluvial)
setwd("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20241229/Figure/Figure4")
ELN=read.csv("eln.csv",header=T,row.names = 1)
link <- ELN[,c("ELN","Adjusted_ELN")] 
link $link=1 
link =reshape::melt(link ,id='link')
variable = summary(link $variable)
link$flow=rep(1:variable[1],length(variable))

link$value <- factor(link$value, levels = c(1,2,3,"Low","High"))

colors <- c("1" = "#B9CEE3", "2" = "#FFE5C1", "3" = "#F6BEB9","Low"="#B9CEE3","High"="#F6BEB9")
ggplot(link ,aes(x=variable,
                 y=link,
                 stratum = value ,
                 alluvium = flow,
                 fill = value))+
  geom_stratum()+ #
 scale_fill_manual(values = colors) +  # 设置颜色映射
  geom_flow(aes.flow = 'forward')+
  geom_text(stat = 'stratum', infer.label=T,size=2.5)+
  theme_minimal() +
  geom_alluvium() +
  theme(axis.title.y = element_blank(),  # 去除纵坐标
        axis.text.y = element_blank(),   # 去除纵坐标刻度标签
        axis.ticks.y = element_blank(),  # 去除纵坐标刻度线
        panel.grid = element_blank(),  # 去除背景线
        legend.position = 'none')  # 去除图例

ELN=read.csv("eln.csv",header=T,row.names = 1)
library(dplyr)

ELN$futime=ELN$futime/30
# 生存曲线绘制
survfit(Surv(futime, fustat) ~ ELN, data = ELN) %>%
  ggsurvplot(
    data =ELN,
    conf.int = T,
    pval = TRUE,
    pval.size = 6,
    legend = "none",
    pval.coord=c(30,0.6),
    xlab = "Time(Months)",
    break.time.by = 10,
    legend.title = "ELN",
    legend.labs = c("Favorable","Intermedite","Adverse"), 
    surv.median.line = "hv",
    risk.table = T,
    palette = c("#00599F","orange","#d80700"),
    cumevents = F,
    risk.table.height = 0.3
  )

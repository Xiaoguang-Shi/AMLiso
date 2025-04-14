library(ggplot2)
increase_data=read.csv("F2B/increase_number.csv",header=T,check.names = F)

ggplot(data = increase_data, mapping = aes(x = id, y = `number`, color = `classification`)) + 
  geom_point() + 
  #geom_line() +
  geom_smooth(se=FALSE)+
  scale_color_manual(values=c("novel"='red',"known"='blue'))+
  labs(x = "Number of samples", y = "Increase numbers of isoforms",fill="") +
  theme_bw() + 
  theme(panel.grid=element_blank())+
  theme(legend.position =c(0.9,0.5))
# Figure2

```{r}
library(ggplot2)

isoform=read.csv("number.csv",check.names = F)
ggplot(isoform, aes(x=id, y=number, fill=factor(type, c("novel", "known")))) +
  geom_bar(stat = "identity",position = 'fill') +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Samples", y = "Numbers of isoforms", fill = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = c(0.5,0.999), 
        #legend.justification = c("center", "top"), 
        legend.box = "horizontal",legend.background = element_rect(fill = "transparent"), 
        legend.margin = margin(t = 0, b = 0, unit = "cm")) +
  guides(fill = guide_legend(keywidth = unit(1, "cm"), keyheight = unit(0.15, "cm"), nrow = 1))
```

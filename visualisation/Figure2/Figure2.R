#devtools::install_github("dzhang32/ggtranscript")
library(ggplot2)
library(ggpubr)
library(patchwork) 
library("ggthemes")
library(magrittr)
library(ggtranscript)
library(patchwork)  

setwd("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20241229/Figure/Figure2")
isoform=read.csv("Figure2A//number.csv",check.names = F)
p1=ggplot(isoform, aes(x=id, y=number, fill=factor(type, c("novel", "known")))) +
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
  


novel=read.csv("figure2/F2C/novel.txt",header=F)
rownames(novel)=novel$V1
known=read.csv("figure2/F2C/known.txt",header=F)
rownames(known)=known$V1

com_novel=intersect(rownames(novel),rownames(gtf_total))
com_known=intersect(rownames(known),rownames(gtf_total))

gtf_known=read.table("Figure2B/known_statistics.tsv",header=T,row.names = 1)
gtf_novel=read.table("Figure2B/novel_statistics.tsv",header=T,row.names = 1)

aml_novel_length=log(gtf_novel$transcript_length)
aml_known_length=log(gtf_known$transcript_length)
# 创建数据框
df <- data.frame(length = c(aml_novel_length, aml_known_length),
                 type = c(rep("novel", length(aml_novel_length)),
                          rep("known", length(aml_known_length))))

# 绘制密度分布图
p2.1=
  ggplot(df, aes(x = length, fill=factor(type,c("novel","known")))) +
  geom_density(alpha=0.5) +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "log ( Length of isoform )", y = "Density")+
  labs(fill = "isoform type") +   theme_bw() + theme(panel.grid=element_blank())


aml_novel_orf_length=log(gtf_novel$orf_Length)
aml_known_orf_length=log(gtf_known$orf_Length)
# 创建数据框
df2 <- data.frame(length = c(aml_novel_orf_length, aml_known_orf_length),
                 type = c(rep("novel", length(aml_novel_orf_length)),
                          rep("known", length(aml_known_orf_length))))

# 绘制密度分布图
p3.2=
  ggplot(df2, aes(x = length, fill=factor(type,c("novel","known")))) +
  geom_density(alpha=0.5) +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "log ( Length of ORF )", y = "Density")+
  labs(fill = "isoform type")  +   theme_bw() + theme(panel.grid=element_blank())




aml_novel_exons=log(gtf_novel$Exon_Count)
aml_known_exons=log(gtf_known$Exon_Count)
# 创建数据框
df3 <- data.frame(length = c(aml_novel_exons, aml_known_exons),
                  type = c(rep("novel", length(aml_novel_exons)),
                           rep("known", length(aml_known_exons))))

# 绘制密度分布图
p3.3=
  ggplot(df3, aes(x = length, fill=factor(type,c("novel","known")))) +
  geom_density(alpha=0.5) +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "log ( number of exons )", y = "Density")+
  labs(fill = "isoform type")  +   theme_bw() + theme(panel.grid=element_blank())



library(ggseqlogo)

a=read.csv("F2D/aml_known_5.txt",header=F)
b=read.csv("F2D/aml_known_3.txt",header=F)
c=read.csv("F2D/aml_novel_5.txt",header = F)
d=read.csv("F2D/aml_novel_3.txt",header = F)


p4.1=
  ggplot()+
  geom_logo(b$V1,method = "bits",col_scheme = "nucleotide")+
  theme_logo()+
  annotate('rect', xmin = 0.5, xmax =1.5, 
           ymin = 0, ymax = 2, 
           alpha = .1, col='black', fill='yellow',
           linewidth = 0.5) +
  annotate('segment', x = 2.5, xend=7, y=0.9, yend=0.9, size=1) +
  annotate('text', x=4.75, y=1, label="5' end of the intron of the known isoforms") 

p4.2=
  ggplot()+
  geom_logo(a$V1,method = "bits",col_scheme = "nucleotide")+
  theme_logo()+
  annotate('rect', xmin = 8.5, xmax = 10.5, 
           ymin = 0, ymax = 2.01, 
           alpha = .1, col='black', fill='yellow',
           linewidth = 0.5) +
  annotate('segment', x = 2, xend=6, y=0.8, yend=0.8, size=1) +
  annotate('text', x=4, y=0.9, label="3' end of the intron of the known isoforms") 


p4.3=
  ggplot()+
  geom_logo(d$V1,method = "bits",col_scheme = "nucleotide")+
  theme_logo()+
  annotate('rect', xmin = 0.5, xmax =2.5, 
           ymin = 0, ymax = 4.25, 
           alpha = .1, col='black', fill='yellow',
           linewidth = 0.5) +
  annotate('segment', x = 4.3, xend=8.6, y=3.25, yend=3.25, size=1) +
  annotate('text', x=6.5, y=3.4, label="5' end of the intron of the novel isoforms") 
p4.4=
  ggplot()+
  geom_logo(c$V1,method = "bits",col_scheme = "nucleotide")+
  theme_logo()+
  annotate('rect', xmin = 8.5, xmax = 10.5, 
           ymin = 0, ymax = 1.96, 
           alpha = .1, col='black',
           fill='yellow',linewidth = 0.5) +
  annotate('segment', x = 2.5, xend=7, y=0.8, yend=0.8, size=1) +
  annotate('text', x=4.75, y=0.9, label="3' end of the intron of the novel isoforms") 




ms=read.csv("F2E/numbers2.csv",check.names = F)


ms$id <- factor(ms$id, levels = c("Jayavelu AK,et al.","11 AML"))

# plotting graph
p5=ggplot(ms, aes(fill = `isoform type`,
               y = number, x = id,bcolor = `isoform type` ))+
  geom_bar(position = "stack", stat = "identity")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "", y = "Numbers of isoforms") +
  theme_bw() +
  theme(panel.grid=element_blank())+
  theme(legend.position =c(0.9,0.75))+coord_flip()



# P6
merge=read.csv("F2F/ggtranscript.csv",header=T)
exons <- merge %>% dplyr::filter(gene_name=="CEBPA",type=="exon")
#exons$width=exons$width-1
# 创建 type 列为因子变量，并定义顺序
exons$classification<- factor(exons$classification, levels = c( "known","novel"))  # 逆序定义
exons$transcript_name
rescaled <- shorten_gaps(
  exons, 
  to_intron(exons, "transcript_name"), 
  group_var = "transcript_name"
)


# 绘图
#p6=
exons  %>% dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = reorder(transcript_id, as.numeric(type)),  
    fill = classification
  )) +
  geom_range() +
  scale_fill_manual(values = c("known" = "blue","novel" = "red")) +
  theme_classic() + 
  ylab(NULL)+
  theme(legend.position = "top",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

  



  
p3=
  ggarrange(p3.1,p3.2,p3.3,ncol = 3, common.legend = TRUE,legend="right")

p4=
  ggarrange(p4.1,p4.2,p4.3,p4.4,ncol = 2)

(p1 | p2) / p3 / (p4 | p5) / p6+plot_annotation(tag_levels = 'A')

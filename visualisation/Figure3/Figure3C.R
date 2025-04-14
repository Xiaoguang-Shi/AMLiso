library(circlize)
library(ComplexHeatmap)
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
setwd("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20241229/Figure/Figure3")


phenotype = read.csv("aml-related-gene.csv",header=T,row.names = 1)

exp=read.csv("tmp.tpm.csv",header=T)
com=intersect(rownames(phenotype),colnames(exp))
exp=exp[,com]

exp_scale=t(scale(t(exp), center=TRUE, scale=TRUE))



age_col=colorRamp2(c(17,50,75), c("blue", "white", "red"))
BMBlast_col=colorRamp2(c(0,60,100), c("blue", "white", "red"))
WBC_col=colorRamp2(c(0.5,12,20,250), c("blue", "white", "red","red"))
PTL_col=colorRamp2(c(0,50,100,1000), c("blue", "white", "red","red"))
RBC_col=colorRamp2(c(1,3,5), c("blue", "white", "red"))
Hb_col=colorRamp2(c(40,100,150), c("blue", "white", "red"))
sv_col = structure(names = c("0","1"), c( "white", "black"))
gender_col = structure(names = c("male", "female"), c( "red", "blue"))
ENL_col = structure(names = c("NA","1","2","3"), c( "grey", "green","yellow","red"))
FAB_col = structure(names = c("AML", "M1","M2","M3","M4","M5","M6"), c( "grey", "red","purple","blue","turquoise","yellowgreen","orange"))
FLT3ITD_col=structure(names = c("NO","YES","NO" ), c( "white", "red","white"))
fusion_col=structure(names = c("0","1"), c( "white", "red"))
mu_col= structure(names = c("In_Frame_Ins",
                            "In_Frame_Del",
                            "Frame_Shift_Ins",
                            "Frame_Shift_Del",
                            "Nonsense_Mutation",
                            "Missense_Mutation",
                            "Splice_Site",
                            "Translation_Start_Site",
                            "NULL",
                            "Multi_Hit"), 
                  c( "brown",
                     "orange", 
                     "purple",
                     "steelblue2",
                     "orangered2",
                     "darkgreen",
                     "yellow",
                     "brown",
                     "white",
                     "black"))
mu_number_col = structure(names = c("0","1","2","3","4"), c( "white", "darkslategray2","darkslategray3","darkslategray4","darkslategray"))



top_annotation = HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = c("#CE332A","#4A7DB3","#68AD57","#8E529E","#EE8632","#F4D03F","#9B5A33","#D262B9","#9A999A","#060964")),
    labels = c("C1","C2","C3","C4","C5","C6","C7","C8"),
    labels_gp = gpar(col = "white", fontsize = 12,fontface = "bold")))
#通用代码

ha = HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = c("#E64B10","#4E3184","#88C63C","#FFC20D","#FF7F00","#347DFF","#549CF8","#8781FF")),
    labels = c("C1","C2","C3","C4","C5","C6","C7","C8"),
    labels_gp = gpar(col = "white", fontsize = 12,fontface = "bold")),
  age = phenotype[[2]], 
  gender_cluster = phenotype[[1]],
  BMBlast_cluster = phenotype[[3]],
  WBC_cluster = phenotype[[4]],
  RBC_cluster = phenotype[[5]],
  Hb_cluster = phenotype[[6]],
  PTL_cluster = phenotype[[7]],
  FAB_cluster=phenotype[[13]],
  ELN_cluster=phenotype[[12]],
  sv_8=phenotype[[8]],
  sv_5=phenotype[[9]],
  sv_7=phenotype[[10]],
  sv_17=phenotype[[11]],
  FLT3ITD=phenotype[[60]],
  PMLRARA=phenotype[[14]],
  RUNX1RUNX1T1=phenotype[[15]],
  CBFBMYH11=phenotype[[16]],
  NUP98=phenotype[[17]],
  KMT2A=phenotype[[18]],
  CEBPA=phenotype[[21]],
  NPM1=phenotype[[22]],
  NRAS=phenotype[[23]],
  TTN=phenotype[[24]],
  TET2=phenotype[[25]],
  DNMT3A=phenotype[[26]],
  FLT3=phenotype[[27]],
  GATA2=phenotype[[28]],
  WT1=phenotype[[29]],
  IDH2=phenotype[[30]],
  SZT2=phenotype[[31]],
  IDH1=phenotype[[32]],
  PTPN11=phenotype[[33]],
  SMC1A=phenotype[[34]],
  DHX15=phenotype[[35]],
  SRCAP=phenotype[[36]],
  RUNX1=phenotype[[37]],
  BCOR=phenotype[[38]],
  BCORL1=phenotype[[39]],
  IKZF1=phenotype[[40]],
  SERINC2=phenotype[[41]],
  KMT2C=phenotype[[42]],
  PKD1=phenotype[[43]],
  USP9X=phenotype[[44]],
  ASXL2=phenotype[[45]],
  CREBBP=phenotype[[46]],
  MADD=phenotype[[47]],
  MAP4=phenotype[[48]],
  STAG2=phenotype[[49]],
  PHRF1=phenotype[[50]],
  RAD21=phenotype[[51]],
  SETD2=phenotype[[52]],
  STAG2=phenotype[[53]],
  SMC3=phenotype[[54]],
  MAPK8IP3=phenotype[[55]],
  SOS1=phenotype[[56]],
  ANKLE1=phenotype[[57]],
  CHD9=phenotype[[58]],
  VariantClassification=phenotype[[20]],
  col = list( age=age_col,
              gender_cluster =gender_col ,
              BMBlast_cluster = BMBlast_col ,
              WBC_cluster = WBC_col,
              RBC_cluster = RBC_col,
              PTL_cluster = PTL_col,
              ELN_cluster=ENL_col ,
              Hb_cluster=Hb_col ,
              sv_8=sv_col,
              sv_5=sv_col,
              sv_7=sv_col,
              sv_17=sv_col,
              FAB_cluster=FAB_col,
              PMLRARA=fusion_col,
              RUNX1RUNX1T1=fusion_col,
              CBFBMYH11=fusion_col,
              NUP98=fusion_col,
              KMT2A=fusion_col,
              fusion=fusion_col,
              CEBPA=mu_col,
              NPM1=mu_col,
              GATA2=mu_col,
              TET2=mu_col,
              FLT3=mu_col,
              FLT3ITD=FLT3ITD_col,
              WT1=mu_col,
              DNMT3A=mu_col,
              RUNX1=mu_col,
              BCOR=mu_col,
              PTPN11=mu_col,
              SMC1A=mu_col,
              IDH2=mu_col,
              NRAS=mu_col,
              TTN=mu_col,
              USP9X=mu_col,
              ANKLE1=mu_col,
              BCORL1=mu_col,
              DHX15=mu_col,
              IKZF1=mu_col,
              KMT2C=mu_col,
              STAG2=mu_col,
              SZT2=mu_col,
              ASXL2=mu_col,
              CHD9=mu_col,
              CREBBP=mu_col,
              GOLGA8R=mu_col,
              IDH1=mu_col,
              MADD=mu_col,
              MAP4=mu_col,
              MAPK8IP3=mu_col,
              PHRF1=mu_col,
              PKD1=mu_col,
              RAD21=mu_col,
              SERINC2=mu_col,
              SETD2=mu_col,
              SMC3=mu_col,
              SOS1=mu_col,
              SRCAP=mu_col,
              WASHC1=mu_col,
              VariantClassification=mu_col),
  na_col = "grey", border = F,
  show_legend = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
  show_annotation_name = T,
  annotation_legend_param = list(
    age = list(title = "age at diagnosis"),
    gender_cluster = list(title = "Gender"),
    BMBlast_cluster = list(title = "BM Blast(%)"),
    WBC_cluster = list(title = "WBC (10^9/L)"),
    RBC_cluster = list(title = "RBC (10^12/L)"),
    Hb_cluster = list(title = "Hb (g/L)"),
    PTL_cluster = list(title = "PTL (10^9/L)"),
    sv_8=list(title = "SV"),
    ELN_cluster=list(title = "ELN"),
    FAB_cluster=list(title ="FAB subtype"),
    VariantClassification=list(title ="Mut"),
    FLT3ITD=list(title ="FLT3ITD")
  ))

col_fun = colorRamp2(c(-2,0,2), c("blue", "white", "red"))

ht_list=
  Heatmap(exp,
          cluster_rows = F,
          cluster_columns = F,
          col= col_fun,
          column_split = phenotype$NMF,
          show_row_names = F,
          row_names_gp = gpar(fontsize = 10),
          show_column_names =F,
          show_heatmap_legend = T,
          row_title = NULL,
          column_title = NULL,
          top_annotation = ha,
          border = TRUE,
          column_gap = unit(0.5, "mm"),
          row_gap = unit(0.5, "mm"),
          row_title_gp = gpar(col = "#FFFFFF00"))



draw(ht_list, 
     annotation_legend_side = "right", heatmap_legend_side = "right")





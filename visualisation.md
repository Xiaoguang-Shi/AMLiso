## visualisation

## Figure2

```{r}
library(ggplot2)
# read data

isoform=read.csv("number.csv",check.names = F)

# Create the plot

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

## Figure3
### Figure 3A
```{r}
# Load required libraries
library(NMF)
library(bigmemory)

# Read gene expression data (genes in rows, samples in columns)
data <- read.table("aml_up_tpm.tsv", header=T, sep='\t', row.names=1)

# Feature selection: Select top 5000 genes with highest variance
row_variances <- apply(data, 1, var)
top_indices <- order(row_variances, decreasing = TRUE)[1:5000]
data <- data[top_indices, ]

# Data preprocessing: Apply log2 transformation (add pseudocount to avoid zeros)
data <- log2(data + 1.0001)

# Run NMF with rank range 2-15 for model selection
rank <- 2:15
nmf_res <- nmf(data, rank=rank, nrun=10000, method="brunet", seed=10073, .options='P40')

# Save workspace and plot cophenetic correlation
save.image(file="mads_rank2_15_nrun10000_normal.Rdata")
pdf(file="cophenetic_2_15_run10000_normal.pdf", width=8, height=7, onefile=F)
plot(nmf_res)
dev.off()

# Final model with selected rank (8)
rank <- 8
nmf_res <- nmf(data, rank=rank, nrun=10000, method="brunet", seed=123456, .options='P40')

# Save clustering results
Cluster <- predict(nmf_res)
Cluster <- as.data.frame(Cluster)
Cluster$Cluster <- paste0("C", Cluster$Cluster)
clusterOut <- rbind(ID=colnames(Cluster), Cluster)
write.table(clusterOut, file="cluster_rank8-normal-5k.txt", sep="\t", quote=F, col.names=F)

# Generate consensus matrix heatmap
pdf(file="heatmap_rank8.pdf", width=6, height=6, onefile=F)
consensusmap(nmf_res,
             annRow=NA,
             annCol=NA,
             main="Consensus matrix",
             info=FALSE)
dev.off()
```

## Figure 3B

```{r}
# Load required libraries
library(ggplot2)
library(ggalluvial)



# Read cluster data with header and row names
cluster <- read.csv("cluster.csv", 
                    header = TRUE, 
                    row.names = 1)

# Prepare data for alluvial plot: select relevant columns and create link variable
link <- cluster[, c("NMF", "WHO", "class_max")]  # "class_max" kept as comment for reference
link$link <- 1  # Create dummy variable for melting

# Reshape data to long format using reshape package
link <- reshape::melt(link, id = 'link')

# Create flow variable for alluvial diagram
variable <- summary(link$variable)
link$flow <- rep(1:variable[1], length(variable))

# Define factor levels for ordered display
link$value <- factor(link$value, levels = c(
  "PML::RARA", "RUNX1::RUNX1T1", "CEBPA_mutation","NPM1_mutation","CBFB::MYH11",
  "NUP98_rearrangement","KMT2A_rearrangement","DEK::NUP214","BCR::ABL1","MR","others",
  "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
  "G1","G3","G4","G7","G8","G2","G6","G5"))

# Define color scheme for categories
colors <- c(
  "PML::RARA" = "#E64B10", "C1" = "#E64B10", "G1" = "#E64B10",
  "RUNX1::RUNX1T1" = "#4E3184", "C2" = "#4E3184", "G3" = "#4E3184",
  "CEBPA_mutation" = "#88C63C", "C3" = "#88C63C", "G4" = "#88C63C",
  "NPM1_mutation" = "#ffc20d", "C4" = "#ffc20d", "G7" = "#ffc20d", "G8" = "#ffc20d",
  "C6" = "#347DFF", "C7" = "#549CF8", "G5" = "#549CF8",
  "G6" = "#6981ff", "C8" = "#8781FF",
  "CBFB::MYH11" = "#DD4D42", "C5" = "#FF7F00", "MR" = "#105E7F",
  "others" = "#295EA4", "NUP98_rearrangement" = "#BBDD8E",
  "KMT2A_rearrangement" = "#F09D9E", "DEK::NUP214" = "#E0D6E4",
  "BCR::ABL1" = "#360A47"
)

# Create alluvial plot
ggplot(link, aes(x = variable,
                 y = link,
                 stratum = value,
                 alluvium = flow,
                 fill = value)) +
  geom_stratum() +  # Create stratum rectangles
  scale_fill_manual(values = colors) +  # Apply color mapping
  geom_flow(aes.flow = 'forward') +  # Create flow curves
  geom_text(stat = 'stratum',  # Add labels
            infer.label = TRUE,
            size = 2.5) +
  theme_minimal() +  # Use minimal theme
  geom_alluvium() +  # Add alluvial flows
  theme(
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove background grid
    legend.position = 'none'         # Remove legend
  )
```

## Figure 3C

```{r}
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

```


### Figure 3D,F
```{r}
# Load required packages
library(survival)       # Survival analysis
library(survminer)      # Survival visualization
library(forestmodel)    # Forest plot visualization
library(dplyr)          # Data manipulation

# Set working directory and load data

surv_data <- read.csv("surv_data.csv", header = TRUE, row.names = 1)

# Convert survival time from days to months
surv_data$futime <- surv_data$futime / 30

# Define color palette for clusters
cluster_palette <- c("#E64B10", "#4E3184", "#88C63C", "#FFC20D",
                     "#FF7F00", "#347DFF", "#549CF8", "#8781FF")

#-----------------------------
# Kaplan-Meier Survival Analysis
#-----------------------------

# Fit survival curve for all clusters
surv_fit_all <- survfit(Surv(futime, fustat) ~ NMF, data = surv_data)

# Generate survival plot for all data
ggsurvplot(
  surv_fit_all,
  data = surv_data,
  conf.int = FALSE,
  pval = TRUE,                     # Display log-rank p-value
  pval.size = 6,                   # P-value text size
  legend = "none",                 # Hide legend
  palette = cluster_palette,       # Custom color scheme
  xlab = "Time (Months)",          # X-axis label
  break.time.by = 10,              # Time axis breaks
  surv.median.line = "hv",         # Show median survival lines
  risk.table = TRUE,               # Display risk table
  risk.table.height = 0.3          # Risk table height ratio
)

#-----------------------------
# Cox Proportional Hazards Model
#-----------------------------

# Prepare data for Cox regression
cox_data <- surv_data %>%
  filter(NMF != "C7") %>%          # Exclude cluster C7
  mutate(
    age = factor(age),
    gender = factor(gender),
    BMT = factor(BMT),
    WBC = as.numeric(WBC),
    RBC = as.numeric(RBC),
    PLT = as.numeric(PTL),
    BMBlast = as.numeric(BMBlast),
    ELN = factor(ELN),
    NMF = relevel(factor(NMF), ref = "C1")  # Set reference cluster
  )

# Fit Cox proportional hazards model
cox_model <- coxph(
  Surv(futime, fustat) ~ 
    age + gender + BMBlast + WBC + RBC + PLT + NMF,
  data = cox_data
)

# Visualize results with forest plot
forest_model(
  cox_model,
  format_options = forest_model_format_options(
    colour = "black",     # Font color
    point_size = 3,       # Hazard ratio point size
    banded = TRUE         # Alternate row colors
  ),
  factor_separate_line = FALSE
)

#-----------------------------
# Subgroup Analysis (BMT=0)
#-----------------------------

# Filter data for BMT=0 subgroup
bmt0_data <- surv_data %>% 
  filter(BMT == "0", NMF != "C7", WHO_subtype != "APL") %>%
  mutate(
    NMF = droplevels(factor(NMF)),  # Remove unused factor levels
    NMF = relevel(NMF, ref = "C4")   # New reference cluster
  )

# Fit subgroup Cox model
cox_bmt0 <- coxph(
  Surv(futime, fustat) ~ 
    age + gender + BMBlast + WBC + RBC + PLT + ELN + NMF,
  data = bmt0_data
)

# Generate subgroup forest plot
forest_model(
  cox_bmt0,
  format_options = forest_model_format_options(
    colour = "black",
    point_size = 3,
    banded = TRUE
  ),
  factor_separate_line = FALSE
)
```


### Figure 3F 

```{r}
library(circlize)
library(ComplexHeatmap)

# Set working directory and load data
setwd("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20241229/Figure/Figure3")

# Load and preprocess phenotype data
phenotype <- read.csv("aml-related-gene.csv", header = TRUE, row.names = 1)

# Load and process clinical expression data
exp <- read.csv("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20241229/clinical/clinical_20250315.csv", 
                header = TRUE)
rownames(exp) <- exp$住院号

# Data alignment and preprocessing
common_samples <- intersect(rownames(phenotype), rownames(exp))
exp <- exp[common_samples, ]

# Prepare data matrices
prepare_matrix <- function(data, cols, scale_type = "column") {
  mat <- as.matrix(data[, cols])
  mat <- t(scale(mat, center = TRUE, scale = TRUE))
  return(mat)
}

BM_matrix <- prepare_matrix(exp, 4:31)
CIBERSORT_matrix <- prepare_matrix(exp, 32:54)

# Annotation configuration
create_top_annotation <- function() {
  HeatmapAnnotation(
    cluster = anno_block(
      gp = gpar(fill = c("#E64B10","#4E3184","#88C63C","#FFC20D",
                         "#FF7F00","#347DFF","#549CF8","#8781FF")),
      labels = c("C1","C2","C3","C4","C5","C6","C7","C8"),
      labels_gp = gpar(col = "white", fontsize = 12, fontface = "bold")
    )
  )
}

# Heatmap visualization function
create_heatmap <- function(matrix, title = "") {
  Heatmap(
    matrix,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("dodgerblue", "white", "red")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = phenotype$NMF,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
    top_annotation = create_top_annotation(),
    column_title = title,
    border = TRUE,
    column_gap = unit(0.5, "mm"),
    row_gap = unit(0.5, "mm"),
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      legend_height = unit(4, "cm")
    )
  )
}

# Generate and save heatmaps
pdf("combined_heatmaps.pdf", width = 16, height = 10)
draw(create_heatmap(BM_matrix, "Bone Marrow Markers") %v% 
       create_heatmap(CIBERSORT_matrix, "Immune Cell Composition"))
dev.off()
```

## Figure 4

```{r}
#-------------------------
# 1. Load Required Libraries
#-------------------------
library(dplyr)           # Data manipulation
library(ComplexHeatmap)  # Advanced heatmap visualization
library(ggplot2)         # Base plotting system
library(ggpubr)          # Publication-ready ggplot extensions
library(glmnet)          # LASSO regression implementation
library(survminer)       # Survival analysis visualization
library(survival)        # Survival analysis core functions
library(tidyverse)       # Data science ecosystem
library(RColorBrewer)    # Color palette management
library(ROCR)            # ROC curve analysis
library(forestmodel)     # Forest plot visualization

#-------------------------
# 2. Data Preparation
#-------------------------
setwd("/Volumes/sandisk/neoantigen/analysis/RNA_seq/featureCounts/175/20240531/lasso")

# Load expression data (genes as rows, samples as columns)
aml <- read.csv("sigGenes_aml_total.csv", header = TRUE, row.names = 1)

# Load clinical survival data
surv <- read.csv("time_nonAPL.csv", header = TRUE, row.names = 1)

# Align datasets by sample IDs
aml <- aml[rownames(surv), ]

#-------------------------
# 3. LASSO-Cox Modeling
#-------------------------
# Prepare input matrices
x <- as.matrix(aml)  # Gene expression matrix
y <- Surv(surv$futime, surv$fustat)  # Survival object

# Fit base LASSO model
set.seed(1234)
fit <- glmnet(x, y, 
              family = "cox",
              maxit = 1e+05)

# Cross-validated model tuning
set.seed(1234)
cvfit <- cv.glmnet(x, y, family = "cox")

# Final model with selected lambda
set.seed(1234)
lasso_best <- glmnet(x = x, y = y,
                     lambda = 0.24,  # Manually selected value
                     family = "cox",
                     maxit = 1e+05)

#-------------------------
# 4. Feature Extraction
#-------------------------
# Extract non-zero coefficients
gene.coef <- as.matrix(round(coef(lasso_best), 4))
coef.lasso <- gene.coef[gene.coef[, 1] != 0, ]
coef.lasso <- data.frame(genes = names(coef.lasso),
                         coefficient = as.numeric(coef.lasso)) %>%
  arrange(-coefficient)

# Save coefficient results
write.csv(coef.lasso, file = "coef.lasso.csv")

#-------------------------
# 5. Risk Score Calculation
#-------------------------
y1 <- as.data.frame(y) %>%
  mutate(LassoScore = predict(lasso_best, newx = x),
         ID = rownames(y)) %>%
  mutate(LASSO.level = ifelse(LassoScore > median(LassoScore),
                              "High", "Low"),
         LASSO.level = factor(LASSO.level, levels = c("Low", "High")))

#-------------------------
# 6. Visualization Modules
#-------------------------
# 6A. Feature Coefficient Plot
ggplot(coef.lasso, aes(reorder(genes, coefficient), coefficient, 
                       fill = coefficient > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  coord_flip() +
  theme_minimal() +
  labs(y = "Regression Coefficient", x = "")

# 6B. Risk Score Distribution
ggscatter(y1, x = "ID", y = "LassoScore", 
          color = "LASSO.level",
          palette = c("blue", "red")) +
  geom_hline(yintercept = median(y1$LassoScore), linetype = 2) +
  theme(axis.text.x = element_blank())

#-------------------------
# 7. Heatmap Visualization
#-------------------------
# Prepare heatmap data
heat.data <- t(scale(t(x[rev(coef.lasso$genes), rownames(y1)])))

# Cluster-based heatmap
Heatmap(heat.data,
        column_split = y1$NMF_cluster,
        top_annotation = HeatmapAnnotation(
          cluster = anno_block(
            gp = gpar(fill = c("#F5B041","#58D68D","#5DADE2",
                               "#8E44AD","#192466","#A2CB12","pink"))
          )
        ))

# Risk score heatmap
Heatmap(heat.data,
        column_split = y1$LASSO.level,
        top_annotation = HeatmapAnnotation(
          cluster = anno_block(
            gp = gpar(fill = c("blue", "red"))
          )
        ))

#-------------------------
# 8. Model Evaluation
#-------------------------
# ROC Curve Analysis
lasso.prob <- predict(cvfit, newx = x, s = c(cvfit$lambda.min, cvfit$lambda.1se))
pred_min <- prediction(lasso.prob[,1], y1$status)
perf_min <- performance(pred_min, "tpr", "fpr")

# Plot ROC curves
plot(perf_min, col = "red", main = "ROC Curves")
plot(performance(prediction(lasso.prob[,2], y1$status), "tpr", "fpr"), 
     col = "blue", add = TRUE)

#-------------------------
# 9. Survival Analysis
#-------------------------
# Survival curve plotting
fit <- survfit(Surv(futime, fustat) ~ LASSO.level, data = y1)
ggsurvplot(fit, data = y1,
           palette = c("#00599F", "#d80700"),
           risk.table = TRUE,
           pval = TRUE)

#-------------------------
# 10. TCGA Validation
#-------------------------
# Load external validation data
tcga <- read.csv("TCGA_tpm.csv")  # TCGA expression data
tcga_clinical <- read.csv("TCGA_clinical.csv")  # TCGA survival data

# Predict risk scores
tcga_risk <- predict(lasso_best, newx = as.matrix(tcga), s = 0.24)

# Survival analysis
fit_tcga <- survfit(Surv(futime_months, fustat) ~ LASSO.level,
                    data = tcga_clinical)
ggsurvplot(fit_tcga, palette = c("#00599F", "#d80700"))


# -------------------------------
# 11. ELN Risk Transition Visualization
# -------------------------------
library(ggplot2)
library(ggalluvial)

# Load and preprocess ELN data
eln_data <- read.csv("eln.csv", row.names = 1, header = TRUE)

# Prepare alluvial plot data structure
eln_flow <- eln_data %>% 
  select(ELN, Adjusted_ELN) %>%
  mutate(PatientID = row.names(.)) %>%
  pivot_longer(-PatientID, names_to = "Assessment", values_to = "Risk")

# Create categorical mapping
risk_levels <- c("1" = "Favorable", "2" = "Intermediate", "3" = "Adverse",
                 "Low" = "Low Risk", "High" = "High Risk")

# Visualization parameters
eln_colors <- c(
  "Favorable" = "#B9CEE3",       # Blue
  "Intermediate" = "#FFE5C1",    # Beige
  "Adverse" = "#F6BEB9",         # Red
  "Low Risk" = "#B9CEE3",        # Matching blue
  "High Risk" = "#F6BEB9"        # Matching red
)

# Generate alluvial plot
ggplot(eln_flow, aes(x = Assessment, 
                     stratum = Risk,
                     alluvium = PatientID,
                     fill = Risk,
                     label = Risk)) +
  geom_flow(alpha = 0.6, 
            curve_type = "sigmoid",
            aes(fill = after_scale(colorspace::lighten(fill, 0.8)))) +
  geom_stratum(width = 0.4) +
  geom_text(stat = "stratum", size = 3, color = "black") +
  scale_fill_manual(values = eln_colors,
                    labels = risk_levels) +
  labs(title = "ELN Risk Category Transitions",
       x = "Assessment Stage",
       y = "Patient Count") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# -------------------------------
# 12. Survival Analysis Visualization
# -------------------------------
library(survival)
library(survminer)

# Prepare survival data
eln_surv <- eln_data %>%
  mutate(
    futime = futime / 30,  # Convert days to months
    ELN = factor(ELN, 
                 levels = c(1:3),
                 labels = c("Favorable", "Intermediate", "Adverse"))
  )

# Fit survival curves
surv_object <- survfit(Surv(futime, fustat) ~ ELN, data = eln_surv)

# Generate survival plot
ggsurvplot(
  surv_object,
  data = eln_surv,
  conf.int = TRUE,
  pval = TRUE,
  pval.size = 5,
  xlab = "Time (Months)",
  break.time.by = 6,
  legend.title = "ELN Risk",
  legend.labs = c("Favorable", "Intermediate", "Adverse"),
  surv.median.line = "hv",
  risk.table = TRUE,
  risk.table.height = 0.3,
  palette = c("#00599F", "#F5B041", "#D80700"),
  ggtheme = theme_minimal()
)

```


## Figure 5

```{r}
# ================================================================================
# CWGCNA Analysis Pipeline - Enhanced Version with Robust Error Handling
# ================================================================================

# ------------------------------
# 1. Environment Configuration
# ------------------------------
library(CWGCNA)
library(WGCNA)


# ------------------------------
# 2. Data Loading & Preprocessing
# ------------------------------
# Read differential expression data
dif <- read.table("AML_vs_CD34_up-normal.tsv", header=TRUE, row.names=1)

# Load normalized expression matrix (VST transformed)
exp <- read.table("exp.tsv",
                  header=TRUE, sep='\t', row.names=1, check.names=FALSE)
exp <- exp[rownames(dif), ]

# Load clinical metadata
clinical <- read.csv("clinical.csv", header=TRUE, row.names=1, check.names=FALSE)
com <- intersect(rownames(clinical), colnames(exp))

# Align datasets
clinical <- clinical[com, ]
exp <- exp[, com]

# ------------------------------
# 3. Sample Anonymization
# ------------------------------
colnames(exp) <- paste0("sample", 1:175)
rownames(clinical) <- paste0("sample", 1:175)
clinical$sampleid <- rownames(clinical)
clinical <- clinical[, c(ncol(clinical), 1:(ncol(clinical)-1))]  # Reorder columns

# ------------------------------
# 4. Initial Variance Analysis
# ------------------------------
anovares <- featuresampling(
  betas = exp,
  pddat = clinical,
  anova = TRUE,
  plotannovares = TRUE,
  featuretype = "probe",
  plottitlesuffix = "Bone marrow",
  threads = 8
)

# ------------------------------
# 5. Batch CWGCNA Analysis
# ------------------------------
bm <- colnames(clinical)[4:ncol(clinical)]

# Batch processing with error logging
analyze_phenotype <- function(id) {
  tryCatch({
    pdf(file=paste0(id,".pdf"), width=10, height=8)
    
    result <- diffwgcna(
      dat = exp,
      pddat = clinical,
      responsevarname = id,
      confoundings = c("age", "gender"),
      featuretype = "gene",
      topvaricancetype = "sd",
      topvaricance = 20000,
      powers = seq(1, 20, 1),
      rsqcutline = 0.85,
      minclustersize = 30,
      mediation = TRUE,
      mergecutheight = 0.25,
      maxblocksize = 5000,
      topn = 1,
      plot = TRUE,
      titleprefix = "AML",
      labelnum = 8,
      titlesize = 17,
      textsize = 16,
      annotextsize = 5,
      pvalcolname = "adj.P.Val",
      pvalcutoff = 0.05,
      isbetaval = TRUE,
      absxcutoff = 0,
      diffcutoff = 0,
      threads = 8
    )
    
    # Save outputs
    save_results(id, result)
    dev.off()
    return(TRUE)
  }, error = function(e) {
    message(paste("Error processing", id, ":", e$message))
    dev.off()
    return(FALSE)
  })
}

save_results <- function(id, res) {
  write.csv(res$mediationres, paste0(id, ".mediationres.csv"))
  write.csv(res$limmares, paste0(id, ".limmares.csv"))
  write.csv(res$wgcnares$tom, paste0(id, ".wgcnares.csv"))
  write.csv(res$sftpowers$fitIndices, paste0(id, ".fitIndices.csv"))
  write.csv(res$melimmares, paste0(id, ".melimmares.csv"))
  save.image(paste0(id,".Rdata"))
}

# Execute analysis
sapply(bm[1], analyze_phenotype)  # Initial test with first variable
sapply(bm[30:50], analyze_phenotype)  # Batch processing

# ------------------------------
# 6. Module-Trait Analysis
# ------------------------------
load("BM.myeloblasts.Rdata")
Mes <- cwgcnares$wgcnares$mergedmelist
nSamples <- nrow(t(exp))
clinical_2 <- clinical[, -c(1:3)]

# Correlation analysis
moduleTraitCor <- cor(Mes, clinical_2, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Matrix reorganization
ordered_index <- order(as.numeric(gsub("ME", "", rownames(moduleTraitCor))))
moduleTraitCor <- moduleTraitCor[ordered_index, ]
moduleTraitPvalue <- moduleTraitPvalue[ordered_index, ]

# Significance filtering
significant_modules <- moduleTraitPvalue <= 0.05
filteredCor <- moduleTraitCor[rowSums(significant_modules) > 0, 
                              colSums(significant_modules) > 0]

# Annotation generation
create_annotations <- function(pvals) {
  cut(pvals, 
      breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
      labels=c("***", "**", "*", ""))
}

# Visualization
pdf("10.Module-trait_relationships_CIBERSORT_pvalue.pdf", height=12, width=15)
par(mar=c(25, 10, 3, 3))
labeledHeatmap(
  Matrix = filteredCor,
  xLabels = colnames(filteredCor),
  yLabels = rownames(filteredCor),
  ySymbols = rownames(filteredCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(100),
  textMatrix = create_annotations(moduleTraitPvalue[row.names(filteredCor), 
                                                    colnames(filteredCor)]),
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-0.5, 0.5),
  main = "Module-Trait Relationships with Significance Annotation"
)
dev.off()
```


## Figure 6

```{r}
# ==============================================================================
# Mutation-Expression Correlation Analysis Pipeline
# ==============================================================================

# ------------------------------
# 1. Load Required Libraries
# ------------------------------
library(dplyr)    # Data manipulation
library(tidyr)    # Data reshaping
library(ggplot2)  # Data visualization
library(reshape2) # Matrix operations
library(corrplot) # Correlation visualization

# ------------------------------
# 2. Data Loading & Initialization
# ------------------------------


# Load mutation data (Variant Allele Frequency)
mutation_data <- read.table("mutation.txt", header = TRUE)
colnames(mutation_data) <- c("sample_id", "gene", "vaf")

# ------------------------------
# 3. Mutation Data Processing
# ------------------------------
# Calculate mean VAF per gene-sample pair
mean_vaf_matrix <- mutation_data %>%
  group_by(sample_id, gene) %>%
  summarise(vaf = mean(vaf), .groups = "keep")

# Calculate max VAF per gene-sample pair
max_vaf_matrix <- mutation_data %>%
  group_by(sample_id, gene) %>%
  summarise(vaf = max(vaf), .groups = "keep")

# Convert to wide format matrices
mean_vaf_matrix <- mean_vaf_matrix %>%
  pivot_wider(names_from = sample_id, values_from = vaf, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

max_vaf_matrix <- max_vaf_matrix %>%
  pivot_wider(names_from = sample_id, values_from = vaf, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

# ------------------------------
# 4. Gene Frequency Filtering
# ------------------------------
gene_frequency <- read.table("gene.stat.tsv", header = FALSE, sep = "\t")
frequent_genes <- gene_frequency[gene_frequency$V1 > 4, ]

# Filter matrices to keep only frequent genes
mean_vaf_matrix <- mean_vaf_matrix[frequent_genes$V2, ]
max_vaf_matrix <- max_vaf_matrix[frequent_genes$V2, ]

# ------------------------------
# 5. Expression Data Processing
# ------------------------------
expression_matrix <- read.table("me.tsv", header = TRUE, row.names = 1)

# Format expression matrix column names
colnames(expression_matrix) <- gsub("ME(\\d+)", "\\1", colnames(expression_matrix)) %>%
  sprintf("ME%02d", as.integer(.) + 1)

expression_matrix <- t(expression_matrix)

# ------------------------------
# 6. Data Alignment
# ------------------------------
common_samples <- intersect(colnames(mean_vaf_matrix), colnames(expression_matrix))
expression_matrix <- expression_matrix[, common_samples]

# Transpose matrices for correlation analysis
expression_samples_as_rows <- t(expression_matrix) %>% as.data.frame()
mutation_samples_as_rows <- t(mean_vaf_matrix) %>% as.data.frame()

# ------------------------------
# 7. Correlation Analysis
# ------------------------------
n_genes <- ncol(mutation_samples_as_rows)
n_modules <- ncol(expression_samples_as_rows)

# Initialize result matrices
correlation_matrix <- matrix(0, nrow = n_modules, ncol = n_genes)
p_value_matrix <- matrix(NA, nrow = n_modules, ncol = n_genes)

# Calculate pairwise correlations
for (module_idx in 1:n_modules) {
  for (gene_idx in 1:n_genes) {
    test_result <- cor.test(expression_samples_as_rows[, module_idx], 
                            mutation_samples_as_rows[, gene_idx],
                            method = "pearson")
    correlation_matrix[module_idx, gene_idx] <- test_result$estimate
    p_value_matrix[module_idx, gene_idx] <- test_result$p.value
  }
}

# ------------------------------
# 8. Multiple Testing Correction
# ------------------------------
fdr_matrix <- apply(p_value_matrix, 2, p.adjust, method = "bonferroni")

# ------------------------------
# 9. Significant Results Extraction
# ------------------------------
significant_results <- which(fdr_matrix < 0.05, arr.ind = TRUE)
output_data <- data.frame(
  Mutation = rownames(mean_vaf_matrix)[significant_results[, 2]],
  Transcript = rownames(expression_matrix)[significant_results[, 1]],
  Correlation = correlation_matrix[significant_results],
  P_value = p_value_matrix[significant_results],
  FDR = fdr_matrix[significant_results]
)

write.csv(output_data, "mutation_mean_cor_fdr0.05.csv")

# ------------------------------
# 10. Heatmap Visualization
# ------------------------------
correlation_heatmap_data <- acast(output_data, Mutation ~ Transcript, 
                                  value.var = "Correlation")
correlation_heatmap_data[is.na(correlation_heatmap_data)] <- 0

ggplot(melt(correlation_heatmap_data), aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = "Gene Expression Modules", 
       y = "Mutation Genes",
       title = "Mutation-Expression Correlation Landscape") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------
# 11. Correlation Matrix Plot
# ------------------------------
correlation_network <- xtabs(Correlation ~ Mutation + Transcript, 
                             data = output_data)

corrplot(correlation_network, 
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black",
         tl.srt = 45,
         title = "Significant Mutation-Expression Correlations",
         mar = c(2, 2, 3, 2))


```

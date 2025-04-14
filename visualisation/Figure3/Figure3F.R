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

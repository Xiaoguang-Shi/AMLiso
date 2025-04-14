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
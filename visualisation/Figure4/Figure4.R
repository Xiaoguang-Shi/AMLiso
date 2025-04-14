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

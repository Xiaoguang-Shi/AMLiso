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
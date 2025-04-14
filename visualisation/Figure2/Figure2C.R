# Visualization of Transcript Feature Comparisons
# -------------------------------------------------
# Load required packages
library(ggplot2)       # Core plotting framework
library(scales)        # Axis scaling and formatting
library(dplyr)         # Data manipulation
library(ggpubr)        # Publication-ready plots
library(gghalves)      # Half-violin plots (original implementation)

# Data Preparation
# -------------------------------------------------
# Import dataset containing transcript features
known <- read.delim("known.stat.txt", na.strings = "NA", check.names = FALSE)
novel <- read.delim("novel.stat.txt", na.strings = "NA", check.names = FALSE)

# Data cleaning pipeline
combined <- bind_rows(
  known %>% filter(`CDS Total Length` > 0) %>% mutate(Group = "Known"),
  novel %>% filter(`CDS Total Length` > 0) %>% mutate(Group = "Novel")
) %>% mutate(
  Group = factor(Group, levels = c("Known", "Novel")),
  across(c(`Exon_Total_Length`, `CDS Total Length`), ~ .x + 1) # Avoid log(0)
)

# Visualization Function
# -------------------------------------------------
#' Generate comparative violin-boxplot visualization
#' @param data Combined dataset
#' @param y_var Feature variable to visualize
#' @param title Plot title
#' @param y_lab Y-axis label
create_violin_plot <- function(data, y_var, title, y_lab) {
  
  # Non-parametric hypothesis testing
  test_res <- wilcox.test(reformulate("Group", y_var), data)
  
  ggplot(data, aes(x = Group, y = .data[[y_var]], fill = Group)) +
    ggpubr::geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    scale_fill_manual(values = c("#866EFE", "#FF7259")) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(
      title = title,
      subtitle = paste("Mann-Whitney U Test, P =", 
                       format.pval(test_res$p.value, digits = 3)),
      y = y_lab,
      x = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(c("Known", "Novel")),
      method = "wilcox.test",
      label = "p.format"
    )
}

# Generate Comparative Plots
# -------------------------------------------------
plot_list <- list(
  create_violin_plot(combined, "Exon_Total_Length", 
                     "Exon Length Comparison", "Exon Length (bp)"),
  create_violin_plot(combined, "Exon Count", 
                     "Exon Count Comparison", "Number of Exons"),
  create_violin_plot(combined, "CDS Total Length", 
                     "CDS Length Comparison", "CDS Length (bp)"),
  create_violin_plot(combined, "3'UTR Length", 
                     "3'UTR Length Comparison", "3'UTR Length (bp)")
)

# Arrange final composite figure
ggarrange(plotlist = plot_list, ncol = 4)
  
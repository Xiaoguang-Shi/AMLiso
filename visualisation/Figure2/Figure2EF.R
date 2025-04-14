# Gene Structure Visualization Pipeline
# ------------------------------------------
# Required Packages
library(ggplot2)         # Core visualization
library(ggtranscript)    # Transcript structure visualization 
library(dplyr)           # Data manipulation
library(gggenes)         # Protein domain visualization

# Data Preparation
# ------------------------------------------
# Set working directory containing annotation files
setwd("/path/to/your/annotation/files")

# Load exon/CDS coordinates data
exons <- read.table("CAMP.tsv", header = TRUE) 

# Data Processing
# ------------------------------------------
# Separate CDS and exon regions
camp_cds <- exons %>% filter(type == "CDS")
camp_exon <- exons %>% filter(type == "exon")

# Transcript Ordering
# ------------------------------------------
# Add prefix for proper transcript ordering
camp_exon$transcript_id <- paste0("1_", camp_exon$transcript_id)
camp_cds$transcript_id <- paste0("1_", camp_cds$transcript_id)

# Convert to factor for reverse ordering in visualization
camp_exon$transcript_id <- factor(camp_exon$transcript_id, 
                                  levels = rev(unique(camp_exon$transcript_id)))
camp_cds$transcript_id <- factor(camp_cds$transcript_id,
                                 levels = rev(unique(camp_cds$transcript_id)))

# Transcript Structure Visualization
# ------------------------------------------
p1 <- camp_exon %>%
  ggplot(aes(
    xstart = start,  # Genomic start position
    xend = end,      # Genomic end position
    y = transcript_id, # Transcript identifier
    fill = classification # Novel/Known classification
  )) +
  # White background ranges for exons
  geom_range(
    fill = "white",
    height = 0.25    # Control visual thickness
  ) +
  # CDS regions overlay
  geom_range(
    data = camp_cds  # CDS-specific data
  ) +
  # Intron connections with strand awareness
  geom_intron(
    data = to_intron(camp_exon, "transcript_id"),
    aes(strand = strand)
  ) +
  # Color scheme for novel/known isoforms
  scale_fill_manual(values = c("known" = "blue", "novel" = "red")) +
  # Clean theme setup
  theme_classic() +
  ylab(NULL) +  # Remove y-axis label
  theme(
    legend.position = "top",
    axis.line.x = element_blank(),  # Remove axis lines
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),   # Remove tick marks
    axis.text.x = element_blank()   # Hide x-axis text
  )

# Protein Domain Visualization
# ------------------------------------------
# Load protein domain annotations
inter <- read.table("camp.interpro.tsv", sep = "\t", 
                    header = TRUE, quote = "")

# Domain Color Scheme
domain_colors <- c(
  "#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#984EA3",
  "#40E0D0", "#FFC0CB", "#00BFFF", "#FFDEAD", "#EE82EE", 
  "#00FFFF", "grey"
)

# Filter relevant domain databases
df <- inter %>% filter(
  dataset %in% c("MobiDBLite", "ProSitePatterns", "Pfam")
)

# Domain Plot Construction
# ------------------------------------------
df$id <- factor(df$id, levels = rev(unique(df$id)))

p2 <- df %>%
  ggplot(aes(
    xmin = start,    # Domain start position
    xmax = end,      # Domain end position
    y = id,          # Protein identifier
    fill = domain1   # Domain type
  )) +
  # Custom color mapping
  scale_fill_manual(values = domain_colors) +
  # Domain representation as arrows
  geom_gene_arrow(
    arrowhead_height = unit(0, "mm"),  # Disable arrowheads
    arrowhead_width = unit(0, "mm"),
    arrow_body_height = unit(3, "mm")  # Domain bar height
  ) +
  # Protein length indicators
  geom_segment(
    aes(y = id, yend = id, x = 0, xend = length),
    linetype = 1
  ) +
  # Length markers
  geom_point(aes(x = length), shape = "|", show.legend = FALSE) +
  geom_point(aes(x = 0), shape = "|", show.legend = FALSE) +
  # Theme customization
  theme_genes() +
  ylab(NULL) +
  xlab("Protein Length (aa)") +
  labs(fill = "Domain Type") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Composite Figure Assembly
# ------------------------------------------
library(ggpubr)
final_plot <- ggarrange(
  p1, p2,
  ncol = 2,          # Side-by-side layout
  widths = c(1, 1.2) # Adjust panel ratio
)

# Export Options
# ------------------------------------------
# ggsave("gene_structure.pdf", final_plot, width = 12, height = 8)
# ggsave("gene_structure.png", final_plot, dpi = 300, width = 12, height = 8)
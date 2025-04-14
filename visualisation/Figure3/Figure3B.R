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

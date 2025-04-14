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




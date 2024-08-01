# Set the working directory to the specified path
setwd("C:/Users/achir/Desktop/Rayca Bio/New Folder")

# Load necessary libraries
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
library(dplyr)

# Read the clinical data CSV file into a dataframe
data1 = read.csv("brca_tcga_clinical_data.csv")
View(data1)  # Display the data in a spreadsheet-like format
dim(data1)   # Print the dimensions of the dataframe

# Select relevant columns for survival analysis
survival_data = data1 %>%
  select('Patient.ID', 'Overall.Survival..Months.', 'Overall.Survival.Status',
         'Race.Category', 'Sex', 'Subtype', 'Cancer.Type.Detailed')

View(survival_data)  # Display the selected columns

# Create a new column 'event' indicating if the patient is deceased (1) or alive (0)
survival_data$event = ifelse(survival_data$Overall.Survival.Status == "1:DECEASED", 1, 0)
head(survival_data)  # Display the first few rows of the modified dataframe

# Convert overall survival from months to years
survival_data["year"] = survival_data["Overall.Survival..Months."] * 0.083333333
survival_data  # Display the dataframe with the new 'year' column

# Save the modified dataframe to a new CSV file
write.csv(survival_data, file = "Clinical data with events.csv")

# Fit a survival curve using the Cox proportional hazards model with 'Subtype' as the strata
input_data = survfit(Surv(as.numeric(year), as.numeric(event)) ~ Subtype, data = survival_data)
input_data  # Print the fitted survival curve object

# Load the survminer library for advanced survival analysis visualization
library(survminer)

# Plot the survival curves using ggsurvplot
ggsurvplot(
  input_data,
  data = survival_data,          # Use the survival data for the plot
  surv.median.line = 'hv',       # Add horizontal and vertical lines at the median survival
  pval = TRUE,                   # Display the p-value of the log-rank test
  risk.table = TRUE,             # Include a risk table below the survival plot
  conf.int = TRUE,               # Display the confidence intervals for the survival curves
  tables.theme = theme_cleantable(), # Use a clean table theme for the risk table
  break.time.by = 2,             # Break the time axis in increments of 2 years
  ggtheme = theme_bw(),          # Use the black-and-white theme for the plot
  legend.title = "Subtype",      # Title for the legend
  legend.labs = c("Basal", "HER2+", "Luminal A", "Luminal B", "BRCA_Normal"), # Labels for the legend
  xlab = "Time (Years)",         # Label for the x-axis
  ylab = "Survival Probability", # Label for the y-axis
  title = "Survival Curves by Breast Cancer Subtype" # Title of the plot
)



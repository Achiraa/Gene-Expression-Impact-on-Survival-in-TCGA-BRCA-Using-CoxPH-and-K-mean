# Set the working directory to the specified path
setwd("C:/Users/achir/Desktop/Rayca Bio")

# Load necessary libraries
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# Read the Cox data CSV file into a dataframe
cox_data = read.csv("cox file.csv")
View(cox_data)  
head(cox_data)  

# Convert the 'year' column to an integer (floor values)
cox_data$year <- floor(cox_data$year)

# Create a survival object
sv = with(cox_data, Surv(year, event))
sv  

# Fit survival curves without grouping and by subtype
sv_fit <- survfit(Surv(year, event) ~ 1, data = cox_data)
sv_fit1 <- survfit(Surv(year, event) ~ Subtype, data = cox_data)

# Plot the survival curves
autoplot(sv_fit)
autoplot(sv_fit1)

# Create a formula with all 500 genes for the Cox model
gene_columns <- colnames(cox_data)[6:ncol(cox_data)]  
gene_columns
formula <- as.formula(paste("Surv(year, event) ~", paste(gene_columns, collapse = " + ")))
formula  

# Fit the Cox proportional hazards model
cox_model <- coxph(formula, data = cox_data, control = coxph.control(iter.max = 1000))
summary(cox_model)  

# Get the summary of the model
cox_summary <- summary(cox_model)

# Extract the p-values of the genes
p_values <- cox_summary$coefficients[, "Pr(>|z|)"]

# Select the top 100 genes based on the smallest p-values
top_genes <- names(sort(p_values))[1:100]
top_genes  

# Extract data for the top 100 genes
top_gene_data <- cox_data[, top_genes]
write.csv(top_gene_data, file = "Top gene 100.csv")  

# Create a new formula with the top 100 genes
top_genes_formula <- as.formula(paste("Surv(year, event) ~", paste(top_genes, collapse = " + ")))

# Fit the Cox model with the top 100 genes
top_genes_cox_model <- coxph(top_genes_formula, data = cox_data)

# Get the summary of the new model
summary(top_genes_cox_model)

# Plot the survival curve for the Cox model with the top 100 genes
coxplot2 = survfit(top_genes_cox_model)
autoplot(coxplot2)


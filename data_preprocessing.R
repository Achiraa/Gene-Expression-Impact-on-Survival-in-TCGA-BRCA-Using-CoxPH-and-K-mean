# Set the working directory to the specified path
setwd("C:/Users/achir/Desktop/Rayca Bio/New Folder")

# Load the dplyr library for data manipulation
library(dplyr)

# Read the RSEM data CSV file into a dataframe
data_preprocessed = read.csv("data_rsem.csv")
dim(data_preprocessed)  # Print the dimensions of the dataframe
View(data_preprocessed) 

# Create unique row names using the second column of the dataframe
unique_row_names <- make.unique(as.character(data_preprocessed[, 2]))
rownames(data_preprocessed) = unique_row_names
View(data_preprocessed) 

# Remove the second column as it's now used as row names
countfile = data_preprocessed[, -2]
View(countfile)  

# Calculate the total counts for each row and add it as a new column
countfile$total = rowSums(countfile[, -1])
head(countfile)  # Display the first few rows of the dataframe

# Select the top 500 rows based on the highest total counts
top_500_data <- countfile %>% arrange(desc(total)) %>% slice(1:500)
View(top_500_data) # Display the top 500 rows
dim(top_500_data)  # Print the dimensions of the top 500 rows dataframe

# Save the top 500 data to a new CSV file
write.csv(top_500_data, file = "count_data_filtered500.csv")









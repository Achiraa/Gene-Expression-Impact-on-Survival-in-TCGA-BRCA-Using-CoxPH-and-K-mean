# Importing necessary libraries
import random
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

# Load the list of genes
gene_list_df = pd.read_csv('list of gene 100.csv')
print("CSV file loaded successfully.")
print(gene_list_df.head())  

# Convert the gene list dataframe to a list of gene names
gene_list = gene_list_df.iloc[:, 0].tolist()

# Randomly pick a gene from the gene list
selected_gene = random.choice(gene_list)
print(f'Selected Gene: {selected_gene}')

# Load the dataset including cluster information and selected gene expression
data_set = pd.read_csv('cluster file edit.csv')
print(data_set.dtypes)

# Define the number of optimal clusters
optimal_clusters = 4

# Perform Cox Proportional Hazards model analysis for each cluster
for cluster in range(optimal_clusters):
    cluster_data = data_set[data_set['Cluster'] == cluster] 
    
    # Prepare data for the CoxPH model
    cox_data = cluster_data[['year', 'event', selected_gene]]
    
    # Initialize the CoxPHFitter
    cph = CoxPHFitter()
    
    # Fit the CoxPH model
    cph.fit(cox_data, duration_col='year', event_col='event')
    
    # Print the summary of the CoxPH model for the current cluster
    print(f'\nCluster {cluster} CoxPH Model Summary:')
    print(cph.summary[['coef', 'exp(coef)', 'p']])

# Perform Kaplan-Meier Survival Analysis for each cluster
for cluster in range(optimal_clusters):
    cluster_data = data_set[data_set['Cluster'] == cluster]  
    
    # Split the data into high and low expression groups based on the median expression of the selected gene
    median_expression = cluster_data[selected_gene].median()
    high_expr = cluster_data[cluster_data[selected_gene] > median_expression]
    low_expr = cluster_data[cluster_data[selected_gene] <= median_expression]
    
    # Initialize the Kaplan-Meier fitter for high and low expression groups
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    # Fit the Kaplan-Meier model for high expression group
    kmf_high.fit(durations=high_expr['year'], event_observed=high_expr['event'], label='High Expression')
    
    # Fit the Kaplan-Meier model for low expression group
    kmf_low.fit(durations=low_expr['year'], event_observed=low_expr['event'], label='Low Expression')
    
    # Plot the survival curves for high and low expression groups
    plt.figure()
    kmf_high.plot_survival_function()
    kmf_low.plot_survival_function()
    plt.title(f'Kaplan-Meier Survival Curve for {selected_gene} in Cluster {cluster}')
    plt.xlabel('Time (years)')
    plt.ylabel('Survival Probability')
    plt.legend()
    plt.show()
    
    # Perform log-rank test to compare survival distributions between high and low expression groups
    logrank_test_result = logrank_test(high_expr['year'], low_expr['year'], event_observed_A=high_expr['event'], event_observed_B=low_expr['event'])
    
    print(f'\nCluster {cluster} Log-Rank Test p-value:', logrank_test_result.p_value)

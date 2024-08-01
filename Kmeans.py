# Importing libraries
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import random 
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE

# Set a random seed for reproducibility
random.seed(123)

# Read the CSV file into a dataframe
data_set = pd.read_csv('cluster file.csv')

# Extract gene expression data from the dataset 
gene_expression_data = data_set.iloc[:, 5:].values
print(gene_expression_data)  

# Standardize the gene expression data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(gene_expression_data)

# Initializing the list to store the WCSS for different cluster counts
wcss_list = [] 

# Calculate WCSS for different numbers of clusters to use the Elbow Method
for i in range(1, 11):
    kmeans = KMeans(n_clusters=i, init='k-means++', random_state=16)
    kmeans.fit(scaled_data)
    wcss_list.append(kmeans.inertia_)

# Plot the Elbow Method graph
plt.plot(range(1, 11), wcss_list)
plt.title('The Elbow Method')
plt.xlabel('Number of clusters (k)')
plt.ylabel('WCSS')
plt.show() 

# Training the K-means model on the dataset with the optimal number of clusters
optimal_clusters = 4  # Set the optimal number of clusters based on the Elbow Method
kmeans = KMeans(n_clusters=optimal_clusters, init='k-means++', random_state=16) 
y_predict = kmeans.fit_predict(scaled_data)

# Assign each patient to a cluster
data_set['Cluster'] = y_predict
labels = data_set['Cluster']

# Calculate the silhouette score to evaluate the clustering
silhouette_avg = silhouette_score(scaled_data, labels)
print(f'Silhouette Score: {silhouette_avg}')


print(data_set[['Patient samples', 'Cluster']].head())
data_set[['Patient samples', 'Cluster']].to_csv('patient_clusters.csv', index=False)

# Visualize the clusters using PCA (Principal Component Analysis)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)

# Add PCA results to the dataframe
data_set['PCA1'] = pca_result[:, 0]
data_set['PCA2'] = pca_result[:, 1]

# Plot the PCA results with cluster labels
plt.figure(figsize=(10, 8))
sns.scatterplot(x='PCA1', y='PCA2', hue='Cluster', palette='viridis', data=data_set)
plt.title('PCA of Gene Expression Data')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.legend()
plt.show()

# t-SNE Visualization
tsne = TSNE(n_components=3, random_state=42)
tsne_result = tsne.fit_transform(scaled_data)

# Add t-SNE results to the dataframe
data_set['TSNE1'] = tsne_result[:, 0]
data_set['TSNE2'] = tsne_result[:, 1]

# Plot the t-SNE results with cluster labels
plt.figure(figsize=(10, 8))
sns.scatterplot(x='TSNE1', y='TSNE2', hue='Cluster', palette='viridis', data=data_set)
plt.title('t-SNE of Gene Expression Data')
plt.xlabel('TSNE1')
plt.ylabel('TSNE2')
plt.legend()
plt.show()



# # Function to calculate silhouette score for different random states
# def get_best_random_state(scaled_data, n_clusters, random_states):
#     best_score = -1
#     best_state = None
#     best_labels = None
#     for state in random_states:
#         kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=state)
#         y_predict = kmeans.fit_predict(scaled_data)
#         score = silhouette_score(scaled_data, y_predict)
#         if score > best_score:
#             best_score = score
#             best_state = state
#             best_labels = y_predict
#     return best_state, best_score, best_labels

# Check silhouette scores for a range of random states
# random_states = range(1, 101)  # You can adjust this range as needed
# optimal_clusters = 4
# best_random_state, best_silhouette_score, best_labels = get_best_random_state(scaled_data, optimal_clusters, random_states)

# print(f'Best Random State: {best_random_state}')
# print(f'Best Silhouette Score: {best_silhouette_score}')
#%%
from matplotlib import legend
import pandas as pd
import numpy as np
from scipy.stats import chi2
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import statsmodels.api as sm
import scanpy as sc


#%%
############################################################################
# Load data
data_df = pd.read_csv("./input/tcga/TCGA-PAAD.csv",  index_col=0)
clinical = pd.read_csv("./input/tcga/TCGA-PAAD_clinical.csv",index_col=0)
data_df.columns = data_df.columns.str.replace('.', '-')
ppi_df = pd.read_excel('./input/ppi.xlsx')
genes = data_df.index

# data pre-process
control_samples = clinical.query("sample_type=='Solid Tissue Normal'").index.to_list()
disease_samples = clinical.query("sample_type=='Primary Tumor'").index.to_list()
n_control = len(control_samples)
n_disease = len(disease_samples)
# Reorder columns to have control samples first, followed by disease samples
data_df = data_df[control_samples + disease_samples]
data = data_df.groupby(level=0).mean()
co = data[control_samples]
di = data[disease_samples]


#%%
############################################################################
plt.figure(figsize=(10,6))
sns.boxplot(data=data)
plt.title('Boxplot of Counts Across Samples')
plt.ylabel('Raw Counts') 
plt.xticks([])
plt.show()


#%%
############################################################################
melted_df = data.melt(var_name='Sample', value_name='Expression')
plt.figure(figsize=(12, 6))
sns.boxplot(x='Sample', y='Expression', data=melted_df)
plt.xticks(rotation=90)
plt.show()


#%%
############################################################################
plt.figure(figsize=(10,6))
for sample in data.columns:
    data[sample].plot(kind='kde')
plt.title('Density Plot of Raw Counts Across Samples')
plt.xlabel('Raw Counts')
plt.ylabel('Density')
plt.show()


#%%
############################################################################
# Clear separation between conditions (e.g., health vs. disease) 
# with minimal technical-driven clusters indicates good normalization.
############################################################################
pca = PCA(n_components=2)
pca_result = pca.fit_transform(data.T)

# Plot PCA result
plt.figure(figsize=(8, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1])
plt.title('PCA Plot of Samples')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.show()




#%%
############################################################################
# DEG
############################################################################
gene = "A1BG"
group1 = co.loc[gene]
group2 = di.loc[gene]
# Perform t-test
tstat, pval, dfreedom = sm.stats.ttest_ind(group1, group2)    
print(f"Differential expression for {gene}: tstat {tstat} pval {pval}")




#%%
############################################################################
# Uniform color distribution and absence of distinct blocks corresponding to 
# technical batches suggest effective normalization.
############################################################################
top_genes = data.var(axis=1).sort_values(ascending=False).head(100).index
top_df = data.loc[top_genes]
plt.figure(figsize=(12, 10))
sns.heatmap(top_df, cmap='viridis', xticklabels=True, yticklabels=False)
plt.title('Heatmap of Top 500 Variable Genes (Normalized)')
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.show()






#%%
############################################################################
# MA plots should show balanced M values around zero across all A levels, 
# indicating no systematic bias.
############################################################################
co_mean = co.mean(axis=1)
di_mean = di.mean(axis=1)
M = np.log2(di_mean / co_mean)
A = 0.5 * (np.log2(di_mean) + np.log2(co_mean))

plt.figure(figsize=(10, 6))
sns.scatterplot(x=A, y=M, alpha=0.5)
plt.axhline(0, color='red', linestyle='--')
plt.title('MA Plot')
plt.xlabel('A (Mean Average)')
plt.ylabel('M (Log Ratio)')
plt.show()




#%%
############################################################################
# similar means and standard deviations across samples indicate effective normalization.
############################################################################
summary_stats = data.describe().T

plt.figure(figsize=(10, 6))
sns.boxplot(data=summary_stats[['mean', '50%', 'std']])
plt.title('Summary Statistics After Normalization')
plt.show()




#%%
############################################################################

adata = sc.AnnData(X=data.T.values)
adata.var_names = data.index
adata.obs_names = data.columns
# Normalize counts to counts per million (CPM) or similar scale
#sc.pp.normalize_total(adata, target_sum=1e6)
# Log-transform the data
#sc.pp.log1p(adata)

#%%
adata.obs['Condition'] = clinical.sample_type

#%%
# Perform PCA
sc.pp.pca(adata)
sc.pl.pca(adata, color="Condition") 

#%%
# Compute neighbors for UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color='Condition')  # Use metadata for coloring


#%%
# Compute dendrogram
sc.tl.dendrogram(adata, groupby='Condition')
sc.pl.heatmap(adata, var_names=adata.var_names[:50], groupby='Condition', cmap='viridis')


#%%
# Set your comparison group
# Perform differential expression
sc.tl.rank_genes_groups(adata, groupby='Condition', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=10)
top_genes = adata.uns['rank_genes_groups']['names'][:10]
print(top_genes)

#%%
import gseapy as gp
# Example of running GSEA
gsea_results = gp.prerank(rnk="path_to_ranked_gene_list.rnk", gene_sets="KEGG_2019_Human")


#%%
# Apply batch correction
sc.pp.combat(adata, key='Batch')  

#%%
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=100)
sc.pl.heatmap(adata, var_names=adata.var.highly_variable, groupby='Condition', cmap='viridis')


#%%
# Compute correlation matrix
correlation_matrix = pd.DataFrame(adata.X).T.corr()
sns.heatmap(correlation_matrix, cmap='coolwarm', square=True)
plt.show()


#%%
adata_subset = adata[:, "A1BG"]
sc.pl.pca(adata_subset)

#%%
adata_genes = adata
sc.pp.scale(adata_genes, max_value=10)
sc.pp.neighbors(adata_genes, n_neighbors=10, use_rep='X')  # Adjust `n_neighbors` if needed
sc.tl.leiden(adata_genes, resolution=0.2)  # Adjust resolution for more or fewer clusters
# Compute UMAP for visualization
sc.tl.umap(adata_genes)
# Plot UMAP with Leiden clusters
sc.pl.umap(adata_genes, color='leiden')








#%%
###############################################################

data = pd.read_parquet("./../orthoclust/input/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.p")
clinical = pd.read_parquet("./../orthoclust/input/depmap/cohorts.p")
sample_list = clinical[["slug","sample_list"]].explode("sample_list").set_index("sample_list")
data = data.loc[sample_list.index.values]
#%%
# Get common samples between data and sample_list
common_samples = data.columns.intersection(sample_list.index)
sample_list = sample_list.loc[sample_list.index.isin(common_samples)]
data = data[common_samples]

#%%
# Create AnnData object
adata = sc.AnnData(X=data.T.values)
adata.obs_names = data.columns  # Set cell lines as observations
adata.var_names = data.index  # Set genes as variables
adata.obs = sample_list
# %%

# Scale the data (center and scale genes)
sc.pp.scale(adata, max_value=10)

#%%
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.leiden(adata, resolution=1.0)
#%%
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')
#%%

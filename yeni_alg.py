#%%
import pandas as pd
import numpy as np
import functions as f


#%%
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
data_df = data_df.groupby(level=0).mean()


#%%
NCData, NDData = f.binarize_expression(data_df, n_control)
NC, ND = f.count_interactions(NCData, NDData, genes, ppi_df)
df = f.calculate_entropy_for_each_state_columns(NC, ND, n_control, n_disease, normalize=True)
df.index = list(zip(ppi_df.Gene1,ppi_df.Gene2))

#%%
subset_df = df[['State_11_c', 'State_11_d', 'State_11_p', 'State_11_adjusted_p',
                'State_11_h', 'State_11_significant',
                'State_-1-1_c', 'State_-1-1_d', 'State_-1-1_p',
                'State_-1-1_adjusted_p', 'State_-1-1_h', 'State_-1-1_significant']]

subset_df.to_csv('TCGA_subset_result.csv')
df.to_csv('TCGA_result.csv')


#%%
for uniq_state in ["State_11_significant", "State_-1-1_significant"]:
    significant_state = subset_df[subset_df['State_11_significant'] == True].index.tolist() 
    unique_proteins = list(set([protein for pair in significant_state for protein in pair]))
    np.savetxt(f'TCGA_{uniq_state}_genes.txt', unique_proteins, fmt="%s")


#%%

# Read the TCGA and NCBI results
tcga_df = pd.read_csv('result/tcga/TCGA_result.csv', index_col=0)
ncbi_df = pd.read_csv('result/ncbi/NCBI_result.csv', index_col=0)

# Get significant pairs for State_11 from both datasets
tcga_sig = tcga_df[tcga_df['State_-1-1_significant'] == True].index.tolist()
ncbi_sig = ncbi_df[ncbi_df['State_-1-1_significant'] == True].index.tolist()

# Convert string tuples to actual tuples
tcga_sig = [eval(pair) for pair in tcga_sig]
ncbi_sig = [eval(pair) for pair in ncbi_sig]

# Get unique genes from each dataset
tcga_genes = set([gene for pair in tcga_sig for gene in pair])
ncbi_genes = set([gene for pair in ncbi_sig for gene in pair])

# Find common genes
common_genes = tcga_genes & ncbi_genes

print(f"Number of unique genes in TCGA: {len(tcga_genes)}")
print(f"Number of unique genes in NCBI: {len(ncbi_genes)}")
print(f"Number of common genes: {len(common_genes)}")


# Find common pairs
common_pairs = set(tcga_sig) & set(ncbi_sig)

print(f"Number of significant pairs in TCGA: {len(tcga_sig)}")
print(f"Number of significant pairs in NCBI: {len(ncbi_sig)}")
print(f"Number of common significant pairs: {len(common_pairs)}")

# Save common pairs to txt file
common_pairs_list = list(common_pairs)
np.savetxt('State_-1-1_common_significant_pairs.txt', common_pairs_list, fmt='%s')

# Save common genes to txt file 
common_genes_list = list(common_genes)
np.savetxt('State_-1-1_common_genes.txt', common_genes_list, fmt='%s')




# %%

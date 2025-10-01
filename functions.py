import pandas as pd
import numpy as np
from scipy.stats import chi2
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns

column_names = ['11', '00', '10', '01', '-10', '0-1', '-11', '1-1', '-1-1']


def binarize_expression(data, n_control):

    # Convert pandas DataFrame to numpy array if needed
    if isinstance(data, pd.DataFrame):
        data = data.values

    print(f"data shape: {data.shape}")
   
    # Binarize expression levels
    avg_exp = np.mean(data, axis=1)
    std_dev = np.std(data, axis=1)
    NData = np.zeros(data.shape)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            upper_threshold = avg_exp[j] + 1.96 * std_dev[j]
            lower_threshold = avg_exp[j] - 1.96 * std_dev[j]
            if data[i, j] >= upper_threshold:
                NData[i, j] = 1
            elif data[i, j] <= lower_threshold:
                NData[i, j] = -1
            else:
                NData[i, j] = 0
    
    # Split the data
    NCData = NData[:, :n_control]
    NDData = NData[:, n_control:]
    
    return NCData, NDData



def count_interactions(NCData, NDData, genes, ppi_df):
    interaction_counts_control = np.zeros((ppi_df.shape[0], 9))  # Control interaction counts
    interaction_counts_disease = np.zeros((ppi_df.shape[0], 9))  # Disease interaction counts

    states = {
        (1, 1): 0,
        (0, 0): 1,
        (1, 0): 2,
        (0, 1): 3,
        (-1, 0): 4,
        (0, -1): 5,
        (-1, 1): 6,
        (1, -1): 7,
        (-1, -1): 8
    }
    
    for z, row in ppi_df.iterrows():
        gene1, gene2 = row['Gene1'], row['Gene2']

        if gene1 in genes and gene2 in genes:
            idx1 = np.where(genes == gene1)[0][0]
            idx2 = np.where(genes == gene2)[0][0]

            for t in range(NCData.shape[1]):  # Control samples
                state = (NCData[idx1, t], NCData[idx2, t])
                if state in states:
                    interaction_counts_control[z, states[state]] += 1

            for t in range(NDData.shape[1]):  # Disease samples
                state = (NDData[idx1, t], NDData[idx2, t])
                if state in states:
                    interaction_counts_disease[z, states[state]] += 1
                
    return interaction_counts_control, interaction_counts_disease

def calculate_entropy_for_each_state_columns(NC, ND, n_control, n_disease, normalize=True):

    # Normalize counts if specified
    if normalize:
        NC = NC / n_control
        ND = ND / n_disease

    column_names = ['11', '00', '10', '01', '-10', '0-1', '-11', '1-1', '-1-1']
    result_dict = {}

    # Loop over each state (column index)
    for i in range(NC.shape[1]):
        c = NC[:, i]
        d = ND[:, i]

        # Calculate the probability P = D / (D + C)
        p = d / (d + c)

        # Small constant for smoothing
        epsilon = 1e-5
        # Adjusted probability calculation with smoothing
        adjusted_p = (d + epsilon) / (d + c + 2 * epsilon)

        # Calculate entropy with the adjusted probabilities
        h = np.array([entropy([p_val, 1 - p_val], base=2) for p_val in adjusted_p])

        # Store each result in the dictionary with formatted column names
        result_dict[f'State_{column_names[i]}_c'] = c
        result_dict[f'State_{column_names[i]}_d'] = d
        result_dict[f'State_{column_names[i]}_p'] = p
        result_dict[f'State_{column_names[i]}_adjusted_p'] = adjusted_p
        result_dict[f'State_{column_names[i]}_h'] = h

        # Calculate percentiles for significance thresholds
        p_95th = np.percentile(adjusted_p, 99)  # 95th percentile for P
        h_5th = np.percentile(h, 1)    # 5th percentile for H

        # Determine significance for each interaction
        significant = (adjusted_p >= p_95th) & (h <= h_5th)

        # Add significance to the result dictionary
        result_dict[f'State_{column_names[i]}_significant'] = significant

    # Convert dictionary to DataFrame
    result_df = pd.DataFrame(result_dict)
    return result_df
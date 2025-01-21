import numpy as np

from .calculations import apply_hypothesis_test

from .z_score import calculate_z_score

from .ranking import calculate_ranks_for_sample


def NMD_test_statistic(group):
    """
    Calculate test statistic for NMD.
    This function calculates a unique test statistic for each sample in the group.
    """
    results = []
    for _, row in group.iterrows():
        # Calculate the test statistic for this specific sample
        ratio = (row['Cyclo_TPM'] + 1) / (row['Noncyclo_TPM'] + 1)
        test_statistic = np.log2(ratio) * np.log2(row['cyclo_count'] + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def Noncyclo_Expression_Outlier_LOE(group):
    """Calculate test statistic for Noncyclo Expression Outlier - Loss of Expression (LOE)."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        min_noncyclo_TPM = other_samples['Noncyclo_TPM'].min()
        median_noncyclo_TPM = other_samples['Noncyclo_TPM'].median()
        
        # Calculate the test statistic
        ratio = (min_noncyclo_TPM + 1) / (row['Noncyclo_TPM'] + 1)
        test_statistic = np.log2(ratio) * np.log2(median_noncyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def Noncyclo_Expression_Outlier_GOE(group):
    """Calculate test statistic for Noncyclo Expression Outlier - Gain of Expression (GOE)."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        max_noncyclo_TPM = other_samples['Noncyclo_TPM'].max()
        median_noncyclo_TPM = other_samples['Noncyclo_TPM'].median()
        
        # Calculate the test statistic
        ratio = (row['Noncyclo_TPM'] + 1) / (max_noncyclo_TPM + 1)
        test_statistic = np.log2(ratio) * np.log2(median_noncyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def Cyclo_Expression_Outlier_LOE(group):
    """Calculate test statistic for Cyclo Expression Outlier - Loss of Expression (LOE)."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        min_cyclo_TPM = other_samples['Cyclo_TPM'].min()
        median_cyclo_TPM = other_samples['Cyclo_TPM'].median()
        
        # Calculate the test statistic
        ratio = (min_cyclo_TPM + 1) / (row['Cyclo_TPM'] + 1)
        test_statistic = np.log2(ratio) * np.log2(median_cyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def Cyclo_Expression_Outlier_GOE(group):
    """Calculate test statistic for Cyclo Expression Outlier - Gain of Expression (GOE)."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        max_cyclo_TPM = other_samples['Cyclo_TPM'].max()
        median_cyclo_TPM = other_samples['Cyclo_TPM'].median()
        
        # Calculate the test statistic
        ratio = (row['Cyclo_TPM'] + 1) / (max_cyclo_TPM + 1)
        test_statistic = np.log2(ratio) * np.log2(median_cyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def NMD_rare_steady_state_transcript(group):
    """
    Calculate test statistic for NMD Rare Steady State Transcript.
    This function calculates a unique test statistic for each sample in the group.
    """
    results = []
    for _, row in group.iterrows():
        # Calculate the test statistic for this specific sample
        cyclo_term = row['Total_bin_cyclo_count_Bin1_le'] / (row['cyclo_count'] + 1)
        noncyclo_term = row['Total_bin_noncyclo_count_Bin1_le'] / (row['noncyclo_count'] + 1)
        test_statistic = (
            (cyclo_term - noncyclo_term) *
            np.log2(row['Total_bin_cyclo_count_Bin1_le'] + 1) *
            np.log2(row['Total_bin_noncyclo_count_Bin2_g'] + 1)
        )
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def process_hypothesis_test(filtered_data, group_col, test_statistic_func):
    """
    Combine hypothesis testing, z-score calculation, and ranking into a single function.
    
    Parameters:
    - filtered_data (pd.DataFrame): The filtered data to process.
    - group_col (str): The column to group by (e.g., 'Isoform_PBid').
    - test_statistic_func (function): The hypothesis test function to apply.
    
    Returns:
    - pd.DataFrame: Ranked data after applying hypothesis test, calculating z-scores, and ranks.
    """
    # Apply hypothesis test
    tested_data = apply_hypothesis_test(filtered_data, group_col, test_statistic_func)
    
    # Calculate z-scores
    z_scored_data = calculate_z_score(tested_data, group_col=group_col, stat_col='test_statistic')
    
    # Calculate ranks
    ranked_data = calculate_ranks_for_sample(z_scored_data, group_col=group_col)
    
    return ranked_data

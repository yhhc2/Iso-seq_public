import numpy as np

percentiles = [99.9, 99.5, 99, 98, 95]

def calculate_ranks_for_sample(df):
    """
    Calculate ranks for each sample at specified percentiles of z-scores.

    Parameters:
    - df: DataFrame containing columns 'Sample', 'Isoform_PBid', 'z_score', and 'test_statistic'.

    Returns:
    - DataFrame with new rank columns added for each percentile.
    """
    # Filter out rows with NaN z-scores
    filtered_df = df[df['z_score'].notna()].copy()

    # Initialize all rank columns globally
    for percentile in percentiles:
        if percentile == 99.5:
            rank_column_name = 'rank_top_99_5_percentile'
        elif percentile == 99.9:
            rank_column_name = 'rank_top_99_9_percentile'
        else:
            rank_column_name = f'rank_top_{percentile}_percentile'
        filtered_df[rank_column_name] = np.nan

    # Process each sample independently
    for sample in filtered_df['Sample'].unique():
        sample_df = filtered_df[filtered_df['Sample'] == sample]

        for percentile in percentiles:
            # Calculate percentile value
            percentile_value = np.percentile(sample_df['z_score'], percentile)

            # Define rank column name
            if percentile == 99.5:
                rank_column_name = 'rank_top_99_5_percentile'
            elif percentile == 99.9:
                rank_column_name = 'rank_top_99_9_percentile'
            else:
                rank_column_name = f'rank_top_{percentile}_percentile'

            # Get top percentile subset and rank
            top_percentile_df = sample_df[sample_df['z_score'] > percentile_value]
            sorted_df = top_percentile_df.sort_values(by='test_statistic', ascending=False).reset_index(drop=True)
            rank_mapping = {isoform: rank + 1 for rank, isoform in enumerate(sorted_df['Isoform_PBid'])}

            # Update the main DataFrame with ranks
            filtered_df.loc[
                (filtered_df['Sample'] == sample) & 
                (filtered_df['Isoform_PBid'].isin(rank_mapping.keys())),
                rank_column_name
            ] = filtered_df['Isoform_PBid'].map(rank_mapping)

    return filtered_df

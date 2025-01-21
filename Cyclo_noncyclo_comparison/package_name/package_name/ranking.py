import numpy as np

percentiles = [99.9, 99.5, 99, 98, 95]

def calculate_ranks_for_sample(df, group_col='Sample', isoform_col='Isoform_PBid', percentiles=[95, 99, 99.5, 99.9]):
    """
    Calculate ranks for each sample at specified percentiles of z-scores.

    Parameters:
    - df: DataFrame containing columns 'Sample', isoform_col, 'z_score', and 'test_statistic'.
    - group_col: Column to group by (default is 'Sample').
    - isoform_col: Column containing isoform identifiers (default is 'Isoform_PBid').
    - percentiles: List of percentiles to calculate ranks for (default is [95, 99, 99.5, 99.9]).

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

    # Process each group independently
    for sample in filtered_df[group_col].unique():
        sample_df = filtered_df[filtered_df[group_col] == sample]

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
            rank_mapping = {isoform: rank + 1 for rank, isoform in enumerate(sorted_df[isoform_col])}

            # Update the main DataFrame with ranks
            filtered_df.loc[
                (filtered_df[group_col] == sample) & 
                (filtered_df[isoform_col].isin(rank_mapping.keys())),
                rank_column_name
            ] = filtered_df[isoform_col].map(rank_mapping)

    return filtered_df

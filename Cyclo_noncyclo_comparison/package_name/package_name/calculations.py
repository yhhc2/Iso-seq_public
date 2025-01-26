def apply_hypothesis_test(df, group_col, test_statistic_func):
    """
    Apply a hypothesis test statistic function to grouped data.
    
    Parameters:
    - df: Input DataFrame.
    - group_col: Column to group by (e.g., Isoform_PBid).
    - test_statistic_func: Function to calculate the test statistic for each group.
    """
    return df.groupby(group_col, group_keys=False).apply(test_statistic_func).reset_index(drop=True)

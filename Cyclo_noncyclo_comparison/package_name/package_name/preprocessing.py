def filter_isoforms(df, count_threshold=10):
    """
    Filter isoforms based on count thresholds.
    Groups by Isoform_PBid and checks if any sample for that isoform meets the threshold.
    """
    isoforms_to_keep = df.groupby('Isoform_PBid').apply(
        lambda group: any(group['cyclo_count'] >= count_threshold) or any(group['noncyclo_count'] >= count_threshold)
    )
    isoforms_to_keep = isoforms_to_keep[isoforms_to_keep].index.tolist()

    # Filter the DataFrame to include only the isoforms to keep
    return df[df['Isoform_PBid'].isin(isoforms_to_keep)]

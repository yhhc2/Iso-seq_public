import pandas as pd
from collections import defaultdict

def parse_read_stats(file_path):
    """
    Parse the read_stats.txt file to extract sample and PB identifiers.

    Parameters:
    - file_path: Path to the read_stats.txt file.

    Returns:
    - A nested dictionary with PB identifiers as keys and sample counts as values.
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    with open(file_path, 'r') as f:
        for line in f:
            # Split the line into components
            read, pb_id = line.strip().split()
            
            # Extract the sample name (before the first "_")
            sample = read.split('_')[0]
            
            # Increment the count for the PB identifier in the sample
            counts[pb_id][sample] += 1
    
    return counts

def create_expression_matrix(file_path, output_file=None):
    """
    Create an expression matrix from a read_stats.txt file.

    Parameters:
    - file_path: Path to the read_stats.txt file.
    - output_file: Path to save the resulting expression matrix as a CSV (optional).

    Returns:
    - A pandas DataFrame representing the expression matrix.
    """
    # Parse the file and aggregate counts
    counts = parse_read_stats(file_path)
    
    # Create a DataFrame from the nested dictionary
    df = pd.DataFrame.from_dict(counts, orient='index').fillna(0).astype(int)
    
    # Sort rows and columns for better readability
    df.sort_index(inplace=True)
    df.sort_index(axis=1, inplace=True)
    
    # Save to a CSV file if specified
    if output_file:
        df.to_csv(output_file)
    
    return df

#!/usr/bin/env python3
"""
Example showing how to use conditional imports for Google Colab compatibility.

This script demonstrates how to import alphagenome modules with Google Colab 
compatibility using the updated colab_utils module.
"""

# Import the conditional imports from colab_utils
from alphagenome.colab_utils import data_table, files, IN_COLAB

# Standard imports
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import pandas as pd

# The data_table and files objects are now safely imported
# They will work in both Colab and local environments
print(f"Running in Google Colab: {IN_COLAB}")

# Set max rows for data display
if IN_COLAB:
    data_table.DataTable.max_rows = 100_000
    data_table.enable_dataframe_formatter()
else:
    pd.set_option('display.max_rows', 100_000)

# Example usage:
def demonstrate_functionality():
    """Demonstrate the functionality works in both environments."""
    
    # Create sample data
    sample_data = pd.DataFrame({
        'chromosome': ['chr1', 'chr2', 'chr3'],
        'position': [1000, 2000, 3000],
        'score': [0.8, 0.6, 0.9]
    })
    
    print("Sample data:")
    print(sample_data)
    
    # File operations (will behave differently in Colab vs local)
    if IN_COLAB:
        # In Colab, this would actually download the file
        print("\nIn Colab environment - file operations are available")
    else:
        # In local environment, this just prints a message
        print("\nIn local environment - file operations are mocked")
        files.download("sample_data.csv")

if __name__ == "__main__":
    demonstrate_functionality()
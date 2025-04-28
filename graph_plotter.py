import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import seaborn as sns

# Function to parse the input file with flexible format
def parse_data(filename):
    with open(filename, 'r') as file:
        content = file.read()
    
    # Split by dataset sections
    sections = content.split("---------------------------------------------")
    sections = [s.strip() for s in sections if s.strip()]
    
    data = []
    for i, section in enumerate(sections):
        entry = {}
        
        # Try to extract filename if present
        filename_match = re.search(r"File: (.+)", section)
        if filename_match:
            path = filename_match.group(1)
            entry['filename'] = path.split('\\')[-1]
        else:
            entry['filename'] = f"graph_{i+1}"
        
        # Extract nodes and edges
        nodes_match = re.search(r"Nodes: (\d+)", section)
        edges_match = re.search(r"Edges: (\d+)", section)
        if nodes_match:
            entry['nodes'] = int(nodes_match.group(1))
        if edges_match:
            entry['edges'] = int(edges_match.group(1))
        
        # Extract H and Alpha
        h_match = re.search(r"H: (\d+)", section)
        alpha_match = re.search(r"Alpha: (\d+\.\d+)", section)
        if h_match:
            entry['h'] = int(h_match.group(1))
        if alpha_match:
            entry['alpha'] = float(alpha_match.group(1))
        
        # Extract CDS size
        cds_size_match = re.search(r"Nodes in CDS: (\d+)", section)
        if cds_size_match:
            entry['cds_size'] = int(cds_size_match.group(1))
        
        # Extract time
        time_match = re.search(r"Time taken: (\d+) ms", section)
        if time_match:
            entry['time_ms'] = int(time_match.group(1))
        
        if entry:
            data.append(entry)
    
    df = pd.DataFrame(data)

    # Ensure all expected columns exist
    required_columns = ['filename', 'nodes', 'edges', 'h', 'alpha', 'cds_size', 'time_ms']
    for col in required_columns:
        if col not in df.columns:
            df[col] = np.nan  # Fill missing columns with NaN

    print(f"Parsed {len(df)} entries with columns: {df.columns.tolist()}")
    return df

# Create plots
def create_plots(df):
    # Calculate edge density
    df['edge_density'] = df['edges'] / (df['nodes'] * (df['nodes'] - 1) / 2)
    
    # Print data summary
    print("Data summary:")
    if 'time_ms' in df.columns:
        print(df[['filename', 'nodes', 'edges', 'time_ms', 'edge_density']].sort_values('nodes'))
    else:
        print(df[['filename', 'nodes', 'edges', 'edge_density']].sort_values('nodes'))
    
    # Set style for all plots
    sns.set(style="whitegrid")
    
    # 1. Execution Time vs Graph Size
    df_valid = df.dropna(subset=['time_ms'])
    
    if not df_valid.empty:
        plt.figure(figsize=(10, 6))
        plt.scatter(df_valid['nodes'], df_valid['time_ms'], alpha=0.7, s=100)
        plt.xlabel('Number of Nodes in Graph', fontsize=12)
        plt.ylabel('Execution Time (ms)', fontsize=12)
        plt.title('Execution Time vs Graph Size', fontsize=14)
        plt.yscale('log')
        plt.xscale('log')
        
        # Add regression line
        if len(df_valid) > 1:
            x = np.log10(df_valid['nodes'])
            y = np.log10(df_valid['time_ms'])
            mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isinf(x) & ~np.isinf(y) & (df_valid['time_ms'] > 0)
            if sum(mask) > 1:
                z = np.polyfit(x[mask], y[mask], 1)
                p = np.poly1d(z)
                x_range = np.linspace(min(x[mask]), max(x[mask]), 100)
                plt.plot(10**x_range, 10**p(x_range), "r--", alpha=0.7)
                plt.annotate(f"Slope: {z[0]:.2f}", 
                             xy=(0.05, 0.95), xycoords='axes fraction',
                             fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        
        # Add dataset labels
        for i, txt in enumerate(df_valid['filename']):
            plt.annotate(txt, (df_valid['nodes'].iloc[i], df_valid['time_ms'].iloc[i]), 
                         fontsize=9, alpha=0.8)
        
        plt.tight_layout()
        plt.savefig('execution_time_vs_graph_size.png', dpi=300, bbox_inches='tight')
    else:
        print("No valid time_ms data available to plot Execution Time vs Graph Size.")
    
    # 2. Histogram for density distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(df['edge_density'].dropna(), bins=10, kde=True)
    plt.xlabel('Edge Density (Ratio of Edges to Possible Edges)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Distribution of Edge Density Across Datasets', fontsize=14)
    
    # Add vertical line for mean
    mean_density = df['edge_density'].mean()
    plt.axvline(x=mean_density, color='r', linestyle='--', alpha=0.7)
    plt.annotate(f"Mean: {mean_density:.4f}", 
                 xy=(mean_density, plt.ylim()[1]*0.9),
                 xytext=(10, 0), textcoords='offset points',
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('edge_density_histogram.png', dpi=300, bbox_inches='tight')
    
    print("Plots created successfully!")

# Main execution
if __name__ == "__main__":
    df = parse_data("output_733.txt")
    
    # Filter out datasets with 0 or missing nodes
    df = df[df['nodes'].fillna(0) > 0]
    
    create_plots(df)

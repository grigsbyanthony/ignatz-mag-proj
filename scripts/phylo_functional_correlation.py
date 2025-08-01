#!/usr/bin/env python3
"""
Correlate phylogenetic relationships with functional differences
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import warnings
warnings.filterwarnings('ignore')

def load_phylogenetic_data(phylo_dir):
    """Load phylogenetic distance matrix"""
    distance_file = phylo_dir / "distance_matrix.csv"
    if distance_file.exists():
        distance_df = pd.read_csv(distance_file, index_col=0)
        return distance_df
    else:
        print(f"Phylogenetic distance matrix not found: {distance_file}")
        return None

def load_functional_data(results_dir):
    """Load all functional comparison data"""
    comp_dir = results_dir / "comparative_analysis"
    
    functional_data = {}
    
    # Load different functional datasets
    datasets = {
        'cog_comparison': 'COG categories',
        'specialized_genes_comparison': 'Specialized functions',
        'interaction_mechanisms_comparison': 'Interaction mechanisms',
        'nutrient_biosynthesis_comparison': 'Nutrient biosynthesis',
        'terpene_genes_comparison': 'Terpene genes'
    }
    
    for filename, description in datasets.items():
        file_path = comp_dir / f"{filename}.csv"
        if file_path.exists():
            df = pd.read_csv(file_path, index_col=0)
            functional_data[description] = df
            print(f"Loaded {description}: {df.shape}")
        else:
            print(f"File not found: {file_path}")
    
    return functional_data

def calculate_functional_distances(functional_df, metric='euclidean'):
    """Calculate functional distance matrix"""
    # Normalize data
    normalized_df = functional_df.div(functional_df.sum(axis=1), axis=0).fillna(0)
    
    # Calculate pairwise distances
    distances = pdist(normalized_df.values, metric=metric)
    distance_matrix = squareform(distances)
    
    # Convert to DataFrame
    distance_df = pd.DataFrame(distance_matrix, 
                              index=normalized_df.index, 
                              columns=normalized_df.index)
    
    return distance_df

def correlate_phylo_functional(phylo_distances, func_distances, mag_order):
    """Correlate phylogenetic and functional distances"""
    # Ensure same ordering
    common_mags = [mag for mag in mag_order if mag in phylo_distances.index and mag in func_distances.index]
    
    if len(common_mags) < 3:
        return None, None, None
    
    # Extract distance vectors (upper triangle)
    phylo_matrix = phylo_distances.loc[common_mags, common_mags].values
    func_matrix = func_distances.loc[common_mags, common_mags].values
    
    # Get upper triangle indices
    triu_indices = np.triu_indices_from(phylo_matrix, k=1)
    phylo_vector = phylo_matrix[triu_indices]
    func_vector = func_matrix[triu_indices]
    
    # Calculate correlations
    pearson_r, pearson_p = pearsonr(phylo_vector, func_vector)
    spearman_r, spearman_p = spearmanr(phylo_vector, func_vector)
    
    return (pearson_r, pearson_p), (spearman_r, spearman_p), (phylo_vector, func_vector)

def create_correlation_plots(correlations, output_dir):
    """Create correlation plots"""
    n_comparisons = len(correlations)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for i, (func_type, data) in enumerate(correlations.items()):
        if i >= len(axes):
            break
            
        pearson_stats, spearman_stats, vectors = data
        phylo_vector, func_vector = vectors
        
        ax = axes[i]
        ax.scatter(phylo_vector, func_vector, alpha=0.6, s=50)
        ax.set_xlabel('Phylogenetic Distance')
        ax.set_ylabel('Functional Distance')
        ax.set_title(f'{func_type}\nPearson r={pearson_stats[0]:.3f} (p={pearson_stats[1]:.3f})')
        
        # Add trend line
        z = np.polyfit(phylo_vector, func_vector, 1)
        p = np.poly1d(z)
        ax.plot(phylo_vector, p(phylo_vector), "r--", alpha=0.8)
    
    # Remove empty subplots
    for i in range(len(correlations), len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    plt.savefig(output_dir / "phylo_functional_correlations.png", dpi=300, bbox_inches='tight')
    plt.close()

def analyze_functional_clusters(phylo_distances, functional_data, output_dir):
    """Analyze functional differences between phylogenetic clusters"""
    # Get phylogenetic clusters (using 20% divergence threshold)
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    
    # Convert distance matrix to condensed form
    phylo_matrix = phylo_distances.values
    condensed_phylo = squareform(phylo_matrix)
    
    # Perform clustering
    phylo_linkage = linkage(condensed_phylo, method='average')
    clusters = fcluster(phylo_linkage, 0.2, criterion='distance')
    
    # Map clusters to MAGs
    mag_clusters = dict(zip(phylo_distances.index, clusters))
    
    print(f"\nPhylogenetic clusters (20% threshold):")
    cluster_members = {}
    for mag, cluster in mag_clusters.items():
        if cluster not in cluster_members:
            cluster_members[cluster] = []
        cluster_members[cluster].append(mag)
    
    for cluster, members in cluster_members.items():
        print(f"  Cluster {cluster}: {members}")
    
    # Analyze functional differences between clusters
    functional_differences = {}
    
    for func_type, func_df in functional_data.items():
        if func_df.empty:
            continue
            
        cluster_means = {}
        for cluster, members in cluster_members.items():
            if len(members) > 1:  # Only analyze clusters with multiple members
                common_members = [m for m in members if m in func_df.columns]
                if len(common_members) > 0:
                    cluster_mean = func_df[common_members].mean(axis=1)
                    cluster_means[f"Cluster_{cluster}"] = cluster_mean
        
        if len(cluster_means) > 1:
            cluster_df = pd.DataFrame(cluster_means)
            functional_differences[func_type] = cluster_df
    
    return mag_clusters, functional_differences

def create_functional_cluster_heatmaps(functional_differences, output_dir):
    """Create heatmaps showing functional differences between clusters"""
    n_plots = len(functional_differences)
    if n_plots == 0:
        return
    
    fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4*n_plots))
    if n_plots == 1:
        axes = [axes]
    
    for i, (func_type, cluster_df) in enumerate(functional_differences.items()):
        # Select top variable functions
        if len(cluster_df) > 20:
            func_variance = cluster_df.var(axis=1).sort_values(ascending=False)
            top_functions = func_variance.head(20).index
            plot_df = cluster_df.loc[top_functions]
        else:
            plot_df = cluster_df
        
        sns.heatmap(plot_df, annot=True, fmt='.1f', cmap='RdYlBu_r', 
                   ax=axes[i], cbar_kws={'label': 'Gene Count'})
        axes[i].set_title(f'{func_type} - Cluster Differences')
        axes[i].set_ylabel('Function')
    
    plt.tight_layout()
    plt.savefig(output_dir / "functional_cluster_differences.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    phylo_dir = results_dir / "phylogenetic_analysis"
    output_dir = phylo_dir
    
    print("Correlating phylogeny with functional differences...")
    
    # Load phylogenetic distances
    print("\n1. Loading phylogenetic distances...")
    phylo_distances = load_phylogenetic_data(phylo_dir)
    if phylo_distances is None:
        return
    
    # Load functional data
    print("\n2. Loading functional comparison data...")
    functional_data = load_functional_data(results_dir)
    
    if not functional_data:
        print("No functional data found!")
        return
    
    # Calculate functional distances and correlations
    print("\n3. Calculating phylogenetic-functional correlations...")
    correlations = {}
    mag_order = list(phylo_distances.index)
    
    for func_type, func_df in functional_data.items():
        if func_df.empty:
            continue
            
        # Calculate functional distances
        func_distances = calculate_functional_distances(func_df.T)  # Transpose so MAGs are rows
        
        # Correlate with phylogenetic distances
        correlation_results = correlate_phylo_functional(phylo_distances, func_distances, mag_order)
        
        if correlation_results[0] is not None:
            correlations[func_type] = correlation_results
            pearson_stats, spearman_stats, _ = correlation_results
            print(f"  {func_type}:")
            print(f"    Pearson r = {pearson_stats[0]:.3f} (p = {pearson_stats[1]:.3f})")
            print(f"    Spearman r = {spearman_stats[0]:.3f} (p = {spearman_stats[1]:.3f})")
    
    # Create correlation plots
    print("\n4. Creating correlation plots...")
    if correlations:
        create_correlation_plots(correlations, output_dir)
    
    # Analyze functional differences between phylogenetic clusters
    print("\n5. Analyzing functional differences between phylogenetic clusters...")
    mag_clusters, functional_differences = analyze_functional_clusters(phylo_distances, functional_data, output_dir)
    
    # Create cluster comparison heatmaps
    print("\n6. Creating functional cluster comparison plots...")
    if functional_differences:
        create_functional_cluster_heatmaps(functional_differences, output_dir)
    
    # Summary statistics
    print("\n=== PHYLOGENY-FUNCTION CORRELATION SUMMARY ===")
    
    if correlations:
        significant_correlations = []
        for func_type, (pearson_stats, spearman_stats, _) in correlations.items():
            if pearson_stats[1] < 0.05:  # Significant correlation
                significant_correlations.append((func_type, pearson_stats[0], pearson_stats[1]))
        
        if significant_correlations:
            print(f"\nSignificant correlations (p < 0.05):")
            for func_type, r, p in significant_correlations:
                print(f"  {func_type}: r = {r:.3f}, p = {p:.3f}")
        else:
            print(f"\nNo significant correlations found (p < 0.05)")
        
        # Best correlations regardless of significance
        print(f"\nStrongest correlations:")
        correlation_strengths = [(func_type, abs(stats[0][0])) for func_type, stats in correlations.items()]
        correlation_strengths.sort(key=lambda x: x[1], reverse=True)
        
        for func_type, abs_r in correlation_strengths[:3]:
            pearson_stats, _, _ = correlations[func_type]
            print(f"  {func_type}: |r| = {abs_r:.3f}, p = {pearson_stats[1]:.3f}")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - phylo_functional_correlations.png")
    print("  - functional_cluster_differences.png")

if __name__ == "__main__":
    main()
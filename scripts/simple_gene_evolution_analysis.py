#!/usr/bin/env python3
"""
Simplified Gene Evolution Analysis for Ignatzschineria MAGs
Focus on gain/loss patterns without complex tree reconstruction
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

def load_pangenome_data(pangenome_dir):
    """Load pan-genome classification data"""
    categories_file = pangenome_dir / "pangenome_categories.csv"
    pa_file = pangenome_dir / "presence_absence_matrix.csv"
    
    if not categories_file.exists() or not pa_file.exists():
        print("Required pan-genome files not found!")
        return None, None
    
    categories_df = pd.read_csv(categories_file)
    pa_matrix = pd.read_csv(pa_file, index_col=0)
    
    return categories_df, pa_matrix

def load_phylogenetic_relationships(phylo_dir):
    """Load phylogenetic distance matrix"""
    distance_file = phylo_dir / "distance_matrix.csv"
    if not distance_file.exists():
        print("Phylogenetic distance matrix not found!")
        return None
    
    distance_df = pd.read_csv(distance_file, index_col=0)
    return distance_df

def load_functional_annotations(results_dir):
    """Load functional annotations"""
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    annotations = {}
    
    for file_path in annotation_files:
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#', 
                           names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                                 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                                 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                                 'BiGG_Reaction', 'PFAMs'])
            
            for _, row in df.iterrows():
                og_groups = str(row['eggNOG_OGs'])
                
                if og_groups != 'nan' and og_groups != '-':
                    primary_og = og_groups.split(',')[0].strip().split('@')[0]
                    
                    if primary_og not in annotations:
                        annotations[primary_og] = {
                            'description': str(row['Description']),
                            'cog_category': str(row['COG_category']),
                            'kegg_pathway': str(row['KEGG_Pathway']),
                            'preferred_name': str(row['Preferred_name'])
                        }
                        
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return annotations

def analyze_evolutionary_patterns(categories_df, pa_matrix, distance_df):
    """Analyze evolutionary patterns based on gene categories and phylogeny"""
    
    # Find most closely related strain pairs
    distance_matrix = distance_df.values
    np.fill_diagonal(distance_matrix, np.inf)  # Ignore self-comparisons
    
    # Get closest pairs
    min_distance_idx = np.unravel_index(np.argmin(distance_matrix), distance_matrix.shape)
    closest_pair = (distance_df.index[min_distance_idx[0]], distance_df.index[min_distance_idx[1]])
    closest_distance = distance_matrix[min_distance_idx]
    
    print(f"Closest strain pair: {closest_pair[0]} <-> {closest_pair[1]} (distance: {closest_distance:.3f})")
    
    # Analyze gene content differences between closest strains
    strain1_genes = set(pa_matrix.loc[pa_matrix[closest_pair[0]] == 1].index)
    strain2_genes = set(pa_matrix.loc[pa_matrix[closest_pair[1]] == 1].index)
    
    shared_genes = strain1_genes & strain2_genes
    strain1_specific = strain1_genes - strain2_genes
    strain2_specific = strain2_genes - strain1_genes
    
    evolutionary_events = {
        'closest_pair': closest_pair,
        'distance': closest_distance,
        'shared_genes': len(shared_genes),
        'strain1_specific': len(strain1_specific),
        'strain2_specific': len(strain2_specific),
        'strain1_unique_genes': list(strain1_specific),
        'strain2_unique_genes': list(strain2_specific)
    }
    
    return evolutionary_events

def analyze_gene_category_evolution(categories_df, pa_matrix):
    """Analyze how different gene categories evolve"""
    
    category_evolution = {}
    
    for category in ['core', 'soft_core', 'shell', 'unique']:
        category_genes = categories_df[categories_df['category'] == category]
        
        if len(category_genes) == 0:
            continue
        
        # Calculate presence patterns
        presence_patterns = []
        for _, row in category_genes.iterrows():
            gene_family = row['ortholog_group']
            if gene_family in pa_matrix.index:
                presence_pattern = pa_matrix.loc[gene_family].values
                presence_patterns.append(presence_pattern)
        
        if presence_patterns:
            presence_array = np.array(presence_patterns)
            
            category_evolution[category] = {
                'gene_count': len(category_genes),
                'mean_presence': np.mean(presence_array),
                'presence_variance': np.var(presence_array),
                'evolutionary_volatility': np.std(presence_array, axis=0).mean()
            }
    
    return category_evolution

def analyze_strain_evolutionary_profiles(pa_matrix, distance_df):
    """Create evolutionary profile for each strain"""
    
    strain_profiles = {}
    
    for strain in pa_matrix.columns:
        # Count genes present in this strain
        present_genes = pa_matrix.loc[pa_matrix[strain] == 1].index.tolist()
        
        # Calculate average distance to other strains
        avg_distance = distance_df.loc[strain].mean()
        
        # Count unique genes (present only in this strain)
        unique_genes = []
        for gene in present_genes:
            if pa_matrix.loc[gene].sum() == 1:  # Only in this strain
                unique_genes.append(gene)
        
        # Count rare genes (present in ≤2 strains)
        rare_genes = []
        for gene in present_genes:
            if pa_matrix.loc[gene].sum() <= 2:
                rare_genes.append(gene)
        
        strain_profiles[strain] = {
            'total_genes': len(present_genes),
            'unique_genes': len(unique_genes),
            'rare_genes': len(rare_genes),
            'avg_phylo_distance': avg_distance,
            'unique_gene_list': unique_genes,
            'rare_gene_list': rare_genes
        }
    
    return strain_profiles

def analyze_functional_evolution(strain_profiles, annotations):
    """Analyze functional categories of evolved genes"""
    
    functional_evolution = {}
    
    for strain, profile in strain_profiles.items():
        unique_functions = Counter()
        rare_functions = Counter()
        
        # Analyze unique genes
        for gene in profile['unique_gene_list']:
            if gene in annotations:
                cog = annotations[gene]['cog_category']
                if cog != 'nan' and cog != '-':
                    for cog_cat in cog:
                        if cog_cat.isalpha():
                            unique_functions[cog_cat] += 1
        
        # Analyze rare genes
        for gene in profile['rare_gene_list']:
            if gene in annotations:
                cog = annotations[gene]['cog_category']
                if cog != 'nan' and cog != '-':
                    for cog_cat in cog:
                        if cog_cat.isalpha():
                            rare_functions[cog_cat] += 1
        
        functional_evolution[strain] = {
            'unique_cog_categories': dict(unique_functions),
            'rare_cog_categories': dict(rare_functions)
        }
    
    return functional_evolution

def create_evolution_visualizations(strain_profiles, category_evolution, functional_evolution, output_dir):
    """Create visualizations for evolutionary analysis"""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Strain evolutionary metrics
    strains = list(strain_profiles.keys())
    display_names = [s.split('_')[1] if '_' in s else s[:8] for s in strains]
    
    unique_counts = [strain_profiles[s]['unique_genes'] for s in strains]
    rare_counts = [strain_profiles[s]['rare_genes'] for s in strains]
    
    x = np.arange(len(strains))
    width = 0.35
    
    axes[0, 0].bar(x - width/2, unique_counts, width, label='Unique genes', color='red', alpha=0.7)
    axes[0, 0].bar(x + width/2, rare_counts, width, label='Rare genes (≤2 strains)', color='orange', alpha=0.7)
    axes[0, 0].set_xlabel('Strain')
    axes[0, 0].set_ylabel('Number of Genes')
    axes[0, 0].set_title('Strain-Specific and Rare Gene Content')
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(display_names, rotation=45, ha='right')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Genome size vs evolutionary distance
    genome_sizes = [strain_profiles[s]['total_genes'] for s in strains]
    avg_distances = [strain_profiles[s]['avg_phylo_distance'] for s in strains]
    
    axes[0, 1].scatter(avg_distances, genome_sizes, alpha=0.7, s=100)
    for i, strain in enumerate(display_names):
        axes[0, 1].annotate(strain, (avg_distances[i], genome_sizes[i]), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    axes[0, 1].set_xlabel('Average Phylogenetic Distance')
    axes[0, 1].set_ylabel('Genome Size (genes)')
    axes[0, 1].set_title('Genome Size vs Phylogenetic Distance')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Gene category evolutionary volatility
    if category_evolution:
        categories = list(category_evolution.keys())
        volatilities = [category_evolution[cat]['evolutionary_volatility'] for cat in categories]
        
        axes[1, 0].bar(categories, volatilities, color='lightblue', alpha=0.7)
        axes[1, 0].set_xlabel('Gene Category')
        axes[1, 0].set_ylabel('Evolutionary Volatility')
        axes[1, 0].set_title('Evolutionary Volatility by Gene Category')
        axes[1, 0].tick_params(axis='x', rotation=45)
    
    # 4. Functional categories of unique genes
    all_unique_cogs = Counter()
    for strain, funcs in functional_evolution.items():
        for cog, count in funcs['unique_cog_categories'].items():
            all_unique_cogs[cog] += count
    
    if all_unique_cogs:
        top_cogs = dict(all_unique_cogs.most_common(8))
        axes[1, 1].bar(top_cogs.keys(), top_cogs.values(), color='lightgreen', alpha=0.7)
        axes[1, 1].set_xlabel('COG Category')
        axes[1, 1].set_ylabel('Number of Unique Genes')
        axes[1, 1].set_title('Functional Categories of Strain-Specific Genes')
        axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / "evolutionary_patterns.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    phylo_dir = results_dir / "phylogenetic_analysis"
    pangenome_dir = results_dir / "pangenome_analysis"
    output_dir = results_dir / "evolutionary_dynamics"
    output_dir.mkdir(exist_ok=True)
    
    print("Starting simplified gene evolution analysis...")
    
    # Step 1: Load data
    print("\n1. Loading pan-genome and phylogenetic data...")
    categories_df, pa_matrix = load_pangenome_data(pangenome_dir)
    if categories_df is None:
        return
    
    distance_df = load_phylogenetic_relationships(phylo_dir)
    if distance_df is None:
        return
    
    print(f"Loaded data for {len(pa_matrix.columns)} strains and {len(pa_matrix)} gene families")
    
    # Step 2: Load annotations
    print("\n2. Loading functional annotations...")
    annotations = load_functional_annotations(results_dir)
    print(f"Loaded annotations for {len(annotations)} gene families")
    
    # Step 3: Analyze evolutionary patterns
    print("\n3. Analyzing evolutionary patterns...")
    evolutionary_events = analyze_evolutionary_patterns(categories_df, pa_matrix, distance_df)
    
    # Step 4: Analyze gene category evolution
    print("\n4. Analyzing gene category evolutionary patterns...")
    category_evolution = analyze_gene_category_evolution(categories_df, pa_matrix)
    
    # Step 5: Create strain evolutionary profiles
    print("\n5. Creating strain evolutionary profiles...")
    strain_profiles = analyze_strain_evolutionary_profiles(pa_matrix, distance_df)
    
    # Step 6: Analyze functional evolution
    print("\n6. Analyzing functional evolution patterns...")
    functional_evolution = analyze_functional_evolution(strain_profiles, annotations)
    
    # Step 7: Create visualizations
    print("\n7. Creating evolutionary visualizations...")
    create_evolution_visualizations(strain_profiles, category_evolution, functional_evolution, output_dir)
    
    # Print results
    print("\n=== EVOLUTIONARY DYNAMICS ANALYSIS ===")
    
    print(f"\nClosest Strain Pair Analysis:")
    print(f"  Pair: {evolutionary_events['closest_pair'][0]} <-> {evolutionary_events['closest_pair'][1]}")
    print(f"  Phylogenetic distance: {evolutionary_events['distance']:.3f}")
    print(f"  Shared genes: {evolutionary_events['shared_genes']}")
    print(f"  {evolutionary_events['closest_pair'][0]} specific: {evolutionary_events['strain1_specific']}")
    print(f"  {evolutionary_events['closest_pair'][1]} specific: {evolutionary_events['strain2_specific']}")
    
    print(f"\nGene Category Evolutionary Patterns:")
    for category, data in category_evolution.items():
        print(f"  {category.title()}:")
        print(f"    Gene families: {data['gene_count']}")
        print(f"    Mean presence: {data['mean_presence']:.2f}")
        print(f"    Evolutionary volatility: {data['evolutionary_volatility']:.3f}")
    
    print(f"\nStrain Evolutionary Profiles:")
    # Sort by unique gene content
    sorted_strains = sorted(strain_profiles.items(), key=lambda x: x[1]['unique_genes'], reverse=True)
    
    for strain, profile in sorted_strains:
        print(f"  {strain}:")
        print(f"    Total genes: {profile['total_genes']}")
        print(f"    Unique genes: {profile['unique_genes']}")
        print(f"    Rare genes: {profile['rare_genes']}")
        print(f"    Avg phylo distance: {profile['avg_phylo_distance']:.3f}")
    
    print(f"\nFunctional Evolution Summary:")
    total_unique_genes = sum(profile['unique_genes'] for profile in strain_profiles.values())
    total_rare_genes = sum(profile['rare_genes'] for profile in strain_profiles.values())
    
    print(f"  Total strain-specific genes: {total_unique_genes}")
    print(f"  Total rare genes (≤2 strains): {total_rare_genes}")
    
    # Most common functional categories in unique genes
    all_unique_cogs = Counter()
    for strain, funcs in functional_evolution.items():
        for cog, count in funcs['unique_cog_categories'].items():
            all_unique_cogs[cog] += count
    
    print(f"  Top COG categories in unique genes: {dict(all_unique_cogs.most_common(5))}")
    
    # Save results
    strain_summary = []
    for strain, profile in strain_profiles.items():
        strain_summary.append({
            'strain': strain,
            'total_genes': profile['total_genes'],
            'unique_genes': profile['unique_genes'],
            'rare_genes': profile['rare_genes'],
            'avg_phylo_distance': profile['avg_phylo_distance']
        })
    
    strain_df = pd.DataFrame(strain_summary)
    strain_df.to_csv(output_dir / "strain_evolutionary_profiles.csv", index=False)
    
    # Gene evolution summary
    gene_evolution_summary = []
    for gene_family in pa_matrix.index:
        presence_count = pa_matrix.loc[gene_family].sum()
        
        gene_evolution_summary.append({
            'gene_family': gene_family,
            'presence_count': presence_count,
            'category': 'unique' if presence_count == 1 else 'rare' if presence_count <= 2 else 'common',
            'annotation': annotations.get(gene_family, {}).get('description', 'Unknown')
        })
    
    gene_df = pd.DataFrame(gene_evolution_summary)
    gene_df.to_csv(output_dir / "gene_evolution_summary.csv", index=False)
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - evolutionary_patterns.png: Evolutionary pattern visualizations")
    print("  - strain_evolutionary_profiles.csv: Per-strain evolutionary metrics")
    print("  - gene_evolution_summary.csv: Per-gene family evolutionary patterns")
    
    return strain_profiles, category_evolution, functional_evolution

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
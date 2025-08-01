#!/usr/bin/env python3
"""
Specialized Functional Analysis of Ignatzschineria MAGs
- Decomposition genes (CAZy)
- Nitrogen cycling genes
- Insect interaction genes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter

def analyze_specialized_genes(results_dir):
    """Analyze specialized gene categories across MAGs"""
    
    # Get all MAG identifiers
    decomp_dir = results_dir / "decomposition_analysis/decomposition_genes"
    nitrogen_dir = results_dir / "decomposition_analysis/nitrogen_genes"
    insect_dir = results_dir / "decomposition_analysis/insect_genes"
    
    # Get MAG IDs from decomposition files
    decomp_files = list(decomp_dir.glob("*_decomposition_genes.tsv"))
    mag_ids = [f.name.split('_genomic_decomposition_genes.tsv')[0] for f in decomp_files]
    
    print(f"Analyzing specialized genes for {len(mag_ids)} MAGs")
    
    # Initialize results
    results = {
        'decomposition': {},
        'nitrogen': {},
        'insect': {}
    }
    
    # Analyze each category
    for mag_id in mag_ids:
        print(f"Processing {mag_id}...")
        
        # Decomposition genes
        decomp_file = decomp_dir / f"{mag_id}_genomic_decomposition_genes.tsv"
        if decomp_file.exists():
            decomp_df = pd.read_csv(decomp_file, sep='\t', header=None)
            results['decomposition'][mag_id] = len(decomp_df)
        else:
            results['decomposition'][mag_id] = 0
            
        # Nitrogen cycling genes
        nitrogen_file = nitrogen_dir / f"{mag_id}_genomic_nitrogen_cycling_genes.tsv"
        if nitrogen_file.exists():
            nitrogen_df = pd.read_csv(nitrogen_file, sep='\t', header=None)
            results['nitrogen'][mag_id] = len(nitrogen_df)
        else:
            results['nitrogen'][mag_id] = 0
            
        # Insect interaction genes
        insect_file = insect_dir / f"{mag_id}_genomic_insect_interaction_genes.tsv"
        if insect_file.exists():
            insect_df = pd.read_csv(insect_file, sep='\t', header=None)
            results['insect'][mag_id] = len(insect_df)
        else:
            results['insect'][mag_id] = 0
    
    return results, mag_ids

def analyze_cazy_detailed(results_dir, mag_ids):
    """Detailed analysis of CAZy enzyme families"""
    cazy_dir = results_dir / "decomposition_analysis/decomposition_genes"
    
    cazy_families = defaultdict(dict)
    
    for mag_id in mag_ids:
        cazy_file = cazy_dir / f"{mag_id}_genomic_cazy_genes.tsv"
        if cazy_file.exists():
            cazy_df = pd.read_csv(cazy_file, sep='\t', header=None)
            
            # Extract CAZy families from annotations (assuming similar format to eggNOG)
            family_counts = Counter()
            for _, row in cazy_df.iterrows():
                # This assumes CAZy info is in a specific column - adjust as needed
                if len(row) > 18:  # CAZy column should be around position 18
                    cazy_info = str(row[18]) if pd.notna(row[18]) else ""
                    if cazy_info and cazy_info != '-':
                        for family in cazy_info.split(','):
                            family = family.strip()
                            if family:
                                family_counts[family] += 1
            
            for family, count in family_counts.items():
                cazy_families[family][mag_id] = count
    
    return cazy_families

def create_visualizations(results, mag_ids, output_dir):
    """Create comparative visualizations"""
    
    plt.style.use('seaborn-v0_8')
    
    # Specialized gene categories comparison
    categories = ['decomposition', 'nitrogen', 'insect']
    data_matrix = []
    
    for category in categories:
        row = [results[category].get(mag_id, 0) for mag_id in mag_ids]
        data_matrix.append(row)
    
    # Create DataFrame for heatmap
    spec_df = pd.DataFrame(data_matrix, columns=mag_ids, index=categories)
    
    # Specialized genes heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(spec_df, annot=True, fmt='.0f', cmap='Greens', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Specialized Gene Categories Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Gene Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "specialized_genes_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Stacked bar chart
    plt.figure(figsize=(12, 8))
    bottom = np.zeros(len(mag_ids))
    colors = ['#2E8B57', '#4682B4', '#CD853F']  # Green, Blue, Brown
    
    for i, category in enumerate(categories):
        values = [results[category].get(mag_id, 0) for mag_id in mag_ids]
        plt.bar(range(len(mag_ids)), values, bottom=bottom, 
                label=category.title(), color=colors[i], alpha=0.8)
        bottom += values
    
    plt.xlabel('MAG')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Specialized Genes Across Ignatzschineria MAGs')
    plt.xticks(range(len(mag_ids)), mag_ids, rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "specialized_genes_stacked.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    return spec_df

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    output_dir = results_dir / "comparative_analysis"
    output_dir.mkdir(exist_ok=True)
    
    # Analyze specialized genes
    results, mag_ids = analyze_specialized_genes(results_dir)
    
    # Analyze CAZy families in detail
    cazy_families = analyze_cazy_detailed(results_dir, mag_ids)
    
    print("\n=== SPECIALIZED FUNCTIONAL ANALYSIS ===")
    
    # Summary statistics
    for category in ['decomposition', 'nitrogen', 'insect']:
        total_genes = sum(results[category].values())
        avg_genes = np.mean(list(results[category].values()))
        std_genes = np.std(list(results[category].values()))
        
        print(f"\n{category.title()} genes:")
        print(f"  Total across all MAGs: {total_genes}")
        print(f"  Average per MAG: {avg_genes:.1f} ± {std_genes:.1f}")
        
        # MAG with most genes in this category
        max_mag = max(results[category], key=results[category].get)
        max_count = results[category][max_mag]
        print(f"  Highest count: {max_mag} ({max_count} genes)")
    
    # Create visualizations
    spec_df = create_visualizations(results, mag_ids, output_dir)
    
    # Save detailed results
    spec_df.to_csv(output_dir / "specialized_genes_comparison.csv")
    
    # CAZy analysis if data available
    if cazy_families:
        cazy_df = pd.DataFrame(cazy_families).fillna(0).T
        if len(cazy_df) > 0:
            cazy_df.to_csv(output_dir / "detailed_cazy_families.csv")
            
            print(f"\nCAZy Enzyme Families:")
            print(f"  Total unique families identified: {len(cazy_df)}")
            
            # Most abundant families
            family_totals = cazy_df.sum(axis=1).sort_values(ascending=False)
            print(f"  Top CAZy families:")
            for family in family_totals.head(5).index:
                print(f"    {family}: {family_totals[family]:.0f} total occurrences")
    
    # Functional capacity analysis
    print(f"\n=== FUNCTIONAL CAPACITY COMPARISON ===")
    
    # Calculate total specialized genes per MAG
    total_specialized = {}
    for mag_id in mag_ids:
        total = (results['decomposition'].get(mag_id, 0) + 
                results['nitrogen'].get(mag_id, 0) + 
                results['insect'].get(mag_id, 0))
        total_specialized[mag_id] = total
    
    # Rank MAGs by specialized function
    ranked_mags = sorted(total_specialized.items(), key=lambda x: x[1], reverse=True)
    
    print("MAGs ranked by total specialized genes:")
    for i, (mag_id, count) in enumerate(ranked_mags, 1):
        print(f"  {i}. {mag_id}: {count} specialized genes")
        decomp = results['decomposition'].get(mag_id, 0)
        nitrogen = results['nitrogen'].get(mag_id, 0) 
        insect = results['insect'].get(mag_id, 0)
        print(f"     Decomposition: {decomp}, Nitrogen: {nitrogen}, Insect: {insect}")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - specialized_genes_comparison.csv")
    print("  - specialized_genes_heatmap.png")
    print("  - specialized_genes_stacked.png")
    if cazy_families:
        print("  - detailed_cazy_families.csv")

if __name__ == "__main__":
    main()
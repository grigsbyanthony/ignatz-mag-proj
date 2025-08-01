#!/usr/bin/env python3
"""
Comparative Functional Analysis of Ignatzschineria MAGs
Author: Claude Code Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict, Counter

def parse_eggnog_annotations(file_path):
    """Parse eggNOG annotation files and extract functional information"""
    df = pd.read_csv(file_path, sep='\t', comment='#', 
                     names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                           'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                           'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                           'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                           'BiGG_Reaction', 'PFAMs'])
    return df

def extract_mag_id(file_path):
    """Extract MAG identifier from file path"""
    filename = Path(file_path).name
    return filename.split('_genomic_eggnog')[0]

def analyze_cog_categories(df):
    """Analyze COG functional categories"""
    cog_counts = Counter()
    for cog_cat in df['COG_category'].dropna():
        if cog_cat != '-':
            for cat in cog_cat:
                cog_counts[cat] += 1
    return cog_counts

def analyze_kegg_pathways(df):
    """Analyze KEGG pathway distribution"""
    pathway_counts = Counter()
    for pathways in df['KEGG_Pathway'].dropna():
        if pathways != '-':
            for pathway in pathways.split(','):
                pathway = pathway.strip()
                if pathway.startswith('ko'):
                    pathway_counts[pathway] += 1
    return pathway_counts

def analyze_cazy_families(df):
    """Analyze CAZy enzyme families"""
    cazy_counts = Counter()
    for cazy in df['CAZy'].dropna():
        if cazy != '-':
            for family in cazy.split(','):
                family = family.strip()
                cazy_counts[family] += 1
    return cazy_counts

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    
    # Get all annotation files
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    print(f"Found {len(annotation_files)} MAG annotation files")
    
    # Initialize data structures
    mag_data = {}
    all_cog_categories = defaultdict(dict)
    all_kegg_pathways = defaultdict(dict)
    all_cazy_families = defaultdict(dict)
    
    # Process each MAG
    for file_path in annotation_files:
        mag_id = extract_mag_id(file_path)
        print(f"Processing {mag_id}...")
        
        # Parse annotations
        df = parse_eggnog_annotations(file_path)
        mag_data[mag_id] = df
        
        # Analyze functional categories
        cog_counts = analyze_cog_categories(df)
        kegg_counts = analyze_kegg_pathways(df)
        cazy_counts = analyze_cazy_families(df)
        
        # Store results
        for cog, count in cog_counts.items():
            all_cog_categories[cog][mag_id] = count
        
        for pathway, count in kegg_counts.items():
            all_kegg_pathways[pathway][mag_id] = count
            
        for family, count in cazy_counts.items():
            all_cazy_families[family][mag_id] = count
    
    # Create comparison DataFrames
    print("\nCreating comparison matrices...")
    
    # COG categories comparison
    cog_df = pd.DataFrame(all_cog_categories).fillna(0).T
    cog_df.index.name = 'COG_Category'
    
    # KEGG pathways comparison (top 20 most common)
    kegg_df = pd.DataFrame(all_kegg_pathways).fillna(0).T
    kegg_totals = kegg_df.sum(axis=1).sort_values(ascending=False)
    kegg_df = kegg_df.loc[kegg_totals.head(20).index]
    kegg_df.index.name = 'KEGG_Pathway'
    
    # CAZy families comparison
    cazy_df = pd.DataFrame(all_cazy_families).fillna(0).T
    cazy_df.index.name = 'CAZy_Family'
    
    # Calculate basic statistics
    print("\n=== COMPARATIVE FUNCTIONAL ANALYSIS RESULTS ===")
    print(f"\nTotal MAGs analyzed: {len(mag_data)}")
    
    # Gene counts per MAG
    gene_counts = {mag_id: len(df) for mag_id, df in mag_data.items()}
    print(f"\nGene counts per MAG:")
    for mag_id, count in sorted(gene_counts.items()):
        print(f"  {mag_id}: {count:,} genes")
    
    # COG category analysis
    print(f"\nCOG Functional Categories (present in all MAGs):")
    core_cogs = cog_df.columns[(cog_df > 0).all()]
    if len(core_cogs) > 0:
        print(f"  Core COG categories: {', '.join(core_cogs)}")
    else:
        print("  No COG categories present in all MAGs")
    
    # Most variable COG categories
    cog_cv = (cog_df.std(axis=1) / cog_df.mean(axis=1)).sort_values(ascending=False)
    print(f"\nMost variable COG categories:")
    for cog in cog_cv.head(5).index:
        print(f"  {cog}: CV = {cog_cv[cog]:.2f}")
    
    # CAZy analysis
    if len(cazy_df) > 0:
        print(f"\nCAZy Enzyme Families:")
        print(f"  Total unique CAZy families: {len(cazy_df)}")
        cazy_totals = cazy_df.sum(axis=1).sort_values(ascending=False)
        print(f"  Most abundant CAZy families:")
        for family in cazy_totals.head(5).index:
            print(f"    {family}: {cazy_totals[family]} total occurrences")
    
    # Save results
    output_dir = results_dir / "comparative_analysis"
    output_dir.mkdir(exist_ok=True)
    
    cog_df.to_csv(output_dir / "cog_comparison.csv")
    kegg_df.to_csv(output_dir / "kegg_comparison.csv")
    if len(cazy_df) > 0:
        cazy_df.to_csv(output_dir / "cazy_comparison.csv")
    
    # Create visualizations
    plt.style.use('seaborn-v0_8')
    
    # COG categories heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(cog_df, annot=True, fmt='.0f', cmap='Blues', cbar_kws={'label': 'Gene Count'})
    plt.title('COG Functional Categories Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('COG Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "cog_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gene count comparison
    plt.figure(figsize=(10, 6))
    mag_ids = list(gene_counts.keys())
    counts = list(gene_counts.values())
    plt.bar(range(len(mag_ids)), counts, color='skyblue', edgecolor='navy', alpha=0.7)
    plt.xlabel('MAG')
    plt.ylabel('Number of Genes')
    plt.title('Gene Count Comparison Across Ignatzschineria MAGs')
    plt.xticks(range(len(mag_ids)), mag_ids, rotation=45, ha='right')
    for i, count in enumerate(counts):
        plt.text(i, count + max(counts)*0.01, f'{count:,}', ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(output_dir / "gene_counts.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # CAZy families heatmap (if data exists)
    if len(cazy_df) > 0 and len(cazy_df) <= 30:  # Only plot if reasonable number of families
        plt.figure(figsize=(12, max(6, len(cazy_df)*0.3)))
        sns.heatmap(cazy_df, annot=True, fmt='.0f', cmap='Oranges', cbar_kws={'label': 'Gene Count'})
        plt.title('CAZy Enzyme Families Across Ignatzschineria MAGs')
        plt.xlabel('MAG')
        plt.ylabel('CAZy Family')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / "cazy_heatmap.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - cog_comparison.csv: COG category counts per MAG")
    print("  - kegg_comparison.csv: KEGG pathway counts per MAG") 
    print("  - cazy_comparison.csv: CAZy family counts per MAG")
    print("  - cog_heatmap.png: COG categories visualization")
    print("  - gene_counts.png: Gene count comparison")
    if len(cazy_df) > 0:
        print("  - cazy_heatmap.png: CAZy families visualization")

if __name__ == "__main__":
    main()
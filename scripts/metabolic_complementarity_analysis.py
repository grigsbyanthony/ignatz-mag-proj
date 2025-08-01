#!/usr/bin/env python3
"""
Metabolic Complementarity Analysis for Ignatzschineria MAGs
Focus on essential nutrients and metabolic support for insect hosts
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import re

def analyze_essential_nutrients(df):
    """Analyze genes involved in essential nutrient biosynthesis"""
    nutrients = {
        'amino_acids': [],
        'vitamins': [],
        'cofactors': [],
        'fatty_acids': [],
        'nucleotides': []
    }
    
    for idx, row in df.iterrows():
        description = str(row['Description']).lower()
        kegg_ko = str(row['KEGG_ko']).lower()
        kegg_pathway = str(row['KEGG_Pathway']).lower()
        
        # Amino acid biosynthesis
        if any(term in description for term in ['amino acid', 'arginine', 'lysine', 'methionine', 'tryptophan', 
                                               'phenylalanine', 'tyrosine', 'histidine', 'threonine', 'valine',
                                               'leucine', 'isoleucine', 'aromatic amino']):
            nutrients['amino_acids'].append(row)
        
        # Vitamin biosynthesis
        elif any(term in description for term in ['biotin', 'thiamine', 'folate', 'pantothen', 'riboflavin', 
                                                'pyridoxal', 'cobalamin', 'vitamin', 'folic acid']):
            nutrients['vitamins'].append(row)
        
        # Cofactor biosynthesis
        elif any(term in description for term in ['molybdopterin', 'heme', 'coenzyme', 'nad', 'fad', 'cofactor']):
            nutrients['cofactors'].append(row)
        
        # Fatty acid biosynthesis
        elif any(term in description for term in ['fatty acid', 'acetyl-coa carboxylase', 'malonyl-coa']):
            nutrients['fatty_acids'].append(row)
        
        # Nucleotide biosynthesis  
        elif any(term in description for term in ['purine', 'pyrimidine', 'nucleotide', 'gmp', 'amp', 'cmp', 'ump']):
            nutrients['nucleotides'].append(row)
    
    return nutrients

def analyze_metabolic_pathways(df):
    """Analyze complete metabolic pathways"""
    pathways = {
        'central_metabolism': [],
        'biosynthesis': [],
        'degradation': [],
        'energy_metabolism': []
    }
    
    for idx, row in df.iterrows():
        kegg_pathway = str(row['KEGG_Pathway']).lower()
        description = str(row['Description']).lower()
        
        # Central metabolism
        if any(pathway in kegg_pathway for pathway in ['map01100', 'map00010', 'map00020', 'map00030']):
            pathways['central_metabolism'].append(row)
        
        # Biosynthetic pathways
        elif any(pathway in kegg_pathway for pathway in ['map00220', 'map00230', 'map00240', 'map00250',
                                                        'map00260', 'map00270', 'map00290', 'map00300']):
            pathways['biosynthesis'].append(row)
        
        # Degradation pathways
        elif any(pathway in kegg_pathway for pathway in ['map00380', 'map00410', 'map00620', 'map00650']):
            pathways['degradation'].append(row)
        
        # Energy metabolism
        elif any(pathway in kegg_pathway for pathway in ['map00190', 'map00710', 'map00720']):
            pathways['energy_metabolism'].append(row)
    
    return pathways

def analyze_host_beneficial_functions(df):
    """Analyze functions that could benefit insect hosts"""
    beneficial = {
        'detoxification': [],
        'immune_support': [],
        'digestion_aid': [],
        'antimicrobial': []
    }
    
    for idx, row in df.iterrows():
        description = str(row['Description']).lower()
        pfams = str(row['PFAMs']).lower()
        
        # Detoxification
        if any(term in description for term in ['oxidase', 'reductase', 'dehydrogenase', 'detox', 'xenobiotic']):
            beneficial['detoxification'].append(row)
        
        # Immune support (antioxidants, etc.)
        elif any(term in description for term in ['catalase', 'superoxide', 'peroxidase', 'glutathione']):
            beneficial['immune_support'].append(row)
        
        # Digestion aid (enzymes)
        elif any(term in description for term in ['peptidase', 'protease', 'lipase', 'amylase', 'cellulase']):
            beneficial['digestion_aid'].append(row)
        
        # Antimicrobial compounds
        elif any(term in description for term in ['antibiotic', 'antimicrobial', 'bacteriocin']):
            beneficial['antimicrobial'].append(row)
    
    return beneficial

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    output_dir = results_dir / "comparative_analysis"
    output_dir.mkdir(exist_ok=True)
    
    # Get all annotation files
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    mag_ids = [f.name.split('_genomic_eggnog.emapper.annotations')[0] for f in annotation_files]
    
    print(f"Analyzing metabolic complementarity for {len(mag_ids)} MAGs")
    
    # Initialize data structures
    all_nutrients = {}
    all_pathways = {}
    all_beneficial = {}
    nutrient_summary = defaultdict(dict)
    pathway_summary = defaultdict(dict)
    beneficial_summary = defaultdict(dict)
    
    # Process each MAG
    for file_path in annotation_files:
        mag_id = file_path.name.split('_genomic_eggnog.emapper.annotations')[0]
        print(f"Processing {mag_id}...")
        
        # Parse annotations
        df = pd.read_csv(file_path, sep='\t', comment='#', 
                        names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                               'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                               'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                               'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                               'BiGG_Reaction', 'PFAMs'])
        
        # Analyze different aspects
        nutrients = analyze_essential_nutrients(df)
        pathways = analyze_metabolic_pathways(df)
        beneficial = analyze_host_beneficial_functions(df)
        
        all_nutrients[mag_id] = nutrients
        all_pathways[mag_id] = pathways
        all_beneficial[mag_id] = beneficial
        
        # Count for summaries
        for nutrient_type, genes in nutrients.items():
            nutrient_summary[nutrient_type][mag_id] = len(genes)
        
        for pathway_type, genes in pathways.items():
            pathway_summary[pathway_type][mag_id] = len(genes)
            
        for beneficial_type, genes in beneficial.items():
            beneficial_summary[beneficial_type][mag_id] = len(genes)
    
    print("\n=== METABOLIC COMPLEMENTARITY ANALYSIS ===")
    
    # Essential nutrient biosynthesis
    print(f"\nEssential Nutrient Biosynthesis:")
    for nutrient_type in nutrient_summary:
        total_genes = sum(nutrient_summary[nutrient_type].values())
        avg_per_mag = np.mean(list(nutrient_summary[nutrient_type].values()))
        print(f"  {nutrient_type.replace('_', ' ').title()}: {total_genes} total ({avg_per_mag:.1f} avg per MAG)")
    
    # Metabolic pathway completeness
    print(f"\nMetabolic Pathway Categories:")
    for pathway_type in pathway_summary:
        total_genes = sum(pathway_summary[pathway_type].values())
        avg_per_mag = np.mean(list(pathway_summary[pathway_type].values()))
        print(f"  {pathway_type.replace('_', ' ').title()}: {total_genes} total ({avg_per_mag:.1f} avg per MAG)")
    
    # Host-beneficial functions
    print(f"\nHost-Beneficial Functions:")
    for beneficial_type in beneficial_summary:
        total_genes = sum(beneficial_summary[beneficial_type].values())
        avg_per_mag = np.mean(list(beneficial_summary[beneficial_type].values()))
        print(f"  {beneficial_type.replace('_', ' ').title()}: {total_genes} total ({avg_per_mag:.1f} avg per MAG)")
    
    # Key findings
    print(f"\n=== KEY METABOLIC CONTRIBUTIONS ===")
    
    # Amino acid biosynthesis leaders
    aa_counts = {mag: len(all_nutrients[mag].get('amino_acids', [])) for mag in mag_ids}
    print(f"\nAmino Acid Biosynthesis Leaders:")
    for mag_id, count in sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:3]:
        print(f"  {mag_id}: {count} amino acid biosynthesis genes")
    
    # Vitamin biosynthesis
    vitamin_counts = {mag: len(all_nutrients[mag].get('vitamins', [])) for mag in mag_ids}
    print(f"\nVitamin Biosynthesis:")
    for mag_id, count in sorted(vitamin_counts.items(), key=lambda x: x[1], reverse=True)[:3]:
        print(f"  {mag_id}: {count} vitamin biosynthesis genes")
    
    # Create comparison matrices
    nutrients_df = pd.DataFrame(nutrient_summary).fillna(0).T
    pathways_df = pd.DataFrame(pathway_summary).fillna(0).T
    beneficial_df = pd.DataFrame(beneficial_summary).fillna(0).T
    
    # Save results
    nutrients_df.to_csv(output_dir / "nutrient_biosynthesis_comparison.csv")
    pathways_df.to_csv(output_dir / "metabolic_pathways_comparison.csv")
    beneficial_df.to_csv(output_dir / "host_beneficial_functions_comparison.csv")
    
    # Create visualizations
    plt.style.use('seaborn-v0_8')
    
    # Essential nutrients heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(nutrients_df, annot=True, fmt='.0f', cmap='Greens', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Essential Nutrient Biosynthesis Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Nutrient Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "nutrient_biosynthesis_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Host-beneficial functions heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(beneficial_df, annot=True, fmt='.0f', cmap='Oranges', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Host-Beneficial Functions Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Beneficial Function')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "host_beneficial_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Metabolic contribution summary
    plt.figure(figsize=(14, 8))
    
    # Calculate total metabolic contribution per MAG
    total_contribution = {}
    for mag_id in mag_ids:
        nutrients_total = sum(len(genes) for genes in all_nutrients[mag_id].values())
        beneficial_total = sum(len(genes) for genes in all_beneficial[mag_id].values())
        total_contribution[mag_id] = nutrients_total + beneficial_total
    
    # Create stacked bar chart
    nutrients_totals = [sum(len(genes) for genes in all_nutrients[mag].values()) for mag in mag_ids]
    beneficial_totals = [sum(len(genes) for genes in all_beneficial[mag].values()) for mag in mag_ids]
    
    x = range(len(mag_ids))
    plt.bar(x, nutrients_totals, label='Essential Nutrients', color='lightblue', alpha=0.8)
    plt.bar(x, beneficial_totals, bottom=nutrients_totals, label='Host-Beneficial Functions', 
            color='lightcoral', alpha=0.8)
    
    plt.xlabel('MAG')
    plt.ylabel('Number of Genes')
    plt.title('Metabolic Contribution to Host Across Ignatzschineria MAGs')
    plt.xticks(x, mag_ids, rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "metabolic_contribution_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Detailed analysis
    print(f"\n=== DETAILED METABOLIC ANALYSIS ===")
    
    # Most common biosynthetic functions
    all_functions = []
    for mag_id in mag_ids:
        for nutrient_type, genes in all_nutrients[mag_id].items():
            for gene in genes:
                desc = str(gene['Description'])
                if desc != 'nan' and desc != '-':
                    all_functions.append(desc)
    
    function_counts = Counter(all_functions)
    print(f"\nMost common biosynthetic functions:")
    for func, count in function_counts.most_common(10):
        print(f"  {func[:80]}{'...' if len(func) > 80 else ''}: {count}")
    
    # MAG rankings by metabolic contribution
    print(f"\nMAGs ranked by total metabolic contribution:")
    for i, (mag_id, count) in enumerate(sorted(total_contribution.items(), 
                                              key=lambda x: x[1], reverse=True), 1):
        nutrients_count = sum(len(genes) for genes in all_nutrients[mag_id].values())
        beneficial_count = sum(len(genes) for genes in all_beneficial[mag_id].values())
        print(f"  {i}. {mag_id}: {count} total genes")
        print(f"     Nutrient biosynthesis: {nutrients_count}, Host-beneficial: {beneficial_count}")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - nutrient_biosynthesis_comparison.csv")
    print("  - metabolic_pathways_comparison.csv")
    print("  - host_beneficial_functions_comparison.csv")
    print("  - nutrient_biosynthesis_heatmap.png")
    print("  - host_beneficial_heatmap.png")
    print("  - metabolic_contribution_summary.png")

if __name__ == "__main__":
    main()
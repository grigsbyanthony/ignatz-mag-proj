#!/usr/bin/env python3
"""
Host-Microbe Interaction Mechanism Analysis for Ignatzschineria MAGs
Focus on secretion systems, toxin-antitoxin systems, and host adaptation strategies
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import re

def parse_interaction_genes(file_path):
    """Parse insect interaction gene files"""
    df = pd.read_csv(file_path, sep='\t', header=None, 
                     names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                           'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                           'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                           'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                           'BiGG_Reaction', 'PFAMs'])
    return df

def classify_interaction_mechanisms(df):
    """Classify genes into interaction mechanism categories"""
    mechanisms = {
        'toxin_antitoxin': [],
        'secretion_transport': [],
        'resistance_systems': [],
        'membrane_systems': [],
        'stress_response': [],
        'metabolic_support': [],
        'cell_envelope': []
    }
    
    for idx, row in df.iterrows():
        description = str(row['Description']).lower()
        preferred_name = str(row['Preferred_name']).lower()
        pfams = str(row['PFAMs']).lower()
        
        # Toxin-Antitoxin Systems
        if any(term in description for term in ['toxin', 'antitoxin', 'abieii', 'yoeb', 'mazf']):
            mechanisms['toxin_antitoxin'].append(row)
        
        # Secretion and Transport Systems
        elif any(term in description for term in ['abc-type transport', 'secretion', 'efflux', 'transporter']):
            mechanisms['secretion_transport'].append(row)
        
        # Resistance Systems
        elif any(term in description for term in ['resistance', 'tellurite', 'copper', 'multidrug']):
            mechanisms['resistance_systems'].append(row)
        
        # Membrane Systems (including Mla system)
        elif any(term in description for term in ['membrane', 'lipoprotein', 'periplasmic', 'mla']):
            mechanisms['membrane_systems'].append(row)
        
        # Stress Response
        elif any(term in description for term in ['stress', 'oxidase', 'peroxidase', 'catalase']):
            mechanisms['stress_response'].append(row)
        
        # Cell Envelope Modification
        elif any(term in description for term in ['lipopolysaccharide', 'peptidoglycan', 'murein', 'glmu']):
            mechanisms['cell_envelope'].append(row)
        
        # Metabolic Support (amino acid, vitamin synthesis)
        elif any(term in description for term in ['amino acid', 'vitamin', 'cofactor', 'thiamine', 'biosynthetic']):
            mechanisms['metabolic_support'].append(row)
    
    return mechanisms

def analyze_secretion_systems(df):
    """Specific analysis of secretion systems"""
    secretion_systems = {
        'Type_I': [],
        'Type_II': [],
        'Type_III': [],
        'Type_IV': [],
        'Type_VI': [],
        'ABC_transporters': [],
        'RND_efflux': [],
        'Mla_system': []
    }
    
    for idx, row in df.iterrows():
        description = str(row['Description']).lower()
        pfams = str(row['PFAMs']).lower()
        
        if 'type iv ta system' in description:
            secretion_systems['Type_IV'].append(row)
        elif 'abc-type transport' in description:
            secretion_systems['ABC_transporters'].append(row)
        elif 'rnd' in description or 'acr_tran' in pfams:
            secretion_systems['RND_efflux'].append(row)
        elif any(mla in str(row['Preferred_name']).lower() for mla in ['mlad', 'mlae', 'mlac']):
            secretion_systems['Mla_system'].append(row)
    
    return secretion_systems

def analyze_host_adaptation_strategies(all_mechanisms):
    """Analyze host adaptation strategies across MAGs"""
    strategies = {
        'immune_evasion': ['toxin_antitoxin', 'membrane_systems', 'cell_envelope'],
        'stress_tolerance': ['resistance_systems', 'stress_response'],
        'metabolic_cooperation': ['metabolic_support', 'secretion_transport'],
        'persistence': ['toxin_antitoxin', 'resistance_systems']
    }
    
    mag_strategies = {}
    for mag_id, mechanisms in all_mechanisms.items():
        mag_strategies[mag_id] = {}
        for strategy, mechanism_types in strategies.items():
            count = sum(len(mechanisms.get(mech_type, [])) for mech_type in mechanism_types)
            mag_strategies[mag_id][strategy] = count
    
    return mag_strategies

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    insect_dir = results_dir / "decomposition_analysis/insect_genes"
    output_dir = results_dir / "comparative_analysis"
    output_dir.mkdir(exist_ok=True)
    
    # Get all insect interaction files
    insect_files = list(insect_dir.glob("*_insect_interaction_genes.tsv"))
    mag_ids = [f.name.split('_genomic_insect_interaction_genes.tsv')[0] for f in insect_files]
    
    print(f"Analyzing host-microbe interaction mechanisms for {len(mag_ids)} MAGs")
    
    # Initialize data structures
    all_mechanisms = {}
    all_secretion_systems = {}
    interaction_summary = defaultdict(dict)
    
    # Process each MAG
    for file_path in insect_files:
        mag_id = file_path.name.split('_genomic_insect_interaction_genes.tsv')[0]
        print(f"Processing {mag_id}...")
        
        # Parse interaction genes
        df = parse_interaction_genes(file_path)
        
        # Classify interaction mechanisms
        mechanisms = classify_interaction_mechanisms(df)
        all_mechanisms[mag_id] = mechanisms
        
        # Analyze secretion systems
        secretion_systems = analyze_secretion_systems(df)
        all_secretion_systems[mag_id] = secretion_systems
        
        # Count mechanisms
        for mech_type, genes in mechanisms.items():
            interaction_summary[mech_type][mag_id] = len(genes)
    
    # Analyze host adaptation strategies
    adaptation_strategies = analyze_host_adaptation_strategies(all_mechanisms)
    
    print("\n=== HOST-MICROBE INTERACTION ANALYSIS ===")
    
    # Summary statistics
    print(f"\nInteraction Mechanism Categories:")
    for mech_type in interaction_summary:
        total_genes = sum(interaction_summary[mech_type].values())
        avg_per_mag = np.mean(list(interaction_summary[mech_type].values()))
        print(f"  {mech_type.replace('_', ' ').title()}: {total_genes} total ({avg_per_mag:.1f} avg per MAG)")
    
    # Key findings
    print(f"\n=== KEY HOST-MICROBE INTERACTION MECHANISMS ===")
    
    # Toxin-Antitoxin Systems
    ta_counts = {mag: len(mechanisms.get('toxin_antitoxin', [])) for mag, mechanisms in all_mechanisms.items()}
    print(f"\nToxin-Antitoxin Systems:")
    for mag_id, count in sorted(ta_counts.items(), key=lambda x: x[1], reverse=True)[:3]:
        print(f"  {mag_id}: {count} TA systems")
    
    # Secretion Systems Analysis
    print(f"\nSecretion Systems:")
    for system_type in ['ABC_transporters', 'RND_efflux', 'Mla_system', 'Type_IV']:
        counts = [len(all_secretion_systems[mag].get(system_type, [])) for mag in mag_ids]
        total = sum(counts)
        if total > 0:
            print(f"  {system_type.replace('_', ' ')}: {total} total across all MAGs")
    
    # Host Adaptation Strategies
    print(f"\n=== HOST ADAPTATION STRATEGIES ===")
    strategy_totals = defaultdict(int)
    for mag_id, strategies in adaptation_strategies.items():
        for strategy, count in strategies.items():
            strategy_totals[strategy] += count
    
    for strategy, total in sorted(strategy_totals.items(), key=lambda x: x[1], reverse=True):
        avg_per_mag = total / len(mag_ids)
        print(f"  {strategy.replace('_', ' ').title()}: {total} genes ({avg_per_mag:.1f} avg per MAG)")
    
    # Create comparison matrices
    mechanisms_df = pd.DataFrame(interaction_summary).fillna(0).T
    strategies_df = pd.DataFrame(adaptation_strategies).fillna(0).T
    
    # Save detailed results
    mechanisms_df.to_csv(output_dir / "interaction_mechanisms_comparison.csv")
    strategies_df.to_csv(output_dir / "adaptation_strategies_comparison.csv")
    
    # Create visualizations
    plt.style.use('seaborn-v0_8')
    
    # Interaction mechanisms heatmap
    plt.figure(figsize=(14, 8))
    sns.heatmap(mechanisms_df, annot=True, fmt='.0f', cmap='Reds', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Host-Microbe Interaction Mechanisms Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Interaction Mechanism')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "interaction_mechanisms_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Adaptation strategies heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(strategies_df, annot=True, fmt='.0f', cmap='Blues', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Host Adaptation Strategies Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Adaptation Strategy')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "adaptation_strategies_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Stacked bar chart for mechanisms
    plt.figure(figsize=(14, 8))
    mechanisms_df_plot = mechanisms_df.T  # Transpose for plotting
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(mechanisms_df_plot.columns)))
    bottom = np.zeros(len(mechanisms_df_plot))
    
    for i, mechanism in enumerate(mechanisms_df_plot.columns):
        plt.bar(range(len(mechanisms_df_plot)), mechanisms_df_plot[mechanism], 
                bottom=bottom, label=mechanism.replace('_', ' ').title(), 
                color=colors[i], alpha=0.8)
        bottom += mechanisms_df_plot[mechanism]
    
    plt.xlabel('MAG')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Host-Microbe Interaction Mechanisms')
    plt.xticks(range(len(mechanisms_df_plot)), mechanisms_df_plot.index, rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / "interaction_mechanisms_stacked.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Detailed mechanism analysis
    print(f"\n=== DETAILED MECHANISM ANALYSIS ===")
    
    # Most common specific mechanisms
    all_descriptions = []
    for mag_id, mechanisms in all_mechanisms.items():
        for mech_type, genes in mechanisms.items():
            for gene in genes:
                desc = str(gene['Description'])
                if desc != 'nan' and desc != '-':
                    all_descriptions.append(desc)
    
    desc_counts = Counter(all_descriptions)
    print(f"\nMost common interaction mechanisms:")
    for desc, count in desc_counts.most_common(10):
        print(f"  {desc[:80]}{'...' if len(desc) > 80 else ''}: {count}")
    
    # MAG rankings by interaction complexity
    total_interactions = {mag: sum(len(genes) for genes in mechanisms.values()) 
                         for mag, mechanisms in all_mechanisms.items()}
    
    print(f"\nMAGs ranked by total interaction genes:")
    for i, (mag_id, count) in enumerate(sorted(total_interactions.items(), 
                                              key=lambda x: x[1], reverse=True), 1):
        print(f"  {i}. {mag_id}: {count} interaction genes")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - interaction_mechanisms_comparison.csv")
    print("  - adaptation_strategies_comparison.csv") 
    print("  - interaction_mechanisms_heatmap.png")
    print("  - adaptation_strategies_heatmap.png")
    print("  - interaction_mechanisms_stacked.png")

if __name__ == "__main__":
    main()
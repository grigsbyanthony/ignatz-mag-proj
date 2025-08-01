#!/usr/bin/env python3
"""
Detailed analysis of accessory genome functional content
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter

def load_pangenome_data(pangenome_dir):
    """Load pan-genome classification results"""
    categories_file = pangenome_dir / "pangenome_categories.csv"
    if not categories_file.exists():
        print(f"Pan-genome categories file not found: {categories_file}")
        return None
    
    pangenome_df = pd.read_csv(categories_file)
    return pangenome_df

def load_functional_annotations(results_dir):
    """Load functional annotations for all genes"""
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    all_annotations = {}
    
    for file_path in annotation_files:
        mag_id = file_path.name.split('_genomic_eggnog.emapper.annotations')[0]
        
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#', 
                           names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                                 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                                 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                                 'BiGG_Reaction', 'PFAMs'])
            
            for _, row in df.iterrows():
                gene_id = row['query']
                og_groups = str(row['eggNOG_OGs'])
                
                if og_groups != 'nan' and og_groups != '-':
                    primary_og = og_groups.split(',')[0].strip().split('@')[0]
                else:
                    primary_og = f"UNIQUE_{mag_id}_{gene_id}"
                
                all_annotations[primary_og] = {
                    'description': str(row['Description']),
                    'cog_category': str(row['COG_category']),
                    'kegg_ko': str(row['KEGG_ko']),
                    'kegg_pathway': str(row['KEGG_Pathway']),
                    'pfams': str(row['PFAMs']),
                    'preferred_name': str(row['Preferred_name'])
                }
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return all_annotations

def analyze_accessory_functions(pangenome_df, annotations):
    """Analyze functional content of accessory genome components"""
    
    # Separate by pan-genome category
    categories = ['shell', 'unique']  # Focus on truly accessory genes
    
    functional_analysis = {}
    
    for category in categories:
        category_genes = pangenome_df[pangenome_df['category'] == category]
        
        cog_counts = Counter()
        kegg_counts = Counter()
        pfam_counts = Counter()
        descriptions = []
        pathway_counts = Counter()
        
        for _, row in category_genes.iterrows():
            og_id = row['ortholog_group']
            
            if og_id in annotations:
                annotation = annotations[og_id]
                
                # COG categories
                cog = annotation['cog_category']
                if cog != 'nan' and cog != '-':
                    for cog_cat in cog:
                        if cog_cat.isalpha():
                            cog_counts[cog_cat] += 1
                
                # KEGG KOs
                kegg_ko = annotation['kegg_ko']
                if kegg_ko != 'nan' and kegg_ko != '-':
                    for ko in kegg_ko.split(','):
                        ko = ko.strip()
                        if ko.startswith('ko:'):
                            kegg_counts[ko] += 1
                
                # KEGG Pathways
                kegg_pathway = annotation['kegg_pathway']
                if kegg_pathway != 'nan' and kegg_pathway != '-':
                    pathways = kegg_pathway.split(',')
                    for pathway in pathways:
                        pathway = pathway.strip()
                        if pathway.startswith('map'):
                            pathway_counts[pathway] += 1
                
                # Pfam domains
                pfams = annotation['pfams']
                if pfams != 'nan' and pfams != '-':
                    for pfam in pfams.split(','):
                        pfam = pfam.strip()
                        if pfam:
                            pfam_counts[pfam] += 1
                
                # Descriptions
                desc = annotation['description']
                if desc != 'nan' and desc != '-' and len(desc) > 20:
                    descriptions.append(desc)
        
        functional_analysis[category] = {
            'total_genes': len(category_genes),
            'cog_categories': dict(cog_counts.most_common(15)),
            'kegg_kos': dict(kegg_counts.most_common(15)),
            'kegg_pathways': dict(pathway_counts.most_common(15)),
            'pfam_domains': dict(pfam_counts.most_common(15)),
            'example_functions': descriptions[:20]
        }
    
    return functional_analysis

def identify_strain_specific_functions(pangenome_df, annotations):
    """Identify strain-specific functional capabilities"""
    
    # Focus on unique genes
    unique_genes = pangenome_df[pangenome_df['category'] == 'unique']
    
    strain_functions = defaultdict(lambda: {
        'functions': [],
        'cog_categories': Counter(),
        'pathways': Counter()
    })
    
    for _, row in unique_genes.iterrows():
        og_id = row['ortholog_group']
        mags_present = row['present_in_mags'].split(',')
        
        if len(mags_present) == 1:  # Truly unique to one strain
            mag_id = mags_present[0]
            
            if og_id in annotations:
                annotation = annotations[og_id]
                
                # Store function description
                desc = annotation['description']
                if desc != 'nan' and desc != '-':
                    strain_functions[mag_id]['functions'].append(desc)
                
                # COG categories
                cog = annotation['cog_category']
                if cog != 'nan' and cog != '-':
                    for cog_cat in cog:
                        if cog_cat.isalpha():
                            strain_functions[mag_id]['cog_categories'][cog_cat] += 1
                
                # KEGG pathways
                kegg_pathway = annotation['kegg_pathway']
                if kegg_pathway != 'nan' and kegg_pathway != '-':
                    pathways = kegg_pathway.split(',')
                    for pathway in pathways:
                        pathway = pathway.strip()
                        if pathway.startswith('map'):
                            strain_functions[mag_id]['pathways'][pathway] += 1
    
    return strain_functions

def analyze_functional_specialization(pangenome_df, annotations):
    """Analyze functional specialization patterns"""
    
    # Identify genes present in few strains (2-4 strains)
    shell_genes = pangenome_df[pangenome_df['category'] == 'shell']
    
    specialization_patterns = {}
    
    # Group by presence count
    for presence_count in range(2, 5):  # 2-4 strains
        subset = shell_genes[shell_genes['presence_count'] == presence_count]
        
        if len(subset) > 0:
            cog_counts = Counter()
            pathway_counts = Counter()
            functions = []
            
            for _, row in subset.iterrows():
                og_id = row['ortholog_group']
                
                if og_id in annotations:
                    annotation = annotations[og_id]
                    
                    # COG categories
                    cog = annotation['cog_category']
                    if cog != 'nan' and cog != '-':
                        for cog_cat in cog:
                            if cog_cat.isalpha():
                                cog_counts[cog_cat] += 1
                    
                    # KEGG pathways
                    kegg_pathway = annotation['kegg_pathway']
                    if kegg_pathway != 'nan' and kegg_pathway != '-':
                        pathways = kegg_pathway.split(',')
                        for pathway in pathways:
                            pathway = pathway.strip()
                            if pathway.startswith('map'):
                                pathway_counts[pathway] += 1
                    
                    # Functions
                    desc = annotation['description']
                    if desc != 'nan' and desc != '-':
                        functions.append(desc)
            
            specialization_patterns[f"{presence_count}_strains"] = {
                'gene_count': len(subset),
                'cog_categories': dict(cog_counts.most_common(10)),
                'pathways': dict(pathway_counts.most_common(10)),
                'example_functions': functions[:10]
            }
    
    return specialization_patterns

def create_accessory_visualizations(functional_analysis, strain_functions, specialization_patterns, output_dir):
    """Create visualizations for accessory genome analysis"""
    
    # 1. COG category distribution in accessory genome
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Shell genes COG distribution
    if 'shell' in functional_analysis:
        shell_cogs = functional_analysis['shell']['cog_categories']
        if shell_cogs:
            cog_categories = list(shell_cogs.keys())[:10]
            cog_counts = [shell_cogs[cog] for cog in cog_categories]
            
            axes[0, 0].bar(cog_categories, cog_counts, color='lightblue', alpha=0.7)
            axes[0, 0].set_title('COG Categories in Shell Genes')
            axes[0, 0].set_ylabel('Number of Genes')
            axes[0, 0].tick_params(axis='x', rotation=45)
    
    # Unique genes COG distribution
    if 'unique' in functional_analysis:
        unique_cogs = functional_analysis['unique']['cog_categories']
        if unique_cogs:
            cog_categories = list(unique_cogs.keys())[:10]
            cog_counts = [unique_cogs[cog] for cog in cog_categories]
            
            axes[0, 1].bar(cog_categories, cog_counts, color='lightcoral', alpha=0.7)
            axes[0, 1].set_title('COG Categories in Unique Genes')
            axes[0, 1].set_ylabel('Number of Genes')
            axes[0, 1].tick_params(axis='x', rotation=45)
    
    # Strain-specific gene counts
    strain_ids = list(strain_functions.keys())
    strain_unique_counts = [len(strain_functions[strain]['functions']) for strain in strain_ids]
    
    if strain_ids:
        # Truncate long strain names for display
        display_names = [strain.split('_')[1] if '_' in strain else strain[:8] for strain in strain_ids]
        
        axes[1, 0].bar(display_names, strain_unique_counts, color='lightgreen', alpha=0.7)
        axes[1, 0].set_title('Strain-Specific Gene Counts')
        axes[1, 0].set_ylabel('Number of Unique Genes')
        axes[1, 0].tick_params(axis='x', rotation=45)
    
    # Specialization by presence count
    if specialization_patterns:
        presence_counts = list(specialization_patterns.keys())
        gene_counts = [specialization_patterns[pc]['gene_count'] for pc in presence_counts]
        
        axes[1, 1].bar(presence_counts, gene_counts, color='lightyellow', alpha=0.7)
        axes[1, 1].set_title('Gene Distribution by Presence Count')
        axes[1, 1].set_ylabel('Number of Gene Families')
        axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / "accessory_genome_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    pangenome_dir = results_dir / "pangenome_analysis"
    output_dir = pangenome_dir
    
    print("Starting accessory genome functional analysis...")
    
    # Load pan-genome data
    print("\n1. Loading pan-genome classification...")
    pangenome_df = load_pangenome_data(pangenome_dir)
    if pangenome_df is None:
        return
    
    # Load functional annotations
    print("\n2. Loading functional annotations...")
    annotations = load_functional_annotations(results_dir)
    print(f"Loaded annotations for {len(annotations)} ortholog groups")
    
    # Analyze accessory genome functions
    print("\n3. Analyzing accessory genome functions...")
    functional_analysis = analyze_accessory_functions(pangenome_df, annotations)
    
    # Identify strain-specific functions
    print("\n4. Identifying strain-specific functions...")
    strain_functions = identify_strain_specific_functions(pangenome_df, annotations)
    
    # Analyze functional specialization
    print("\n5. Analyzing functional specialization patterns...")
    specialization_patterns = analyze_functional_specialization(pangenome_df, annotations)
    
    # Create visualizations
    print("\n6. Creating visualizations...")
    create_accessory_visualizations(functional_analysis, strain_functions, specialization_patterns, output_dir)
    
    # Print results
    print("\n=== ACCESSORY GENOME FUNCTIONAL ANALYSIS ===")
    
    for category, analysis in functional_analysis.items():
        print(f"\n{category.upper()} GENES ({analysis['total_genes']} gene families):")
        
        print(f"  Top COG categories:")
        for cog, count in list(analysis['cog_categories'].items())[:5]:
            print(f"    {cog}: {count} genes")
        
        print(f"  Top KEGG pathways:")
        for pathway, count in list(analysis['kegg_pathways'].items())[:3]:
            print(f"    {pathway}: {count} genes")
        
        print(f"  Example functions:")
        for func in analysis['example_functions'][:3]:
            print(f"    {func[:80]}{'...' if len(func) > 80 else ''}")
    
    print(f"\n=== STRAIN-SPECIFIC FUNCTIONS ===")
    for strain, data in strain_functions.items():
        if data['functions']:
            print(f"\n{strain}: {len(data['functions'])} unique functions")
            print(f"  Top COG categories: {dict(data['cog_categories'].most_common(3))}")
            print(f"  Example functions:")
            for func in data['functions'][:2]:
                print(f"    {func[:80]}{'...' if len(func) > 80 else ''}")
    
    print(f"\n=== FUNCTIONAL SPECIALIZATION PATTERNS ===")
    for pattern, data in specialization_patterns.items():
        print(f"\n{pattern}: {data['gene_count']} gene families")
        if data['cog_categories']:
            print(f"  Top COG categories: {dict(list(data['cog_categories'].items())[:3])}")
    
    # Save detailed results
    # Strain-specific functions summary
    strain_summary = []
    for strain, data in strain_functions.items():
        strain_summary.append({
            'strain': strain,
            'unique_functions': len(data['functions']),
            'top_cog_categories': dict(data['cog_categories'].most_common(5)),
            'top_pathways': dict(data['pathways'].most_common(5))
        })
    
    strain_df = pd.DataFrame(strain_summary)
    strain_df.to_csv(output_dir / "strain_specific_functions.csv", index=False)
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - accessory_genome_analysis.png: Accessory genome functional overview")
    print("  - strain_specific_functions.csv: Strain-specific functional capabilities")
    
    return functional_analysis, strain_functions, specialization_patterns

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
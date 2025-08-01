#!/usr/bin/env python3
"""
Pan-genome Analysis for Ignatzschineria MAGs
Identify core, accessory, and unique genes across all strains
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
import itertools
import warnings
warnings.filterwarnings('ignore')

def load_all_proteins(results_dir):
    """Load all protein sequences from all MAGs"""
    gene_prediction_dir = results_dir / "gene_prediction"
    protein_files = list(gene_prediction_dir.glob("*_proteins.faa"))
    
    all_proteins = {}
    mag_gene_counts = {}
    
    for file_path in protein_files:
        mag_id = file_path.name.split('_genomic_proteins.faa')[0]
        proteins = {}
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                gene_id = record.id
                sequence = str(record.seq)
                proteins[gene_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'mag_id': mag_id
                }
            
            all_proteins[mag_id] = proteins
            mag_gene_counts[mag_id] = len(proteins)
            print(f"Loaded {len(proteins)} proteins from {mag_id}")
            
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
    
    return all_proteins, mag_gene_counts

def load_ortholog_assignments(results_dir):
    """Load ortholog assignments from eggNOG annotations"""
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    gene_to_og = {}  # gene_id -> ortholog_group
    og_to_genes = defaultdict(dict)  # og_id -> {mag_id: [gene_ids]}
    gene_annotations = {}  # gene_id -> annotation_data
    
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
                
                # Store annotation data
                gene_annotations[f"{mag_id}:{gene_id}"] = {
                    'description': str(row['Description']),
                    'cog_category': str(row['COG_category']),
                    'kegg_ko': str(row['KEGG_ko']),
                    'kegg_pathway': str(row['KEGG_Pathway'])
                }
                
                if og_groups != 'nan' and og_groups != '-':
                    # Use the first (most specific) ortholog group
                    ogs = [og.strip() for og in og_groups.split(',')]
                    primary_og = ogs[0].split('@')[0] if '@' in ogs[0] else ogs[0]
                    
                    gene_to_og[f"{mag_id}:{gene_id}"] = primary_og
                    
                    if mag_id not in og_to_genes[primary_og]:
                        og_to_genes[primary_og][mag_id] = []
                    og_to_genes[primary_og][mag_id].append(gene_id)
                else:
                    # Genes without ortholog assignment - create unique ID
                    unique_og = f"UNIQUE_{mag_id}_{gene_id}"
                    gene_to_og[f"{mag_id}:{gene_id}"] = unique_og
                    og_to_genes[unique_og][mag_id] = [gene_id]
                    
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return gene_to_og, og_to_genes, gene_annotations

def classify_pangenome_categories(og_to_genes, mag_ids):
    """Classify ortholog groups into core, accessory, and unique categories"""
    total_mags = len(mag_ids)
    
    core_genes = []        # Present in ≥95% of MAGs
    soft_core_genes = []   # Present in 80-94% of MAGs  
    shell_genes = []       # Present in 15-79% of MAGs
    cloud_genes = []       # Present in <15% of MAGs
    unique_genes = []      # Present in only 1 MAG
    
    category_thresholds = {
        'core': 0.95,
        'soft_core': 0.80,
        'shell': 0.15
    }
    
    for og_id, mag_genes in og_to_genes.items():
        presence_count = len(mag_genes)
        presence_fraction = presence_count / total_mags
        
        if presence_fraction >= category_thresholds['core']:
            core_genes.append((og_id, presence_count, mag_genes))
        elif presence_fraction >= category_thresholds['soft_core']:
            soft_core_genes.append((og_id, presence_count, mag_genes))
        elif presence_fraction >= category_thresholds['shell']:
            shell_genes.append((og_id, presence_count, mag_genes))
        elif presence_count == 1:
            unique_genes.append((og_id, presence_count, mag_genes))
        else:
            cloud_genes.append((og_id, presence_count, mag_genes))
    
    pangenome_categories = {
        'core': core_genes,
        'soft_core': soft_core_genes,
        'shell': shell_genes,
        'cloud': cloud_genes,
        'unique': unique_genes
    }
    
    return pangenome_categories

def calculate_pangenome_statistics(pangenome_categories, mag_ids):
    """Calculate pan-genome size and statistics"""
    stats = {}
    
    # Count genes in each category
    for category, genes in pangenome_categories.items():
        stats[f"{category}_count"] = len(genes)
        
        # Calculate total gene instances
        total_instances = sum(len(mag_genes) for _, _, mag_genes in genes)
        stats[f"{category}_instances"] = total_instances
    
    # Total pan-genome size
    stats['pangenome_size'] = sum(stats[f"{cat}_count"] for cat in pangenome_categories.keys())
    
    # Core genome size (strict + soft core)
    stats['core_genome_size'] = stats['core_count'] + stats['soft_core_count']
    
    # Accessory genome size
    stats['accessory_genome_size'] = stats['shell_count'] + stats['cloud_count'] + stats['unique_count']
    
    # Average genome size
    total_gene_instances = sum(stats[f"{cat}_instances"] for cat in pangenome_categories.keys())
    stats['average_genome_size'] = total_gene_instances / len(mag_ids)
    
    return stats

def analyze_pangenome_accumulation(pangenome_categories, mag_ids):
    """Analyze pan-genome accumulation curve"""
    # Create presence/absence matrix
    all_ogs = []
    for category_genes in pangenome_categories.values():
        all_ogs.extend([og_id for og_id, _, _ in category_genes])
    
    presence_matrix = pd.DataFrame(index=all_ogs, columns=mag_ids, dtype=int)
    presence_matrix.fillna(0, inplace=True)
    
    for category_genes in pangenome_categories.values():
        for og_id, _, mag_genes in category_genes:
            for mag_id in mag_genes.keys():
                presence_matrix.loc[og_id, mag_id] = 1
    
    # Calculate accumulation curves
    accumulation_data = []
    
    # Generate random orders for robust estimation
    n_permutations = 100
    for _ in range(n_permutations):
        mag_order = np.random.permutation(mag_ids)
        
        core_accumulation = []
        pan_accumulation = []
        
        for i in range(1, len(mag_order) + 1):
            current_mags = mag_order[:i]
            
            # Pan-genome: genes present in at least one MAG
            pan_genes = (presence_matrix[current_mags].sum(axis=1) > 0).sum()
            pan_accumulation.append(pan_genes)
            
            # Core genome: genes present in all MAGs
            core_genes = (presence_matrix[current_mags].sum(axis=1) == i).sum()
            core_accumulation.append(core_genes)
        
        for i, (pan, core) in enumerate(zip(pan_accumulation, core_accumulation)):
            accumulation_data.append({
                'n_genomes': i + 1,
                'pan_size': pan,
                'core_size': core,
                'permutation': _
            })
    
    accumulation_df = pd.DataFrame(accumulation_data)
    
    # Calculate means and confidence intervals
    summary_stats = accumulation_df.groupby('n_genomes').agg({
        'pan_size': ['mean', 'std', 'min', 'max'],
        'core_size': ['mean', 'std', 'min', 'max']
    }).round(1)
    
    return accumulation_df, summary_stats, presence_matrix

def analyze_functional_categories(pangenome_categories, gene_annotations):
    """Analyze functional categories in different pan-genome components"""
    functional_analysis = {}
    
    for category, genes in pangenome_categories.items():
        cog_counts = Counter()
        kegg_counts = Counter()
        descriptions = []
        
        for og_id, _, mag_genes in genes:
            # Get annotations for genes in this ortholog group
            for mag_id, gene_list in mag_genes.items():
                for gene_id in gene_list:
                    gene_key = f"{mag_id}:{gene_id}"
                    if gene_key in gene_annotations:
                        annotation = gene_annotations[gene_key]
                        
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
                        
                        # Descriptions
                        desc = annotation['description']
                        if desc != 'nan' and desc != '-' and len(desc) > 10:
                            descriptions.append(desc)
        
        functional_analysis[category] = {
            'cog_categories': dict(cog_counts.most_common(10)),
            'kegg_kos': dict(kegg_counts.most_common(10)),
            'example_functions': descriptions[:10]
        }
    
    return functional_analysis

def create_pangenome_visualizations(pangenome_categories, stats, accumulation_df, output_dir):
    """Create pan-genome visualizations"""
    
    # 1. Pan-genome size distribution pie chart
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Pie chart of gene categories
    categories = ['core', 'soft_core', 'shell', 'cloud', 'unique']
    sizes = [stats[f"{cat}_count"] for cat in categories]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    ax1.pie(sizes, labels=categories, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title('Pan-genome Gene Categories')
    
    # Bar chart of gene counts
    ax2.bar(categories, sizes, color=colors, alpha=0.7)
    ax2.set_ylabel('Number of Gene Families')
    ax2.set_title('Gene Family Counts by Category')
    ax2.tick_params(axis='x', rotation=45)
    
    # Pan-genome accumulation curve
    summary_stats = accumulation_df.groupby('n_genomes').agg({
        'pan_size': ['mean', 'std'],
        'core_size': ['mean', 'std']
    })
    
    n_genomes = summary_stats.index
    pan_mean = summary_stats[('pan_size', 'mean')]
    pan_std = summary_stats[('pan_size', 'std')]
    core_mean = summary_stats[('core_size', 'mean')]
    core_std = summary_stats[('core_size', 'std')]
    
    ax3.fill_between(n_genomes, pan_mean - pan_std, pan_mean + pan_std, alpha=0.3, color='blue')
    ax3.plot(n_genomes, pan_mean, 'b-', linewidth=2, label='Pan-genome')
    
    ax3.fill_between(n_genomes, core_mean - core_std, core_mean + core_std, alpha=0.3, color='red')
    ax3.plot(n_genomes, core_mean, 'r-', linewidth=2, label='Core genome')
    
    ax3.set_xlabel('Number of Genomes')
    ax3.set_ylabel('Number of Gene Families')
    ax3.set_title('Pan-genome Accumulation Curve')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Genome size distribution
    genome_sizes = []
    mag_names = []
    
    # Calculate individual genome sizes from presence matrix would go here
    # For now, show category contributions
    categories_subset = ['core', 'soft_core', 'shell', 'cloud', 'unique']
    bottoms = np.zeros(len(categories_subset))
    
    for i, cat in enumerate(categories_subset):
        count = stats[f"{cat}_count"]
        ax4.bar(cat, count, bottom=bottoms[i], color=colors[i], alpha=0.7)
    
    ax4.set_ylabel('Gene Families')
    ax4.set_title('Pan-genome Composition')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / "pangenome_overview.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    output_dir = results_dir / "pangenome_analysis"
    output_dir.mkdir(exist_ok=True)
    
    print("Starting pan-genome analysis...")
    
    # Step 1: Load all protein sequences
    print("\n1. Loading protein sequences...")
    all_proteins, mag_gene_counts = load_all_proteins(results_dir)
    mag_ids = list(all_proteins.keys())
    
    print(f"Loaded proteins from {len(mag_ids)} MAGs:")
    for mag_id, count in mag_gene_counts.items():
        print(f"  {mag_id}: {count} genes")
    
    # Step 2: Load ortholog assignments
    print("\n2. Loading ortholog assignments...")
    gene_to_og, og_to_genes, gene_annotations = load_ortholog_assignments(results_dir)
    
    print(f"Total ortholog groups: {len(og_to_genes)}")
    print(f"Genes with annotations: {len(gene_annotations)}")
    
    # Step 3: Classify pan-genome categories
    print("\n3. Classifying pan-genome categories...")
    pangenome_categories = classify_pangenome_categories(og_to_genes, mag_ids)
    
    # Step 4: Calculate statistics
    print("\n4. Calculating pan-genome statistics...")
    stats = calculate_pangenome_statistics(pangenome_categories, mag_ids)
    
    # Step 5: Analyze accumulation
    print("\n5. Analyzing pan-genome accumulation...")
    accumulation_df, summary_stats, presence_matrix = analyze_pangenome_accumulation(pangenome_categories, mag_ids)
    
    # Step 6: Functional analysis
    print("\n6. Analyzing functional categories...")
    functional_analysis = analyze_functional_categories(pangenome_categories, gene_annotations)
    
    # Print results
    print("\n=== PAN-GENOME ANALYSIS RESULTS ===")
    
    print(f"\nGenome Statistics:")
    print(f"  Total MAGs: {len(mag_ids)}")
    print(f"  Average genome size: {stats['average_genome_size']:.0f} genes")
    
    print(f"\nPan-genome Size: {stats['pangenome_size']:,} gene families")
    print(f"  Core genome: {stats['core_count']:,} families ({stats['core_count']/stats['pangenome_size']*100:.1f}%)")
    print(f"  Soft-core: {stats['soft_core_count']:,} families ({stats['soft_core_count']/stats['pangenome_size']*100:.1f}%)")
    print(f"  Shell: {stats['shell_count']:,} families ({stats['shell_count']/stats['pangenome_size']*100:.1f}%)")
    print(f"  Cloud: {stats['cloud_count']:,} families ({stats['cloud_count']/stats['pangenome_size']*100:.1f}%)")
    print(f"  Unique: {stats['unique_count']:,} families ({stats['unique_count']/stats['pangenome_size']*100:.1f}%)")
    
    print(f"\nCore vs Accessory:")
    print(f"  Core genome (≥80%): {stats['core_genome_size']:,} families")
    print(f"  Accessory genome: {stats['accessory_genome_size']:,} families")
    print(f"  Core/Pan ratio: {stats['core_genome_size']/stats['pangenome_size']:.2f}")
    
    # Pan-genome openness assessment
    final_pan_size = summary_stats.loc[len(mag_ids), ('pan_size', 'mean')]
    final_core_size = summary_stats.loc[len(mag_ids), ('core_size', 'mean')]
    
    print(f"\nPan-genome Properties:")
    print(f"  Pan-genome appears {'open' if stats['unique_count'] > stats['core_count']*0.1 else 'closed'}")
    print(f"  Unique genes per MAG: {stats['unique_count']/len(mag_ids):.1f} on average")
    
    # Create visualizations
    print("\n7. Creating visualizations...")
    create_pangenome_visualizations(pangenome_categories, stats, accumulation_df, output_dir)
    
    # Save detailed results
    # Pan-genome categories
    pangenome_summary = []
    for category, genes in pangenome_categories.items():
        for og_id, presence_count, mag_genes in genes:
            pangenome_summary.append({
                'ortholog_group': og_id,
                'category': category,
                'presence_count': presence_count,
                'present_in_mags': ','.join(mag_genes.keys())
            })
    
    pangenome_df = pd.DataFrame(pangenome_summary)
    pangenome_df.to_csv(output_dir / "pangenome_categories.csv", index=False)
    
    # Presence/absence matrix
    presence_matrix.to_csv(output_dir / "presence_absence_matrix.csv")
    
    # Accumulation statistics
    accumulation_summary = summary_stats.round(1)
    accumulation_summary.to_csv(output_dir / "accumulation_statistics.csv")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - pangenome_overview.png: Overview visualizations")
    print("  - pangenome_categories.csv: Detailed gene family classifications")
    print("  - presence_absence_matrix.csv: Gene presence/absence across MAGs")
    print("  - accumulation_statistics.csv: Pan-genome accumulation data")
    
    return stats, pangenome_categories, functional_analysis

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
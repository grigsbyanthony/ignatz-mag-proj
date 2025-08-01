#!/usr/bin/env python3
"""
Gene Gain/Loss Analysis for Ignatzschineria MAGs
Analyze evolutionary dynamics of gene families across phylogenetic tree
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import warnings
warnings.filterwarnings('ignore')

def load_phylogenetic_tree(phylo_dir):
    """Load phylogenetic distance matrix and create tree structure"""
    distance_file = phylo_dir / "distance_matrix.csv"
    if not distance_file.exists():
        print(f"Distance matrix not found: {distance_file}")
        return None, None
    
    distance_df = pd.read_csv(distance_file, index_col=0)
    
    # Convert to condensed distance matrix for clustering
    distance_matrix = distance_df.values
    condensed_distances = squareform(distance_matrix)
    
    # Create hierarchical clustering
    linkage_matrix = linkage(condensed_distances, method='average')
    
    # Convert to tree object for easier traversal
    tree = to_tree(linkage_matrix, rd=False)
    
    return tree, distance_df.index.tolist()

def load_presence_absence_matrix(pangenome_dir):
    """Load gene presence/absence matrix"""
    pa_file = pangenome_dir / "presence_absence_matrix.csv"
    if not pa_file.exists():
        print(f"Presence/absence matrix not found: {pa_file}")
        return None
    
    pa_matrix = pd.read_csv(pa_file, index_col=0)
    return pa_matrix

def load_functional_annotations(results_dir):
    """Load functional annotations for gene families"""
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

def create_tree_node_mapping(tree, mag_ids):
    """Create mapping between tree nodes and MAG IDs"""
    # Get leaf nodes
    leaves = tree.get_leaves()
    
    # Map leaf indices to MAG IDs
    leaf_to_mag = {}
    for i, leaf in enumerate(leaves):
        if i < len(mag_ids):
            leaf_to_mag[leaf.id] = mag_ids[i]
    
    # Create node mapping for all nodes
    node_mapping = {}
    node_counter = [len(mag_ids)]  # Start numbering internal nodes after leaves
    
    def assign_node_ids(node):
        if node.is_leaf():
            # Leaf node - use original ID
            node_mapping[node.id] = leaf_to_mag.get(node.id, f"Unknown_{node.id}")
        else:
            # Internal node - assign new ID
            node_mapping[node.id] = f"Node_{node_counter[0]}"
            node_counter[0] += 1
        
        # Recursively assign IDs to children
        if hasattr(node, 'left') and node.left:
            assign_node_ids(node.left)
        if hasattr(node, 'right') and node.right:
            assign_node_ids(node.right)
    
    assign_node_ids(tree)
    return node_mapping, leaf_to_mag

def reconstruct_ancestral_states(tree, pa_matrix, mag_ids):
    """Reconstruct ancestral gene presence states using parsimony"""
    
    # Create mapping from MAG IDs to leaf indices
    mag_to_index = {mag: i for i, mag in enumerate(mag_ids)}
    
    # Initialize ancestral states dictionary
    ancestral_states = {}
    
    # Process each gene family
    for gene_family in pa_matrix.index:
        gene_states = {}
        
        # Set leaf states
        for mag_id in mag_ids:
            if mag_id in pa_matrix.columns:
                gene_states[mag_id] = pa_matrix.loc[gene_family, mag_id]
            else:
                gene_states[mag_id] = 0
        
        # Reconstruct internal node states using simple parsimony
        def reconstruct_node(node, node_mapping):
            if node.is_leaf():
                # Leaf node - use observed state
                mag_id = node_mapping.get(node.id, f"Unknown_{node.id}")
                return gene_states.get(mag_id, 0)
            else:
                # Internal node - use parsimony
                left_state = reconstruct_node(node.left, node_mapping) if node.left else 0
                right_state = reconstruct_node(node.right, node_mapping) if node.right else 0
                
                # Simple parsimony: if both children have same state, use that
                # If different, choose presence (1) to minimize losses
                if left_state == right_state:
                    return left_state
                else:
                    return 1  # Favor presence to minimize evolutionary events
        
        # Create node mapping
        node_mapping, _ = create_tree_node_mapping(tree, mag_ids)
        
        # Reconstruct all node states
        def store_states(node):
            node_id = node_mapping.get(node.id, f"Node_{node.id}")
            state = reconstruct_node(node, node_mapping)
            gene_states[node_id] = state
            
            if hasattr(node, 'left') and node.left:
                store_states(node.left)
            if hasattr(node, 'right') and node.right:
                store_states(node.right)
        
        store_states(tree)
        ancestral_states[gene_family] = gene_states
    
    return ancestral_states

def identify_gain_loss_events(tree, ancestral_states, mag_ids):
    """Identify gene gain and loss events on each branch"""
    
    node_mapping, _ = create_tree_node_mapping(tree, mag_ids)
    
    gain_events = defaultdict(list)  # node -> [gained_genes]
    loss_events = defaultdict(list)  # node -> [lost_genes]
    
    def analyze_branch(node):
        if node.is_leaf():
            return
        
        # Get node IDs
        node_id = node_mapping.get(node.id, f"Node_{node.id}")
        
        # Analyze left child
        if node.left:
            left_id = node_mapping.get(node.left.id, f"Node_{node.left.id}")
            
            for gene_family, states in ancestral_states.items():
                parent_state = states.get(node_id, 0)
                child_state = states.get(left_id, 0)
                
                if parent_state == 0 and child_state == 1:
                    gain_events[left_id].append(gene_family)
                elif parent_state == 1 and child_state == 0:
                    loss_events[left_id].append(gene_family)
            
            analyze_branch(node.left)
        
        # Analyze right child
        if node.right:
            right_id = node_mapping.get(node.right.id, f"Node_{node.right.id}")
            
            for gene_family, states in ancestral_states.items():
                parent_state = states.get(node_id, 0)
                child_state = states.get(right_id, 0)
                
                if parent_state == 0 and child_state == 1:
                    gain_events[right_id].append(gene_family)
                elif parent_state == 1 and child_state == 0:
                    loss_events[right_id].append(gene_family)
            
            analyze_branch(node.right)
    
    analyze_branch(tree)
    
    return gain_events, loss_events

def analyze_functional_dynamics(gain_events, loss_events, annotations):
    """Analyze functional categories of gained and lost genes"""
    
    functional_analysis = {}
    
    # Analyze gains
    all_gained_genes = []
    for node, genes in gain_events.items():
        all_gained_genes.extend(genes)
    
    gained_cogs = Counter()
    gained_pathways = Counter()
    gained_functions = []
    
    for gene in all_gained_genes:
        if gene in annotations:
            annotation = annotations[gene]
            
            # COG categories
            cog = annotation['cog_category']
            if cog != 'nan' and cog != '-':
                for cog_cat in cog:
                    if cog_cat.isalpha():
                        gained_cogs[cog_cat] += 1
            
            # KEGG pathways
            pathway = annotation['kegg_pathway']
            if pathway != 'nan' and pathway != '-':
                for p in pathway.split(','):
                    p = p.strip()
                    if p.startswith('map'):
                        gained_pathways[p] += 1
            
            # Functions
            desc = annotation['description']
            if desc != 'nan' and desc != '-':
                gained_functions.append(desc)
    
    # Analyze losses
    all_lost_genes = []
    for node, genes in loss_events.items():
        all_lost_genes.extend(genes)
    
    lost_cogs = Counter()
    lost_pathways = Counter()
    lost_functions = []
    
    for gene in all_lost_genes:
        if gene in annotations:
            annotation = annotations[gene]
            
            # COG categories
            cog = annotation['cog_category']
            if cog != 'nan' and cog != '-':
                for cog_cat in cog:
                    if cog_cat.isalpha():
                        lost_cogs[cog_cat] += 1
            
            # KEGG pathways
            pathway = annotation['kegg_pathway']
            if pathway != 'nan' and pathway != '-':
                for p in pathway.split(','):
                    p = p.strip()
                    if p.startswith('map'):
                        lost_pathways[p] += 1
            
            # Functions
            desc = annotation['description']
            if desc != 'nan' and desc != '-':
                lost_functions.append(desc)
    
    functional_analysis = {
        'gains': {
            'total_events': len(all_gained_genes),
            'cog_categories': dict(gained_cogs.most_common(10)),
            'kegg_pathways': dict(gained_pathways.most_common(10)),
            'example_functions': gained_functions[:10]
        },
        'losses': {
            'total_events': len(all_lost_genes),
            'cog_categories': dict(lost_cogs.most_common(10)),
            'kegg_pathways': dict(lost_pathways.most_common(10)),
            'example_functions': lost_functions[:10]
        }
    }
    
    return functional_analysis

def analyze_strain_specific_evolution(gain_events, loss_events, mag_ids):
    """Analyze evolutionary patterns for each strain"""
    
    strain_evolution = {}
    
    for mag_id in mag_ids:
        gains = gain_events.get(mag_id, [])
        losses = loss_events.get(mag_id, [])
        
        strain_evolution[mag_id] = {
            'gains': len(gains),
            'losses': len(losses),
            'net_change': len(gains) - len(losses),
            'gained_genes': gains,
            'lost_genes': losses
        }
    
    return strain_evolution

def create_evolutionary_visualizations(strain_evolution, functional_analysis, output_dir):
    """Create visualizations for evolutionary dynamics"""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Gains vs Losses by strain
    strains = list(strain_evolution.keys())
    gains = [strain_evolution[s]['gains'] for s in strains]
    losses = [strain_evolution[s]['losses'] for s in strains]
    
    # Truncate strain names for display
    display_names = [s.split('_')[1] if '_' in s else s[:8] for s in strains]
    
    x = np.arange(len(strains))
    width = 0.35
    
    axes[0, 0].bar(x - width/2, gains, width, label='Gains', color='green', alpha=0.7)
    axes[0, 0].bar(x + width/2, losses, width, label='Losses', color='red', alpha=0.7)
    axes[0, 0].set_xlabel('Strain')
    axes[0, 0].set_ylabel('Number of Gene Events')
    axes[0, 0].set_title('Gene Gains vs Losses by Strain')
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(display_names, rotation=45, ha='right')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Net evolutionary change
    net_changes = [strain_evolution[s]['net_change'] for s in strains]
    colors = ['green' if x > 0 else 'red' if x < 0 else 'gray' for x in net_changes]
    
    axes[0, 1].bar(display_names, net_changes, color=colors, alpha=0.7)
    axes[0, 1].set_xlabel('Strain')
    axes[0, 1].set_ylabel('Net Gene Change (Gains - Losses)')
    axes[0, 1].set_title('Net Evolutionary Change by Strain')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].axhline(y=0, color='black', linestyle='-', alpha=0.5)
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Functional categories of gained genes
    if functional_analysis['gains']['cog_categories']:
        gained_cogs = functional_analysis['gains']['cog_categories']
        cog_cats = list(gained_cogs.keys())[:8]
        cog_gains = [gained_cogs[cog] for cog in cog_cats]
        
        axes[1, 0].bar(cog_cats, cog_gains, color='lightgreen', alpha=0.7)
        axes[1, 0].set_xlabel('COG Category')
        axes[1, 0].set_ylabel('Number of Gained Genes')
        axes[1, 0].set_title('Functional Categories of Gained Genes')
        axes[1, 0].tick_params(axis='x', rotation=45)
    
    # 4. Functional categories of lost genes
    if functional_analysis['losses']['cog_categories']:
        lost_cogs = functional_analysis['losses']['cog_categories']
        cog_cats = list(lost_cogs.keys())[:8]
        cog_losses = [lost_cogs[cog] for cog in cog_cats]
        
        axes[1, 1].bar(cog_cats, cog_losses, color='lightcoral', alpha=0.7)
        axes[1, 1].set_xlabel('COG Category')
        axes[1, 1].set_ylabel('Number of Lost Genes')
        axes[1, 1].set_title('Functional Categories of Lost Genes')
        axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / "evolutionary_dynamics.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    phylo_dir = results_dir / "phylogenetic_analysis"
    pangenome_dir = results_dir / "pangenome_analysis"
    output_dir = results_dir / "evolutionary_dynamics"
    output_dir.mkdir(exist_ok=True)
    
    print("Starting gene gain/loss analysis...")
    
    # Step 1: Load phylogenetic tree
    print("\n1. Loading phylogenetic tree...")
    tree, mag_ids = load_phylogenetic_tree(phylo_dir)
    if tree is None:
        print("Could not load phylogenetic tree!")
        return
    
    print(f"Loaded tree with {len(mag_ids)} taxa")
    
    # Step 2: Load presence/absence matrix
    print("\n2. Loading gene presence/absence matrix...")
    pa_matrix = load_presence_absence_matrix(pangenome_dir)
    if pa_matrix is None:
        print("Could not load presence/absence matrix!")
        return
    
    print(f"Loaded presence/absence data for {len(pa_matrix)} gene families")
    
    # Step 3: Load functional annotations
    print("\n3. Loading functional annotations...")
    annotations = load_functional_annotations(results_dir)
    print(f"Loaded annotations for {len(annotations)} gene families")
    
    # Step 4: Reconstruct ancestral states
    print("\n4. Reconstructing ancestral gene content...")
    ancestral_states = reconstruct_ancestral_states(tree, pa_matrix, mag_ids)
    print(f"Reconstructed ancestral states for {len(ancestral_states)} gene families")
    
    # Step 5: Identify gain/loss events
    print("\n5. Identifying gene gain and loss events...")
    gain_events, loss_events = identify_gain_loss_events(tree, ancestral_states, mag_ids)
    
    total_gains = sum(len(genes) for genes in gain_events.values())
    total_losses = sum(len(genes) for genes in loss_events.values())
    
    print(f"Identified {total_gains} gene gain events and {total_losses} gene loss events")
    
    # Step 6: Analyze functional categories
    print("\n6. Analyzing functional categories of evolutionary events...")
    functional_analysis = analyze_functional_dynamics(gain_events, loss_events, annotations)
    
    # Step 7: Analyze strain-specific evolution
    print("\n7. Analyzing strain-specific evolutionary patterns...")
    strain_evolution = analyze_strain_specific_evolution(gain_events, loss_events, mag_ids)
    
    # Step 8: Create visualizations
    print("\n8. Creating evolutionary dynamics visualizations...")
    create_evolutionary_visualizations(strain_evolution, functional_analysis, output_dir)
    
    # Print results
    print("\n=== EVOLUTIONARY DYNAMICS ANALYSIS ===")
    
    print(f"\nOverall Statistics:")
    print(f"  Total gene gain events: {total_gains}")
    print(f"  Total gene loss events: {total_losses}")
    print(f"  Net evolutionary change: {total_gains - total_losses}")
    print(f"  Gain/Loss ratio: {total_gains/max(total_losses, 1):.2f}")
    
    print(f"\nStrain-specific Evolution:")
    sorted_strains = sorted(strain_evolution.items(), key=lambda x: x[1]['net_change'], reverse=True)
    
    for strain, data in sorted_strains:
        print(f"  {strain}:")
        print(f"    Gains: {data['gains']}, Losses: {data['losses']}")
        print(f"    Net change: {data['net_change']:+d}")
    
    print(f"\nFunctional Analysis of Evolutionary Events:")
    
    print(f"\n  GAINED GENES ({functional_analysis['gains']['total_events']} events):")
    print(f"    Top COG categories: {dict(list(functional_analysis['gains']['cog_categories'].items())[:5])}")
    if functional_analysis['gains']['example_functions']:
        print(f"    Example functions:")
        for func in functional_analysis['gains']['example_functions'][:3]:
            print(f"      {func[:70]}{'...' if len(func) > 70 else ''}")
    
    print(f"\n  LOST GENES ({functional_analysis['losses']['total_events']} events):")
    print(f"    Top COG categories: {dict(list(functional_analysis['losses']['cog_categories'].items())[:5])}")
    if functional_analysis['losses']['example_functions']:
        print(f"    Example functions:")
        for func in functional_analysis['losses']['example_functions'][:3]:
            print(f"      {func[:70]}{'...' if len(func) > 70 else ''}")
    
    # Save detailed results
    # Strain evolution summary
    strain_summary = []
    for strain, data in strain_evolution.items():
        strain_summary.append({
            'strain': strain,
            'gains': data['gains'],
            'losses': data['losses'],
            'net_change': data['net_change']
        })
    
    strain_df = pd.DataFrame(strain_summary)
    strain_df.to_csv(output_dir / "strain_evolutionary_dynamics.csv", index=False)
    
    # Gain/loss events by gene family
    gain_loss_summary = []
    all_genes = set()
    all_genes.update(*[genes for genes in gain_events.values()])
    all_genes.update(*[genes for genes in loss_events.values()])
    
    for gene in all_genes:
        gain_count = sum(1 for genes in gain_events.values() if gene in genes)
        loss_count = sum(1 for genes in loss_events.values() if gene in genes)
        
        gain_loss_summary.append({
            'gene_family': gene,
            'total_gains': gain_count,
            'total_losses': loss_count,
            'net_evolution': gain_count - loss_count,
            'annotation': annotations.get(gene, {}).get('description', 'Unknown')
        })
    
    gain_loss_df = pd.DataFrame(gain_loss_summary)
    gain_loss_df.to_csv(output_dir / "gene_family_evolution.csv", index=False)
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - evolutionary_dynamics.png: Overview of evolutionary patterns")
    print("  - strain_evolutionary_dynamics.csv: Per-strain gain/loss summary")
    print("  - gene_family_evolution.csv: Per-gene family evolutionary events")
    
    return strain_evolution, functional_analysis, gain_events, loss_events

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
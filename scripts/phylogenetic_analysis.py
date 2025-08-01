#!/usr/bin/env python3
"""
Core Genome Phylogenetic Analysis for Ignatzschineria MAGs
Identify core genes, create alignments, and build phylogenetic relationships
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import matplotlib.pyplot as plt
import subprocess
import os

def parse_protein_files(results_dir):
    """Parse protein sequence files for all MAGs"""
    gene_prediction_dir = results_dir / "gene_prediction"
    protein_files = list(gene_prediction_dir.glob("*_proteins.faa"))
    
    mag_proteins = {}
    for file_path in protein_files:
        mag_id = file_path.name.split('_genomic_proteins.faa')[0]
        proteins = {}
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                # Extract gene ID from header
                gene_id = record.id
                proteins[gene_id] = str(record.seq)
            
            mag_proteins[mag_id] = proteins
            print(f"Loaded {len(proteins)} proteins from {mag_id}")
        except Exception as e:
            print(f"Error parsing {file_path}: {e}")
    
    return mag_proteins

def identify_orthologs_from_eggnog(results_dir):
    """Identify orthologous groups from eggNOG annotations"""
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    ortholog_groups = defaultdict(dict)  # og_id -> {mag_id: [gene_ids]}
    gene_to_og = {}  # gene_id -> og_id
    
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
                    # Parse ortholog groups (can be multiple)
                    ogs = [og.strip() for og in og_groups.split(',')]
                    for og in ogs:
                        if '@' in og:  # Valid OG format
                            og_id = og.split('@')[0]
                            if mag_id not in ortholog_groups[og_id]:
                                ortholog_groups[og_id][mag_id] = []
                            ortholog_groups[og_id][mag_id].append(gene_id)
                            gene_to_og[f"{mag_id}:{gene_id}"] = og_id
            
            print(f"Processed annotations for {mag_id}: {len(df)} genes")
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return ortholog_groups, gene_to_og

def identify_core_genes(ortholog_groups, mag_ids, min_presence=0.9):
    """Identify core genes present in most MAGs"""
    core_genes = []
    total_mags = len(mag_ids)
    min_mags = int(total_mags * min_presence)
    
    for og_id, mag_genes in ortholog_groups.items():
        present_mags = len(mag_genes)
        if present_mags >= min_mags:
            # Check if single copy in each MAG (avoid paralogs)
            single_copy = all(len(genes) == 1 for genes in mag_genes.values())
            if single_copy:
                core_genes.append((og_id, present_mags, mag_genes))
    
    # Sort by presence (most universal first)
    core_genes.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\nCore gene analysis:")
    print(f"Total ortholog groups: {len(ortholog_groups)}")
    print(f"Core genes (≥{min_presence*100}% MAGs, single-copy): {len(core_genes)}")
    
    return core_genes

def extract_core_gene_sequences(core_genes, mag_proteins, output_dir):
    """Extract sequences for core genes and prepare for alignment"""
    core_seq_dir = output_dir / "core_gene_sequences"
    core_seq_dir.mkdir(exist_ok=True)
    
    valid_core_genes = []
    
    for i, (og_id, presence, mag_genes) in enumerate(core_genes[:100]):  # Limit to top 100
        sequences = {}
        
        # Extract sequences for this core gene
        for mag_id, gene_list in mag_genes.items():
            if mag_id in mag_proteins:
                gene_id = gene_list[0]  # Single copy gene
                if gene_id in mag_proteins[mag_id]:
                    seq = mag_proteins[mag_id][gene_id]
                    sequences[mag_id] = seq
        
        # Only keep if we have sequences for most MAGs
        if len(sequences) >= len(mag_genes):
            valid_core_genes.append((og_id, sequences))
            
            # Write sequences to file
            output_file = core_seq_dir / f"{og_id}.faa"
            with open(output_file, 'w') as f:
                for mag_id, seq in sequences.items():
                    f.write(f">{mag_id}\n{seq}\n")
    
    print(f"Extracted sequences for {len(valid_core_genes)} core genes")
    return valid_core_genes

def create_core_genome_alignment(valid_core_genes, output_dir, max_genes=50):
    """Create concatenated alignment of core genes"""
    alignment_dir = output_dir / "alignments"
    alignment_dir.mkdir(exist_ok=True)
    
    # Limit to manageable number of genes
    selected_genes = valid_core_genes[:max_genes]
    
    # Get MAG IDs from first gene
    mag_ids = list(selected_genes[0][1].keys())
    
    # Initialize concatenated sequences
    concatenated_seqs = {mag_id: "" for mag_id in mag_ids}
    gene_info = []
    
    print(f"Creating alignment from {len(selected_genes)} core genes...")
    
    for og_id, sequences in selected_genes:
        if len(sequences) < len(mag_ids):
            continue  # Skip if missing sequences
        
        # Simple alignment by padding to same length (basic approach)
        max_len = max(len(seq) for seq in sequences.values())
        
        for mag_id in mag_ids:
            if mag_id in sequences:
                seq = sequences[mag_id]
                # Pad sequence to max length
                padded_seq = seq + '-' * (max_len - len(seq))
                concatenated_seqs[mag_id] += padded_seq
            else:
                # Fill with gaps if missing
                concatenated_seqs[mag_id] += '-' * max_len
        
        gene_info.append((og_id, max_len))
    
    # Write concatenated alignment
    concat_file = alignment_dir / "core_genome_alignment.fasta"
    with open(concat_file, 'w') as f:
        for mag_id, seq in concatenated_seqs.items():
            f.write(f">{mag_id}\n{seq}\n")
    
    # Write gene information
    info_file = alignment_dir / "gene_positions.txt"
    with open(info_file, 'w') as f:
        f.write("Gene_ID\tLength\tStart_Position\tEnd_Position\n")
        pos = 0
        for og_id, length in gene_info:
            f.write(f"{og_id}\t{length}\t{pos}\t{pos + length - 1}\n")
            pos += length
    
    print(f"Created concatenated alignment: {len(concatenated_seqs[mag_ids[0]])} positions")
    print(f"Alignment file: {concat_file}")
    
    return concat_file, len(concatenated_seqs[mag_ids[0]])

def calculate_distance_matrix(alignment_file):
    """Calculate pairwise distances between sequences"""
    sequences = {}
    
    # Read alignment
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences[record.id] = str(record.seq)
    
    mag_ids = list(sequences.keys())
    n_seqs = len(mag_ids)
    distance_matrix = np.zeros((n_seqs, n_seqs))
    
    # Calculate pairwise distances
    for i, mag1 in enumerate(mag_ids):
        for j, mag2 in enumerate(mag_ids):
            if i <= j:
                seq1 = sequences[mag1]
                seq2 = sequences[mag2]
                
                # Calculate Hamming distance (proportion of different positions)
                differences = sum(c1 != c2 for c1, c2 in zip(seq1, seq2) if c1 != '-' and c2 != '-')
                valid_positions = sum(1 for c1, c2 in zip(seq1, seq2) if c1 != '-' and c2 != '-')
                
                if valid_positions > 0:
                    distance = differences / valid_positions
                else:
                    distance = 1.0
                
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
    
    return distance_matrix, mag_ids

def create_simple_tree(distance_matrix, mag_ids, output_dir):
    """Create simple phylogenetic tree using hierarchical clustering"""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    
    # Convert to condensed distance matrix
    condensed_distances = squareform(distance_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distances, method='average')
    
    # Create dendrogram
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, labels=mag_ids, orientation='left', leaf_font_size=10)
    plt.title('Ignatzschineria Core Genome Phylogeny')
    plt.xlabel('Evolutionary Distance')
    plt.tight_layout()
    
    tree_file = output_dir / "core_genome_tree.png"
    plt.savefig(tree_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return linkage_matrix, tree_file

def analyze_phylogenetic_groups(linkage_matrix, mag_ids, distance_matrix):
    """Analyze phylogenetic groupings and relationships"""
    from scipy.cluster.hierarchy import fcluster
    
    # Create clusters at different distance thresholds
    clusters_02 = fcluster(linkage_matrix, 0.2, criterion='distance')
    clusters_01 = fcluster(linkage_matrix, 0.1, criterion='distance')
    clusters_005 = fcluster(linkage_matrix, 0.05, criterion='distance')
    
    print("\n=== PHYLOGENETIC ANALYSIS RESULTS ===")
    
    # Distance statistics
    distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    print(f"\nEvolutionary distances:")
    print(f"  Mean distance: {np.mean(distances):.4f}")
    print(f"  Min distance: {np.min(distances):.4f}")
    print(f"  Max distance: {np.max(distances):.4f}")
    print(f"  Std deviation: {np.std(distances):.4f}")
    
    # Clustering results
    print(f"\nPhylogenetic clusters:")
    print(f"  At 20% divergence: {max(clusters_02)} clusters")
    print(f"  At 10% divergence: {max(clusters_01)} clusters")
    print(f"  At 5% divergence: {max(clusters_005)} clusters")
    
    # Show closest and most distant pairs
    min_idx = np.unravel_index(np.argmin(distance_matrix + np.eye(len(mag_ids))), distance_matrix.shape)
    max_idx = np.unravel_index(np.argmax(distance_matrix), distance_matrix.shape)
    
    print(f"\nClosest pair:")
    print(f"  {mag_ids[min_idx[0]]} <-> {mag_ids[min_idx[1]]}: {distance_matrix[min_idx]:.4f}")
    
    print(f"\nMost distant pair:")
    print(f"  {mag_ids[max_idx[0]]} <-> {mag_ids[max_idx[1]]}: {distance_matrix[max_idx]:.4f}")
    
    return {
        'clusters_02': dict(zip(mag_ids, clusters_02)),
        'clusters_01': dict(zip(mag_ids, clusters_01)),
        'clusters_005': dict(zip(mag_ids, clusters_005)),
        'distances': distances
    }

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    output_dir = results_dir / "phylogenetic_analysis"
    output_dir.mkdir(exist_ok=True)
    
    print("Starting core genome phylogenetic analysis...")
    
    # Step 1: Parse protein sequences
    print("\n1. Loading protein sequences...")
    mag_proteins = parse_protein_files(results_dir)
    mag_ids = list(mag_proteins.keys())
    print(f"Loaded proteins from {len(mag_ids)} MAGs")
    
    # Step 2: Identify orthologous groups
    print("\n2. Identifying orthologous groups...")
    ortholog_groups, gene_to_og = identify_orthologs_from_eggnog(results_dir)
    
    # Step 3: Identify core genes
    print("\n3. Identifying core genes...")
    core_genes = identify_core_genes(ortholog_groups, mag_ids, min_presence=0.8)
    
    if len(core_genes) < 10:
        print("WARNING: Very few core genes found. Lowering stringency...")
        core_genes = identify_core_genes(ortholog_groups, mag_ids, min_presence=0.6)
    
    # Step 4: Extract core gene sequences
    print("\n4. Extracting core gene sequences...")
    valid_core_genes = extract_core_gene_sequences(core_genes, mag_proteins, output_dir)
    
    if len(valid_core_genes) == 0:
        print("ERROR: No valid core genes found with sequences!")
        return
    
    # Step 5: Create alignment
    print("\n5. Creating core genome alignment...")
    alignment_file, alignment_length = create_core_genome_alignment(valid_core_genes, output_dir)
    
    # Step 6: Calculate distances
    print("\n6. Calculating evolutionary distances...")
    distance_matrix, ordered_mag_ids = calculate_distance_matrix(alignment_file)
    
    # Step 7: Build tree
    print("\n7. Building phylogenetic tree...")
    linkage_matrix, tree_file = create_simple_tree(distance_matrix, ordered_mag_ids, output_dir)
    
    # Step 8: Analyze results
    print("\n8. Analyzing phylogenetic relationships...")
    phylo_results = analyze_phylogenetic_groups(linkage_matrix, ordered_mag_ids, distance_matrix)
    
    # Save results
    results_summary = {
        'total_core_genes': len(valid_core_genes),
        'alignment_length': alignment_length,
        'mag_count': len(mag_ids),
        'mean_distance': np.mean(phylo_results['distances']),
        'phylogenetic_groups': phylo_results
    }
    
    # Save distance matrix
    distance_df = pd.DataFrame(distance_matrix, index=ordered_mag_ids, columns=ordered_mag_ids)
    distance_df.to_csv(output_dir / "distance_matrix.csv")
    
    # Save core gene list
    core_genes_df = pd.DataFrame([
        {'ortholog_group': og_id, 'presence': presence, 'mags_present': len(mag_genes)}
        for og_id, presence, mag_genes in core_genes[:len(valid_core_genes)]
    ])
    core_genes_df.to_csv(output_dir / "core_genes_list.csv", index=False)
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - core_genome_tree.png: Phylogenetic tree")
    print("  - distance_matrix.csv: Pairwise evolutionary distances")
    print("  - core_genes_list.csv: List of core genes used")
    print("  - core_genome_alignment.fasta: Concatenated alignment")
    
    return results_summary

if __name__ == "__main__":
    try:
        results = main()
    except ImportError as e:
        print(f"Missing required package: {e}")
        print("Please install: pip install scipy biopython")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
#!/usr/bin/env python3
"""
Terpene Production Analysis for Ignatzschineria MAGs
Focus on isoprenoid biosynthesis pathways and terpene synthases
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import re

def parse_annotations_for_terpenes(file_path):
    """Parse eggNOG annotations for terpene-related genes"""
    df = pd.read_csv(file_path, sep='\t', comment='#', 
                     names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                           'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                           'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                           'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                           'BiGG_Reaction', 'PFAMs'])
    return df

def classify_terpene_genes(df):
    """Classify genes into terpene biosynthesis categories"""
    terpene_categories = {
        'mep_pathway': [],           # MEP/DOXP pathway for IPP/DMAPP
        'prenyl_transferases': [],   # Prenyltransferases
        'terpene_synthases': [],     # Actual terpene synthases
        'quinone_biosynthesis': [],  # Ubiquinone/menaquinone (prenylated)
        'sterol_related': [],        # Sterol biosynthesis
        'carotenoid_related': []     # Carotenoid biosynthesis
    }
    
    for idx, row in df.iterrows():
        description = str(row['Description']).lower()
        kegg_ko = str(row['KEGG_ko']).lower()
        preferred_name = str(row['Preferred_name']).lower()
        pfams = str(row['PFAMs']).lower()
        kegg_pathway = str(row['KEGG_Pathway']).lower()
        
        # MEP/DOXP pathway (main bacterial isoprenoid pathway)
        if any(term in description for term in ['isopentenyl diphosphate', 'dimethylallyl', 'isoprenoid']) or \
           any(ko in kegg_ko for ko in ['k01770', 'k03527', 'k00919', 'k01662', 'k00099', 'k03526']):
            terpene_categories['mep_pathway'].append(row)
        
        # Prenyltransferases (attach isoprenoid chains)
        elif any(term in description for term in ['prenyltransferase', 'prenyl', 'polyprenyl', 'undecaprenyl']):
            terpene_categories['prenyl_transferases'].append(row)
        
        # Terpene synthases (create terpene skeletons)
        elif any(term in description for term in ['terpene synthase', 'sesquiterpene', 'monoterpene', 'diterpene']):
            terpene_categories['terpene_synthases'].append(row)
        
        # Quinone biosynthesis (uses prenyl chains)
        elif any(term in description for term in ['ubiquinone', 'menaquinone', 'octaprenyl', 'phylloquinone']):
            terpene_categories['quinone_biosynthesis'].append(row)
        
        # Sterol-related
        elif any(term in description for term in ['sterol', 'cholesterol', 'squalene']):
            terpene_categories['sterol_related'].append(row)
        
        # Carotenoid-related
        elif any(term in description for term in ['carotenoid', 'phytoene', 'lycopene', 'beta-carotene']):
            terpene_categories['carotenoid_related'].append(row)
    
    return terpene_categories

def analyze_mep_pathway_completeness(terpene_genes):
    """Analyze completeness of MEP pathway"""
    mep_enzymes = {
        'dxs': 'K01662',      # 1-deoxy-D-xylulose-5-phosphate synthase
        'dxr': 'K00099',      # 1-deoxy-D-xylulose-5-phosphate reductoisomerase
        'ispD': 'K00991',     # 2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase
        'ispE': 'K00919',     # 4-diphosphocytidyl-2-C-methyl-D-erythritol kinase
        'ispF': 'K01770',     # 2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase
        'ispG': 'K03526',     # (E)-4-hydroxy-3-methylbut-2-enyl-diphosphate synthase
        'ispH': 'K03527',     # 4-hydroxy-3-methylbut-2-enyl diphosphate reductase
        'ispA': 'K00795',     # farnesyl diphosphate synthase
        'ispB': 'K02523'      # octaprenyl diphosphate synthase
    }
    
    pathway_presence = {enzyme: False for enzyme in mep_enzymes}
    
    for gene in terpene_genes['mep_pathway']:
        kegg_ko = str(gene['KEGG_ko']).lower()
        for enzyme, ko_id in mep_enzymes.items():
            if ko_id.lower() in kegg_ko:
                pathway_presence[enzyme] = True
    
    return pathway_presence

def identify_specific_terpene_products(terpene_genes):
    """Identify specific terpene products that could be synthesized"""
    products = {
        'ubiquinone': False,
        'menaquinone': False,
        'hopenes': False,
        'carotenoids': False,
        'other_terpenes': False
    }
    
    all_genes = []
    for category, genes in terpene_genes.items():
        all_genes.extend(genes)
    
    for gene in all_genes:
        description = str(gene['Description']).lower()
        kegg_pathway = str(gene['KEGG_Pathway']).lower()
        
        if 'ubiquinone' in description or 'map00130' in kegg_pathway:
            products['ubiquinone'] = True
        if 'menaquinone' in description or 'map00130' in kegg_pathway:
            products['menaquinone'] = True
        if 'hopene' in description or 'squalene' in description:
            products['hopenes'] = True
        if 'carotenoid' in description or 'map00906' in kegg_pathway:
            products['carotenoids'] = True
        if any(term in description for term in ['terpene', 'sesquiterpene', 'monoterpene']):
            products['other_terpenes'] = True
    
    return products

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    output_dir = results_dir / "comparative_analysis"
    output_dir.mkdir(exist_ok=True)
    
    # Get all annotation files
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    mag_ids = [f.name.split('_genomic_eggnog.emapper.annotations')[0] for f in annotation_files]
    
    print(f"Analyzing terpene production potential for {len(mag_ids)} MAGs")
    
    # Initialize data structures
    all_terpene_genes = {}
    terpene_summary = defaultdict(dict)
    mep_pathway_completeness = {}
    terpene_products = {}
    
    # Process each MAG
    for file_path in annotation_files:
        mag_id = file_path.name.split('_genomic_eggnog.emapper.annotations')[0]
        print(f"Processing {mag_id}...")
        
        # Parse annotations
        df = parse_annotations_for_terpenes(file_path)
        
        # Classify terpene genes
        terpene_genes = classify_terpene_genes(df)
        all_terpene_genes[mag_id] = terpene_genes
        
        # Analyze MEP pathway completeness
        mep_completeness = analyze_mep_pathway_completeness(terpene_genes)
        mep_pathway_completeness[mag_id] = mep_completeness
        
        # Identify potential products
        products = identify_specific_terpene_products(terpene_genes)
        terpene_products[mag_id] = products
        
        # Count genes for summary
        for category, genes in terpene_genes.items():
            terpene_summary[category][mag_id] = len(genes)
    
    print("\n=== TERPENE PRODUCTION ANALYSIS ===")
    
    # Summary statistics
    print(f"\nTerpene Biosynthesis Gene Categories:")
    for category in terpene_summary:
        total_genes = sum(terpene_summary[category].values())
        avg_per_mag = np.mean(list(terpene_summary[category].values()))
        print(f"  {category.replace('_', ' ').title()}: {total_genes} total ({avg_per_mag:.1f} avg per MAG)")
    
    # MEP pathway analysis
    print(f"\n=== MEP PATHWAY COMPLETENESS ===")
    mep_enzymes = ['dxs', 'dxr', 'ispD', 'ispE', 'ispF', 'ispG', 'ispH', 'ispA', 'ispB']
    
    for mag_id in mag_ids:
        completeness = mep_pathway_completeness[mag_id]
        present_enzymes = sum(completeness.values())
        total_enzymes = len(completeness)
        percentage = (present_enzymes / total_enzymes) * 100
        
        print(f"\n{mag_id}: {present_enzymes}/{total_enzymes} enzymes ({percentage:.1f}%)")
        missing_enzymes = [enzyme for enzyme, present in completeness.items() if not present]
        if missing_enzymes:
            print(f"  Missing: {', '.join(missing_enzymes)}")
        else:
            print("  Complete MEP pathway!")
    
    # Terpene product analysis
    print(f"\n=== POTENTIAL TERPENE PRODUCTS ===")
    product_summary = defaultdict(int)
    
    for mag_id, products in terpene_products.items():
        print(f"\n{mag_id}:")
        for product, present in products.items():
            if present:
                product_summary[product] += 1
                print(f"  ✓ {product.replace('_', ' ').title()}")
    
    print(f"\nProduct Summary (MAGs with capacity):")
    for product, count in sorted(product_summary.items(), key=lambda x: x[1], reverse=True):
        print(f"  {product.replace('_', ' ').title()}: {count}/{len(mag_ids)} MAGs")
    
    # Key findings
    print(f"\n=== KEY FINDINGS ===")
    
    # Top performers
    total_terpene_genes = {mag: sum(len(genes) for genes in all_terpene_genes[mag].values()) 
                          for mag in mag_ids}
    
    print(f"\nMAGs with highest terpene gene content:")
    for mag_id, count in sorted(total_terpene_genes.items(), key=lambda x: x[1], reverse=True)[:3]:
        print(f"  {mag_id}: {count} terpene-related genes")
    
    # MEP pathway leaders
    complete_pathways = []
    for mag_id, completeness in mep_pathway_completeness.items():
        if all(completeness.values()):
            complete_pathways.append(mag_id)
    
    if complete_pathways:
        print(f"\nMAGs with complete MEP pathways: {len(complete_pathways)}")
        for mag_id in complete_pathways:
            print(f"  {mag_id}")
    else:
        print(f"\nNo MAGs with completely intact MEP pathways")
    
    # Create comparison matrices
    terpene_df = pd.DataFrame(terpene_summary).fillna(0).T
    
    # MEP pathway completeness matrix
    mep_df = pd.DataFrame({mag_id: {enzyme: 1 if present else 0 
                                   for enzyme, present in completeness.items()}
                          for mag_id, completeness in mep_pathway_completeness.items()}).T
    
    # Product capability matrix
    products_df = pd.DataFrame({mag_id: {product: 1 if present else 0 
                                        for product, present in products.items()}
                               for mag_id, products in terpene_products.items()}).T
    
    # Save results
    terpene_df.to_csv(output_dir / "terpene_genes_comparison.csv")
    mep_df.to_csv(output_dir / "mep_pathway_completeness.csv")
    products_df.to_csv(output_dir / "terpene_products_capability.csv")
    
    # Create visualizations
    plt.style.use('seaborn-v0_8')
    
    # Terpene gene categories heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(terpene_df, annot=True, fmt='.0f', cmap='Purples', 
                cbar_kws={'label': 'Number of Genes'})
    plt.title('Terpene Biosynthesis Genes Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('Gene Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "terpene_genes_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # MEP pathway completeness heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(mep_df, annot=True, fmt='.0f', cmap='RdYlGn', 
                cbar_kws={'label': 'Enzyme Present (1=Yes, 0=No)'})
    plt.title('MEP Pathway Completeness Across Ignatzschineria MAGs')
    plt.xlabel('MAG')
    plt.ylabel('MEP Pathway Enzyme')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "mep_pathway_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Terpene production potential summary
    plt.figure(figsize=(12, 8))
    
    # Calculate total terpene potential per MAG
    total_scores = []
    for mag_id in mag_ids:
        score = sum(len(genes) for genes in all_terpene_genes[mag_id].values())
        total_scores.append(score)
    
    # Create bar plot
    colors = plt.cm.viridis(np.linspace(0, 1, len(mag_ids)))
    bars = plt.bar(range(len(mag_ids)), total_scores, color=colors, alpha=0.8)
    
    plt.xlabel('MAG')
    plt.ylabel('Total Terpene-Related Genes')
    plt.title('Terpene Production Potential Across Ignatzschineria MAGs')
    plt.xticks(range(len(mag_ids)), mag_ids, rotation=45, ha='right')
    
    # Add value labels on bars
    for i, (bar, score) in enumerate(zip(bars, total_scores)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(total_scores)*0.01, 
                str(score), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(output_dir / "terpene_potential_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Detailed gene analysis
    print(f"\n=== SPECIFIC GENE ANALYSIS ===")
    
    # Most common terpene-related functions
    all_descriptions = []
    for mag_id in mag_ids:
        for category, genes in all_terpene_genes[mag_id].items():
            for gene in genes:
                desc = str(gene['Description'])
                if desc != 'nan' and desc != '-':
                    all_descriptions.append(desc)
    
    desc_counts = Counter(all_descriptions)
    print(f"\nMost common terpene-related functions:")
    for desc, count in desc_counts.most_common(8):
        print(f"  {desc[:70]}{'...' if len(desc) > 70 else ''}: {count}")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - terpene_genes_comparison.csv")
    print("  - mep_pathway_completeness.csv")
    print("  - terpene_products_capability.csv")
    print("  - terpene_genes_heatmap.png")
    print("  - mep_pathway_heatmap.png")
    print("  - terpene_potential_summary.png")

if __name__ == "__main__":
    main()
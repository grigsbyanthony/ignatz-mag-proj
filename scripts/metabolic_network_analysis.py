#!/usr/bin/env python3
"""
Metabolic Network Reconstruction for Ignatzschineria MAGs
Analyze pathway completeness and metabolic capabilities
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# KEGG pathway and module definitions for key metabolic processes
KEGG_PATHWAYS = {
    'map00010': 'Glycolysis / Gluconeogenesis',
    'map00020': 'Citrate cycle (TCA cycle)',
    'map00030': 'Pentose phosphate pathway',
    'map00040': 'Pentose and glucuronate interconversions',
    'map00051': 'Fructose and mannose metabolism',
    'map00052': 'Galactose metabolism',
    'map00053': 'Ascorbate and aldarate metabolism',
    'map00061': 'Fatty acid biosynthesis',
    'map00071': 'Fatty acid degradation',
    'map00190': 'Oxidative phosphorylation',
    'map00220': 'Arginine biosynthesis',
    'map00230': 'Purine metabolism',
    'map00240': 'Pyrimidine metabolism',
    'map00250': 'Alanine, aspartate and glutamate metabolism',
    'map00260': 'Glycine, serine and threonine metabolism',
    'map00270': 'Cysteine and methionine metabolism',
    'map00280': 'Valine, leucine and isoleucine degradation',
    'map00290': 'Valine, leucine and isoleucine biosynthesis',
    'map00300': 'Lysine biosynthesis',
    'map00310': 'Lysine degradation',
    'map00330': 'Arginine and proline metabolism',
    'map00340': 'Histidine metabolism',
    'map00350': 'Tyrosine metabolism',
    'map00360': 'Phenylalanine metabolism',
    'map00380': 'Tryptophan metabolism',
    'map00400': 'Phenylalanine, tyrosine and tryptophan biosynthesis',
    'map00450': 'Selenocompound metabolism',
    'map00460': 'Cyanoamino acid metabolism',
    'map00471': 'D-Glutamine and D-glutamate metabolism',
    'map00480': 'Glutathione metabolism',
    'map00500': 'Starch and sucrose metabolism',
    'map00520': 'Amino sugar and nucleotide sugar metabolism',
    'map00540': 'Lipopolysaccharide biosynthesis',
    'map00550': 'Peptidoglycan biosynthesis',
    'map00561': 'Glycerolipid metabolism',
    'map00562': 'Inositol phosphate metabolism',
    'map00564': 'Glycerophospholipid metabolism',
    'map00620': 'Pyruvate metabolism',
    'map00630': 'Glyoxylate and dicarboxylate metabolism',
    'map00640': 'Propanoate metabolism',
    'map00650': 'Butanoate metabolism',
    'map00660': 'C5-Branched dibasic acid metabolism',
    'map00670': 'One carbon pool by folate',
    'map00680': 'Methane metabolism',
    'map00710': 'Carbon fixation in photosynthetic organisms',
    'map00720': 'Carbon fixation pathways in prokaryotes',
    'map00730': 'Thiamine metabolism',
    'map00740': 'Riboflavin metabolism',
    'map00750': 'Vitamin B6 metabolism',
    'map00760': 'Nicotinate and nicotinamide metabolism',
    'map00770': 'Pantothenate and CoA biosynthesis',
    'map00780': 'Biotin metabolism',
    'map00785': 'Lipoic acid metabolism',
    'map00790': 'Folate biosynthesis',
    'map00860': 'Porphyrin and chlorophyll metabolism',
    'map00900': 'Terpenoid backbone biosynthesis',
    'map00910': 'Nitrogen metabolism',
    'map00920': 'Sulfur metabolism'
}

# Key metabolic modules with their essential genes (using simplified pathways)
METABOLIC_MODULES = {
    'Glycolysis': {'name': 'Glycolysis (core pathway)', 'genes': ['K00844', 'K01810', 'K00850', 'K01623', 'K01803', 'K00134']},
    'TCA_cycle': {'name': 'Citrate cycle (TCA cycle)', 'genes': ['K01647', 'K01681', 'K01902', 'K01903', 'K00164', 'K00658']},
    'PPP_oxidative': {'name': 'Pentose phosphate pathway (oxidative)', 'genes': ['K00036', 'K00033', 'K01057']},
    'PPP_non_oxidative': {'name': 'Pentose phosphate pathway (non-oxidative)', 'genes': ['K00615', 'K01783', 'K01807']},
    'Fatty_acid_biosynthesis': {'name': 'Fatty acid biosynthesis', 'genes': ['K00059', 'K00648', 'K02078', 'K00668']},
    'Amino_acid_biosynthesis': {'name': 'Amino acid biosynthesis (core)', 'genes': ['K01586', 'K00928', 'K01915', 'K00133']},
    'Purine_biosynthesis': {'name': 'Purine biosynthesis', 'genes': ['K00764', 'K01945', 'K01589', 'K00602']},
    'Pyrimidine_biosynthesis': {'name': 'Pyrimidine biosynthesis', 'genes': ['K01465', 'K00609', 'K00762', 'K00940']},
    'Oxidative_phosphorylation': {'name': 'Oxidative phosphorylation (Complex I)', 'genes': ['K00330', 'K00331', 'K00332', 'K00333']},
    'ATP_synthase': {'name': 'ATP synthase', 'genes': ['K02111', 'K02112', 'K02115', 'K02108']},
    'Riboflavin_biosynthesis': {'name': 'Riboflavin biosynthesis', 'genes': ['K00794', 'K00793', 'K00082', 'K00861']},
    'Folate_biosynthesis': {'name': 'Folate biosynthesis', 'genes': ['K01495', 'K00796', 'K00287', 'K00950']},
    'Biotin_biosynthesis': {'name': 'Biotin biosynthesis', 'genes': ['K00652', 'K00833', 'K01935', 'K00060']},
    'Cobalamin_biosynthesis': {'name': 'Cobalamin biosynthesis', 'genes': ['K02189', 'K02190', 'K02191', 'K00595']},
    'Nitrogen_metabolism': {'name': 'Nitrogen metabolism', 'genes': ['K00362', 'K00363', 'K01915', 'K01916']},
    'Sulfur_metabolism': {'name': 'Sulfur metabolism', 'genes': ['K00392', 'K00380', 'K00381', 'K01011']},
    'Cell_wall_biosynthesis': {'name': 'Peptidoglycan biosynthesis', 'genes': ['K01000', 'K02563', 'K01921', 'K05364']},
    'DNA_replication': {'name': 'DNA replication', 'genes': ['K02314', 'K02315', 'K03469', 'K02342']},
    'Transcription': {'name': 'RNA polymerase', 'genes': ['K03040', 'K03041', 'K03043', 'K03046']},
}

def load_kegg_annotations(results_dir):
    """Load KEGG annotations from eggNOG mapper results"""
    eggnog_dir = results_dir / "functional_annotation/eggnog"
    annotation_files = list(eggnog_dir.glob("*_eggnog.emapper.annotations"))
    
    mag_kegg_data = {}
    
    for file_path in annotation_files:
        mag_id = file_path.name.split('_genomic_eggnog.emapper.annotations')[0]
        
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#', 
                           names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                                 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                                 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                                 'BiGG_Reaction', 'PFAMs'])
            
            kegg_data = {
                'genes': [],
                'kos': set(),
                'pathways': set(),
                'modules': set(),
                'reactions': set(),
                'ec_numbers': set()
            }
            
            for _, row in df.iterrows():
                gene_id = row['query']
                
                # KEGG KO numbers
                kegg_ko = str(row['KEGG_ko'])
                if kegg_ko != 'nan' and kegg_ko != '-':
                    kos = [ko.strip() for ko in kegg_ko.split(',') if ko.strip().startswith('ko:')]
                    for ko in kos:
                        ko_id = ko.replace('ko:', '')
                        kegg_data['kos'].add(ko_id)
                        kegg_data['genes'].append({
                            'gene_id': gene_id,
                            'ko': ko_id,
                            'description': str(row['Description'])
                        })
                
                # KEGG Pathways
                kegg_pathway = str(row['KEGG_Pathway'])
                if kegg_pathway != 'nan' and kegg_pathway != '-':
                    pathways = [p.strip() for p in kegg_pathway.split(',') if p.strip().startswith('map')]
                    kegg_data['pathways'].update(pathways)
                
                # KEGG Modules
                kegg_module = str(row['KEGG_Module'])
                if kegg_module != 'nan' and kegg_module != '-':
                    modules = [m.strip() for m in kegg_module.split(',') if m.strip().startswith('M')]
                    kegg_data['modules'].update(modules)
                
                # KEGG Reactions
                kegg_reaction = str(row['KEGG_Reaction'])
                if kegg_reaction != 'nan' and kegg_reaction != '-':
                    reactions = [r.strip() for r in kegg_reaction.split(',') if r.strip().startswith('R')]
                    kegg_data['reactions'].update(reactions)
                
                # EC numbers
                ec = str(row['EC'])
                if ec != 'nan' and ec != '-':
                    ecs = [e.strip() for e in ec.split(',') if e.strip()]
                    kegg_data['ec_numbers'].update(ecs)
            
            mag_kegg_data[mag_id] = kegg_data
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return mag_kegg_data

def calculate_pathway_completeness(mag_kegg_data):
    """Calculate pathway completeness for each MAG"""
    pathway_completeness = {}
    
    for mag_id, kegg_data in mag_kegg_data.items():
        mag_kos = kegg_data['kos']
        
        pathway_scores = {}
        
        # Calculate completeness for each metabolic module
        for module_id, module_info in METABOLIC_MODULES.items():
            required_genes = set(module_info['genes'])
            present_genes = required_genes & mag_kos
            
            completeness = len(present_genes) / len(required_genes) if required_genes else 0
            
            pathway_scores[module_id] = {
                'name': module_info['name'],
                'completeness': completeness,
                'present_genes': list(present_genes),
                'missing_genes': list(required_genes - present_genes),
                'total_genes': len(required_genes)
            }
        
        pathway_completeness[mag_id] = pathway_scores
    
    return pathway_completeness

def analyze_metabolic_capabilities(mag_kegg_data, pathway_completeness):
    """Analyze metabolic capabilities and specializations"""
    
    metabolic_profiles = {}
    
    for mag_id, kegg_data in mag_kegg_data.items():
        
        # Basic metabolic statistics
        total_kos = len(kegg_data['kos'])
        total_pathways = len(kegg_data['pathways'])
        total_modules = len(kegg_data['modules'])
        
        # Analyze pathway completeness
        complete_pathways = []
        partial_pathways = []
        missing_pathways = []
        
        for module_id, module_data in pathway_completeness[mag_id].items():
            if module_data['completeness'] >= 0.8:
                complete_pathways.append((module_id, module_data['name'], module_data['completeness']))
            elif module_data['completeness'] >= 0.3:
                partial_pathways.append((module_id, module_data['name'], module_data['completeness']))
            else:
                missing_pathways.append((module_id, module_data['name'], module_data['completeness']))
        
        # Categorize metabolic capabilities
        metabolic_categories = {
            'energy_metabolism': [],
            'carbon_metabolism': [],
            'amino_acid_metabolism': [],
            'cofactor_biosynthesis': [],
            'nucleotide_metabolism': []
        }
        
        # Map pathways to pathways present in the MAG
        for pathway in kegg_data['pathways']:
            if pathway in KEGG_PATHWAYS:
                pathway_name = KEGG_PATHWAYS[pathway]
                
                if any(term in pathway_name.lower() for term in ['oxidative phosphorylation', 'citrate cycle', 'glycolysis']):
                    metabolic_categories['energy_metabolism'].append(pathway_name)
                elif any(term in pathway_name.lower() for term in ['carbon', 'pentose', 'pyruvate', 'propanoate']):
                    metabolic_categories['carbon_metabolism'].append(pathway_name)
                elif any(term in pathway_name.lower() for term in ['amino', 'alanine', 'glycine', 'valine', 'lysine', 'arginine']):
                    metabolic_categories['amino_acid_metabolism'].append(pathway_name)
                elif any(term in pathway_name.lower() for term in ['vitamin', 'thiamine', 'riboflavin', 'biotin', 'folate']):
                    metabolic_categories['cofactor_biosynthesis'].append(pathway_name)
                elif any(term in pathway_name.lower() for term in ['purine', 'pyrimidine', 'nucleotide']):
                    metabolic_categories['nucleotide_metabolism'].append(pathway_name)
        
        metabolic_profiles[mag_id] = {
            'total_kos': total_kos,
            'total_pathways': total_pathways,
            'total_modules': total_modules,
            'complete_pathways': complete_pathways,
            'partial_pathways': partial_pathways,
            'missing_pathways': missing_pathways,
            'metabolic_categories': metabolic_categories,
            'metabolic_versatility': len(complete_pathways) + 0.5 * len(partial_pathways)
        }
    
    return metabolic_profiles

def identify_metabolic_specializations(metabolic_profiles):
    """Identify strain-specific metabolic specializations"""
    
    # Find pathways that are complete in few strains
    pathway_distribution = defaultdict(list)
    
    for mag_id, profile in metabolic_profiles.items():
        for module_id, name, completeness in profile['complete_pathways']:
            pathway_distribution[module_id].append((mag_id, completeness))
    
    specializations = {}
    
    for module_id, strain_data in pathway_distribution.items():
        if len(strain_data) <= 3:  # Present in 3 or fewer strains
            module_name = METABOLIC_MODULES.get(module_id, {}).get('name', 'Unknown')
            specializations[module_id] = {
                'name': module_name,
                'strain_count': len(strain_data),
                'strains': strain_data
            }
    
    return specializations

def analyze_metabolic_complementarity(metabolic_profiles):
    """Analyze metabolic complementarity between strains"""
    
    mag_ids = list(metabolic_profiles.keys())
    complementarity_matrix = np.zeros((len(mag_ids), len(mag_ids)))
    
    for i, mag1 in enumerate(mag_ids):
        for j, mag2 in enumerate(mag_ids):
            if i != j:
                # Get complete pathways for each strain
                mag1_pathways = set([module_id for module_id, _, _ in metabolic_profiles[mag1]['complete_pathways']])
                mag2_pathways = set([module_id for module_id, _, _ in metabolic_profiles[mag2]['complete_pathways']])
                
                # Calculate complementarity: unique pathways in mag2 that mag1 lacks
                unique_in_mag2 = mag2_pathways - mag1_pathways
                complementarity = len(unique_in_mag2)
                
                complementarity_matrix[i, j] = complementarity
    
    return complementarity_matrix, mag_ids

def create_metabolic_visualizations(pathway_completeness, metabolic_profiles, specializations, 
                                  complementarity_matrix, mag_ids, output_dir):
    """Create visualizations for metabolic analysis"""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Pathway completeness heatmap
    # Create matrix of pathway completeness scores
    modules = list(METABOLIC_MODULES.keys())
    completeness_matrix = np.zeros((len(mag_ids), len(modules)))
    
    for i, mag_id in enumerate(mag_ids):
        for j, module_id in enumerate(modules):
            if module_id in pathway_completeness[mag_id]:
                completeness_matrix[i, j] = pathway_completeness[mag_id][module_id]['completeness']
    
    # Truncate MAG names for display
    display_mag_names = [mag.split('_')[1] if '_' in mag else mag[:8] for mag in mag_ids]
    module_names = [METABOLIC_MODULES[m]['name'][:25] + '...' if len(METABOLIC_MODULES[m]['name']) > 25 
                   else METABOLIC_MODULES[m]['name'] for m in modules]
    
    im1 = axes[0, 0].imshow(completeness_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
    axes[0, 0].set_xticks(range(len(modules)))
    axes[0, 0].set_xticklabels(module_names, rotation=45, ha='right', fontsize=8)
    axes[0, 0].set_yticks(range(len(mag_ids)))
    axes[0, 0].set_yticklabels(display_mag_names, fontsize=8)
    axes[0, 0].set_title('Metabolic Pathway Completeness')
    plt.colorbar(im1, ax=axes[0, 0], shrink=0.8)
    
    # 2. Metabolic versatility comparison
    versatility_scores = [profile['metabolic_versatility'] for profile in metabolic_profiles.values()]
    
    axes[0, 1].bar(display_mag_names, versatility_scores, color='lightblue', alpha=0.7)
    axes[0, 1].set_ylabel('Metabolic Versatility Score')
    axes[0, 1].set_title('Metabolic Versatility by Strain')
    axes[0, 1].tick_params(axis='x', rotation=45)
    
    # 3. Metabolic specializations
    if specializations:
        spec_names = [spec['name'][:20] + '...' if len(spec['name']) > 20 else spec['name'] 
                     for spec in specializations.values()]
        spec_counts = [spec['strain_count'] for spec in specializations.values()]
        
        axes[1, 0].bar(range(len(spec_names)), spec_counts, color='lightcoral', alpha=0.7)
        axes[1, 0].set_xticks(range(len(spec_names)))
        axes[1, 0].set_xticklabels(spec_names, rotation=45, ha='right', fontsize=8)
        axes[1, 0].set_ylabel('Number of Strains')
        axes[1, 0].set_title('Metabolic Specializations (Rare Pathways)')
    
    # 4. Metabolic complementarity heatmap
    im2 = axes[1, 1].imshow(complementarity_matrix, cmap='Blues', aspect='auto')
    axes[1, 1].set_xticks(range(len(mag_ids)))
    axes[1, 1].set_xticklabels(display_mag_names, rotation=45, ha='right', fontsize=8)
    axes[1, 1].set_yticks(range(len(mag_ids)))
    axes[1, 1].set_yticklabels(display_mag_names, fontsize=8)
    axes[1, 1].set_title('Metabolic Complementarity Matrix')
    plt.colorbar(im2, ax=axes[1, 1], shrink=0.8)
    
    plt.tight_layout()
    plt.savefig(output_dir / "metabolic_network_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    output_dir = results_dir / "metabolic_analysis"
    output_dir.mkdir(exist_ok=True)
    
    print("Starting metabolic network reconstruction...")
    
    # Step 1: Load KEGG annotations
    print("\n1. Loading KEGG annotations...")
    mag_kegg_data = load_kegg_annotations(results_dir)
    
    if not mag_kegg_data:
        print("No KEGG data loaded!")
        return
    
    print(f"Loaded KEGG data for {len(mag_kegg_data)} MAGs")
    for mag_id, kegg_data in mag_kegg_data.items():
        print(f"  {mag_id}: {len(kegg_data['kos'])} KOs, {len(kegg_data['pathways'])} pathways")
    
    # Step 2: Calculate pathway completeness
    print("\n2. Calculating pathway completeness...")
    pathway_completeness = calculate_pathway_completeness(mag_kegg_data)
    print(f"Calculated completeness for {len(METABOLIC_MODULES)} metabolic modules")
    
    # Step 3: Analyze metabolic capabilities
    print("\n3. Analyzing metabolic capabilities...")
    metabolic_profiles = analyze_metabolic_capabilities(mag_kegg_data, pathway_completeness)
    
    # Step 4: Identify metabolic specializations
    print("\n4. Identifying metabolic specializations...")
    specializations = identify_metabolic_specializations(metabolic_profiles)
    print(f"Identified {len(specializations)} metabolic specializations")
    
    # Step 5: Analyze metabolic complementarity
    print("\n5. Analyzing metabolic complementarity...")
    complementarity_matrix, mag_ids = analyze_metabolic_complementarity(metabolic_profiles)
    
    # Step 6: Create visualizations
    print("\n6. Creating metabolic network visualizations...")
    create_metabolic_visualizations(pathway_completeness, metabolic_profiles, specializations,
                                  complementarity_matrix, mag_ids, output_dir)
    
    # Print results
    print("\n=== METABOLIC NETWORK ANALYSIS RESULTS ===")
    
    print(f"\nOverall Metabolic Statistics:")
    total_kos = sum(len(data['kos']) for data in mag_kegg_data.values())
    avg_kos_per_strain = total_kos / len(mag_kegg_data)
    print(f"  Total unique KOs across all strains: {len(set().union(*[data['kos'] for data in mag_kegg_data.values()]))}")
    print(f"  Average KOs per strain: {avg_kos_per_strain:.1f}")
    
    print(f"\nMetabolic Versatility Ranking:")
    versatility_ranking = sorted(metabolic_profiles.items(), 
                               key=lambda x: x[1]['metabolic_versatility'], reverse=True)
    
    for i, (mag_id, profile) in enumerate(versatility_ranking[:5]):
        print(f"  {i+1}. {mag_id}: {profile['metabolic_versatility']:.1f} points")
        print(f"     Complete pathways: {len(profile['complete_pathways'])}")
        print(f"     Partial pathways: {len(profile['partial_pathways'])}")
    
    print(f"\nMetabolic Specializations (Rare Pathways):")
    for module_id, spec in specializations.items():
        print(f"  {spec['name']} ({module_id}): {spec['strain_count']} strains")
        for strain, completeness in spec['strains']:
            print(f"    {strain}: {completeness:.2f}")
    
    print(f"\nKey Metabolic Capabilities by Category:")
    category_summary = defaultdict(Counter)
    
    for mag_id, profile in metabolic_profiles.items():
        for category, pathways in profile['metabolic_categories'].items():
            category_summary[category][len(pathways)] += 1
    
    for category, counts in category_summary.items():
        avg_pathways = sum(k * v for k, v in counts.items()) / sum(counts.values())
        print(f"  {category.replace('_', ' ').title()}: {avg_pathways:.1f} pathways per strain on average")
    
    # Find most complementary strain pairs
    print(f"\nMost Complementary Strain Pairs:")
    max_complementarity = 0
    best_pairs = []
    
    for i, mag1 in enumerate(mag_ids):
        for j, mag2 in enumerate(mag_ids):
            if i < j:  # Avoid duplicates
                comp_score = complementarity_matrix[i, j] + complementarity_matrix[j, i]
                if comp_score > max_complementarity:
                    max_complementarity = comp_score
                    best_pairs = [(mag1, mag2, comp_score)]
                elif comp_score == max_complementarity:
                    best_pairs.append((mag1, mag2, comp_score))
    
    for mag1, mag2, score in best_pairs[:3]:
        print(f"  {mag1} <-> {mag2}: {score} complementary pathways")
    
    # Save detailed results
    # Pathway completeness matrix
    completeness_data = []
    for mag_id in mag_ids:
        for module_id, module_data in pathway_completeness[mag_id].items():
            completeness_data.append({
                'mag_id': mag_id,
                'module_id': module_id,
                'module_name': module_data['name'],
                'completeness': module_data['completeness'],
                'present_genes': ','.join(module_data['present_genes']),
                'missing_genes': ','.join(module_data['missing_genes'])
            })
    
    completeness_df = pd.DataFrame(completeness_data)
    completeness_df.to_csv(output_dir / "pathway_completeness.csv", index=False)
    
    # Metabolic profiles summary
    profile_data = []
    for mag_id, profile in metabolic_profiles.items():
        profile_data.append({
            'mag_id': mag_id,
            'total_kos': profile['total_kos'],
            'total_pathways': profile['total_pathways'],
            'complete_pathways': len(profile['complete_pathways']),
            'partial_pathways': len(profile['partial_pathways']),
            'metabolic_versatility': profile['metabolic_versatility']
        })
    
    profile_df = pd.DataFrame(profile_data)
    profile_df.to_csv(output_dir / "metabolic_profiles.csv", index=False)
    
    # Complementarity matrix
    comp_df = pd.DataFrame(complementarity_matrix, index=mag_ids, columns=mag_ids)
    comp_df.to_csv(output_dir / "metabolic_complementarity.csv")
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - metabolic_network_analysis.png: Overview of metabolic capabilities")
    print("  - pathway_completeness.csv: Detailed pathway completeness scores")
    print("  - metabolic_profiles.csv: Metabolic capability summary per strain")
    print("  - metabolic_complementarity.csv: Strain complementarity matrix")
    
    return pathway_completeness, metabolic_profiles, specializations

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
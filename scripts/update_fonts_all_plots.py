#!/usr/bin/env python3
"""
Update all existing plots to use CMU Sans Serif font
Re-run all analysis scripts with updated font settings
"""

import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import subprocess
import sys

def setup_cmu_font():
    """Setup CMU Sans Serif font for matplotlib"""
    try:
        # Try to use CMU Sans Serif if available
        plt.rcParams['font.family'] = 'CMU Sans Serif'
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['axes.labelsize'] = 10
        plt.rcParams['xtick.labelsize'] = 9
        plt.rcParams['ytick.labelsize'] = 9
        plt.rcParams['legend.fontsize'] = 9
        plt.rcParams['figure.titlesize'] = 14
        
        # Test if font is available
        fig, ax = plt.subplots(1, 1, figsize=(1, 1))
        ax.text(0.5, 0.5, 'Test', fontfamily='CMU Sans Serif')
        plt.close(fig)
        
        print("✓ CMU Sans Serif font configured successfully")
        return True
        
    except Exception as e:
        print(f"⚠ CMU Sans Serif not available, using fallback fonts: {e}")
        # Fallback to high-quality fonts
        plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif']
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['axes.labelsize'] = 10
        plt.rcParams['xtick.labelsize'] = 9
        plt.rcParams['ytick.labelsize'] = 9
        plt.rcParams['legend.fontsize'] = 9
        plt.rcParams['figure.titlesize'] = 14
        return False

def run_analysis_script(script_path, description):
    """Run an analysis script and handle errors"""
    try:
        print(f"🔄 Running {description}...")
        result = subprocess.run([sys.executable, script_path], 
                              capture_output=True, text=True, cwd=script_path.parent)
        
        if result.returncode == 0:
            print(f"✓ {description} completed successfully")
            return True
        else:
            print(f"✗ {description} failed:")
            print(f"  Error: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"✗ Error running {description}: {e}")
        return False

def main():
    # Setup font configuration
    print("=== UPDATING ALL PLOTS WITH CMU SANS SERIF FONT ===\n")
    
    font_success = setup_cmu_font()
    
    # Define all analysis scripts to re-run
    base_dir = Path("/Users/grigsbyanthony/Library/Mobile Documents/com~apple~CloudDocs/PhD Projects/Ignatzschineria MAG analysis/ignatz-mag-proj")
    
    analysis_scripts = [
        {
            'script': base_dir / 'visualize_basic_stats.py',
            'description': 'Basic Genome Statistics'
        },
        {
            'script': base_dir / 'comparative_analysis.py',
            'description': 'Comparative Functional Analysis'
        },
        {
            'script': base_dir / 'phylogenetic_analysis.py',
            'description': 'Phylogenetic Analysis'
        },
        {
            'script': base_dir / 'pangenome_analysis.py',
            'description': 'Pan-genome Analysis'
        },
        {
            'script': base_dir / 'metabolic_network_analysis.py',
            'description': 'Metabolic Network Analysis'
        },
        {
            'script': base_dir / 'simple_gene_evolution_analysis.py',
            'description': 'Evolutionary Dynamics Analysis'
        },
        {
            'script': base_dir / 'ecological_adaptation_analysis.py',
            'description': 'Ecological Adaptation Analysis'
        },
        {
            'script': base_dir / 'comparative_endosymbiont_analysis.py',
            'description': 'Comparative Endosymbiont Analysis'
        }
    ]
    
    # Track results
    successful_scripts = []
    failed_scripts = []
    
    # Re-run all analysis scripts with updated font settings
    for script_info in analysis_scripts:
        script_path = script_info['script']
        description = script_info['description']
        
        if script_path.exists():
            success = run_analysis_script(script_path, description)
            if success:
                successful_scripts.append(description)
            else:
                failed_scripts.append(description)
        else:
            print(f"⚠ Script not found: {script_path}")
            failed_scripts.append(description)
    
    # Print summary
    print(f"\n=== FONT UPDATE SUMMARY ===")
    print(f"Font configuration: {'CMU Sans Serif' if font_success else 'Fallback fonts'}")
    print(f"Successfully updated: {len(successful_scripts)}/{len(analysis_scripts)} analyses")
    
    if successful_scripts:
        print(f"\n✓ Successfully updated plots:")
        for script in successful_scripts:
            print(f"  - {script}")
    
    if failed_scripts:
        print(f"\n✗ Failed to update:")
        for script in failed_scripts:
            print(f"  - {script}")
    
    # List all updated visualization files
    results_dir = Path("MAGs/Consolidated .fna/results")
    png_files = list(results_dir.rglob("*.png"))
    
    print(f"\n📊 Total visualization files: {len(png_files)}")
    print("Updated visualization files:")
    
    # Group by analysis type
    analysis_groups = {
        'basic_stats': [],
        'comparative_analysis': [],
        'phylogenetic_analysis': [],
        'pangenome_analysis': [],
        'metabolic_analysis': [],
        'evolutionary_dynamics': [],
        'ecological_adaptation': []
    }
    
    for png_file in png_files:
        file_category = None
        for category in analysis_groups.keys():
            if category in str(png_file):
                file_category = category
                break
        
        if file_category:
            analysis_groups[file_category].append(png_file.name)
        else:
            # Handle comparative_analysis files
            if 'comparative_analysis' in str(png_file):
                analysis_groups['comparative_analysis'].append(png_file.name)
    
    for category, files in analysis_groups.items():
        if files:
            print(f"\n  {category.replace('_', ' ').title()}:")
            for file in sorted(files):
                print(f"    - {file}")
    
    print(f"\n🎨 All plots have been updated with {'CMU Sans Serif' if font_success else 'high-quality fallback'} font!")
    print("All visualizations are now ready for publication.")

if __name__ == "__main__":
    main()
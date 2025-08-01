#!/usr/bin/env python3
"""
Visualize Basic Statistics for Ignatzschineria MAGs
Extract and visualize genome assembly statistics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re

def load_basic_stats(results_dir):
    """Load basic statistics from detailed stats files"""
    basic_stats_dir = results_dir / "basic_stats"
    stats_files = list(basic_stats_dir.glob("*_detailed_stats.txt"))
    
    mag_stats = []
    
    for file_path in stats_files:
        # Extract MAG ID from filename
        mag_id = file_path.name.replace('_genomic_detailed_stats.txt', '')
        
        try:
            # Read the stats file
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Parse the data line (skip header)
            if len(lines) >= 2:
                data_line = lines[1].strip()
                # The file is actually one long line with spaces, need to parse carefully
                # Split by multiple spaces to get proper fields
                import re
                fields = re.split(r'\s+', data_line)
                
                if len(fields) >= 18:  # Ensure we have enough fields
                    # Extract key statistics based on the actual format
                    stats = {
                        'mag_id': mag_id,
                        'num_contigs': int(fields[3]),
                        'total_length': int(fields[4].replace(',', '')),
                        'min_contig_length': int(fields[5].replace(',', '')),
                        'avg_contig_length': float(fields[6].replace(',', '')),
                        'max_contig_length': int(fields[7].replace(',', '')),
                        'n50': int(fields[12].replace(',', '')),
                        'n50_count': int(fields[13]),
                        'gc_content': float(fields[-2])  # GC% is second to last field
                    }
                else:
                    print(f"Insufficient fields in {file_path}: {len(fields)} fields found")
                    print(f"Fields: {fields}")
                    continue
                
                mag_stats.append(stats)
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return pd.DataFrame(mag_stats)

def create_basic_stats_visualizations(stats_df, output_dir):
    """Create comprehensive visualizations of basic genome statistics"""
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Create display names
    display_names = [mag.replace('GCA_', '').replace('_ASM', '_ASM').split('_')[0] + '_' + 
                    mag.replace('GCA_', '').replace('_ASM', '_ASM').split('_')[1][:7] 
                    for mag in stats_df['mag_id']]
    
    # 1. Genome Size Distribution
    axes[0, 0].bar(range(len(stats_df)), stats_df['total_length'] / 1e6, 
                   color='lightblue', alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('MAG')
    axes[0, 0].set_ylabel('Genome Size (Mb)')
    axes[0, 0].set_title('Genome Size Distribution')
    axes[0, 0].set_xticks(range(len(stats_df)))
    axes[0, 0].set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, v in enumerate(stats_df['total_length'] / 1e6):
        axes[0, 0].text(i, v + 0.01, f'{v:.2f}', ha='center', va='bottom', fontsize=8)
    
    # 2. Number of Contigs
    axes[0, 1].bar(range(len(stats_df)), stats_df['num_contigs'], 
                   color='lightcoral', alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('MAG')
    axes[0, 1].set_ylabel('Number of Contigs')
    axes[0, 1].set_title('Assembly Fragmentation')
    axes[0, 1].set_xticks(range(len(stats_df)))
    axes[0, 1].set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(stats_df['num_contigs']):
        axes[0, 1].text(i, v + 2, str(v), ha='center', va='bottom', fontsize=8)
    
    # 3. N50 Values
    axes[0, 2].bar(range(len(stats_df)), stats_df['n50'] / 1000, 
                   color='lightgreen', alpha=0.7, edgecolor='black')
    axes[0, 2].set_xlabel('MAG')
    axes[0, 2].set_ylabel('N50 (kb)')
    axes[0, 2].set_title('Assembly Contiguity (N50)')
    axes[0, 2].set_xticks(range(len(stats_df)))
    axes[0, 2].set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
    axes[0, 2].grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(stats_df['n50'] / 1000):
        axes[0, 2].text(i, v + 0.2, f'{v:.1f}', ha='center', va='bottom', fontsize=8)
    
    # 4. GC Content
    axes[1, 0].bar(range(len(stats_df)), stats_df['gc_content'], 
                   color='lightyellow', alpha=0.7, edgecolor='black')
    axes[1, 0].set_xlabel('MAG')
    axes[1, 0].set_ylabel('GC Content (%)')
    axes[1, 0].set_title('GC Content Distribution')
    axes[1, 0].set_xticks(range(len(stats_df)))
    axes[1, 0].set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
    axes[1, 0].grid(True, alpha=0.3)
    
    # Add horizontal line for average
    avg_gc = stats_df['gc_content'].mean()
    axes[1, 0].axhline(y=avg_gc, color='red', linestyle='--', alpha=0.7, 
                      label=f'Average: {avg_gc:.1f}%')
    axes[1, 0].legend()
    
    # Add value labels
    for i, v in enumerate(stats_df['gc_content']):
        axes[1, 0].text(i, v + 0.3, f'{v:.1f}', ha='center', va='bottom', fontsize=8)
    
    # 5. Average Contig Length
    axes[1, 1].bar(range(len(stats_df)), stats_df['avg_contig_length'] / 1000, 
                   color='lightpink', alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('MAG')
    axes[1, 1].set_ylabel('Average Contig Length (kb)')
    axes[1, 1].set_title('Average Contig Size')
    axes[1, 1].set_xticks(range(len(stats_df)))
    axes[1, 1].set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(stats_df['avg_contig_length'] / 1000):
        axes[1, 1].text(i, v + 0.1, f'{v:.1f}', ha='center', va='bottom', fontsize=8)
    
    # 6. Assembly Quality Scatter Plot (N50 vs Number of Contigs)
    scatter = axes[1, 2].scatter(stats_df['num_contigs'], stats_df['n50'] / 1000, 
                                c=stats_df['total_length'] / 1e6, cmap='viridis', 
                                s=100, alpha=0.7, edgecolor='black')
    
    axes[1, 2].set_xlabel('Number of Contigs')
    axes[1, 2].set_ylabel('N50 (kb)')
    axes[1, 2].set_title('Assembly Quality Assessment')
    axes[1, 2].grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=axes[1, 2])
    cbar.set_label('Genome Size (Mb)')
    
    # Add strain labels
    for i, (x, y, name) in enumerate(zip(stats_df['num_contigs'], 
                                        stats_df['n50'] / 1000, 
                                        display_names)):
        axes[1, 2].annotate(name.split('_')[0], (x, y), 
                           xytext=(5, 5), textcoords='offset points', 
                           fontsize=7, alpha=0.8)
    
    plt.tight_layout()
    plt.savefig(output_dir / "basic_genome_statistics.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_table(stats_df, output_dir):
    """Create a summary table of basic statistics"""
    
    # Calculate summary statistics
    summary_stats = {
        'Metric': [
            'Total Length (Mb)',
            'Number of Contigs', 
            'N50 (kb)',
            'Max Contig (kb)',
            'Average Contig (kb)',
            'GC Content (%)'
        ],
        'Min': [
            f"{stats_df['total_length'].min() / 1e6:.2f}",
            f"{stats_df['num_contigs'].min()}",
            f"{stats_df['n50'].min() / 1000:.1f}",
            f"{stats_df['max_contig_length'].min() / 1000:.1f}",
            f"{stats_df['avg_contig_length'].min() / 1000:.1f}",
            f"{stats_df['gc_content'].min():.1f}"
        ],
        'Max': [
            f"{stats_df['total_length'].max() / 1e6:.2f}",
            f"{stats_df['num_contigs'].max()}",
            f"{stats_df['n50'].max() / 1000:.1f}",
            f"{stats_df['max_contig_length'].max() / 1000:.1f}",
            f"{stats_df['avg_contig_length'].max() / 1000:.1f}",
            f"{stats_df['gc_content'].max():.1f}"
        ],
        'Mean': [
            f"{stats_df['total_length'].mean() / 1e6:.2f}",
            f"{stats_df['num_contigs'].mean():.1f}",
            f"{stats_df['n50'].mean() / 1000:.1f}",
            f"{stats_df['max_contig_length'].mean() / 1000:.1f}",
            f"{stats_df['avg_contig_length'].mean() / 1000:.1f}",
            f"{stats_df['gc_content'].mean():.1f}"
        ],
        'Std Dev': [
            f"{stats_df['total_length'].std() / 1e6:.2f}",
            f"{stats_df['num_contigs'].std():.1f}",
            f"{stats_df['n50'].std() / 1000:.1f}",
            f"{stats_df['max_contig_length'].std() / 1000:.1f}",
            f"{stats_df['avg_contig_length'].std() / 1000:.1f}",
            f"{stats_df['gc_content'].std():.1f}"
        ]
    }
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Create table visualization
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('tight')
    ax.axis('off')
    
    table = ax.table(cellText=summary_df.values,
                    colLabels=summary_df.columns,
                    cellLoc='center',
                    loc='center',
                    bbox=[0, 0, 1, 1])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 2)
    
    # Style the table
    for i in range(len(summary_df.columns)):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    for i in range(1, len(summary_df) + 1):
        for j in range(len(summary_df.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')
    
    plt.title('Basic Genome Statistics Summary', fontsize=16, fontweight='bold', pad=20)
    plt.savefig(output_dir / "genome_statistics_summary_table.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save as CSV
    summary_df.to_csv(output_dir / "genome_statistics_summary.csv", index=False)
    
    return summary_df

def print_detailed_stats(stats_df):
    """Print detailed statistics to console"""
    
    print("=== IGNATZSCHINERIA MAG BASIC STATISTICS ===\n")
    
    print("Individual MAG Statistics:")
    print("-" * 80)
    
    for _, row in stats_df.iterrows():
        mag_short = row['mag_id'].replace('GCA_', '').split('_')[0]
        print(f"{mag_short}:")
        print(f"  Genome Size: {row['total_length'] / 1e6:.2f} Mb")
        print(f"  Contigs: {row['num_contigs']}")
        print(f"  N50: {row['n50'] / 1000:.1f} kb")
        print(f"  Max Contig: {row['max_contig_length'] / 1000:.1f} kb")
        print(f"  GC Content: {row['gc_content']:.1f}%")
        print()
    
    print("Summary Statistics:")
    print("-" * 40)
    print(f"Total MAGs analyzed: {len(stats_df)}")
    print(f"Genome size range: {stats_df['total_length'].min() / 1e6:.2f} - {stats_df['total_length'].max() / 1e6:.2f} Mb")
    print(f"Average genome size: {stats_df['total_length'].mean() / 1e6:.2f} ± {stats_df['total_length'].std() / 1e6:.2f} Mb")
    print(f"Contig count range: {stats_df['num_contigs'].min()} - {stats_df['num_contigs'].max()}")
    print(f"Average contig count: {stats_df['num_contigs'].mean():.1f} ± {stats_df['num_contigs'].std():.1f}")
    print(f"N50 range: {stats_df['n50'].min() / 1000:.1f} - {stats_df['n50'].max() / 1000:.1f} kb")
    print(f"Average N50: {stats_df['n50'].mean() / 1000:.1f} ± {stats_df['n50'].std() / 1000:.1f} kb")
    print(f"GC content range: {stats_df['gc_content'].min():.1f} - {stats_df['gc_content'].max():.1f}%")
    print(f"Average GC content: {stats_df['gc_content'].mean():.1f} ± {stats_df['gc_content'].std():.1f}%")
    
    # Assembly quality assessment
    print(f"\nAssembly Quality Assessment:")
    print("-" * 40)
    
    # High quality: N50 > 10kb and contigs < 200
    high_quality = stats_df[(stats_df['n50'] > 10000) & (stats_df['num_contigs'] < 200)]
    print(f"High quality assemblies (N50 >10kb, contigs <200): {len(high_quality)}")
    
    # Medium quality: N50 > 5kb or contigs < 300
    medium_quality = stats_df[((stats_df['n50'] > 5000) | (stats_df['num_contigs'] < 300)) & 
                             ~((stats_df['n50'] > 10000) & (stats_df['num_contigs'] < 200))]
    print(f"Medium quality assemblies: {len(medium_quality)}")
    
    # Low quality: remaining
    low_quality = len(stats_df) - len(high_quality) - len(medium_quality)
    print(f"Lower quality assemblies: {low_quality}")

def main():
    # Define paths
    results_dir = Path("MAGs/Consolidated .fna/results")
    output_dir = results_dir / "basic_stats"
    
    print("Loading basic genome statistics...")
    
    # Load basic statistics
    stats_df = load_basic_stats(results_dir)
    
    if stats_df.empty:
        print("No basic statistics found!")
        return
    
    print(f"Loaded statistics for {len(stats_df)} MAGs")
    
    # Create visualizations
    print("\nCreating visualizations...")
    create_basic_stats_visualizations(stats_df, output_dir)
    
    # Create summary table
    print("Creating summary table...")
    summary_df = create_summary_table(stats_df, output_dir)
    
    # Print detailed statistics
    print_detailed_stats(stats_df)
    
    print(f"\nResults saved to: {output_dir}")
    print("Generated files:")
    print("  - basic_genome_statistics.png: Comprehensive genome statistics visualization")
    print("  - genome_statistics_summary_table.png: Summary statistics table")
    print("  - genome_statistics_summary.csv: Summary statistics in CSV format")
    
    return stats_df, summary_df

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
"""
Publication Figure Generation: GC Content Thermodynamics
Generates a bar chart comparing the GC percentage of specific spacer architectures 
against the global promoter baseline.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

def main():
    # Architectural Categories
    labels = ['24 bp Dof-Dof\n(Flexible Hinge)', '20 bp bZIP-Dof\n(Rigid Tether)']
    
    # Exact GC percentages established in your previous analysis
    gc_percentages = [30.95, 34.52]
    global_baseline = 31.70
    
    x = np.arange(len(labels))
    width = 0.5
    
    # Create the high-resolution figure
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    # Color bars based on relation to baseline (Blue for flexible/melted, Red for rigid)
    colors = ['#2c7bb6' if val < global_baseline else '#d7191c' for val in gc_percentages]
    
    bars = ax.bar(x, gc_percentages, width, color=colors, edgecolor='black', linewidth=1.5)
    
    # Draw the global baseline
    ax.axhline(global_baseline, color='black', linestyle='--', linewidth=2, label=f'Global Promoter Baseline ({global_baseline}%)')
    
    # Formatting and labels
    ax.set_ylabel('Spacer GC Content (%)', fontsize=12, fontweight='bold')
    ax.set_title('Thermodynamic Divergence of Regulatory Spacer Elements', fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11, fontweight='bold')
    ax.set_ylim(25, 45) # Zoom in to highlight the divergence
    
    # Add exact percentages on top of the bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
                
    ax.legend(fontsize=11, loc='upper left')
    
    plt.tight_layout()
    output_path = os.path.expanduser("~/promoter_project/Figure_GC_Content.pdf")
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    
    print(f"GC Content figure successfully generated and saved to: {output_path}")

if __name__ == "__main__":
    main()

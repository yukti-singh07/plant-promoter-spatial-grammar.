"""
Publication Figure Generation: Position Probability Matrix (PPM) Heatmap
Visualizes the nucleotide frequencies across the entire length of the 
regulatory complexes, proving the degenerate nature of the spacers.
"""

import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter

def calculate_ppm(sequences, length):
    """Calculates the Position Probability Matrix (A, C, G, T frequencies)."""
    if not sequences:
        return None
        
    total_seqs = len(sequences)
    # Initialize dataframe with 0s
    ppm = pd.DataFrame(0.0, index=['A', 'C', 'G', 'T'], columns=range(1, length + 1))
    
    for seq in sequences:
        for i, nt in enumerate(seq):
            if nt in ppm.index:
                ppm.loc[nt, i + 1] += 1
                
    # Normalize to probabilities (0.0 to 1.0)
    ppm = ppm / total_seqs
    return ppm

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    
    targets = {
        "20 bp bZIP-Dof (Rigid Tether)": {
            "pattern": r"(ACGT.{20}AAAG|AAAG.{20}ACGT|ACGT.{20}CTTT|CTTT.{20}ACGT)",
            "length": 28
        },
        "24 bp Dof-Dof (Flexible Hinge)": {
            "pattern": r"(AAAG.{24}AAAG|CTTT.{24}CTTT)",
            "length": 32
        }
    }
    
    extracted_sequences = {name: [] for name in targets.keys()}
    
    print("Extracting sequences and generating PPM heatmaps...")
    
    for crop in crops:
        fasta_file = os.path.join(base_dir, crop, "true_promoters.fasta")
        if not os.path.exists(fasta_file):
            continue
            
        with open(fasta_file, 'r') as f:
            content = f.read().split('>')[1:]
            
        for entry in content:
            lines = entry.strip().split('\n')
            if not lines: continue
            seq = "".join(lines[1:]).upper()
            
            for name, data in targets.items():
                for match in re.finditer(data["pattern"], seq):
                    extracted_sequences[name].append(match.group(0))

    # Create a 2-panel figure (one above the other)
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), dpi=300, gridspec_kw={'height_ratios': [1, 1]})
    
    for i, (name, data) in enumerate(targets.items()):
        seqs = extracted_sequences[name]
        length = data["length"]
        
        ppm = calculate_ppm(seqs, length)
        
        if ppm is not None:
            # Generate heatmap
            sns.heatmap(
                ppm, 
                ax=axes[i],
                cmap="Blues",       # Dark blue means high conservation, white means low
                cbar_kws={'label': 'Frequency'},
                linewidths=0.5,
                linecolor='lightgrey',
                vmin=0, vmax=1      # Lock scale from 0% to 100%
            )
            
            axes[i].set_title(f'{name} Nucleotide Probability Matrix (n={len(seqs)})', fontsize=12, fontweight='bold')
            axes[i].set_xlabel('Nucleotide Position', fontsize=10)
            axes[i].set_ylabel('Base', fontsize=10, fontweight='bold')
            # Rotate Y-axis labels so letters are upright
            axes[i].tick_params(axis='y', rotation=0)
            
    plt.tight_layout()
    output_path = os.path.expanduser("~/promoter_project/Figure_Motif_PPM_Heatmap.pdf")
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    
    print(f"PPM Heatmap generated and saved to: {output_path}")

if __name__ == "__main__":
    main()

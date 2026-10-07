"""
Publication Figure Generation: TSS Anchoring
Maps coordinates and generates a density distribution plot of motif distances.
"""

import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    
    targets = {
        "20 bp bZIP-Dof (Rigid Tether)": r"(ACGT.{20}AAAG|AAAG.{20}ACGT|ACGT.{20}CTTT|CTTT.{20}ACGT)",
        "24 bp Dof-Dof (Flexible Hinge)": r"(AAAG.{24}AAAG|CTTT.{24}CTTT)" 
    }
    
    data = []
    
    print("Extracting coordinates and generating density plot...")
    
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
            seq_length = len(seq)
            
            for name, pattern in targets.items():
                for match in re.finditer(pattern, seq):
                    dist_to_tss = seq_length - match.start()
                    # Limit to typical promoter bounds for clean visualization
                    if dist_to_tss <= 2000: 
                        data.append({"Complex": name, "Distance to TSS (bp)": dist_to_tss})
                        
    df = pd.DataFrame(data)
    
    # Create the high-resolution figure
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Plot Kernel Density Estimate with fill
    sns.kdeplot(
        data=df, 
        x="Distance to TSS (bp)", 
        hue="Complex", 
        fill=True, 
        common_norm=False, 
        palette=["#d7191c", "#2c7bb6"],
        alpha=0.5, 
        linewidth=2
    )
    
    # Add vertical lines at the medians
    plt.axvline(475, color='#d7191c', linestyle='--', label='20 bp Median (475 bp)')
    plt.axvline(502, color='#2c7bb6', linestyle='--', label='24 bp Median (502 bp)')
    
    # Formatting - X-axis now tightly cropped to the data limit
    plt.xlim(0, 1050) 
    plt.xlabel('Distance Upstream from Transcription Start Site (bp)', fontsize=12, fontweight='bold')
    plt.ylabel('Relative Frequency / Motif Density', fontsize=12, fontweight='bold')
    plt.title('Proximal Promoter Anchoring of Stress-Responsive Complexes', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    
    plt.tight_layout()
    output_path = os.path.expanduser("~/promoter_project/Figure_TSS_Anchoring.pdf")
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    
    print(f"Density plot successfully generated and saved to: {output_path}")

if __name__ == "__main__":
    main()

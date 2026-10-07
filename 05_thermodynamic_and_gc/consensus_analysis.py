"""
Motif Consensus Analysis
Extracts the exact sequence of the motif complexes and calculates 
the consensus nucleotide at 50%, 70%, and 90% conservation cutoffs.
"""

import os
import re
from collections import Counter

def calculate_consensus(sequences, length, cutoffs=[0.5, 0.7, 0.9]):
    total_seqs = len(sequences)
    if total_seqs == 0:
        return {}

    # Initialize a list of Counters for each position
    position_counts = [Counter() for _ in range(length)]
    
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            position_counts[i][nucleotide] += 1

    consensus_results = {}
    
    for cutoff in cutoffs:
        consensus_string = []
        for counts in position_counts:
            # Find the most common nucleotide at this position
            most_common_nt, count = counts.most_common(1)[0]
            frequency = count / total_seqs
            
            # If the frequency meets or exceeds the cutoff, use the nucleotide; otherwise use 'N'
            if frequency >= cutoff:
                consensus_string.append(most_common_nt)
            else:
                consensus_string.append('N')
                
        consensus_results[int(cutoff * 100)] = "".join(consensus_string)
        
    return consensus_results

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    
    # Define patterns (4bp motif + spacer + 4bp motif)
    # 20bp complex total length = 28; 24bp complex total length = 32
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
    
    print("Extracting full sequences and calculating consensus thresholds...\n")
    
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

    # Calculate and print consensus for each complex
    for name, data in targets.items():
        seqs = extracted_sequences[name]
        length = data["length"]
        
        print(f"=== {name} ===")
        print(f"Total sequences analyzed: {len(seqs)}")
        
        consensus_strings = calculate_consensus(seqs, length)
        
        if consensus_strings:
            print(f"50% Cutoff (Majority): {consensus_strings[50]}")
            print(f"70% Cutoff (Conserved): {consensus_strings[70]}")
            print(f"90% Cutoff (Strict):    {consensus_strings[90]}\n")
        else:
            print("No sequences found.\n")

if __name__ == "__main__":
    main()

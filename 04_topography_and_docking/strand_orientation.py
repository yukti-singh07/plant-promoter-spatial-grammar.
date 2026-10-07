"""
Strand Orientation (3D Docking Angle) Analysis
Evaluates the directionality (+ vs - strand) of the statistically enriched 
motif pairings (20 bp bZIP-Dof and 24 bp Dof-Dof).
"""

import os
import re

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    
    # We test the enriched geometries using AAAG (forward) and CTTT (reverse)
    targets = [
        {"category": "20 bp Heterotypic (bZIP upstream)", "name": "ACGT ... AAAG (+/+)", "pattern": r"ACGT.{20}AAAG"},
        {"category": "20 bp Heterotypic (bZIP upstream)", "name": "ACGT ... CTTT (+/-)", "pattern": r"ACGT.{20}CTTT"},
        
        {"category": "20 bp Heterotypic (Dof upstream)", "name": "AAAG ... ACGT (+/+)", "pattern": r"AAAG.{20}ACGT"},
        {"category": "20 bp Heterotypic (Dof upstream)", "name": "CTTT ... ACGT (-/+)", "pattern": r"CTTT.{20}ACGT"},
        
        {"category": "24 bp Homotypic (Dof-Dof)", "name": "AAAG ... AAAG (Head-to-Tail, +/+)", "pattern": r"AAAG.{24}AAAG"},
        {"category": "24 bp Homotypic (Dof-Dof)", "name": "AAAG ... CTTT (Head-to-Head, +/-)", "pattern": r"AAAG.{24}CTTT"},
        {"category": "24 bp Homotypic (Dof-Dof)", "name": "CTTT ... AAAG (Tail-to-Tail, -/+)", "pattern": r"CTTT.{24}AAAG"},
        {"category": "24 bp Homotypic (Dof-Dof)", "name": "CTTT ... CTTT (Tail-to-Head, -/-)", "pattern": r"CTTT.{24}CTTT"}
    ]
    
    results = {target["name"]: {"count": 0, "category": target["category"]} for target in targets}
    
    print("Scanning genomic geometries for 3D strand orientation...")
    
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
            
            for target in targets:
                # Using lookahead to catch any overlapping structural constraints
                pattern = f"(?=({target['pattern']}))"
                matches = re.findall(pattern, seq)
                results[target["name"]]["count"] += len(matches)
                
    print("\n=== STRAND ORIENTATION & 3D DOCKING RESULTS ===\n")
    
    current_category = ""
    for name, data in results.items():
        if data["category"] != current_category:
            print(f"[{data['category']}]")
            current_category = data["category"]
        print(f"  {name}: {data['count']} occurrences")
        
    print("\nNote: ACGT is a perfect palindrome, making its innate strand indiscernible.")
    print("Docking polarity is determined entirely by the Dof (AAAG/CTTT) orientation.")

if __name__ == "__main__":
    main()

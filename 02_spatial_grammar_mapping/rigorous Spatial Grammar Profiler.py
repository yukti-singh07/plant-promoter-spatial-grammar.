import os
import re
import csv

crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]
base_dir = os.path.expanduser("~/promoter_project")
output_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")

def parse_fasta(file_path):
    with open(file_path, 'r') as f:
        header, seq = '', []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "Species", "Gene_ID", "Total_ACGT", "Total_AAAG",
        "ACGT_ACGT_Count", "ACGT_ACGT_Spacers",
        "ACGT_AAAG_Count", "ACGT_AAAG_Spacers",
        "AAAG_AAAG_Count", "AAAG_AAAG_Spacers",
        "AAAG_ACGT_Count", "AAAG_ACGT_Spacers"
    ])

    for crop in crops:
        fasta_path = os.path.join(base_dir, crop, "true_promoters.fasta")
        if not os.path.exists(fasta_path):
            continue
            
        print(f"Processing {crop}...")
        
        for header, seq in parse_fasta(fasta_path):
            seq = seq.upper()
            
            # 1. Map every motif coordinate
            acgt_matches = [(m.start(), 'ACGT') for m in re.finditer(r'ACGT', seq)]
            aaag_matches = [(m.start(), 'AAAG') for m in re.finditer(r'AAAG', seq)]
            
            # 2. Sort linearly 5' to 3'
            all_motifs = sorted(acgt_matches + aaag_matches, key=lambda x: x[0])
            
            spacers = {'ACGT_ACGT': [], 'ACGT_AAAG': [], 'AAAG_AAAG': [], 'AAAG_ACGT': []}
            
            # 3. Calculate distances with strict 30bp cutoff
            for i in range(len(all_motifs)):
                for j in range(i + 1, len(all_motifs)):
                    pos1, motif1 = all_motifs[i]
                    pos2, motif2 = all_motifs[j]
                    
                    dist = pos2 - pos1 - 4 # Subtract 4bp motif length
                    
                    # If the next motif is more than 30bp away, stop checking. 
                    # Because they are sorted, everything else will be even further away.
                    if dist > 30:
                        break 
                        
                    if dist >= 0:
                        pair_key = f"{motif1}_{motif2}"
                        spacers[pair_key].append(str(dist))
            
            # 4. Write data
            writer.writerow([
                crop, header, len(acgt_matches), len(aaag_matches),
                len(spacers['ACGT_ACGT']), ";".join(spacers['ACGT_ACGT']) if spacers['ACGT_ACGT'] else "NA",
                len(spacers['ACGT_AAAG']), ";".join(spacers['ACGT_AAAG']) if spacers['ACGT_AAAG'] else "NA",
                len(spacers['AAAG_AAAG']), ";".join(spacers['AAAG_AAAG']) if spacers['AAAG_AAAG'] else "NA",
                len(spacers['AAAG_ACGT']), ";".join(spacers['AAAG_ACGT']) if spacers['AAAG_ACGT'] else "NA"
            ])

print(f"\nAnalysis complete! Results saved to {output_csv}")
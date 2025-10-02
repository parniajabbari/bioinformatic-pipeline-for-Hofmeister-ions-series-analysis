import csv
import os

# File paths
input_csv = '/home/anto/parnia/F/result.csv'
dssp_dir = '/home/anto/parnia/F/dssp_results'
output_csv = '/home/anto/parnia/F/result.csv'

# Residue code mapping
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
dssp_structure_map = {
    'H': 'Alpha helix',
    'B': 'Beta bridge',
    'E': 'Extended strand',
    'G': '3-10 helix',
    'I': 'Pi helix',
    'T': 'Turn',
    'S': 'Bend',
    ' ': 'Coil',
    'PP': 'Polyproline II helix',
    'P><': 'Polyproline II helix',  
    'P': 'Polyproline II helix'
}


# Helper function to parse DSSP
def parse_dssp(file_path):
    dssp_data = {}
    with open(file_path, 'r') as file:
        start_reading = False
        for line in file:
            if line.startswith("  #  RESIDUE AA STRUCTURE"):
                start_reading = True
                continue
            if start_reading:
                residue_id = line[5:10].strip()  # Residue ID
                chain_id = line[11].strip()  # Chain ID
                aa = line[13].strip()  # Amino acid
                secondary_structure = line[16].strip()  # Secondary structure
                dssp_data[(chain_id, residue_id)] = {
                    'aa': aa,
                    'secondary_structure': secondary_structure
                }
    return dssp_data

# Read input CSV
with open(input_csv, 'r') as file:
    reader = csv.DictReader(file)
    input_data = list(reader)

# Process each row and find secondary structure
for row in input_data:
    pdb_id = row['PDB ID']
    residue_name = row['Residue Name']
    residue_id = row['Residue ID'].strip()
    chain_id = row['Chain ID']

    # Skip water molecules (HOH)
    if residue_name == 'HOH':
        row['Secondary Structure'] = 'N/A'
        continue

    if residue_name in three_to_one:
        residue_name_one_letter = three_to_one[residue_name]
    else:
        row['Secondary Structure'] = 'N/A'  # Non-residues will be marked as N/A
        continue

    # Locate corresponding DSSP file
    dssp_file_path = os.path.join(dssp_dir, f"{pdb_id}.dssp")
    if not os.path.exists(dssp_file_path):
        row['Secondary Structure'] = 'DSSP file not found'
        continue

    # Parse DSSP file
    dssp_data = parse_dssp(dssp_file_path)

    # Find matching residue
    key = (chain_id, residue_id)
    if key in dssp_data:
        dssp_info = dssp_data[key]
        if dssp_info['aa'] == residue_name_one_letter:
            secondary_structure = dssp_info['secondary_structure']
            row['Secondary Structure'] = dssp_structure_map.get(secondary_structure, 'Coil')  # Use "Coil" for empty structures
        else:
            row['Secondary Structure'] = f'Residue mismatch in DSSP (Expected: {residue_name_one_letter}, Found: {dssp_info["aa"]})'
    else:
        row['Secondary Structure'] = f'Residue not in DSSP for key: {key}'

# Write output CSV
output_columns = list(input_data[0].keys())
with open(output_csv, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=output_columns)
    writer.writeheader()
    writer.writerows(input_data)

print(f"Secondary structure information saved to {output_csv}.")

import csv
import os
import argparse

# -----------------------------
# 1. Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Add secondary structure information (from DSSP) to the ion result CSV."
)
parser.add_argument(
    "ion",
    type=str,
    help="Ion symbol (e.g., F, CL, SCN, SO4, etc.)"
)
parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Optional output directory (default: ~/Desktop/{ION}/)"
)
args = parser.parse_args()
ION_SYMBOL = args.ion.upper()

# -----------------------------
# 2. Define directories
# -----------------------------
BASE_DESKTOP = os.path.expanduser("~/Desktop")
BASE_DIR = args.output or os.path.join(BASE_DESKTOP, ION_SYMBOL)

input_csv = os.path.join(BASE_DIR, "result.csv")
dssp_dir = os.path.join(BASE_DIR, "dssp_results")
output_csv = os.path.join(BASE_DIR, "result_with_dssp.csv")

# -----------------------------
# 3. Residue mappings
# -----------------------------
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

# -----------------------------
# 4. Helper to parse DSSP files
# -----------------------------
def parse_dssp(file_path):
    dssp_data = {}
    with open(file_path, 'r') as file:
        start_reading = False
        for line in file:
            if line.startswith("  #  RESIDUE AA STRUCTURE"):
                start_reading = True
                continue
            if start_reading:
                residue_id = line[5:10].strip()
                chain_id = line[11].strip()
                aa = line[13].strip()
                secondary_structure = line[16].strip()
                dssp_data[(chain_id, residue_id)] = {
                    'aa': aa,
                    'secondary_structure': secondary_structure
                }
    return dssp_data

# -----------------------------
# 5. Main processing
# -----------------------------
def main():
    if not os.path.exists(input_csv):
        print(f"⚠️ Input CSV not found: {input_csv}")
        return
    if not os.path.exists(dssp_dir):
        print(f"⚠️ DSSP directory not found: {dssp_dir}")
        return

    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        input_data = list(reader)

    for row in input_data:
        pdb_id = row.get('PDB ID', '').strip()
        residue_name = row.get('Residue Name', '').strip()
        residue_id = row.get('Residue ID', '').strip()
        chain_id = row.get('Chain ID', '').strip()

        if residue_name == 'HOH':
            row['Secondary Structure'] = 'N/A'
            continue

        if residue_name not in three_to_one:
            row['Secondary Structure'] = 'N/A'
            continue

        residue_name_one = three_to_one[residue_name]
        dssp_file_path = os.path.join(dssp_dir, f"{pdb_id}.dssp")

        if not os.path.exists(dssp_file_path):
            row['Secondary Structure'] = 'DSSP file not found'
            continue

        dssp_data = parse_dssp(dssp_file_path)
        key = (chain_id, residue_id)

        if key in dssp_data:
            dssp_info = dssp_data[key]
            if dssp_info['aa'] == residue_name_one:
                ss = dssp_info['secondary_structure']
                row['Secondary Structure'] = dssp_structure_map.get(ss, 'Coil')
            else:
                row['Secondary Structure'] = (
                    f"Residue mismatch (Expected: {residue_name_one}, Found: {dssp_info['aa']})"
                )
        else:
            row['Secondary Structure'] = f"Residue not in DSSP ({key})"

    output_columns = list(input_data[0].keys())

    with open(output_csv, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=output_columns)
        writer.writeheader()
        writer.writerows(input_data)

    print(f"✅ Secondary structure data saved to: {output_csv}")

if __name__ == "__main__":
    main()

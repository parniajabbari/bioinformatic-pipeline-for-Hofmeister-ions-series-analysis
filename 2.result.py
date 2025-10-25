import os
import argparse
from Bio.PDB import MMCIFParser
import numpy as np
import csv

# -----------------------------
# 1. Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Analyze Hofmeister ion interactions in CIF structures."
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
OUTPUT_DIRECTORY = args.output or os.path.join(BASE_DESKTOP, ION_SYMBOL)
CIF_DIRECTORY = os.path.join(OUTPUT_DIRECTORY, "cif")

os.makedirs(CIF_DIRECTORY, exist_ok=True)

# ⚠️ Always save to result.csv
output_csv_path = os.path.join(OUTPUT_DIRECTORY, "result.csv")

# -----------------------------
# 3. Core analysis functions
# -----------------------------
def calculate_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def calculate_angle(a, b, c):
    ab = np.array(b) - np.array(a)
    bc = np.array(c) - np.array(b)
    cos_angle = np.dot(ab, bc) / (np.linalg.norm(ab) * np.linalg.norm(bc))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))

def classify_ion_interaction(dist, interaction_type):
    if dist < 4 and interaction_type == "N/A":
        return "Cofactor/ligand"
    return interaction_type

# -----------------------------
# 4. Main analysis routine
# -----------------------------
def main():
    cif_files = [os.path.join(CIF_DIRECTORY, f) for f in os.listdir(CIF_DIRECTORY) if f.endswith('.cif')]

    if not cif_files:
        print(f"⚠️ No CIF files found in {CIF_DIRECTORY}")
        return

    with open(output_csv_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "PDB ID", f"{ION_SYMBOL} Full Name", f"{ION_SYMBOL} Atom Name",
            "Residue Atom Name", "Residue Name", "Residue ID",
            "Distance (Å)", "Angle (°)", "Interaction Type",
            f"{ION_SYMBOL} Atom Coordinates", "Residue Atom Coordinates", "Chain ID"
        ])

        for file_path in cif_files:
            pdb_id = os.path.splitext(os.path.basename(file_path))[0]
            parser = MMCIFParser(QUIET=True)

            try:
                structure = parser.get_structure(pdb_id, file_path)
            except Exception as e:
                print(f"❌ Error processing {pdb_id}: {e}")
                continue

            ion_atoms = []
            other_atoms = []

            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    for residue in chain:
                        for atom in residue:
                            if (residue.resname.strip().upper() == ION_SYMBOL or
                                atom.element.strip().upper() == ION_SYMBOL):
                                ion_atoms.append({
                                    'name': atom.get_name(),
                                    'coord': atom.coord,
                                    'chain': chain_id,
                                    'full_id': f"('{residue.resname}', {residue.id[1]}, '{chain_id}')"
                                })
                            else:
                                other_atoms.append({
                                    'name': atom.get_name(),
                                    'coord': atom.coord,
                                    'resname': residue.resname,
                                    'res_id': residue.id[1],
                                    'chain': chain_id
                                })

            print(f"🔍 {pdb_id}: Found {len(ion_atoms)} {ION_SYMBOL} atoms and {len(other_atoms)} other atoms.")

            results = []

            # --- Ion–Ion interactions ---
            for i, I1 in enumerate(ion_atoms):
                for I2 in ion_atoms[i + 1:]:
                    if I1['full_id'] != I2['full_id']:
                        dist = calculate_distance(I1['coord'], I2['coord'])
                        if dist < 4:
                            interaction_type = "Ion"
                            results.append([
                                I1['full_id'], I1['name'], I2['name'], ION_SYMBOL, "N/A",
                                f"{dist:.2f}", "N/A", interaction_type,
                                str(I1['coord'].tolist()), str(I2['coord'].tolist()), I1['chain']
                            ])

            # --- Ion–Other interactions ---
            for I in ion_atoms:
                for other in other_atoms:
                    dist = calculate_distance(I['coord'], other['coord'])
                    if dist > 8:
                        continue

                    interaction_type = "N/A"
                    angle = None

                    if other['resname'] == 'PEG':
                        interaction_type = classify_ion_interaction(dist, interaction_type)

                    hofmeister_ions = {
                        'CA': 1.00, 'MG': 0.72, 'NA': 1.02, 'K': 1.38, 'SR': 1.24, 'BA': 1.43,
                        'CL': 1.67, 'BR': 1.82, 'I': 2.06, 'SO4': 2.58, 'NO3': 1.90
                    }
                    ion_radius = hofmeister_ions.get(other['resname'].upper(), 0)
                    if ion_radius > 0 and dist <= 4:
                        interaction_type = "Ion"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        metals = ['ZN', 'MG', 'FE', 'CA', 'CU', 'NI', 'MN', 'CO', 'NA', 'K', 'CD',
                                  'AL', 'AG', 'AU', 'PT', 'LI', 'SR', 'BA', 'RB']
                        if other['resname'].upper() in metals:
                            interaction_type = "Interaction with Metal"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        if other['resname'].upper() in ['ASN', 'LEU', 'ALA', 'CYS', 'GLN', 'GLY',
                                                        'ILE', 'MET', 'PRO', 'SER', 'THR', 'VAL', 'OH']:
                            interaction_type = "Aliphatic"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        if other['resname'].upper() in ['PHE', 'TYR', 'TRP']:
                            interaction_type = "Aromatic"

                    if interaction_type == "N/A":
                        ion_charge_dict = {'F': -1, 'CL': -1, 'BR': -1, 'I': -1}
                        ion_overall_charge = ion_charge_dict.get(ION_SYMBOL, None)
                        if ion_overall_charge is not None:
                            specific_neg_atoms = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
                            specific_pos_atoms = {'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ']}
                            resname = other['resname'].upper()
                            atom_name = other['name'].upper()
                            if ion_overall_charge > 0:
                                if resname in specific_neg_atoms and atom_name in specific_neg_atoms[resname]:
                                    if 2.2 <= dist <= 4.0:
                                        interaction_type = "Salt Bridge"
                                elif 3.5 <= dist <= 4.0:
                                    interaction_type = "Salt Bridge"
                            elif ion_overall_charge < 0:
                                if resname in specific_pos_atoms and atom_name in specific_pos_atoms[resname]:
                                    if 2.2 <= dist <= 4.0:
                                        interaction_type = "Salt Bridge"
                                elif 3.5 <= dist <= 4.0:
                                    interaction_type = "Salt Bridge"

                    if interaction_type == "N/A" and 2.2 <= dist <= 3.5:
                        for other2 in other_atoms:
                            if other['coord'].tolist() != other2['coord'].tolist() or other['name'] != other2['name']:
                                angle = calculate_angle(I['coord'], other['coord'], other2['coord'])
                                if angle < 135:
                                    if other['resname'].upper() == "HOH":
                                        interaction_type = "H-bond with Water"
                                    elif other['name'] in ['C', 'N', 'O', 'CA']:
                                        interaction_type = "H-bond with Main Chain"
                                    else:
                                        interaction_type = "H-bond with Side Chain"

                    results.append([
                        I['full_id'], I['name'], other['name'], other['resname'], other['res_id'],
                        f"{dist:.2f}", angle if angle else "N/A",
                        interaction_type, str(I['coord'].tolist()), str(other['coord'].tolist()), I['chain']
                    ])

            for row in results:
                writer.writerow([pdb_id] + row)

            print(f"✅ Processed {pdb_id} with {len(results)} interactions.")

    print(f"\n📁 Results saved to: {output_csv_path}")

if __name__ == "__main__":
    main()

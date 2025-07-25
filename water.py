# --- Import Required Libraries ---
import pandas as pd
import matplotlib.pyplot as plt
from Bio import PDB
import os
import numpy as np

# --- Define File Paths ---
cif_directory = "/home/parnia/parnia/F/cif"  # Directory containing CIF files
F_count_file = "/home/parnia/parnia/F/F_ion_information_filtered.csv"  # CSV with F counts

# --- Load and Prepare F Count Data ---
F_data = pd.read_csv(F_count_file, header=0).iloc[:-1]  # Skip the final 'Total' row
F_data = F_data.iloc[:, [0, -1]]  # Use first (PDB ID) and last (F count) columns
F_data.columns = ['pdb_id', 'F_count']

# --- Parameters ---
distance_cutoff = 4.0     # Maximum distance for counting water
water_molecules_per_F = []  # Stores water counts per F binding site

# --- Utility Functions ---

def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms."""
    return np.linalg.norm(atom1.coord - atom2.coord)

def count_nearby_water_molecules(pdb_id):
    """
    Count water molecules (HOH) within ≤4.0 Å of each F site in the structure.
    """
    file_path = os.path.join(cif_directory, f"{pdb_id}.cif")
    parser = PDB.MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, file_path)
    except Exception as e:
        print(f"[ERROR] Failed to parse {pdb_id}: {e}")
        return []

    F_water_counts = []

    for model in structure:
        for chain in model:
            F_residues = [res for res in chain if res.get_resname() == 'F']
            HOH_residues = [res for res in chain if res.get_resname() == 'HOH']

            for F in F_residues:
                F_water_count = 0

                for hoh in HOH_residues:
                    for F_atom in F:
                        for hoh_atom in hoh:
                            distance = calculate_distance(F_atom, hoh_atom)
                            if distance <= distance_cutoff:
                                F_water_count += 1
                                break  # Count each HOH once per F

                F_water_counts.append(F_water_count)

    return F_water_counts

# --- Process Structures ---

for _, row in F_data.iterrows():
    pdb_id = str(row['pdb_id']).strip()

    try:
        F_count = int(row['F_count'])
    except (ValueError, TypeError):
        print(f"[WARNING] Invalid F count in row: {row}")
        continue

    water_counts = count_nearby_water_molecules(pdb_id)
    water_molecules_per_F.extend(water_counts)

# --- Build Water Distribution Histogram (0–5+ bins) ---

water_distribution = [0] * 6  # Index: 0-5 (5+ goes to index 5)

for count in water_molecules_per_F:
    if count <= 5:
        water_distribution[count] += 1
    else:
        water_distribution[5] += 1

# --- Plot Results ---

plt.figure(figsize=(10, 6))
bars = plt.bar(range(6), water_distribution, edgecolor='black', align='center')

plt.title('Distribution of Water Molecules Around F Binding Sites (≤ 4 Å)')
plt.xlabel('Water Molecules per F Site')
plt.ylabel('Number of F Binding Sites')
plt.xticks(range(6), ['0', '1', '2', '3', '4', '5+'])

# Annotate bars with counts
for bar, count in zip(bars, water_distribution):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), str(count),
             ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.show()



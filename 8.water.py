import pandas as pd
import matplotlib.pyplot as plt
from Bio import PDB
import os
import numpy as np

# ============================
# --- USER PARAMETERS ---
# ============================

ION_SYMBOL = "GAI"  # <-- Change this to your ion symbol, e.g. "GAI", "SO4", "MG", etc.

# ✅ Fixed f-string formatting
cif_directory = f"/Users/respina/desktop/{ION_SYMBOL}/cif"  # Directory containing CIF files
ion_count_file = f"/Users/respina/desktop/{ION_SYMBOL}/{ION_SYMBOL}_ion_information_filtered.csv"  # CSV with ion counts
distance_cutoff = 4.0  # Maximum distance (Å) to count nearby water

# ============================
# --- LOAD ION COUNT DATA ---
# ============================

ion_data = pd.read_csv(ion_count_file, header=0).iloc[:-1]  # Skip the final 'Total' row
ion_data = ion_data.iloc[:, [0, -1]]  # Use first (PDB ID) and last (ion count) columns
ion_data.columns = ['pdb_id', 'ion_count']

# ============================
# --- FUNCTIONS ---
# ============================

def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms."""
    return np.linalg.norm(atom1.coord - atom2.coord)

def count_nearby_water_molecules(pdb_id, ion_symbol):
    """
    Count water molecules (HOH) within ≤4.0 Å of each ion site in the structure.
    """
    file_path = os.path.join(cif_directory, f"{pdb_id}.cif")
    parser = PDB.MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, file_path)
    except Exception as e:
        print(f"[ERROR] Failed to parse {pdb_id}: {e}")
        return []

    ion_water_counts = []

    for model in structure:
        for chain in model:
            ion_residues = [res for res in chain if res.get_resname().strip() == ion_symbol]
            HOH_residues = [res for res in chain if res.get_resname() == 'HOH']

            for ion in ion_residues:
                water_count = 0

                for hoh in HOH_residues:
                    for ion_atom in ion:
                        for hoh_atom in hoh:
                            distance = calculate_distance(ion_atom, hoh_atom)
                            if distance <= distance_cutoff:
                                water_count += 1
                                break  # Count each HOH once per ion site

                ion_water_counts.append(water_count)

    return ion_water_counts

# ============================
# --- PROCESS STRUCTURES ---
# ============================

water_molecules_per_ion = []

for _, row in ion_data.iterrows():
    pdb_id = str(row['pdb_id']).strip()

    try:
        ion_count = int(row['ion_count'])
    except (ValueError, TypeError):
        print(f"[WARNING] Invalid ion count in row: {row}")
        continue

    water_counts = count_nearby_water_molecules(pdb_id, ION_SYMBOL)
    water_molecules_per_ion.extend(water_counts)

# ============================
# --- BUILD WATER DISTRIBUTION ---
# ============================

water_distribution = [0] * 6  # Bins for 0–5+

for count in water_molecules_per_ion:
    if count <= 5:
        water_distribution[count] += 1
    else:
        water_distribution[5] += 1

# ============================
# --- PLOT RESULTS ---
# ============================

plt.figure(figsize=(10, 6))
bars = plt.bar(range(6), water_distribution, edgecolor='black', align='center')

plt.title(f'Distribution of Water Molecules Around {ION_SYMBOL} Binding Sites (≤ {distance_cutoff} Å)')
plt.xlabel(f'Water Molecules per {ION_SYMBOL} Site')
plt.ylabel(f'Number of {ION_SYMBOL} Binding Sites')
plt.xticks(range(6), ['0', '1', '2', '3', '4', '5+'])

# Annotate bars with counts
for bar, count in zip(bars, water_distribution):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
             ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.show()


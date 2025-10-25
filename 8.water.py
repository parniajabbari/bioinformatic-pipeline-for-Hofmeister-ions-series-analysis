import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio import PDB
import numpy as np

# ------------------------
# CLI Arguments
# ------------------------
parser = argparse.ArgumentParser(description="Plot water distribution around ions")
parser.add_argument("ion", help="Ion symbol (e.g., GAI, SO4, MG, LCP)")
parser.add_argument("--output", help="Output directory (default: ~/Desktop/<ION>)")
args = parser.parse_args()

ION_SYMBOL = args.ion.upper()
BASE_DIR = args.output or os.path.join(os.path.expanduser("~"), "Desktop", ION_SYMBOL)
os.makedirs(BASE_DIR, exist_ok=True)

# ------------------------
# User Parameters
# ------------------------
CIF_DIR = os.path.join(BASE_DIR, "cif")
ION_COUNT_FILE = os.path.join(BASE_DIR, f"{ION_SYMBOL}_ion_information_filtered.csv")
DISTANCE_CUTOFF = 4.0
OUTPUT_FILE = os.path.join(BASE_DIR, "water_distribution.png")  # fixed name

if not os.path.exists(ION_COUNT_FILE):
    raise FileNotFoundError(f"Ion count CSV not found: {ION_COUNT_FILE}")

# ------------------------
# Load Ion Count Data
# ------------------------
ion_data = pd.read_csv(ION_COUNT_FILE, header=0).iloc[:-1]  # Skip 'Total' row
ion_data = ion_data.iloc[:, [0, -1]]  # PDB ID and last column
ion_data.columns = ['pdb_id', 'ion_count']

# ------------------------
# Functions
# ------------------------
def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms."""
    return np.linalg.norm(atom1.coord - atom2.coord)

def count_nearby_water_molecules(pdb_id, ion_symbol):
    """Count water molecules (HOH) within DISTANCE_CUTOFF Å of each ion site."""
    file_path = os.path.join(CIF_DIR, f"{pdb_id}.cif")
    parser = PDB.MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, file_path)
    except Exception as e:
        print(f"[ERROR] Failed to parse {pdb_id}: {e}")
        return []

    counts = []
    for model in structure:
        for chain in model:
            ion_residues = [res for res in chain if res.get_resname().strip() == ion_symbol]
            HOH_residues = [res for res in chain if res.get_resname() == 'HOH']

            for ion in ion_residues:
                water_count = 0
                for hoh in HOH_residues:
                    for ion_atom in ion:
                        for hoh_atom in hoh:
                            if calculate_distance(ion_atom, hoh_atom) <= DISTANCE_CUTOFF:
                                water_count += 1
                                break
                counts.append(water_count)
    return counts

# ------------------------
# Process Structures
# ------------------------
water_molecules_per_ion = []

for _, row in ion_data.iterrows():
    pdb_id = str(row['pdb_id']).strip()
    try:
        int(row['ion_count'])
    except (ValueError, TypeError):
        print(f"[WARNING] Invalid ion count in row: {row}")
        continue

    water_counts = count_nearby_water_molecules(pdb_id, ION_SYMBOL)
    water_molecules_per_ion.extend(water_counts)

# ------------------------
# Build Water Distribution
# ------------------------
water_distribution = [0] * 6  # 0-5+
for count in water_molecules_per_ion:
    if count <= 5:
        water_distribution[count] += 1
    else:
        water_distribution[5] += 1

# ------------------------
# Plot Results
# ------------------------
plt.figure(figsize=(10, 6))
bars = plt.bar(range(6), water_distribution, edgecolor='black', align='center')

plt.title(f'Distribution of Water Molecules Around {ION_SYMBOL} Binding Sites (≤ {DISTANCE_CUTOFF} Å)')
plt.xlabel(f'Water Molecules per {ION_SYMBOL} Site')
plt.ylabel(f'Number of {ION_SYMBOL} Binding Sites')
plt.xticks(range(6), ['0', '1', '2', '3', '4', '5+'])

# Annotate bars
for bar, count in zip(bars, water_distribution):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
             ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=300)
plt.show()

print(f"✅ Water distribution plot saved to: {OUTPUT_FILE}")



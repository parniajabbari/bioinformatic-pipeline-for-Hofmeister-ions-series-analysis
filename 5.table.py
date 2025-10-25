import os
import pandas as pd
import argparse

# ------------------------
# CLI Arguments
# ------------------------
parser = argparse.ArgumentParser(description="Calculate interaction statistics per ion")
parser.add_argument("ion", help="Ion symbol (e.g., F, SCN, AL, LCP)")
parser.add_argument("--output", help="Output directory (default: ~/Desktop/<ION>)")
args = parser.parse_args()

ION_SYMBOL = args.ion.upper()
BASE_DIR = args.output or os.path.join(os.path.expanduser("~"), "Desktop", ION_SYMBOL)
os.makedirs(BASE_DIR, exist_ok=True)

# Input and output files
INPUT_FILE = os.path.join(BASE_DIR, "result.csv")
OUTPUT_FILE = os.path.join(BASE_DIR, "interaction_statistics_with_residue_counts.csv")  # fixed name

# ------------------------
# Residue Characteristics
# ------------------------
characteristics = {
    'LEU': 'Hydrophobic Aliphatic', 'ILE': 'Hydrophobic Aliphatic', 'MET': 'Hydrophobic Aliphatic',
    'ALA': 'Hydrophobic Aliphatic (Tiny)', 'GLY': 'Hydrophobic Aliphatic (Tiny)',
    'PRO': 'Hydrophobic Aliphatic (Small)', 'VAL': 'Hydrophobic Aliphatic (Small)', 'CYS': 'Hydrophobic Aliphatic (Small)',
    'PHE': 'Hydrophobic Aromatic', 'TYR': 'Hydrophobic Aromatic', 'TRP': 'Hydrophobic Aromatic',
    'LYS': 'Charged Positive', 'ARG': 'Charged Positive', 'HIS': 'Charged Positive',
    'ASP': 'Charged Negative (Small)', 'GLU': 'Charged Negative',
    'SER': 'Polar (Tiny)', 'GLN': 'Polar', 'ASN': 'Polar (Small)', 'THR': 'Polar (Small)',
    'OH': 'Hydroxyl Group', 'HOH': 'Water',
    'ZN': 'Metal', 'MG': 'Metal', 'FE': 'Metal', 'CA': 'Metal', 'CU': 'Metal', 'NI': 'Metal',
    'SCN': 'Ion', 'IOD': 'Ion', 'CL': 'Ion', 'BR': 'Ion', 'F': 'Ion',
    'SO4': 'Ion', 'PO4': 'Ion', 'CO3': 'Ion', 'HCO3': 'Ion',
    'NA': 'Ion', 'K': 'Ion', 'CS': 'Ion', 'LI': 'Ion', 'BA': 'Ion', 'SR': 'Ion',
    'NH4': 'Ion', 'NO3': 'Ion', 'ClO4': 'Ion', 'LCP': 'Ion', "2HP": "Ion",
    'AL': 'Metal', 'FE2': 'Ion', 'FE3': 'Ion', 'AG': 'Ion'
}

# ------------------------
# Load Input Data
# ------------------------
if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"Input CSV not found: {INPUT_FILE}")

data = pd.read_csv(INPUT_FILE)
data.columns = data.columns.str.strip().str.lower()

if 'distance (å)' not in data.columns:
    raise KeyError("Column 'Distance (Å)' not found in input CSV")

# Filter interactions within 4 Å
filtered_data = data[data['distance (å)'] <= 4.0]

# Residue-specific counts
filtered_data['residue specific count'] = filtered_data.groupby('residue name')['residue name'].transform('count')

# Create unique residue ID
filtered_data['residue_id'] = (
    filtered_data['residue name'].astype(str) + "_" +
    filtered_data['residue id'].astype(str) + "_" +
    filtered_data['chain id'].astype(str)
)

# Count unique residues
unique_residue_counts = (
    filtered_data[['residue name', 'residue_id']]
    .drop_duplicates()
    .groupby('residue name')
    .size()
    .to_dict()
)

# Interaction counts
interaction_counts = filtered_data['interaction type'].value_counts()
total_interactions = filtered_data.shape[0]

results = []

for interaction in interaction_counts.index:
    interaction_data = filtered_data[filtered_data['interaction type'] == interaction]
    grouped_residues = interaction_data.groupby(['interaction type', 'residue name'])
    
    for (interaction_type, residue), group in grouped_residues:
        count = group.shape[0]
        residue_specific_count = group['residue specific count'].iloc[0]
        approximate_percentage = round((count / total_interactions) * 100, 2)
        characteristic = characteristics.get(residue, 'Cofactor/Ligand')
        interaction_type = characteristic if characteristic in ['Cofactor/Ligand', 'Metal', 'Ion'] else interaction
        atom_specific_counts = group['residue atom name'].value_counts().to_dict()
        unique_residue_count = unique_residue_counts.get(residue, 0)

        results.append({
            'Interaction Type': interaction_type,
            'Residue Name': residue,
            'Approx. Percentage in Shell': approximate_percentage,
            'atoms of the residue (column B) in interaction type (column A) / Total atoms of all residues in Shell': f"{count} / {total_interactions}",
            'atoms of the residue (column B) in interaction type (column A) / total number of residue (column B) in the shell': f"{count} / {residue_specific_count}",
            'Characteristic': characteristic,
            'Atom-Specific Counts': atom_specific_counts,
            'Unique Residues in Shell': unique_residue_count
        })

# Aggregate results
results_df = pd.DataFrame(results)
results_df_aggregated = results_df.groupby(['Interaction Type', 'Residue Name']).agg({
    'Approx. Percentage in Shell': 'first',
    'atoms of the residue (column B) in interaction type (column A) / Total atoms of all residues in Shell': 'first',
    'atoms of the residue (column B) in interaction type (column A) / total number of residue (column B) in the shell': 'first',
    'Characteristic': 'first',
    'Atom-Specific Counts': lambda x: {k: sum(d.get(k, 0) for d in x) for k in set(k for d in x for k in d)},
    'Unique Residues in Shell': 'first'
}).reset_index()

# Save output with fixed name
results_df_aggregated.to_csv(OUTPUT_FILE, index=False)
print(f"✅ Aggregated statistics saved to: {OUTPUT_FILE}")


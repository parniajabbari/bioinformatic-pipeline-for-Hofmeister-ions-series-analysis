import pandas as pd

# Define residue characteristics
characteristics = {
    'LEU': 'Hydrophobic Aliphatic', 'ILE': 'Hydrophobic Aliphatic', 'MET': 'Hydrophobic Aliphatic',
    'ALA': 'Hydrophobic Aliphatic (Tiny)', 'GLY': 'Hydrophobic Aliphatic (Tiny)',
    'PRO': 'Hydrophobic Aliphatic (Small)', 'VAL': 'Hydrophobic Aliphatic (Small)', 'CYS': 'Hydrophobic Aliphatic (Small)',
    'PHE': 'Hydrophobic Aromatic', 'TYR': 'Hydrophobic Aromatic', 'TRP': 'Hydrophobic Aromatic',
    'LYS': 'Charged Positive', 'ARG': 'Charged Positive', 'HIS': 'Charged Positive',
    'ASP': 'Charged Negative (Small)', 'GLU': 'Charged Negative',
    'SER': 'Polar (Tiny)', 'GLN': 'Polar', 'ASN': 'Polar (Small)', 'THR': 'Polar (Small)',
    'OH': 'Hydroxyl Group', 'HOH': 'Water',
    'ZN': 'Metal', 'MG': 'Metal', 'FE': 'Metal', 'CA': 'Metal', 'CU': 'Metal', 'NI': 'Metal', 'SCN': 'Ion' , 'IOD': 'Ion', 'CL': 'Ion', 'BR': 'Ion', 'F': 'Ion', 
    'SO4': 'Ion', 'PO4': 'Ion', 'CO3': 'Ion', 'HCO3': 'Ion',
    'NA': 'Ion', 'K': 'Ion', 'CS': 'Ion','LI': 'Ion', 'BA': 'Ion', 'SR': 'Ion' ,  
    'NH4': 'Ion', 'NO3': 'Ion', 'SCN': 'Ion', 'ClO4': 'Ion',  
    'FE2': 'Ion', 'FE3': 'Ion', 'ZN': 'Metal', 'CU': 'Metal', 'MN': 'Metal', 
    'CO': 'Metal', 'NI': 'Metal', 'MO': 'Metal', 'W': 'Metal', 'V': 'Metal',  
    'AG': 'Ion' , 'AL': 'Metal', 'CR': 'Metal', 'TI': 'Metal', 'MN': 'Metal', 'CO': 'Metal', 'RU': 'Metal', 'RH': 'Metal', 'PD': 'Metal',
    'OS': 'Metal', 'IR': 'Metal', 'PT': 'Metal', 'AU': 'Metal', 'HG': 'Metal', 'PB': 'Metal', 'BI': 'Metal', 'LA': 'Metal',
    'CE': 'Metal', 'PR': 'Metal', 'ND': 'Metal', 'SM': 'Metal', 'EU': 'Metal', 'GD': 'Metal', 'TB': 'Metal', 'DY': 'Metal',
    'HO': 'Metal', 'ER': 'Metal', 'TM': 'Metal', 'YB': 'Metal', 'LU': 'Metal', 'AC': 'Metal', 'TH': 'Metal', 'PA': 'Metal',
    'U': 'Metal', 'NP': 'Metal', 'PU': 'Metal', 'AM': 'Metal', 'CM': 'Metal', 'BK': 'Metal', 'CF': 'Metal', 'ES': 'Metal',
    'FM': 'Metal', 'MD': 'Metal', 'NO': 'Metal', 'Lr': 'Metal' ,'LCP' : 'Ion' , "2HP": "Ion"
}

# Load the input data
input_file = '/home/parnia/parnia/F/result.csv'
data = pd.read_csv(input_file)

# Clean and standardize column names (convert to lowercase and strip spaces)
data.columns = data.columns.str.strip().str.lower()

# Ensure the 'distance (å)' column exists
if 'distance (å)' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise KeyError("Column 'Distance (Å)' (case-insensitive) not found in the input file.")

# Filter interactions within the 4 Å shell
filtered_data = data[data['distance (å)'] <= 4.0]

# Add a column for residue-specific counts within the 4 Å shell
filtered_data['residue specific count'] = filtered_data.groupby('residue name')['residue name'].transform('count')

# === NEW: Create unique residue ID for identifying unique residues ===
# Make sure your input CSV has: residue_name, residue_number, chain_id
filtered_data['residue_id'] = (
    filtered_data['residue name'].astype(str) + "_" +
    filtered_data['residue id'].astype(str) + "_" +
    filtered_data['chain id'].astype(str)
)

# === NEW: Count unique residues per residue name ===
unique_residue_counts = (
    filtered_data[['residue name', 'residue_id']]
    .drop_duplicates()
    .groupby('residue name')
    .size()
    .to_dict()
)

# Count occurrences of each interaction type, residue name, and residue atom
interaction_counts = filtered_data['interaction type'].value_counts()
residue_counts = filtered_data['residue name'].value_counts()
atom_counts = filtered_data['residue atom name'].value_counts()

# Total number of interactions in the 4 Å shell
total_interactions_in_shell = filtered_data.shape[0]

# Create a list to hold the results
results = []

# Iterate through each interaction type in the shell
for interaction in interaction_counts.index:
    total_interaction_count = interaction_counts[interaction]
    # Filter data for the specific interaction type
    interaction_data = filtered_data[filtered_data['interaction type'] == interaction]
    # Group by 'residue name' to aggregate atom-specific counts
    grouped_residues = interaction_data.groupby(['interaction type', 'residue name'])
    for (interaction_type, residue), group in grouped_residues:
        count = group.shape[0]
        residue_specific_count = group['residue specific count'].iloc[0]
        approximate_percentage = round((count / total_interactions_in_shell) * 100, 2)
        # Determine the characteristic of the residue
        characteristic = characteristics.get(residue, 'Cofactor/Ligand')
        # Update the interaction type based on the characteristic
        interaction_type = characteristic if characteristic in ['Cofactor/Ligand', 'Metal', 'Ion'] else interaction
        # Aggregate atom-specific counts
        atom_specific_counts = group['residue atom name'].value_counts().to_dict()
        # Get number of unique residues of this type in shell
        unique_residue_count = unique_residue_counts.get(residue, 0)

        # Append the result to the list
        results.append({
            'Interaction Type': interaction_type,
            'Residue Name': residue,
            'Approx. Percentage in Shell': approximate_percentage,
            'res atoms in a interaction/ Total atoms in Shell': f"{count} / {total_interactions_in_shell}",
            'res atoms in a interaction / all the res atoms in shell': f"{count} / {residue_specific_count}",
            'Characteristic': characteristic,
            'Atom-Specific Counts': atom_specific_counts,
            'Unique Residues in Shell': unique_residue_count  # ← NEW column added
        })

# Create a DataFrame from the results
results_df = pd.DataFrame(results)

# Aggregate results for metals and cofactors (same residue name and interaction type)
results_df_aggregated = results_df.groupby(['Interaction Type', 'Residue Name']).agg({
    'Approx. Percentage in Shell': 'first',
    'res atoms in a interaction/ Total atoms in Shell': 'first',
    'res atoms in a interaction / all the res atoms in shell': 'first',
    'Characteristic': 'first',
    'Atom-Specific Counts': lambda x: {key: sum(d.get(key, 0) for d in x) for key in set(k for d in x for k in d)},
    'Unique Residues in Shell': 'first'  # ← also carry this into the final table
}).reset_index()

# Save the results to a CSV file
output_file = '/home/parnia/parnia/F/interaction_statistics_with_residue_counts.csv'
results_df_aggregated.to_csv(output_file, index=False)
print(f"Aggregated statistics within 4 Å shell with residue-specific counts saved to {output_file}")

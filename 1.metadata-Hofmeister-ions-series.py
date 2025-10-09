#In the script,  F( Fluoride) is interested ion, which can be replaced with other ions, then extract information regarding the ion, and then download the CIF files.
import os
import requests
import pandas as pd
from Bio.PDB import MMCIFParser

# Define constants
ION_SYMBOL = "F"  # Ion of interest
SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
SUPFAM_TEMPLATE = "https://supfam.org/SUPERFAMILY/cgi-bin/gene.cgi?genome=up;seqid=%s"
CIF_DOWNLOAD_TEMPLATE = "https://files.rcsb.org/download/{}.cif"
CIF_DIRECTORY = f"/Users/respina/desktop/{ION_SYMBOL}/cif"
OUTPUT_DIRECTORY = f"/Users/respina/desktop/{ION_SYMBOL}"

# Ensure CIF directory exists
os.makedirs(CIF_DIRECTORY, exist_ok=True)

def make_request(search_params, rows=40000):
    search_params['rows'] = rows
    search_params['wt'] = 'json'
    try:
        response = requests.get(SEARCH_URL, params=search_params)
        response.raise_for_status()
        return response.json().get('response', {}).get('docs', [])
    except requests.RequestException as e:
        print(f"Request failed: {e}")
        return []

def format_search_query(filter_fields=None):
    return {
        'q': f'compound_id:{ION_SYMBOL}',
        'fl': ','.join(filter_fields) if filter_fields else None
    }

def generate_supfam_link(uniprot_id):
    if uniprot_id:
        url = SUPFAM_TEMPLATE % uniprot_id
        return f'=HYPERLINK("{url}", "{uniprot_id}")'
    return "N/A"

filter_fields = [
    'pdb_id', 'experimental_method', 'uniprot', 'title', 'resolution',
    'organism_scientific_name', 'interacting_ligands', 'assembly_type',
    'struct_asym_id', 'number_of_protein_chains', 'all_enzyme_names',
    'number_of_bound_molecules', 'molecule_type'
]

pdb_results = make_request(format_search_query(filter_fields))
print(f"Number of PDB entries found: {len(pdb_results)}")
tab1_data = []

for entry in pdb_results:
    pdb_id = entry.get('pdb_id', "").upper()
    experimental_method = entry.get('experimental_method', "")
    resolution = entry.get('resolution', None)
    molecule_type = entry.get('molecule_type', "")

    if isinstance(experimental_method, list):
        experimental_method = experimental_method[0]
    if isinstance(resolution, list):
        resolution = resolution[0]
    if isinstance(molecule_type, list):
        molecule_type = molecule_type[0]

    if isinstance(experimental_method, str) and "x-ray" in experimental_method.lower() and resolution is not None:
        try:
            if float(resolution) > 2.0:
                continue
        except ValueError:
            continue
    else:
        continue

    if isinstance(molecule_type, str) and "protein" not in molecule_type.lower():
        continue

    uniprot_id = entry.get('uniprot', [""])[0].split(':')[0]
    supfam_link = generate_supfam_link(uniprot_id)

    tab1_data.append({
        'pdb_id': pdb_id,
        'experimental_method': experimental_method,
        'uniprot': uniprot_id,
        'title': entry.get('title', ""),
        'resolution': resolution,
        'organism_scientific_name': entry.get('organism_scientific_name', ""),
        'interacting_ligands': entry.get('interacting_ligands', ""),
        'assembly_type': entry.get('assembly_type', ""),
        'struct_asym_id': entry.get('struct_asym_id', ""),
        'number_of_protein_chains': entry.get('number_of_protein_chains', ""),
        'all_enzyme_names': entry.get('all_enzyme_names', "N/A"),
        'number_of_bound_molecules': entry.get('number_of_bound_molecules', ""),
        'molecule_type': molecule_type,
        'Superfamily': supfam_link
    })

tab1_df = pd.DataFrame(tab1_data)

if tab1_df.empty:
    print("No valid PDB entries found. Exiting script.")
    exit()

tab1_df['resolution'] = pd.to_numeric(tab1_df['resolution'], errors='coerce')
filtered_df = tab1_df.loc[tab1_df.groupby('uniprot')['resolution'].idxmin()]

assembly_mapping = {
    'monomer': 1, 'dimer': 2, 'trimer': 3, 'tetramer': 4, 'pentamer': 5, 'hexamer': 6,
    'heptamer': 7, 'octamer': 8, 'nonamer': 9, 'decamer': 10, 'undecamer': 11,
    'dodecamer': 12, 'tridecamer': 13, 'tetradecamer': 14, 'pentadecamer': 15,
    'hexadecamer': 16, 'heptadecamer': 17, 'octadecamer': 18, 'nonadecamer': 19,
    'eicosamer': 20
}

def extract_numeric_assembly(value):
    if isinstance(value, list) and len(value) > 0:
        val = value[0].lower()
        if val.endswith('mer'):
            try:
                return int(val.split('-')[0])
            except:
                return assembly_mapping.get(val, None)
        return assembly_mapping.get(val, None)
    return None

filtered_df['assembly_type_numeric'] = filtered_df['assembly_type'].apply(extract_numeric_assembly)

os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)
csv_file = os.path.join(OUTPUT_DIRECTORY, f"{ION_SYMBOL}_ion_information_filtered.csv")
filtered_df.to_csv(csv_file, index=False)

def download_cif_file(pdb_id):
    cif_url = CIF_DOWNLOAD_TEMPLATE.format(pdb_id)
    cif_path = os.path.join(CIF_DIRECTORY, f"{pdb_id}.cif")
    if not os.path.exists(cif_path):
        response = requests.get(cif_url)
        if response.status_code == 200:
            with open(cif_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded: {pdb_id}.cif")
        else:
            print(f"Failed to download: {pdb_id}.cif")

for pdb_id in filtered_df['pdb_id']:
    download_cif_file(pdb_id)

def get_F_count(cif_file_path, ion_identifier=ION_SYMBOL):
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("PDB", cif_file_path)
        ion_count = sum(
            1 for model in structure
            for chain in model
            for residue in chain
            if residue.get_resname().strip() == ion_identifier
        )
        return ion_count
    except Exception as e:
        print(f"Error processing {cif_file_path}: {e}")
        return 0

def process_cif_files(cif_directory):
    results = []
    for cif_file in os.listdir(cif_directory):
        if cif_file.endswith(".cif"):
            cif_file_path = os.path.join(cif_directory, cif_file)
            f_count = get_F_count(cif_file_path)
            results.append({
                'PDB ID': cif_file.split('.')[0],
                f'{ION_SYMBOL}_count': f_count
            })
    return pd.DataFrame(results)

cif_results_df = process_cif_files(CIF_DIRECTORY)
merged_df = pd.merge(filtered_df, cif_results_df, left_on='pdb_id', right_on='PDB ID', how='outer')
total_f_count = merged_df[f'{ION_SYMBOL}_count'].sum()
total_row = pd.DataFrame([{'PDB ID': 'Total', f'{ION_SYMBOL}_count': total_f_count}])
merged_df = pd.concat([merged_df, total_row], ignore_index=True)

matched_df = merged_df[
    pd.to_numeric(merged_df['assembly_type_numeric'], errors='coerce') ==
    pd.to_numeric(merged_df['number_of_protein_chains'], errors='coerce')
]

match_count = matched_df[matched_df['PDB ID'] != 'Total'].shape[0]
total_entries = merged_df[merged_df['PDB ID'] != 'Total'].shape[0]
percentage_match = (match_count / total_entries) * 100 if total_entries > 0 else 0

print(f"\nNumber of PDB entries where 'assembly_type_numeric' equals 'number_of_protein_chains': {match_count}")
print(f"Percentage: {percentage_match:.2f}%\n")

merged_df.to_csv(csv_file, index=False)
print(f"Final results with {ION_SYMBOL} count saved to {csv_file}")
print(f"Total {ION_SYMBOL} ion count: {total_f_count}")



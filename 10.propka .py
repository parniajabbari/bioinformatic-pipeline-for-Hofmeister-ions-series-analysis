import os
import requests
import propka.run

# 🧪 Ion Symbol Setup — change this once!
ion_symbol = "GAI"   # ← Your new ion symbol (e.g., GAI, LCP, XYZ, etc.)


def fetch_pdb(pdb_id, download_dir):
    """Fetch the PDB file from the RCSB PDB repository."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    local_filename = os.path.join(download_dir, f'{pdb_id}.pdb')
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(local_filename, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded {pdb_id} to {local_filename}")
        return local_filename
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None


def replace_all_ions(pdb_filename, ion_symbol):
    """
    Replace common ion symbols (CA, MG, ZN, NA, CL, MN, FE, CU, K, etc.)
    with the user-specified ion symbol.
    """
    ion_symbols = ["CA", "MG", "ZN", "NA", "CL", "MN", "FE", "CU", "K", "CO", "NI"]
    try:
        with open(pdb_filename, 'r') as file:
            content = file.read()

        for ion in ion_symbols:
            content = content.replace(f" {ion} ", f" {ion_symbol} ")

        with open(pdb_filename, 'w') as file:
            file.write(content)

        print(f"All ions replaced with {ion_symbol} in {pdb_filename}")
    except Exception as e:
        print(f"Error replacing ions in {pdb_filename}: {e}")


def calculate_pka(pdb_filename, output_dir):
    """Run PROPKA calculation and save results in the specified directory."""
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        os.chdir(output_dir)
        propka.run.single(pdb_filename)
        print(f"pKa calculation completed for {pdb_filename}. Results saved in {output_dir}")
    except Exception as e:
        print(f"Error processing {pdb_filename}: {e}")


def process_pdb_ids(pdb_ids_file, download_dir):
    """Process each PDB ID: download, modify ions, and calculate pKa."""
    # Use ion_symbol in output folder name
    output_dir = f'/Users/respina/desktop/mmcif_files/propka-{ion_symbol}'

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    with open(pdb_ids_file, 'r') as file:
        pdb_ids = [line.strip() for line in file if line.strip()]

    for pdb_id in pdb_ids:
        pdb_filename = fetch_pdb(pdb_id, download_dir)
        if pdb_filename:
            replace_all_ions(pdb_filename, ion_symbol)
            calculate_pka(pdb_filename, output_dir)


if __name__ == "__main__":
    pdb_ids_file = '/Users/respina/desktop/mmcif_files/pdb_ids.txt'
    download_dir = '/Users/respina/desktop/mmcif_files/pdb_files'

    process_pdb_ids(pdb_ids_file, download_dir)

import os
import csv
import requests
import propka.run

# 🧪 Ion Symbol Setup — change this once!
ion_symbol = "GAI"  # Example: GAI, LCP, XYZ, etc.


# -----------------------------------------------------------------------------
# 1️⃣ Fetch PDB file
# -----------------------------------------------------------------------------
def fetch_pdb(pdb_id, download_dir):
    """Fetch the PDB file from the RCSB PDB repository."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    local_filename = os.path.join(download_dir, f'{pdb_id}.pdb')
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(local_filename, 'wb') as file:
            file.write(response.content)
        print(f"✅ Downloaded {pdb_id} → {local_filename}")
        return local_filename
    except requests.exceptions.RequestException as e:
        print(f"❌ Error downloading {pdb_id}: {e}")
        return None


# -----------------------------------------------------------------------------
# 2️⃣ Safely replace ion names (preserve PDB column format)
# -----------------------------------------------------------------------------
def replace_all_ions(pdb_filename, ion_symbol):
    """
    Replace ion atom and residue names safely without breaking PDB formatting.
    """
    ion_symbols = ["CA", "MG", "ZN", "NA", "CL", "MN", "FE", "CU", "K", "CO", "NI"]
    try:
        lines = []
        with open(pdb_filename, "r") as file:
            for line in file:
                if line.startswith(("HETATM", "ATOM")):
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    if res_name in ion_symbols or atom_name in ion_symbols:
                        # Replace both atom and residue name, keeping proper padding
                        new_line = (
                            line[:12]
                            + f"{ion_symbol:>4}"  # atom name
                            + line[16:17]
                            + f"{ion_symbol:>3}"  # residue name
                            + line[20:]
                        )
                        line = new_line
                lines.append(line)

        with open(pdb_filename, "w") as f:
            f.writelines(lines)

        print(f"🔁 All ions safely replaced with {ion_symbol} in {pdb_filename}")
    except Exception as e:
        print(f"❌ Error replacing ions in {pdb_filename}: {e}")


# -----------------------------------------------------------------------------
# 3️⃣ Run PROPKA
# -----------------------------------------------------------------------------
def calculate_pka(pdb_filename, output_dir):
    """Run PROPKA calculation and save results in the specified directory."""
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        pdb_base = os.path.splitext(os.path.basename(pdb_filename))[0]
        result_file = os.path.join(output_dir, f"{pdb_base}.pka")

        # Skip if result already exists
        if os.path.exists(result_file):
            print(f"⏩ Skipping {pdb_base} (PROPKA already done)")
            return True

        # Run PROPKA
        os.chdir(output_dir)
        propka.run.single(pdb_filename)
        print(f"🧮 pKa calculation completed for {pdb_filename}. Results saved in {output_dir}")
        return True
    except Exception as e:
        print(f"❌ Error processing {pdb_filename}: {e}")
        return False


# -----------------------------------------------------------------------------
# 4️⃣ Extract PROPKA results
# -----------------------------------------------------------------------------
def extract_propka_results(pdb_filename, output_dir):
    """Extracts and returns PROPKA results from the generated .pka file."""
    pdb_base = os.path.splitext(os.path.basename(pdb_filename))[0]
    propka_file = os.path.join(output_dir, f"{pdb_base}.pka")

    if not os.path.exists(propka_file):
        print(f"⚠️ PROPKA result file not found for {pdb_base}")
        return []

    results = []
    with open(propka_file, "r") as f:
        for line in f:
            # Example line format: " ASP  25 A     3.57        4.00        -0.43"
            if line.strip().startswith(("ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG")):
                parts = line.split()
                if len(parts) >= 5:
                    residue = parts[0]
                    res_num = parts[1]
                    chain = parts[2]
                    pka_value = parts[3]
                    results.append({
                        "PDB_ID": pdb_base,
                        "Residue": residue,
                        "Res_Num": res_num,
                        "Chain": chain,
                        "pKa": pka_value
                    })

    # Print summary
    if results:
        print(f"📊 pKa values for {pdb_base}:")
        for r in results:
            print(f"  {r['Residue']} {r['Res_Num']}{r['Chain']} → pKa = {r['pKa']}")
    else:
        print(f"⚠️ No titratable residues found in {pdb_base}")

    return results


# -----------------------------------------------------------------------------
# 5️⃣ Read CSV of PDB IDs
# -----------------------------------------------------------------------------
def read_pdb_ids_from_csv(csv_file):
    """Read PDB IDs from the first column of the provided CSV file."""
    pdb_ids = []
    try:
        with open(csv_file, 'r', newline='') as file:
            reader = csv.reader(file)
            next(reader, None)  # Skip header if exists
            for row in reader:
                if row and row[0].strip():
                    pdb_ids.append(row[0].strip().upper())
        print(f"📄 Found {len(pdb_ids)} PDB IDs in {csv_file}")
    except Exception as e:
        print(f"❌ Error reading {csv_file}: {e}")
    return pdb_ids


# -----------------------------------------------------------------------------
# 6️⃣ Main pipeline
# -----------------------------------------------------------------------------
def process_pdb_ids(csv_file, download_dir):
    """Process each PDB ID from the CSV file: download, modify ions, run PROPKA, extract results."""
    output_dir = f'/Users/respina/desktop/{ion_symbol}/propka-{ion_symbol}'
    summary_file = os.path.join(output_dir, "propka_summary.csv")
    all_results = []

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    pdb_ids = read_pdb_ids_from_csv(csv_file)

    for pdb_id in pdb_ids:
        pdb_filename = fetch_pdb(pdb_id, download_dir)
        if not pdb_filename:
            continue

        replace_all_ions(pdb_filename, ion_symbol)

        success = calculate_pka(pdb_filename, output_dir)
        if success:
            results = extract_propka_results(pdb_filename, output_dir)
            all_results.extend(results)

    # 🧾 Save all pKa results to a single CSV file
    if all_results:
        with open(summary_file, "w", newline='') as csvfile:
            fieldnames = ["PDB_ID", "Residue", "Res_Num", "Chain", "pKa"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\n💾 All PROPKA results saved to: {summary_file}")
    else:
        print("\n⚠️ No PROPKA results to save.")


# -----------------------------------------------------------------------------
# 🚀 Entry point
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    csv_file = f'/Users/respina/desktop/{ion_symbol}/{ion_symbol_ion_information_filtered.csv'
    download_dir = f'/Users/respina/desktop/{ion_symbol}/pdb_files'

    process_pdb_ids(csv_file, download_dir)

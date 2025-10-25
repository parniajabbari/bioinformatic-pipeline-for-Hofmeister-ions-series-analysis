import os
import csv
import requests
import propka.run
from datetime import datetime
import argparse

# -----------------------------------------------------------------------------
# --- ARGUMENT PARSER ---
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Download PDBs, replace ions, run PROPKA, and extract pKa values.")
parser.add_argument("ion", help="Ion symbol, e.g., GAI, LCP")
parser.add_argument("--output", default=None, help="Output directory for all results (optional)")
args = parser.parse_args()

ION_SYMBOL = args.ion
OUTPUT_DIR = args.output  # Can be None; we will define default later

# -----------------------------------------------------------------------------
# --- LOGGING ---
# -----------------------------------------------------------------------------
def setup_logging(base_dir):
    """Create a log file to record all activity."""
    log_dir = os.path.join(base_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{ION_SYMBOL}_propka_log.txt")
    return log_file

def log_message(message, log_file):
    """Append a timestamped message to the log file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}\n"
    print(message)
    with open(log_file, "a") as f:
        f.write(line)

# -----------------------------------------------------------------------------
# --- PDB DOWNLOAD ---
# -----------------------------------------------------------------------------
def fetch_pdb(pdb_id, download_dir, log_file):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    local_filename = os.path.join(download_dir, f'{pdb_id}.pdb')
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(local_filename, 'wb') as file:
            file.write(response.content)
        log_message(f"✅ Downloaded {pdb_id} → {local_filename}", log_file)
        return local_filename
    except requests.exceptions.RequestException as e:
        log_message(f"❌ Error downloading {pdb_id}: {e}", log_file)
        return None

# -----------------------------------------------------------------------------
# --- REPLACE ION SYMBOLS ---
# -----------------------------------------------------------------------------
def replace_all_ions(pdb_filename, ion_symbol, log_file):
    ion_symbols = ["CA", "MG", "ZN", "NA", "CL", "MN", "FE", "CU", "K", "CO", "NI"]
    try:
        lines = []
        with open(pdb_filename, "r") as file:
            for line in file:
                if line.startswith(("HETATM", "ATOM")):
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    if res_name in ion_symbols or atom_name in ion_symbols:
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

        log_message(f"🔁 All ions safely replaced with {ion_symbol} in {pdb_filename}", log_file)
    except Exception as e:
        log_message(f"❌ Error replacing ions in {pdb_filename}: {e}", log_file)

# -----------------------------------------------------------------------------
# --- RUN PROPKA ---
# -----------------------------------------------------------------------------
def calculate_pka(pdb_filename, output_dir, log_file):
    try:
        os.makedirs(output_dir, exist_ok=True)
        pdb_base = os.path.splitext(os.path.basename(pdb_filename))[0]
        result_file = os.path.join(output_dir, f"{pdb_base}.pka")

        if os.path.exists(result_file):
            log_message(f"⏩ Skipping {pdb_base} (PROPKA already done)", log_file)
            return True

        os.chdir(output_dir)
        propka.run.single(pdb_filename)
        log_message(f"🧮 pKa calculation completed for {pdb_filename}. Results saved in {output_dir}", log_file)
        return True
    except Exception as e:
        log_message(f"❌ Error processing {pdb_filename}: {e}", log_file)
        return False

# -----------------------------------------------------------------------------
# --- EXTRACT PROPKA RESULTS ---
# -----------------------------------------------------------------------------
def extract_propka_results(pdb_filename, output_dir, log_file):
    pdb_base = os.path.splitext(os.path.basename(pdb_filename))[0]
    propka_file = os.path.join(output_dir, f"{pdb_base}.pka")

    if not os.path.exists(propka_file):
        log_message(f"⚠️ PROPKA result file not found for {pdb_base}", log_file)
        return []

    results = []
    with open(propka_file, "r") as f:
        for line in f:
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

    if results:
        log_message(f"📊 Extracted {len(results)} pKa entries for {pdb_base}", log_file)
    else:
        log_message(f"⚠️ No titratable residues found in {pdb_base}", log_file)
    return results

# -----------------------------------------------------------------------------
# --- READ PDB IDS FROM CSV ---
# -----------------------------------------------------------------------------
def read_pdb_ids_from_csv(csv_file, log_file):
    pdb_ids = []
    try:
        with open(csv_file, 'r', newline='') as file:
            reader = csv.reader(file)
            next(reader, None)
            for row in reader:
                if row and row[0].strip():
                    pdb_ids.append(row[0].strip().upper())
        log_message(f"📄 Found {len(pdb_ids)} PDB IDs in {csv_file}", log_file)
    except Exception as e:
        log_message(f"❌ Error reading {csv_file}: {e}", log_file)
    return pdb_ids

# -----------------------------------------------------------------------------
# --- MAIN PIPELINE ---
# -----------------------------------------------------------------------------
def process_pdb_ids(csv_file, download_dir, output_dir):
    os.makedirs(download_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    log_file = setup_logging(output_dir)

    all_results = []
    pdb_ids = read_pdb_ids_from_csv(csv_file, log_file)

    for pdb_id in pdb_ids:
        log_message(f"\n🚀 Starting process for {pdb_id}", log_file)
        pdb_filename = fetch_pdb(pdb_id, download_dir, log_file)
        if not pdb_filename:
            continue

        replace_all_ions(pdb_filename, ION_SYMBOL, log_file)

        success = calculate_pka(pdb_filename, output_dir, log_file)
        if success:
            results = extract_propka_results(pdb_filename, output_dir, log_file)
            all_results.extend(results)

    summary_file = os.path.join(output_dir, "propka_summary.csv")
    if all_results:
        with open(summary_file, "w", newline='') as csvfile:
            fieldnames = ["PDB_ID", "Residue", "Res_Num", "Chain", "pKa"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        log_message(f"\n💾 All PROPKA results saved to: {summary_file}", log_file)
    else:
        log_message("\n⚠️ No PROPKA results to save.", log_file)

# -----------------------------------------------------------------------------
# --- ENTRY POINT ---
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Default directories if --output is not provided
    base_dir = os.path.join("/Users/respina/desktop", ION_SYMBOL)
    if OUTPUT_DIR is None:
        OUTPUT_DIR = base_dir

    csv_file = os.path.join(base_dir, f"{ION_SYMBOL}_ion_information_filtered.csv")
    download_dir = os.path.join(OUTPUT_DIR, "pdb_files")

    process_pdb_ids(csv_file, download_dir, OUTPUT_DIR)

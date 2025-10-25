import os
import subprocess
import argparse

# -----------------------------
# 1. Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Convert CIF files to DSSP format for a given Hofmeister ion."
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
BASE_DIR = args.output or os.path.join(BASE_DESKTOP, ION_SYMBOL)

input_dir = os.path.join(BASE_DIR, "cif")
output_dir = os.path.join(BASE_DIR, "dssp_results")

os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# 3. Conversion function
# -----------------------------
def convert_cif_to_dssp(cif_file, output_file):
    try:
        subprocess.run(['mkdssp', '-i', cif_file, '-o', output_file], check=True)
        print(f"✅ Converted {os.path.basename(cif_file)} → {os.path.basename(output_file)}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error converting {os.path.basename(cif_file)}: {e}")

# -----------------------------
# 4. Main conversion loop
# -----------------------------
def main():
    if not os.path.exists(input_dir):
        print(f"⚠️ CIF directory not found: {input_dir}")
        return

    cif_files = [f for f in os.listdir(input_dir) if f.endswith(".cif")]

    if not cif_files:
        print(f"⚠️ No CIF files found in {input_dir}")
        return

    for cif_file in cif_files:
        cif_path = os.path.join(input_dir, cif_file)
        dssp_path = os.path.join(output_dir, f"{os.path.splitext(cif_file)[0]}.dssp")
        convert_cif_to_dssp(cif_path, dssp_path)

    print(f"\n📁 All conversions complete. DSSP files saved in: {output_dir}")

if __name__ == "__main__":
    main()


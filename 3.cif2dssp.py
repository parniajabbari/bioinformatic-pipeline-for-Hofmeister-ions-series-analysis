#after downloading the cif or Pdb can convert them to dssp 

import os
import subprocess


# Directory containing .cif files
input_dir = "/Users/respina/desktop/GAI/cif"

# Directory to store the output .dssp files
output_dir = "/Users/respina/desktop/GAI/dssp_results"


# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to convert .cif to .dssp using mkdssp
def convert_cif_to_dssp(cif_file, output_file):
    try:
        # Run the mkdssp command
        subprocess.run(['mkdssp', '-i', cif_file, '-o', output_file], check=True)
        print(f"Converted {os.path.basename(cif_file)} to {os.path.basename(output_file)}")
    except subprocess.CalledProcessError as e:
        print(f"Error converting {os.path.basename(cif_file)}: {e}")

# Iterate through all .cif files in the input directory
for cif_file in os.listdir(input_dir):
    if cif_file.endswith(".cif"):
        # Full path of the input .cif file
        cif_file_path = os.path.join(input_dir, cif_file)
        
        # Define the output .dssp file path
        dssp_file_path = os.path.join(output_dir, f"{os.path.splitext(cif_file)[0]}.dssp")
        
        # Convert .cif to .dssp
        convert_cif_to_dssp(cif_file_path, dssp_file_path)

print("Conversion process completed.")


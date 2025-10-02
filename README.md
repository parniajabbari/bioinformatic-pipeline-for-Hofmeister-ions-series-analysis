# bioinformatic-pipeline-for-Hofmeister-ions-series-analysis


The Hofmeister series refers to a classification of ions that influence the stability and structure of proteins and other macromolecules in aqueous solutions. These ions are categorized into kosmotropes, which stabilize proteins by enhancing the structure of water, and chaotropes, which destabilize proteins by disrupting the water structure. This series helps us understand how ions affect protein folding, solubility, and other important biological processes that are fundamental concepts in protein chemistry. For this purpose, the bioinformatic pipelines have designed to automate the extraction of meta-data from diverse datasets (PDB, PDBsum, and Uniport). For the computational analyses in pipelines extensively Python prgramming has utilized. 

<img width="738" height="303" alt="image" src="https://github.com/user-attachments/assets/a7d08e81-66d9-497d-b1f4-57e0a6146b69" />

Workflow for Metadata Extraction and Structural Analysis
To systematically analyze protein-ion interactions involving Hofmeister series ions, we developed and executed an automated workflow within a High-Performance Computing (HPC) environment, using Bash scripting as well as Python programming.


Ion–Protein Interaction Analysis Pipeline This pipeline helps you study how ions (default: F■) interact with proteins using PDB CIF files. You can
replace F with other ions (Cl, Br, I, etc.) easily.
■ Requirements
Before starting, make sure you have installed: - Python 3.8+ - Packages: pip install biopython pandas
matplotlib requests numpy - External tools: mkdssp (for secondary structure analysis)

1■■ Extract Metadata & Download CIF Files
Run: python metadata-Hofmeister-ions-series.py
- Searches PDB for entries with your ion of interest (ION_SYMBOL in script). - Filters for high-resolution protein X-ray structures (≤2.0 Å).
- Downloads CIF files into /ION/cif.
- Creates ION_ion_information_filtered.csv with metadata + ion counts.
  
2■■ Analyze Ion Interactions
Run: python result.py
- Reads CIFs, finds all ion–atom interactions (distance, angle).
- Classifies interactions: Ion–Ion, Salt Bridge, Metal, H-bond, Aliphatic, Aromatic, Cofactor/Ligand, Water.
- Saves detailed results into result.csv.
  
3■■ Convert CIF → DSSP (Secondary Structure)
Run: python cif2dssp.py
- Uses mkdssp to convert CIF → DSSP. - Stores .dssp files in /dssp_results.
  
4■■ Add Secondary Structure Info
Run: python DSSP.column.py
- Merges DSSP secondary structure info with interaction data in result.csv. - Assigns each residue to: α-helix, β-strand, turn, coil, etc.
  
5■■ Generate Interaction Tables
Run: python table.py
- Processes result.csv.
- Summarizes residue types, interaction counts within interested shell, unique residue counts. - Output: interaction_statistics_with_residue_counts.csv.
  
6■■ Plot Interaction Type Distribution (Pie Chart)
Run: python pichart-interaction-type.py
- Creates pie chart of interaction types.
- Saves as interaction_type_distribution.png.
  
7■■ Plot Secondary Structure Distribution (Pie Chart)
Run: python pichart-dssp.py
- Visualizes secondary structure preference of ion interactions. - Saves as secondary_structure_distribution.png.
  
8■■ Analyze Hydration Shell
Run: python water.py
- Counts nearby water molecules (≤4 Å) around each ion. - Produces histogram of hydration shell sizes.








# bioinformatic-pipeline-for-Hofmeister-ions-series-analysis


The Hofmeister series refers to a classification of ions that influence the stability and structure of proteins and other macromolecules in aqueous solutions. These ions are categorized into kosmotropes, which stabilize proteins by enhancing the structure of water, and chaotropes, which destabilize proteins by disrupting the water structure. This series helps us understand how ions affect protein folding, solubility, and other important biological processes that are fundamental concepts in protein chemistry. For this purpose, the bioinformatic pipelines have designed to automate the extraction of meta-data from diverse datasets (PDB, PDBsum, and Uniport). For the computational analyses in pipelines extensively Python prgramming has utilized. 

<img width="738" height="303" alt="image" src="https://github.com/user-attachments/assets/a7d08e81-66d9-497d-b1f4-57e0a6146b69" />

Workflow for Metadata Extraction and Structural Analysis
To systematically analyze protein-ion interactions involving Hofmeister series ions, we developed and executed an automated workflow within a High-Performance Computing (HPC) environment, using Bash scripting as well as Python programming.


This pipeline helps you study how ions (default: F) interact with proteins using PDB CIF files. You can
replace F with other ions (Cl, Br, I, etc.) easily.

■ Requirements
Before starting, make sure you have installed: - Python 3.8+ - Packages: pip install biopython pandas
matplotlib requests numpy - External tools: mkdssp (for secondary structure analysis)

1. Extract Metadata & Download CIF Files
Run: python metadata-Hofmeister-ions-series.py
- Searches PDB for entries with your ion of interest (ION_SYMBOL in script). - Filters for high-resolution protein X-ray structures (≤2.0 Å).
- Downloads CIF files into /ION/cif.
- Creates ION_ion_information_filtered.csv with metadata + ion counts.

- The information are: 

pdb_id , experimental_method , UniProt , title , Resolution , organism_scientific_name , interacting_ligands , assembly_type  , struct_asym_id , number_of_protein_chains , all_enzyme_names , number_of_bound_molecules , SUPFAM , number of chains , lON_count , Sum ION

  
2. Analyze Ion Interactions
Run: python result.py
- Reads CIFs, finds all ion–atom interactions (distance, angle).
- Classifies interactions: Ion–Ion, Salt Bridge, Metal, H-bond, Aliphatic, Aromatic, Cofactor/Ligand, Water.
- Saves detailed results into result.csv.
  
3. Convert CIF → DSSP (Secondary Structure)
Run: python cif2dssp.py
- Uses mkdssp to convert CIF → DSSP. - Stores .dssp files in /dssp_results.
  
4. Add Secondary Structure Info
Run: python DSSP.column.py
- Merges DSSP secondary structure info with interaction data in result.csv. - Assigns each residue to: α-helix, β-strand, turn, coil, etc.
  
5. Generate Interaction Tables
Run: python table.py
- Processes result.csv.
- Summarizes residue types, interaction counts within interested shell, unique residue counts. - Output: interaction_statistics_with_residue_counts.csv.

- the columns in your output CSV based on the following photo:
Interaction Type – Type of interaction (H-bond, Salt Bridge, Metal, Ion, Cofactor/Ligand).
Residue Name – Three-letter code of the interacting residue or ion.
Approx. Percentage in Shell – % of all interactions in the 4 Å shell involving this residue.
Atoms in Interaction / Total Atoms in Shell – Atom-level contribution of this residue relative to all atoms.
Atoms in Interaction / Total Residue Atoms in Shell – Atom-level contribution relative to this residue type.
Characteristic – Chemical/functional class of the residue (Hydrophobic, Polar, Metal, Ion, etc.).
Atom-Specific Counts – Counts of each residue atom involved in interactions.
Unique Residues in Shell – Number of unique residue instances of this type in the shell.

  ![WhatsApp Image 2024-11-04 at 16 29 35](https://github.com/user-attachments/assets/b13d407e-f580-4f84-ab1e-29f31cd82bab)

6. Plot Interaction Type Distribution (Pie Chart)
Run: python pichart-interaction-type.py
- Creates pie chart of interaction types.
- Saves as interaction_type_distribution.png.

- <img width="328" height="265" alt="image" src="https://github.com/user-attachments/assets/637805b7-8760-4462-adcc-ae47d668dcc5" />

7. Plot Secondary Structure Distribution (Pie Chart)
Run: python pichart-dssp.py
- Visualizes secondary structure preference of ion interactions. - Saves as secondary_structure_distribution.png.
- <img width="327" height="263" alt="image" src="https://github.com/user-attachments/assets/bdfabb81-5a05-426a-8bde-0e6581c7eb02" />

8. Analyze Hydration Shell
Run: python water.py
- Counts nearby water molecules (≤4 or 5 Å) around each ion. - Produces histogram of hydration shell sizes.

- <img width="490" height="294" alt="image" src="https://github.com/user-attachments/assets/47fbf0f9-6157-4ddc-bbcc-c8d752cb0ae6" />


9. Can also in detail look at the frequency of the each residues in the entries that contain the ion 

<img width="558" height="295" alt="image" src="https://github.com/user-attachments/assets/ff5db81f-b1ac-4e68-93d6-af7fe84dca64" />


10. PROPKA Calculation of pKa Values
Addionally, after ion replacement, the script runs PROPKA to predict pKa values of ionizable residues in the modified protein structures.This provides insight into how different ions affect local protein charge and stability, a key factor in Hofmeister effects.Results are stored in a folder automatically named after the ion_symbol (propka-GAI), keeping outputs organized for multiple ions.








# bioinformatic-pipeline-for-Hofmeister-ions-series-analysis


Executive Summary
The Hofmeister series refers to a classification of ions that influence the stability and structure of proteins and other macromolecules in aqueous solutions. These ions are categorized into kosmotropes, which stabilize proteins by enhancing the structure of water, and chaotropes, which destabilize proteins by disrupting the water structure. This series helps us understand how ions affect protein folding, solubility, and other important biological processes that are fundamental concepts in protein chemistry. For this purpose, the bioinformatic pipelines have designed to automate the extraction of meta-data from diverse datasets (PDB, PDBsum, and Uniport). For the computational analyses in pipelines extensively Python prgramming has utilized. 

<img width="738" height="303" alt="image" src="https://github.com/user-attachments/assets/a7d08e81-66d9-497d-b1f4-57e0a6146b69" />


Methodology 

1.	Data Acquisition and Filtering

1.1	Identification of PDB Entries containing the ions
All structural entries containing Hofmeister ions as standalone ligands were retrieved from the Protein Data Bank in Europe (PDBe). Ligands were considered "standalone" if they were not covalently bound or part of a larger complex but instead existed as discrete ions interacting non-covalently within the structure. The Protein Data Bank in Europe (PDBe) is a member of the Worldwide Protein Data Bank (wwPDB) that provides access to curated macromolecular structure data, including proteins, nucleic acids, and complexes [1].
1.2	Data Filtering and Quality Control
To ensure that only high-quality and relevant structural data were included in the analysis, several filtering criteria were applied as follows:
1.2.1	Experiment Method: Only structures determined using X-ray crystallography were considered. Structures resolved through NMR spectroscopy, cryo-electron microscopy, powder diffraction, and neutron diffraction were excluded. 
1.2.2	Resolution Threshold: Structures were required to have a crystallographic resolution of 2.0 Å or better (≤ 2.0 Å). High-resolution structures increase the reliability of ion position, reduce coordinate uncertainty, and allow for better interpretation of ligand-protein interactions.
1.2.3	Molecular Type: Only structures classified as proteins were included. Entries containing nucleic acids (DNA/RNA), either alone or in combination with proteins, were excluded to maintain a focused scope on protein-ion interactions.
1.2.4	UniProt ID Filtering: To reduce redundancy and avoid overrepresentation of specific proteins, entries sharing the same UniProt accession number were compared. From each group of entries corresponding to a single UniProt ID, only the structure with the highest resolution was retained for analysis. 


<img width="484" height="238" alt="image" src="https://github.com/user-attachments/assets/71b9b53f-1331-46ab-82e6-877017ef8235" />


<img width="474" height="210" alt="image" src="https://github.com/user-attachments/assets/95c9d379-0ee3-45ef-bd5c-5818e37bad99" />


2.	Workflow for Metadata Extraction and Structural Analysis
To systematically analyze protein-ion interactions involving Hofmeister series ions, we developed and executed an automated workflow within a High-Performance Computing (HPC) environment, using Bash scripting as well as Python programming.

2.1 Extraction of Metadata and Ion Information
Metadata for each PDB entry was extracted using of PDBe API, PDBsum, UniProt, and Subfamily via the requests library, then the pandas library was used for Data manipulation and analysis, also different libraries of Biopython were used. This automated process collected the following information for each structure and then stored it in CSV format.
•	PDB ID
•	Experimental Method
•	Resolution
•	Title
•	Organism Scientific Name
•	UniProt Accession Number
•	Assembly Type and Struct_Asym_ID
•	Number of Protein Chains
•	Number of Total Chains
•	Number of Bound Molecules
•	Interacting Ligands
•	All Enzyme Names
•	Superfamily (SUPFAM) Classification
•	Ion Count and Total Ion Sum






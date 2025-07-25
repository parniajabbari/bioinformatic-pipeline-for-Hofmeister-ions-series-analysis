import os
from Bio.PDB import MMCIFParser
import numpy as np
import csv

def calculate_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def calculate_angle(a, b, c):
    ab = np.array(b) - np.array(a)
    bc = np.array(c) - np.array(b)
    cos_angle = np.dot(ab, bc) / (np.linalg.norm(ab) * np.linalg.norm(bc))
    angle = np.degrees(np.arccos(cos_angle))
    return angle

def classify_F_interaction(dist, interaction_type):
    if dist < 4 and interaction_type == "N/A":
        return "Cofactor/ligand"
    return interaction_type

def main():
    cif_dir = "/home/parnia/parnia/Ca/cif"  # Path remains the same
    output_csv_path = "/home/parnia/parnia/Ca/result.csv"

    cif_files = [os.path.join(cif_dir, f) for f in os.listdir(cif_dir) if f.endswith('.cif')]

    with open(output_csv_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "PDB ID", "F Full Name", "F Atom Name", "Residue Atom Name",
            "Residue Name", "Residue ID", "Distance (\u00c5)",
            "Angle (\u00b0)", "Interaction Type",
            "F Atom Coordinates", "Residue Atom Coordinates", "Chain ID"
        ])

        for file_path in cif_files:
            pdb_id = os.path.splitext(os.path.basename(file_path))[0]
            parser = MMCIFParser()
            try:
                structure = parser.get_structure(pdb_id, file_path)
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")
                continue

            F_atoms = []
            other_atoms = []

            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    for residue in chain:
                        for atom in residue:
                            if atom.element == 'F':
                                F_atoms.append({
                                    'name': atom.get_name(),
                                    'coord': atom.coord,
                                    'chain': chain_id,
                                    'full_id': f"('{residue.resname}', {residue.id[1]}, '{chain_id}')"
                                })
                            else:
                                other_atoms.append({
                                    'name': atom.get_name(),
                                    'coord': atom.coord,
                                    'resname': residue.resname,
                                    'res_id': residue.id[1],
                                    'chain': chain_id
                                })

            results = []
            # F–F interaction check
            for i, F1 in enumerate(F_atoms):
                for F2 in F_atoms[i + 1:]:
                    if F1['full_id'] != F2['full_id']:
                        dist = calculate_distance(F1['coord'], F2['coord'])
                        if dist < 4:
                            interaction_type = "Ion"
                            results.append([
                                F1['full_id'], F1['name'], F2['name'], "F", "N/A",
                                f"{dist:.2f}", "N/A", interaction_type,
                                str(F1['coord'].tolist()), str(F2['coord'].tolist()), F1['chain']
                            ])

            # F interactions with non-F atoms
            for F in F_atoms:
                for other in other_atoms:
                    dist = calculate_distance(F['coord'], other['coord'])
                    interaction_type = "N/A"
                    angle = None

                    if dist > 8:
                        continue

                    if other['resname'] == 'PEG':
                        interaction_type = classify_F_interaction(dist, interaction_type)

                    hofmeister_ions = {
                        'CA': 1.00, 'MG': 0.72, 'NA': 1.02, 'K': 1.38, 'SR': 1.24, 'BA': 1.43,
                        'CL': 1.67, 'BR': 1.82, 'I': 2.06, 'SO4': 2.58, 'NO3': 1.90
                    }
                    ion_radius = hofmeister_ions.get(other['resname'], 0)
                    if ion_radius > 0 and dist <= 4:
                        interaction_type = "Ion"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        metals = ['ZN', 'MG', 'FE', 'CA', 'CU', 'NI', 'MN', 'CO', 'NA', 'K', 'CD', 'AL', 'AG', 'AU', 'PT', 'LI', 'SR', 'BA', 'RB']
                        if other['resname'] in metals:
                            interaction_type = "Interaction with Metal"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        if other['resname'] in ['ASN', 'LEU', 'ALA', 'CYS', 'GLN', 'GLY', 'ILE', 'MET', 'PRO', 'SER', 'THR', 'VAL', 'OH']:
                            interaction_type = "Aliphatic"

                    if interaction_type == "N/A" and 3.5 <= dist <= 4:
                        if other['resname'] in ['PHE', 'TYR', 'TRP']:
                            interaction_type = "Aromatic"

                    if interaction_type == "N/A":
                        ion_charge_dict = {
                            'F': -1  # Fluoride ion is typically -1
                        }
                        ion_overall_charge = ion_charge_dict.get('F', None)
                        if ion_overall_charge is not None:
                            specific_neg_atoms = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
                            specific_pos_atoms = {'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ']}

                            resname = other['resname']
                            atom_name = other['name']

                            if ion_overall_charge > 0:
                                if resname in specific_neg_atoms and atom_name in specific_neg_atoms[resname]:
                                    if 2.2 <= dist <= 4.0:
                                        interaction_type = "Salt Bridge"
                                elif 3.5 <= dist <= 4.0:
                                    interaction_type = "Salt Bridge"
                            elif ion_overall_charge < 0:
                                if resname in specific_pos_atoms and atom_name in specific_pos_atoms[resname]:
                                    if 2.2 <= dist <= 4.0:
                                        interaction_type = "Salt Bridge"
                                elif 3.5 <= dist <= 4.0:
                                    interaction_type = "Salt Bridge"

                    if interaction_type == "N/A" and 2.2 <= dist <= 3.5:
                        for other2 in other_atoms:
                            if other['coord'].tolist() != other2['coord'].tolist() or other['name'] != other2['name']:
                                angle = calculate_angle(F['coord'], other['coord'], other2['coord'])
                                if angle < 135:
                                    if other['resname'] == "HOH":
                                        interaction_type = "H-bond with Water"
                                    elif other['name'] in ['C', 'N', 'O', 'CA']:
                                        interaction_type = "H-bond with Main Chain"
                                    else:
                                        interaction_type = "H-bond with Side Chain"

                    results.append([
                        F['full_id'], F['name'], other['name'],
                        other['resname'], other['res_id'],
                        f"{dist:.2f}", angle if angle else "N/A",
                        interaction_type, str(F['coord'].tolist()), str(other['coord'].tolist()), F['chain']
                    ])

            for row in results:
                writer.writerow([pdb_id] + row)
            print(f"Processed {pdb_id} with {len(results)} interactions.")

    all_rows = []
    unique_f_set = set()
    with open(output_csv_path, mode='r') as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
        for row in reader:
            try:
                distance = float(row["Distance (\u00c5)"])
                if distance <= 5:
                    unique_f_set.add(row["F Full Name"])
            except ValueError:
                pass
            all_rows.append(row)

    new_fieldname = "Unique F in shell (≤ 4 Å)"
    updated_fieldnames = fieldnames + [new_fieldname]

    with open(output_csv_path, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=updated_fieldnames)
        writer.writeheader()
        for row in all_rows:
            row[new_fieldname] = len(unique_f_set)
            writer.writerow(row)

    print(f"Unique F in shell (≤ 4 Å): {len(unique_f_set)}")

if __name__ == "__main__":
    main()



#!/usr/bin/env python

import argparse
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

def get_charge_and_multiplicity(smiles):
    """
    Predict the charge and multiplicity of a molecule given its SMILES string.
    NOTE: This is an RDKit-based approximation. For complex systems or 
          accurate spin states, use a quantum chemical approach.

    Parameters:
        smiles (str): SMILES representation of the molecule.

    Returns:
        tuple: (charge, multiplicity, unpaired_electrons)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Add hydrogens so RDKit correctly assigns radical electrons
    mol = Chem.AddHs(mol)

    # 1) Calculate formal charge
    charge = Chem.GetFormalCharge(mol)

    # 2) Count the total unpaired electrons (radicals)
    unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())

    # 3) Determine multiplicity = N + 1, where N is number of unpaired electrons
    multiplicity = unpaired_electrons + 1 if unpaired_electrons > 0 else 1

    return charge, multiplicity, unpaired_electrons

def generate_3d_structure(smiles):
    """
    Converts a SMILES string to a 3D molecular structure using RDKit.
    """
    print(f"Generating 3D structure for SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol = Chem.AddHs(mol)
    # Attempt 3D embedding
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        raise ValueError(f"Could not generate 3D conformer for SMILES: {smiles}")

    return mol

def conformer_search(mol, num_conformers=10):
    """
    Performs a conformer search using MMFF94 and returns the molecule plus 
    the conformer ID of the lowest energy conformer.
    """
    print(f"Performing conformer search for molecule with up to {num_conformers} conformers using MMFF94.")
    conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=42)
    energies = []

    for conf_id in conformer_ids:
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        if mmff_props is None:
            print(f"MMFF94 force field could not be created for conformer {conf_id}. Skipping.")
            continue

        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
        if ff is None:
            print(f"MMFF94 optimization failed for conformer {conf_id}. Skipping.")
            continue

        ff.Minimize()
        energy = ff.CalcEnergy()
        energies.append((conf_id, energy))

    if not energies:
        raise ValueError("No valid conformers found with MMFF94. Check the input molecule.")

    # Select the lowest energy conformer
    lowest_conf_id, lowest_energy = min(energies, key=lambda x: x[1])
    print(f"Lowest energy conformer ID: {lowest_conf_id}, Energy: {lowest_energy}")
    return mol, lowest_conf_id

def create_gaussian_input(molecule, charge, spin_multiplicity, filename):
    """
    Generates a Gaussian input file for a given molecule with the specified
    charge and spin multiplicity.
    """
    print(f"Creating Gaussian input file: {filename}")
    route_parameters = {
        "opt": "",
        "freq": "NoRaman",
        #"EmpiricalDispersion": "GD3BJ",
    }

    chk_file = filename.replace('.gjf', '.chk')

    gaussian_input = GaussianInput(
        molecule,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        functional="M062X",
        basis_set="6-31++G(d,p)",
        route_parameters=route_parameters,
        link0_parameters={
            "%mem": "64GB",
            "%nprocshared": "16",
            "%Chk": chk_file
        },
        title="Generated Gaussian Input"
    )

    if os.path.exists(filename):
        print(f"Warning: {filename} already exists and will be overwritten.")

    gaussian_input.write_file(filename)
    print(f"Gaussian input file {filename} created successfully.\n")

def process_smiles(smiles, structure_label):
    """
    Given a SMILES string and a structure label (e.g., 'oxidised' or 'reduced'),
    this function:
      1) Validates the SMILES string.
      2) Predicts charge and multiplicity.
      3) Generates a 3D structure and optimizes conformers.
      4) Creates a Pymatgen Molecule from the best conformer.
    Returns: (pmg_molecule, charge, multiplicity)
    """
    if not isinstance(smiles, str) or pd.isna(smiles):
        raise ValueError(f"Invalid or missing SMILES for {structure_label} structure: {smiles}")

    # 1) Predict charge and multiplicity
    charge, multiplicity, _ = get_charge_and_multiplicity(smiles)

    # 2) Generate 3D structure
    rdkit_mol = generate_3d_structure(smiles)
    rdkit_mol, lowest_conf_id = conformer_search(rdkit_mol)

    # 3) Convert RDKit molecule to Pymatgen Molecule
    atoms = [atom.GetSymbol() for atom in rdkit_mol.GetAtoms()]
    coords = rdkit_mol.GetConformer(lowest_conf_id).GetPositions()
    pmg_mol = Molecule(atoms, coords)

    return pmg_mol, charge, multiplicity

def main():
    parser = argparse.ArgumentParser(description="Generate 3D structures and Gaussian input files from SMILES.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the CSV file.')
    args = parser.parse_args()

    file_path = args.file

    print(f"Loading CSV file: {file_path}")
    df = pd.read_csv(file_path)
    print(f"CSV file {file_path} loaded successfully.\n")

    for idx, row in df.iterrows():
        name = row.get('name', f"entry_{idx}")
        oxidised_smiles = row.get('oxidised_smiles', "")
        reduced_smiles = row.get('reduced_smiles', "")

        # Only process if both oxidised and reduced SMILES are present
        if oxidised_smiles and reduced_smiles:
            print("="*60)
            print(f"Processing row {idx}: {name}")
            print(f"  Oxidised SMILES: {oxidised_smiles}")
            print(f"  Reduced SMILES:  {reduced_smiles}")

            # Oxidised neutral
            try:
                pmg_mol_ox, charge_ox, mult_ox = process_smiles(oxidised_smiles, "oxidised")
                oxidised_filename = f"{name}_oxidised_neutral.gjf"
                create_gaussian_input(
                    pmg_mol_ox,
                    charge=charge_ox,
                    spin_multiplicity=mult_ox,
                    filename=oxidised_filename
                )
            except ValueError as e:
                print(f"Error processing oxidised SMILES for {name}: {e}")

            # Reduced neutral
            try:
                pmg_mol_red, charge_red, mult_red = process_smiles(reduced_smiles, "reduced")
                reduced_filename = f"{name}_reduced_neutral.gjf"
                create_gaussian_input(
                    pmg_mol_red,
                    charge=charge_red,
                    spin_multiplicity=mult_red,
                    filename=reduced_filename
                )
            except ValueError as e:
                print(f"Error processing reduced SMILES for {name}: {e}")
        else:
            # Skip rows where either column is missing
            print(f"Skipping {name} - must have both oxidised and reduced SMILES.\n")

if __name__ == "__main__":
    main()
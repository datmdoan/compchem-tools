#!/usr/bin/env python3
"""
smiles_to_structure.py

Converts SMILES strings to 3D molecular structures and generates Gaussian 16
input files (.gjf) for DFT calculations.

Two modes of operation
----------------------
pair   Reads a CSV with ``oxidised_smiles`` and ``reduced_smiles`` columns.
       Auto-detects charge and multiplicity from the SMILES string (via
       RDKit formal-charge / radical-electron counts).  Produces one .gjf
       per oxidation state per row.

ion    Reads a CSV with a ``SMILES`` column and a 1-based row
       index (``-i``).  Generates a neutral (charge 0, multiplicity 1) and
       a cation (charge +1, multiplicity 2) .gjf for that single molecule.

Usage
-----
    # Neutral redox-pair workflow (one row per structure_N subfolder)
    python smiles_to_structure.py pair -f structure_1_oxidised.csv

    # Ion workflow (called from a PBS array job)
    python smiles_to_structure.py ion -f molecules.csv -i 3
"""

import argparse
import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

# ── Gaussian parameters ────────────────────────────────────────────────
FUNCTIONAL = "M062X"
BASIS_SET = "6-31++G(d,p)"
MEMORY = "64GB"
NPROC = "16"
NUM_CONFORMERS = 10


# ── Helper functions ───────────────────────────────────────────────────

def get_charge_and_multiplicity(smiles):
    """Return (charge, multiplicity) for a SMILES string.

    Uses RDKit formal charges and radical-electron counts.  This is an
    approximation — complex open-shell systems may need manual input.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)
    charge = Chem.GetFormalCharge(mol)
    unpaired = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    multiplicity = unpaired + 1 if unpaired > 0 else 1
    return charge, multiplicity


def generate_3d_structure(smiles):
    """Convert a SMILES string to a 3D RDKit molecule."""
    print(f"Generating 3D structure for SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        raise ValueError(f"Could not generate 3D conformer for SMILES: {smiles}")
    return mol


def conformer_search(mol, num_conformers=NUM_CONFORMERS):
    """Run an MMFF94 conformer search and return (mol, best_conf_id)."""
    print(f"Conformer search ({num_conformers} conformers, MMFF94).")
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers,
                                          randomSeed=42)
    energies = []
    for cid in conf_ids:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        if props is None:
            print(f"  MMFF94 properties unavailable for conformer {cid}, skipping.")
            continue
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        if ff is None:
            print(f"  MMFF94 optimisation failed for conformer {cid}, skipping.")
            continue
        ff.Minimize()
        energies.append((cid, ff.CalcEnergy()))

    if not energies:
        raise ValueError("No valid conformers found.  Check the input molecule.")

    best_id, best_e = min(energies, key=lambda x: x[1])
    print(f"  Lowest-energy conformer: {best_id}  E = {best_e:.2f} kcal/mol")
    return mol, best_id


def rdkit_to_pymatgen(mol, conf_id):
    """Convert an RDKit conformer to a pymatgen Molecule."""
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = mol.GetConformer(conf_id).GetPositions()
    return Molecule(atoms, coords)


def write_gaussian_input(pmg_mol, charge, multiplicity, filename):
    """Write a Gaussian 16 .gjf file."""
    print(f"Writing Gaussian input: {filename}")
    chk = filename.replace(".gjf", ".chk")
    route = {"opt": "", "freq": "NoRaman"}

    gi = GaussianInput(
        pmg_mol,
        charge=charge,
        spin_multiplicity=multiplicity,
        functional=FUNCTIONAL,
        basis_set=BASIS_SET,
        route_parameters=route,
        link0_parameters={
            "%mem": MEMORY,
            "%nprocshared": NPROC,
            "%Chk": chk,
        },
        title="Generated Gaussian Input",
    )
    if os.path.exists(filename):
        print(f"  Warning: overwriting {filename}")
    gi.write_file(filename)
    print(f"  {filename} created.\n")


# ── Mode: pair ─────────────────────────────────────────────────────────

def run_pair_mode(csv_path):
    """Process a CSV with oxidised_smiles / reduced_smiles columns.

    Charge and multiplicity are auto-detected from the SMILES string.
    """
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows from {csv_path}\n")

    for idx, row in df.iterrows():
        name = row.get("name", f"entry_{idx}")
        ox_smi = row.get("oxidised_smiles", "")
        red_smi = row.get("reduced_smiles", "")

        if not ox_smi or not red_smi:
            print(f"Skipping {name} — need both oxidised and reduced SMILES.\n")
            continue

        print("=" * 60)
        print(f"Row {idx}: {name}")
        print(f"  Oxidised SMILES: {ox_smi}")
        print(f"  Reduced SMILES:  {red_smi}")

        # ── oxidised ──
        try:
            charge, mult = get_charge_and_multiplicity(ox_smi)
            mol = generate_3d_structure(ox_smi)
            mol, cid = conformer_search(mol)
            pmg = rdkit_to_pymatgen(mol, cid)
            write_gaussian_input(pmg, charge, mult,
                                 f"{name}_oxidised_neutral.gjf")
        except ValueError as exc:
            print(f"  Error (oxidised) for {name}: {exc}")

        # ── reduced ──
        try:
            charge, mult = get_charge_and_multiplicity(red_smi)
            mol = generate_3d_structure(red_smi)
            mol, cid = conformer_search(mol)
            pmg = rdkit_to_pymatgen(mol, cid)
            write_gaussian_input(pmg, charge, mult,
                                 f"{name}_reduced_neutral.gjf")
        except ValueError as exc:
            print(f"  Error (reduced) for {name}: {exc}")


# ── Mode: ion ──────────────────────────────────────────────────────────

def run_ion_mode(csv_path, index):
    """Process a single row from a CSV with a ``SMILES`` column.

    Generates two .gjf files:
      * neutral  — charge 0, multiplicity 1
      * cation   — charge +1, multiplicity 2

    Parameters
    ----------
    csv_path : str
        Path to the master CSV.
    index : int
        1-based row index (matches PBS_ARRAY_INDEX).
    """
    df = pd.read_csv(csv_path)
    idx_0 = index - 1
    if idx_0 < 0 or idx_0 >= len(df):
        raise ValueError(f"Index {index} out of range (CSV has {len(df)} rows).")

    smiles = df.iloc[idx_0]["SMILES"]
    label = f"structure_{index}"
    print(f"Processing {label}: {smiles}\n")

    mol = generate_3d_structure(smiles)
    mol, cid = conformer_search(mol)
    pmg = rdkit_to_pymatgen(mol, cid)

    write_gaussian_input(pmg, charge=0, multiplicity=1,
                         filename=f"neutral_gas_{label}.gjf")
    write_gaussian_input(pmg, charge=1, multiplicity=2,
                         filename=f"cation_gas_{label}.gjf")


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate Gaussian 16 input files from SMILES.",
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    # pair sub-command
    p_pair = sub.add_parser(
        "pair",
        help="Neutral redox-pair workflow (oxidised + reduced columns).",
    )
    p_pair.add_argument("-f", "--file", required=True,
                        help="CSV with oxidised_smiles and reduced_smiles.")

    # ion sub-command
    p_ion = sub.add_parser(
        "ion",
        help="Ion workflow (neutral + cation from a single SMILES).",
    )
    p_ion.add_argument("-f", "--file", required=True,
                       help="CSV with a 'SMILES' column.")
    p_ion.add_argument("-i", "--index", type=int, required=True,
                       help="1-based row index (PBS_ARRAY_INDEX).")

    args = parser.parse_args()

    if args.mode == "pair":
        run_pair_mode(args.file)
    elif args.mode == "ion":
        run_ion_mode(args.file, args.index)


if __name__ == "__main__":
    main()

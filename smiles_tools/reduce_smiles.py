"""
reduce_smiles.py

Generates reduced form SMILES from oxidised form SMILES for phenazines,
alloxazines and quinones. Each molecule class has its own reduction function.

CLI:
    python reduce_smiles.py --input molecules.csv --type phenazine --output reduced.csv
    python reduce_smiles.py --input molecules.csv --type quinone   --output reduced.csv

Input CSV must have a column named 'SMILES'.

"""

import os
import argparse
import pandas as pd
from rdkit import Chem

try:
    import networkx as nx
    _NX_AVAILABLE = True
except ImportError:
    _NX_AVAILABLE = False


def reduce_phenazine(smiles):
    """
    Reduce a phenazine by protonating bare aromatic nitrogens ([n&H0]).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")

    pattern = Chem.MolFromSmarts("[n&H0]")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return smiles

    emol = Chem.EditableMol(mol)
    new_mol = emol.GetMol()

    for match in matches:
        n_atom = new_mol.GetAtomWithIdx(match[0])
        n_atom.SetNumExplicitHs(n_atom.GetNumExplicitHs() + 1)

    try:
        Chem.Kekulize(new_mol, clearAromaticFlags=True)
    except Chem.KekulizeException as e:
        raise ValueError(f"Kekulisation failed for {smiles}: {e}")

    return Chem.MolToSmiles(new_mol, canonical=True)


def reduce_alloxazine(smiles):
    """
    Reduce an alloxazine by converting the ring N=C-C=N motif to NH-C=C-NH.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        raise ValueError(f"Kekulisation failed for SMILES: {smiles}")

    pattern = Chem.MolFromSmarts("[#7;R]=[#6]-[#6]=[#7;R]")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return Chem.MolToSmiles(mol)

    emol = Chem.EditableMol(mol)

    for match in matches:
        n1_idx, c1_idx, c2_idx, n2_idx = match

        emol.RemoveBond(n1_idx, c1_idx)
        emol.RemoveBond(n2_idx, c2_idx)
        emol.AddBond(n1_idx, c1_idx, Chem.BondType.SINGLE)
        emol.AddBond(n2_idx, c2_idx, Chem.BondType.SINGLE)

        emol.RemoveBond(c1_idx, c2_idx)
        emol.AddBond(c1_idx, c2_idx, Chem.BondType.DOUBLE)

        n1_atom = emol.GetMol().GetAtomWithIdx(n1_idx)
        n1_atom.SetNumExplicitHs(n1_atom.GetNumExplicitHs() + 1)
        n2_atom = emol.GetMol().GetAtomWithIdx(n2_idx)
        n2_atom.SetNumExplicitHs(n2_atom.GetNumExplicitHs() + 1)

    result_mol = emol.GetMol()
    Chem.SanitizeMol(result_mol)
    return Chem.MolToSmiles(result_mol, canonical=True)


def _mol_to_networkx(mol):
    """
    Build a NetworkX graph from an RDKit Mol (nodes are atom indices).
    """  
    if not _NX_AVAILABLE:
        raise ImportError("networkx is required for quinone reduction. Install with: pip install networkx")
    G = nx.Graph()
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G


def _try_path_flip(mol, path):
    """
    Flip single <-> double bonds along an atom-index path.
    Returns a sanitised Mol on success, or None if sanitisation fails.
    """
    emol = Chem.EditableMol(mol)
    for i in range(len(path) - 1):
        a1, a2 = path[i], path[i + 1]
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            return None
        bt = bond.GetBondType()
        emol.RemoveBond(a1, a2)
        if bt == Chem.BondType.SINGLE:
            emol.AddBond(a1, a2, Chem.BondType.DOUBLE)
        elif bt == Chem.BondType.DOUBLE:
            emol.AddBond(a1, a2, Chem.BondType.SINGLE)
        else:
            return None

    candidate = emol.GetMol()
    try:
        Chem.SanitizeMol(candidate)
    except Exception:
        return None
    return candidate


def reduce_quinone_like(smiles, max_path_length=10):
    """
    Reduce a quinone by converting both ring C=O groups to C-OH, then
    re-establishing conjugation via bond-alternation along paths between the
    two carbonyl carbons (enumerated with NetworkX).

    Requires exactly two ring carbonyls ([C;R]=O). 
    max_path_length controls how far the path search goes
    Returns the original SMILES if reduction is not achievable.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return Chem.MolToSmiles(mol)

    pattern = Chem.MolFromSmarts("[C;R]=O")
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) != 2:
        return Chem.MolToSmiles(mol)

    (c1_idx, o1_idx), (c2_idx, o2_idx) = matches

    emol = Chem.EditableMol(mol)
    emol.RemoveBond(c1_idx, o1_idx)
    emol.RemoveBond(c2_idx, o2_idx)
    emol.AddBond(c1_idx, o1_idx, Chem.BondType.SINGLE)
    emol.AddBond(c2_idx, o2_idx, Chem.BondType.SINGLE)

    temp_mol = emol.GetMol()
    for o_idx in [o1_idx, o2_idx]:
        temp_mol.GetAtomWithIdx(o_idx).SetNumExplicitHs(
            temp_mol.GetAtomWithIdx(o_idx).GetNumExplicitHs() + 1
        )

    G = _mol_to_networkx(temp_mol)  # raises ImportError if networkx not installed
    try:
        paths = list(nx.all_simple_paths(G, source=c1_idx, target=c2_idx, cutoff=max_path_length))
    except Exception as e:
        if "NetworkXNoPath" in type(e).__name__ or "NodeNotFound" in type(e).__name__:
            return Chem.MolToSmiles(temp_mol)
        raise

    for path in paths:
        candidate = _try_path_flip(temp_mol, path)
        if candidate is not None:
            try:
                Chem.Kekulize(candidate, clearAromaticFlags=True)
            except Chem.KekulizeException:
                continue
            return Chem.MolToSmiles(candidate, canonical=True)


    try:
        Chem.SanitizeMol(temp_mol)
    except Exception:
        pass
    return Chem.MolToSmiles(temp_mol, canonical=True)


SUPPORTED_TYPES = ("phenazine", "alloxazine", "quinone")


def reduce_smiles(smiles, molecule_type):
    """
    molecule_type must be one of: 'phenazine', 'alloxazine', 'quinone'.
    """
    molecule_type = molecule_type.lower()
    if molecule_type == "phenazine":
        return reduce_phenazine(smiles)
    elif molecule_type == "alloxazine":
        return reduce_alloxazine(smiles)
    elif molecule_type == "quinone":
        return reduce_quinone_like(smiles)
    else:
        raise ValueError(f"Unknown molecule type '{molecule_type}'. Choose from: {SUPPORTED_TYPES}")



def process_csv(input_csv, output_csv, molecule_type, smiles_col="SMILES"):
    """Read a CSV, reduce all SMILES, write a new CSV with a 'Reduced SMILES' column."""
    df = pd.read_csv(input_csv)
    if smiles_col not in df.columns:
        raise ValueError(f"Column '{smiles_col}' not found in {input_csv}. "
                         f"Available columns: {list(df.columns)}")

    def _safe_reduce(smi):
        try:
            return reduce_smiles(smi, molecule_type)
        except Exception as e:
            return f"Error: {e}"

    df["Reduced SMILES"] = df[smiles_col].apply(_safe_reduce)

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} reduced SMILES to: {output_csv}")

    n_errors = df["Reduced SMILES"].str.startswith("Error:").sum()
    if n_errors:
        print(f"  Warning: {n_errors} entries failed reduction.")


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Produce reduced (hydrogenated) SMILES from oxidised molecules.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reduce_smiles.py --input quinones.csv --type quinone --output quinones_reduced.csv
  python reduce_smiles.py --input phenazines.csv --type phenazine --output phenazines_reduced.csv
  python reduce_smiles.py --input alloxazines.csv --type alloxazine --output alloxazines_reduced.csv
        """
    )
    parser.add_argument("--input",  required=True, help="Path to input CSV file.")
    parser.add_argument("--output", required=True, help="Path for output CSV file.")
    parser.add_argument("--type",   required=True, choices=SUPPORTED_TYPES,
                        help="Molecule class to reduce.")
    parser.add_argument("--smiles-col", default="SMILES",
                        help="Name of the SMILES column in the CSV (default: 'SMILES').")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    process_csv(
        input_csv=args.input,
        output_csv=args.output,
        molecule_type=args.type,
        smiles_col=args.smiles_col,
    )

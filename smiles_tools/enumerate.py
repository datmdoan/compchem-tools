#!/usr/bin/env python3
"""
enumerate.py

This script generates functionalised derivatives SMILES from backbone SMILES by attaching
functional groups to available ring positions.

Three modes:
  single    one FG per molecule per position
  multi     all combinations of FGs across positions
  targeted  specific FGs at specific positions (position_map dict)

By default, only positions on the unmodified core ring are used
(scaffold_filter=True). Ring N-H positions are handled automatically.

Example usage::

    from smiles_tools.enumerate import enumerate_library

    df = enumerate_library(
        backbones=["c1ccc2c(c1)Nc3ccccc3S2"],
        functional_groups=["Cl", "[OH]", "N(C)(C)", "C#N"],
        mode="single",
    )

"""

import argparse
import os
from datetime import datetime
from itertools import combinations, product

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


# ── SMILES normalisation ──────────────────────────────────────────────

def normalize_smiles(smiles, method="canonical"):
    """Normalise a SMILES string.

    method="canonical" uses RDKit canonical SMILES directly.
    method="inchi" round-trips through InChI 
    Returns None if the input cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if method == "inchi":
        inchi = Chem.MolToInchi(mol)
        if inchi is None:
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        mol2 = Chem.MolFromInchi(inchi)
        if mol2 is None:
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return Chem.MolToSmiles(mol2, isomericSmiles=True, canonical=True)
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


# ── Position-detection helpers ─────────────────────────────────────────

def is_scaffold_atom(mol, atom_idx):
    """True if the atom is in a ring with no exocyclic substituents.

    Any neighbour outside a ring is treated as a pre-existing substituent.
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    if not atom.IsInRing():
        return False
    for neighbour in atom.GetNeighbors():
        if not neighbour.IsInRing():
            return False
    return True


def is_h_bearing_ring_nitrogen(mol, atom_idx):
    """True if the atom is a ring nitrogen with a hydrogen attached."""
    atom = mol.GetAtomWithIdx(atom_idx)
    return (
        atom.GetSymbol() == "N"
        and atom.IsInRing()
        and atom.GetTotalNumHs() > 0
    )


def get_valid_positions(mol, scaffold_filter=True):
    """Return (normal_indices, nh_indices) for a molecule.

    normal_indices: positions for standard FG attachment.
    nh_indices: ring N-H positions that need add_functional_group_replace_h.
    Set scaffold_filter=False to include any atom with available valence.
    """
    nh_indices = [
        i for i in range(mol.GetNumAtoms())
        if is_h_bearing_ring_nitrogen(mol, i)
    ]
    nh_set = set(nh_indices)

    if scaffold_filter:
        normal_indices = [
            i for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetImplicitValence() > 0
            and is_scaffold_atom(mol, i)
            and i not in nh_set
        ]
    else:
        normal_indices = [
            i for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetImplicitValence() > 0
            and i not in nh_set
        ]
    return normal_indices, nh_indices


# ── Functional-group attachment ────────────────────────────────────────

def add_functional_group(mol, fg_smiles, atom_index):
    """Attach fg_smiles to mol at atom_index via a single bond."""
    fg = Chem.MolFromSmiles(fg_smiles)
    if fg is None:
        raise ValueError(f"Invalid functional-group SMILES: {fg_smiles}")
    mol_copy = Chem.RWMol(mol)
    combined = Chem.CombineMols(mol_copy, fg)
    rw = Chem.RWMol(combined)
    fg_start = mol_copy.GetNumAtoms()
    rw.AddBond(atom_index, fg_start, Chem.BondType.SINGLE)
    return rw.GetMol()


def add_functional_group_replace_h(mol, fg_smiles, atom_index):
    """Attach fg_smiles to mol at atom_index, replacing an existing H.

    Used for N-H → N-R substitution on ring nitrogens.
    """
    fg = Chem.MolFromSmiles(fg_smiles)
    if fg is None:
        raise ValueError(f"Invalid functional-group SMILES: {fg_smiles}")
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(atom_index)
    h_count = atom.GetTotalNumHs()
    if h_count > 0:
        atom.SetNumExplicitHs(h_count - 1)
    combined = Chem.CombineMols(rw, fg)
    rw2 = Chem.RWMol(combined)
    fg_start = rw.GetNumAtoms()
    rw2.AddBond(atom_index, fg_start, Chem.BondType.SINGLE)
    result = rw2.GetMol()
    Chem.SanitizeMol(result)
    return result


def _attach_fg(mol, fg_smiles, atom_index, nh_set):
    """Route to the correct attachment function based on atom type."""
    if atom_index in nh_set:
        return add_functional_group_replace_h(mol, fg_smiles, atom_index)
    return add_functional_group(mol, fg_smiles, atom_index)


# ── Enumeration strategies ─────────────────────────────────────────────

def enumerate_single(mol, functional_groups, backbone_smiles,
                     scaffold_filter=True, normalisation="canonical"):
    """Single-substitution enumeration: one FG per molecule per position."""
    normal, nh = get_valid_positions(mol, scaffold_filter)
    nh_set = set(nh)
    all_indices = normal + nh

    rows = [{
        "Backbone": backbone_smiles,
        "SMILES": normalize_smiles(
            Chem.MolToSmiles(mol, isomericSmiles=True), normalisation),
        "Position": "None",
        "Functional Group": "None",
    }]

    for fg in functional_groups:
        for idx in all_indices:
            try:
                deriv = _attach_fg(mol, fg, idx, nh_set)
                smi = Chem.MolToSmiles(deriv, isomericSmiles=True, canonical=True)
                rows.append({
                    "Backbone": backbone_smiles,
                    "SMILES": normalize_smiles(smi, normalisation),
                    "Position": idx,
                    "Functional Group": fg,
                })
            except Exception as exc:
                print(f"  Skipping {backbone_smiles} pos {idx} + {fg}: {exc}")
    return rows


def enumerate_multi(mol, functional_groups, backbone_smiles,
                    atom_indices=None, scaffold_filter=True,
                    normalisation="canonical"):
    """Multi-substitution enumeration: all combinations of FGs across positions.

    Positions are auto-detected if atom_indices is None.
    """
    if atom_indices is not None:
        nh = [
            i for i in atom_indices
            if is_h_bearing_ring_nitrogen(mol, i)
        ]
        normal = [i for i in atom_indices if i not in set(nh)]
        all_indices = normal + nh
    else:
        normal, nh = get_valid_positions(mol, scaffold_filter)
        all_indices = normal + nh

    nh_set = set(nh)

    rows = [{
        "Backbone": backbone_smiles,
        "SMILES": normalize_smiles(
            Chem.MolToSmiles(mol, isomericSmiles=True), normalisation),
        "Position": "None",
        "Functional Group": "None",
    }]

    for r in range(1, len(all_indices) + 1):
        for selected in combinations(all_indices, r):
            for fg_combo in product(functional_groups, repeat=r):
                try:
                    temp = mol
                    for fg, idx in zip(fg_combo, selected):
                        temp = _attach_fg(temp, fg, idx, nh_set)
                    smi = Chem.MolToSmiles(temp, isomericSmiles=True,
                                           canonical=True)
                    rows.append({
                        "Backbone": backbone_smiles,
                        "SMILES": normalize_smiles(smi, normalisation),
                        "Position": ", ".join(map(str, selected)),
                        "Functional Group": ", ".join(fg_combo),
                    })
                except Exception as exc:
                    pos_str = ", ".join(map(str, selected))
                    fg_str = ", ".join(fg_combo)
                    print(f"  Skipping {backbone_smiles} pos [{pos_str}] "
                          f"+ [{fg_str}]: {exc}")
    return rows


def enumerate_targeted(mol, position_map, backbone_smiles,
                       normalisation="canonical"):
    """Enumerate with specific FGs assigned to specific positions.

    position_map is a dict of {atom_index: [fg_smiles, ...]}.
    Every combination across positions is generated (one FG choice per
    position per molecule), e.g.::

        position_map = {
            3: ["Cl", "[OH]"],
            7: ["N(C)(C)", "C#N"],
        }

    Positions can be left out entirely (the backbone is always included).
    N-H ring nitrogens are handled automatically.
    """
    nh_set = set(
        i for i in position_map
        if is_h_bearing_ring_nitrogen(mol, i)
    )

    rows = [{
        "Backbone": backbone_smiles,
        "SMILES": normalize_smiles(
            Chem.MolToSmiles(mol, isomericSmiles=True), normalisation),
        "Position": "None",
        "Functional Group": "None",
    }]

    positions = list(position_map.keys())
    fg_options = [position_map[p] for p in positions]

    # Generate every combination: one FG choice per position
    for fg_combo in product(*fg_options):
        try:
            temp = mol
            for idx, fg in zip(positions, fg_combo):
                temp = _attach_fg(temp, fg, idx, nh_set)
            smi = Chem.MolToSmiles(temp, isomericSmiles=True, canonical=True)
            rows.append({
                "Backbone": backbone_smiles,
                "SMILES": normalize_smiles(smi, normalisation),
                "Position": ", ".join(map(str, positions)),
                "Functional Group": ", ".join(fg_combo),
            })
        except Exception as exc:
            fg_str = ", ".join(fg_combo)
            print(f"  Skipping {backbone_smiles} [{fg_str}]: {exc}")

    return rows



def enumerate_library(backbones, functional_groups=None, mode="single",
                      atom_indices=None, scaffold_filter=True,
                      normalisation="canonical", position_map=None):
    """Build a deduplicated DataFrame of all functionalised derivatives.

    backbones          list of backbone SMILES
    functional_groups  list of FG SMILES (not required for 'targeted' mode)
    mode               'single', 'multi', or 'targeted'
    atom_indices       explicit positions for multi mode (auto if None)
    scaffold_filter    restrict to core ring atoms (True by default)
    normalisation      'canonical' or 'inchi'
    position_map       dict of {atom_index: [fg_smiles, ...]} for targeted mode,
                       e.g. {3: ["Cl", "[OH]"], 7: ["N(C)(C)"]}

    Returns a DataFrame with columns: Backbone, SMILES, Position,
    Functional Group.
    """
    all_rows = []
    for bb in backbones:
        mol = Chem.MolFromSmiles(bb)
        if mol is None:
            print(f"Could not parse backbone SMILES '{bb}' — skipping.")
            continue

        if mode == "single":
            rows = enumerate_single(mol, functional_groups, bb,
                                    scaffold_filter, normalisation)
        elif mode == "multi":
            rows = enumerate_multi(mol, functional_groups, bb,
                                   atom_indices, scaffold_filter,
                                   normalisation)
        elif mode == "targeted":
            if not position_map:
                raise ValueError("'targeted' mode requires a position_map dict.")
            rows = enumerate_targeted(mol, position_map, bb, normalisation)
        else:
            raise ValueError(f"Unknown mode '{mode}'. Use 'single', 'multi', or 'targeted'.")
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df.drop_duplicates(subset=["SMILES"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    print(f"{len(df)} unique structures after deduplication.")
    return df


def display_molecules(smiles_list, mols_per_row=4, img_size=(300, 300),
                      batch_size=50):
    """Draw a grid of 2D structures from a SMILES list (Jupyter only)."""
    from rdkit.Chem import Draw
    from IPython.display import display as ipy_display

    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i : i + batch_size]
        mols = [Chem.MolFromSmiles(s) for s in batch if s]
        img = Draw.MolsToGridImage(
            mols, molsPerRow=mols_per_row, subImgSize=img_size,
            legends=batch,
        )
        ipy_display(img)


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Enumerate functionalised derivatives from backbone SMILES.",
    )
    parser.add_argument("--backbones", nargs="+", required=True,
                        help="Backbone SMILES.")
    parser.add_argument("--fg", nargs="+", default=None,
                        help="Functional-group SMILES (not needed for targeted mode).")
    parser.add_argument("--mode", choices=["single", "multi", "targeted"], default="single",
                        help="Enumeration mode (default: single).")
    parser.add_argument("--indices", nargs="*", type=int, default=None,
                        help="Atom indices for multi mode.")
    parser.add_argument("--position-map", default=None,
                        help='JSON dict for targeted mode, e.g. \'{"3": ["Cl", "[OH]"], "7": ["N(C)(C)"]}\'')
    parser.add_argument("--no-scaffold-filter", action="store_true",
                        help="Allow functionalisation at any available atom.")
    parser.add_argument("--normalisation", choices=["canonical", "inchi"],
                        default="canonical", help="SMILES normalisation method.")
    parser.add_argument("--output", "-o", default=None,
                        help="Output CSV (default: timestamped filename).")
    args = parser.parse_args()

    position_map = None
    if args.position_map:
        import json
        raw = json.loads(args.position_map)
        position_map = {int(k): v for k, v in raw.items()}

    if args.mode == "targeted" and position_map is None:
        parser.error("--position-map is required for targeted mode.")
    if args.mode != "targeted" and args.fg is None:
        parser.error("--fg is required for single and multi modes.")

    df = enumerate_library(
        backbones=args.backbones,
        functional_groups=args.fg,
        mode=args.mode,
        atom_indices=args.indices,
        scaffold_filter=not args.no_scaffold_filter,
        normalisation=args.normalisation,
        position_map=position_map,
    )

    if args.output:
        out_path = args.output
    else:
        ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out_path = f"enumerated_derivatives_{ts}.csv"

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()

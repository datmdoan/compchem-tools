"""
smiles_tools

RDKit cheminformatics utilities:

- **reduce_smiles** — generate reduced (hydrogenated) SMILES from oxidised
  phenazines, alloxazines and quinones.
- **enumerate** — combinatorial enumeration of functionalised derivatives
  from backbone scaffolds.
"""

from .reduce_smiles import (
    reduce_phenazine,
    reduce_alloxazine,
    reduce_quinone_like,
    reduce_smiles,
    process_csv,
)

from .enumerate import (
    enumerate_library,
    enumerate_single,
    enumerate_multi,
    enumerate_targeted,
    normalize_smiles,
    add_functional_group,
    add_functional_group_replace_h,
    get_substitutable_positions,
    is_scaffold_atom,
    is_h_bearing_ring_nitrogen,
    display_molecules,
)

__all__ = [
    # reduce_smiles
    "reduce_phenazine",
    "reduce_alloxazine",
    "reduce_quinone_like",
    "reduce_smiles",
    "process_csv",
    # enumerate
    "enumerate_library",
    "enumerate_single",
    "enumerate_multi",
    "enumerate_targeted",
    "normalize_smiles",
    "add_functional_group",
    "add_functional_group_replace_h",
    "get_substitutable_positions",
    "is_scaffold_atom",
    "is_h_bearing_ring_nitrogen",
    "display_molecules",
]

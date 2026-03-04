#!/usr/bin/env python3
"""
autodft_pairs.py

Launcher for the two-electron PCET reaction DFT workflow on HPC.

For each valid row in the master CSV this script:
  1. Creates ``structure_N/oxidised/`` and ``structure_N/reduced/`` folders
     under RESULTS_DIRECTORY.
  2. Writes a mini one-row CSV and copies ``smiles_to_structure.py`` into
     each subfolder.
  3. Fills placeholders in ``autodft_pairs.pbs`` and submits via ``qsub``.

The PBS array uses odd indices for oxidised and even for reduced so that a
single array covers both oxidation states (2 × N_structures subjobs).

Usage
-----
    python autodft_pairs.py
"""

import os
import shutil
import subprocess

import pandas as pd

# ── User settings ──────────────────────────────────────────────────────
RESULTS_DIRECTORY = "/path/to/results"
MASTER_CSV = "/path/to/molecules.csv"
SMILES_TO_STRUCTURE_SCRIPT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "smiles_to_structure.py"
)
# ───────────────────────────────────────────────────────────────────────


def main():
    if not os.path.isfile(SMILES_TO_STRUCTURE_SCRIPT):
        raise FileNotFoundError(
            f"smiles_to_structure.py not found: {SMILES_TO_STRUCTURE_SCRIPT}"
        )
    if not os.path.isfile(MASTER_CSV):
        raise FileNotFoundError(f"CSV file not found: {MASTER_CSV}")

    df = pd.read_csv(MASTER_CSV)

    os.makedirs(RESULTS_DIRECTORY, exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIRECTORY, "o_files"), exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIRECTORY, "e_files"), exist_ok=True)

    valid_rows = 0
    for _, row in df.iterrows():
        ox = row.get("oxidised_smiles", "")
        red = row.get("reduced_smiles", "")

        if not isinstance(ox, str) or not isinstance(red, str) or not ox or not red:
            continue

        valid_rows += 1
        struct_folder = os.path.join(RESULTS_DIRECTORY, f"structure_{valid_rows}")

        for state, smi_ox, smi_red in [
            ("oxidised", ox, ""),
            ("reduced", "", red),
        ]:
            sub = os.path.join(struct_folder, state)
            os.makedirs(sub, exist_ok=True)

            # Mini CSV for smiles_to_structure.py pair mode
            mini_csv = os.path.join(
                sub, f"structure_{valid_rows}_{state}.csv"
            )
            pd.DataFrame(
                {"oxidised_smiles": [smi_ox], "reduced_smiles": [smi_red]}
            ).to_csv(mini_csv, index=False)

            shutil.copy(SMILES_TO_STRUCTURE_SCRIPT, sub)

    total_subjobs = 2 * valid_rows
    print(f"{valid_rows} valid structures → {total_subjobs} PBS subjobs.")

    # Fill the PBS template
    pbs_template = os.path.join(os.path.dirname(__file__), "autodft_pairs.pbs")
    if not os.path.isfile(pbs_template):
        raise FileNotFoundError(f"PBS template not found: {pbs_template}")

    with open(pbs_template) as fh:
        content = fh.read()

    content = (
        content
        .replace("__TOTAL_SUBJOBS__", str(total_subjobs))
        .replace("__RESULTS_DIRECTORY__", RESULTS_DIRECTORY)
        .replace("__WORKING_DIRECTORY__", os.getcwd())
    )

    submitted = os.path.join(os.getcwd(),
                             "gaussian_pairs_submitted.pbs")
    with open(submitted, "w") as fh:
        fh.write(content)

    print(f"Submitting job: {submitted}")
    subprocess.run(["qsub", submitted])


if __name__ == "__main__":
    main()

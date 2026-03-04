#!/usr/bin/env python3
"""
autodft_ions.py

Launcher for the single-electron electron transfer reaction DFT workflow on HPC.

For each molecule in the input CSV this script:
  1. Creates ``structure_N/`` under RESULTS_DIRECTORY.
  2. Copies ``smiles_to_structure.py`` into each folder.
  3. Fills placeholders in ``autodft_ions.pbs`` and submits via ``qsub``.

The PBS array job then runs Gaussian 16 OPT+FREQ (M062X, B3LYP fallback)
followed by a single-point energy (SPE) using SMD in water.

Usage
-----
    python autodft_ions.py
"""

import os
import shutil
import subprocess

import pandas as pd

# ── User settings ──────────────────────────────────────────────────────
RESULTS_DIRECTORY = "/path/to/results"
CSV_FILE_PATH = "/path/to/molecules.csv"
# ───────────────────────────────────────────────────────────────────────


def main():
    if not os.path.isfile(CSV_FILE_PATH):
        raise FileNotFoundError(f"CSV file not found: {CSV_FILE_PATH}")

    df = pd.read_csv(CSV_FILE_PATH)
    total_molecules = len(df)
    print(f"Found {total_molecules} molecules in {CSV_FILE_PATH}.")

    working_directory = os.getcwd()

    # Ensure smiles_to_structure.py is alongside this script
    s2s_path = os.path.join(working_directory, "smiles_to_structure.py")
    if not os.path.isfile(s2s_path):
        raise FileNotFoundError(
            f"smiles_to_structure.py not found in {working_directory}"
        )

    # Create structure_N folders and copy helper script
    for idx in range(total_molecules):
        folder = os.path.join(RESULTS_DIRECTORY, f"structure_{idx + 1}")
        os.makedirs(folder, exist_ok=True)
        shutil.copy(s2s_path, folder)

    # Load and fill the PBS template
    pbs_template = os.path.join(working_directory, "autodft_ions.pbs")
    if not os.path.isfile(pbs_template):
        raise FileNotFoundError(f"PBS template not found: {pbs_template}")

    with open(pbs_template) as fh:
        content = fh.read()

    content = (
        content
        .replace("__TOTAL_MOLECULES__", str(total_molecules))
        .replace("__RESULTS_DIRECTORY__", RESULTS_DIRECTORY)
        .replace("__WORKING_DIRECTORY__", working_directory)
        .replace("__CSV_FILE_PATH__", CSV_FILE_PATH)
    )

    submitted = os.path.join(working_directory,
                             "gaussian_array_job_submitted.pbs")
    with open(submitted, "w") as fh:
        fh.write(content)

    print(f"Submitting job: {submitted}")
    subprocess.run(["qsub", submitted])


if __name__ == "__main__":
    main()

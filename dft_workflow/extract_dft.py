#!/usr/bin/env python3
"""
extract_dft.py

Parses Gaussian 16 log files produced by the autodft pipeline and extracts a
consolidated CSV of energies, orbital data, spin densities and timings.

Two directory layouts are supported (selected via subcommand):

ion    ``structure_N/*.log``
       Used by the ion (1 e⁻) workflow.  Calculation type is inferred
       from the filename (``spe`` → SPE, otherwise OPT/FREQ).  All .log
       files in each ``structure_N/`` folder are written to the CSV.

pair   ``structure_N/{oxidised,reduced}/*.log``
       Used by the neutral redox-pair (2 e⁻) workflow.  For each
       oxidation state the script picks the *best* log file (proper
       termination, highest resubmission number) per calculation type.

Usage
-----
    # Ion workflow
    python extract_dft.py ion \\
        --root /path/to/outputs \\
        --smiles /path/to/molecules.csv \\
        --output /path/to/extracted.csv

    # Neutral-pair workflow
    python extract_dft.py pair \\
        --root /path/to/outputs \\
        --smiles /path/to/molecules.csv \\
        --output /path/to/extracted.csv
"""

import argparse
import csv
import os
import re

from pymatgen.io.gaussian import GaussianOutput

HARTREE_TO_EV = 27.2114


def read_input_csv(path):
    """Return a list of dicts from *path*, stripping whitespace from keys."""
    with open(path, "r", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        reader.fieldnames = [f.strip() for f in reader.fieldnames]
        return [{k.strip(): v for k, v in row.items()} for row in reader]


def parse_time_to_minutes(text):
    """Convert a Gaussian timing string to minutes."""
    days = hours = minutes = seconds = 0
    m = re.search(r"(\d+)\s+days", text)
    if m:
        days = int(m.group(1))
    m = re.search(r"(\d+)\s+hours", text)
    if m:
        hours = int(m.group(1))
    m = re.search(r"(\d+)\s+minutes", text)
    if m:
        minutes = int(m.group(1))
    m = re.search(r"(\d+\.?\d*)\s+seconds", text)
    if m:
        seconds = float(m.group(1))
    return days * 1440 + hours * 60 + minutes + seconds / 60


def extract_job_times(log_file):
    cpu_time = elapsed_time = None
    with open(log_file) as fh:
        for line in fh:
            if "Job cpu time" in line:
                cpu_time = parse_time_to_minutes(line.split(":")[-1].strip())
            elif "Elapsed time" in line:
                elapsed_time = parse_time_to_minutes(line.split(":")[-1].strip())
    return cpu_time, elapsed_time


def extract_mulliken_spin_densities(log_file):
    densities = {}
    with open(log_file) as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if "Mulliken atomic spin densities" in line or \
           "Mulliken charges and spin densities" in line:
            idx = i + 2
            while idx < len(lines):
                cur = lines[idx].strip()
                if not cur or not re.match(r"^\d+\s", cur):
                    break
                parts = cur.split()
                densities[int(parts[0])] = {
                    "atom_type": parts[1],
                    "spin_density": float(parts[-1]),
                }
                idx += 1
            break
    return densities


def extract_homo_lumo(log_file):
    """Return (alpha_homo, alpha_lumo, beta_homo, beta_lumo) in Hartree."""
    try:
        with open(log_file) as fh:
            lines = fh.readlines()
        alpha_homo = alpha_lumo = beta_homo = beta_lumo = None
        last_alpha_occ = []
        last_beta_occ = []
        collecting_alpha = collecting_beta = False

        for line in lines:
            s = line.strip()
            if "Alpha  occ. eigenvalues" in s:
                collecting_alpha = True
                last_alpha_occ = list(map(float, re.findall(r"(-?\d+\.\d+)", s)))
            elif "Beta  occ. eigenvalues" in s:
                collecting_beta = True
                last_beta_occ = list(map(float, re.findall(r"(-?\d+\.\d+)", s)))
            elif collecting_alpha and "Alpha virt. eigenvalues" in s:
                vals = list(map(float, re.findall(r"(-?\d+\.\d+)", s)))
                alpha_lumo = vals[0] if vals else None
                collecting_alpha = False
            elif collecting_beta and "Beta virt. eigenvalues" in s:
                vals = list(map(float, re.findall(r"(-?\d+\.\d+)", s)))
                beta_lumo = vals[0] if vals else None
                collecting_beta = False

        alpha_homo = last_alpha_occ[-1] if last_alpha_occ else None
        beta_homo = last_beta_occ[-1] if last_beta_occ else None
        return alpha_homo, alpha_lumo, beta_homo, beta_lumo
    except Exception as exc:
        print(f"Error extracting HOMO/LUMO from {log_file}: {exc}")
        return None, None, None, None


def extract_free_energies(log_file):
    """Return the Sum of electronic and thermal Free Energies (Hartree)."""
    try:
        with open(log_file) as fh:
            for line in fh:
                if "Sum of electronic and thermal Free Energies=" in line:
                    m = re.search(r"(-?\d+\.\d+)", line)
                    if m:
                        return float(m.group(1))
    except Exception as exc:
        print(f"Error extracting free energy from {log_file}: {exc}")
    return None


def extract_error_lines(log_file):
    """Return lines containing 'error termination'."""
    errors = []
    with open(log_file) as fh:
        for line in fh:
            if "error termination" in line.lower():
                errors.append(line.strip())
    return "\n".join(errors)


def calc_type_from_filename(filename):
    """Infer calculation type from the Gaussian log filename."""
    low = filename.lower()
    if "spe" in low:
        return "SPE"
    if "freq" in low and "opt" not in low:
        return "FREQ"
    if "opt" in low and "freq" not in low:
        return "OPT"
    return "OPT/FREQ"


def _build_row(log_file, smiles_data, calc_type, structure_id,
               oxidation_state=None, smiles_columns=None):
    """Parse a Gaussian log and return a dict ready for csv.DictWriter."""
    try:
        gout = GaussianOutput(log_file)
    except Exception as exc:
        print(f"Error parsing {log_file}: {exc}")
        return None

    row = {
        "filename": os.path.basename(log_file),
        "structure_id": structure_id,
        "calculation_type": calc_type,
        "charge": gout.charge,
        "spin_multiplicity": gout.spin_multiplicity,
        "functional": gout.functional,
        "basis_set": gout.basis_set,
        "final_energy": gout.final_energy,
        "gibbs_free_energy": gout.corrections.get("Gibbs Free Energy"),
        "sum_electronic_thermal_free_energy": extract_free_energies(log_file),
        "termination_status": gout.properly_terminated,
        "alpha_HOMO": None,
        "alpha_LUMO": None,
        "beta_HOMO": None,
        "beta_LUMO": None,
        "alpha_gap": None,
        "beta_gap": None,
        "alpha_gap_eV": None,
        "beta_gap_eV": None,
        "max_spin_density": None,
        "max_spin_atom_type": None,
        "cpu_time_minutes": None,
        "elapsed_time_minutes": None,
    }

    if oxidation_state is not None:
        row["oxidation_state"] = oxidation_state

    # Propagate input-CSV metadata
    if smiles_columns:
        for col in smiles_columns:
            row[col] = smiles_data.get(col, "")

    ah, al, bh, bl = extract_homo_lumo(log_file)
    row["alpha_HOMO"], row["alpha_LUMO"] = ah, al
    row["beta_HOMO"], row["beta_LUMO"] = bh, bl
    if ah is not None and al is not None:
        row["alpha_gap"] = al - ah
        row["alpha_gap_eV"] = row["alpha_gap"] * HARTREE_TO_EV
    if bh is not None and bl is not None:
        row["beta_gap"] = bl - bh
        row["beta_gap_eV"] = row["beta_gap"] * HARTREE_TO_EV

    spins = extract_mulliken_spin_densities(log_file)
    if spins:
        best = max(spins.items(), key=lambda x: abs(x[1]["spin_density"]))
        row["max_spin_density"] = best[1]["spin_density"]
        row["max_spin_atom_type"] = best[1]["atom_type"]

    cpu, wall = extract_job_times(log_file)
    row["cpu_time_minutes"] = cpu
    row["elapsed_time_minutes"] = wall

    return row



ION_SMILES_COLS = [
    "Backbone", "SMILES", "Position", "Functional Group",
    "Redox Potential", "LogS (AqSolPred v1.0s)", "ETS or HTS", "SA Score",
]


def extract_ion(root, smiles_csv, output_csv):
    """Walk ``structure_N/*.log`` directories (ion workflow)."""
    smiles_list = read_input_csv(smiles_csv)
    id_to_data = {}
    for idx, row in enumerate(smiles_list):
        id_to_data[f"structure_{idx + 1}"] = row

    struct_dirs = sorted(
        [d for d in os.listdir(root)
         if os.path.isdir(os.path.join(root, d)) and d.startswith("structure_")],
        key=lambda s: int(re.match(r"structure_(\d+)", s).group(1)),
    )

    fieldnames = (
        ["filename", "structure_id", "calculation_type"]
        + ION_SMILES_COLS
        + [
            "charge", "spin_multiplicity", "functional", "basis_set",
            "final_energy", "gibbs_free_energy",
            "sum_electronic_thermal_free_energy", "termination_status",
            "alpha_HOMO", "alpha_LUMO", "beta_HOMO", "beta_LUMO",
            "alpha_gap", "beta_gap", "alpha_gap_eV", "beta_gap_eV",
            "max_spin_density", "max_spin_atom_type",
            "cpu_time_minutes", "elapsed_time_minutes",
        ]
    )

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()

        for sd in struct_dirs:
            dirpath = os.path.join(root, sd)
            sdata = id_to_data.get(sd, {})
            logs = [f for f in os.listdir(dirpath) if f.endswith(".log")]

            for fname in sorted(logs):
                fpath = os.path.join(dirpath, fname)
                ctype = calc_type_from_filename(fname)
                row = _build_row(fpath, sdata, ctype, sd,
                                 smiles_columns=ION_SMILES_COLS)
                if row:
                    writer.writerow(row)
                else:
                    print(f"Failed to extract: {fpath}")

    print(f"Results written to {output_csv}")



PAIR_SMILES_COLS = ["ID", "SMILES"]

# Global error collector
_gaussian_errors = []


def _get_grouped_logs(subdir):
    """Return logs grouped by calc type, sorted by resubmission number."""
    logs = [f for f in os.listdir(subdir) if f.endswith(".log")]

    def priority(fn):
        m = re.search(r"resub(\d*)", fn.lower())
        n = int(m.group(1)) if m and m.group(1) else (1 if m else 0)
        return n

    logs.sort(key=priority, reverse=True)
    grouped = {"OPT/FREQ": [], "SPE": []}
    for fn in logs:
        path = os.path.join(subdir, fn)
        if "spe" in fn.lower():
            grouped["SPE"].append(path)
        else:
            grouped["OPT/FREQ"].append(path)
    return grouped


def _select_best(files):
    """Pick the first properly-terminated file (list is pre-sorted)."""
    for fp in files:
        try:
            if GaussianOutput(fp).properly_terminated:
                return fp
        except Exception:
            pass
    return files[0] if files else None


def extract_pair(root, smiles_csv, output_csv):
    """Walk ``structure_N/{oxidised,reduced}/*.log`` directories."""
    smiles_list = read_input_csv(smiles_csv)
    id_to_data = {r.get("ID", ""): r for r in smiles_list if r.get("ID")}

    struct_dirs = sorted(
        [d for d in os.listdir(root)
         if os.path.isdir(os.path.join(root, d)) and d.startswith("structure_")],
        key=lambda s: int(re.match(r"structure_(\d+)", s).group(1)),
    )

    fieldnames = (
        ["filename", "structure_id", "oxidation_state", "calculation_type"]
        + PAIR_SMILES_COLS
        + [
            "charge", "spin_multiplicity", "functional", "basis_set",
            "final_energy", "gibbs_free_energy",
            "sum_electronic_thermal_free_energy", "termination_status",
            "alpha_HOMO", "alpha_LUMO", "beta_HOMO", "beta_LUMO",
            "alpha_gap", "beta_gap", "alpha_gap_eV", "beta_gap_eV",
            "max_spin_density", "max_spin_atom_type",
            "cpu_time_minutes", "elapsed_time_minutes",
        ]
    )

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()

        for sd in struct_dirs:
            dirpath = os.path.join(root, sd)
            sdata = id_to_data.get(sd, {})

            for state in ("oxidised", "reduced"):
                subdir = os.path.join(dirpath, state)
                if not os.path.isdir(subdir):
                    print(f"  {sd}/{state}: directory not found.")
                    continue

                grouped = _get_grouped_logs(subdir)
                if not any(grouped.values()):
                    print(f"  {sd}/{state}: no log files.")
                    continue

                for ctype, files in grouped.items():
                    if not files:
                        continue
                    best = _select_best(files)
                    if best is None:
                        continue

                    row = _build_row(
                        best, sdata, ctype, sd,
                        oxidation_state=state,
                        smiles_columns=PAIR_SMILES_COLS,
                    )
                    if row:
                        writer.writerow(row)
                        if not row.get("termination_status", True):
                            errs = extract_error_lines(best)
                            if errs:
                                hdr = (f"Structure: {sd}, "
                                       f"Oxidation: {state}, "
                                       f"File: {os.path.basename(best)}\n")
                                _gaussian_errors.append(hdr + errs + "\n")

    # Write error summary
    if _gaussian_errors:
        err_path = os.path.join(os.path.dirname(output_csv),
                                "gaussian_errors.txt")
        with open(err_path, "w") as fh:
            fh.write("\n".join(_gaussian_errors))
        print(f"Gaussian errors written to {err_path}")

    print(f"Results written to {output_csv}")


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Extract DFT results from Gaussian log files.",
    )
    sub = parser.add_subparsers(dest="layout", required=True)

    for name, help_text in [
        ("ion", "Ion workflow — structure_N/*.log"),
        ("pair", "Neutral-pair workflow — structure_N/{oxidised,reduced}/*.log"),
    ]:
        sp = sub.add_parser(name, help=help_text)
        sp.add_argument("--root", "-r", required=True,
                        help="Root directory containing structure_N folders.")
        sp.add_argument("--smiles", "-s", required=True,
                        help="Input CSV with SMILES / metadata.")
        sp.add_argument("--output", "-o", required=True,
                        help="Path for the output CSV.")

    args = parser.parse_args()

    if args.layout == "ion":
        extract_ion(args.root, args.smiles, args.output)
    else:
        extract_pair(args.root, args.smiles, args.output)


if __name__ == "__main__":
    main()

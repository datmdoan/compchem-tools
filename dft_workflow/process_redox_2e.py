#!/usr/bin/env python3
"""
process_redox_2e.py

Calculate two-electron PCET redox potentials (E⁰ vs SHE) from DFT data
produced by the neutral redox-pair workflow.

Usage
-----
    python process_redox_2e.py -i extracted.csv -o redox_results.csv
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

# ── Constants ──────────────────────────────────────────────────────────
HARTREE_TO_EV = 27.2114
HARTREE_TO_KJMOL = 2625.5
N_ELECTRONS = 2
FARADAY_KJ = 96.485          # kJ mol⁻¹ V⁻¹
SHE_OFFSET = 4.28            # absolute potential of SHE (V)
G_SOLV_PROTON_KJMOL = -1112.5
ENERGY_DIFF_THRESHOLD = 1e-6  # Hartree — sanity check


def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower()
            for c in re.split(r"(\d+)", s)]


def calculate_redox_potentials(csv_file, output_file=None):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    # Numeric conversions
    for col in ("final_energy", "sum_electronic_thermal_free_energy",
                "gibbs_free_energy"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "calculation_type" not in df.columns:
        print("ERROR: 'calculation_type' column not found.")
        return
    df["calculation_type"] = df["calculation_type"].astype(str).str.strip().str.upper()

    structure_ids = sorted(df["structure_id"].unique(), key=natural_sort_key)
    print(f"Found {len(structure_ids)} unique structures.")

    results = []

    for sid in structure_ids:
        sd = df[df["structure_id"] == sid].copy()

        # Locate the four required rows
        red_gas = _pick(sd, "reduced", "OPT/FREQ")
        ox_gas = _pick(sd, "oxidised", "OPT/FREQ")
        red_solv = _pick(sd, "reduced", "SPE")
        ox_solv = _pick(sd, "oxidised", "SPE")

        if any(x is None for x in (red_gas, ox_gas, red_solv, ox_solv)):
            print(f"{sid}: incomplete data — skipping.")
            results.append(_empty_result(sid, "Incomplete data"))
            continue

        try:
            E_gas_red = red_gas["final_energy"]
            G_gas_red = red_gas["sum_electronic_thermal_free_energy"]
            E_gas_ox = ox_gas["final_energy"]
            G_gas_ox = ox_gas["sum_electronic_thermal_free_energy"]
            E_solv_red = red_solv["final_energy"]
            E_solv_ox = ox_solv["final_energy"]

            for name, val in [("E_gas_red", E_gas_red), ("G_gas_red", G_gas_red),
                              ("E_gas_ox", E_gas_ox), ("G_gas_ox", G_gas_ox),
                              ("E_solv_red", E_solv_red), ("E_solv_ox", E_solv_ox)]:
                if pd.isna(val) or not np.isfinite(val):
                    raise ValueError(f"Invalid energy: {name} = {val}")

            # Sanity check — identical structures?
            if abs(G_gas_red - G_gas_ox) < ENERGY_DIFF_THRESHOLD:
                msg = (f"Gas-phase free-energy difference too small "
                       f"({abs(G_gas_red - G_gas_ox):.2e} Ha)")
                print(f"{sid}: {msg}")
                results.append(_empty_result(sid, msg))
                continue

            # Solution-phase free energies
            G_solv_red = E_solv_red + (G_gas_red - E_gas_red)
            G_solv_ox = E_solv_ox + (G_gas_ox - E_gas_ox)

            # Solvation free energies
            dG_solv_red = G_solv_red - G_gas_red
            dG_solv_ox = G_solv_ox - G_gas_ox

            # Gas-phase reaction free energy
            G_rxn_gas = G_gas_red - G_gas_ox

            # kJ/mol conversions
            dG_solv_red_kj = dG_solv_red * HARTREE_TO_KJMOL
            dG_solv_ox_kj = dG_solv_ox * HARTREE_TO_KJMOL
            G_rxn_gas_kj = G_rxn_gas * HARTREE_TO_KJMOL

            # Full solution-phase redox free energy
            dG_rxn_solv_kj = (
                dG_solv_red_kj
                + G_rxn_gas_kj
                - dG_solv_ox_kj
                - N_ELECTRONS * G_SOLV_PROTON_KJMOL
            )

            E0 = -(dG_rxn_solv_kj / (N_ELECTRONS * FARADAY_KJ)) - SHE_OFFSET

            print(f"{sid}: E⁰ = {E0:.4f} V vs SHE")

            results.append({
                "structure_id": sid,
                "E0_vs_SHE_V": E0,
                "deltaG_redox_solv_kjmol": dG_rxn_solv_kj,
                "deltaG_solv_reduced_kJmol": dG_solv_red_kj,
                "deltaG_solv_oxidised_kJmol": dG_solv_ox_kj,
                "G_rxn_gas_kJmol": G_rxn_gas_kj,
                "reduced_gas_file": red_gas["filename"],
                "oxidised_gas_file": ox_gas["filename"],
                "reduced_solv_file": red_solv["filename"],
                "oxidised_solv_file": ox_solv["filename"],
                "oxidised_alpha_HOMO": ox_gas.get("alpha_HOMO"),
                "oxidised_alpha_LUMO": ox_gas.get("alpha_LUMO"),
                "oxidised_alpha_gap_eV": ox_gas.get("alpha_gap_eV"),
                "reduced_alpha_HOMO": red_gas.get("alpha_HOMO"),
                "reduced_alpha_LUMO": red_gas.get("alpha_LUMO"),
                "reduced_alpha_gap_eV": red_gas.get("alpha_gap_eV"),
                "functional": red_gas.get("functional"),
                "basis_set": red_gas.get("basis_set"),
                "calculation_status": "Success",
            })

        except Exception as exc:
            print(f"{sid}: error — {exc}")
            results.append(_empty_result(sid, str(exc)))

    # Save
    if results:
        out = pd.DataFrame(results)
        out["_sort"] = out["structure_id"].apply(natural_sort_key)
        out = out.sort_values("_sort").drop("_sort", axis=1)

        out_path = output_file or "redox_potential_results.csv"
        out.to_csv(out_path, index=False)

        ok = out["calculation_status"] == "Success"
        print(f"\n{ok.sum()} / {len(out)} structures calculated successfully.")
        if ok.sum() > 0:
            print(f"  Mean E⁰ = {out.loc[ok, 'E0_vs_SHE_V'].mean():.4f} V")
            print(f"  Range   = [{out.loc[ok, 'E0_vs_SHE_V'].min():.4f}, "
                  f"{out.loc[ok, 'E0_vs_SHE_V'].max():.4f}] V")

        # Optional summary CSV
        if ok.sum() > 0:
            summary = out.loc[ok, [
                "structure_id", "E0_vs_SHE_V",
                "oxidised_alpha_gap_eV", "reduced_alpha_gap_eV",
                "oxidised_alpha_HOMO", "reduced_alpha_HOMO",
            ]]
            summary_path = os.path.join(
                os.path.dirname(out_path), "redox_summary.csv"
            )
            summary.to_csv(summary_path, index=False)

        print(f"Results saved to {out_path}")
    else:
        print("No results to save.")


# ── Helpers ────────────────────────────────────────────────────────────

def _pick(df, ox_state, calc_type):
    """Return a dict for the first matching row, or None."""
    mask = (
        (df["calculation_type"] == calc_type)
        & (df["oxidation_state"].str.lower() == ox_state)
    )
    sub = df[mask]
    if sub.empty:
        return None
    row = sub.iloc[0]
    if pd.isna(row["final_energy"]):
        return None
    if calc_type == "OPT/FREQ" and pd.isna(row["sum_electronic_thermal_free_energy"]):
        return None
    return row.to_dict()


def _empty_result(sid, status):
    return {
        "structure_id": sid,
        "E0_vs_SHE_V": None,
        "deltaG_redox_solv_kjmol": None,
        "deltaG_solv_reduced_kJmol": None,
        "deltaG_solv_oxidised_kJmol": None,
        "G_rxn_gas_kJmol": None,
        "calculation_status": status,
    }


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Calculate 2 e⁻ redox potentials from DFT data.",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="CSV from extract_dft.py pair.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output CSV for redox potentials.")
    args = parser.parse_args()
    calculate_redox_potentials(args.input, args.output)


if __name__ == "__main__":
    main()

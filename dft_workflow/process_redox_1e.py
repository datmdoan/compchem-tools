#!/usr/bin/env python3
"""
process_redox_1e.py

Calculates one-electron ET reaction redox potentials (E⁰ vs SHE) from DFT data produced by the ion workflow.

Usage
-----
    python process_redox_1e.py -i extracted.csv -o redox_results.csv
"""

import argparse

import pandas as pd

# ── Constants ──────────────────────────────────────────────────────────
HARTREE_TO_EV = 27.2114
HARTREE_TO_KJMOL = 2625.5
N_ELECTRONS = 1
SHE_OFFSET = 4.28  # absolute potential of SHE (V)


def calculate_redox_potentials(csv_file, output_file):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    # Clean calculation_type
    df["calculation_type"] = (
        df["calculation_type"]
        .astype(str).str.strip().str.upper().str.replace(r"\s+", "", regex=True)
    )

    # Numeric conversions
    df["charge"] = pd.to_numeric(df["charge"], errors="coerce")
    df = df.dropna(subset=["charge"])
    df["charge"] = df["charge"].astype(int)
    df["final_energy"] = pd.to_numeric(df["final_energy"], errors="coerce")
    df["sum_electronic_thermal_free_energy"] = pd.to_numeric(
        df["sum_electronic_thermal_free_energy"], errors="coerce"
    )

    # Drop rows missing required energies
    opt_mask = df["calculation_type"] == "OPT/FREQ"
    df = df[~(opt_mask & df["sum_electronic_thermal_free_energy"].isna())]
    df = df.dropna(subset=["final_energy"])

    results = []
    for structure_id in df["structure_id"].unique():
        sd = df[df["structure_id"] == structure_id]

        # Neutral gas OPT/FREQ
        ng = sd[(sd["calculation_type"] == "OPT/FREQ") & (sd["charge"] == 0)]
        if ng.empty:
            print(f"{structure_id}: missing neutral gas OPT/FREQ — skipping.")
            continue
        ng = ng.iloc[0]

        # Cation gas OPT/FREQ
        cg = sd[(sd["calculation_type"] == "OPT/FREQ") & (sd["charge"] == 1)]
        if cg.empty:
            print(f"{structure_id}: missing cation gas OPT/FREQ — skipping.")
            continue
        cg = cg.iloc[0]

        # Neutral solution SPE
        ns = sd[(sd["calculation_type"] == "SPE") & (sd["charge"] == 0)]
        if ns.empty:
            print(f"{structure_id}: missing neutral solution SPE — skipping.")
            continue
        ns = ns.iloc[0]

        # Cation solution SPE
        cs = sd[(sd["calculation_type"] == "SPE") & (sd["charge"] == 1)]
        if cs.empty:
            print(f"{structure_id}: missing cation solution SPE — skipping.")
            continue
        cs = cs.iloc[0]

        E_gas_neutral = ng["final_energy"]
        G_gas_neutral = ng["sum_electronic_thermal_free_energy"]
        E_gas_cation = cg["final_energy"]
        G_gas_cation = cg["sum_electronic_thermal_free_energy"]
        E_solv_neutral = ns["final_energy"]
        E_solv_cation = cs["final_energy"]

        # Solvated Gibbs free energies
        G_solv_neutral = E_solv_neutral + G_gas_neutral - E_gas_neutral
        G_solv_cation = E_solv_cation + G_gas_cation - E_gas_cation

        # Solvation free energies
        dG_solv_neutral = G_solv_neutral - G_gas_neutral
        dG_solv_cation = G_solv_cation - G_gas_cation

        # Gas-phase reaction free energy
        G_rxn = G_gas_cation - G_gas_neutral

        # Solution-phase reaction free energy
        dG_rxn_solv = (-dG_solv_cation + dG_solv_neutral) - G_rxn

        # Redox potential
        E0 = -(dG_rxn_solv * HARTREE_TO_EV / N_ELECTRONS) - SHE_OFFSET

        results.append({
            "structure_id": structure_id,
            "E0_vs_SHE_V": E0,
            "deltaG_redox_solv_Hartree": dG_rxn_solv,
            "deltaG_solv_neutral_kJmol": dG_solv_neutral * HARTREE_TO_KJMOL,
            "deltaG_solv_cation_kJmol": dG_solv_cation * HARTREE_TO_KJMOL,
        })

    if results:
        out = pd.DataFrame(results)
        out.to_csv(output_file, index=False)
        print(f"Saved {len(out)} results to {output_file}")
    else:
        print("No complete structures found.")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate 1 e⁻ oxidation potentials from DFT data.",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="CSV from extract_dft.py ion.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output CSV for redox potentials.")
    args = parser.parse_args()
    calculate_redox_potentials(args.input, args.output)


if __name__ == "__main__":
    main()

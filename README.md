# compchem-tools

Python tools for high-throughput redox-potential screening of organic redox-active molecules — RDKit utilities and automated Gaussian 16 DFT on HPC.

---

## RDKit tools (`smiles_tools/`)

### Functional group enumeration

Build combinatorial libraries by attaching functional groups to molecular backbones.

**CLI:**
```bash
# single substitution — one FG per position per molecule
python -m smiles_tools.enumerate \
    --backbones "c1ccc2c(c1)Nc3ccccc3S2" \
    --fg Cl "[OH]" "N(C)(C)" "C#N" \
    --mode single \
    --output derivatives.csv

# multi substitution — all combinations across positions
python -m smiles_tools.enumerate \
    --backbones "c1ccc2c(c1)Nc3ccccc3S2" \
    --fg Cl "[OH]" \
    --mode multi \
    --output multi_derivatives.csv

# targeted — specific FGs at specific atom indices
python -m smiles_tools.enumerate \
    --backbones "c1ccc2c(c1)Nc3ccccc3S2" \
    --mode targeted \
    --position-map '{"1": ["Cl", "[OH]"], "6": ["N(C)(C)"]}' \
    --output targeted.csv
```

**Python:**
```python
from smiles_tools import enumerate_library

df = enumerate_library(
    backbones=["c1ccc2c(c1)Nc3ccccc3S2"],
    functional_groups=["Cl", "[OH]", "N(C)(C)", "C#N"],
    mode="single",
)
```

### Reduced SMILES generation

Generates reduced (hydrogenated) SMILES from oxidised forms — supports phenazines, alloxazines and quinones.

**CLI:**
```bash
python -m smiles_tools.reduce_smiles \
    --input molecules.csv \
    --output reduced.csv \
    --type quinone \
    --smiles-col SMILES
```

**Python:**
```python
from smiles_tools import reduce_smiles
reduced = reduce_smiles("c1ccc2nc3ccccc3nc2c1", molecule_type="phenazine")
```

---

## DFT workflow (`dft_workflow/`)

Automated Gaussian 16 pipeline on HPC (PBS scheduler) for computing redox potentials.

### Pipeline

```
autodft launcher → conformer generation (ETKDG/MMFF94) → OPT+FREQ (M062X/6-31++G(d,p))
→ convergence check → resubmission (up to 3×) or B3LYP fallback
→ SPE (M062X/6-31++G(d,p), SMD water) → extract energies → compute E⁰
```

The launcher (`autodft_ions.py` or `autodft_pairs.py`) reads the input CSV, creates structure folders, and submits a PBS array job. Everything from conformer generation through to SPE runs automatically inside each PBS task. Once the jobs finish, run `extract_dft.py` then `process_redox_1e.py` or `process_redox_2e.py` to get redox potentials.

### Ion workflow (1e⁻ electron transfer)

Calculates E⁰ for **M → M⁺ + e⁻** from a CSV with a `SMILES` column.

```bash
# 1. Launch
python autodft_ions.py

# 2. Wait for PBS jobs to finish

# 3. Extract energies from Gaussian log files
python extract_dft.py ion --base-dir ./structures --output extracted.csv

# 4. Compute 1e⁻ redox potentials
python process_redox_1e.py --input extracted.csv --output redox_1e.csv
```

### Neutral-pair workflow (2e⁻ PCET)

Calculates E⁰ for **Ox + 2H⁺ + 2e⁻ → Red** from a CSV with `oxidised_smiles` and `reduced_smiles` columns. Each molecule gets two PBS tasks (odd index = oxidised, even = reduced).

```bash
# 1. Launch
python autodft_pairs.py

# 2. Wait for PBS jobs to finish

# 3. Extract energies from Gaussian log files
python extract_dft.py pair --base-dir ./structures --output extracted.csv

# 4. Compute 2e⁻ redox potentials
python process_redox_2e.py --input extracted.csv --output redox_2e.csv
```

See `example_files/` for sample input CSVs and a worked phenazine example.

---

## Installation

```bash
git clone https://github.com/datmdoan/compchem-tools.git
cd compchem-tools
pip install -r requirements.txt
```

Python ≥ 3.9. Gaussian 16 and a PBS scheduler are needed for the DFT workflow only.

---

## Citation

DOI to be added upon publication.

## Licence

MIT

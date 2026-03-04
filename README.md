# compchem-tools

Python tools for high-throughput redox-potential screening of organic redox-active molecules — cheminformatics utilities and Gaussian 16 DFT automation on HPC.

---

## Cheminformatics (`smiles_tools/`)

### Functional group enumeration

Generates libraries of functionalised derivatives from backbone SMILES.

```python
from smiles_tools import enumerate_library

df = enumerate_library(
    backbones=["c1ccc2c(c1)Nc3ccccc3S2"],
    functional_groups=["Cl", "[OH]", "N(C)(C)", "C#N"],
    mode="single",   # or "multi", "targeted"
)
```

Three enumeration modes:
- **single** — one FG per molecule per position
- **multi** — all combinations of FGs across positions
- **targeted** — specify which FGs go at which positions

```python
# targeted example: position 1 gets Cl or OH, position 6 gets N(C)(C)
df = enumerate_library(
    backbones=["c1ccc2c(c1)Nc3ccccc3S2"],
    mode="targeted",
    position_map={1: ["Cl", "[OH]"], 6: ["N(C)(C)"]},
)
```

Or via the CLI:
```bash
python -m smiles_tools.enumerate --backbones "c1ccc2c(c1)Nc3ccccc3S2" \
    --fg Cl "[OH]" "N(C)(C)" --mode single --output derivatives.csv
```

### Redox SMILES generation

Generate reduced (hydrogenated) SMILES for phenazines, alloxazines and quinones.

```python
from smiles_tools import reduce_smiles

reduced = reduce_smiles("c1ccc2nc3ccccc3nc2c1", molecule_type="phenazine")
```

```bash
python -m smiles_tools.reduce_smiles --input molecules.csv --output reduced.csv \
    --type quinone --smiles-col SMILES
```

---

## DFT workflow (`dft_workflow/`)

Automates Gaussian 16 OPT+FREQ → SPE calculations on the CX3 HPC cluster and computes redox potentials vs SHE.

### Ion workflow (1e⁻ ET)

E⁰ for **M⁺ + e⁻ → M**

```
molecules.csv → smiles_to_structure.py ion → autodft_ions.py/.pbs
→ extract_dft.py ion → process_redox_1e.py → redox_results.csv
```

### Neutral-pair workflow (2e⁻ PCET)

E⁰ for **Ox + 2H⁺ + 2e⁻ → Red**

```
molecules.csv → smiles_to_structure.py pair → autodft_pairs.py/.pbs
→ extract_dft.py pair → process_redox_2e.py → redox_results.csv
```

**Gaussian settings:** M062X/6-31++G(d,p), SMD water

---

## Installation

```bash
git clone https://github.com/your-username/compchem-tools.git
cd compchem-tools
pip install -r requirements.txt
```

Requires Python ≥ 3.9. Gaussian 16 and a PBS scheduler are needed for the DFT workflow only.

---

## Citation

DOI to be added upon publication.

## Licence

MIT

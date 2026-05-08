# SPA Benchmark Data

This repository contains the raw data and code for benchmarking the Separable Pair Ansatz (SPA) on hydrogen chains ($H_n$), alkanes, and small molecules.

## Installation

Clone this repository and install in developer mode:

```bash
git clone https://github.com/lily-barta/spa-benchmark.git
cd spa-benchmark
pip install -e .
```


## Generating Data 

### 1. Hydrogen Chains ($H_n$)

#### Dissociation Curve
Generate a CSV file for a given linear $H_n$ molecule.

```python
from spa_benchmark.main_Hn import run_dissociation
run_dissociation(n=4, max_iter=31, d_min=0.5, d_max=3.5, get_fci=True)
```

- `n` – Number of H atoms
- `max_iter` – Number of points along the dissociation curve
- `d_min`, `d_max` – Minimum and maximum interatomic distances (Å)

Output CSV file: `results_h{n}.csv`

Each row of the generated file contains:
- `distance` – Interatomic distance
- `spa` – Optimized SPA energy
- `fci` – Reference FCI ground-state energy (if `get_fci=True`)
- `fid` – Fidelity between SPA and FCI states (if `get_fci=True`)
- `hf` – HF energy (if `get_hf=True`)
- `hf_fid` - HF/FCI fidelity (if `get_hf=True`)
- `ccsdt` - CCSD(T) energy (if `get_ccsdt=True`)

```python
from spa_benchmark.plotting import plot_dissociation_curves, plot_dissociation_accuracy
plot_dissociation_curves(n=6)
plot_dissociation_accuracy()
```

#### Runtime and Accuracy Scaling
Generate data for increasing $H_n$ chain length for a fixed interatomic distance.

```python
from spa_benchmark.main_Hn import run_scaling
run_scaling(n_min=2, n_max=30, distance=1.0)
```
Output CSV files: 
- `results_vs_n.csv` – SPA/FCI/fidelity vs. n
- `runtime_vs_n.csv` – Total runtime and individual runtime contributions vs. n

```python
from spa_benchmark.plotting import plot_runtime_scaling, plot_accuracy_scaling, plot_spa_scaling
plot_runtime_scaling()
plot_accuracy_scaling()
plot_spa_scaling(quantity="fidelity")
```

#### SPA/FCI Fidelity Spectrum
Generate eigenenergies and corresponding SPA/FCI fidelity for a given linear $H_n$ molecule and interatomic distance.

```python
from spa_benchmark.main_Hn import run_fidelity_spectrum
run_fidelity_spectrum(n=6, distance=1.0, nroots=100)
```

- `n` – Number of H atoms
- `distance` - Interatomic distance
- `nroots` - Number of FCI eigenstates to compute

Output CSV file: `h{n}_spectrum_{distance}.csv`

```python
from spa_benchmark.plotting import plot_fidelity_spectrum
plot_fidelity_spectrum(filename="h6_spectrum_0.9.csv")
```

---

### 2. Alkanes and Small Molecules

Single-point SPA benchmarks for molecules defined by `.xyz` geometry files. Geometry files are stored in the `geometries/` directory at the root of the repository.

```python
from spa_benchmark.main import generate_data_point
generate_data_point(
    geometry="geometries/ch4.xyz",
    get_ccsdt=True,
    get_fci=True,
    verbose=True,
    filename="alkanes.csv",
    filename_t="alkanes_t.csv"
)
```

To benchmark a series of molecules, iterate over geometry files:
```python
# Alkane series
molecules = ['ch4', 'c2h6', 'c3h8', 'c4h10', 'c5h12', 'c6h14', 'c7h16', 'c8h18']
filename, filename_t = "alkanes.csv", "alkanes_t.csv"

# Small molecules
# molecules = ['h2o', 'nh3', 'ch4', 'hf', 'h2co', 'c2h4', 'c2h2', 'ch3oh']
# filename, filename_t = "small_molecs.csv", "small_molecs_t.csv"

for mol in molecules:
    generate_data_point(
        geometry=f"geometries/{mol}.xyz",
        get_ccsdt=True,
        get_fci=True,
        verbose=True,
        filename=filename,
        filename_t=filename_t,
    )
```

Each call appends one row to the output CSVs:
- `{filename}.csv` – SPA, FCI, and CCSD(T) energies and SPA/FCI fidelity per molecule
- `{filename_t}.csv` – Runtime contributions

---

All output CSV files are written to the current working directory.
Pre-generated example data is available in `data_hn/`.

## Notes on Computational Cost

- **FCI calculations** become computationally demanding for larger systems.  
  For example, a single-point FCI calculation for $H_{14}$ may take approximately **30 minutes to 1 hour**, depending on hardware.  
  For larger chains, FCI becomes prohibitively expensive and is therefore not recommended.
  In the provided scaling script, FCI calculations are automatically disabled for larger systems to avoid excessive runtimes.

- **Near-degeneracies at large interatomic distances:**
  As the $H_n$ chains dissociate, the FCI ground state becomes degenerate.  
  To obtain a meaningful fidelity in this regime, multiple FCI roots must be included.  
  The required number of roots grows rapidly with system size, and it becomes prohibitively expensive beyond $H_{10}$.  
  For this reason, fidelity calculations at large separations ($d > 2.5$ Å) is disabled for larger systems.


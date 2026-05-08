# SPA Benchmark Data

This repository contains the raw data and code for benchmarking the Separable Pair Ansatz (SPA) method on hydrogen chains ($H_n$).

## Installation

Clone this repository and install in developer mode:

```bash
git clone https://github.com/lily-barta/spa-benchmark.git
cd spa-benchmark
pip install -e .
```


## Generating Data 

### 1. Dissociation Curve
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

### 2. Runtime and Accuracy Scaling
Generate data for increasing $H_n$ chain length for a fixed interatomic distance.

```python
from spa_benchmark.main_Hn import run_scaling

run_scaling(n_min=2, n_max=30, distance=1.0)
```
Output CSV files: 
- `results_vs_n.csv` – SPA/FCI/fidelity vs. n
- `runtime_vs_n.csv` – Total runtime and individual runtime contributions vs. n

### 3. SPA/FCI Fidelity Spectrum
Generate a CSV file containing eigenenergies and corresponding SPA/FCI fidelity for a given linear $H_n$ molecule and interatomic distance.

```python
from spa_benchmark.main_Hn import run_fidelity_spectrum

run_fidelity_spectrum(n=6, distance=1.0, nroots=100)
```

- `n` – Number of H atoms
- `distance` - Interatomic distance
- `nroots` - Number of FCI eigenstates to compute

Output CSV file: `h{n}_spectrum_{distance}.csv`

Output CSV files are written to the current working directory.
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

## Reproducing Plots

Plots can be generated directly from the CSV files using the provided plotting utilities.

### Example

```python
from plotting import plot_timing, plot_accuracy, plot_dissociation,

plot_timing()
plot_accuracy()
plot_dissociation(n=6)
```

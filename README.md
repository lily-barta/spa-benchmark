# SPA Benchmark Data

This repository contains the raw data and code for benchmarking the Separable Pair Ansatz (SPA) method on hydrogen chains ($H_n$).

## Installation

Clone this repository and install in developer mode:

```bash
git clone https://github.com/lily-barta/spa-benchmark.git
cd spa-benchmark
pip install -e .
```

You will also need to install Sunrise quantum chemistry package
pip install git+https://github.com/tequilahub/sunrise.git@devel
```

## Generating Data 

### 1. Dissociation Curve
Generate a CSV file for a given linear $H_n$ molecule.

```python
from spa_benchmark import run_dissociation

run_dissociation(n=4, max_iter=31, d_min=0.5, d_max=3.5, get_fci=True, get_var=True)
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
- `var` – Energy variance of the SPA state (if `get_var=True`)

### 2. Single-Point Scaling
Generate data for increasing $H_n$ chain length at a fixed interatomic distance.

```python
from spa_benchmark import run_scaling

run_scaling(n_min=2, n_max=30, distance=1.0)
```
Output CSV files: 
- `results_vs_n.csv` – SPA/FCI/fidelity/variance vs. n
- `timing_vs_n.csv` – Total runtime and individual timing contributions vs. n

Output CSV files are written to the current working directory.
Pre-generated example data is available in `data_hn/`.

## Notes on Computational Cost

- **FCI calculations** become computationally demanding for larger systems.  
  For example, a single-point FCI calculation for $H_{14}$ may take approximately **30 minutes to 1 hour**, depending on hardware.  
  For larger chains, FCI becomes prohibitively expensive and is therefore not recommended.

- **Energy variance calculations** are also computationally intensive.  
  Even for $H_{10}$, the variance evaluation is very slow (**3/4 hours**) and is therefore also not recommended beyond this system size.

- In the provided scaling script, FCI and variance calculations are automatically disabled for larger systems to avoid excessive runtimes.

- **Near-degeneracies at large interatomic distances:**  
As the $H_n$ chains dissociate, the FCI ground state becomes degenerate.  
To obtain a meaningful fidelity in this regime, multiple FCI roots must be included.  
The required number of roots grows rapidly with system size, and for $H_{10}$ it becomes prohibitively expensive.  
For this reason, fidelity calculations at large separations ($d > 2.5$ Å) are not recommended for larger systems.

## Reproducing Plots

Plots can be generated directly from the CSV files using the provided plotting utilities.

### Example

```python
from plotting import plot_timing, plot_accuracy, plot_dissociation,

plot_timing()
plot_accuracy()
plot_dissociation(n=6)
```

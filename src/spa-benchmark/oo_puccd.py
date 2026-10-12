import tequila as tq
import sunrise as sun
from sunrise.measurement import get_hcb_part
import numpy
import csv
import time
import sys, os
# print(sys.getrecursionlimit())
sys.setrecursionlimit(30000)

def generate_geometry(n, distance):
    geom = ""
    for i in range(n):
        geom += f"H 0.0 0.0 {i * distance}\n"
    return geom

def get_edges_and_guess(n):
    edges = [(i, i + 1) for i in range(0, n, 2)]
    guess = numpy.eye(n)
    for edge in edges:
        guess[edge[0]][edge[1]] = 1
        guess[edge[1]][edge[0]] = -1
    return edges, guess.T

def run_spa(n, distance):
    
    geometry = generate_geometry(n, distance)
    mol = sun.Molecule(geometry=geometry, basis_set="sto-3g", transformation="reordered-jordan-wigner").use_native_orbitals()

    # Reorder spa orbitals to match upccd's convention
    # Put the first orbital of each adjacent pair in UpCCD's occupied block.
    _, guess = get_edges_and_guess(n)
    order = list(range(0, n, 2)) + list(range(1, n, 2))
    mol = mol.transform_orbitals(numpy.eye(n)[:, order])
    guess = guess[numpy.ix_(order, order)]
    edges = [(i, i + n // 2) for i in range(n // 2)]

    # Get spa optimised orbitals
    spa_circuit = mol.make_spa_ansatz(edges=edges, hcb=True)
    solver = sun.SPAFP.SPASolver(decompose=True, grouping=edges)
    spa = sun.optimize_orbitals(
        mol, spa_circuit, use_hcb=True, vqe_solver=solver,
        initial_guess=guess, silent=True)
    return spa.molecule


def run_oo_puccd_fci_scaling(n_min=2, n_max=10, distance=1.0,
                            filename="results_oo_puccd_fci_vs_n.csv",
                            verbose=False):
    
    optimizer_options = {"method": "BFGS", "initial_values": 0.0, "silent": True}
    orbital_options = {"max_cycle_macro": 50, "max_cycle_micro": 3}
    columns = ["n", "oo_puccd", "fci", "fid", "oo_puccd_t", "vqe_t"]
    rows = []

    with open(filename, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        for n in range(n_min, n_max + 1, 2):
            oo_start = time.perf_counter()

            mol = run_spa(n, distance)
            U = mol.make_ansatz(name="HCB-UpCCD")
            opt = sun.optimize_orbitals(
                mol, U, use_hcb=True, return_mcscf=True, silent=True,
                vqe_solver_arguments={"optimizer_arguments": dict(
                    optimizer_options)}, pyscf_arguments=orbital_options.copy(),
            )
            oo_puccd_t = time.perf_counter() - oo_start
            mol = opt.molecule

            H = mol.make_hardcore_boson_hamiltonian()
            E = tq.ExpectationValue(H=H, U=U)
            start = time.perf_counter()
            vqe = tq.minimize(E, **optimizer_options)
            vqe_t = time.perf_counter() - start

            wfn = tq.simulate(U, variables=vqe.variables)

            fci_energie, fci_wfn = mol.compute_energy(
                "fci", get_wfn=True, use_hcb=True)
            
            fidelity = abs(wfn.inner(fci_wfn))**2

            row = {"n": n, "oo_puccd": vqe.energy, "fci": fci_energie,
                   "fid": fidelity, "oo_puccd_t": oo_puccd_t, "vqe_t": vqe_t}
            rows.append(row)
            writer.writerow(row)
            file.flush()
            if verbose:
                print(f"H{n}: OO-UpCCD = {row['oo_puccd']:.10f}, "
                      f"FCI = {row['fci']:.10f}, fidelity = {fidelity:.8f}; "
                      f"OO total = {oo_puccd_t:.3f}s, VQE = {vqe_t:.3f}s, ")

    return rows


if __name__ == "__main__":
    rows = run_oo_puccd_fci_scaling(
        n_min=2, n_max=30, distance=1.0,
        filename="test.csv",
        verbose=True,
    )

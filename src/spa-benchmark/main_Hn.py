import sunrise as sun
import numpy
import csv
import time
import sys
sys.setrecursionlimit(30000)

def generate_geometry(n, iter, max_iter, d_min=0.5, d_max=4.0):
    if max_iter == 1:
        R = d_min
    else:
        R = d_min + (d_max - d_min) * iter / (max_iter - 1)

    geom = ""
    for i in range(n):
        geom += f"H 0.0 0.0 {i * R}\n"
    return geom, R

def get_edges_and_guess(n):
    edges = [(i, i + 1) for i in range(0, n, 2)]
    guess = numpy.eye(n)
    for edge in edges:
        guess[edge[0]][edge[1]] = 1
        guess[edge[1]][edge[0]] = -1
    return edges, guess.T

def generate_data_point(n, iter, max_iter, d_min=0.5, d_max=4.0, nroots=1, get_fci=False, get_hf=False, get_ccsdt=False, verbose=False):
    start=time.time()

    # Build tequila/sunrise molecule
    geometry, distance = generate_geometry(n, iter, max_iter, d_min=d_min, d_max=d_max)
    mol = sun.Molecule(geometry=geometry,basis_set='sto-3g',transformation='reordered-jordan-wigner')

    if verbose:
        print(f"\nIteration : {iter + 1} / {max_iter}")
        print(f"Interatomic distance = {distance:.5f}")
        
    results = {"n": n, "distance": f"{distance:.3f}"}
    results_t = {"n": n}

    if get_hf:
        # Compute HF energy
        hf = mol.compute_energy("hf")
        results["hf"] = f"{hf:.10f}"

        # Compute HF/FCI fidelity
        U_hf = mol.prepare_reference()
        wfn_hf = sun.simulate(U_hf)
        nroots_map = {4: 6, 6: 20, 8: 70, 10: 252} # should be 252 for h10
        if distance >= 2.5 and n > 10:
            print(f"\n!!! Warning !!! \nFCI for H{n} at distance {distance:.2f} Å is degenerate. \n"
                "Fidelity requires many FCI roots and is not computed due to high cost.\n")
        elif distance >= 2.5 and n in nroots_map:
            nroots = nroots_map[n]
        ci0 = wfn_hf if n > 6 else None
        fci, wfn_fci = mol.compute_energy("fci", get_wfn=True, nroots=nroots, ci0=ci0)
        fci0 = fci if nroots == 1 else fci[0]  

        if nroots == 1:
            fidelity = abs(wfn_hf.inner(wfn_fci))**2
        else:
            fidelity = 0.0
            for i in range(nroots):
                if abs(fci[i]-fci0) <= 0.0016: 
                    fidelity += abs(wfn_hf.inner(wfn_fci[i]))**2
        results["hf_fid"] = f"{fidelity:.6f}"

        if verbose:
            print(f"HF      : {hf:.10f}")
            print(f"FCI      : {fci0:.10f}")
            print(f"HF error    : {hf-fci0:.10f}")
            print(f"fidelity : {fidelity:.6f}")
    
    if get_ccsdt:
        # Compute CCSD(T) energy
        ccsdt = mol.compute_energy('CCSD(T)')
        results["ccsdt"] = f"{ccsdt:.10f}"
        if verbose:
            print(f"CCSD(T)  : {ccsdt:.10f}")

    mol = mol.use_native_orbitals()

    # Optimise orbitals
    oo_start=time.time()
    edges, guess = get_edges_and_guess(n)
    opt = sun.SPAFP.run_spa(mol=mol, edges=edges, initial_guess=guess)
    mol = opt.molecule
    results_t["orb_opt"] = f"{time.time() - oo_start:.6f}"
    if verbose:
        print(f"Orb opt took {time.time() - oo_start:.6f}s") 

    # Fast HCB-SPA VQE
    spa_start=time.time()
    U = mol.make_spa_ansatz(edges=edges, hcb=True)
    H = mol.make_hardcore_boson_hamiltonian()
    grouping = sun.SPAFP.make_decomposed_clusters(U)
    vqe_solver = sun.SPAFP.SPASolver(decompose=True,grouping=grouping)
    spa = vqe_solver(H=H, circuit=U, molecule=mol)
    results["spa"] = f"{spa.energy:.10f}"
    results_t["spa"] = f"{time.time() - spa_start:.6f}"
    if verbose:
        print(f"VQE SPA  : {spa.energy:.10f}")
        print(f"SPA energy took {time.time() - spa_start:.6f}s") 

    if get_fci:
        # Get HCB-SPA wavefunction
        spa_start = time.time()
        wfn_spa_hcb = sun.simulate(U, variables=spa.variables)
        results_t["spa_wfn"] = f"{time.time() - spa_start:.6f}"
        if verbose:
            print(f"SPA wfn took {time.time() - spa_start:.6f}s") 
    
        # Compute FCI energy and wavefunction
        fci_start = time.time()
        nroots_map = {4: 6, 6: 20, 8: 70, 10: 252}
        if distance >= 2.5 and n > 10:
            print(f"\n!!! Warning !!! \nFCI for H{n} at distance {distance:.2f} Å is degenerate. \n"
                "Fidelity requires many FCI roots and is not computed due to high cost.\n")
        elif distance >= 2.5 and n in nroots_map:
            nroots = nroots_map[n]
        ci0 = wfn_spa_hcb if n > 6 and distance > 1.0 else None
        fci, wfn_fci = mol.compute_energy("fci", get_wfn=True, nroots=nroots, ci0=ci0, use_hcb=True)
        fci0 = fci if nroots == 1 else fci[0]     
        results["fci"] = f"{fci0:.10f}"
        results_t["fci"] = f"{time.time() - fci_start:.6f}"
        if verbose:
            print(f"FCI      : {fci0:.10f}")
            print(f"Spa error    : {spa.energy-fci0:.10f}")
            print(f"FCI energy+wfn took {time.time() - fci_start:.6f}s")
    
        # Compute fidelity (SPA/FCI overlap)
        fid_start = time.time()
        if nroots == 1:
            fidelity = abs(wfn_spa_hcb.inner(wfn_fci))**2
        else:
            fidelity = 0.0
            for i in range(nroots):
                if abs(fci[i]-fci[0]) < 0.0016: 
                    fidelity += abs(wfn_spa_hcb.inner(wfn_fci[i]))**2
        results["fid"] = f"{fidelity:.6f}"
        results_t["fid"] = f"{time.time() - fid_start:.6f}"
        if verbose:
            print(f"fidelity : {fidelity:.6f}")
            print(f"Fidelity took {time.time() - fid_start:.6f}s")

    results_t["total"] = f"{time.time() - start:.6f}"
    if verbose:
            print(f"Data point took {time.time()-start:.6f}s")
    return results, results_t


def run_dissociation(n, max_iter, d_min=0.5, d_max=4.0, nroots=1, filename=None, get_fci=False, get_hf=False, get_ccsdt=False, verbose=False):
    """
    :param n: Number of H atoms in H_n chain.
    :param max_iter: Total number of dissociation points.
    :param d_min: Minimum interatomic distance (in angstroms).
    :param d_max: Maximum interatomic distance (in angstroms).
    :param nroots: Number of FCI eigenstates used as reference. Use nroots > 1 in the presence of near-degeneracies.
    :param filename: Name of the output CSV file. If None, default is f"results_h{n}.csv".
    :param get_fci: If True, compute FCI energy and SPA/FCI fidelity.
    :param get_hf: If True, compute HF energy and HF/FCI fidelity.
    :param get_ccsdt: If True, compute CCSD(T) energy.
    :param verbose: If True, prints results
    
    :return: Generates a CSV file with one row per dissociation point.
    """
    if filename is None:
        filename = f"results_h{n}.csv"
    
    with open(filename, "w", newline="") as file:
        columns = ["distance", "spa"]
        if get_fci:
            columns += ["fci", "fid"]
        if get_hf:
            columns += ["hf", "hf_fid"]
        if get_ccsdt:
            columns.append("ccsdt")
        writer = csv.DictWriter(file, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for iter in range(max_iter):
            data, data_t = generate_data_point(n, iter, max_iter, d_min=d_min, d_max=d_max, nroots=nroots, get_fci=get_fci, get_ccsdt=get_ccsdt, get_hf=get_hf, verbose=verbose)
            writer.writerow(data)

def run_single_point(n, distance, nroots=1, get_fci=False, get_hf=False, get_ccsdt=False, verbose=True):
    """
    :param n: Number of H atoms in H_n chain.
    :param distance: Interatomic distance (in angstroms).
    :param nroots: Number of FCI eigenstates used as reference. Use nroots > 1 in the presence of near-degeneracies. 
    :param get_fci: If True, compute FCI energy and SPA/FCI fidelity.
    :param get_hf: If True, compute HF energy and HF/FCI fidelity.
    :param get_ccsdt: If True, compute CCSD(T) energy.
    :param verbose: If True, prints results

    :return: 2 dictionaries, one containing the computed energies/fidelities, and one containing the computed runtimes.
    """
    return generate_data_point(n=n, iter=0, max_iter=1, d_min=distance, nroots=nroots, get_fci=get_fci, get_ccsdt=get_ccsdt, get_hf=get_hf, verbose=verbose)

def run_scaling(n_min=2, n_max=10, distance=1.0, nroots=1, filename_t="timing_vs_n.csv", filename="results_vs_n.csv", verbose=False):
    """
    :param n_min: Minimum number of H atoms in H_n chain.
    :param n_max: Maximum number of H atoms in H_n chain.
    :param distance: Interatomic distance (in angstroms).
    :param nroots: Number of FCI eigenstates used as reference. Use nroots > 1 in the presence of near-degeneracies.
    :param filename_t: Name of the output CSV file containing the runtimes (default "timing_vs_n.csv").
    :param filename: Name of the output CSV file for SPA/FCI/fidelity results (default "results_vs_n.csv"). 
    
    :return Generates two CSV files, one row per H_n chain: one for results and one for runtimes.
    """
    columns = ["n", "spa", "fci", "fid"]
    columns_t = ["n", "orb_opt", "spa", "spa_wfn", "fci", "fid", "total"]

    with open(filename, "w", newline="") as file, \
         open(filename_t, "w", newline="") as file_t:

        writer = csv.DictWriter(file, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer_t = csv.DictWriter(file_t, fieldnames=columns_t, extrasaction='ignore')
        writer_t.writeheader()

        for n in range(n_min, n_max + 1, 2):
            get_fci = True if n <= 14 else False
            data, data_t = run_single_point(n=n, distance=distance, nroots=nroots, get_fci=get_fci, verbose=verbose)
            writer.writerow(data)
            writer_t.writerow(data_t)


if __name__ == "__main__":

    # run_dissociation(n=4, max_iter=11, d_max=3.0, d_min=1.0, verbose=True, get_fci=True)
    run_scaling(n_max=30)

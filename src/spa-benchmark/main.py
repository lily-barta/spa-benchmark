import sunrise as sun
import numpy 
import csv
import time
import os


def generate_data_point(geometry, get_hf=False, get_ccsdt=False, get_fci=False, verbose=True, filename=None, filename_t=None, save_to_csv=True):
    """
    :param geometry: Geometry string or path to .xyz file.
    :param get_hf: Bool, if True compute HF energy and HF/FCI fidelity.
    :param get_ccsdt: Bool, if True compute CCSD(T) energy.
    :param get_fci: Bool, if True compute FCI energy and SPA/FCI fidelity.
    :param save_to_csv: Bool, if True save results in a .csv file.
    :param filename: Name of the generated .csv file containing energies and fidelities. Default is results_name.
    :param filename_t: Name of the generated .csv file containing runtimes. Default is results_name_t.
    """
    start = time.time()

    # Build tequila/sunrise molecule
    mol = sun.Molecule(geometry=geometry,basis_set='sto-3g',transformation='reordered-jordan-wigner')
    if os.path.isfile(geometry):
        name = os.path.splitext(os.path.basename(geometry))[0] 
    else:
        name = mol.parameters.name
    results = {"name": name}
    results_t = {"name": name}
    if verbose:
        print(f"Results for {name}:")

    if get_hf:
        # Compute HF energy
        hf = mol.compute_energy("hf")
        results["hf"] = f"{hf:.10f}"
        # Compute HF/FCI fidelity
        U_hf = mol.prepare_reference()
        wfn_hf = sun.simulate(U_hf)
        fci, fci_wfn = mol.compute_energy("fci", get_wfn=True)
        fidelity = abs(wfn_hf.inner(fci_wfn))**2
        results["hf_fid"] = f"{fidelity:.6f}"
        if verbose:
            print(f"HF      : {hf:.10f}")
            print(f"HF error: {hf-fci:.10f}")
            print(f"HF/FCI fidelity : {fidelity:.6f}")
        
    if get_ccsdt:
        # Compute CCSD(T) energy
        ccsdt_start = time.time()
        ccsdt = mol.compute_energy('CCSD(T)')
        results["ccsdt"] = f"{ccsdt:.10f}"
        results_t["ccsdt"] = f"{time.time() - ccsdt_start:.6f}"
        if verbose:
            print(f"CCSD(T)  : {ccsdt:.10f}")
            print(f"CCSD(T) took {time.time() - ccsdt_start:.6f}s")
    
    # Localise orbitals
    ol_start = time.time()
    mol, edges = sun.CLPO.generate_CLPO_molecule_edges(mol)
    results_t["orb_loc"] = f"{time.time() - ol_start:.6f}"
    if verbose:
        print(f"Orb loc took {time.time() - ol_start:.6f}s")

    # Optimise orbitals
    oo_start = time.time()
    guess = numpy.eye(mol.n_orbitals)
    opt = sun.SPAFP.run_spa(mol=mol, edges=edges, initial_guess=guess)
    mol = opt.molecule
    results_t["orb_opt"] = f"{time.time() - oo_start:.6f}"
    if verbose:
        print(f"Orb opt took {time.time() - oo_start:.6f}s")

    # Fast HCB-SPA VQE
    U = mol.make_spa_ansatz(edges=edges, hcb=True)
    H = mol.make_hardcore_boson_hamiltonian()
    grouping = sun.SPAFP.make_decomposed_clusters(U)
    spa_start=time.time()
    vqe_solver = sun.SPAFP.SPASolver(decompose=True, grouping=grouping)
    spa = vqe_solver(H=H, circuit=U, molecule=mol)
    results["spa"] = f"{spa.energy:.10f}"
    results_t["spa"] = f"{time.time() - spa_start:.6f}"
    if verbose:
        print(f"SPA energy took {time.time() - spa_start:.6f}s")
        print(f"VQE SPA  : {spa.energy:.10f}")

    if get_fci:
        if mol.n_orbitals > 14:
            print(f"get_fci={get_fci}, but molecule is too big (n_orbitals = {mol.n_orbitals} > 14). Skipping FCI calculation.")
        else:
            # Get HCB-SPA wavefunction
            spa_start = time.time()
            wfn_spa = sun.simulate(U, variables=spa.variables)
            results_t["spa_wfn"] = f"{time.time() - spa_start:.6f}"
            if verbose:
                print(f"SPA wfn took {time.time() - spa_start:.6f}s")

            # Compute FCI energy and wavefunction
            fci_start = time.time()
            fci, fci_wfn = mol.compute_energy("fci", get_wfn=True, use_hcb=True)
            results["fci"] = f"{fci:.10f}"
            results_t["fci"] = f"{time.time() - fci_start:.6f}"
            if verbose:
                print(f"FCI      : {fci:.10f}")
                print(f"SPA error: {spa.energy-fci:.10f}")
                print(f"FCI energy+wfn took {time.time()-fci_start:.6f}s")

            # Compute fidelity (SPA/FCI overlap)
            fid_start = time.time()
            fidelity = abs(wfn_spa.inner(fci_wfn))**2
            results["fid"] = f"{fidelity:.6f}"
            results_t["fid"] = f"{time.time() - fid_start:.6f}"
            if verbose:
                print(f"fidelity : {fidelity:.6f}")
                print(f"Fidelity took {time.time()-fid_start:.6f}s")

    results_t["total"] = f"{time.time() - start:.6f}"
    if verbose:
        print(f"Data point took {time.time()-start:.6f}s")

    # Store results
    if save_to_csv:
        if filename is None:
            filename = f"results_{name}.csv"
        if filename_t is None:
            filename_t = f"results_{name}_t.csv"
        
        file_exists = os.path.exists(filename)
        with open(filename, mode="a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(results)

        file_exists = os.path.exists(filename_t)
        with open(filename_t, mode="a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results_t.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(results_t)

if __name__ == "__main__":
    molecules = ['ch4', 'c2h6', 'c3h8', 'c4h10', 'c5h12', 'c6h14', 'c7h16', 'c8h18']
    filename = "alkanes.csv"
    filename_t = "alkanes_t.csv"

    # molecules = ['h2o', 'nh3', 'ch4', 'hf', 'h2co', 'c2h4', 'c2h2', 'ch3oh']
    # filename = "small_molecs.csv"
    # filename_t = "small_molecs_t.csv"

    for mol in molecules:
        generate_data_point(
            geometry=f"geometries/{mol}.xyz",
            get_ccsdt=True,
            get_fci=True,
            verbose=True,
            filename=filename,
            filename_t=filename_t
        )

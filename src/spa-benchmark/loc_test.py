import sunrise as sun
import tequila as tq
import numpy 
import csv
import time
import os
import sys
sys.setrecursionlimit(10000000)


def _match_spa_edges(occupied, virtual, exchange):
    """Maximize total exchange with at most one virtual per occupied orbital.

    Rows/columns of exchange follow occupied/virtual, respectively. Unmatched
    occupied orbitals remain singletons; every virtual must be assigned.
    """
    from scipy.optimize import linear_sum_assignment

    if len(virtual) > len(occupied):
        raise ValueError(
            "One-to-one SPA pairing cannot cover more virtual than occupied "
            "orbitals. Use a smaller active space or supply explicit edges "
            "allowing multiple virtuals per occupied orbital."
        )
    edges = [[i] for i in occupied]
    if virtual:
        rows, columns = linear_sum_assignment(-numpy.asarray(exchange))
        for row, column in zip(rows, columns):
            edges[row].append(virtual[column])
    return edges


def localize_spa_orbitals(mol, localization="clpo", edges=None):
    """Localize orbitals and return disjoint SPA edges in active-orbital indices.

    Boys/PM rotate active occupied and virtual spaces separately, keeping frozen
    orbitals fixed. With no explicit edges, maximize the total exchange (ia|ai)
    using one-to-one occupied/virtual matching. Unmatched occupied orbitals
    remain singletons (e.g. lone pairs); their chemical identity is not imposed.
    Pairing is separate from the Boys/PM localization criterion.
    """
    method = localization.lower().replace("_", "-").replace(" ", "-")
    if method in ("pm", "pipek", "pipekmezey"):
        method = "pipek-mezey"
    if method not in ("clpo", "boys", "pipek-mezey"):
        raise ValueError("localization must be 'clpo', 'boys', or 'pipek-mezey'")
    if method == "clpo":
        return sun.CLPO.generate_CLPO_molecule_edges(mol, edges=edges)

    from pyscf import lo

    occupied = [orb.idx for orb in mol.reference_orbitals]
    virtual = [i for i in range(mol.n_orbitals) if i not in occupied]
    active = [orb.idx_total for orb in mol.integral_manager.active_orbitals]
    coefficients = mol.integral_manager.orbital_coefficients[:, active]
    localized = coefficients.copy()
    for block in (occupied, virtual):
        if len(block) > 1:
            if method == "boys":
                localizer = lo.Boys(mol.pyscf_molecule, coefficients[:, block])
            else:
                localizer = lo.PM(mol.pyscf_molecule, coefficients[:, block])
            localizer.verbose = 0
            localized[:, block] = localizer.kernel()
    # transform_orbitals expects coefficients relative to the current basis,
    # whereas PySCF returns coefficients in the AO basis.
    overlap = mol.integral_manager.overlap_integrals
    rotation = coefficients.T @ overlap @ localized
    localized_mol = mol.transform_orbitals(rotation)
    if edges is None:
        _, _, integrals = localized_mol.get_integrals(two_body_ordering="chem")
        exchange = numpy.array([[integrals.elems[i, a, a, i] for a in virtual]
                                for i in occupied])
        edges = _match_spa_edges(occupied, virtual, exchange)
    else:
        edges = [list(edge) for edge in edges]
    return localized_mol, edges


def generate_data_point(geometry, verbose=True, filename=None, save_to_csv=True,
                        localization="clpo", edges=None, plot_orbitals=False):
    """
    :param geometry: Geometry string or path to .xyz file.
    :param save_to_csv: Bool, if True save results in a .csv file.
    :param filename: CSV; default results_<name>_<localization>.csv.
    :param localization: 'clpo' (default), 'boys', or 'pipek-mezey' ('pm').
    :param edges: Optional SPA groups in localized active-orbital indices.
        Boys/PM default to global one-to-one exchange matching, with unmatched
        occupied orbitals left as singletons; see localize_spa_orbitals.
    :param plot_orbitals: Write orbital plots before and after optimization.
    """
    start = time.time()

    mol = tq.Molecule(geometry=geometry,basis_set='sto-3g',transformation='reordered-jordan-wigner',backend='pyscf')
    
    if os.path.isfile(geometry):
        name = os.path.splitext(os.path.basename(geometry))[0] 
    else:
        name = mol.parameters.name

    localization = localization.lower().replace("_", "-").replace(" ", "-")
    if localization in ("pm", "pipek", "pipekmezey"):
        localization = "pipek-mezey"
    results = {"name": name, "localization": localization}
    if verbose:
        print(f"Results for {name}:")

    
    # Localise orbitals
    ol_start = time.time()
    mol, edges = localize_spa_orbitals(mol, localization, edges)
    ol_time = time.time() - ol_start

    # edges = [(0,16),(1,15),(2,21),(3,13),(4,11),(5,18),(6,17),(7,12),(8,14),(9,19),(10,20)] #c4h6 boys
    # edges = [(0,15),(1,22),(2,16),(3,14),(4,25),(5,19),(6,17),(7,18),(8,20),(9,24),(10,13),(11,21),(12,23)] #c4h10 boys
    # edges = [(0,15),(1,20),(2,16),(3,21),(4,14),(5,23),(6,19),(7,18),(8,22),(9,24),(10,13),(11,25),(12,17)] #c4h10 pm

    if plot_orbitals:
        sun.plot_MO(molecule=mol, filename=localization)
    if verbose:
        print(f"Orb loc took {time.time() - ol_start:.6f}s")

    # Compute non opt spa
    U = mol.make_spa_ansatz(edges=edges, hcb=True)
    H = mol.make_hardcore_boson_hamiltonian()
    grouping = sun.SPAFP.make_decomposed_clusters(U)
    spa_start=time.time()
    vqe_solver = sun.SPAFP.SPASolver(decompose=True, grouping=grouping)
    spa = vqe_solver(H=H, circuit=U, molecule=mol)
    results["spa_localized"] = f"{spa.energy:.10f}"
    if verbose:
        print(f"SPA energy took {time.time() - spa_start:.6f}s")
        print(f"VQE SPA  : {spa.energy:.10f}")

    # Optimise orbitals and compute spa
    oo_start = time.time()
    guess = numpy.eye(mol.n_orbitals)
    opt = sun.SPAFP.run_spa(mol=mol, edges=edges, initial_guess=guess,pyscf_arguments={"max_cycle_macro": 50, "max_cycle_micro": 3},return_mcscf=True)
    mol = opt.molecule
    macroiterations = (opt.mcscf_local_data or {}).get("imacro", "unavailable")
    results["spa_opt"] = f"{opt.energy:.10f}"
    results["loc_time"] = f"{ol_time:.6f}"
    results["opt_time"] = f"{time.time() - oo_start:.6f}"
    results["converged"] = opt.mcscf_object.converged
    results["macroiterations"] = macroiterations
    results["edges"] = edges
    if plot_orbitals:
        sun.plot_MO(molecule=mol, filename=f"{localization}_opt")
    if verbose:
        print(f"OO macroiterations: {macroiterations}; converged: {opt.mcscf_object.converged}")
        print(f"Orb opt took {time.time() - oo_start:.6f}s")
        print(f"Data point took {time.time()-start:.6f}s")

    if save_to_csv:
        if filename is None:
            filename = f"results_{name}.csv"
    
        file_exists = os.path.exists(filename) and os.path.getsize(filename) > 0
        if file_exists:
            with open(filename, newline="") as f:
                if next(csv.reader(f)) != list(results):
                    raise ValueError(f"CSV columns do not match: {filename}; use a new file")
        with open(filename, mode="a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(results)


if __name__ == "__main__":

    edges = [{0, 15}, {1, 22}, {16, 2}, {3, 14}, {25, 4}, {19, 5}, {17, 6}, {18, 7}, {8, 20}, {24, 9}, {10, 13}, {11, 21}, {12, 23}]
    result = generate_data_point(
        geometry="../geometries/c4h10.xyz",
        localization="boys",
        plot_orbitals=False,
        edges=edges
    )
    
    # for method in ["clpo","boys","pm"]:
    #     result = generate_data_point(
    #         geometry="../geometries/alkenes/c4h6.xyz",
    #         localization=method,
    #         plot_orbitals=False
    #     )
    # for method in ["clpo","boys","pm"]:
    #     result = generate_data_point(
    #         geometry="../geometries/c4h10.xyz",
    #         localization=method,
    #         plot_orbitals=False
    #     )

    

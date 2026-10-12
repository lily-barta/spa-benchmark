from pyblock2.driver.core import DMRGDriver, SymmetryTypes
import tequila as tq
import csv, os

class ResultsTable:
    def __init__(self, filename="results_dmrg.csv"):
        self.filename = filename
        self.rows = []
        # Load existing data if file existsP    
        if os.path.exists(self.filename):
            with open(self.filename, mode='r', newline='') as file:
                reader = csv.DictReader(file)
                self.rows = list(reader)

    def add_entry(self, name, energy):
        """Add a name and energy to the table and save to CSV immediately."""
        new_entry = {"name": name, "dmrg": energy}
        self.rows.append(new_entry)
        self._save_to_csv()

    def _save_to_csv(self):
        """Save current data to CSV file."""
        with open(self.filename, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["name", "dmrg"])
            writer.writeheader()
            writer.writerows(self.rows)

def generate_geometry(n, distance):
    geom = ""
    for i in range(n):
        geom += f"H 0.0 0.0 {i * distance}\n"
    return geom

def run_dmrg(n, dim, distance=1, save_data=True):
    geometry = generate_geometry(n, distance)
    mol = tq.Molecule(geometry=geometry, basis_set='sto-3g').use_native_orbitals()
    c, h1, g = mol.get_integrals(two_body_ordering="chem")

    driver = DMRGDriver(scratch="./tmp", symm_type=SymmetryTypes.SU2, n_threads=4)
    driver.initialize_system(n_sites=n, n_elec=n, spin=0)
    mpo = driver.get_qc_mpo(h1e=h1, g2e=g.elems, ecore=c, iprint=0)
    ket = driver.get_random_mps(tag="KET", bond_dim=dim, nroots=1)
    energies = driver.dmrg(mpo, ket, n_sweeps=30, bond_dims=[dim], noises=[1e-5] * 4 + [0], thrds=[1e-10] * 8, iprint=1, tol=1e-5)

    if save_data:
        table = ResultsTable(filename="results_dmrg_atomics.csv")
        table.add_entry(f"h{n}", energies)


if __name__ == "__main__":
    for n in range(2,32,2):
        run_dmrg(n,dim=200)
    # run_dmrg(16,dim=200)
    # for n in range(24,28,2):
    #     run_dmrg(n,dim=800)
    # for n in range(28,32,2):
    #     run_dmrg(n,dim=600)
    # run_dmrg(10, save_data=False)
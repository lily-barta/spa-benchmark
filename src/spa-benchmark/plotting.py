import pandas as pd
import matplotlib.pyplot as plt


def plot_runtime_scaling(filename="runtime_vs_n.csv"):
    data = pd.read_csv(filename)

    for col in data.columns:
        if col not in ["n", "total"]:
            plt.plot(data["n"], data[col], marker="o", label=col)

    plt.yscale("log")
    plt.xlabel("Number of hydrogens (n)")
    plt.ylabel("Time (s)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("runtime_scaling.pdf")
    plt.close()


def plot_accuracy_scaling(filename="results_vs_n.csv"):
    data = pd.read_csv(filename)
    data = data[data["n"] <= 14]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
    ax1.plot(data["n"], data["spa"]-data["fci"], "o-")
    ax2.plot(data["n"], data["fid"], "o-")
    
    ax1.set_ylabel("SPA error (Hartree)")
    ax2.set_xlabel("Number of hydrogen atoms)")
    ax2.set_ylabel("Fidelity")
    
    plt.tight_layout()
    plt.savefig("accuracy_scaling.pdf")
    plt.close()


def plot_spa_scaling(filename="results_vs_n.csv"):
    data = pd.read_csv(filename)
    plt.plot(data["n"], data["spa"], "o-")
    
    plt.ylabel("SPA energy (Hartree)")
    plt.xlabel("Number of hydrogen atoms)")
    
    plt.tight_layout()
    plt.savefig("spa_scaling.pdf")
    plt.close()


def plot_dissociation(n=6, filename=None):
    if filename is None:
      filename = f"results_h{n}.csv"
    data = pd.read_csv(filename)

    # Energy curve
    plt.plot(data["distance"], data["spa"], "o-", label="SPA")
    plt.plot(data["distance"], data["fci"], "o-", label="FCI")
    plt.xlabel("Interatomic distance (Å)")
    plt.ylabel("Energy (eH)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"dissociation_h{n}.pdf")
    plt.close()

    # Error curve
    plt.plot(data["distance"], data["spa"] - data["fci"], "o-")
    plt.xlabel("Interatomic distance (Å)")
    plt.ylabel("Error (eH)")
    plt.tight_layout()
    plt.savefig(f"error_h{n}.pdf")
    plt.close()

    # Fidelity curve
    plt.plot(data["distance"], data["fid"], "o-")
    plt.xlabel("Interatomic distance (Å)")
    plt.ylabel("Fidelity")
    plt.tight_layout()
    plt.savefig(f"fidelity_h{n}.pdf")
    plt.close()

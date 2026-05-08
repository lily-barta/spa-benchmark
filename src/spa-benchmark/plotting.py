import pandas as pd
import matplotlib.pyplot as plt


def plot_timing(filename="timing_vs_n.csv"):
    data = pd.read_csv(filename)

    for col in data.columns:
        if col not in ["n", "total_t"]:
            plt.plot(data["n"], data[col], marker="o", label=col)

    plt.yscale("log")
    plt.xlabel("Number of hydrogens (n)")
    plt.ylabel("Time (s)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("timing_vs_n.pdf")
    plt.close()


def plot_accuracy(filename="results_vs_n.csv"):
    data = pd.read_csv(filename)
    data = data[data["n"] <= 14]

    plt.plot(data["n"], 1 - data["fid"], "o-", label="1 - Fidelity")
    plt.plot(data["n"], data["spa"] - data["fci"], "o-", label="Error (eH)")

    plt.xlabel("Number of hydrogens (n)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("results_vs_n.pdf")
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

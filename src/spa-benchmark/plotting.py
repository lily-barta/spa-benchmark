import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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
    
    ax1.set_ylabel("SPA error (eH)")
    ax2.set_xlabel("Number of hydrogen atoms)")
    ax2.set_ylabel("Fidelity")
    
    plt.tight_layout()
    plt.savefig("accuracy_scaling.pdf")
    plt.close()


def plot_spa_scaling(filename="results_vs_n.csv"):
    data = pd.read_csv(filename)
    plt.plot(data["n"], data["spa"], "o-")
    
    plt.ylabel("SPA energy (eH)")
    plt.xlabel("Number of hydrogen atoms)")
    
    plt.tight_layout()
    plt.savefig("spa_scaling.pdf")
    plt.close()


def plot_dissociation_curves(n=6, filename=None):
    if filename is None:
      filename = f"results_h{n}.csv"
    data = pd.read_csv(filename)

    plt.plot(data["distance"], data["spa"], "o-", label="SPA", markersize=4)
    plt.plot(data["distance"], data["hf"], "v-", label="HF", markersize=4)
    plt.plot(data["distance"], data["ccsdt"], "d-", label="CCSD(T)", markersize=4)
    plt.plot(data["distance"], data["fci"], "s-", label="FCI", markersize=4)
    plt.xlabel("Interatomic distance (Å)")
    plt.ylabel("Energy (Hartree)")
    # plt.xlim(0.5,3.5)
    # plt.ylim(-3.3,-2.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"dissociation_energies_h{n}.pdf")
    plt.close()

def plot_dissociation_accuracy():
    files = ["results_h4.csv", "results_h6.csv", "results_h8.csv", "results_h10.csv"]
    dataframes = [pd.read_csv(f) for f in files]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:4]
    markers = ['o', 's', 'v', 'd']
    labels = ["$H_4$", "$H_6$", "$H_8$", "$H_{10}$"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
    for i, df in enumerate(dataframes):
        ax1.plot(df["distance"], df["spa"] - df["fci"], linestyle='-', marker=markers[i], color=colors[i], markersize=4)
        ax2.plot(df["distance"], df["fid"], linestyle='-', color=colors[i], marker=markers[i], markersize=4)
        ax2.plot(df["distance"], df["hf_fid"], linestyle='--', color=colors[i])
    
    # Axis formatting
    ax1.set_ylabel("SPA error (eH)")
    ax1.set_yticks([0.0, 0.02, 0.04, 0.06, 0.08])
    ax1.set_ylim(0.0, 0.08)
    ax1.set_xlim(0.5, 3.5)
    
    ax2.set_xlabel("Interatomic distance (Å)")
    ax2.set_ylabel("Fidelity")
    ax2.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax2.set_ylim(0.0, 1.0)
    ax2.set_xlim(0.5, 3.5)
    
    # Legend 1: line styles
    style_handles = [
        Line2D([0], [0], color='black', linestyle='-', label='SPA'),
        Line2D([0], [0], color='black', linestyle='--', label='HF')
    ]
    legend1 = ax1.legend(handles=style_handles, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    ax1.add_artist(legend1)
    
    # Legend 2: H_n series
    series_handles = [Line2D([0], [0], color=color, marker=marker, linestyle='-', markersize=5, label=label)
        for color, marker, label in zip(colors, markers, labels)
    ]
    ax1.legend(handles=series_handles, loc='upper right', bbox_to_anchor=(1.0, 0.85))
    
    plt.tight_layout()
    plt.savefig("dissociation_accuracy.pdf")
    plt.close()

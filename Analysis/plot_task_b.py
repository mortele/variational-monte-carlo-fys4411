import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess

def plot_simple_HO(filename="simple_HO", D=1, omega=1.0, alpha_range=(0.1,1.1, 11), save=False):
    alphas = np.linspace(*alpha_range)
    Ns = [1,10,100,500]
    stepLengths = [0.1, 0.1, 0.5, 1.0]

    filename = f"{filename}_{D}D.txt"
    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        k = 1
        for N, stepLength in zip(Ns, stepLengths):
            for i, alpha in enumerate(alphas):
                cpp_utils.vmcRun(D=D, filename=filename, N=N, alpha=alpha, stepLength=stepLength)
                print(f"Done {k}/{total}...")
                k += 1

    df = cpp_utils.vmcLoad(filename=filename)

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_N = df[ df.Particles == N ]
        E, E_var, alpha, MCCs = df_N.Energy.to_numpy(), df_N.Energy_var.to_numpy(), df_N["WF1"].to_numpy(), df_N["Metro-steps"].to_numpy()
        Error = np.sqrt(E_var/MCCs)
        ax.plot(alpha, E, c=c[i], label=f"{N =}", marker="o")
        ax.plot(alpha, E+Error, c=c[i], ls="--")
        ax.plot(alpha, E-Error, c=c[i], ls="--")
        print(f"Minimum energy at alpha = {alpha[np.argmin(E)]}")
        
    
    ax.legend(ncol=4, bbox_to_anchor=(1.05, 1.15))
    ax.set(xlabel=r"$\alpha$", ylabel="E $[\hbar \omega]$")
    ax.set_yscale("log")
    ylim = ax.get_ylim()
    ax.vlines(x=0.5, ymin=ylim[0], ymax=ylim[1], color="gray", ls="--", alpha=0.4)
    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()


if __name__ == "__main__":
    plot_simple_HO(D=1, save=True)
    plot_simple_HO(D=2, save=True)
    plot_simple_HO(D=3, save=True)
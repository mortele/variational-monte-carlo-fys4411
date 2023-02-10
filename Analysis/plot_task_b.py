import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess

def plot_simple_HO(filename="simple_HO", omega=1.0, alpha_range=(0.1,1.1, 11), save=False):
    alphas = np.linspace(*alpha_range)
    Ns = [1,10,100,500]
    stepLengths = [0.1, 0.1, 0.5, 1.0]
    D = 2

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
    for i, N in enumerate(Ns):
        fig, ax = plt.subplots()
        df_N = df[ df.Particles == N ]
        E, E_var, alpha, MCCs = df_N.Energy.to_numpy(), df_N.Energy_var.to_numpy(), df_N["WF1"].to_numpy(), df_N["Metro-steps"].to_numpy()
        Error = np.sqrt(E_var/MCCs)
        ax.plot(alpha, E, c=c[i], label=f"{N =}", marker="o")
        ax.plot(alpha, E+Error, c=c[i], ls="--")
        ax.plot(alpha, E-Error, c=c[i], ls="--")
        ax.legend()
        ax.set(xlabel=r"$\alpha$", ylabel="E $[\hbar \omega]$")
        print(f"Minimum energy at alpha = {alpha[np.argmin(E)]}")
        
        if save:
            plot_utils.save(filename.replace(".txt",f"{N}_plot"))
        plt.show()


if __name__ == "__main__":
    plot_simple_HO()
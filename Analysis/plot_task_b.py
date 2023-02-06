import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess

def plot_simple_HO(filename="simple_HO.txt", omega=1.0, alpha_range=(0.1,1.1, 11)):
    alphas = np.linspace(*alpha_range)
    Ns = [1,10,100,500]

    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        k = 1
        for N in Ns:
            for i, alpha in enumerate(alphas):
                cpp_utils.vmcRun(filename=filename, N=N, alpha=alpha)
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
        ax.plot(alpha, E+np.sqrt(E_var), c=c[i], ls="--", label=r"$\pm\sqrt{\sigma^2}$")
        ax.plot(alpha, E-np.sqrt(E_var), c=c[i], ls="--")
        ax.legend()

        y_lim_old = ax.get_ylim()
        y_lim_new = (0.75*np.min(E), y_lim_old[1])
        ax.set_ylim(y_lim_new)

        plot_utils.save(filename.replace(".txt",f"{N}_plot"))
        plt.show()

if __name__ == "__main__":
    plot_simple_HO()
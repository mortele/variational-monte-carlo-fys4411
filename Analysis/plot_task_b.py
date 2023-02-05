import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess

def plot_simple_HO(filename="simple_HO.txt", omega=1.0, alpha_range=(0.5,1.5, 11)):
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

    df = cpp_utils.vmcLoad(filename="simple_HO.dat")

    fig, ax = plt.subplots()
    c = plot_utils.colors
    for i, N in enumerate(Ns):
        df_N = df[ df.Particles == N ]
        E, E_std, alpha = df_N.Energy.to_numpy(), df_N.Energy_std.to_numpy(), df_N["WF1"].to_numpy()

        ax.plot(alpha, E, c=c[i], label=f"{N =}")
        ax.plot(alpha, E+E_std, c=c[i], ls="--")
        ax.plot(alpha, E-E_std, c=c[i])

    ax.set_yscale("log")
    ax.legned()
    plot_utils.save(filename.replace(".txt","_plot"))
    plt.show()

if __name__ == "__main__":
    plot_simple_HO()
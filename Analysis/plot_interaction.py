import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils

def plot_energy_per_particle(filename="energy_per_particle", D=3, save=False):
    Ns = [1, 10, 50, 100]
    alphas = np.linspace(0.3, 0.7, 11)

    deltaT = 0.5

    filename = f"{filename}_{D}D.txt"
    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        for i, N in enumerate(Ns):
            for j, alpha in enumerate(alphas):
                cpp_utils.vmcRun(D=D, filename=filename, N=N, alpha=alpha, stepLength=deltaT, importance=True, logMet=15, logEq=10)
                print(f"Done {i*len(alpha)+j+1}/{total}...")

    df = cpp_utils.vmcLoad(filename=filename)

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_N = df[ df.Particles == N ]
        E, E_var, alpha, MCCs = df_N.Energy.to_numpy(), df_N.Energy_var.to_numpy(), df_N["WF1"].to_numpy(), df_N["Metro-steps"].to_numpy()
        Error = np.sqrt(E_var/MCCs)
        ax.errorbar(alpha, E/N, Error/N, c=c[i], label=f"{N =}", marker="o")
        print(f"Minimum energy at alpha = {alpha[np.argmin(E)]}")

    ax.legend(ncol=4, bbox_to_anchor=(1.05, 1.15))
    ax.set(xlabel=r"$\alpha$", ylabel="E $[\hbar \omega]$")
    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()

if __name__ == "__main__":
    plot_energy_per_particle()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils
import seaborn as sns
import time

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
cmap = plot_utils.cmap 


def plot_alpha_search(filename="gradientSearch", D=3, omega=1.0, alpha_range=(0.2,1.2, 11), save=False):
    alphas = np.linspace(*alpha_range)
    Ns = [10] # Number of particles
    stepLengths = [1.5]#, 0.1, 0.5, 1.0]
    epsilon = 0.01
    logMet = 16 # 2 ^ 16 = 65536
    logEq = 14 # 2 ^ 14 = 16384

    filename = f"{filename}_{D}D.txt" # prolly a good idea to add the dimension
    if not cpp_utils.dataPath(filename).exists(): # If the file does not exist, run the gradient search code
        total = len(Ns)*len(alphas)
        for N, stepLength in zip(Ns, stepLengths):
            for i, alpha in enumerate(alphas):
                cpp_utils.gradientRun(logMet=logMet, logEq=logEq, D=D, N=N, alpha=alpha, stepLength=stepLength, epsilon=epsilon, filename=filename)
                print(f"Gradient run done for alpha = {alpha} and N = {N}...")

    info = f"_D={D}_N={Ns}_stepLength={stepLengths[0]}_met={logMet}_eq={logEq}_eps={epsilon}"

    df = cpp_utils.gradientLoad(filename=f"../Data/{filename}") # Load the data
    df_detailed = cpp_utils.gradientLoad(filename=f"../Data/detailed_{filename}")

    ax = sns.lineplot(data=df_detailed, x="Epoch", y="Alpha", hue="Alpha_0", legend="full", palette=cmap)
    ax.get_legend().remove()
    plt.xlabel("Epoch")
    plt.ylabel("Alpha")

    if save:
        plot_utils.save("alpha_search" + info)
    plt.show()

    return df, df_detailed, info


def plot_energy_var(df_detailed, info, save=False):
    # Plot the energy variance
    ax = sns.lineplot(data=df_detailed, x="Alpha", y="Energy_var", hue="Alpha_0", legend=False, palette=cmap)

    norm = plt.Normalize(df_detailed["Alpha_0"].min(), df_detailed["Alpha_0"].max())

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label="$Alpha_0$", cmap=cmap)

    plt.xlabel("Alpha")
    plt.ylabel("$Var(E)[(\hbar \omega)^2]$")


    if save:
        plot_utils.save("alpha_search_energy_var" + info)
    plt.show()


def plot_energy_per_particle(filename="GD_energy_per_particle", D=3, interacting=False, save=False):
    Ns = np.array([1, 10, 30, 50])
    alphas = [0.51]
    epsilon = 0.01 # this does not need to be super small. This is a tolerance with respect to the gradient, but notice 
                     # that this gets multiplied by the learning rate, so the paramaeter update at the end is what matters
    lr = 0.01 # this is scaled by the number of particles as epsilon = 1/(ln(N) + 1)
    stepLength = 1.5
    logMet = 2^19 # 2^19 = 524288
    logEq = 2^14 # 2 ^ 14 = 16384

    # start time measurement
    start = time.time()
    filename = f"{filename}_{D}D_interac{interacting}.txt"
    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        for i, N in enumerate(Ns):
            for j, alpha in enumerate(alphas):
                cpp_utils.gradientRun(logMet=logMet, logEq=logEq, D=D, epsilon=epsilon, filename=filename, N=N, alpha=alpha, stepLength=stepLength, interacting=interacting, lr = lr)
                print(f"Done N = {N}, alpha = {alpha} {i*len(alphas)+j+1}/{total}...")
                print(f"Time elapsed: {time.time() - start}")
                print("Expected time if linear: ", (time.time() - start)/(i*len(alphas)+j+1)*total / 60, "min")

    df_detailed = cpp_utils.gradientLoad(filename=f"../Data/detailed_{filename}")

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_detailed_N = df_detailed[ df_detailed.Particles == N ]
        E, E_std, alpha, MCCs = df_detailed_N["Energy"].to_numpy(), df_detailed_N["Energy_var"].to_numpy(), df_detailed_N["WF1"].to_numpy(), df_detailed_N["Metro-steps"].to_numpy()
        ax.errorbar(alpha, E/N, np.sqrt(E_std)/(MCCs), c=c[i], label=f"{N =}")
        #print(f"Minimum energy at alpha = {alphas[np.argmin(E)]}")

    ax.legend(ncol=4, bbox_to_anchor=(1.05, 1.15))
    ax.set(xlabel=r"$\alpha$", ylabel="$expectation E_L [\hbar \omega]$")

    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()

if __name__ == "__main__":
    df, df_detailed, info = plot_alpha_search(filename="GD",D=3, save=True)

    plot_energy_var(df_detailed, info, save=True)

    plot_energy_per_particle(filename="GD_energy_per_particle", D=3, interacting=True, save=False)



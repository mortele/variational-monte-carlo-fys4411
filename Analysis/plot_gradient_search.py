import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils
import seaborn as sns

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
cmap = plot_utils.cmap 


def plot_alpha_serch(filename="gradientSearch", D=3, omega=1.0, alpha_range=(0.1,1.1, 11), save=False):
    alphas = np.linspace(*alpha_range)
    Ns = [10] # Number of particles
    stepLengths = [0.1]#, 0.1, 0.5, 1.0]
    epsilon = 0.0001
    logMet = 4
    logEq = 3

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


### the following function might not be usefull for plotting, but it helped me catch some potential problems, so i will leave it now so i can test things later.
"""
def plot_energy_search_comparison(filename="gradientSearch_comparison", D=3, omega=1.0, alpha_range=(0.1,1.1, 11), save=False):
    #select alpha randomly from the range
    Ns = [1,10,100,500]

    # make alphas be a list of length len(Ns) with random alphas between alpha_range
    alphas = [np.random.choice(np.linspace(*alpha_range)) for _ in range(len(Ns))]
    stepLengths = [0.05]#, 0.1, 0.5, 1.0]
    stepLength = stepLengths[0]
    epsilon = 0.0001
    logMet = 6
    logEq = 5

    filename = f"{filename}_{D}D.txt" # prolly a good idea to add the dimension
    if not cpp_utils.dataPath(filename).exists(): # If the file does not exist, run the gradient search code
        for N in Ns:
            # select alpha according to N index
            alpha = alphas[Ns.index(N)]
            cpp_utils.gradientRun(logMet=logMet, logEq=logEq, lr=0.05, D=D, N=N, alpha=alpha, stepLength=stepLength, epsilon=epsilon, filename=filename)
            print(f"Gradient run done for alpha = {alpha} and N = {N}...")

    info = f"_D={D}_N={Ns}_stepLength={stepLengths[0]}_met={logMet}_eq={logEq}_eps={epsilon}"

    df_detailed = cpp_utils.gradientLoad(filename=f"../Data/detailed_{filename}")

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_N = df_detailed[ df_detailed.Particles == N ]
        E, E_var, alpha, MCCs = df_N.Energy.to_numpy(), df_N.Energy_var.to_numpy(), df_N["WF1"].to_numpy(), df_N["Metro-steps"].to_numpy()
        print("E",E)
        exit()
        Error = np.sqrt(E_var/MCCs)
        ax.plot(alpha, E, c=c[i], label=f"{N =}", marker="o")
        ax.plot(alpha, E+Error, c=c[i], ls="--")
        ax.plot(alpha, E-Error, c=c[i], ls="--")
        print(f"Minimum energy at alpha = {alpha[np.argmin(E)]}")

    if save:
        plot_utils.save("minimize_energy_search" + info)
    plt.show()

    return df, df_detailed, info    
"""

if __name__ == "__main__":
    df, df_detailed, info = plot_alpha_serch(D=3, save=True)

    plot_energy_var(df_detailed, info, save=True)




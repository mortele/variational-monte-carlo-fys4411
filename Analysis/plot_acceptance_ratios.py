import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils

def plot_acceptance(filename="acceptance", stepLength_range=(0.1, 10, 10), save=False, impo=False):
    N = 3
    n = 100
    stepLengths = np.linspace(*stepLength_range)
    alphas = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8]
    xlabel = {
        False: r"$\delta$",
        True: r"$\Delta t$",
    }
    extension = {
        False: "met",
        True: "impo"
    }

    filename = f"{filename}_{extension[impo]}.txt"
    total = len(stepLengths)*len(alphas)
    if not cpp_utils.dataPath(filename).exists():
        for i, stepLength in enumerate(stepLengths):
            for j, alpha in enumerate(alphas):
                cpp_utils.vmcRun(logMet=15, logEq=10, D=3, filename=filename, stepLength=stepLength, alpha=alpha, importance=impo)
                print(f"{i*len(alphas)+j+1}/{total}")

    df = cpp_utils.vmcLoad(filename)

    fig, ax = plt.subplots()
    for alpha in alphas:
        A = df[df["WF1"]==alpha]["Accept_ratio"]
        ax.plot(stepLengths, A, label=rf"$\alpha$ = {alpha}")

    ax.set(xlabel=xlabel[impo], ylabel="$A_r$")
    ax.legend()
    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()

if __name__ == '__main__':
    # plot_acceptance(stepLength_range=(0.1, 10, 100), impo=False, save=True)
    plot_acceptance(stepLength_range=(0.1, 10, 100), impo=True, save=True)
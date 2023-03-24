import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils

def number_of_mccs(filename="increase_mccs", D=3, omega=1.0, mccs_range=(5,22, 50), save=False):
    mccs = np.linspace(*mccs_range)

    N = 3
    stepLength = 1.5
    deltaT = 0.5
    filename = f"{filename}_{D}.txt"
    n = len(mccs)
    if not cpp_utils.dataPath(filename).exists():
        for i, mcc in enumerate(mccs):
            cpp_utils.vmcRun(logMet=mcc, logEq=15, D=D, filename=filename, N=N, alpha=0.51, stepLength=deltaT, importance=True)
            cpp_utils.vmcRun(logMet=mcc, logEq=15, D=D, filename=filename, N=N, alpha=0.51, stepLength=stepLength, importance=False)
            print(f"{i+1}/{n}")
    

    df = cpp_utils.vmcLoad(filename=filename)

    df["Metro-steps"] = np.log2(df["Metro-steps"])
    df_impo = df[df["Imposampling"] == 1]
    df_brute = df[df["Imposampling"] == 0]

    fig, ax = plt.subplots()
    ax.plot(df_impo["Metro-steps"], df_impo["Energy_var"], label="Metropolis-Hastings")
    ax.plot(df_brute["Metro-steps"], df_brute["Energy_var"], label="Metropolis")
    ax.set(xlabel=r"$log_{2}(M)$", ylabel=r"$Var(E_L) [\hbar^2 \omega^2]$")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    number_of_mccs(filename="increase_mccs")
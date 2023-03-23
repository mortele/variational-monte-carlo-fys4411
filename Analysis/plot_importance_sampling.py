import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils

def number_of_mccs(filename="increase_mccs", D=3, omega=1.0, mccs_range=(0,6, 10), save=False):
    mccs = np.linspace(0,6.5,20)

    N = 3
    stepLength = 0.4
    deltaT = 0.01
    filename_noeq = f"{filename}_{D}_noeq.txt"
    filename_eq = f"{filename}_{D}_eq.txt"
    n = 2*len(mccs)
    if not cpp_utils.dataPath(filename_noeq).exists():
        for i, mcc in enumerate(mccs):
            cpp_utils.vmcRun(logMet=mcc, logEq=0, D=D, filename=filename_noeq, N=N, alpha=0.6, stepLength=deltaT, importance=True)
            cpp_utils.vmcRun(logMet=mcc, logEq=0, D=D, filename=filename_noeq, N=N, alpha=0.6, stepLength=stepLength, importance=False)
            print(f"{i+1}/{n}")
    if not cpp_utils.dataPath(filename_eq).exists():
        for i, mcc in enumerate(mccs):
            cpp_utils.vmcRun(logMet=mcc, logEq=4, D=D, filename=filename_eq, N=N, alpha=0.6, stepLength=deltaT, importance=True)
            cpp_utils.vmcRun(logMet=mcc, logEq=4, D=D, filename=filename_eq, N=N, alpha=0.6, stepLength=stepLength, importance=False)
            print(f"{int(n/2)+i+1}/{n}")
    


    df = cpp_utils.vmcLoad(filename=filename_noeq)
    fig, ax = plt.subplots()

    df["Metro-steps"] = np.log10(df["Metro-steps"])
    df_impo = df[df["Imposampling"] == 1]
    df_brute = df[df["Imposampling"] == 0]

    c = plot_utils.colors
    ax.plot(df_impo["Metro-steps"], df_impo["Energy"], ls="--", c=c[0])
    ax.plot(df_brute["Metro-steps"], df_brute["Energy"], ls="--", c=c[1])
   
    df = cpp_utils.vmcLoad(filename=filename_eq)

    df["Metro-steps"] = np.log10(df["Metro-steps"])
    df_impo = df[df["Imposampling"] == 1]
    df_brute = df[df["Imposampling"] == 0]

    ax.plot(df_impo["Metro-steps"], df_impo["Energy"], label="Metropolis-Hastings", c=c[0])
    ax.plot(df_brute["Metro-steps"], df_brute["Energy"], label="Metropolis", c=c[1])
    ax.hlines(y=0.5*N*D, xmin=mccs_range[0], xmax=mccs_range[1], ls="-.", color="gray", alpha=.5)
    ax.set(xlabel=r"$log_{10}(M)$", ylabel=r"$\langle E_L \rangle [\hbar \omega]$",ylim=(0.5*N*D-0.5, 5.5))
    ax.legend()
    plt.show()
if __name__ == "__main__":
    number_of_mccs()
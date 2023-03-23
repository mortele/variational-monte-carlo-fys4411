import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

import plot_utils
import cpp_utils


def plot_timing(filename="", D=1, Ns = np.arange(2, 10), trials = 2, save=False):
    filename = f"{filename}_{D}D.txt"
    MCC0 = 2**13
    if not cpp_utils.dataPath(filename).exists():
        for i, N in enumerate(Ns):
            MCC = np.log2(MCC0*N)
            for _ in range(trials):
                cpp_utils.timingRun(D=D, filename=filename, N=N, alpha=0.43, stepLength=0.2, analytical=1, logMet=MCC, logEq=4)
                cpp_utils.timingRun(D=D, filename=filename, N=N, alpha=0.43, stepLength=0.2, analytical=0, logMet=MCC)
            print(f"Done {i+1}/{len(Ns)}... {N = }")

    df = cpp_utils.vmcLoad(filename=filename)
    df["Time"] = df["Time"]*1e-6
    

    # Extract time in seconds
    df_analytical = df[df["Analytical"] == 1]
    df_numerical = df[df["Analytical"] == 0]

    group_analytical = df_analytical.groupby("Particles", as_index=False)
    group_numerical = df_numerical.groupby("Particles", as_index=False)

    time_analytical = group_analytical["Time"].mean().to_numpy()[:,1]
    std_analytical = group_analytical["Time"].std().to_numpy()[:,1]
    time_numerical = group_numerical["Time"].mean().to_numpy()[:,1]
    std_numerical = group_numerical["Time"].std().to_numpy()[:,1]
    
    cols = ["Energy", "Energy_std"]
    E_analytical = group_analytical[cols].mean().to_numpy()[:,1:]
    E_numerical = group_numerical[cols].mean().to_numpy()[:,1:]
    energy_analytical, energy_std_analytical = E_analytical[:,0], E_analytical[:,1]
    energy_numerical, energy_std_numerical = E_numerical[:,0], E_numerical[:,1]

    fig, ax = plt.subplots()
    ax.errorbar(Ns, time_analytical, std_analytical, label="Analytical")
    ax.errorbar(Ns, time_numerical, std_numerical, label="Numerical")
    ax.legend()
    ax.set(xlabel="N", ylabel="Time [s]", title=f"Derivative time comparison D = {D}")
    print("Energy per particle ana vs num")
    print(energy_analytical/Ns)
    print(energy_numerical/Ns)
    print("\n\n\n")
    print("Std per particle ana vs num")
    print(energy_std_analytical/Ns)
    print(energy_std_numerical/Ns)

    if save: 
        plot_utils.save(filename.replace(".txt", ""))
    plt.show()

if __name__ == "__main__":
    Ns = np.array([1,2,3,4,5,10,15,20,25,30,40,50])
    plot_timing("derivativeTiming", Ns=Ns, trials=10, D=3, save=True)

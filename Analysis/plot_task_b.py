import numpy as np
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess

def plot_simple_HO(omega=1.0, alpha_range=(0.5,1.5, 11)):
    alphas = np.linspace(*alpha_range)
    Ns = [1,10,100,500]

    for N in Ns:
        
        E = np.zeros_like(alphas)
        E_std = np.zeros_like(alphas)
        for i, alpha in enumerate(alphas):
            subprocess.run()
            # Read in E and E_std 
            pass

if __name__ == "__main__":
    # plot_simple_HO()
    cpp_utils.runVMC()
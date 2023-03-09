import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess
import plot_task_b

N_max=100
steps=5
n=np.linspace(1,N_max,steps)
TimeAnalytical=np.zeros(steps)
for i in range(1,steps):
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], timing=True, filename="timeana.txt")
    

print(n)
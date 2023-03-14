import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils

import subprocess
import plot_task_b

N_max=10 #probably more needed
steps=3
n=np.linspace(2,N_max,steps)
print(n)
TimeAnalytical=np.zeros(steps) #make sure to vary the amount of  metropolis steps with the amount of particles
for i in range(0,steps):
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=True, timing=True, filename="1timeana.txt")
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=True, timing=True, filename="2timeana.txt")
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=True, timing=True, filename="3timeana.txt")
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=False, timing=True, filename="1timenum.txt")
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=False, timing=True, filename="2timenum.txt")
    TimeAnalytical[i]=cpp_utils.vmcRun(N=n[i], analytical=False, timing=True, filename="3timenum.txt")

print(n)



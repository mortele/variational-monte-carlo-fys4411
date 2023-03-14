import numpy as np 
import pandas as pd 
import plot_utils
import matplotlib.pyplot as plt 
import cpp_utils as cput

A1=cput.vmcLoad("../Data/1timeana.txt")

A2=cput.vmcLoad("../Data/2timeana.txt")
A3=cput.vmcLoad("../Data/3timeana.txt")
N1=cput.vmcLoad("../Data/1timenum.txt")
N2=cput.vmcLoad("../Data/2timenum.txt")
N3=cput.vmcLoad("../Data/3timenum.txt")

"""
A1=pd.read_csv("../Data/1timeana.txt")
A2=pd.read_csv("../Data/2timeana.txt")
A3=pd.read_csv("../Data/3timeana.txt")
N1=pd.read_csv("../Data/1timenum.txt")
N2=pd.read_csv("../Data/2timenum.txt")
N3=pd.read_csv("../Data/3timenum.txt")
"""
print(A1.columns)
print(A1)
print(A1.loc[1].at["Time"])
Aav=np.zeros(len(A1["Time"]))
Nav=np.zeros(len(N1["Time"]))
PAav=np.zeros(len(A1["Particles"]))
PNav=np.zeros(len(N1["Particles"]))
for i in range(0, len(A1["Time"])):
    Aav[i]=(A1.loc[i].at["Time"]+A2.loc[i].at["Time"]+A3.loc[i].at["Time"])/3
    PAav[i]=(A1.loc[i].at["Particles"])
    Nav[i]=(N1.loc[i].at["Time"]+N2.loc[i].at["Time"]+N3.loc[i].at["Time"])/3
    PNav[i]=(N1.loc[i].at["Particles"])
print(Aav)
print(PAav)
plt.plot(PAav, Aav, label="Analytical")
plt.plot(PNav, Nav, label="Numerical")
plt.legend()
plt.savefig("testplottiming.pdf")
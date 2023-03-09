import numpy as np 
import pandas as pd 
import plot_utils
import matplotlib.pyplot as plt 


df = pd.read_csv("../Data/time.txt")

df = df[['Time', 'Particles']]

plt.plot(df)
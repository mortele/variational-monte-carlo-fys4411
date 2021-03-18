#!/usr/bin/env python3

import os 

import numpy as np
from numpy.random import randint
from time import time

#where to save data
DATA_ID = "build/results"

def data_path(data_id): 
    return os.path.join(DATA_ID, data_id)

def tsboot(data, statistic, R, l):
    t = np.zeros(R)
    n = len(data)
    k = int(np.ceil(float(n)/l))
    inds = np.arange(n)
    t0 = time()

    #time series bootstrap
    for i in range(R):
        #bootstrap sample from k cunks of data
        #chunksize l
        _data = np.concatenate([data[j:j+l] for j in randint(0, n-1, k)])[0:n]
        t[i] = statistic(_data)

    #analysis
    print(  "Runtime: %g sec" % (time()-t0))
    print(  "Bootstrap statistics: ")
    print(  "original   bias    std.error")
    print(  "%7g    %4g %14g" % (statistic(data), \
                np.mean(t) - statistic(data), \
                np.std(t)   ))
    return t

def stat(data):
    return np.mean(data)


infile = open(data_path("energies.csv"), "r")

#read in data
X = np.genfromtxt(infile, delimiter=";")

t = tsboot(X, stat, 2**12, 2**10)

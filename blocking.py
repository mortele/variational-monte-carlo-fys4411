#!/usr/bin/env python3
"""
first attempt at the blocking code from Marius Jonsson's paper (Standard error
estimation by an automated blocking
method)[https://www.duo.uio.no/bitstream/handle/10852/68360/PhysRevE.98.043304.pdf?sequence=2&isAllowed=y]
"""

from numpy import log2, floor, array, zeros, mean, arange, cumsum, var, loadtxt
import numpy.linalg as lin
import os
import csv
import re


DATA_ID = "results"

def data_path(data_id): 
    return os.path.join(DATA_ID, data_id)

def block(X):
    #preliminary
    d = log2(len(X))
    if ( d - floor(d) != 0 ):
        print( "Warning: Data size: %g, not a power of 2." % floor(2**d) )
        print( "Truncating data to %d" % 2**floor(d) )
        X = X[:2**int(floor(d))]
    d = int(floor(d))
    #n = 2**d #This seems to entirely be overwritten later without use? 
    s, gamma = zeros(d), zeros(d)
    mu = mean(X)

    #Estimate autocovariance and variances
    #per blocking trasnform
    for i in arange(0, d):
        n = len(X)
        #estimate autocovariance of X
        gamma[i] = (n)**(-1) * sum( ( X[0:(n-1)] - mu )*( X[1:n] - mu ) )
        #Variance of X
        s[i] = var(X)
        #Blocking transform
        X = .5 * ( X[0::2] + X[1::2] )

    #Here ends what I mostly get. 

    #Generate test observator M_k from theorem 
    M = ( cumsum(\
        ( (gamma/s)**2 * 2**arange(1, d+1)[::-1] )[::-1] ) )[::-1]
    #List of magic numbers
    q =array([6.634897,  9.210340,  11.344867, 13.276704, 15.086272, 
              16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
              24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 
              31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
              38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 
              45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    #Determine blocking stop with magic
    for k in arange(0, d):
        if (M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    return s[k]/2**(d-k)

#X = loadtxt("../block/resources/data.txt")
#print("standard error = %g" %block(X)**.5)
"""
runexample: 
    $ standard error = 0.00392796
"""
#read in data
#1c
print("1c) ")
infile = open(data_path("energies_1c.csv"), "r")
csv_reader = csv.reader(infile, delimiter=";")
i = 0
for row in csv_reader:
    print("*"*50)
    row = list(filter(None, row))
    X = array(row, dtype=float)
    print("standard error = %g" %block(X)**.5)
infile.close()

#1d
print("1d) ")
infile = open(data_path("energies_1d.csv"), "r")
csv_reader = csv.reader(infile, delimiter=";")
for row in csv_reader:
    print("*"*50)
    row = list(filter(None, row))
    X = array(row, dtype=float)
    print("standard error = %g" %block(X)**.5)
    #break
infile.close()

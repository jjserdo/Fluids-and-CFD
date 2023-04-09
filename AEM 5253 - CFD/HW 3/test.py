# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:05:46 2022

@author: jjser
"""
import numpy as np

t = np.arange(2)
x = np.arange(5)
D = 1
dx = 1

#C = np.zeros((len(t),len(x)))
#C[0] = cti
for n in range(len(t)-1):
    A = np.zeros((len(x),len(x)))
    B = np.zeros((len(x),1))
    print(A)
    print(B)
    #d = D/(2*dx**2)
    d = 1
    B[0] = 0
    B[-1] = 0
    A[0][0] = 1
    A[-1][-1] = 1
    for i in np.arange(1,len(B)-1):
        A[i,i-1] = -d
        A[i,i] = 1+2*d
        A[i,i+1] = -d
        B[i-1] = d
        B[i] = 1-2*d
        B[i+1] = d
    print(A)
    print(B)
    #C[n+1] = np.linalg.solve(A,B)
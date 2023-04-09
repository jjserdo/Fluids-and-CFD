# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 19:18:21 2023

@author: jjser
Justine John A. Serdoncillo 
AEM 8202 Fluids II
Homework 3
"""

import numpy as np
import matplotlib.pyplot as plt

# %%
"""
    Personal Created Functions
    - FalknerSkan()
"""

def Blasius(deta,nf,fpp0):
    N = int(nf/deta)
    etas = np.linspace(0,nf,N,True)
    f = np.zeros(len(etas))
    fp = np.zeros(len(etas))
    fpp = np.zeros(len(etas))
    fpp[0] = fpp0
    for space in range(len(etas)-1):
        etas[space+1] = etas[space] + deta
        f[space+1] = f[space] + deta * fp[space]
        fp[space+1] = fp[space] + deta * fpp[space]
        fpp[space+1] = fpp[space] - deta/2 * f[space] * fpp[space]
    return etas, f, fp, fpp

def FalknerSkan(deta,nf,fpp0,alpha):
    N = int(nf/deta)
    etas = np.linspace(0,nf,N,True)
    f = np.zeros(len(etas))
    fp = np.zeros(len(etas))
    fpp = np.zeros(len(etas))
    fpp[0] = fpp0
    for space in range(len(etas)-1):
        etas[space+1] = etas[space] + deta
        f[space+1] = f[space] + deta * fp[space]
        fp[space+1] = fp[space] + deta * fpp[space]
        sheesh = (f[space+1] + f[space])/2
        sheesh2 = (fp[space+1] + fp[space])/2
        slope = (alpha+1)/2 * sheesh * fpp[space] + alpha * sheesh2**2 - alpha
        fpp[space+1] = fpp[space] - deta * slope
    return etas, f, fp, fpp

def shootFalknerSkan(deta,nf,alpha):
    diff = 10
    TOL = 1E-8
    final = 0.01
    MAXITE = 10000
    fpp0 = 0
    ite = 0
    bestfpp0 = 0
    bestdiff = 10
    
    
    while diff > TOL and ite < MAXITE:
        ite += 1
        etas, f, fp, fpp = FalknerSkan(deta,nf,fpp0,alpha)
        diff = np.abs(fp[-1] - 1)
        if diff < bestdiff:
            bestdiff = diff
            bestfpp0 = fpp0
        fpp0 += final/MAXITE
    return bestfpp0, ite, bestdiff
    
# %% 
"""
    Problem Specific Functions
"""

def testBlasius():
    etas, f, fp, fpp = Blasius(0.02,7,0.332)
    fig, ax = plt.subplots()
    ax.plot(etas, fp)

    return etas, f, fp, fpp

def testFalknerSkan(case):
    if case == 'A':
        etas, f, fp, fpp = FalknerSkan(0.02,5,0.332,0)
        fig, ax = plt.subplots()
        ax.plot(etas, fp)
        return etas, f, fp, fpp
    elif case == 'B':
        alpha = 1/9
        fpp0, ite, diff = shootFalknerSkan(0.01, 100, alpha)
        print(fpp0)
        print(ite)
        print(diff)
        etas, f, fp, fpp = FalknerSkan(0.01,100, fpp0, alpha)
        fig, ax = plt.subplots()
        ax.plot(etas, fp)
        ax.grid()
        ax.plot(np.linspace(0,100,100,True), np.ones(100), '--')
        
def problem2():
    """
        Falkner-Skan Similarity Solutions of the Laminar Boundary-Layer Equations
        - fppp + (n+1)/2 * f * fpp - n * fp**2 + n = 0
            - Given: f(0) = 0, fp(0) = 0
            - Shoot with: fpp(0) = idk
            - Target: fp(inf) = 1
    """
    ALPHAS = [-0.0904, 0, 1/3, 1/9, 1]
    
    fig, ax = plt.subplots()
    ax.set_title('Falkner-Skan Solutions')
    ax.set_xlabel('$ \\eta $')
    ax.set_ylabel('$ fp(\\eta) $')
    ax.set_xlim([0,10])
    ax.set_ylim([0,2])
    ax.plot(np.linspace(0,100,100,True), np.ones(100), '--')
    for alpha in ALPHAS:
        fpp0, ite, bestdiff = shootFalknerSkan(0.01, 10, alpha)
        if alpha == ALPHAS[0]:
            fpp0 = 0
        print('$f^{pp}(0) $= '+str(np.around(fpp0,4)))
        etas, f, fp, fpp = FalknerSkan(0.01, 10, fpp0, alpha)
        label = "$\\alpha$ = " + str(np.around(alpha,4))
        ax.plot(etas, fp, label=label)
    ax.grid()
    ax.legend(loc='lower right')    
    
# %% 
"""
    Main Function
"""

if __name__ == "__main__":
    #etas, f, fp, fpp = testFalknerSkan('A')
    testFalknerSkan('B') # works for 0
    #etas, f, fp, fpp = testBlasius()
    #problem2()
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 19:18:21 2023

@author: jjser
Justine John A. Serdoncillo 
AEM 8202 Fluids II
Homework 3 v2 Updates:
    - shooting is just finding the best fpp0 by screening
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
    gg = True
    for space in range(len(etas)-1):
        etas[space+1] = etas[space] + deta
        f[space+1] = f[space] + deta * fp[space]
        fp[space+1] = fp[space] + deta * fpp[space]
        sheesh = (f[space+1] + f[space])/2
        slope = (alpha+1)/2 * sheesh * fpp[space] + alpha * fp[space+1]**2 - alpha
        fpp[space+1] = fpp[space] - deta * slope
        if fp[space+1] > 100:
            gg = False
            break
    return etas, f, fp, fpp, gg

def shootFalknerSkan(deta,nf,alpha):
    diff = 10
    TOL = 1E-4
    MAXITE = 1000
    fpp0 = 0
    ite = 0
    bestfpp0 = 0
    bestdiff = 10
    while diff > TOL and ite < MAXITE:
        ite += 1
        etas, f, fp, fpp, gg = FalknerSkan(deta,nf,fpp0,alpha)
        if gg == True:
            #print(fp[-1])
            diff = np.abs(fp[-1] - 1)
            if diff < bestdiff:
                bestdiff = diff
                bestfpp0 = fpp0
            #print('ite = '+ str(ite))
            #print('fpp0 = '+ str(fpp0))
            #print('diff = '+ str(diff))
        else:
            diff = 10
            fpp0 += 0.001
        
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
        alpha = -0.0904
        fpp0, ite, diff = shootFalknerSkan(0.02, 10, alpha)
        etas, f, fp, fpp = FalknerSkan(0.020, 10, fpp0, alpha)
        fig, ax = plt.subplots()
        ax.plot(etas, fp)
        
def problem2():
    """
        Falkner-Skan Similarity Solutions of the Laminar Boundary-Layer Equations
        - fppp + (n+1)/2 * f * fpp - n * fp**2 + n = 0
            - Given: f(0) = 0, fp(0) = 0
            - Shoot with: fpp(0) = idk
            - Target: fp(inf) = 1
    """
    ALPHAS = [-0.0904, 0, 1/3, 1/9, 1]
    ITES = np.zeros(len(ALPHAS))
    DIFFS = np.zeros(len(ALPHAS))
    FPP0S = np.zeros(len(ALPHAS))
    
    fig, ax = plt.subplots()
    ax.set_title('Falkner-Skan Solutions')
    ax.set_xlabel('$ \\eta $')
    ax.set_ylabel('$ fp(\\eta) $')
    ax.set_xlim([0,10])
    ax.set_ylim([0,2])
    ax.plot(np.linspace(0,10,100,True), np.ones(100), '--')
    for a in range(len(ALPHAS)):
        fpp0, ite, diff = shootFalknerSkan(0.1, 10, ALPHAS[a])
        FPP0S[a] = fpp0
        ITES[a] = ite
        DIFFS[a] = diff
        etas, f, fp, fpp, gg = FalknerSkan(0.1, 10, fpp0, ALPHAS[a])
        label = "$\\alpha$ = " + str(np.around(ALPHAS[a],4))
        ax.plot(etas, fp, label=label)
    ax.grid()
    ax.legend(loc='lower right')
    
    fig1, ax1 = plt.subplots()
    ax1.plot(ALPHAS, ITES,'*')
    ax1.set_title('Iteration counts')
    
    fig2, ax2 = plt.subplots()
    ax2.plot(ALPHAS, DIFFS,'*')
    ax2.set_title('Differences counts')
    return ALPHAS, ITES, DIFFS, FPP0S
    
# %% 
"""
    Main Function
"""

a,b,c,d = problem2()
print(np.around(a,4))
print(b)
print(np.around(c,4))
print(np.around(d,4))

#testFalknerSkan('B')
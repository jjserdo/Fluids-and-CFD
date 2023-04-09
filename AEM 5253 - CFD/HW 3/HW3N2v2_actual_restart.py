# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:33:14 2022

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# Initializing Everything
ni = 21 # 20
nj = 21 # 20
x = np.zeros((ni,nj))
y = np.zeros((ni,nj))

dx = 2/(ni-1)
dy = 1/(nj-1)

def bound(arr1, arr2):
    ni,nj = arr1.shape
    # Create top and bottom to be dy/dx = 0
    for i in range(ni):
        arr1[i,0] = arr2[i,1]
    for i in range(ni):
        arr1[i,-1] = arr2[i,-2]
    # Initialize boundary condition f(0) = 0
    arr1[0,:] = 0
    # Initialize boundary condition f(y) = y
    for j in range(nj):
        arr1[-1,j] = y[-1,j]
    return arr1

def grid():
    # Initializing the Grid x and y
    for i in np.arange(0,ni): # goes from 1 to ni
        for j in np.arange(0,nj): # goes from 1 to nj
            x[i,j] = (i) * 2.0 / (ni-1)
            y[i,j] = (j) * 1.0 / (nj-1)

def anal():
    u = np.zeros((ni,nj))
    for i in np.arange(0,ni): # goes from 1 to ni
        for j in np.arange(0,nj): # goes from 1 to nj
            u[i,j] = x[i,j] / 4
            for n in (np.arange(10) * 2 + 1):
                u[i,j] -= 4 * (np.sinh(n*math.pi*x[i,j]) * np.cos(n*math.pi*y[i,j])) / ( (n*math.pi)**2 * np.sinh(2*n*math.pi ) )
    return u                   

def jaco():
    TOL = 1e-6
    res = 1
    MAXITE = 3000
    ite = 0
    resp = []
    
    #u = np.zeros((ni,nj))
    u = np.random.rand(ni,nj)
    u = bound(u,u)
    uk = u
    while res > TOL and ite <= MAXITE:
        ite += 1
        res = 0
        
        # Compute the new u(k+1)
        uk1 = np.zeros((ni,nj))
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                uk1[i,j] = ( (uk[i-1,j] + uk[i+1,j])/dx**2 + (uk[i,j-1] + uk[i,j+1])/dy**2 ) * (1 / (2/dx**2 + 2/dy**2) )
        
        # Compute Residual       
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                res += (uk1[i,j]-uk[i,j])**2
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        
        u = bound(u,uk1)
        resp.append(res)
        uk = uk1
    
    ufin = bound(uk1,uk1)
    
    return ufin, ite, resp
    
grid()
#### Analytical Solution
u = anal()
levels = np.arange(0.0,1.0,0.05)
cs = plt.contourf(x,y,u, levels = levels)
plt.colorbar(cs)
plt.title('Analytical Solution')
plt.show()

#### Jacobi
u, ite, resp = jaco()
print(' #### Jacobi ####')
print('iterations:', str(ite))
print('final residual: ', str(np.round(resp[-1],10)) )
levels = np.arange(0.0,1.0,0.05)
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Jacobi Solution')
plt.show()

#### residual plot
plt.plot(np.arange(1,ite+1), resp)
plt.title('Jacobi Solution')
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.show()


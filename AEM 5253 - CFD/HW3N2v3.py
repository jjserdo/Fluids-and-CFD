# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:33:14 2022

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# Initializing Everything
ni = 5 # 20
nj = 5 # 20
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
    u = np.zeros((ni,nj))
    #u = np.random.rand(ni,nj)
    u = bound(u,u)
    uk = u.copy()
    uk1 = u.copy()
    while res > TOL and ite < MAXITE:
        ite += 1
        res = 0
        # Compute the new u(k+1)
        uk1 = np.zeros((ni,nj))
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                uk1[i,j] = ( (uk[i-1,j] + uk[i+1,j])/dx**2 + (uk[i,j-1] + uk[i,j+1])/dy**2 ) * (1 / (2/dx**2 + 2/dy**2) )
        #print(uk)
        #print(uk1)
        # Compute Residual       
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                res += (uk1[i,j]-uk[i,j])**2
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        resp.append(res)
        uk = bound(uk1,uk1)
    ufin = bound(uk1,uk1) 
    return ufin, ite, resp
def gase():
    TOL = 1e-6
    res = 1
    MAXITE = 10
    ite = 0
    resp = []
    u = np.zeros((ni,nj))
    #u = np.random.rand(ni,nj)
    u = bound(u,u)
    u = bound(u,u)
    uk = u.copy()
    uk1 = u.copy()
    while res > TOL and ite < MAXITE:
        ite += 1
        res = 0
        for j in np.arange(1,nj-1): # goes from 2 to nj-1
            for i in np.arange(1,ni-1): # goes from 2 to ni-1
                #print(( (uk1[i-1,j] + uk1[i+1,j])/dx**2 + (uk1[i,j-1] + uk1[i,j+1])/dy**2 ) * (1 / (2/dx**2 + 2/dy**2) ))
                #print(uk1[i,j])
                uk1[i,j] = ( (uk1[i-1,j] + uk1[i+1,j])/dx**2 + (uk1[i,j-1] + uk1[i,j+1])/dy**2 ) * (1 / (2/dx**2 + 2/dy**2) )
                #print(uk1[i,j])
                #print('\n')
        # Compute Residual       
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                res += (uk1[i,j]-uk[i,j])**2
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        print(res)
        resp.append(res)
        uk = bound(uk1,uk1) 
    ufin = bound(uk1,uk1) 
    return ufin, ite, resp
def sor(w):
    TOL = 1e-6
    res = 1
    MAXITE = 3000
    ite = 0
    resp = []
    #u = np.zeros((ni,nj))
    u = np.random.rand(ni,nj)
    u = bound(u,u)
    uk1 = u
    while res > TOL and ite <= MAXITE:
        ite += 1
        res = 0
        uk = bound(uk1,uk1)
        uk1 = bound(uk1,uk1)
        # Compute the new u(k+1)
        for j in np.arange(1,nj-1): # goes from 2 to nj-1
            for i in np.arange(1,ni-1): # goes from 2 to ni-1
                uk1[i,j] = ( (uk[i-1,j] + uk[i+1,j])/dx**2 + (uk[i,j-1] + uk1[i,j+1])/dy**2 ) * (1 / (2/dx**2 + 2/dy**2) ) + (1-w)*uk[i,j]
        
        # Compute Residual       
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                res += (uk1[i,j]-uk[i,j])**2
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        resp.append(res)
    ufin = bound(uk1,uk1) 
    return ufin, ite, resp

grid()
#### Analytical Solution
u = anal()
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Analytical Solution')
plt.show()

'''
#### Jacobi ####
u, ite, resp = jaco()
print(' #### Jacobi ####')
print('Iterations to Convergence:', str(ite))
print('Final residual: ', str(np.round(resp[-1],10)) )
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Jacobi Solution')
plt.show()
plt.plot(np.arange(1,ite+1), resp)
plt.title('Jacobi Residuals')
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.yscale('log')
plt.show()
'''

#### Gauss-Seidel ####
u, ite, resp = gase()
print(' #### Gauss-Seidel ####')
print('Iterations to Convergence:', str(ite))
print('Final residual: ', str(np.round(resp[-1],10)) )
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Gauss-Seidel Solution')
plt.show()

#plt.plot(np.arange(1,ite+1), resp)
#plt.title('Gauss-Seidel Residuals')
#plt.xlabel('Iterations')
#plt.ylabel('Residual')
#plt.yscale('log')
#plt.show()

'''
#### Successive-Over Relaxation ###
wlist = np.arange(0.1,2,0.01)
itesor = np.zeros(len(wlist))
for i, w in enumerate(wlist):
    u, ite, resp = sor(w)
    itesor[i] = ite
plt.title('SOR Iterations')
plt.xlabel('w')
plt.ylabel('Iterations')
plt.show()
u, ite, resp = sor(1.95)
print(' #### Gauss-Seidel ####')
print('Iterations to Convergence:', str(ite))
print('Final residual: ', str(np.round(resp[-1],10)) )
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Gauss-Seidel Solution')
plt.show()
plt.plot(np.arange(1,ite+1), resp)
plt.title('Gauss-Seidel Residuals')
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.yscale('log')
plt.show()
'''

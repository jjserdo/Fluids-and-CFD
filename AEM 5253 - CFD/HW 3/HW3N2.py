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
    #print(ni,nj)
    for i in range(ni):
        arr1[i,0] = arr2[i,1]
    for j in range(nj):
        arr1[i,-1] = arr2[i,-2]
    arr1[0,:] = 0
    arr1[-1,:] = 1
    return arr1

def grid():
    # Initializing the Grid
    for i in np.arange(0,ni): # goes from 1 to ni
        for j in np.arange(0,nj): # goes from 1 to nj
            x[i,j] = (i) * 2.0 / (ni-1)
            y[i,j] = (j) * 1.0 / (nj-1)
#grid()
#print(y)
#print(y) 
#y = bound(y)
#print(y)

def anal():
    u = np.zeros((ni,nj))
    for i in np.arange(0,ni): # goes from 1 to ni
        for j in np.arange(0,nj): # goes from 1 to nj
            u[i,j] = x[i,j] / 4
            for n in (np.arange(10) * 2 + 1):
                u[i,j] += 4 * (np.sinh(n*math.pi*x[i,j]) * np.cos(n*math.pi*y[i,j])) / ( (n*math.pi)**2 * np.sinh(2*n*math.pi ) )
    return u                   

def jaco():
    u = np.random.rand(ni,nj)
    u = bound(u,u)
    uold = u
    
    tol = 1e-3
    res = 1
    maxite = 1000
    ite = 0
    resp = []
    
    while res > tol and ite <= maxite:
        ite += 1
        res = 0
        
        unew = np.zeros((ni,nj))
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                unew[i,j] = (u[i-1,j] - u[i+1,j])/dx**2 + (u[i,j-1] - u[i,j+1])/dy**2
                unew[i,j] = unew[i,j] * (1/(2/dx**2 + 2/dy**2))
                
                res += ( unew[i,j]-uold[i,j] )**2
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        
        #u = bound(u,unew)
        resp.append(res)
        u = unew
        uold = unew
    
    return unew, ite, resp
'''
def gase():
    u = np.random.rand(ni,nj)
    u = bound(u,u)
    
    tol = 1e-6
    res = 1
    maxite = 1000
    ite = 0
    resp = []
    
    while res > tol and ite <= maxite:
        ite += 1
        res = 0
        unew = np.zeros((ni,nj))
        for i in np.arange(1,ni-1): # goes from 2 to ni-1
            for j in np.arange(1,nj-1): # goes from 2 to nj-1
                unew[i,j] = (u[i-1,j] - u[i+1,j])/dx**2 + (u[i,j-1] - u[i,j+1])/dy**2
                unew[i,j] = unew[i,j] * (1/(2/dx**2 + 2/dy**2))
                
                res += (unew[i,j]-u[i,j])**2
                
        res = math.sqrt(res / ((ni-2)*(nj-2)) )
        u = bound(unew, unew)
        resp.append(res)
    return unew, ite, resp

def sory():
    pass
'''   
grid()
#print(x)
#print(y)

#### Analytical Solution
u = anal()
#print(u)
levels = np.arange(0.0,1.0,0.05)
cs = plt.contourf(x,y,u, levels = levels)
#plt.contour(x,y,u, levels = levels)
#plt.colorbar(cs)
plt.title('Analytical Solution')
plt.show()

#### Jacobi
u, ite, resp = jaco()
print(' #### Jacobi ####')
print('iterations:', str(ite))
print('final residual: ', str(np.round(resp[-1],10)) )
levels = np.arange(0.0,1.0,0.05)
cs = plt.contourf(x,y,u)
#plt.colorbar(cs)
plt.title('Jacobi Solution')
plt.show()

#### residual plot
plt.plot(np.arange(1,ite+1), resp)
plt.show()

'''
#### Gauss-Seidel
u, ite, resp = gase()
print(' #### Gauss-Seidel ####')
print('iterations:', str(ite))
print('final residual: ', str(np.round(resp[-1],10)) )
cs = plt.contour(x,y,u)
#plt.colorbar(cs)
plt.title('Gauss-Seidel Solution')
plt.show()

#### residual plot
plt.plot(np.arange(1,ite+1), resp)
plt.show()

#### SOR
wnum = np.linspace(0,2,10)
itee = np.zeros(len(wnum))
for ite, w in enumerate(wnum):
    #u, ite, resp = sorr(w)
    #itee[i] = ite
#print(u)
cs = plt.contourf(x,y,u)
plt.colorbar(cs)
plt.title('Jacobi Solution')
plt.show()

#### residual plot
#plt.plot(wnum, itee)
#plt.show()
'''





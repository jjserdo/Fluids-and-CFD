# -*- coding: utf-8 -*-
"""
Created using prof's notes in Tuesday 12/6 class

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# Initializing Everything
ni = 50 # number of cells in i
nj = 50 # number of cells in j
g = 1 # number of ghost cells
deta = (3.72-0.255) / (ni-1)
dxi =  (2*np.pi) / (nj-1)

eta = np.zeros((ni+2*g,nj+2*g))
xi = np.zeros((ni+2*g,nj+2*g))
x = np.zeros((ni+2*g,nj+2*g))
y = np.zeros((ni+2*g,nj+2*g))

Jdet = np.zeros((ni+2*g,nj+2*g))

dxidx = np.zeros((ni+2*g,nj+2*g))
dxidy = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))

# Notes: ghost cells are treated as 0 index
ib = g
jb = g
ie = ib + ni - 1
je = jb + nj - 1

# Initializing the Grid
for i in np.arange(ib,ie+1): # goes from ib to ie
    for j in np.arange(jb,je+1): # goes from jb to je
        eta[i,j] = (i-ib) * (3.72-0.255) / (ni-1) + 0.255
        xi[i,j] = (j-jb) * 2 * np.pi / (nj-1) + 0
for i in np.arange(ib,ie+1):
    for j in np.arange(jb,je+1):
        x[i,j] = np.cosh(eta[i,j]) * np.cos(xi[i,j])
        y[i,j] = np.sinh(eta[i,j]) * np.sin(xi[i,j])        

# Calculating the derivatives
for i in np.arange(ib,ie):
    for j in np.arange(jb,je):
        if i >=2 and i<= ni-1:
            # central
            dx_dxi = (x[i+1,j]-x[i-1,j]) / (2 * dxi)
            dy_dxi = (y[i+1,j]-y[i-1,j]) / (2 * dxi)
        elif i == 1:
            # forward
            dx_dxi = (x[i+1,j]-x[i,j]) / (dxi)
            dy_dxi = (y[i+1,j]-y[i,j]) / (dxi)
        elif i == ni:
            # forward
            dx_dxi = (x[i,j]-x[i-1,j]) / (dxi)
            dy_dxi = (y[i,j]-y[i-1,j]) / (dxi)
        
        if j >=2 and j<= nj-1:
            # central
            dx_deta = (x[i,j+1]-x[i,j-1]) / (2 * deta)
            dy_deta = (y[i,j+1]-y[i,j-1]) / (2 * deta)
        elif j == 1:
            # forward
            dx_deta = (x[i,j+1]-x[i,j]) / (deta)
            dy_deta = (y[i,j+1]-y[i,j]) / (deta)
        elif j == nj:
            # forward
            dx_deta = (x[i,j]-x[i,j-1]) / (deta)
            dy_deta = (y[i,j]-y[i,j-1]) / (deta)   
        Jdet[i,j] = 1/(dx_dxi*dy_deta-dx_deta*dy_dxi)
        dxidx[i,j] = dy_deta * Jdet[i,j]
        dxidy[i,j] = -dx_deta * Jdet[i,j]
        detadx[i,j] = -dy_dxi * Jdet[i,j]
        detady[i,j] = dx_dxi * Jdet[i,j]
        
#levels = [-1,0,2]
plt.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], Jdet[ib:ie+1,jb:je+1])
plt.xlim([-1.25,1.25])
plt.ylim([-1.25,1.25])
plt.show()

# compute dt 
dt = np.zeros((ni+2*g,nj+2*g))
for i in np.arange(ib,ie):
    for j in np.arange(jb,je):
        dfirst = np.abs(dxidx[i,j] * 0.4) + math.sqrt(dxidx[i,j]**2+dxidy[i,j]**2)
        dsecond = np.abs(detadx[i,j] * 0.4) + math.sqrt(detadx[i,j]**2+detady[i,j]**2)
        first = dxi / dfirst
        second = deta / dsecond
        dt[i,j] = min(first, second)
        
dT = np.min(dt[ib:ie,jb:je])
print(dT)
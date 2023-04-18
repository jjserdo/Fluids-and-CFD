# -*- coding: utf-8 -*-
"""
Created using prof's notes in Tuesday 12/6 class

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt
import math

#### Initializing Everything
# Parameters
ni = 100 # number of cells in i
nj = 50 # number of cells in j
g = 1 # number of ghost cells
dxi =  (2*np.pi) / (ni-1)
deta = (3.72-0.255) / (nj-1)

#dt = 0.014538054419394451 # calculations based on MeshGen2testJ
dt = 0.0144
 
# Grid
eta = np.zeros((ni+2*g,nj+2*g))
xi = np.zeros((ni+2*g,nj+2*g))
x = np.zeros((ni+2*g,nj+2*g))
y = np.zeros((ni+2*g,nj+2*g))

# Fluid Quantities
r = np.zeros((ni+2*g,nj+2*g))
ru = np.zeros((ni+2*g,nj+2*g))
rv = np.zeros((ni+2*g,nj+2*g))
re = np.zeros((ni+2*g,nj+2*g))
 
# Derivatives and Jacobian
dxidx = np.zeros((ni+2*g,nj+2*g))
dxidy = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))
Jdet = np.zeros((ni+2*g,nj+2*g))

# RHS 
rhs = np.zeros((ni+2*g,nj+2*g))

#### Initializing the Grid
# Notes: ghost cells are treated as 0 index
ib = g
jb = g
ie = ib + ni - 1
je = jb + nj - 1

# Initializing the Grid
for i in np.arange(ib,ie+1): # goes from ib to ie
    for j in np.arange(jb,je+1): # goes from jb to je
        xi[i,j] = (i-ib) * 2 * np.pi / (ni-1) 
        eta[i,j] = (j-jb) * (3.72-0.255) / (nj-1) + 0.255
        
for i in np.arange(ib,ie+1):
    for j in np.arange(jb,je+1):
        x[i,j] = np.cosh(eta[i,j]) * np.cos(xi[i,j])
        y[i,j] = np.sinh(eta[i,j]) * np.sin(xi[i,j])        
 
#### Calculating Derivatives and Jacobian
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
        
#### Initial Condition
# Set Initial Condition
r0 = 1
#u0 = [0.4, 0.85, 1.5, 6]
u0 = 0.4
v0 = 0
p0 = 1
g = 1.4
for i in np.arange(ib,ie+1): # goes from ib to ie
    for j in np.arange(jb,je+1): # goes from jb to je
        r[i,j] = r0
        ru[i,j] = r0 * u0
        rv[i,j] = r0 * v0
        q = (u0**2 + v0**2)/2
        re[i,j] = r0 * q + p0/(g-1)

#### Boundary Conditions
def bound():
    # eta last should be equal to eta first because periodic
        # so last equal first
        # r[:,-g-1] = r[:,1] ?
    # no slip flow but not sure how to implement
    # 
    return

def solvebound():
    # solves for the u and v values at the ellipse using the boundary conditions
    pass

#### Functions for Flux-Vector splitting

#### Solver, residual
TOL = 1e-6
res = 1
MAXITE = 3000
ite = 0
#while res > TOL and ite < MAXITE:
#    ite += 1
#    res = 0
    
    # Compute New Values of U
    #for i in np.arange(ib,ie): # I think I'm actually iterating from ib+1 to ie+1
    #    for j in np.arange(jb,je):
    #        rhs = 0
    #        rhs = fplus():
    #        dfdxi = fplus(u[i,j]) - fplus([u[i-1,j]])
    #        dfdxi += fminu(u[i+1,j]) - fminu([u[i,j]])
    #        dgdeta = gplus(u[i,j]) - gplus([u[i,j-1]])
    #        dgdeta += gminu(u[i,j+1]) - gminu([u[i,j]])
    #        rhs = dfdxi + dgdeta 
    #        unew = uold + dt * dgdeta
    
    # Compute Residual       
    #for i in np.arange(ib,ie):
    #    for j in np.arange(jb,je):
    #        res += (rk1[i,j]-r1[i,j])**2
    #res = math.sqrt(res / ((ni-2)*(nj-2)) )
    #resp.append(res)
    #rk = bound(rk1).copy() # Does this even make sense 

#### Post-processing, plotting contours
# Calculating umag
#umag  = np.zeros((ni+2*g,nj+2*g))
#for i in np.arange(ib,ie+1): # goes from ib to ie
#    for j in np.arange(jb,je+1): # goes from jb to je
#        um = u[i,j]**2 + v[i,j]**2
#        umag[i,j] = math.sqrt(um)

# Plotting the contours of Jacobian of the grid transformation
plt.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], Jdet[ib:ie+1,jb:je+1])
plt.xlim([-1.25,1.25])
plt.ylim([-1.25,1.25])
plt.show()


# Plotting the contours of the Velocity Magnitude
#plt.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1])
#plt.xlim([-1.25,1.25])
#plt.ylim([-1.25,1.25])
#plt.show()

##### Compute Drag 
# Manually note the drag coeff from each of the mach numbers from 0.4 to 1.5
A = 2 * x[1,1]

def computedrag():
    #idk how to compute the xcomp of pressure
    fd = 0
    for i in np.arange(ib,ie+1): 
        for j in np.arange(jb,je+1):
            pass 
            # How to calculate the x direction of pressure only
            # pressure is only perpedicular to the surface of the airfoil and can
            # get maybe an x comp from that,
            # but now how about anything that's not touching the airfoil?
            #fd = p * dx?
    cd = fd / (0.5 * r0 * u0**2 * A)
    return cd 

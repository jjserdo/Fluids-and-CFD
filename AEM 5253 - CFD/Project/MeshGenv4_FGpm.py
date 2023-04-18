# -*- coding: utf-8 -*-
"""
Created using prof's notes in Tuesday 12/6 class

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la

#### Initializing Everything
# Parameters
ni = 100 # number of cells in i
nj = 50 # number of cells in j
g = 1 # number of ghost cells
dxi =  (2*np.pi) / (ni-1)
deta = (3.72-0.255) / (nj-1)
dt = 0.0144
 
# Grid
eta = np.zeros((ni+2*g,nj+2*g))
xi = np.zeros((ni+2*g,nj+2*g))
x = np.zeros((ni+2*g,nj+2*g))
y = np.zeros((ni+2*g,nj+2*g))

# Fluid Quantities
U = np.zeros((ni+2*g,nj+2*g,4))    # rho, rhou, rhov, rhoE
ruvp = np.zeros((ni+2*g,nj+2*g,4)) # rho, u, v and pressure
 
# Derivatives and Jacobian
dxidx = np.zeros((ni+2*g,nj+2*g))
dxidy = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))
J = np.zeros((ni+2*g,nj+2*g))

# RHS 
rhs = np.zeros((ni+2*g,nj+2*g))

#### Initializing the Grid
# Notes: ghost cells are treated as 0 index
ib = g
jb = g
ie = ib + ni - 1
je = jb + nj - 1
# Initializing the Grid
def ell_grid():
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            xi[i,j] = (i-ib) * 2 * np.pi / (ni-1) 
            eta[i,j] = (j-jb) * (3.72-0.255) / (nj-1) + 0.255
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,je+1):
            x[i,j] = np.cosh(eta[i,j]) * np.cos(xi[i,j])
            y[i,j] = np.sinh(eta[i,j]) * np.sin(xi[i,j])        
 
#### Calculating Derivatives and Jacobian
def deri_jaco():
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,je+1):
            if i > ib and i < ie:
                # central
                dx_dxi = (x[i+1,j]-x[i-1,j]) / (2 * dxi)
                dy_dxi = (y[i+1,j]-y[i-1,j]) / (2 * dxi)
            elif i == ib:
                # forward
                dx_dxi = (x[i+1,j]-x[i,j]) / (dxi)
                dy_dxi = (y[i+1,j]-y[i,j]) / (dxi)
            elif i == ie:
                # backward
                dx_dxi = (x[i,j]-x[i-1,j]) / (dxi)
                dy_dxi = (y[i,j]-y[i-1,j]) / (dxi)
            
            if j > jb and j < je:
                # central
                dx_deta = (x[i,j+1]-x[i,j-1]) / (2 * deta)
                dy_deta = (y[i,j+1]-y[i,j-1]) / (2 * deta)
            elif j == jb:
                # forward
                dx_deta = (x[i,j+1]-x[i,j]) / (deta)
                dy_deta = (y[i,j+1]-y[i,j]) / (deta)
            elif j == je:
                # backward
                dx_deta = (x[i,j]-x[i,j-1]) / (deta)
                dy_deta = (y[i,j]-y[i,j-1]) / (deta)    
            J[i,j] = 1 / (dx_dxi*dy_deta-dx_deta*dy_dxi)
            dxidx[i,j] = dy_deta * J[i,j]
            dxidy[i,j] = -dx_deta * J[i,j]
            detadx[i,j] = -dy_dxi * J[i,j]
            detady[i,j] = dx_dxi * J[i,j]
        
#### Initial Condition
# Set Initial Condition
ga = 1.4
r0 = 1.0
#u0 = [0.4, 0.85, 1.5, 6]
u0 = 0.4
v0 = 0.0
p0 = 1.0 / ga
def ell_init():
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
        
            U[i,j,0] = r0
            U[i,j,1] = r0 * u0
            U[i,j,2] = r0 * v0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            E0 = r0q + p0/(ga-1)
            U[i,j,3] = E0
    ell_bound()
    get_actual()
    
# Function to get actual values
def get_actual():
    ruvp[:,:,0] = U[:,:,0]
    ruvp[:,:,1] = U[:,:,1] / U[:,:,0]
    ruvp[:,:,2] = U[:,:,2] / U[:,:,0]
    rq = 0.5 * (U[:,:,1]**2 + U[:,:,2]**2) / U[:,:,0]
    ruvp[:,:,3] = (U[:,:,3] - rq)*(ga-1)

#### Boundary Conditions
def ell_bound():
    U[ie,:,:] = U[ib,:,:]
    for i in np.arange(ib,ie+1):
        j = je
        U[i,j,0] = r0
        U[i,j,1] = r0 * u0
        U[i,j,2] = r0 * v0
        r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
        E0 = r0q + p0/(ga-1)
        U[i,je,3] = E0

def solvebound(i,j):
    # solves for the u and v values at the ellipse using the boundary conditions 
    A = np.zeros((2,2))
    B = np.zeros(2)
    A[0,0] = dxidx[i,j] / deta
    A[0,1] = dxidy[i,j] / deta
    A[1,0] = detadx[i,j]
    A[1,1] = detady[i,j]
    C = la.solve(A,B)
    return C
    
#### Functions for Flux-Vector splitting
# Compute Eigenvalues
def computelambdas(i,j,k1,k2):
    r = ruvp[i,j,0]
    u = ruvp[i,j,1] 
    v = ruvp[i,j,2]
    p = ruvp[i,j,3]
    c = math.sqrt(ga*p/r)
    
    l = np.zeros(4)
    l1 = k1 * u + k2 * v
    l2 = l1
    l3 = l1 + c * math.sqrt(k1**2+k2**2)
    l4 = l1 - c * math.sqrt(k1**2+k2**2)
    l[0] = l1
    l[1] = l2
    l[2] = l3
    l[3] = l4
    return l
# Return lambda plus and lambda minus
def lpm(i,j,k1,k2):
    l = computelambdas(i,j,k1,k2)
    ep = 0.01
    lp = np.zeros(4)
    lm = np.zeros(4)
    for i in range(4):
        lp[i] = 0.5 * (l[i] + math.sqrt(l[i]**2 + ep ** 2))
        lm[i] = 0.5 * (l[i] - math.sqrt(l[i]**2 + ep ** 2))
    return lp, lm, l
# F and G plus and minus
def FGpm(l,i,j,k1,k2):
    l1 = l[0]
    l3 = l[2]
    l4 = l[3]
    
    r = ruvp[i,j,0]
    u = ruvp[i,j,1] 
    v = ruvp[i,j,2]
    p = ruvp[i,j,3]
    
    c = math.sqrt(ga*p/r)
    
    k1b = k1/math.sqrt(k1**2+k2**2)
    k2b = k2/math.sqrt(k1**2+k2**2)
    
    FGpm = np.zeros(4)
    FGpm[0] = (0.5*r/ga) * (2*(ga-1)*l1 + l3 + l4)
    FGpm[1] = (0.5*r/ga) * (2*(ga-1)*l1*u + l3*(u+c*k1b) + l4*(u-c*k1b))
    FGpm[2] = (0.5*r/ga) * (2*(ga-1)*l1*v + l3*(v+c*k2b) + l4*(v-c*k2b))
    W = ( (3-ga)*(l3+l4)*c**2 )/( 2*(ga-1) )
    FGpm[3] = (0.5*r/ga) * ((ga-1)*l1*(u**2+v**2) + 0.5*l3*((u+c*k1b)**2+(v+c*k2b)**2) + 0.5*l4*((u-c*k1b)**2+(v-c*k2b)**2) + W)
    return FGpm

#### Solver, residual
ti = 0
tf = 100
itei = 0
maxite = 500
resi = 1
tol = 1e-6
def solve():
    t = ti
    ite = itei
    res = resi
    get_actual()
    while t <= tf and ite < maxite and res > tol: # and res > tol
        t += dt
        ite += 1
        ell_bound()
        rhoold = U[:,:,0].copy()
        for i in np.arange(ib,ie+1): # ie
            for j in np.arange(jb,je+1): # je
                # Do Explicit Euler to find update u vector
                rhs = np.zeros(4)
                # Evaluate df/dxi
                k1 = dxidx[i,j] / J[i,j]
                k2 = dxidy[i,j] / J[i,j]
                lp, lm, l = lpm(i,j,k1,k2)
                rhs += FGpm(lp,i  ,j,k1,k2) / dxi
                rhs -= FGpm(lp,i-1,j,k1,k2) / dxi
                rhs += FGpm(lm,i+1,j,k1,k2) / dxi
                rhs -= FGpm(lm,i  ,j,k1,k2) / dxi
                
                # Evaluate dg/deta
                k1 = detadx[i,j] / J[i,j]
                k2 = detady[i,j] / J[i,j]
                lp, lm, l = lpm(i,j,k1,k2)
                rhs += FGpm(lp,i,j  ,k1,k2) / deta
                rhs -= FGpm(lp,i,j-1,k1,k2) / deta
                rhs += FGpm(lm,i,j+1,k1,k2) / deta
                rhs -= FGpm(lm,i,j  ,k1,k2) / deta
                U[i,j] -= rhs * dt * J[i,j]
        # Compute Residual
        res = 0
        for i in np.arange(ib,ie):
            for j in np.arange(jb,je):
                res += (U[i,j,0]-rhoold[i,j])**2
        res = math.sqrt( res / (ni*nj) )
        get_actual()
    
    #print('Solution Solved lmao')
    print('Final time at', str(np.round(t,5)))
    print('Final iteration at', str(ite))
    print('Final Residual at', str(np.round(res,5)))

#### Post-processing, plotting contours
def get_umag():
    umag  = np.zeros((ni+2*g,nj+2*g))
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            um = ruvp[i,j,1]**2 + ruvp[i,j,2]**2
            umag[i,j] = math.sqrt(um)
    return umag
    

# Plotting the contours of Jacobian of the grid transformation
def plot_jaco():
    fig, ax = plt.subplots()
    ax.set_title('Jacobian Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-1.25,1.25])
    ax.set_ylim([-1.25,1.25])
    #levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], J[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    plt.show() 

# Plotting the contours of the Velocity Magnitude
def plot_umag():
    umag = get_umag()
    fig, ax = plt.subplots()
    ax.set_title('Velocity Magnitude Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-1.25,1.25])
    ax.set_ylim([-1.25,1.25])
    #levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    plt.show() 

##### Compute Drag 
def compute_area():
    A = 0
    return A
def compute_drag():
    A = compute_area()
    #idk how to compute the xcomp of pressure
    fd = 0
    for i in np.arange(ib,ie+1): 
        for j in np.arange(jb,je+1):
            pass 
            # How to calculate the x direction of pressure only
            # pressure is only perpedicular to the surface of the airfoil and can
            # get maybe an x comp from that
            # but now how about anything that's not touching the airfoil?
            #fd = p * dx?
    cd = fd / (0.5 * r0 * u0**2 * A)
    return cd 

##############################################
def main():
    ell_grid()
    deri_jaco()
    #ell_init()
    plot_jaco()
    
    #solve()
    plot_umag()
    

if __name__ == "__main__":
    main()
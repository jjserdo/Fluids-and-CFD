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
ni = 56 # number of cells in i
nj = 2 # number of cells in j
g = 1 # number of ghost cells
#dxi =  (4.5-0) / (ni-1)
#deta = (1-0) / (nj-1)
dxi =  (4.5-0) / (ni-1)
deta = (1-0) / (nj-1)
dt = 0.032

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
Jdet = np.zeros((ni+2*g,nj+2*g))

# RHS 
rhs = np.zeros((ni+2*g,nj+2*g))

#### Initializing the Grid
# Notes: ghost cells are treated as 0 index
ib = g
jb = g
ie = ib + ni - 1
je = jb + nj - 1
def sod_grid():
    # Set Initial Condition
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            xi[i,j] = (i-ib) * (4.5-0) / (ni-1) + 0
            eta[i,j] = (j-jb) * (1-0) / (nj-1) + 0        
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            x[i,j] = xi[i,j]
            y[i,j] = eta[i,j]
 
#### Calculating Derivatives and Jacobian
def deri_jaco():
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
ga = 1.4
r0 = 1.0
rf = 0.125
u0 = 0.0
v0 = 0.0
p0 = 1.0 / ga
pf = 10.0 * p0
c = 1
def sod_init():
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            U[i,j,1] = r0 * u0
            U[i,j,2] = r0 * v0

            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            rfq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            
            E0 = r0q + p0/(ga-1)
            Ef = rfq + pf/(ga-1)
            
            if i*dxi <= 1.95:
                U[i,j,0] = r0
                U[i,j,3] = E0
            else:
                U[i,j,0] = rf
                U[i,j,3] = Ef
    sod_bound()
    get_actual()
                
#### Boundary Conditions
def sod_bound():
    U[ :, 0,:] = U[ :,jb,:]
    U[ :,-1,:] = U[ :,je,:]
    U[ 0, :,:] = U[ib, :,:]
    U[-1, :,:] = U[ie, :,:]
def noslip():
    pass

#### Functions for upwinding
# Eigenvalues 
def get_actual():
    ruvp[:,:,0] = U[:,:,0]
    ruvp[:,:,1] = U[:,:,1] / U[:,:,0]
    ruvp[:,:,2] = U[:,:,2] / U[:,:,0]
    q = 0.5 * (ruvp[:,:,1]**2 + ruvp[:,:,2]**2)
    ruvp[:,:,3] = (U[:,:,3] - U[:,:,0] * q)*(ga-1)

# F and G plus and minus
def FGpm(l,i,j,k1,k2):
    l1 = l[0]
    l3 = l[2]
    l4 = l[3]
    u = U[i,j,1]
    v = U[i,j,2]
    k1b = k1/math.sqrt(k1**2+k2**2)
    k2b = k2/math.sqrt(k1**2+k2**2)
    FGpm = np.zeros(4)
    FGpm[0] = 2*(ga-1)*l1 + l3 + l4
    FGpm[1] = 2*(ga-1)*l1*u + l3*(u+c*k1b) + l4*(u-c*k1b)
    FGpm[2] = 2*(ga-1)*l1*v + l3*(u+c*k2b) + l4*(u-c*k2b)
    W = ((3-ga)*(l3+l4)*c**2)/(2*(ga-1))
    FGpm[3] = (ga-1)*l1*(u**2+v**2) + 0.5*l3*((u+c*k1b)**2+(v+c*k2b)**2) + 0.5*l4*((u-c*k1b)**2+(v-c*k2b)**2) + W
    return FGpm

def computelambdas(i,j,k1,k2):
    l = np.zeros(4)
    l1 = k1 * ruvp[i,j,1] + k2 * ruvp[i,j,2]
    l2 = l1
    l3 = l1 + c * math.sqrt(k1**2+k2**2)
    l4 = l1 - c * math.sqrt(k1**2+k2**2)
    l[0] = l1
    l[1] = l2
    l[2] = l3
    l[3] = l4
    return l

def lpm(i,j,k1,k2):
    l = computelambdas(i,j,k1,k2)
    ep = 0.01
    lp = np.zeros(4)
    lm = np.zeros(4)
    for i in range(4):
        lp[i] = 0.5 * (l[i] + math.sqrt(l[i]**2 + ep ** 2))
        lm[i] = 0.5 * (l[i] + math.sqrt(l[i]**2 + ep ** 2))
    return lp, lm
    
#### Solver, residual
ti = 0
tf = 1.0
itei = 0
maxite = 1
resi = 1
tol = 1e-6
def solve():
    t = ti
    ite = itei
    res = resi
    while t <= tf and ite < maxite and res > tol:
        t += dt
        ite += 1
        sod_bound()
        rhoold = U[:,:,0]
        #print(ruvp[:,:,0])
        for i in np.arange(ib,ie):
            for j in np.arange(jb,je):
                # Do Explicit Euler to find update u vector
                rhs = np.zeros(4)
    
                # Evaluate df/dxi
                lp, lm = lpm(i,j,1,0)
                rhs += FGpm(lp,i  ,j,1,0) / dxi
                rhs -= FGpm(lp,i-1,j,1,0) / dxi
                rhs += FGpm(lm,i+1,j,1,0) / dxi
                rhs -= FGpm(lm,i  ,j,1,0) / dxi
                # Evaluate dg/deta
                lp, lm = lpm(i,j,0,1)
                rhs += FGpm(lp,i  ,j,0,1) / dxi
                rhs -= FGpm(lp,i-1,j,0,1) / dxi
                rhs += FGpm(lm,i+1,j,0,1) / dxi
                rhs -= FGpm(lm,i  ,j,0,1) / dxi
                #print(rhs) # got 0 wat
                
                U[i,j] += rhs * dt * Jdet[i,j]

        get_actual()
        #print(ruvp[:,:,0])
        
        # Compute Residual
        res = 0
        for i in np.arange(ib,ie):
            for j in np.arange(jb,je):
                res += (ruvp[i,j,0]-rhoold[i,j])**2
        res = math.sqrt(res / (ni*nj) )
    
    print('Solution Solved lmao')
    print('Final time at', str(t))
    print('Final iteration at', str(ite))
    print('Final Residual at', str(res))

#### Post-processing, plotting contours
        
# Plot of density distribution distribution inside the shock tube
def plot_density():
    fig, ax = plt.subplots()
    ax.set_title('rho/rho0 vs x')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([0,4.5])
    ax.set_ylim([0,1.25])
    ax.plot(x[ib:ie+1,je], ruvp[ib:ie+1,je,0],'bo')

# Contour Plot of Jacobian
def plot_jaco():
    fig, ax = plt.subplots()
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], Jdet[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    plt.show()  # Jdet.min() = -1, Jdet.max() = 0

def main():
    sod_grid()
    deri_jaco()
    sod_init()
    plot_jaco()
    solve()
    

if __name__ == "__main__":
    main()
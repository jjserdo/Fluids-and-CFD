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
J = np.zeros((ni+2*g,nj+2*g))

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
    Bound(dxidx)
    Bound(dxidy)
    Bound(detadx)
    Bound(detady)
    Bound(J)
    
        
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
    Bound(U)
    get_actual()
                
#### Boundary Conditions
def Bound(arr):
    arr[ :, 0] = arr[ :,jb]
    arr[ :,-1] = arr[ :,je]
    arr[ 0, :] = arr[ib, :]
    arr[-1, :] = arr[ie, :]
    return arr
def sod_bound():
    U[ :, 0,:] = U[ :,jb,:]
    U[ :,-1,:] = U[ :,je,:]
    U[ 0, :,:] = U[ib, :,:]
    U[-1, :,:] = U[ie, :,:]
def noslip():
    pass
def wallbound(i,j):
    U[i,j,2] = 0 
    

#### Functions for upwinding
def get_actual():
    ruvp[:,:,0] = U[:,:,0]
    ruvp[:,:,1] = U[:,:,1] / U[:,:,0]
    ruvp[:,:,2] = U[:,:,2] / U[:,:,0]
    q = 0.5 * (ruvp[:,:,1]**2 + ruvp[:,:,2]**2)
    ruvp[:,:,3] = (U[:,:,3] - U[:,:,0] * q)*(ga-1)

# Compute Eigenvalues
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
        lm[i] = 0.5 * (l[i] - math.sqrt(l[i]**2 + ep ** 2))
    return lp, lm, l

# F and G plus and minus
def FGpm(i,j,FG,PM):
    if FG == 'F':
        k1 = dxidx[i,j] / J[i,j]
        k2 = dxidy[i,j] / J[i,j]
    elif FG == 'G':
        k1 = detadx[i,j] / J[i,j]
        k2 = detady[i,j] / J[i,j]
        
    lp, lm, L = lpm(i,j,k1,k2)
    
    if PM == 'p':
        l = lp
    elif PM == 'm':
        l = lm

    l1 = l[0]
    l3 = l[2]
    l4 = l[3]
    r = U[i,j,0]
    u = U[i,j,1]
    v = U[i,j,2]
    
    k1b = k1/math.sqrt(k1**2+k2**2)
    k2b = k2/math.sqrt(k1**2+k2**2)

    FGpm = np.zeros(4)
    FGpm[0] = (0.5*r/ga) * (2*(ga-1)*l1 + l3 + l4)
    FGpm[1] = (0.5*r/ga) * (2*(ga-1)*l1*u + l3*(u+c*k1b) + l4*(u-c*k1b))
    FGpm[2] = (0.5*r/ga) * (2*(ga-1)*l1*v + l3*(u+c*k2b) + l4*(u-c*k2b))
    W = ((3-ga)*(l3+l4)*c**2)/(2*(ga-1))
    FGpm[3] = (0.5*r/ga) * ((ga-1)*l1*(u**2+v**2) + 0.5*l3*((u+c*k1b)**2+(v+c*k2b)**2) + 0.5*l4*((u-c*k1b)**2+(v-c*k2b)**2) + W)
    return FGpm

#### Solver, residual
ti = 0
tf = 1.0
itei = 0
maxite = 10
resi = 1
tol = 1e-6
def solve():
    t = ti
    ite = itei
    res = resi
    while t <= tf and ite < maxite: # and res > tol
        t += dt
        ite += 1
        sod_bound()
        rhoold = U[:,:,0].copy()
        for i in np.arange(ib,ie+1): # ie
            for j in np.arange(jb,je+1): # je
                # Do Explicit Euler to find update u vector
                rhs = np.zeros(4)
                
                # Evaluate df/dxi
                rhs += FGpm(i  ,j,'F','p') / dxi 
                rhs -= FGpm(i-1,j,'F','p') / dxi 
                rhs += FGpm(i+1,j,'F','m') / dxi 
                rhs -= FGpm(i  ,j,'F','m') / dxi
                
                # Evaluate dg/deta
                rhs += FGpm(i,j  ,'G','p') / deta
                rhs -= FGpm(i,j-1,'G','p') / deta
                rhs += FGpm(i,j+1,'G','m') / deta
                rhs -= FGpm(i,j  ,'G','m') / deta
                
                U[i,j] += rhs * dt * J[i,j]
                #if i==24 and j==jb:
                    #print(U[i,j])

        #get_actual()
        #print(ruvp[:,:,0])
        
        # Compute Residual
        res = 0
        for i in np.arange(ib,ie):
            for j in np.arange(jb,je):
                res += (U[i,j,0]-rhoold[i,j])**2
        res = math.sqrt(res / (ni*nj) )
        #print(rhoold[:,:])
        #print(U[:,:,0])
        #print(res)
    
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
    ax.set_ylim([0,1])
    ax.plot(x[ib:ie+1,je], ruvp[ib:ie+1,je,0],'bo')

# Contour Plot of Jacobian
def plot_jaco():
    fig, ax = plt.subplots()
    ax.set_title('Jacobian Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([0,4.5])
    ax.set_ylim([0,1])
    levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], J[ib:ie+1,jb:je+1], levels=levels)
    cbar = fig.colorbar(cf)
    plt.show()  # Jdet.min() = -1, Jdet.max() = 0

##############################################
def main():
    sod_grid()
    deri_jaco()
    sod_init()
    solve()
    get_actual()
    plot_density()
    plot_jaco()
    #plot_umag() # umag looks suspicious
    #get_actual()
    

if __name__ == "__main__":
    main()
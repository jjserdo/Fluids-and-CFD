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
ni = 51 #number of cells in i
nj = 2 # number of cells in j
g = 1 # number of ghost cells
dxi =  (1-0) / (ni-1)
deta = (1-0) / (nj-1)
dt = dxi * 0.4
print('dx is', str(np.round(dxi,5)))
print('dt is', str(np.round(dt ,5)))

# Grid
eta = np.zeros((ni+2*g,nj+2*g))
xi  = np.zeros((ni+2*g,nj+2*g))
x   = np.zeros((ni+2*g,nj+2*g))
y   = np.zeros((ni+2*g,nj+2*g))

# Fluid Quantities
U    = np.zeros((ni+2*g,nj+2*g,4)) # rho, rhou, rhov, rhoE
ruvp = np.zeros((ni+2*g,nj+2*g,4)) # rho, u, v and pressure
 
# Derivatives and Jacobian
dxidx  = np.zeros((ni+2*g,nj+2*g))
dxidy  = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))
J      = np.zeros((ni+2*g,nj+2*g))

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
            xi[i,j]  = (i-ib) * (1-0) / (ni-1) + 0
            eta[i,j] = (j-jb) * (1-0) / (nj-1) + 0        
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            x[i,j] =  xi[i,j]
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
            J[i,j] = 1 / (dx_dxi*dy_deta-dy_dxi*dx_deta)
            dxidx[i,j] = dy_deta * J[i,j]
            dxidy[i,j] = -dx_deta * J[i,j]
            detadx[i,j] = -dy_dxi * J[i,j]
            detady[i,j] = dx_dxi * J[i,j]
    Bound(J)
    Bound(dxidx)
    Bound(dxidy)
    Bound(detadx)
    Bound(detady)        
#### Initial Condition
# Set Initial Condition
ga = 1.4
r0 = 1.0
rf = 0.125
p0 = 1.0
pf = 0.1 * p0
u0 = 0.0
v0 = 0.0
def sod_init():
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,je+1): #np.arange(jb,je+1): # goes from jb to je
            U[i,j,1] = r0 * u0
            U[i,j,2] = r0 * v0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            rfq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0         
            E0 = r0q + p0/(ga-1)
            Ef = rfq + pf/(ga-1)
            if i*dxi <= 0.5:
                U[i,j,0] = r0
                U[i,j,3] = E0
            else:
                U[i,j,0] = rf
                U[i,j,3] = Ef
    Ubound_up()
    Ubound_down()
    Ubound_left()
    Ubound_right()
    get_actual()
#### Boundary Conditions
def Bound(arr):
    arr[ :, 0] = arr[ :,jb]
    arr[ :,-1] = arr[ :,je]
    arr[ 0, :] = arr[ib, :]
    arr[-1, :] = arr[ie, :]
    return arr 
def Ubound_up():
    U[ :,-1,:] = U[ :,je,:]
def Ubound_down():
    U[ :, 0,:] = U[ :,jb,:]
def Ubound_left():
    U[ 0, :,:] = U[ib, :,:]
def Ubound_right():
    U[-1, :,:] = U[ie, :,:] 

#### Functions for upwinding
def get_actual():
    ruvp[:,:,0] = U[:,:,0]
    ruvp[:,:,1] = U[:,:,1] / U[:,:,0]
    ruvp[:,:,2] = U[:,:,2] / U[:,:,0]
    rq = 0.5 * ( U[:,:,1]**2 + U[:,:,2]**2 ) / U[:,:,0]
    ruvp[:,:,3] = ( U[:,:,3] - rq ) * ( ga - 1 )

# Compute Eigenvalues
def lpm(i,j,k1,k2):
    r = ruvp[i,j,0]
    u = ruvp[i,j,1] 
    v = ruvp[i,j,2]
    p = ruvp[i,j,3]
    c = math.sqrt(ga*p/r)
    l1 = k1 * u + k2 * v
    l2 = l1
    l3 = l1 + c * math.sqrt(k1**2+k2**2)
    l4 = l1 - c * math.sqrt(k1**2+k2**2)
    l = np.zeros(4)
    l[0] = l1
    l[1] = l2
    l[2] = l3
    l[3] = l4
    ep = 0.1 * c * math.sqrt(k1**2+k2**2)
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
    W = ((3-ga)*(l3+l4)*c**2)/(2*(ga-1))
    FGpm[3] = (0.5*r/ga) * ((ga-1)*l1*(u**2+v**2) + 0.5*l3*((u+c*k1b)**2+(v+c*k2b)**2) + 0.5*l4*((u-c*k1b)**2+(v-c*k2b)**2) + W)
    return FGpm

#### Solver, residual
ti = 0
tf = 0.2
itei = 0
maxite = 1000
resi = 1
tol = 1e-6
def solve():
    t = ti
    ite = itei
    res = resi
    get_actual()
    while t <= tf and ite < maxite: # and res > tol
        t += dt
        ite += 1
        Ubound_up()
        Ubound_down()
        rhoold = U[:,:,0].copy()
        rhs = np.zeros((ni+2*g,nj+2*g,4))
        for i in np.arange(ib,ie+1): # ie
            for j in np.arange(jb,je+1): # je
                # Do Explicit Euler to find update u vector
                rhs[i,j] = np.zeros(4)
                # Evaluate df/dxi
                rhs[i,j] += FGpm(i  ,j,'F','p') / dxi 
                rhs[i,j] -= FGpm(i-1,j,'F','p') / dxi 
                rhs[i,j] += FGpm(i+1,j,'F','m') / dxi 
                rhs[i,j] -= FGpm(i  ,j,'F','m') / dxi
                # Evaluate dg/deta
                rhs[i,j] += FGpm(i,j  ,'G','p') / deta
                rhs[i,j] -= FGpm(i,j-1,'G','p') / deta
                rhs[i,j] += FGpm(i,j+1,'G','m') / deta
                rhs[i,j] -= FGpm(i,j  ,'G','m') / deta
        for i in np.arange(ib,ie+1): # ie
            for j in np.arange(jb,je+1): # je
                U[i,j] -= rhs[i,j] * dt * J[i,j]       
        # Compute Residual
        res = 0
        for i in np.arange(ib,ie):
            for j in np.arange(jb,je):
                res += (U[i,j,0]-rhoold[i,j])**2
        res = math.sqrt(res / ((ni)*nj) )
        get_actual()
    
    #print('Solution Solved lmao')
    print('Final time at', str(np.round(t,5)))
    print('Final iteration at', str(ite))
    print('Final Residual at', str(np.round(res,5)))

#### Post-processing, plotting contours       
# Plot of density distribution distribution inside the shock tube
def plot_density():
    fig, ax = plt.subplots()
    ax.set_title('rho/rho0 vs x')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_xlim([0,4.5])
    #ax.set_ylim([0,1])
    ax.plot(x[ib:ie+1,je], ruvp[ib:ie+1,je,0],'bo')
def plot_pressure():
    fig, ax = plt.subplots()
    ax.set_title('pressure vs x')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_xlim([0,4.5])
    #ax.set_ylim([0,1])
    ax.plot(x[ib:ie+1,je], ruvp[ib:ie+1,je,3],'bo')

# Contour Plot of Jacobian
def plot_jaco():
    fig, ax = plt.subplots()
    ax.set_title('Jacobian Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_xlim([0,4.5])
    #ax.set_ylim([0,1])
    levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], J[ib:ie+1,jb:je+1], levels=levels)
    cbar = fig.colorbar(cf)
    plt.show()  
    
# Plotting the contours of the Velocity Magnitude
def get_umag():
    umag  = np.zeros((ni+2*g,nj+2*g))
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            um = ruvp[i,j,1]**2 + ruvp[i,j,2]**2
            umag[i,j] = math.sqrt(um)
    return umag
def plot_umag():
    umag = get_umag()
    fig, ax = plt.subplots()
    ax.set_title('Velocity Magnitude Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    plt.show() 

##############################################
def main():
    sod_grid()
    deri_jaco()
    sod_init()
    plot_jaco()
    
    solve()
    plot_density()
    plot_pressure()
    plot_umag()
    

if __name__ == "__main__":
    main()
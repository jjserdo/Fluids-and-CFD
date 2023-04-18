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
u = np.zeros((ni+2*g,nj+2*g,4))
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
p0 = 1.0/ ga
pf = 10.0 * p0
c = 1
for i in np.arange(ib,ie+1): # goes from ib to ie
    for j in np.arange(jb,je+1): # goes from jb to je
        u[i,j,1] = r0 * u0
        u[i,j,2] = r0 * v0
        
        ruvp[i,j,0] = u[i,j,0]
        ruvp[i,j,1] = u[i,j,1] / u[i,j,0]
        ruvp[i,j,2] = u[i,j,2] / u[i,j,0]
        q = 0.5 * (ruvp[i,j,1]**2 + ruvp[i,j,2]**2)
        
        E0 = r0*q + p0/(ga-1)
        Ef = rf*q + pf/(ga-1)
        
        if i*dxi >= 1.95:
            u[i,j,0] = r0
            u[i,j,3] = E0
        else:
            u[i,j,0] = rf
            u[i,j,3] = Ef
            
        ruvp[i,j,3] = (u[i,j,3] - u[i,j,0] * q)*(ga-1)

#### Boundary Conditions
# tbh idk 

#### Functions for upwinding
# Eigenvalues 

ruvp[i,j,0] = u[i,j,0]
ruvp[i,j,1] = u[i,j,1] / u[i,j,0]
ruvp[i,j,2] = u[i,j,2] / u[i,j,0]
q = 0.5 * (ruvp[i,j,1]**2 + ruvp[i,j,2]**2)
ruvp[i,j,3] = (u[i,j,3] - u[i,j,0] * q)*(ga-1)

l1 = dxidx[i,j] * ruvp[i,j,1] + dxidy[i,j] * ruvp[i,j,2]
l2 = l1
l3 = l1 + c * math.sqrt(l1**2+l2**2)
l4 = l1 - c * math.sqrt(l1**2+l2**2)
# F and G plus and minus
def FGpm(l1,l2,l3,l4,u,v,k1,k2):
    k1b = k1/math.sqrt(k1**2+k2**2)
    k2b = k2/math.sqrt(k1**2+k2**2)
    FGpm = np.zeros(4)
    FGpm[0] = 2*(ga-1)*l1 + l3 + l4
    FGpm[1] = 2*(ga-1)*l1*u + l3*(u+c*k1b) + l4*(u-c*k1b)
    FGpm[2] = 2*(ga-1)*l1*v + l3*(u+c*k2b) + l4*(u-c*k2b)
    W = ((3-ga)*(l3+l4)*c**2)/(2*(ga-1))
    FGpm[3] = (ga-1)*l1*(u**2+v**2) + 0.5*l3*((u+c*k1b)**2+(v+c*k2b)**2) + 0.5*l4*((u-c*k1b)**2+(v-c*k2b)**2) + W

def computelambdas(i,j):
    u = ruvp[i,j,1]
    v = ruvp[i,j,2]
    l1 = dxidx[i,j] * ruvp[i,j,1] + dxidy[i,j] * ruvp[i,j,2]
    l2 = l1
    l3 = l1 + c * math.sqrt(l1**2+l2**2)
    l4 = l1 - c * math.sqrt(l1**2+l2**2)
    return l1, l2, l3, l4

def getuv(i,j):
    u = u[i,j,1]
    v = u[i,j,2]
    return u,v 

def lpm(l):
    ep = 0.01
    lp = 0.5 * (l + math.sqrt(l**2 + ep ** 2))
    lm = 0.5 * (l + math.sqrt(l**2 + ep ** 2))
    return lp, lm

def bound():
    pass
    
#### Solver, residual
t = 0
tf = 1.0
while t <= tf:
    t += dt
    bound()
    rhoold = u[i,j,0]
    for i in np.arange(ib,ie):
        for j in np.arange(jb,je):
            # Do Explicit Euler to find update u vector
            rhs = np.zeros(4)
            # what is u and v here for now?
            #rhs = FGpm(l1,l2,l3,l4,u,v,1,0) # generic F and G
            # set l1,l2,l3 and l4, u and v
            l1,l2,l3,l4,u,v = computelambdas(i,j)
            rhs += FGpm(l1,l2,l3,l4,u,v,1,0) # Fplus
            l1,l2,l3,l4,u,v = computelambdas(i-1,j)
            rhs -= FGpm(l1,l2,l3,l4,u,v,1,0) # Fplus 
            l1,l2,l3,l4,u,v = computelambdas(i+1,j)
            rhs += FGpm(l1,l2,l3,l4,u,v,1,0) # Fminus
            l1,l2,l3,l4,u,v = computelambdas(i,j)
            rhs -= FGpm(l1,l2,l3,l4,u,v,1,0) # Fminus
            l1,l2,l3,l4,u,v = computelambdas(i,j)
            rhs += FGpm(l1,l2,l3,l4,u,v,0,1) # Gplus
            l1,l2,l3,l4,u,v = computelambdas(i,j-1)
            rhs -= FGpm(l1,l2,l3,l4,u,v,0,1) # Gplus 
            l1,l2,l3,l4,u,v = computelambdas(i,j+1)
            rhs += FGpm(l1,l2,l3,l4,u,v,0,1) # Gminus
            l1,l2,l3,l4,u,v = computelambdas(i,j)
            rhs -= FGpm(l1,l2,l3,l4,u,v,0,1) # Gminus
            u[i,j] += rhs * dt * Jdet[i,j]
            
            # Store rho, uvel, vvel, pres from u vector
            ruvp[i,j,0] = u[i,j,0]
            ruvp[i,j,1] = u[i,j,1] / u[i,j,0]
            ruvp[i,j,2] = u[i,j,2] / u[i,j,0]
            q = 0.5 * (ruvp[i,j,1]**2 + ruvp[i,j,2]**2)
            ruvp[i,j,3] = (u[i,j,3] - u[i,j,0] * q)*(ga-1)
            
    # Compute Residual
    res = 0
    for i in np.arange(ib,ie):
        for j in np.arange(jb,je):
            res += (ruvp[i,j,0]-rhoold[i,j])**2
    res = math.sqrt(res / ((ni-2)*(nj-2)) )
    
    
print('Solution Solved lmao')
print('Final time at', str(t))

#### Post-processing, plotting contours
# Computing rho, u vel, v vel, vel mag and pressure
rho = np.zeros((ni+2*g,nj+2*g))
uvel = np.zeros((ni+2*g,nj+2*g))
vvel = np.zeros((ni+2*g,nj+2*g))
velmag = np.zeros((ni+2*g,nj+2*g))
pres = np.zeros((ni+2*g,nj+2*g))
for i in np.arange(ib,ie+1): # goes from ib to ie
    for j in np.arange(jb,je+1): # goes from jb to je
        rho[i,j] = u[i,j,0]
        uvel[i,j] = u[i,j,1] / u[i,j,0]
        vvel[i,j] = u[i,j,2] / u[i,j,0]
        velmag[i,j] = math.sqrt(uvel[i,j]**2 + vvel[i,j]**2)
        quack = 0.5 * rho[i,j] * (uvel[i,j]**2 + vvel[i,j]**2)
        vvel[i,j] = ( u[i,j,3] - quack ) * (ga-1)
        
# Plot of density distribution distribution inside the shock tube
fig, ax = plt.subplots()
ax.set_title('rho/rho0 vs x')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim([0,4.5])
ax.set_ylim([0,1.25])
ax.plot(x[ib:ie+1,je], rho[ib:ie+1,je],'bo')

# Contour Plot of Jacobian
fig, ax = plt.subplots()
cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], Jdet[ib:ie+1,jb:je+1])
cbar = fig.colorbar(cf)
plt.show()  # Jdet.min() = -1, Jdet.max() = 0

'''
def main():
    grid()

if __name__ = '__main__':
    main():
'''
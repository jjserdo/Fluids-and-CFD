import os
os.chdir('C:/Users/jjser/OneDrive - Syracuse University/Fall 2022/AEM 5253/Project')

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la
from numpy import asarray
import pandas
import time

#### Initializing Everything
# Parameters
ni = 51 # number of cells in i, make this odd so Cd is accurate
nj = 51 # number of cells in j, because this makes the ellipse symmetric
g = 1 # number of ghost cells
xii = 0*np.pi
xif = 2*np.pi
etai = 0.255
etaf = 4.41 #3.72 is 10 chord lengths away, 4.41 is 20
dxi =  (xif  - xii ) / (ni-1)
deta = (etaf - etai) / (nj-1) 
dt = 0.020 #0.0144 # comptime()
 
# Grid
eta = np.zeros((ni+2*g,nj+2*g))
xi = np.zeros((ni+2*g,nj+2*g))
x = np.zeros((ni+2*g,nj+2*g))
y = np.zeros((ni+2*g,nj+2*g))

# Fluid Quantities
###########################
def review(rev=0):
    if rev:
        u0new = 6
        Uload = np.loadtxt('data/U_'+str(int(u0new*100))+'.txt')
        Unew = Uload.reshape(
            Uload.shape[0], Uload.shape[1] // 4, 4)
        U = Unew.copy()
    else:
        U = np.zeros((ni+2*g,nj+2*g,4)) # rho, rhou, rhov, rhoE
    return U
rev = 0
U = review(rev)       
############################
ruvp = np.zeros((ni+2*g,nj+2*g,4)) # rho, u, v and pressure
 
# Derivatives and Jacobian
dxidx  = np.zeros((ni+2*g,nj+2*g))
dxidy  = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))
J      = np.zeros((ni+2*g,nj+2*g))

#### Initializing the Grid
# Notes: ghost cells are treated as 0 index
ib = g # 1
jb = g # 1
ie = ib + ni - 1 # ni
je = jb + nj - 1 # nj
# Initializing the Grid
def ell_grid():
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            xi[i,j]  = (i-ib) * (xif  - xii ) / (ni-1) + xii 
            #xi[i,j] = -1 * xi[i,j]
            eta[i,j] = (j-jb) * (etaf - etai) / (nj-1) + etai
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
            dxidx[i,j] = dy_deta / (dx_dxi*dy_deta-dy_dxi*dx_deta)
            dxidy[i,j] = -dx_deta / (dx_dxi*dy_deta-dy_dxi*dx_deta)
            detadx[i,j] = -dy_dxi / (dx_dxi*dy_deta-dy_dxi*dx_deta)
            detady[i,j] = dx_dxi / (dx_dxi*dy_deta-dy_dxi*dx_deta)
            J[i,j] = dxidx[i,j]*detady[i,j] - detadx[i,j]*dxidy[i,j]
    Bound(J)
    Bound(dxidx)
    Bound(dxidy)
    Bound(detadx)
    Bound(detady)            
#### Initial Condition
# Set Initial Condition
ga = 1.4
r0 = 1.0
p0 = 1.0 / ga
c0 = math.sqrt(p0*ga/r0)
u0 = 6 * c0
#u0 = [0.4, 0.85, 1.5, 6]
v0 = 0.0
def ell_init():
    # All of the other initial conditions
    # Far Field Initial Conditions, only at the outer layer
    for i in np.arange(ib,ie+1):
        for j in np.arange(je+1,je+2): # from je+1 only
            U[i,j,0] = r0  
            U[i,j,1] = r0 * u0 #u0 
            U[i,j,2] = r0 * v0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            E0 = r0q + p0/(ga-1)
            U[i,j,3] = E0
    # Initial Conditions Not including the boundary
    for i in np.arange(ib,ie+1): 
        for j in np.arange(jb+2,je+1):  # from jb+1 to je
            U[i,j,0] = r0
            U[i,j,1] = 0
            U[i,j,2] = 0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            E0 = r0q + p0/(ga-1)
            U[i,j,3] = E0
    # Periodic Boundary Conditions 
    ell_bound()
    # Testing Drag Solver
    get_actual()
    
# Compute time step
''' hmm idk how to implement this one '''
def comp_time():
    get_actual()
    ell_init()
    '''
    T = np.zeros((ni+2*g,2))
    for i in np.arange(ib,ie):
        j = je+1
        c = 1
        lmi = abs(dxidx[i,j] *ruvp[i,j,1]) + math.sqrt(dxidx[i,j]**2+dxidy[i,j]**2) * c
        lmj = abs(detadx[i,j]*ruvp[i,j,1]) + math.sqrt(detadx[i,j]**2+detady[i,j]**2) * c
        #print(lmi)
        #print(lmj)
        T[i,0] = dxi/lmi
        T[i,1] = deta/lmj
    T = T[ib:ie,:].copy()
    pos = np.argmin(T)
    Pos = np.unravel_index(pos, T.shape)
    dt = T[Pos]
    print(dt)
    '''
    
    T = np.zeros((ni+2*g,nj+2*g,2))
    c = np.zeros((ni+2*g,nj+2*g))
    Rangei = np.arange(ib,ie) # ib to ie-1
    Rangej = np.arange(jb+1,je+1+1) # from jb+1 to je
    for i in Rangei: 
        for j in Rangej: 
            c[i,j] = math.sqrt(ga*ruvp[i,j,3]/ruvp[i,j,0])
            #lmi = abs(dxidx[i,j]*u0+dxidy[i,j]*v0) + math.sqrt(dxidx[i,j]**2+dxidy[i,j]**2) * c[i,j]
            #lmj = abs(detadx[i,j]*u0+detady[i,j]*v0) + math.sqrt(detadx[i,j]**2+detady[i,j]**2) * c[i,j]
            lmi = abs(dxidx[i,j]*ruvp[i,j,1]+dxidy[i,j]*ruvp[i,j,2]) + math.sqrt(dxidx[i,j]**2+dxidy[i,j]**2) * c[i,j]
            lmj = abs(detadx[i,j]*ruvp[i,j,1]+detady[i,j]*ruvp[i,j,2]) + math.sqrt(detadx[i,j]**2+detady[i,j]**2) * c[i,j]
            T[i,j,0] = dxi/lmi
            T[i,j,1] = deta/lmj
    T = T[ib:ie,jb+1:je+2,:].copy()
    pos = np.argmin(T)
    Pos = np.unravel_index(pos, T.shape)
    dt = T[Pos]
    print(Pos)
    print(dt)
    '''
    T = np.zeros((ni+2*g,nj+2*g,2))
    Rangei = np.arange(ib,ie) # ib to ie-1
    for i in Rangei: 
        j = je+1
        lmi = abs(dxidx[i,j]*1.5) + math.sqrt(dxidx[i,j]**2+dxidy[i,j]**2)
        lmj = abs(detadx[i,j]*1.5) + math.sqrt(detadx[i,j]**2+detady[i,j]**2) 
        T[i,j,0] = dxi/lmi
        T[i,j,1] = deta/lmj
    T = T[ib:ie,je+1:je+2,:].copy()
    pos = np.argmin(T)
    Pos = np.unravel_index(pos, T.shape)
    dt = T[Pos]
    print(Pos)
    print(dt)
    '''
            
# Function to get actual values
def get_actual():
    for i in np.arange(ib-1,ie+1+1): # goes from ib-1 to ie+1
        for j in np.arange(jb,je+1+1): # goes from jb to je+1
            ruvp[i,j,0] = U[i,j,0]
            ruvp[i,j,1] = U[i,j,1] / U[i,j,0]
            ruvp[i,j,2] = U[i,j,2] / U[i,j,0]
            rq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / U[i,j,0]
            ruvp[i,j,3] = (U[i,j,3] - rq)*(ga-1)

#### Boundary Conditions
def Bound(arr):
    arr[ :, 0] = arr[ :,jb]
    arr[ :,-1] = arr[ :,je]
    arr[ 0, :] = arr[ib, :]
    arr[-1, :] = arr[ie, :]
    return arr 
def per(arr):
    arr[ie,:] = arr[ib,:]   
    arr[ie+1,:] = arr[ib+1,:]   
    arr[0,:] = arr[ie-1,:]  
    return arr
def ell_bound():
    ell_per()
    ell_neu()
    get_actual()
    ell_nos()
def ell_per():
    U[ie  ,:,:] = U[ib  ,:,:]   
    U[ie+1,:,:] = U[ib+1,:,:]
    U[0   ,:,:] = U[ie-1,:,:]   
def ell_neu():
    U[:,jb,0] = U[:,jb+1,0]  
    U[:,jb,3] = U[:,jb+1,3] 
def ell_nos():
    get_actual()
    # Solve for No Slip Boundary Conditions
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,jb+1): # only on jb
            C = ell_solve(i,j)
            U[i,j,1] = C[0]
            U[i,j,2] = C[1] 
    get_actual()
    ell_per()
    get_actual()
def ell_solve(i,j):
    # solves for the u and v values at the ellipse using the boundary conditions 
    A = np.zeros((2,2))
    B = np.zeros(2)
    ubar = dxidx[i,j+1] * ruvp[i,j+1,1] + dxidy[i,j+1] * ruvp[i,j+1,2]
    B[0] = ubar
    B[1] = 0 
    A[0,0] = dxidx[i,j] 
    A[0,1] = dxidy[i,j] 
    A[1,0] = detadx[i,j]
    A[1,1] = detady[i,j]
    C = la.solve(A,B)
    C = C * ruvp[i,j,0]
    return C
    
#### Functions for Flux-Vector splitting
# Compute Eigenvalues
def computelambdas(i,j,k1,k2):
    r = ruvp[i,j,0]
    u = ruvp[i,j,1] 
    v = ruvp[i,j,2]
    p = ruvp[i,j,3]
    #print('r is ', str(r), 'at ', str(i), ' ', str(j))
    #print('p is ', str(p))
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
def FGpm(i,j,FG,PM):
    if FG == 'F':
        k1 = dxidx[i,j] #/ J[i,j]
        #k1 = abs(k1)
        k2 = dxidy[i,j] #/ J[i,j]
        #k2 = abs(k2)
    elif FG == 'G':
        k1 = detadx[i,j] #/ J[i,j]
        #k1 = abs(k1)
        k2 = detady[i,j] #/ J[i,j]
        #k2 = abs(k2)
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
itei = 0
resi = 1
def solve(tf=100.0, maxite=1000, tol=1e-6):
    print('u0 is', str(np.round(u0 ,5)))
    print('dt is', str(np.round(dt ,5)))
    t = ti
    ite = itei
    res = resi
    get_actual()
    flag = False
    while t <= tf and ite < maxite and res > tol and flag==False: # and res > tol
        t += dt
        ite += 1
        ell_bound()
        rhoold = U[:,:,0].copy()
        rhs = np.zeros((ni+2*g,nj+2*g,4))
        
        Rangei = np.arange(ib,ie) # ib to ie-1
        Rangej = np.arange(jb+1,je+1) # from jb+1 to je
        for i in Rangei: 
            for j in Rangej: 
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
        #print(rhs[ib,jb+1:jb+4])
        for i in Rangei: 
            for j in Rangej: 
                U[i,j] -= rhs[i,j] * dt * J[i,j]
                
        # Compute Residual
        res = 0
        for i in Rangei:
            for j in Rangej:
                res += (U[i,j,0]-rhoold[i,j])**2
        res = math.sqrt(res / ( len(Rangei)*len(Rangej) ) )
        ell_per()
        #print(ruvp[ib,jb+1:jb+4,:])
        get_actual()
        #print(ruvp[ib,jb+1:jb+4,:])
        #print(ite)
        if np.any(ruvp[ib:ie,jb+1:je+1,0] < 0) or np.any(ruvp[ib:ie,jb+1:je+1,3] < 0):
            flag = True
            print('Failed lmao')
            A = np.where(ruvp[ib:ie,jb+1:je+1,0] < 0)
            B = np.where(ruvp[ib:ie,jb+1:je+1,3] < 0)
            print(A)
            print(B)
        
    ell_bound()
    #ell_per()
    get_actual()
    
    print('Final time at', str(np.round(t,5)))
    print('Final iteration at', str(ite))
    if res != 1:
        print('Final Residual at', str(np.round(res,5)))

#### Post-processing, plotting contours
def plot_stru():
    fig, ax = plt.subplots()
    # Draws Vertical Lines
    for i in np.arange(ib,ie+1): 
        ax.plot(xi[i][jb:je+1], eta[i][jb:je+1], c='blue') 
    # Draws Horizontal Lines
    #for j in np.arange(jb,je+1): 
        #ax.plot(xi.T[j][ib:ie+1], eta.T[j][ib:ie+1], c='blue')
    ax.plot(xi[4,3],eta[4,3],'r*')
    ax.set_title('Structured Grid')
    plt.show()
def plot_unst(lim=None):
    fig, ax = plt.subplots()
    # Draws Radial Lines
    for i in np.arange(ib,ie+1): 
        ax.plot(x[i][jb:je+1], y[i][jb:je+1], c='blue') 
    # Draws Circular Lines
    for j in np.arange(jb,je+1): 
        ax.plot(x.T[j][ib:ie+1], y.T[j][ib:ie+1], c='blue')
    ax.set_title('Unstructured Grid')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    plt.show()    
def get_umag():
    get_actual()
    umag  = np.zeros((ni+2*g,nj+2*g))
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            um = ruvp[i,j,1]**2 + ruvp[i,j,2]**2
            umag[i,j] = math.sqrt(um)
    return umag   
# Plotting the contours of Jacobian of the grid transformation
def plot_jaco(lim=None):
    fig, ax = plt.subplots()
    ax.set_title('Jacobian Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    #levels = [-1.5,-0.5,0.5,1.5]
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], J[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    plt.show() 

# Plotting the contours of the Velocity Magnitude
def plot_umag(lim=None,save=False,levels=False):
    umag = get_umag()
    #print(ruvp[ib,:,1])
    #print(umag[ib,:])
    fig, ax = plt.subplots()
    ax.set_title('Velocity Magnitude Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    #levels = [-1.5,-0.5,0.5,1.5]
    if levels == True:
        l = [0,0.15,0.3,0.7,0.85,1,3,5,10]
        cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1], levels=l)
    else:
        cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1])
    cbar = fig.colorbar(cf)
    if save == True:
        if lim == None:
            plt.savefig('images/umag_'+str(int(u0*10))+'.png')
        else:
            plt.savefig('images/umag'+str(int(lim*100))+'_'+str(int(u0*10))+'.png')
    plt.show() 
def plot_vel(lim=None,quant='u'):
    umag = get_umag()
    fig, ax = plt.subplots()
    ax.set_title('Velocity Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    #levels = [-1.5,-0.5,0.5,1.5]
    if quant == 'u':
        uv = 1
    elif quant == 'v':
        uv = 2
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], ruvp[ib:ie+1,jb:je+1,uv])
    cbar = fig.colorbar(cf)
    plt.show() 
def plot_density(lim=None):
    get_actual()
    fig, ax = plt.subplots()
    ax.set_title('Density Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], ruvp[ib:ie+1,jb:je+1,0])
    cbar = fig.colorbar(cf)
    plt.show() 
def plot_pressure(lim=None):
    get_actual()
    fig, ax = plt.subplots()
    ax.set_title('pressure vs x')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], ruvp[ib:ie+1,jb:je+1,3])
    cbar = fig.colorbar(cf)
    plt.show()
def plot_flow(lim=None):
    get_actual()
    fig, ax = plt.subplots()
    ax.set_title('Velocity Flow')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
    ax.quiver(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], ruvp[ib:ie+1,jb:je+1,1], ruvp[ib:ie+1,jb:je+1,2], scale=7,scale_units='inches')
    # Draws Radial Lines
    for i in np.arange(ib,ie+1): 
        ax.plot(x[i][jb:je+1], y[i][jb:je+1], c='blue',alpha=0.4) 
    # Draws Circular Lines
    for j in np.arange(jb,je+1): 
        ax.plot(x.T[j][ib:ie+1], y.T[j][ib:ie+1], c='blue',alpha=0.4)
    plt.show()

##### Compute Drag 
def compute_area():
    A = 0
    X =  np.zeros(ni)
    Y =  np.zeros(ni)
    dl =  np.zeros(ni)
    for i in np.arange(ib,ie):
        for j in np.arange(jb,jb+1):
            X[i-ib] = 0.5 * (x[i+1,j] + x[i,j])
            Y[i-ib] = 0.5 * (y[i+1,j] + y[i,j])
    for i in np.arange(ib,ie): 
        for j in np.arange(jb,jb+1): # on jb only
            dx = x[i+1,j] - x[i,j]
            dy = y[i+1,j] - y[i,j]
            dl[i] = math.sqrt(dx**2+dy**2)
            A += dl[i]
    return A, X, Y, dl
def compute_cdrag():
    A,X,Y,dl = compute_area()
    print('Area is', str(np.round(A,5)))
    #idk how to compute the xcomp of pressure
    fd = 0
    U =  np.zeros(ni)
    V =  np.zeros(ni)
    xdir = [1,0]
    ydir = [0,1]
    for i in np.arange(ib,ie): 
        for j in np.arange(jb,jb+1): # on jb only
            p = 0.5 * (ruvp[i,j,3]+ruvp[i+1,j,3])
            normal = [0,0]
            normal[0] = +1 * (y[i+1,j] - y[i,j])
            normal[1] = -1 * (x[i+1,j] - x[i,j])
            normal = normal / la.norm(normal)
            px = p * np.dot(normal,xdir)
            py = p * np.dot(normal,ydir)
            fd += px * dl[i]
            U[i-ib] = px
            V[i-ib] = py
            #print(px)
    #print(fd)
    #print(0.5 * r0 * u0**2 * A)
    cd = fd / (0.5 * r0 * u0**2 * A)
    print('Drag Coefficient is', str(np.round(cd,5)))
    # Draw a quiver just for fun
    # Just for fun lmao
    fig, ax = plt.subplots()
    # Draws Radial Lines
    for i in np.arange(ib,ie+1):
        ax.plot(x[i][jb:je+1], y[i][jb:je+1], c='blue',alpha=0.4) 
    # Draws Circular Lines
    for j in np.arange(jb,je+1): 
        ax.plot(x.T[j][ib:ie+1], y.T[j][ib:ie+1], c='blue',alpha=0.4)
    ax.set_xlim([-1.5,1.5])
    ax.set_ylim([-1.5,1.5])
    ax.set_title('Pressure Acting on the Ellipse')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    oldx =  np.zeros(ni)
    oldy =  np.zeros(ni)
    for i in range(ni):
        oldx[i] = X[i]-0.025*U[i]*ruvp[i+ib,jb,3]/10
        oldy[i] = Y[i]-0.025*V[i]*ruvp[i+ib,jb,3]/10
        ax.annotate("", xytext=(oldx[i],oldy[i]),xy=(X[i], Y[i]),arrowprops=dict(arrowstyle="->"))
    return cd 

##############################################
def main():
    
    ell_grid() 
    deri_jaco() 
    if rev == 0:
        ell_init() 
    
    comp_time() # broken
    if rev == 0:
        pass
        #solve(1000.0, 10, 1e-6) # comment at review
    
    ###### Plotting Stuff ####
    #plot_stru() 
    #plot_unst(1.5) 
    plot_jaco(1.25) 
    plot_vel(1.25,'u')
    plot_vel(1.25,'v')
    #plot_umag()
    plot_umag(1.25)
    #plot_umag(levels=False)
    #plot_umag(lim=5,levels=False)
    #plot_flow(1.5)
    #plot_pressure(1.5)
    ####### Computations
    #compute_cdrag()
    save = 0
    if rev == 0 and save == 1:
        U_reshaped = U.reshape(U.shape[0], -1)
        np.savetxt('data/U_'+str(int(u0*100))+'v2.txt', U_reshaped)
    
if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %2s seconds ---" % np.round((time.time() - start_time),4))
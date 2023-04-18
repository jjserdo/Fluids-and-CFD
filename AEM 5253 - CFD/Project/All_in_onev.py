import os
os.chdir('C:/Users/jjser/OneDrive - Syracuse University/Fall 2022/AEM 5253/Project')

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la
from numpy import asarray
import pandas
import time

def xiieta_grid(ni,nj,xiii,xiif,etai,etaf):
    ## ## Creating xi and eta ## ##
    xii = np.zeros((ni+2,nj+2))
    eta = np.zeros((ni+2,nj+2))
    dxii = (xiif - xiii) / (ni-1)
    deta = (etaf - etai) / (nj-1)
    xiii2 = xiii - dxii
    etai2 = etai - deta
    for i in np.arange(ni+2): 
        for j in np.arange(nj+2): 
            xii[i,j] = (i) * (xiif - xiii) / (ni-1) + xiii2 
            eta[i,j] = (j) * (etaf - etai) / (nj-1) + etai2
    return xii,eta
def sodsod_grid(xii,eta):
    Ni = xii.shape[0]
    Nj = xii.shape[1]
    x   = np.zeros((Ni,Nj))
    y   = np.zeros((Ni,Nj))
    ## ## Creating x and y ## ##
    for i in np.arange(Ni):
        for j in np.arange(Nj):
            x[i,j] = xii[i,j]
            y[i,j] = eta[i,j]
    return x,y
def ellell_grid(xii,eta):
    ni = xii.shape[0]
    nj = xii.shape[1]
    x   = np.zeros((ni,nj))
    y   = np.zeros((ni,nj))
    ## ## Creating x and y ## ##
    for i in np.arange(ni):
        for j in np.arange(nj):
            x[i,j] = np.cosh(eta[i,j]) * np.cos(xii[i,j])
            y[i,j] = np.sinh(eta[i,j]) * np.sin(xii[i,j]) 
    return x,y
#### Calculating Derivatives and Jacobian
def deri_jaco(xii,eta,x,y):
    ni = xii.shape[0]
    nj = xii.shape[1]
    dxii = xii[1,0]-xii[0,0]
    deta = eta[0,1]-eta[0,0]
    dxiidx = np.zeros((ni,nj))
    dxiidy = np.zeros((ni,nj))
    detadx = np.zeros((ni,nj))
    detady = np.zeros((ni,nj))
    J      = np.zeros((ni,nj))
    for i in np.arange(ni):
        for j in np.arange(nj):
            if i > 0 and i < ni-1:
                # central
                dxdxii = (x[i+1,j]-x[i-1,j]) / (2*dxii)
                dydxii = (y[i+1,j]-y[i-1,j]) / (2*dxii)
            elif i == 0:
                # forward
                dxdxii = (x[i+1,j]-x[i,j]) / (dxii)
                dydxii = (y[i+1,j]-y[i,j]) / (dxii)
            elif i == ni-1:
                # backward
                dxdxii = (x[i,j]-x[i-1,j]) / (dxii)
                dydxii = (y[i,j]-y[i-1,j]) / (dxii)
            if j > 0 and j < nj-1:
                # central
                dxdeta = (x[i,j+1]-x[i,j-1]) / (2*deta)
                dydeta = (y[i,j+1]-y[i,j-1]) / (2*deta)
            elif j == 0:
                # forward
                dxdeta = (x[i,j+1]-x[i,j]) / (deta)
                dydeta = (y[i,j+1]-y[i,j]) / (deta)
            elif j == nj-1:
                # backward
                dxdeta = (x[i,j]-x[i,j-1]) / (deta)
                dydeta = (y[i,j]-y[i,j-1]) / (deta)    
            J[i,j]      = 1 / (dxdxii*dydeta-dydxii*dxdeta)
            dxiidx[i,j] =  dydeta * J[i,j]
            dxiidy[i,j] = -dxdeta * J[i,j]
            detadx[i,j] = -dydxii * J[i,j]
            detady[i,j] =  dxdxii * J[i,j]   
    return dxiidx,dxiidy,detadx,detady,J
def sodsod_init(x,y,ga,r0,rf,p0,pf,u0,v0,cut):
    ni = x.shape[0]
    nj = x.shape[1]  
    dx = x[1,0]-x[0,0]
    U = np.zeros((ni,nj,4))
    for i in np.arange(ni):
        for j in np.arange(nj): 
            U[i,j,1] = r0 * u0
            U[i,j,2] = r0 * v0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            rfq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0         
            E0 = r0q + p0/(ga-1)
            Ef = rfq + pf/(ga-1)
            if i*dx <= cut:
                U[i,j,0] = r0
                U[i,j,3] = E0
            else:
                U[i,j,0] = rf
                U[i,j,3] = Ef
    U = sodsod_bound(U)
    ruvp = get_actual(U,ga)
    return U,ruvp
def sodsod_bound(U):
    U[:, 0,:] = U[:, 1,:]
    U[:,-1,:] = U[:,-2,:]
    return U
def ellell_init(x,y,ga,r0,p0,M,v0,dxiidx,dxiidy,detadx,detady):
    c0 = math.sqrt(p0*ga/r0)
    u0 = M * c0
    ni = x.shape[0]
    nj = x.shape[1]
    U = np.zeros((ni,nj,4))
    # All of the other initial conditions
    # Far Field Initial Conditions, only at the outer layer
    for i in np.arange(ni):
        for j in np.arange(2,nj-1): # from je+1 only
            U[i,j,0] = r0  
            U[i,j,1] = r0 * u0 
            U[i,j,2] = r0 * v0
            r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
            E0 = r0q + p0/(ga-1)
            U[i,j,3] = E0
    # Initial Conditions Not including the boundary
    for i in np.arange(ni): 
        j = 2
        U[i,j,0] = r0
        U[i,j,1] = 0
        U[i,j,2] = 0
        r0q = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / r0
        E0 = r0q + p0/(ga-1)
        U[i,j,3] = E0
    # Periodic Boundary Conditions 
    U = ellell_bound(U,ga,dxiidx,dxiidy,detadx,detady)
    # Testing Drag Solver
    ruvp = get_actual(U,ga)
    return U,ruvp
def get_actual(U, ga):
    ni = U.shape[0]
    nj = U.shape[1] 
    ruvp = U.copy()
    for i in np.arange(ni): # goes from ib-1 to ie+1
        for j in np.arange(nj): # goes from jb to je+1
            ruvp[i,j,0] = U[i,j,0]
            ruvp[i,j,1] = U[i,j,1] / U[i,j,0]
            ruvp[i,j,2] = U[i,j,2] / U[i,j,0]
            rq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / U[i,j,0]
            ruvp[i,j,3] = (U[i,j,3] - rq)*(ga-1)
    return ruvp
def ellell_bound(U,ga,dxiidx,dxiidy,detadx,detady):
    U = ellell_per(U)
    U = ellell_neu(U)
    U = ellell_nos(U,ga,dxiidx,dxiidy,detadx,detady)
    return U
def ellell_per(U):
    U[-1  ,:,:] = U[0  ,:,:]   
    return U
def ellell_neu(U):
    U[:,0,0] = U[:,1,0]  
    U[:,0,3] = U[:,1,3] 
    return U
def ellell_nos(U,ga,dxiidx,dxiidy,detadx,detady):
    ruvp = get_actual(U,ga)
    ni = U.shape[0]
    # Solve for No Slip Boundary Conditions
    for i in np.arange(ni):
        C = ellell_solve(i,2,ruvp,dxiidx,dxiidy,detadx,detady)
        U[i,2,1] = C[0]
        U[i,2,2] = C[1] 
    U = ellell_per(U)
    ruvp = get_actual(U,ga)
def ellell_solve(i,j,ruvp,dxiidx,dxiidy,detadx,detady):
    # solves for the u and v values at the ellipse using the boundary conditions 
    A = np.zeros((2,2))
    B = np.zeros(2)
    ubar = dxiidx[i,j+1] * ruvp[i,j+1,1] + dxiidy[i,j+1] * ruvp[i,j+1,2]
    B[0] = ubar
    B[1] = 0 
    A[0,0] = dxiidx[i,j] 
    A[0,1] = dxiidy[i,j] 
    A[1,0] = detadx[i,j]
    A[1,1] = detady[i,j]
    C = la.solve(A,B)
    C = C * ruvp[i,j,0]
    return C
    
#### Functions for Flux-Vector splitting
# Compute Eigenvalues
def computelambdas(i,j,k1,k2,ruvp,ga):
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
def lpm(i,j,k1,k2,ruvp,ga):
    l = computelambdas(i,j,k1,k2,ruvp,ga)
    r = ruvp[i,j,0]
    p = ruvp[i,j,3]
    c = math.sqrt(ga*p/r)
    ep = 0.01 * c
    lp = np.zeros(4)
    lm = np.zeros(4)
    for i in range(4):
        lp[i] = 0.5 * (l[i] + math.sqrt(l[i]**2 + ep ** 2))
        lm[i] = 0.5 * (l[i] - math.sqrt(l[i]**2 + ep ** 2))
    return lp, lm, l
# F and G plus and minus
def FGpm(i,j,FG,PM,dxiidx,dxiidy,detadx,detady,J,ruvp,ga):
    if FG == 'F':
        k1 = dxiidx[i,j] / J[i,j]
        k2 = dxiidy[i,j] / J[i,j]
    elif FG == 'G':
        k1 = detadx[i,j] / J[i,j]
        k2 = detady[i,j] / J[i,j]
    lp, lm, L = lpm(i,j,k1,k2,ruvp,ga)
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
def solve(U,ruvp,xii,eta,dxiidx,dxiidy,detadx,detady,J,dt,tf,maxite,typee,ga,ti=0,itei=0,resi=1):
    t = ti
    ite = itei
    res = resi
    tol = 1e-6
    U = get_actual(U,ga)
    flag = False
    dxii = xii[1,0]-xii[0,0]
    deta = eta[0,1]-eta[0,0]
    ni = U.shape[0]
    nj = U.shape[1]
    while t <= tf and ite < maxite and res > tol and flag==False: # and res > tol
        t += dt
        ite += 1
        if typee == 'sod':
            U = sodsod_bound(U)
            A = np.arange(1,ni-1)
            B = np.arange(1,nj-1)
        elif typee == 'ell':
            U = ellell_bound(U,ga,dxiidx,dxiidy,detadx,detady)
            A = np.arange(1,ni-1)
            B = np.arange(2,nj-1)
        rhoold = U[:,:,0].copy()
        for i in A: 
            for j in B: 
                # Do Explicit Euler to find update u vector
                rhs = np.zeros(4)
                # Evaluate df/dxi
                rhs += FGpm(i  ,j,'F','p',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / dxii
                rhs -= FGpm(i-1,j,'F','p',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / dxii 
                rhs += FGpm(i+1,j,'F','m',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / dxii
                rhs -= FGpm(i  ,j,'F','m',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / dxii
                # Evaluate dg/deta
                rhs += FGpm(i,j  ,'G','p',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / deta
                rhs -= FGpm(i,j-1,'G','p',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / deta
                rhs += FGpm(i,j+1,'G','m',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / deta
                rhs -= FGpm(i,j  ,'G','m',dxiidx,dxiidy,detadx,detady,J,ruvp,ga) / deta   
                U[i,j] -= rhs * dt * J[i,j]
        # Compute Residual
        res = 0
        for i in np.arange(ni):
            for j in np.arange(nj):
                res += (U[i,j,0]-rhoold[i,j])**2
        res = math.sqrt(res / ( len(A)*len(B) ) )
        ruvp = get_actual(U,ga)
        if np.any(ruvp[:,:,0] < 0) or np.any(ruvp[:,:,3] < 0):
            flag = True
            print('Failed lmao')
            A = np.where(ruvp[:,:,0] < 0)
            B = np.where(ruvp[:,:,3] < 0)
            print(A)
            print(B) 
    if typee == 'sod':
        U = sodsod_bound(U)
    elif typee == 'ell':
        U = ellell_bound(U,ga,dxiidx,dxiidy,detadx,detady)
    ruvp = get_actual(U,ga)
    print('Final time at', str(np.round(t,5)))
    print('Final iteration at', str(ite))
    if res != 1:
        print('Final Residual at', str(np.round(res,5)))
    return U,ruvp,t,ite,res

######################################################### Plotting ############
def plot_stru(xii,eta,fig,ax):
    #fig, ax = plt.subplots()
    ni = xii.shape[0]
    nj = xii.shape[1]
    # Draws Horizontal Lines
    for i in np.arange(1,ni-1): 
        ax.plot(  xii[i,1:-1],   eta[i,1:-1], c='blue') 
    # Draws Vertical Lines
    for j in np.arange(1,nj-1): 
        ax.plot(xii.T[j,1:-1], eta.T[j,1:-1], c='blue')
    ax.set_title('Structured Grid')
    xlabel = "$\\xi $"
    ylabel = "$\\eta $"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig, ax
def plot_unst(x,y,fig,ax,lim=None,alf=1):
    ni = x.shape[0]
    nj = x.shape[1]
    # Draws Radial Lines
    for i in np.arange(1,ni-1): 
        ax.plot(  x[i,1:-1],   y[i,1:-1], c='blue', alpha=alf) 
    # Draws Circular Lines
    for j in np.arange(nj): 
        ax.plot(x.T[j,1:-1], y.T[j,1:-1], c='blue', alpha=alf)
    ax.set_title('Unstructured Grid')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    return fig,ax      
# Plotting the contours of Jacobian of the grid transformation
def plot_jaco(x,y,J,fig,ax,lim=None):
    fig, ax = plt.subplots()
    ax.set_title('Jacobian Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    J = np.round(J,1)
    l = [-1.5,0,1.5]
    cf = ax.contourf(x[1:-1,1:-1], y[1:-1,1:-1], J[1:-1,1:-1],levels=l)
    cbar = fig.colorbar(cf)
    return fig, ax
def plot_sod_quant(x,y,ruvp,quant,fig,ax,lim=None):
    if quant == 0:
        title = 'Density vs x'
    if quant == 1:
        title = 'U-Velocity vs x'
    if quant == 2:
        title = 'V-Density vs x'
    if quant == 3:
        title = 'Pressure vs x'
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    ax.plot(x[1:-1,-2], ruvp[1:-1,-1,quant],'bo')
    return fig, ax
# Plotting the contours of the Velocity Magnitude
def plot_umag(x,y,ruvp,fig,ax,lim=None,save=False,levels=False):
    ni = ruvp.shape[0]
    nj = ruvp.shape[1]
    umag  = np.zeros((ni,nj))
    for i in np.arange(ni): 
        for j in np.arange(nj): 
            umag[i,j] = math.sqrt(ruvp[i,j,1]**2 + ruvp[i,j,2]**2)
    ax.set_title('Velocity Magnitude Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    if levels == True:
        l = [0,0.15,0.3,0.7,0.85,1,3,5,10]
        cf = ax.contourf(x[1:-1,1:-1], y[1:-1,1:-1], umag[1:-1,1:-1], levels=l)
    else:
        cf = ax.contourf(x[0:-1,0:-1,1], y[0:-1,0:-1,1], umag[0:-1,0:-1,1])
        #cf = ax.contourf(x[1:-1,1:-1], y[1:-1,1:-1], umag[1:-1,1:-1])
    cbar = fig.colorbar(cf)
    return fig, ax
def plot_quant(x,y,ruvp,quant,fig,ax,lim=None):
    if quant == 0:
        title = 'Density Contour Plot'
    if quant == 1:
        title = 'U-Velocity Contour Plot'
    if quant == 2:
        title = 'V-Density Contour Plot'
    if quant == 3:
        title = 'Pressure Contour Plot'
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    cf = ax.contourf(x[0:-1,0:-1,1], y[0:-1,0:-1,1], ruvp[0:-1,0:-1,quant])
    cbar = fig.colorbar(cf)
    return fig, ax
def plot_flow(x,y,ruvp,fig,ax,lim=None):
    if lim != None:
        ax.set_xlim([lim[0],lim[1]])
        ax.set_ylim([lim[2],lim[3]])
    ax.quiver(x[0:-1,0:-1], y[0:-1,0:-1], ruvp[0:-1,0:-1,1], ruvp[0:-1,0:-1,2], scale=7,scale_units='inches')
    fig, ax = plot_unst(x,y,fig,ax,alf=0.4)
    ax.set_title('Velocity Flow')
    return fig,ax
def plot_presell(x,y,pxx,pyy,fig,ax):
    fig,ax = plot_unst(x, y, fig, ax)
    ax.set_title('Pressure Acting on the Ellipse')
    ni = x.shape[0]
    oldx =  np.zeros(ni)
    oldy =  np.zeros(ni)
    A, xmid, ymid, dl = compute_area(x, y)
    for i in range(ni):
        oldx[i] = xmid[i]-0.025*pxx[i]
        oldy[i] = ymid[i]-0.025*pyy[i]
        ax.annotate("", xytext=(oldx[i],oldy[i]),xy=(xmid[i], ymid[i]),arrowprops=dict(arrowstyle="->"))
    return fig,ax 
######################################################## Computations #########
# Compute time step
''' Make it realize between ghosts and not based on array sizes '''
def ellell_comp_time(U,ga,xii,eta,dxiidx,dxiidy,detadx,detady,typee):
    ruvp = get_actual(U,ga)
    ni = dxiidx.shape[0]
    nj = dxiidx.shape[1]
    dxii = xii[1,0]-xii[0,0]
    deta = eta[0,1]-eta[0,0]
    T = np.zeros((ni,nj,2))
    for i in np.arange(ni): 
        for j in np.arange(nj): 
            c = math.sqrt(ga*ruvp[i,j,3]/ruvp[i,j,0])
            lmi = abs(dxiidx[i,j]*ruvp[i,j,1]+dxiidy[i,j]*ruvp[i,j,2]) + math.sqrt(dxiidx[i,j]**2+dxiidy[i,j]**2) * c
            lmj = abs(detadx[i,j]*ruvp[i,j,1]+detady[i,j]*ruvp[i,j,2]) + math.sqrt(detadx[i,j]**2+detady[i,j]**2) * c
            T[i,j,0] = dxii/lmi
            T[i,j,1] = deta/lmj
    pos = np.argmin(T)
    Pos = np.unravel_index(pos, T.shape)
    dt = T[Pos]
    return dt, Pos
def compute_area(x,y):
    A = 0
    ni = x.shape[0]
    xmid = np.zeros(ni)
    ymid = np.zeros(ni)
    dl   = np.zeros(ni)
    for i in np.arange(ni-1):
        xmid[i] = 0.5 * (x[i+1,0] + x[i,0])
        ymid[i] = 0.5 * (y[i+1,0] + y[i,0])
        dx = x[i+1,0] - x[i,0]
        dy = y[i+1,0] - y[i,0]
        dl[i] = math.sqrt(dx**2+dy**2)
        A += dl[i]
    return A, xmid, ymid, dl
def compute_cdrag(x,y,ruvp,r0,u0):
    ni = x.shape[0]
    A, xmid, ymid, dl = compute_area(x, y)
    fd = 0
    pxx,pyy = calc_presell(x, y, ruvp)
    for i in np.arange(ni-1): 
        fd += pxx[i] * dl[i]
    cd = fd / (0.5 * r0 * u0**2 * A)
    return cd 
def calc_presell(x,y,ruvp):
    ni = x.shape[0]
    pxx =  np.zeros(ni)
    pyy =  np.zeros(ni)
    xdir = [1,0]
    ydir = [0,1]
    for i in np.arange(ni-1): 
        p = 0.5 * (ruvp[i,0,3]+ruvp[i+1,0,3])
        normal = [0,0]
        normal[0] = +1 * (y[i+1,0] - y[i,0])
        normal[1] = -1 * (x[i+1,0] - x[i,0])
        normal = normal / la.norm(normal)
        pxx[i] = p * np.dot(normal,xdir)
        pyy[i] = p * np.dot(normal,ydir)
    return pxx,pyy

################################################# Main ########################
def save(u0,U,suffix):
    U_reshaped = U.reshape(U.shape[0], -1)
    np.savetxt('data/U_'+str(int(u0*100))+suffix+'.txt', U_reshaped)
def main(typee):
    if typee == 'sod':
        ########### Initialization #####
        #xii,eta = xiieta_grid(51,2,0,1,0,1) # ni,nj,xii,xif,etai,etaf
        xii,eta = xiieta_grid(51,2,0,4.5,0,1)
        x,y     = sodsod_grid(xii,eta)
        dxiidx,dxiidy,detadx,detady,J = deri_jaco(xii,eta,x,y)
        #print(J.shape)
        #print(dxiidx.shape)
        ########### Solving ############

        U,ruvp = sodsod_init(x,y,1.4,1,0.125,1,0.1,0,0,1.95)
        U,ruvp,t,ite,res = solve(U,ruvp,xii,eta,dxiidx,dxiidy,detadx,detady,J,0.008,1,1000,'sod',1.4,ti=0,itei=0,resi=1)
        #U,ruvp = solve(U,ruvp)

        #U,ruvp = sodsod_init(x,y,1.4,1,0.125,1,0.1,0,0,0.5) #(x,y,ga,r0,rf,p0,pf,u0,v0,cut):
        ########### Plots ###############
        fig1, ax1 = plt.subplots()
        plot_stru(xii,eta,fig1,ax1)
        fig1.savefig('images/sod_struv2.png')
        
        fig2, ax2 = plt.subplots()
        plot_unst(x,y,fig2,ax2)
        fig2.savefig('images/sod_unstv2.png')
        
        fig3, ax3 = plt.subplots()
        plot_jaco(x,y,J,fig3,ax3)
        fig3.savefig('images/sod_jaco_45v2.png')
        
        fig4, ax4 = plt.subplots()
        plot_sod_quant(x,y,ruvp,0,fig4,ax4)
        fig4.savefig('images/sod_rho_45v2.png')
        
        fig5, ax5 = plt.subplots()
        plot_sod_quant(x,y,ruvp,3   ,fig5,ax5)
        fig5.savefig('images/sod_pres_45v2.png')
        
        fig6, ax6 = plt.subplots()
        plot_umag(x,y,ruvp,fig6,ax6)
        fig6.savefig('images/sod_umag_45v2.png')
        
    elif typee == 'ell':
        xii,eta = xiieta_grid(61,50,0,2*np.pi,0.255,4.41)
        x,y     = ellell_grid(xii,eta)
        dxiidx,dxiidy,detadx,detady,J = deri_jaco(xii,eta,x,y)
        
        U,ruvp = ellell_init(x,y,1.4,1,1/1.4,0.4,0,dxiidx,dxiidy,detadx,detady)
        #U,ruvp,t,ite,res = solve(U,ruvp,xii,eta,dxiidx,dxiidy,detadx,detady,J,0.008,1,1000,'sod',1.4,ti=0,itei=0,resi=1)
        
        fig1, ax1 = plt.subplots()
        plot_unst(x,y,fig1,ax1)
        fig1.savefig('images/ell_unstv2.png')
        
        fig2, ax2 = plt.subplots()
        plot_unst(x,y,fig2,ax2,[-1.25,1.25,-1.25,1.25])
        fig2.savefig('images/ell_unst2v2.png')
        
        fig3, ax3 = plt.subplots()
        plot_jaco(x,y,J,fig3,ax3,[-1.25,1.25,-1.25,1.25])
        fig3.savefig('images/ell_jacov2.png')
        
        fig4, ax4 = plt.subplots()
        plot_flow(x,y,ruvp,fig4,ax4,[-1.25,1.25,-1.25,1.25])
        fig4.savefig('images/ell_flowv2.png')
        
        fig5, ax5 = plt.subplots()
        plot_umag(x,y,ruvp,fig6,ax6)
        fig5.savefig('images/U_402.png')
        
        fig6, ax6 = plt.subplots()
        plot_umag(x,y,ruvp,fig6,ax6,[-1.25,1.25,-1.25,1.25])
        fig6.savefig('images/U_40v2.png')
        
        fig7, ax7 = plt.subplots()
        plot_umag(x,y,ruvp,fig7,ax7,[-1.25,1.25,-1.25,1.25])
        fig7.savefig('images/ell_presell.png')
        
    return U
def load(filename):
    Uload = np.loadtxt('data/'+filename+'.txt')
    Unew  = Uload.reshape(Uload.shape[0], Uload.shape[1] // 4, 4)
    U = Unew.copy()

################################################## On Run #####################
if __name__ == "__main__":
    start_time = time.time()
    #U = load('U_150v3')
    main('ell')
    #U = main('ell')
    #save(1.5,U,v3)
    print("--- %2s seconds ---" % np.round((time.time() - start_time),4))
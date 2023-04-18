import os
os.chdir('C:/Users/jjser/OneDrive - Syracuse University/Fall 2022/AEM 5253/Project')

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la
from numpy import asarray
import pandas
import time

ga = 1.4
r0 = 1
u0 = 1.5
Uload = np.loadtxt('data/U_'+str(int(u0*100))+'v2.txt')
Unew = Uload.reshape(
    Uload.shape[0], Uload.shape[1] // 4, 4)
U = Unew.copy()
g = 1
ni = U.shape[0]-2*g
nj = U.shape[0]-2*g
xii = 0*np.pi
xif = 2*np.pi
etai = 0.255
etaf = 4.41 
dxi =  (xif  - xii ) / (ni-1)
deta = (etaf - etai) / (nj-1) 
eta = np.zeros((ni+2*g,nj+2*g))
xi = np.zeros((ni+2*g,nj+2*g))
x = np.zeros((ni+2*g,nj+2*g))
y = np.zeros((ni+2*g,nj+2*g))
ruvp = np.zeros((ni+2*g,nj+2*g,4)) 
dxidx  = np.zeros((ni+2*g,nj+2*g))
dxidy  = np.zeros((ni+2*g,nj+2*g))
detadx = np.zeros((ni+2*g,nj+2*g))
detady = np.zeros((ni+2*g,nj+2*g))
J      = np.zeros((ni+2*g,nj+2*g))
ib = g # 1
jb = g # 1
ie = ib + ni - 1 # ni
je = jb + nj - 1 # nj
def ell_grid():
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            xi[i,j]  = (i-ib) * (xif  - xii ) / (ni-1) + xii 
            xi[i,j] = -1 * xi[i,j]
            eta[i,j] = (j-jb) * (etaf - etai) / (nj-1) + etai
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,je+1):
            x[i,j] = np.cosh(eta[i,j]) * np.cos(xi[i,j])
            y[i,j] = np.sinh(eta[i,j]) * np.sin(xi[i,j])        
def deri_jaco():
    for i in np.arange(ib,ie+1):
        for j in np.arange(jb,je+1):
            if i > ib and i < ie:
                dx_dxi = (x[i+1,j]-x[i-1,j]) / (2 * dxi)
                dy_dxi = (y[i+1,j]-y[i-1,j]) / (2 * dxi)
            elif i == ib:
                dx_dxi = (x[i+1,j]-x[i,j]) / (dxi)
                dy_dxi = (y[i+1,j]-y[i,j]) / (dxi)
            elif i == ie:
                dx_dxi = (x[i,j]-x[i-1,j]) / (dxi)
                dy_dxi = (y[i,j]-y[i-1,j]) / (dxi)
            if j > jb and j < je:
                dx_deta = (x[i,j+1]-x[i,j-1]) / (2 * deta)
                dy_deta = (y[i,j+1]-y[i,j-1]) / (2 * deta)
            elif j == jb:
                dx_deta = (x[i,j+1]-x[i,j]) / (deta)
                dy_deta = (y[i,j+1]-y[i,j]) / (deta)
            elif j == je:
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
def get_actual():
    for i in np.arange(ib-1,ie+1+1): # goes from ib-1 to ie+1
        for j in np.arange(jb,je+1+1): # goes from jb to je+1
            ruvp[i,j,0] = U[i,j,0]
            ruvp[i,j,1] = U[i,j,1] / U[i,j,0]
            ruvp[i,j,2] = U[i,j,2] / U[i,j,0]
            rq = 0.5 * (U[i,j,1]**2 + U[i,j,2]**2) / U[i,j,0]
            ruvp[i,j,3] = (U[i,j,3] - rq)*(ga-1)
def Bound(arr):
    arr[ :, 0] = arr[ :,jb]
    arr[ :,-1] = arr[ :,je]
    arr[ 0, :] = arr[ib, :]
    arr[-1, :] = arr[ie, :]
    return arr 
def plot_stru():
    fig, ax = plt.subplots()
    # Draws Horizontal Lines
    for i in np.arange(ib,ie+1): 
        ax.plot(xi[i][jb:je+1], eta[i][jb:je+1], c='blue') 
    # Draws Vertical Lines
    for j in np.arange(jb,je+1): 
        ax.plot(xi.T[j][ib:ie+1], eta.T[j][ib:ie+1], c='blue')
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
        #l = [0,0.15,0.3,0.7,0.85,1,3,5,10]
        c = ('r', 'orange', 'y', 'g', 'b')
        l = [0,0.15,0.3,0.85,1,3]
        cf = ax.contourf(x[ib:ie+1,jb:je+1], y[ib:ie+1,jb:je+1], umag[ib:ie+1,jb:je+1], levels=l,colors=c) #, cmap='gray') #
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
    ########### 
    #plot_stru() 
    #plot_unst(1.5) 
    plot_jaco(1.25) 
    #plot_vel(5,'u')
    #plot_vel(5,'v')
    #plot_umag()
    #plot_umag(5)
    #plot_umag(levels=False)
    plot_umag(lim=5,levels=True,save=False)
    plot_flow(1.5)
    #plot_pressure(1.5)
    #compute_cdrag()
      
if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %2s seconds ---" % np.round((time.time() - start_time),4))
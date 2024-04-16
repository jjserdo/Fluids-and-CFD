# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 19:35:19 2022

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt

# Initialize values
L = 1
a = 1
N = 50
Tp = L/a
Tf = 10 * Tp
dx = L / N
CFL = 0.5
dt = CFL * dx / a

# Create x array
xarray0 = - 2 * dx
xarray1 = L + dx
x = np.linspace(xarray0, xarray1, N+4)

print(x.shape)
print(x[0:3])
print(x[-1])
print(x[-2])

print(dt)
T = Tf / dt
t = np.linspace(0, Tf, int(T)+1)
Tc = int(len(t)/10)
print(t[0])
print(t[Tc])
print(t[-1])

''' Initial Condition Functions '''
def ghost(f):
    f[1]   = f[N+1]
    f[0]   = f[N]
    f[N+2] = f[2]
    f[N+3] = f[3]
    return f
    
def f0(x, shape):
    if shape == 0:
        f0 = np.exp(-100*(x-0.5)**2)
    if shape == 1:
        f0 = np.zeros(len(x))
        for i in range(len(f0)):
            if x[i] > 0.25 and x[i] <= 0.75:
                f0[i] = 1
    f0 = ghost(f0)
    return f0

''' Right Hand Side Functions '''
def s1(i, j):
    dfdx = (f[i][j] - f[i][j-1]) / (dx)
    return dfdx

def s2(i, j):
    dfdx = (f[i][j+1] - f[i][j-1]) / (2*dx)
    return dfdx

def s3(i, j):
    dfdx = (f[i][j-2]-6*f[i][j-1]+3*f[i][j]+2*f[i][j+1]) / (12*dx)
    return dfdx

def s4(i, j):
    dfdx = (-f[i][j+2]+8*f[i][j+1]-8*f[i][j-1]+f[i][j-2]) / (12*dx)
    return dfdx

''' Left Hand Side Functions '''
def exp(i, j, scheme):
    if scheme == 1:
        df = -a * s1(i, j)
    if scheme == 2:
        df = -a * s2(i, j)
    if scheme == 3:
        df = -a * s3(i, j)
    if scheme == 4:
        df = -a * s4(i, j)
        
    fi_1 = f[i][j] + dt * df
    return fi_1

def rk3(i, j, scheme):
    if scheme == 1:
        k0 = s1(i,j)
        k1 = s1(i,j) + dt/2 * k0
        k2 = s1(i,j) + 2 * dt * k1 - dt*k0
    if scheme == 2:
        k0 = s2(i,j)
        k1 = s2(i,j) + dt/2 * k0
        k2 = s2(i,j) + 2 * dt * k1 - dt*k0
    if scheme == 3:
        k0 = s3(i,j)
        k1 = s3(i,j) + dt/2 * k0
        k2 = s3(i,j) + 2 * dt * k1 - dt*k0
    if scheme == 4:
        k0 = s4(i,j)
        k1 = s4(i,j) + dt/2 * k0
        k2 = s4(i,j) + 2 * dt * k1 - dt*k0
        
    df  = -a * (k0 + 4*k1 + k2) / 6
    fi_1 = f[i][j] + dt * df
    return fi_1

''' Basic Code '''
f = np.zeros((len(t),len(x)))

''' Gaussian '''
# Initialize Function
count = 1
for shape in np.arange(2):
    f[0] = f0(x, shape)
    for k in np.arange(4)+1:
        # Calculate values of f for all time steps
        for i in np.arange(len(t)-1):
            # Calculate values of f for all distance
            for j in np.arange(N)+2:
                f[i+1][j] = rk3(i, j, k)
            # Update boundary conditions
            f[i+1] = ghost(f[i+1])
        '''
        # Plot f in different times
        num = 10
        start = 200
        for i in np.arange(num)+start:
            plt.plot(x[2:N+2],f[i][2:N+2])  
        '''
        num = 3 # max 6 at 10 Tp
        for n in np.arange(num):
            label = 't = ' + str(2*n) + 'Tp'
            plt.plot(x[2:N+2],f[2*Tc*n][2:N+2],label=label)
            
        if shape == 0:
            tshape = 'Gaussian'
        if shape == 1:
            tshape = 'Square'
        title = 'Advection for ' + tshape + ' with Rk3 and S' + str(k) + ' scheme'
        plt.title(title)
        plt.xlabel('x')
        plt.ylabel('f')
        plt.legend()
        
        figname = 'images/ad' + str(count) + '.jpg'
        #plt.savefig(figname)
        plt.show()
        count += 1


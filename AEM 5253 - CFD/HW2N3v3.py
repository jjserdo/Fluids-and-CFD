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

t = 1
dx = L / N
xarray0 = - 2 * dx
xarray1 = L + dx
x = np.linspace(xarray0, xarray1, N+4)
CFL = 0.5
dt = CFL * dx / a
T = Tf / dt
t = np.linspace(0, Tf, int(T)+1)
Tc = int(len(t)/10)

''' Functions '''
def f0_gauss(x):
    f0 = np.exp(-100*(x-0.5)**2)
    return f0

def f0_square(x):
    f0 = np.zeros(len(x))
    for i in range(len(f0)):
        if x[i] > 0.25 and x[i] <= 0.75:
            f0[i] = 1
    return f0

def g1(i, j):
    g = -a * (f[i][j+1] - f[i][j-1]) / (2*dx)
    return g

def g2(x):
    g = - a * (100 - 200 * x) * np.exp( -100*(x-0.5)**2 ) 
    return g

def rk3(i, j):
    k0 = g2(x[j])
    k1 = g2(x[j] + dt/2 * k0)
    k2 = g2(x[j] + dt/2 * k1 - dt*k0)
    fi_1 = f[i][j] + dt/6 * (k0 + 4*k1 + k2)
    return fi_1

def exp1(i, j):
    df = g1(i, j)
    fi_1 = f[i][j] + dt * df
    return fi_1

def exp2(i, j):
    df = g2(x[j])
    fi_1 = f[i][j] + dt * df
    return fi_1

''' Code for Gaussian '''
f = np.zeros((len(t),len(x))) # (81,8)
f[0] = f0_gauss(x)
f[0][0]   = f[0][N-1]
f[0][1]   = f[0][N]
f[0][N+2] = f[0][2]
f[0][N+3] = f[0][3]

label = 't = 0'
plt.plot(x[2:N+2],f[0][2:N+2],label=label)

for i in np.arange(len(t)-1):
    for j in np.arange(N+3):
        f[i+1][j] = exp1(i, j)
    f[i+1][2]   = f[i+1][N+2]
    f[i+1][1]   = f[i+1][N+1]
    f[i+1][0]   = f[i+1][N]
    f[0][N+2] = f[i+1][2]
    f[0][N+3] = f[i+1][3]
      
for i in np.arange(10)+1:
    plt.plot(x[2:N+2],f[i][2:N+2])  
plt.title('testing')
plt.legend()
plt.show()

label = 't = 0'
plt.plot(x[2:N+2],f[0][2:N+2],label=label)

for n in np.arange(1)+1:
    label = 't = ' + str(2*n) + 'Tp'
    plt.plot(x[2:N+2],f[2*Tc*n][2:N+2],label=label)

plt.title('testing')
plt.legend()
plt.show()

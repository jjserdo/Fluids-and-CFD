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

def g(i, j):
    g = -a * (f[i][j+1] - f[i][j-1]) / (2*dx)
    return g

def u11(i,j):
    u1 = f[i][j] - a / 2 * (f[i][j+1]-f[i][j-1])
    return u1

def u21(i,j):
    u2 = 0.75*f[i][j] + 0.25*u11(i,j) - a/8*(u11(i,j+1)-u11(i,j-1))
    return u2

def rk3(i, j):
    fi_1 = 1/3*f[i][j] + u21(i,j) -a/3 * (u21(i,j+1)-u21(i,j-1))
    return fi_1
def u12(i,j):
    u1 = f[i][j] + dt * g(i,j)
    return u1
def u21(i,j):
    u2 = 0.75*f[i][j] + 0.25*u12(i,j) + 0.25*dt*g(i,j)
    return u2
def rk3(i, j):
    fi_1 = 1/3*f[i][j]+u21(i,j)-a/3*(u21(i,j+1)-u21(i,j-1))
    return fi_1


'''
def rk3(fi,ti):
    k0 = g(i, j)
    k1 = g(fi+dt/2*k0,ti+dt/2)
    k2 = g(fi+dt/2*k1,ti+dt/2)
    k3 = g(fi+dt  *k2,ti+dt  )
    fi_1 = fi + dt/6 * (k0+2*k1+2*k2+k3)
    return fi_1
'''

def exp(i, j):
    df = g(i, j)
    fi_1 = f[i][j] + df * dt
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
    for j in np.arange(N)+2:
        f[i+1][j] = exp(i, j)
        
for i in np.arange(30)+1:
    plt.plot(x[2:N+2],f[i][2:N+2])        
'''
for n in np.arange(6)+1:
    label = 't = ' + str(2*n) + 'Tp'
    plt.plot(x[2:N+2],f[2*Tc*n][2:N+2],label=label)
'''
plt.title('testing')
plt.legend()
plt.show()

''' Code for Square '''
f = np.zeros((len(t),len(x))) # (81,8)
f[0] = f0_square(x)
f[0][0]   = f[0][N-1]
f[0][1]   = f[0][N]
f[0][N+2] = f[0][2]
f[0][N+3] = f[0][3]

label = 't = 0'
plt.plot(x[2:N+2],f[0][2:N+2],label=label)

for i in np.arange(len(t)-1):
    for j in np.arange(N)+2:
        f[i+1][j] = exp(i, j)
        
for i in np.arange(10)+1:
    plt.plot(x[2:N+2],f[i][2:N+2])        
'''
for n in np.arange(6)+1:
    label = 't = ' + str(2*n) + 'Tp'
    plt.plot(x[2:N+2],f[2*Tc*n][2:N+2],label=label)
'''
plt.title('testing')
plt.legend()
plt.show()

''' Code for Gaussian rk3 '''
f = np.zeros((len(t),len(x))) # (81,8)
f[0] = f0_gauss(x)
f[0][0]   = f[0][N-1]
f[0][1]   = f[0][N]
f[0][N+2] = f[0][2]
f[0][N+3] = f[0][3]

label = 't = 0'
plt.plot(x[2:N+2],f[0][2:N+2],label=label)

for i in np.arange(len(t)-1):
    for j in np.arange(N)+1:
        f[i+1][j] = rk3(i, j)
        
for i in np.arange(10)+1:
    plt.plot(x[2:N+2],f[i][2:N+2])        
'''
for n in np.arange(6)+1:
    label = 't = ' + str(2*n) + 'Tp'
    plt.plot(x[2:N+2],f[2*Tc*n][2:N+2],label=label)
'''
plt.title('testing')
plt.legend()
plt.show()




'''

f1 = np.zeros(len(f))
print(range(4))
for i in range(len(f1)):
    f1[i+1] = rk3()

def g(f,i):
    g = -a * (f[i+1]-f[i-1]) / (2 * dx)
    return g

def rk4(ti,tf,yi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    y[0] = yi
    for i in range(len(t)-1):
        k0 = g(y[i],t[i])
        k1 = g(y[i]+dt/2*k0,t[i]+dt/2)
        k2 = g(y[i]+dt/2*k1,t[i]+dt/2)
        k3 = g(y[i]+dt  *k2,t[i]+dt  )
        y[i+1] = y[i] + dt/6 * (k0+2*k1+2*k2+k3)
    return t, y

def rk3(fi,ti):
    k0 = g(fi,ti)
    k1 = g(fi+dt/2*k0,ti+dt/2)
    k2 = g(fi+dt/2*k1,ti+dt/2)
    k3 = g(fi+dt  *k2,ti+dt  )
    fi_1 = fi + dt/6 * (k0+2*k1+2*k2+k3)
    return fi_1
        
# Advection Initialize Equation Solver
def advec(L,a,N,Tp,Tf):
    
    pass
    return True

# Plots
fig, ax = plt.subplots()

for i in linspace(0,):
    t, y = imp(0,15,4,i)
    label = 'Implicit Euler dt=' + str(i) + 's'
    ax1.plot(t,y,'-o',label=label)
    
#ax.set_title('Implicit Euler Method with Analytical Solution')
#ax.set_xlabel('t')
#ax.set_ylabel('y')
#ax.set_xlim([0,15])
#ax1.set_ylim([0, 4])
#ax.legend()
#plt.savefig('images/imp1a.jpg')
plt.show()
'''
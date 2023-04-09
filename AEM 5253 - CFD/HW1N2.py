"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1 Number 2
"""

import numpy as np
import matplotlib.pyplot as plt

dt = [0.1,0.5,1]

# solvers
def g(y, t):
    return - (2+0.01*t**2) * y;
def anl(ti,tf,yi):
    t = np.linspace(ti,tf,1000)
    y = np.zeros(len(t))
    y[0] = yi
    yi * np.exp(-2*t - 0.01/3*t**3)
    return t, y
def imp(ti,tf,yi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    y[0] = yi
    for i in range(len(t)-1):
        nume = y[i] 
        denu = 1 + dt * (2+0.01*t[i]**2)
        y[i+1] = nume / denu
######
def exp(ti,tf,yi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    y[0] = yi
    for i in range(len(t)-1):
        y[i+1] = y[i] + dt * g(y[i],t[i])
    return t, y
def rk2(ti,tf,yi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    y[0] = yi
    for i in range(len(t)-1):
        k0 = g(y[i],t[i])
        k1 = g(y[i]+dt/2*k0,t[i]+dt/2)
        y[i+1] = y[i] + dt * k1
    return t, y
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

# Analytical
t_anal, y_anal = anl(0,15,4)

# Explicit  
fig, ax = plt.subplots()
ax.plot(t_anal,y_anal,'-r',label='Analytical Solution') 
for i in dt:
    t, y = exp(0,15,4,i)
    label = 'Explicit Euler dt=' + str(i) + 's'
    ax.plot(t,y,'-o',label=label)
ax.set_title('Explicit Euler Method with Analytical Solution')
ax.set_xlabel('t')
ax.set_ylabel('y')
ax.set_xlim([0,15])
ax.legend()
plt.savefig('images/exp1a.jpg')
plt.show()

# Implicit
fig, ax = plt.subplots()
ax.plot(t_anal,y_anal,'-r',label='Analytical Solution') 
for i in dt:
    t, y = imp(0,15,4,i)
    label = 'Implicit Euler dt=' + str(i) + 's'
    ax.plot(t,y,'-o',label=label)
ax.set_title('Implicit Euler Method with Analytical Solution')
ax.set_xlabel('t')
ax.set_ylabel('y')
ax.set_xlim([0,15])
#ax.set_ylim([0, 4])
ax.legend()
plt.savefig('images/imp1a.jpg')
plt.show()

# RK2
fig, ax = plt.subplots()
ax.plot(t_anal,y_anal,'-r',label='Analytical Solution') 
for i in dt:
    t, y = rk2(0,15,4,i)
    label = 'RK2 dt=' + str(i) + 's'
    ax.plot(t,y,'-o',label=label)
ax.set_title('Runge-Kutta 2nd Order with Analytical Solution')
ax.set_xlabel('t')
ax.set_ylabel('y')
ax.set_xlim([0,15])
ax.legend()
plt.savefig('images/rk21a.jpg')
plt.show()

# RK4
fig, ax = plt.subplots()
ax.plot(t_anal,y_anal,'-r',label='Analytical Solution') 
for i in dt:
    t, y = rk4(0,15,4,i)
    label = 'RK4 dt=' + str(i) + 's'
    ax.plot(t,y,'-o',label=label)
ax.set_title('Runge-Kutta 4th Order with Analytical Solution')
ax.set_xlabel('t')
ax.set_ylabel('y')
ax.set_xlim([0,15])
ax.legend()
plt.savefig('images/rk41a.jpg')
plt.show()


 


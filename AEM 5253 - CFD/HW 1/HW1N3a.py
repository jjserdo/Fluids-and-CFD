"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1 Number 3 Part B
"""
import numpy as np
import matplotlib.pyplot as plt
import math

dt = [0.15, 0.5, 1]
grav = 9.81
Leng = 0.6

# g(y(t), t)
def g1(theta, dtheta, t):
    return dtheta;
def g2(theta, dtheta, t):
    return -grav/Leng * theta;

# solvers
def anl(ti,tf,yi,vi):
    t = np.linspace(ti,tf,1000)
    y = np.zeros(len(t))
    y[0] = yi
    for i in range(len(t)):    
        y[i] = yi * math.cos(t[i] * math.sqrt(grav/Leng))
    return t, y
def exp(ti,tf,yi,vi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] = yi
    v[0] = vi
    for i in range(len(t)-1):
        y[i+1] = y[i] + dt * g1(y[i], v[i], t[i]);
        v[i+1] = v[i] + dt * g2(y[i], v[i], t[i]);
    return t, y
def imp(ti,tf,yi,vi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] = yi
    v[0] = vi
    for i in range(len(t)-1):
        # version 1; based on coupled equations
        nume = y[i] + v[i] * dt
        denu = 1 + grav/Leng * dt**2
        y[i+1] = nume / denu
        v[i+1] = (y[i+1] - y[i]) / dt;
    return t, y
def rk4(ti,tf,yi,vi,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] = yi
    v[0] = vi
    for i in range(len(t)-1):
        k01 = g1(y[i], v[i], t[i])
        k02 = g2(y[i], v[i], t[i])
        k11 = g1(y[i]+dt/2*k01, v[i]+dt/2*k02, t[i]+dt/2)
        k12 = g2(y[i]+dt/2*k01, v[i]+dt/2*k02, t[i]+dt/2)
        k21 = g1(y[i]+dt/2*k11, v[i]+dt/2*k12, t[i]+dt/2)
        k22 = g2(y[i]+dt/2*k11, v[i]+dt/2*k12, t[i]+dt/2)
        k31 = g1(y[i]+dt*k21, v[i]+dt*k22, t[i]+dt)
        k32 = g2(y[i]+dt*k21, v[i]+dt*k22, t[i]+dt)
        y[i+1] = y[i] + dt/6 * (k01+2*k11+2*k21+k31)
        v[i+1] = v[i] + dt/6 * (k02+2*k12+2*k22+k32)
    return t, y

# Analytical
t_anal, y_anal = anl(0,6,math.radians(10),0)
    
# Explicit  
fig, ax = plt.subplots(2,2)
fig.suptitle('Explicit Euler Method with Analytical Solution')
ax[0,0].plot(t_anal,y_anal,'-r') 
ax[0,1].plot(t_anal,y_anal,'-r') 
ax[1,0].plot(t_anal,y_anal,'-r') 
ax[0,0].set_title("dt = 0.15s")
ax[0,1].set_title("dt = 0.5s")
ax[1,0].set_title("dt = 1s")
t, y = exp(0,6,math.radians(10),0,0.15)
ax[0,0].plot(t,y,'-o')
t, y = exp(0,6,math.radians(10),0,0.5)
ax[0,1].plot(t,y,'-o')
t, y = exp(0,6,math.radians(10),0,1)
ax[1,0].plot(t,y,'-o')

ax[0,0].set_xlabel('t')
ax[0,0].set_ylabel('theta')
ax[0,1].set_xlabel('t')
ax[0,1].set_ylabel('theta')
ax[1,0].set_xlabel('t')
ax[1,0].set_ylabel('theta')

fig.tight_layout()
plt.savefig('images/exp3a.jpg')
plt.show()

# Implicit
fig, ax = plt.subplots(2,2)
fig.suptitle('Implicit Euler Method with Analytical Solution')
ax[0,0].plot(t_anal,y_anal,'-r') 
ax[0,1].plot(t_anal,y_anal,'-r') 
ax[1,0].plot(t_anal,y_anal,'-r') 
ax[0,0].set_title("dt = 0.15s")
ax[0,1].set_title("dt = 0.5s")
ax[1,0].set_title("dt = 1s")
t, y = imp(0,6,math.radians(10),0,0.15)
ax[0,0].plot(t,y,'-o')
t, y = imp(0,6,math.radians(10),0,0.5)
ax[0,1].plot(t,y,'-o')
t, y = imp(0,6,math.radians(10),0,1)
ax[1,0].plot(t,y,'-o')

ax[0,0].set_xlabel('t')
ax[0,0].set_ylabel('theta')
ax[0,1].set_xlabel('t')
ax[0,1].set_ylabel('theta')
ax[1,0].set_xlabel('t')
ax[1,0].set_ylabel('theta')

fig.tight_layout()
plt.savefig('images/imp3a.jpg')
plt.show()

# Runge-Kutta 4
fig, ax = plt.subplots(2,2)
fig.suptitle('Runge-Kutta 4th Order with Analytical Solution')
ax[0,0].plot(t_anal,y_anal,'-r') 
ax[0,1].plot(t_anal,y_anal,'-r') 
ax[1,0].plot(t_anal,y_anal,'-r') 
ax[0,0].set_title("dt = 0.15s")
ax[0,1].set_title("dt = 0.5s")
ax[1,0].set_title("dt = 1s")
t, y = rk4(0,6,math.radians(10),0,0.15)
ax[0,0].plot(t,y,'-o')
t, y = rk4(0,6,math.radians(10),0,0.5)
ax[0,1].plot(t,y,'-o')
t, y = rk4(0,6,math.radians(10),0,1)
ax[1,0].plot(t,y,'-o')

ax[0,0].set_xlabel('t')
ax[0,0].set_ylabel('theta')
ax[0,1].set_xlabel('t')
ax[0,1].set_ylabel('theta')
ax[1,0].set_xlabel('t')
ax[1,0].set_ylabel('theta')

fig.tight_layout()
plt.savefig('images/rk43a.jpg')
plt.show()


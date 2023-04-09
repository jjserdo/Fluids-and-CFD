
"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1 Number 3
"""
import numpy as np
import matplotlib.pyplot as plt
import math

### Number 1 ###
ti = 0
tf = 6
thetai = math.radians(10)
dthetai = 0
# dt = 0.15, 0.5, 1
grav = 9.81
Leng = 0.6

# g(y(t), t)
def g1(theta, dtheta, t):
    return dtheta;
def g2(theta, dtheta, t):
    return -grav/Leng * theta;
# solvers
def anl()
def exp(dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] =  thetai
    v[0] = dthetai
    for i in range(len(t)-1):
        y[i+1] = y[i] + dt * g1(y[i], v[i], t[i]);
        v[i+1] = v[i] + dt * g2(y[i], v[i], t[i]);
    return t, y
def imp(dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] =  thetai
    v[0] = dthetai
    for i in range(len(t)-1):
        # version 1; based on coupled equations
        nume = y[i] + v[i] * dt
        denu = 1 + grav/Leng * dt**2
        y[i+1] = nume / denu
        v[i+1] = (y[i+1] - y[i]) / dt;
    return t, y
def rk4(dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] =  thetai
    v[0] = dthetai
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
    

# Analytical Solution solution
t_anal = np.linspace(ti,tf,1000)
theta_anal = np.zeros(len(t_anal))
for i in range(len(t_anal)):    
    theta_anal[i] = thetai * math.cos(t_anal[i] * math.sqrt(grav/Leng))

fig1, ax1 = plt.subplots()
ax1.plot(t_anal,theta_anal,'-r',label='Analytical Solution') 
t, theta = imp(0.005) # imp(0.01) is good but hm
ax1.plot(t,theta,'-bo',label='Explicit Euler')
ax1.set_title('Solution with Analytical Solution')
ax1.set_xlabel('Time')
ax1.set_ylabel('y')
ax1.set_xlim([ti,tf])
#ax1.set_ylim([0, 4])
plt.show()



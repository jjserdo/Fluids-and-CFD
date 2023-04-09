"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1 Number 4
"""
import numpy as np
import matplotlib.pyplot as plt
import math

### Number 1 ###
sigm = 10
b = 8/3.

def gx(x, y, z):
    return sigm*(y-x);
def gy(x, y, z):
    return r*x - y  - x*z;
def gz(x, y, z):
    return x*y - b*z;
    
# solvers
def rk4(ti,tf,xi,yi,zi,r,dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    x = np.zeros(len(t))
    y = np.zeros(len(t))
    z = np.zeros(len(t))
    x[0] = xi
    y[0] = yi
    z[0] = zi
    for i in range(len(t)-1):
         k0x = gx(x[i], y[i], z[i])
         k0y = gy(x[i], y[i], z[i])
         k0z = gz(x[i], y[i], z[i])
         k1x = gx(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z)
         k1y = gy(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z)
         k1z = gz(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z)
         k2x = gx(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z)
         k2y = gy(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z)
         k2z = gz(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z)
         k3x = gx(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z)
         k3y = gy(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z)
         k3z = gz(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z)
         x[i+1] = x[i] + dt/6 * (k0x+2*k1x+2*k2x+k3x)
         y[i+1] = y[i] + dt/6 * (k0y+2*k1y+2*k2y+k3y)
         z[i+1] = z[i] + dt/6 * (k0z+2*k1z+2*k2z+k3z)
    return t, x, y, z

# 4 a
fig, ax = plt.subplots(2,2)
#fig.suptitle('Runge-Kutta 4th Order with Analytical Solution')
ax[0,0].set_title("X Y")
ax[0,1].set_title("X Z")
ax[1,0].set_title("Y Z")
t, y = rk4(0,6,math.radians(10),0,0.15)
ax[0,0].plot(t,y,'-o')
t, y = rk4(0,6,math.radians(10),0,0.5)
ax[0,1].plot(t,y,'-o')
t, y = rk4(0,6,math.radians(10),0,1)
ax[1,0].plot(t,y,'-o')
fig.tight_layout()
plt.savefig('images/rk43b.jpg')
plt.show()


t, x, y, z = function4(1, 1, 1)
fig1, ax1 = plt.subplots()
ax1.plot(x, y,'-r',label='x y') 
ax1.set_title('X Y')
ax1.set_xlabel('X')
ax1.set_ylabel('y')
#ax1.set_xlim([0,6])
#ax1.set_ylim([0, 4])
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(x, z,'-r',label='x z') 
ax2.set_title('X Z')
ax2.set_xlabel('X')
ax2.set_ylabel('Z')
#ax1.set_xlim([0,6])
#ax1.set_ylim([0, 4])
plt.show()

fig3, ax3 = plt.subplots()
ax3.plot(y, z,'-r',label='y z') 
ax3.set_title('Y Z')
ax3.set_xlabel('Y')
ax3.set_ylabel('Z')
#ax1.set_xlim([0,6])
#ax1.set_ylim([0, 4])
plt.show()

##########
t1, x1, y1, z1 = function4(6, 6, 6)
t2, x2, y2, z2 = function4(6, 6.01, 6)

fig1c, ax1c = plt.subplots()
ax1c.plot(x1, y1,'-r',label='x y 6') 
ax1c.plot(x2, y2,'-b',label='x y 6.01') 
ax1c.set_title('X Y')
ax1c.set_xlabel('X')
ax1c.set_ylabel('y')
ax1c.legend()
plt.show()

fig2c, ax2c = plt.subplots()
ax2c.plot(x1, z1,'-r',label='x z 6') 
ax2c.plot(x2, z2,'-b',label='x z 6.01') 
ax2c.set_title('X Z')
ax2c.set_xlabel('X')
ax2c.set_ylabel('Z')
ax2c.legend()
plt.show()

fig3c, ax3c = plt.subplots()
ax3c.plot(y1, z1,'-r',label='y z 6') 
ax3c.plot(y2, z2,'-b',label='y z 6.01')
ax3c.set_title('Y Z')
ax3c.set_xlabel('Y')
ax3c.set_ylabel('Z')
ax3c.legend()
plt.show()

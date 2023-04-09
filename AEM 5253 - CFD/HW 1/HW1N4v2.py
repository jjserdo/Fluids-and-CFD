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

def gx(x, y, z, r):
    return sigm*(y-x);
def gy(x, y, z, r):
    return r*x - y  - x*z;
def gz(x, y, z, r):
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
         k0x = gx(x[i], y[i], z[i], r)
         k0y = gy(x[i], y[i], z[i], r)
         k0z = gz(x[i], y[i], z[i], r)
         k1x = gx(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z, r)
         k1y = gy(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z, r)
         k1z = gz(x[i]+dt/2*k0x, y[i]+dt/2*k0y, z[i]+dt/2*k0z, r)
         k2x = gx(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z, r)
         k2y = gy(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z, r)
         k2z = gz(x[i]+dt/2*k1x, y[i]+dt/2*k1y, z[i]+dt/2*k1z, r)
         k3x = gx(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z, r)
         k3y = gy(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z, r)
         k3z = gz(x[i]+dt*k2x, y[i]+dt*k2y, z[i]+dt*k2z, r)
         x[i+1] = x[i] + dt/6 * (k0x+2*k1x+2*k2x+k3x)
         y[i+1] = y[i] + dt/6 * (k0y+2*k1y+2*k2y+k3y)
         z[i+1] = z[i] + dt/6 * (k0z+2*k1z+2*k2z+k3z)
    return t, x, y, z

# 4 a
fig, ax = plt.subplots(2,2)
fig.suptitle('Lorenze Attractor for r = 20')
ax[0,0].set_title("X Y")
ax[0,1].set_title("X Z")
ax[1,0].set_title("Y Z")
t, x, y, z = rk4(0,25,1,1,1,20,0.005)
ax[0,0].plot(x,y,'-')
ax[0,1].plot(x,z,'-')
ax[1,0].plot(y,z,'-')

ax[0,0].set_xlabel('X')
ax[0,0].set_ylabel('Y')
ax[0,1].set_xlabel('X')
ax[0,1].set_ylabel('Z')
ax[1,0].set_xlabel('Y')
ax[1,0].set_ylabel('Z')

fig.tight_layout()
plt.savefig('images/rk44a.jpg')
plt.show()

# 4 b
fig, ax = plt.subplots(2,2)
fig.suptitle('Lorenze Attractor for r = 28')
ax[0,0].set_title("X Y")
ax[0,1].set_title("X Z")
ax[1,0].set_title("Y Z")
t, x, y, z = rk4(0,25,1,1,1,28,0.005)
ax[0,0].plot(x,y,'-')
ax[0,1].plot(x,z,'-')
ax[1,0].plot(y,z,'-')

ax[0,0].set_xlabel('X')
ax[0,0].set_ylabel('Y')
ax[0,1].set_xlabel('X')
ax[0,1].set_ylabel('Z')
ax[1,0].set_xlabel('Y')
ax[1,0].set_ylabel('Z')

fig.tight_layout()
plt.savefig('images/rk44b.jpg')
plt.show()

# 4 c
fig, ax = plt.subplots(2,2)
fig.suptitle('Lorenze Attractor for r = 28 for barely different starting points')
ax[0,0].set_title("X Y")
ax[0,1].set_title("X Z")
ax[1,0].set_title("Y Z")
t, x, y, z = rk4(0,25,6,6,6,28,0.005)
ax[0,0].plot(x,y,'-', label='y=6')
ax[0,1].plot(x,z,'-', label='y=6')
ax[1,0].plot(y,z,'-', label='y=6')
t, x, y, z = rk4(0,25,6,6.01,6,28,0.005)
ax[0,0].plot(x,y,'-', label='y=6.01')
ax[0,1].plot(x,z,'-', label='y=6.01')
ax[1,0].plot(y,z,'-', label='y=6.01')

ax[0,0].set_xlabel('X')
ax[0,0].set_ylabel('Y')
ax[0,1].set_xlabel('X')
ax[0,1].set_ylabel('Z')
ax[1,0].set_xlabel('Y')
ax[1,0].set_ylabel('Z')

ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
fig.tight_layout()
plt.savefig('images/rk44c.jpg')
plt.show()
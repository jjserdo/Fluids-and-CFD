# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 14:28:18 2022

@author: jjser
"""
### first method, does not work
             k0 = gx(x[i],y[i])
             k1 = gx(x[i]+dt/2*k0, y[i]+dt/2*k0)
             k2 = gx(x[i]+dt/2*k1, y[i]+dt/2*k1)
             k3 = gx(x[i]+dt  *k2, y[i]+dt  *k2)
             x[i+1] = x[i] + dt/6 * (k0+2*k1+2*k2+k3)
             k0 = gy(x[i],y[i],z[i])
             k1 = gy(x[i]+dt/2*k0, y[i]+dt/2*k0, z[i]+dt/2*k0)
             k2 = gy(x[i]+dt/2*k1, y[i]+dt/2*k1, z[i]+dt/2*k1)
             k3 = gy(x[i]+dt  *k2, y[i]+dt  *k2, z[i]+dt  *k2)
             y[i+1] = x[i] + dt/6 * (k0+2*k1+2*k2+k3)
             k0 = gz(x[i],y[i],z[i])
             k1 = gz(x[i]+dt/2*k0, y[i]+dt/2*k0, z[i]+dt/2*k0)
             k2 = gz(x[i]+dt/2*k1, y[i]+dt/2*k1, z[i]+dt/2*k1)
             k3 = gz(x[i]+dt  *k2, y[i]+dt  *k2, z[i]+dt  *k2)
             z[i+1] = z[i] + dt/6 * (k0+2*k1+2*k2+k3)
         
         
# third method, best guess that it will work      
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
         
    
 # second method, not sure yet if it works        
            k0 = gx(x[i], y[i], z[i])
            k1 = gx(x[i]+dt/2*k0, y[i], z[i])
            k2 = gx(x[i]+dt/2*k1, y[i], z[i])
            k3 = gx(x[i]+dt  *k2, y[i], z[i])
            x[i+1] = x[i] + dt/6 * (k0+2*k1+2*k2+k3)
            k0 = gy(x[i], y[i], z[i])
            k1 = gy(x[i], y[i]+dt/2*k0, z[i])
            k2 = gy(x[i], y[i]+dt/2*k1, z[i])
            k3 = gy(x[i], y[i]+dt  *k2, z[i])
            y[i+1] = x[i] + dt/6 * (k0+2*k1+2*k2+k3)
            k0 = gz(x[i], y[i], z[i])
            k1 = gz(x[i], y[i], z[i]+dt/2*k0)
            k2 = gz(x[i], y[i], z[i]+dt/2*k1)
            k3 = gz(x[i], y[i], z[i]+dt  *k2)
            z[i+1] = z[i] + dt/6 * (k0+2*k1+2*k2+k3)
            
            
def rk4(dt):
    N = int((tf-ti)/dt)
    t = np.linspace(ti,tf,N+1)
    y = np.zeros(len(t))
    v = np.zeros(len(t))
    y[0] = thetai
    for i in range(len(t)-1):
        k0 = g(y[i],t[i])
        k1 = g(y[i]+dt/2*k0,t[i]+dt/2)
        k2 = g(y[i]+dt/2*k1,t[i]+dt/2)
        k3 = g(y[i]+dt  *k2,t[i]+dt  )
        y[i+1] = y[i] + dt/6 * (k0+2*k1+2*k2+k3)
    return t, y
    
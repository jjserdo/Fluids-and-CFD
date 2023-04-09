# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 10:06:00 2023

@author: jjser
"""
# Imported from CFD Project
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la
from numpy import asarray
import pandas
import time
##################
start_time = time.time()
func = 'z'
s = 10
g = 500

x = np.linspace(-s,s,g)
y = np.linspace(-s,s,g)
X,Y = np.meshgrid(x,y)
R = np.arctan2(Y,X)
T = np.sqrt(X**2+Y**2)

## F(z) = z = x + iy
if func == 'z': # 0.2049 seconds
    phi = X
    psi = Y
elif func == 'z**2'
    phi = 
    psi = 

fig, ax = plt.subplots()
ax.set_title('Complex Plots')
ax.set_xlabel('Real')
ax.set_ylabel('Imaginary')
ax.set_xlim([-s,s])
ax.set_ylim([-s,s])
levels = np.linspace(np.min(phi),np.max(phi),20)
ax.contour(X, Y, phi, levels=levels)
levels = np.linspace(np.min(phi),np.max(phi),20)
ax.contour(X, Y, psi, levels=levels)

print("--- %10s seconds ---" % np.round((time.time() - start_time),4))

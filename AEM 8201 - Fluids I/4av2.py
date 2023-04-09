# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:22:05 2022

@author: jjser
"""

import numpy as np
from matplotlib import pyplot as plt
import math

A = 1
rf = 10
zf = 10

x = np.linspace(-20,20,100)
y = np.linspace(-20,20,100)

a = 10

xg, yg = np.meshgrid(x,y)

psi = np.zeros((len(x),len(y)))

for i in range(len(x)):
    for j in range(len(y)):
        R2 = xg[i,j]**2 + yg[i,j]**2
        psi[i,j] = A*a**4/10 * (R2/a**2) * (1-(R2/a**2))
        
fig, ax = plt.subplots()
ax.contour(xg, yg, psi) 
ax.set_title('Streamlines of the flow')
ax.set_xlabel('x')
ax.set_ylabel('y')
        
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:22:05 2022

@author: jjser
"""

import numpy as np
from matplotlib import pyplot as plt
import math

A = 1
rf = 15
zf = 15
R = np.linspace(-rf,rf,100)
z = np.linspace(-zf,zf,100)
a = 5 #math.sqrt(rf**2 + zf**2)
rg, zg = np.meshgrid(R,z)

psi = np.zeros((len(R),len(z)))

for i in range(len(R)):
    for j in range(len(z)):
        psi[i,j] = A*a**4/10 * (rg[i,j]**2/a**2) * (1-(rg[i,j]**2/a**2)-zg[i,j]**2/a**2)
        
fig, ax = plt.subplots()
#levels = np.logspace(0, 15, num=10)
levels = [-500,-100,-10,-5,-0.5,-0.1,0,0.5,5,10]
ax.contour(rg, zg, psi, levels = levels) 
ax.set_title('Streamlines of the flow')
ax.set_xlabel('R')
ax.set_ylabel('z')
        
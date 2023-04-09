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
R = np.linspace(0,rf,100)
z = np.linspace(-zf,zf,100)
a = math.sqrt(rf**2 + zf**2)
rg, zg = np.meshgrid(R,z)

psi = np.zeros((len(R),len(z)))

for i in range(len(R)):
    for j in range(len(z)):
        psi[i,j] = A*a**4/10 * (rg[i,j]**2/a**2) * (1-(rg[i,j]**2/a**2)-zg[i,j]**2/a**2)
        
fig, ax = plt.subplots()
ax.contour(rg, zg, psi) 
ax.set_title('Streamlines of the flow')
ax.set_xlabel('R')
ax.set_ylabel('z')
        
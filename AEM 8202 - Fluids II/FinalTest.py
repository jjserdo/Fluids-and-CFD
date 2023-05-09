# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:52:50 2023

@author: jjser
Updates:
    - testing the z, eta, F(z) method from Nichols
"""

# Imported from CFD Project
import numpy as np
import time
from scipy.stats import norm 
import matplotlib.pyplot as plt

s = 3
X = np.linspace(-s,s,1000)
Y = np.linspace(-s,s,1000)
x,y = np.meshgrid(X,Y)
r = np.arctan2(y,x)
t = np.sqrt(x**2+y**2)
z = x + 1j*y

m = 0.1
b = 1
a = 1 + m 
num = 10

l = 1/2 * (z + np.sqrt(z+2) * np.sqrt(z-2) )
Fz = (l+m) + b**2/(l+m)

psi = Fz.imag

fig, ax = plt.subplots()
mean = np.mean(psi)
std = np.std(psi)
#levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
ax.contour(x, y, psi)#, levels=levels)